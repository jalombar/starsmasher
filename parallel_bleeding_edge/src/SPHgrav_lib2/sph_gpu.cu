// The SPH of StarSmasher on the GPUs (input option ngpusph, see init.f).
//
// This file is compiled with --fmad=false (SPHgrav_lib2/Makefile): every floating-point operation that decides what the CPU
// code of rho_and_h (advance.f90) decides is made here as one multiplication or addition at a time, in the same order, so that
// the results are the bits the CPU produces.  The gravity files are not compiled that way.
//
// The neighbour lists.  The squared distance of two particles is the sum (dx*dx+dy*dy)+dz*dz in double precision, as in the
// kd tree of the CPU (kdtree2.f90) and in the sums over the neighbours.  A particle's list is every particle at a squared
// distance below r2 = 4*h*h, in the order of the particle index: whatever finds the candidates, that set in that order is
// what the CPU makes, so the search only has to be complete.  The search is made on a tree of boxes built here for each call
// (it is not the tree of the CPU, which is not built when the GPU does the work):
//   1. the particles are sorted along a Morton curve (a radix sort of 63-bit keys: exact and the same on every run),
//   2. the leaves hold LEAF consecutive particles of that order, a node holds the boxes of BR nodes of the level below,
//   3. a thread takes a particle of the rank, walks the tree with a test of the boxes and of the points in single precision
//      with a margin (it can only let too many through, never drop one), and tests what is left with the exact double
//      precision arithmetic.  The candidates are put in the order of their index (a sort in the thread).
// A particle with more than CAP candidates is left to the CPU (flag), by a rule that depends on the particles only.
//
// The particles are done in chunks (SPHGPU_CHUNK, 131072 by default): while the kernel of a chunk runs, the results of the chunk before
// are copied to the CPU, the lists one after the other into list(), and spread over the arrays of the CPU.
//
// The smoothing length.  The thread then solves the smoothing length of its particle with the candidates in its hands: the same
// sequence of trial values, bracketing, Newton and bisection steps, limits, tolerance and sums as solve_particle, bfun,
// bonetsum_d2 and bonetsumwithphi of advance.f90, in the same order of the arithmetic (the sums run over the candidates in the
// order of the index).  The candidates are searched again, further out, when a trial radius goes past what was searched.  The
// tables of the kernel are those of the CPU (kernels.f), read with the same linear interpolation.
//
// Nothing here changes the values of any array of the CPU but the ones the caller names; a call that does not return 0 has
// changed nothing, and the caller does the work on the CPU.
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <cmath>
#include <algorithm>
#include <vector>
#include <utility>
#include <chrono>
#include <stdexcept>
#include <cuda_runtime.h>
#include <cub/cub.cuh>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace
{

// (double4 is deprecated in CUDA 13: a vector type of our own)
struct __align__(16) D4
{
  double x, y, z, w;
};
__host__ __device__ inline D4 make_d4(double x, double y, double z, double w)
{
  D4 r;
  r.x = x; r.y = y; r.z = z; r.w = w;
  return r;
}

// -DSPH_CLOCKS: the cycles spent in the parts of the kernels (lane 0 of each warp), counted in g_clk, printed with SPHGPU_PROF=2
#ifdef SPH_CLOCKS
__device__ unsigned long long g_clk[16];
#define CLK(v) const long long v = clock64()
#define CLKADD(i, a, b) do { if ((threadIdx.x & 31) == 0) atomicAdd(&g_clk[i], (unsigned long long)((b) - (a))); } while (0)
#else
#define CLK(v)
#define CLKADD(i, a, b) do { } while (0)
#endif

constexpr int CAP = 128;      // most candidates of one particle on the GPU
#ifndef SPH_LEAF
#define SPH_LEAF 8
#endif
constexpr int LEAF = SPH_LEAF;   // particles of a leaf
constexpr int BR = 4;         // nodes of a level under a node
constexpr int MAXLEV = 24;
constexpr int NTHR = 128;     // threads of a block of the solve kernel
constexpr int NBLK = 4;
constexpr double GROW = 1.04;      // the first search reaches this factor past the radius asked for (the next ones too)
constexpr int MAXEVAL = 4000;      // evaluations of the function of the solve after which the particle is left to the CPU
constexpr int LSH = 26;       // a stack entry is level << LSH | index
constexpr int LMASK = (1 << LSH) - 1;

// SPHGPU_PROF=1 in the environment: the time of the phases of every call, on stderr (it synchronizes the device between them)
int prof_level()
{
  static const int p = getenv("SPHGPU_PROF") ? atoi(getenv("SPHGPU_PROF")) : 0;
  return p;
}

// SPHGPU_PROF=2: events between the kernels of the first device, printed at the end of the call
struct Mark
{
  const char *name;
  cudaEvent_t ev;
};
std::vector<Mark> g_marks;
void mark(const char *name, cudaStream_t st)
{
  if (prof_level() < 2) return;
  Mark m;
  m.name = name;
  cudaEventCreate(&m.ev);
  cudaEventRecord(m.ev, st);
  g_marks.push_back(m);
}
void print_marks(const char *title)
{
  if (g_marks.empty()) return;
  cudaEventSynchronize(g_marks.back().ev);
  fprintf(stderr, "SPHGPU %s:", title);
  for (size_t k = 1; k < g_marks.size(); k++)
  {
    float ms = 0;
    cudaEventElapsedTime(&ms, g_marks[k - 1].ev, g_marks[k].ev);
    fprintf(stderr, " %s %.3f", g_marks[k].name, ms);
  }
  fprintf(stderr, "\n");
  for (size_t k = 0; k < g_marks.size(); k++) cudaEventDestroy(g_marks[k].ev);
  g_marks.clear();
}
double now_ms()
{
  return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now().time_since_epoch()).count();
}

struct SphError : public std::runtime_error
{
  int code;
  SphError(int c, const std::string &m) : std::runtime_error(m), code(c) {}
};

#define CK(call)                                                                                    \
  do {                                                                                              \
    cudaError_t e_ = (call);                                                                        \
    if (e_ != cudaSuccess) {                                                                        \
      char b_[512];                                                                                 \
      snprintf(b_, sizeof b_, "CUDA error %s at %s:%d: %s", #call, __FILE__, __LINE__, cudaGetErrorString(e_)); \
      throw SphError(e_ == cudaErrorMemoryAllocation ? 4 : 3, b_);                                  \
    }                                                                                               \
  } while (0)

// what a thread writes for each of its particles
struct SphOut
{
  double hp, rho, omega, om0, psi, wn;   // the smoothing length and the sums (om0 is bonet_0mega)
  int idx;    // the particle, from 0
  int nn;     // length of its list
  int first;  // where the list starts in the list of the device
  int flag;   // 1: the GPU did not do it (the CPU must)
};

// what the kernels read (by value)
struct View
{
  int n, nlow, nup;               // particles; the particles of the rank are nlow..nup-1 (from 0)
  const D4 *rec;             // x y z mass of the particles in the order of the Morton curve
  const float4 *frec;             // x y z (single precision, from the corner of the box) and 1 for gas, 0 else
  const int *perm;                // the particle at each place of the curve
  const float4 *box;              // two float4 per node: lower and upper corner
  int nlev;
  int loff[MAXLEV], lcnt[MAXLEV];
  const double *bb;               // corner of the box of the particles (0..2), its largest side (6)
  const double *hp;               // smoothing lengths of the particles of the rank (index i-nlow)
  const double *divv;             // div v of the particles of the rank (index i-nlow)
  const double *am;               // masses (index of the particle)
  const double2 *gd;              // the kernel tables: g and dg/dh
  const D4 *wdp;             // w, dw/dh and dphi/dh
  double ctab, c4, dt, dxfrac, hfloor, nnopt;   // c4 is 4*ctab
  int nselfg;
  int tbase, cn;                  // the chunk: its first particle (of the slice of this device), how many
  int flagmod;                    // test only (SPHGPU_FLAGMOD): leave to the CPU every particle whose index is a multiple of this (0: none)
  const int *tgt;                 // the places of the curve of the particles of the rank
  const int *ntgt;
  int dev, ndev;                  // this device is number dev of ndev that share the particles of the rank
  SphOut *out;
  int *list;                      // where the lists of the chunk go (the places are relative to it)
  int *ltot;                      // the entries written
};

// ---------------------------------------------------------------------------------------------------- the box of the particles
__global__ void k_bbox(int n, const double *x, const double *y, const double *z, double *part)
{
  __shared__ double s[6][256];
  const int t = threadIdx.x;
  double lo[3] = {1e300, 1e300, 1e300}, hi[3] = {-1e300, -1e300, -1e300};
  for (int i = blockIdx.x * blockDim.x + t; i < n; i += gridDim.x * blockDim.x)
  {
    const double p[3] = {x[i], y[i], z[i]};
    for (int k = 0; k < 3; k++)
    {
      lo[k] = p[k] < lo[k] ? p[k] : lo[k];
      hi[k] = p[k] > hi[k] ? p[k] : hi[k];
    }
  }
  for (int k = 0; k < 3; k++) { s[k][t] = lo[k]; s[3 + k][t] = hi[k]; }
  __syncthreads();
  for (int h = 128; h > 0; h >>= 1)
  {
    if (t < h)
      for (int k = 0; k < 3; k++)
      {
        s[k][t] = s[k][t + h] < s[k][t] ? s[k][t + h] : s[k][t];
        s[3 + k][t] = s[3 + k][t + h] > s[3 + k][t] ? s[3 + k][t + h] : s[3 + k][t];
      }
    __syncthreads();
  }
  if (t == 0)
    for (int k = 0; k < 6; k++) part[blockIdx.x * 6 + k] = s[k][0];
}

// bb[0..2] lower corner, bb[3..5] upper corner, bb[6] the largest side, bb[7] the scale of the Morton cells
__global__ void k_bbox_final(int nb, const double *part, double *bb)
{
  __shared__ double s[6][256];
  const int t = threadIdx.x;
  double v[6];
  for (int k = 0; k < 3; k++) { v[k] = 1e300; v[3 + k] = -1e300; }
  for (int b = t; b < nb; b += blockDim.x)
    for (int k = 0; k < 3; k++)
    {
      v[k] = part[b * 6 + k] < v[k] ? part[b * 6 + k] : v[k];
      v[3 + k] = part[b * 6 + 3 + k] > v[3 + k] ? part[b * 6 + 3 + k] : v[3 + k];
    }
  for (int k = 0; k < 6; k++) s[k][t] = v[k];
  __syncthreads();
  for (int h = 128; h > 0; h >>= 1)
  {
    if (t < h)
      for (int k = 0; k < 3; k++)
      {
        s[k][t] = s[k][t + h] < s[k][t] ? s[k][t + h] : s[k][t];
        s[3 + k][t] = s[3 + k][t + h] > s[3 + k][t] ? s[3 + k][t + h] : s[3 + k][t];
      }
    __syncthreads();
  }
  if (t == 0)
  {
    for (int k = 0; k < 6; k++) bb[k] = s[k][0];
    double e = bb[3] - bb[0];
    e = (bb[4] - bb[1]) > e ? (bb[4] - bb[1]) : e;
    e = (bb[5] - bb[2]) > e ? (bb[5] - bb[2]) : e;
    bb[6] = e;
    bb[7] = e > 0.0 ? 2097151.0 / e : 0.0;
  }
}

// ---------------------------------------------------------------------------------------------------- the Morton curve
__device__ __forceinline__ unsigned long long spread21(unsigned long long v)
{
  v &= 0x1fffffULL;
  v = (v | (v << 32)) & 0x1f00000000ffffULL;
  v = (v | (v << 16)) & 0x1f0000ff0000ffULL;
  v = (v | (v << 8)) & 0x100f00f00f00f00fULL;
  v = (v | (v << 4)) & 0x10c30c30c30c30c3ULL;
  v = (v | (v << 2)) & 0x1249249249249249ULL;
  return v;
}

__global__ void k_morton(int n, const double *x, const double *y, const double *z, const double *bb, unsigned long long *key, int *idx)
{
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= n) return;
  const double s = bb[7];
  unsigned long long c[3];
  const double p[3] = {x[i] - bb[0], y[i] - bb[1], z[i] - bb[2]};
  for (int k = 0; k < 3; k++)
  {
    double v = p[k] * s;
    v = v < 0.0 ? 0.0 : (v > 2097151.0 ? 2097151.0 : v);
    c[k] = (unsigned long long)v;
  }
  key[i] = spread21(c[0]) | (spread21(c[1]) << 1) | (spread21(c[2]) << 2);
  idx[i] = i;
}

// the particles in the order of the curve: exact positions, and the single precision copy for the first tests
__global__ void k_gather(int n, const int *perm, const double *x, const double *y, const double *z, const double *am, const double *u,
                         const double *bb, D4 *rec, float4 *frec)
{
  const int k = blockIdx.x * blockDim.x + threadIdx.x;
  if (k >= n) return;
  const int j = perm[k];
  rec[k] = make_d4(x[j], y[j], z[j], am[j]);
  frec[k] = make_float4((float)(x[j] - bb[0]), (float)(y[j] - bb[1]), (float)(z[j] - bb[2]), u[j] != 0.0 ? 1.0f : 0.0f);
}

// ---------------------------------------------------------------------------------------------------- the boxes
__global__ void k_leafbox(int n, int nleaf, const float4 *frec, float4 *box)
{
  const int b = blockIdx.x * blockDim.x + threadIdx.x;
  if (b >= nleaf) return;
  const int p0 = b * LEAF, p1 = min(p0 + LEAF, n);
  float4 lo = frec[p0], hi = frec[p0];
  for (int p = p0 + 1; p < p1; p++)
  {
    const float4 f = frec[p];
    lo.x = fminf(lo.x, f.x); lo.y = fminf(lo.y, f.y); lo.z = fminf(lo.z, f.z);
    hi.x = fmaxf(hi.x, f.x); hi.y = fmaxf(hi.y, f.y); hi.z = fmaxf(hi.z, f.z);
  }
  box[2 * b] = lo;
  box[2 * b + 1] = hi;
}

// the nodes cnt of the level l from the nodes of the level l-1 (cntbelow of them, from offbelow)
__global__ void k_upbox(int cnt, int cntbelow, int offbelow, int off, float4 *box)
{
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= cnt) return;
  const int c0 = BR * i, c1 = min(c0 + BR, cntbelow);
  float4 lo = box[2 * (offbelow + c0)], hi = box[2 * (offbelow + c0) + 1];
  for (int c = c0 + 1; c < c1; c++)
  {
    const float4 l = box[2 * (offbelow + c)], h = box[2 * (offbelow + c) + 1];
    lo.x = fminf(lo.x, l.x); lo.y = fminf(lo.y, l.y); lo.z = fminf(lo.z, l.z);
    hi.x = fmaxf(hi.x, h.x); hi.y = fmaxf(hi.y, h.y); hi.z = fmaxf(hi.z, h.z);
  }
  box[2 * (off + i)] = lo;
  box[2 * (off + i) + 1] = hi;
}

__global__ void k_ownflag(int n, int nlow, int nup, const int *perm, int *flag)
{
  const int k = blockIdx.x * blockDim.x + threadIdx.x;
  if (k >= n) return;
  const int j = perm[k];
  flag[k] = (j >= nlow && j < nup) ? 1 : 0;
}

// the places of the curve of the particles of the rank, in the order of the curve, and their number (scan: the exclusive sum of flag)
__global__ void k_ownlist(int n, const int *flag, const int *scan, int *tgt, int *cnt)
{
  const int k = blockIdx.x * blockDim.x + threadIdx.x;
  if (k >= n) return;
  if (flag[k]) tgt[scan[k]] = k;
  if (k == n - 1) cnt[0] = scan[k] + flag[k];
}

// ---------------------------------------------------------------------------------------------------- the search of one particle
// The candidates of the particle at (xi,yi,zi) found from the root of the tree: every particle at a squared distance not above Rc2,
// with the squared distance of the CPU.  They go to ck (the index of the particle times 2, plus 1 for a gas particle) and cd, element
// k at [k*st], in the order of the index.  Returns their number, or -1 if there are more than CAP.
__device__ int collect(const View &V, double xi, double yi, double zi, double Rc2, unsigned *ck, double *cd, int st)
{
  if (!(Rc2 < 1e300)) return -1;
  const double R = sqrt(Rc2);
  // the single precision tests can only let too many through: the radius has a margin for the rounding of the positions
  // (relative to the largest side of the box) and of the arithmetic
  const float Rf = (float)(R * (1.0 + 2.0e-4) + 1.0e-6 * V.bb[6]);
  const float Rf2 = Rf * Rf;
  const float qx = (float)(xi - V.bb[0]), qy = (float)(yi - V.bb[1]), qz = (float)(zi - V.bb[2]);
  int stk[3 * MAXLEV + 4];
  int sp = 0, m = 0;
  stk[sp++] = (V.nlev - 1) << LSH;
  while (sp > 0)
  {
    const int e = stk[--sp];
    const int l = e >> LSH, i = e & LMASK;
    if (l == 0)
    {
      const int p0 = i * LEAF, p1 = min(p0 + LEAF, V.n);
      for (int p = p0; p < p1; p++)
      {
        const float4 f = V.frec[p];
        const float dx = f.x - qx, dy = f.y - qy, dz = f.z - qz;
        if ((dx * dx + dy * dy) + dz * dz <= Rf2)
        {
          const D4 r = V.rec[p];
          const double ex = xi - r.x, ey = yi - r.y, ez = zi - r.z;
          const double d2 = (ex * ex + ey * ey) + ez * ez;
          if (d2 <= Rc2)
          {
            if (m == CAP) return -1;
            ck[m * st] = ((unsigned)V.perm[p] << 1) | (f.w > 0.5f ? 1u : 0u);
            cd[m * st] = d2;
            m++;
          }
        }
      }
    }
    else
    {
      const int c0 = BR * i, c1 = min(c0 + BR, V.lcnt[l - 1]);
      for (int c = c0; c < c1; c++)
      {
        const int nd = V.loff[l - 1] + c;
        const float4 lo = V.box[2 * nd], hi = V.box[2 * nd + 1];
        const float gx = fmaxf(fmaxf(lo.x - qx, qx - hi.x), 0.0f);
        const float gy = fmaxf(fmaxf(lo.y - qy, qy - hi.y), 0.0f);
        const float gz = fmaxf(fmaxf(lo.z - qz, qz - hi.z), 0.0f);
        if ((gx * gx + gy * gy) + gz * gz <= Rf2) stk[sp++] = ((l - 1) << LSH) | c;
      }
    }
  }
  // By index (insertion: the lists have a few tens of entries).  Other ways were tried: the rank of each key found by comparing it
  // with all the others, and a walk of the tree by the 32 threads of a warp together for a group of particles, the points sorted once
  // and each thread taking its own from the sorted list; neither was faster than this in the time of the whole kernel.
  for (int a = 1; a < m; a++)
  {
    const unsigned k = ck[a * st];
    const double d = cd[a * st];
    int b = a - 1;
    while (b >= 0 && ck[b * st] > k) { ck[(b + 1) * st] = ck[b * st]; cd[(b + 1) * st] = cd[b * st]; b--; }
    ck[(b + 1) * st] = k;
    cd[(b + 1) * st] = d;
  }
  return m;
}

// ---------------------------------------------------------------------------------------------------- the smoothing length
struct Cand
{
  unsigned *ck;    // the candidates in the order of the index: 2*index + 1 if the particle is gas
  double *cd;      // their squared distances; element k at [k*st]
  int st;
  int m;
  double Rc2;      // the candidates are every particle at a squared distance up to this
  int nrec, nev;   // searches again from the root, evaluations of the function (for the counters of -DSPH_CLOCKS)
};

// the candidates again from the root of the tree, out to a squared radius r2 (with the margin); false if there are too many
__device__ bool recollect(const View &V, double xi, double yi, double zi, double r2, Cand &C)
{
  const double Rn = r2 * (GROW * GROW);
  const int m = collect(V, xi, yi, zi, Rn, C.ck, C.cd, C.st);
  C.nrec++;
  if (m < 0) return false;
  C.m = m;
  C.Rc2 = Rn;
  return true;
}

// bfun of advance.f90 (bonetsum_d2 on the candidates): the function of the solve, wn-nnopt, and its derivative, 0mega, for the
// smoothing length h with 4h^2=r2.  False if the candidates could not be made.
__device__ bool bfun(const View &V, double xi, double yi, double zi, Cand &C, double r2, double h, double &f, double &df)
{
  if (r2 > C.Rc2 && !recollect(V, xi, yi, zi, r2, C)) return false;
  C.nev++;
  const double h2 = h * h;
  const double ctoh2 = V.ctab / h2;
  double wn = 0.0, om0 = 0.0;
  for (int k = 0; k < C.m; k++)
  {
    const double d2 = C.cd[k * C.st];
    if (d2 < r2 && (C.ck[k * C.st] & 1u))
    {
      double xq = ctoh2 * d2;
      xq = xq < V.c4 ? xq : V.c4;
      const int it = (int)xq;
      const double fq = xq - (double)it;
      const double2 g0 = V.gd[it], g1 = V.gd[it + 1];
      wn = wn + (g0.x + fq * (g1.x - g0.x));
      om0 = om0 - (g0.y + fq * (g1.y - g0.y));
    }
  }
  om0 = om0 / h;
  f = wn - V.nnopt;
  df = om0;
  return true;
}

// The root of the function of the solve, starting from the guess of the CPU code, with the steps of solve_particle: the bracketing
// by drifting (label 10), then Newton steps kept inside the bracket by bisection (label 12).  The smoothing length without the
// floor hfloor is returned in h.  False if the candidates could not be made or the solve went on too long.
__device__ bool solve_h(const View &V, double xi, double yi, double zi, Cand &C, double hin, double divvi, double &hsol)
{
  const double h0 = hin - V.hfloor;
  const double xacc = 2.0e-6 * h0;
  const double a1 = h0 * (1.0 + (divvi * V.dt) / 3.0);
  const double hpguess = a1 > xacc ? a1 : xacc;
  const double b1 = V.dxfrac * h0, b2 = fabs((h0 * divvi) * V.dt);
  const double b3 = b1 > b2 ? b1 : b2, b4 = 0.5 * hpguess;
  double dxmax = b3 < b4 ? b3 : b4;
  double xl = hpguess;
  dxmax = fabs(dxmax);
  double h = xl, r2, fl, fh, ff, df, dx, ezr, xh, xtmp;
  int neval = 0;
  r2 = 4.0 * (h * h);
  if (!bfun(V, xi, yi, zi, C, r2, h, fl, df)) return false;
  neval++;
  dx = fl / df;
  if (fabs(dx) < xacc) { hsol = xl - dx; return true; }
  if (fl == 0.0) { hsol = xl; return true; }
  {
    const double c1 = fabs(dx) < dxmax ? fabs(dx) : dxmax;
    const double c2 = fabs(0.5 * xl);
    dx = copysign(c1 < c2 ? c1 : c2, dx);
  }
  // label 10: the root is not yet bracketed, drift a little more
  for (;;)
  {
    xh = xl - dx;
    h = xh;
    r2 = 4.0 * (h * h);
    if (!bfun(V, xi, yi, zi, C, r2, h, fh, df)) return false;
    if (++neval > MAXEVAL) return false;
    dx = fh / df;
    if (fabs(dx) < xacc) { hsol = xh - dx; return true; }
    if (fh == 0.0) { hsol = xh; return true; }
    dx = copysign(fabs(dx) < dxmax ? fabs(dx) : dxmax, dx);
    ezr = xh - dx;
    if (fl * fh > 0.0)
    {
      dx = copysign(dx, fh);
      xl = xh;
      continue;
    }
    break;
  }
  // the root is bracketed: xl must give a negative value
  if (fl > 0.0) { xtmp = xh; xh = xl; xl = xtmp; }
  // label 12
  for (;;)
  {
    for (int j = 1; j <= 50; j++)
    {
      if ((ezr - xh) * (ezr - xl) >= 0.0)
      {
        dx = 0.5 * (xh - xl);
        ezr = xl + dx;
      }
      if (fabs(dx) < xacc) { hsol = ezr; return true; }
      h = ezr;
      r2 = 4.0 * (h * h);
      if (!bfun(V, xi, yi, zi, C, r2, h, ff, df)) return false;
      if (++neval > MAXEVAL) return false;
      if (ff < 0.0) xl = ezr; else xh = ezr;
      dx = ff / df;
      ezr = ezr - dx;
    }
    // the maximum number of iterations: the CPU code writes a warning and starts again from the middle of the bracket
    dx = 0.5 * (xh - xl);
    ezr = xl + dx;
  }
}

// The smoothing length of every particle of the rank, its density and the other sums over the neighbours (bonetsumwithphi), and its
// neighbour list, as solve_particle makes them.  A point particle (u=0, dynhco=0) keeps its smoothing
// length and gets the list and the density.
__global__ void __launch_bounds__(NTHR, NBLK) k_solve(View V)
{
  // the places of the lists are found by the warp, not by the block: a warp then waits only for its own threads (the iterations of the
  // solve differ a lot from one particle to the next)
  typedef cub::WarpScan<int> WS;
  __shared__ typename WS::TempStorage ws[NTHR / 32];
  const int nt = *V.ntgt;
  const int t0 = (int)(((long long)nt * V.dev) / V.ndev);
  const int tl = blockIdx.x * NTHR + threadIdx.x;      // the particle in the chunk
  const int tid = t0 + V.tbase + tl;
  const bool valid = tl < V.cn;
  int nn = 0, flag = 0, i = -1;
  double hfin = 0.0, r2f = 0.0, rho = 0.0, omg = 0.0, om0 = 0.0, psi = 0.0, wn = 0.0;
  // the candidates of the thread (local memory: the compiler lays it out so that the threads of a warp read neighbouring words)
  unsigned ck[CAP];
  double cd[CAP];
  Cand C;
  C.ck = ck; C.cd = cd; C.st = 1; C.m = 0; C.Rc2 = -1.0; C.nrec = 0; C.nev = 0;
  CLK(s0);
  if (valid)
  {
    const int k = V.tgt[tid];
    i = V.perm[k];
    const D4 r = V.rec[k];
    const bool gas = V.frec[k].w > 0.5f;
    const double hin = V.hp[i - V.nlow];
    bool ok = !(V.flagmod > 0 && (i + 1) % V.flagmod == 0);
    if (!gas)
    {
      // the list at the smoothing length it has
      hfin = hin;
      r2f = 4.0 * (hfin * hfin);
      if (ok) ok = recollect(V, r.x, r.y, r.z, r2f, C);
    }
    else
    {
      double hs = 0.0;
      // the first search reaches past the guess by the margin, the later ones past the radius asked for
      const double hg = hin - V.hfloor;
      const double hguess0 = hg * (1.0 + (V.divv[i - V.nlow] * V.dt) / 3.0);
      const double hstart = hguess0 > 2.0e-6 * hg ? hguess0 : 2.0e-6 * hg;
      if (ok) ok = recollect(V, r.x, r.y, r.z, 4.0 * (hstart * hstart), C);
      if (ok) ok = solve_h(V, r.x, r.y, r.z, C, hin, V.divv[i - V.nlow], hs);
      CLK(s1);
      CLKADD(5, s0, s1);
      if (ok)
      {
        hfin = hs + V.hfloor;
        r2f = 4.0 * (hfin * hfin);
        if (r2f > C.Rc2) ok = recollect(V, r.x, r.y, r.z, r2f, C);
      }
    }
    if (!ok)
    {
      flag = 1;
      hfin = hin;
    }
    else
    {
      for (int a = 0; a < C.m; a++)
        if (cd[(a)] < r2f) nn++;
      const double hpi = hfin;
      const double h2 = hpi * hpi;
      const double ctoh2 = V.ctab / h2;
      if (!gas)
      {
        // bonetsumwithphi without the sums of the smoothing length (point particle, dynhco=0)
        for (int a = 0; a < C.m; a++)
        {
          const double r2p = cd[(a)];
          if (!(r2p < r2f)) continue;
          if (ck[(a)] & 1u)
          {
            double xq = ctoh2 * r2p;
            xq = xq < V.c4 ? xq : V.c4;
            const int it = (int)xq;
            const double fq = xq - (double)it;
            const D4 w0 = V.wdp[it], w1 = V.wdp[it + 1];
            rho = rho + V.am[ck[(a)] >> 1] * (w0.x + fq * (w1.x - w0.x));
          }
        }
        rho = rho / ((hpi * hpi) * hpi);
        omg = 0.0;
        psi = 0.0;
        om0 = 1.0e30;
      }
      else
      {
        const double hpt = hpi - V.hfloor;
        const double h2t = hpt * hpt;
        const double ctoh2t = V.ctab / h2t;
        const bool samet = V.hfloor == 0.0;
        const double h3 = hpi * h2, h4 = h2 * h2;
        for (int a = 0; a < C.m; a++)
        {
          const double r2p = cd[(a)];
          if (!(r2p < r2f)) continue;
          const unsigned jj = ck[(a)] >> 1;
          const bool gj = (ck[(a)] & 1u) != 0u;
          const double amj = V.am[jj];
          const double xqr = ctoh2 * r2p;
          const double xq = xqr < V.c4 ? xqr : V.c4;
          const int it = (int)xq;
          const double fq = xq - (double)it;
          const D4 w0 = V.wdp[it], w1 = V.wdp[it + 1];
          if (r2p < 4.0 * h2)
          {
            if ((V.nselfg == 1 || !gj) && true)
              psi = psi + amj * (w0.z + fq * (w1.z - w0.z));
          }
          if (gj)
          {
            if (r2p < 4.0 * h2t)
            {
              // without a floor h2t is h2, the index and the fraction of the smoothing length without the floor are the same ones, unless the
              // table was clamped
              int itt = it;
              double fqt = fq;
              if (!(samet && xqr <= V.c4))
              {
                const double xqt = ctoh2t * r2p;
                itt = (int)xqt;
                fqt = xqt - (double)itt;
              }
              const double2 g0 = V.gd[itt], g1 = V.gd[itt + 1];
              wn = wn + (g0.x + fqt * (g1.x - g0.x));
              om0 = om0 - (g0.y + fqt * (g1.y - g0.y));
            }
            if (r2p < 4.0 * h2)
            {
              rho = rho + amj * (w0.x + fq * (w1.x - w0.x));
              omg = omg + amj * (w0.y + fq * (w1.y - w0.y));
            }
          }
        }
        rho = rho / h3;
        omg = omg / h4;
        psi = psi / h2;
        om0 = om0 / hpt;
        if (om0 == 0.0)
        {
          om0 = 1.0e-4;
          omg = 0.0;
          psi = 0.0;
        }
      }
    }
  }
  CLK(s2);
  CLKADD(6, s0, s2);
#ifdef SPH_CLOCKS
  {
    // searches again from the root (beyond the first), evaluations of the function, and the warps in which some thread searched again
    const unsigned any = __ballot_sync(0xffffffffu, C.nrec > 1);
    if ((threadIdx.x & 31) == 0 && valid) { atomicAdd(&g_clk[11], 1ull); if (any) atomicAdd(&g_clk[12], 1ull); }
    if (valid) { atomicAdd(&g_clk[13], (unsigned long long)(C.nrec > 1 ? C.nrec - 1 : 0)); atomicAdd(&g_clk[14], (unsigned long long)C.nev); }
  }
#endif
  int off, total;
  const int mine = (valid && !flag) ? nn : 0;
  WS(ws[threadIdx.x >> 5]).ExclusiveSum(mine, off, total);
  int wbase = 0;
  if ((threadIdx.x & 31) == 0 && total > 0) wbase = atomicAdd(V.ltot, total);
  wbase = __shfl_sync(0xffffffffu, wbase, 0);
  if (valid)
  {
    const int base = wbase + off;
    if (!flag)
    {
      int c = 0;
      for (int a = 0; a < C.m; a++)
        if (cd[(a)] < r2f) V.list[base + c++] = (int)(ck[(a)] >> 1) + 1;   // (the lists of the CPU hold the index from 1, as Fortran numbers the particles)
    }
    SphOut o;
    o.hp = hfin; o.rho = rho; o.omega = omg; o.om0 = om0; o.psi = psi; o.wn = wn;
    o.idx = i; o.nn = flag ? 0 : nn; o.first = base; o.flag = flag;
    V.out[V.tbase + tl] = o;
  }
  CLK(s3);
  CLKADD(7, s2, s3);
}


// ---------------------------------------------------------------------------------------------------- one device
struct Dev
{
  int id = -1;
  cudaStream_t st = 0, st2 = 0;   // st: the build of the tree of boxes; st2: the copies of the results
  cudaStream_t sc[4] = {0, 0, 0, 0};   // the kernels of the chunks, on different streams so that the tail of one is filled by the next
  cudaEvent_t evB = 0;            // the tree of boxes is built
  int ncap = 0, nowncap = 0, chcap = 0;
  double *x = 0, *y = 0, *z = 0, *am = 0, *u = 0, *hp = 0, *divv = 0, *bbpart = 0, *bb = 0;
  double2 *gd = 0;             // the tables of the kernel
  D4 *wdp = 0;
  double ctab = 0;
  int ntab = 0;
  unsigned long long *keya = 0, *keyb = 0;
  int *ida = 0, *idb = 0;
  D4 *rec = 0;
  float4 *frec = 0, *box = 0;
  int *oflag = 0, *oscan = 0;
  int *tgt = 0, *cnt = 0;      // cnt[0]: number of particles of the rank
  SphOut *out = 0;
  int *list = 0, *ccnt = 0;    // the lists of chunk c start at c*chunk*CAP, ccnt[c] entries
  void *tmp = 0;
  size_t tmpbytes = 0;
  int nlev = 0, loff[MAXLEV], lcnt[MAXLEV], nbox = 0;
  SphOut *hout = 0;            // pinned
  int *hcc = 0;                // pinned: the entries of the lists of each chunk
  int houtcap = 0;
  std::vector<cudaEvent_t> evS, evC;   // the kernels of chunk c done / the copies of chunk c done
};

std::vector<Dev> g_dev;
bool g_ready = false;

template <class T> void dfree(T *&p) { if (p) { cudaFree(p); p = 0; } }

void freedev(Dev &d)
{
  cudaSetDevice(d.id);
  dfree(d.x); dfree(d.y); dfree(d.z); dfree(d.am); dfree(d.u); dfree(d.hp); dfree(d.divv); dfree(d.bbpart); dfree(d.bb);
  dfree(d.keya); dfree(d.keyb); dfree(d.ida); dfree(d.idb); dfree(d.rec); dfree(d.frec); dfree(d.box);
  dfree(d.oflag); dfree(d.oscan); dfree(d.tgt); dfree(d.cnt); dfree(d.out); dfree(d.list); dfree(d.ccnt); dfree(d.tmp);
  if (d.hout) { cudaFreeHost(d.hout); d.hout = 0; }
  if (d.hcc) { cudaFreeHost(d.hcc); d.hcc = 0; }
  d.ncap = d.nowncap = d.chcap = 0; d.tmpbytes = 0; d.houtcap = 0;
}

template <class T> void dmalloc(T *&p, size_t count)
{
  CK(cudaMalloc((void **)&p, std::max<size_t>(count, 1) * sizeof(T)));
}

// the particles of a chunk (SPHGPU_CHUNK in the environment changes it)
int chunk_size()
{
  static const int c = getenv("SPHGPU_CHUNK") ? std::max(1024, atoi(getenv("SPHGPU_CHUNK"))) : 131072;
  return c;
}

// room for n particles of which nown are the rank's, in chunks of ch
void reserve(Dev &d, int n, int nown, int ch)
{
  CK(cudaSetDevice(d.id));
  if (n > d.ncap || nown > d.nowncap || ch > d.chcap)
  {
    cudaStreamSynchronize(d.st);
    cudaStreamSynchronize(d.st2);
    for (int q = 0; q < 4; q++) if (d.sc[q]) cudaStreamSynchronize(d.sc[q]);
    const int ncap = std::max(n, d.ncap) + std::max(n, d.ncap) / 16 + 1024;
    const int ocap = std::max(nown, d.nowncap) + std::max(nown, d.nowncap) / 16 + 1024;
    const int chcap = std::max(ch, d.chcap);
    freedev(d);
    CK(cudaSetDevice(d.id));
    dmalloc(d.x, ncap); dmalloc(d.y, ncap); dmalloc(d.z, ncap); dmalloc(d.am, ncap); dmalloc(d.u, ncap);
    dmalloc(d.hp, ocap); dmalloc(d.divv, ocap);
    dmalloc(d.bbpart, 6 * 256); dmalloc(d.bb, 8);
    dmalloc(d.keya, ncap); dmalloc(d.keyb, ncap); dmalloc(d.ida, ncap); dmalloc(d.idb, ncap);
    dmalloc(d.rec, ncap); dmalloc(d.frec, ncap);
    // the boxes: the levels sum to less than n/LEAF*(1+1/BR+...) + a few nodes per level
    dmalloc(d.box, 2 * ((size_t)ncap / LEAF * 4 / 3 + 4 * MAXLEV));
    dmalloc(d.oflag, ncap); dmalloc(d.oscan, ncap); dmalloc(d.tgt, ncap); dmalloc(d.cnt, 4);
    dmalloc(d.out, ocap);
    dmalloc(d.list, ((size_t)ocap + chcap) * CAP);
    dmalloc(d.ccnt, 256);
    size_t b1 = 0, b2 = 0;
    CK(cub::DeviceRadixSort::SortPairs(NULL, b1, (unsigned long long *)0, (unsigned long long *)0, (int *)0, (int *)0, ncap, 0, 63));
    CK(cub::DeviceScan::ExclusiveSum(NULL, b2, (int *)0, (int *)0, ncap));
    d.tmpbytes = std::max(b1, b2) + 256;
    CK(cudaMalloc(&d.tmp, d.tmpbytes));
    CK(cudaMallocHost((void **)&d.hout, (size_t)ocap * sizeof(SphOut)));
    CK(cudaMallocHost((void **)&d.hcc, 256 * sizeof(int)));
    d.houtcap = ocap;
    d.ncap = ncap; d.nowncap = ocap; d.chcap = chcap;
  }
  if (!d.evB) CK(cudaEventCreateWithFlags(&d.evB, cudaEventDisableTiming));
  for (int k = 0; k < 4; k++)
    if (!d.sc[k]) CK(cudaStreamCreateWithFlags(&d.sc[k], cudaStreamNonBlocking));
  while (d.evS.size() < 256)
  {
    cudaEvent_t a, b;
    CK(cudaEventCreateWithFlags(&a, cudaEventDisableTiming));
    CK(cudaEventCreateWithFlags(&b, cudaEventDisableTiming));
    d.evS.push_back(a);
    d.evC.push_back(b);
  }
}

// Pin the arrays of the caller so that the copies are fast; harmless if it cannot be done.  The arrays lie one after the other in
// common blocks, so two of them can share a page, and a page can be registered once: the ranges are made of whole pages, merged
// where they touch, and each merged range is registered the first time it is asked for.
std::vector<std::pair<uintptr_t, uintptr_t> > g_pinned;

void pin_all(const std::vector<std::pair<const void *, size_t> > &arrays)
{
  const uintptr_t page = 4096;
  std::vector<std::pair<uintptr_t, uintptr_t> > r;
  for (size_t k = 0; k < arrays.size(); k++)
  {
    const uintptr_t a = (uintptr_t)arrays[k].first, b = a + arrays[k].second;
    r.push_back(std::make_pair(a / page * page, (b + page - 1) / page * page));
  }
  std::sort(r.begin(), r.end());
  std::vector<std::pair<uintptr_t, uintptr_t> > m;
  for (size_t k = 0; k < r.size(); k++)
  {
    if (!m.empty() && r[k].first <= m.back().second) m.back().second = std::max(m.back().second, r[k].second);
    else m.push_back(r[k]);
  }
  for (size_t k = 0; k < m.size(); k++)
  {
    bool known = false;
    for (size_t j = 0; j < g_pinned.size(); j++)
      if (m[k].first < g_pinned[j].second && g_pinned[j].first < m[k].second) known = true;  // done (or overlapping: leave it)
    if (known) continue;
    g_pinned.push_back(m[k]);
    const cudaError_t e = cudaHostRegister((void *)m[k].first, m[k].second - m[k].first, cudaHostRegisterPortable);
    if (e != cudaSuccess) cudaGetLastError();
  }
}

inline int nblocks(int n, int t) { return (n + t - 1) / t; }

// Copies the particles to the device and builds the tree of boxes there, up to the places of the rank's particles in the
// order of the curve (d.tgt, d.cnt[0]): asynchronous.  hp: the smoothing lengths of the rank's particles (nown of them).
void upload_and_build(Dev &d, int n, int nlow, int nup, const double *x, const double *y, const double *z, const double *am,
                      const double *u, const double *hp, const double *divv)
{
  CK(cudaSetDevice(d.id));
  const int nown = nup - nlow;
  const bool mk = (d.id == g_dev[0].id);
  if (mk) mark("start", d.st);
  CK(cudaMemsetAsync(d.cnt, 0, 4 * sizeof(int), d.st));
  CK(cudaMemcpyAsync(d.x, x, (size_t)n * 8, cudaMemcpyHostToDevice, d.st));
  CK(cudaMemcpyAsync(d.y, y, (size_t)n * 8, cudaMemcpyHostToDevice, d.st));
  CK(cudaMemcpyAsync(d.z, z, (size_t)n * 8, cudaMemcpyHostToDevice, d.st));
  CK(cudaMemcpyAsync(d.am, am, (size_t)n * 8, cudaMemcpyHostToDevice, d.st));
  CK(cudaMemcpyAsync(d.u, u, (size_t)n * 8, cudaMemcpyHostToDevice, d.st));
  CK(cudaMemcpyAsync(d.hp, hp, (size_t)nown * 8, cudaMemcpyHostToDevice, d.st));
  CK(cudaMemcpyAsync(d.divv, divv, (size_t)nown * 8, cudaMemcpyHostToDevice, d.st));
  if (mk) mark("h2d", d.st);
  // the box of the particles and the curve
  const int nbb = std::min(256, nblocks(n, 256));
  k_bbox<<<nbb, 256, 0, d.st>>>(n, d.x, d.y, d.z, d.bbpart);
  k_bbox_final<<<1, 256, 0, d.st>>>(nbb, d.bbpart, d.bb);
  if (mk) mark("bbox", d.st);
  k_morton<<<nblocks(n, 256), 256, 0, d.st>>>(n, d.x, d.y, d.z, d.bb, d.keya, d.ida);
  if (mk) mark("morton", d.st);
  size_t bytes = d.tmpbytes;
  CK(cub::DeviceRadixSort::SortPairs(d.tmp, bytes, d.keya, d.keyb, d.ida, d.idb, n, 0, 63, d.st));
  if (mk) mark("sort", d.st);
  k_gather<<<nblocks(n, 256), 256, 0, d.st>>>(n, d.idb, d.x, d.y, d.z, d.am, d.u, d.bb, d.rec, d.frec);
  if (mk) mark("gather", d.st);
  // the boxes
  const int nl = (n + LEAF - 1) / LEAF;
  int lev = 0, off = 0;
  d.lcnt[0] = nl;
  d.loff[0] = 0;
  while (d.lcnt[lev] > 1)
  {
    off += d.lcnt[lev];
    lev++;
    if (lev >= MAXLEV) throw SphError(3, "too many levels in the tree of boxes");
    d.lcnt[lev] = (d.lcnt[lev - 1] + BR - 1) / BR;
    d.loff[lev] = off;
  }
  d.nlev = lev + 1;
  d.nbox = off + 1;
  k_leafbox<<<nblocks(nl, 128), 128, 0, d.st>>>(n, nl, d.frec, d.box);
  for (int l = 1; l < d.nlev; l++)
    k_upbox<<<nblocks(d.lcnt[l], 128), 128, 0, d.st>>>(d.lcnt[l], d.lcnt[l - 1], d.loff[l - 1], d.loff[l], d.box);
  if (mk) mark("boxes", d.st);
  // the particles of the rank, in the order of the curve
  k_ownflag<<<nblocks(n, 256), 256, 0, d.st>>>(n, nlow, nup, d.idb, d.oflag);
  bytes = d.tmpbytes;
  CK(cub::DeviceScan::ExclusiveSum(d.tmp, bytes, d.oflag, d.oscan, n, d.st));
  k_ownlist<<<nblocks(n, 256), 256, 0, d.st>>>(n, d.oflag, d.oscan, d.tgt, d.cnt);
  if (mk) mark("own", d.st);
  CK(cudaGetLastError());
}

struct Phys { double dt, hfloor, nnopt; int nselfg; };

View makeview(const Dev &d, int ndev, int idev, int n, int nlow, int nup, const Phys &ph)
{
  View V;
  memset(&V, 0, sizeof V);
  V.n = n; V.nlow = nlow; V.nup = nup;
  V.rec = d.rec; V.frec = d.frec; V.perm = d.idb; V.box = d.box;
  V.nlev = d.nlev;
  for (int l = 0; l < d.nlev; l++) { V.loff[l] = d.loff[l]; V.lcnt[l] = d.lcnt[l]; }
  V.hp = d.hp; V.divv = d.divv; V.am = d.am; V.gd = d.gd; V.wdp = d.wdp;
  V.ctab = d.ctab; V.c4 = 4.0 * d.ctab; V.dt = ph.dt; V.dxfrac = ph.dt > 0.0 ? 2.0e-5 : 0.5; V.hfloor = ph.hfloor; V.nnopt = ph.nnopt; V.nselfg = ph.nselfg;
  V.tgt = d.tgt; V.ntgt = d.cnt; V.bb = d.bb;
  V.flagmod = getenv("SPHGPU_FLAGMOD") ? atoi(getenv("SPHGPU_FLAGMOD")) : 0;
  V.dev = idev; V.ndev = ndev;
  V.out = d.out;
  return V;
}

// the particles of the rank that device k of ndev takes (all the rank's particles are targets)
inline int slice_of(int nown, int ndev, int k) { return (int)(((long long)nown * (k + 1)) / ndev - ((long long)nown * k) / ndev); }

// The arrays of the caller that the results go to (Fortran arrays, index from 0 here).
struct Host
{
  double *hp, *rho, *omega, *om0, *psi, *wn, *zeta;
  const double *u;
  int *nn, *first, *sflag;   // nn and first indexed by the particle; sflag by the particle of the rank
  int *list;
  size_t listcap;
};

// Launches the work of one device: for each chunk of its particles the kernel (on one of four streams, so that the tail of one chunk is
// filled by the next), then an event.  Asynchronous.  ev: with SPHGPU_PROF=2 on the first device, the events between the kernels.
void launch_chunks(Dev &d, const View &V0, int slice, int ch, std::vector<cudaEvent_t> *ev)
{
  CK(cudaSetDevice(d.id));
  const int nch = (slice + ch - 1) / ch;
  if (nch > 250) throw SphError(3, "too many chunks");
  CK(cudaMemsetAsync(d.ccnt, 0, 256 * sizeof(int), d.st));
  CK(cudaEventRecord(d.evB, d.st));
  for (int c = 0; c < nch; c++)
  {
    cudaStream_t sc = d.sc[c % 4];
    CK(cudaStreamWaitEvent(sc, d.evB, 0));
    View V = V0;
    V.tbase = c * ch;
    V.cn = std::min(ch, slice - c * ch);
    V.list = d.list + (size_t)c * ch * CAP;
    V.ltot = d.ccnt + c;
    if (ev)
    {
      cudaEvent_t e;
      cudaEventCreate(&e); cudaEventRecord(e, sc); ev->push_back(e);
    }
    k_solve<<<nblocks(V.cn, NTHR), NTHR, 0, sc>>>(V);
    CK(cudaGetLastError());
    if (ev)
    {
      cudaEvent_t e;
      cudaEventCreate(&e); cudaEventRecord(e, sc); ev->push_back(e);
    }
    CK(cudaEventRecord(d.evS[c], sc));
  }
}

// The results, chunk by chunk as the kernels finish them: the number of entries of the lists of the chunk is read, the records and the
// lists are copied to the caller (the lists one after the other in list(1:), whatever the device), and the records of the chunk before
// are spread over the arrays of the caller while the copies and the kernels go on.  The arrays of the caller change as the records are
// spread: on a return with ierr=2 (the lists do not fit) the caller has to put back what it needs of them (the guess of h).
void drain(int nown, int nlow, int ch, const Host &H, int *ltot, int *ierr, long *nflagged)
{
  const int ndev = (int)g_dev.size();
  std::vector<int> slice(ndev), nch(ndev);
  int nchmax = 0;
  for (int k = 0; k < ndev; k++)
  {
    slice[k] = slice_of(nown, ndev, k);
    nch[k] = (slice[k] + ch - 1) / ch;
    nchmax = std::max(nchmax, nch[k]);
  }
  long long hb = 0;
  long nfl = 0;
  bool full = false;
  std::vector<std::vector<long long> > base(ndev, std::vector<long long>(nchmax + 1, 0));
  auto spread = [&](int c) {
    // the records of chunk c of every device, to the arrays of the caller
    for (int k = 0; k < ndev; k++)
    {
      if (c >= nch[k]) continue;
      Dev &d = g_dev[k];
      CK(cudaSetDevice(d.id));
      CK(cudaEventSynchronize(d.evC[c]));
      const int o0 = c * ch, o1 = std::min(slice[k], (c + 1) * ch);
      const long long b0 = base[k][c];
#pragma omp parallel for schedule(static) reduction(+ : nfl)
      for (int o = o0; o < o1; o++)
      {
        const SphOut s = d.hout[o];
        const int i = s.idx;
        H.sflag[i - nlow] = s.flag;
        nfl += s.flag;
        if (s.flag) continue;
        H.hp[i] = s.hp;
        H.rho[i] = s.rho;
        H.omega[i] = s.omega;
        H.om0[i] = s.om0;
        H.psi[i] = s.psi;
        if (H.u[i] != 0.0) H.wn[i] = s.wn;
        H.zeta[i] = 0.0;
        H.nn[i] = s.nn;
        H.first[i] = (int)(b0 + s.first);
      }
    }
  };
  for (int c = 0; c < nchmax; c++)
  {
    for (int k = 0; k < ndev; k++)
    {
      if (c >= nch[k]) continue;
      Dev &d = g_dev[k];
      CK(cudaSetDevice(d.id));
      CK(cudaEventSynchronize(d.evS[c]));
      CK(cudaMemcpyAsync(&d.hcc[c], d.ccnt + c, sizeof(int), cudaMemcpyDeviceToHost, d.st2));
      CK(cudaStreamSynchronize(d.st2));
      const int tot = d.hcc[c];
      const int o0 = c * ch, o1 = std::min(slice[k], (c + 1) * ch);
      base[k][c] = hb;
      if ((size_t)(hb + tot) > H.listcap) full = true;
      if (!full)
      {
        CK(cudaMemcpyAsync(d.hout + o0, d.out + o0, (size_t)(o1 - o0) * sizeof(SphOut), cudaMemcpyDeviceToHost, d.st2));
        CK(cudaMemcpyAsync(H.list + hb, d.list + (size_t)c * ch * CAP, (size_t)tot * 4, cudaMemcpyDeviceToHost, d.st2));
      }
      CK(cudaEventRecord(d.evC[c], d.st2));
      hb += tot;
    }
    if (!full && c > 0) spread(c - 1);
  }
  for (int k = 0; k < ndev; k++)
  {
    CK(cudaSetDevice(g_dev[k].id));
    CK(cudaStreamSynchronize(g_dev[k].st2));
    CK(cudaStreamSynchronize(g_dev[k].st));
    for (int q = 0; q < 4; q++) CK(cudaStreamSynchronize(g_dev[k].sc[q]));
  }
  if (full) { *ierr = 2; return; }
  if (nchmax > 0) spread(nchmax - 1);
  *ltot = (int)hb;
  if (nflagged) *nflagged = nfl;
}

}  // namespace

extern "C"
{
  // 1 if this library can do the SPH on the GPU, 2 if it can but the CPU code of this build fuses multiplications and additions
  // (MARCH, or a processor that does it by default), so that the CPU results would not be the ones of the GPU
  void gpu_sph_available_(int *ok)
  {
#ifdef SS_CPU_FMA
    *ok = 2;
#else
    *ok = 1;
#endif
  }

  // The devices of this process: ngpuperrank of them from localrank*ngpuperrank on (as the gravity does), or if ngpuperrank<=0
  // the devices of the node divided among ngravprocs processes.  ierr: 0 ok, 1 the node has not the devices, 3 error of CUDA.
  void gpu_sph_setup_(int *localrank, int *ngpuperrank, int *ngravprocs, int *ierr)
  {
    *ierr = 0;
    if (g_ready) return;
    try
    {
      int ndevice = 0;
      if (cudaGetDeviceCount(&ndevice) != cudaSuccess || ndevice < 1) { cudaGetLastError(); *ierr = 1; return; }
      int nd = *ngpuperrank;
      if (nd <= 0) nd = std::max(1, ndevice / std::max(1, *ngravprocs));
      const int first = *localrank * nd;
      if (first + nd > ndevice) { *ierr = 1; return; }
      g_dev.resize(nd);
      for (int k = 0; k < nd; k++)
      {
        g_dev[k].id = first + k;
        CK(cudaSetDevice(first + k));
        CK(cudaStreamCreateWithFlags(&g_dev[k].st, cudaStreamNonBlocking));
      }
      CK(cudaSetDevice(first));   // the device the gravity of this process expects to be the current one of the thread
      g_ready = true;
    }
    catch (const SphError &e)
    {
      fprintf(stderr, "StarSmasher (SPH on the GPU): %s\n", e.what());
      *ierr = e.code;
    }
  }

  // number of devices of this process
  void gpu_sph_ndev_(int *nd) { *nd = (int)g_dev.size(); }

  // The tables of the kernel (kernels.f): g, dg/dh, w, dw/dh and dphi/dh with the guard entry (ntab+1 entries each), and ctab.  Copied to
  // every device of the process, once.  ierr: 0 ok, 3 error of CUDA, 4 no memory.
  void gpu_sph_tables_(int *ntab, double *gtab, double *dgdhtab, double *wtab, double *dwdhtab, double *dphidhtab, double *ctab, int *ierr)
  {
    *ierr = 0;
    try
    {
      if (!g_ready) throw SphError(3, "gpu_sph_tables_ called before gpu_sph_setup_");
      const int nt = *ntab + 1;
      std::vector<double2> gd(nt);
      std::vector<D4> wdp(nt);
      for (int k = 0; k < nt; k++)
      {
        gd[k] = make_double2(gtab[k], dgdhtab[k]);
        wdp[k] = make_d4(wtab[k], dwdhtab[k], dphidhtab[k], 0.0);
      }
      for (size_t k = 0; k < g_dev.size(); k++)
      {
        Dev &d = g_dev[k];
        CK(cudaSetDevice(d.id));
        dfree(d.gd); dfree(d.wdp);
        dmalloc(d.gd, nt); dmalloc(d.wdp, nt);
        CK(cudaMemcpy(d.gd, gd.data(), nt * sizeof(double2), cudaMemcpyHostToDevice));
        CK(cudaMemcpy(d.wdp, wdp.data(), nt * sizeof(D4), cudaMemcpyHostToDevice));
        d.ctab = *ctab;
        d.ntab = *ntab;
      }
      CK(cudaSetDevice(g_dev[0].id));
    }
    catch (const SphError &e)
    {
      fprintf(stderr, "StarSmasher (SPH on the GPU): %s\n", e.what());
      *ierr = e.code;
    }
  }

  // The smoothing lengths, densities and the other sums of the particles nlow..nup (from 1) and their neighbour lists, solved on the
  // GPU as rho_and_h does it on the CPU: hp (the guess on entry), rho, omega, 0mega, psi, wn, zeta, nn and first are set for every
  // particle that the GPU did (not wn of a point particle), the lists are written to list(1:ltot).  sflag(i-nlow+1) is 1 for a
  // particle that the GPU did not do, which is left as it was.  ierr: 0 ok, 2 the lists do not fit in list(1:listcap), 3 error of
  // CUDA, 4 no memory.  With ierr 2 the arrays hp... may have been changed (the caller keeps the guesses of h).
  void gpu_sph_rhoh_(int *n_, int *nlow_, int *nup_, int *nmaxarr, double *x, double *y, double *z, double *am, double *u, double *hp,
                     double *divv, double *rho, double *omega, double *om0, double *psi, double *wn, double *zeta, int *nn, int *first,
                     int *sflag, int *list, int *listcap_, int *ltot, double *dt, double *hfloor, int *nnopt, int *nselfg, int *ierr)
  {
    *ierr = 0;
    *ltot = 0;
    const int n = *n_, nlow = *nlow_ - 1, nup = *nup_, nown = nup - nlow;
    const size_t listcap = *listcap_;
    try
    {
      if (!g_ready) throw SphError(3, "gpu_sph_rhoh_ called before gpu_sph_setup_");
      const int ndev = (int)g_dev.size();
      CK(cudaSetDevice(g_dev[0].id));
      const size_t nm = (size_t)*nmaxarr * 8;
      std::vector<std::pair<const void *, size_t> > arr;
      const double *dp[13] = {x, y, z, am, u, hp, divv, rho, omega, om0, psi, wn, zeta};
      for (int k = 0; k < 13; k++) arr.push_back(std::make_pair((const void *)dp[k], nm));
      arr.push_back(std::make_pair((const void *)nn, (size_t)*nmaxarr * 4));
      arr.push_back(std::make_pair((const void *)first, (size_t)*nmaxarr * 4));
      arr.push_back(std::make_pair((const void *)list, listcap * 4));
      pin_all(arr);
      Phys ph;
      ph.dt = *dt; ph.hfloor = *hfloor; ph.nnopt = (double)*nnopt; ph.nselfg = *nselfg;
      const int plev = prof_level();
#ifdef SPH_CLOCKS
      if (plev >= 2)
      {
        unsigned long long z[16] = {0};
        CK(cudaSetDevice(g_dev[0].id));
        CK(cudaMemcpyToSymbol(g_clk, z, sizeof z));
      }
#endif
      const double t0 = now_ms();
      const int ch = std::min(chunk_size(), std::max(1, (nown + ndev - 1) / ndev));
      for (int k = 0; k < ndev; k++)
      {
        if (!g_dev[k].gd) throw SphError(3, "the tables of the kernel were not sent to the device");
        reserve(g_dev[k], n, nown, ch);
      }
      for (int k = 0; k < ndev; k++) upload_and_build(g_dev[k], n, nlow, nup, x, y, z, am, u, hp + nlow, divv + nlow);
      double t1 = 0;
      if (plev) { for (int k = 0; k < ndev; k++) { CK(cudaSetDevice(g_dev[k].id)); CK(cudaStreamSynchronize(g_dev[k].st)); } t1 = now_ms(); }
      std::vector<cudaEvent_t> ev;
      for (int k = 0; k < ndev; k++)
      {
        Dev &d = g_dev[k];
        CK(cudaSetDevice(d.id));
        const View V = makeview(d, ndev, k, n, nlow, nup, ph);
        const int ns = slice_of(nown, ndev, k);
        if (k == 0) mark("(before the kernels)", d.st);
        if (ns > 0) launch_chunks(d, V, ns, ch, (k == 0 && plev >= 2) ? &ev : 0);
      }
      Host H;
      H.hp = hp; H.rho = rho; H.omega = omega; H.om0 = om0; H.psi = psi; H.wn = wn; H.zeta = zeta; H.u = u;
      H.nn = nn; H.first = first; H.sflag = sflag; H.list = list; H.listcap = listcap;
      long nfl = 0;
      drain(nown, nlow, ch, H, ltot, ierr, &nfl);
      CK(cudaSetDevice(g_dev[0].id));
      if (plev)
      {
        double tg = 0, ts = 0;
        for (size_t e = 0; e + 1 < ev.size(); e += 2)
        {
          float b = 0;
          cudaEventElapsedTime(&b, ev[e], ev[e + 1]);
          ts += b;
        }
        for (size_t e = 0; e < ev.size(); e++) cudaEventDestroy(ev[e]);
        fprintf(stderr, "SPHGPU solve: upload+build %.2f ms, kernels+copies %.2f ms, all %.2f ms (%d particles of the rank, %d entries, %d devices, %ld left to the CPU, chunk %d)",
                t1 - t0, now_ms() - t1, now_ms() - t0, nown, *ltot, ndev, nfl, ch);
        if (plev >= 2) fprintf(stderr, "; kernel %.2f ms", ts + tg);
        fprintf(stderr, "\n");
#ifdef SPH_CLOCKS
        if (plev >= 2)
        {
          unsigned long long z[16];
          CK(cudaSetDevice(g_dev[0].id));
          CK(cudaMemcpyFromSymbol(z, g_clk, sizeof z));
          const double gr = z[8] > 0 ? (double)z[8] : 1.0;
          const double wn = (double)((nown + 31) / 32), wsn = (double)((nown + NTHR - 1) / NTHR) * (NTHR / 32);
          fprintf(stderr, "SPHGPU cycles per warp (k_solve): search and solve_h %.0f, all up to the lists %.0f, lists %.0f; searches again per particle %.4f, warps with one %.3f, evaluations per particle %.3f\n",
                  z[5] / wsn, z[6] / wsn, z[7] / wsn, (double)z[13] / nown, z[11] ? (double)z[12] / z[11] : 0.0, (double)z[14] / nown);
        }
#endif
        if (plev >= 2) print_marks("build phases (ms)");
      }
    }
    catch (const SphError &e)
    {
      fprintf(stderr, "StarSmasher (SPH on the GPU): %s\n", e.what());
      *ierr = e.code;
      if (!g_dev.empty()) cudaSetDevice(g_dev[0].id);
    }
  }
}
