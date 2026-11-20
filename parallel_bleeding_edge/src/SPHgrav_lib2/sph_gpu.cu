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
//   3. a thread takes a particle of the rank and walks the tree with a test of the boxes and of the points in single precision
//      with a margin (it can only let too many through, never drop one).  The places that pass are put in the order of the
//      index (a merge sort in the thread) and then tested, in that order, with the exact double precision arithmetic.
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
// The tests of non-negative squared distances are made on the integer bits of the doubles (dlt, dle): the same answers, and the double precision unit, which is slow on most GPUs, does less.
// The solve makes one evaluation of the function at a time (SolveSt); when only a few threads of a warp are not done after their first evaluations, the 32 threads of the warp make the
// evaluations of each of them together (coop_run), adding the terms of the sums in the order of the index.
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
constexpr int VMINLIFE = 12;  // calls a set of Verlet lists has to serve (the build included) for them to be worth building
#ifndef SPH_LEAF
#define SPH_LEAF 8
#endif
constexpr int LEAF = SPH_LEAF;   // particles of a leaf
constexpr int BR = 4;         // nodes of a level under a node
constexpr int MAXLEV = 24;
constexpr int NTHR = 128;     // threads of a block of the solve kernel
constexpr int NBLK = 4;
constexpr double GROW = 1.04;      // the first search reaches this factor past the radius asked for (the next ones too)
#ifndef SPH_RESCUEK
#define SPH_RESCUEK 3
#endif
#ifndef SPH_RESCUEMAX
#define SPH_RESCUEMAX 3
#endif
constexpr int RESCUEK = SPH_RESCUEK;      // evaluations of the solve a thread makes by itself before the warp looks at who is not done
constexpr int RESCUEMAX = SPH_RESCUEMAX;  // the most threads of a warp that the warp solves by working together
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
  // Verlet lists (bit 4 of ngpusph): vmode 0 off; 1 the candidates of the first search of a particle are read from the list of the call that built it when the
  // radius it covers, less twice the largest displacement since then, reaches the radius asked for; 2 this call builds the lists (a margin of vskin of the radius)
  int vmode;
  const double *px, *py, *pz, *pu;   // positions and u in the order of the index
  int *vl_len;                       // for each particle of the rank (index - nlow): its entries (-1: none)
  int *vl_idx;                       // CAP indices of the particle of each (in the order of the index)
  double *vl_R;                      // the radius they cover, from the position at the build
  const double *dmaxp;               // the largest displacement of a particle since the build
  double vskin;
  int *nfall;                        // particles with a list that did not reach: their first search went to the tree
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

// ---------------------------------------------------------------------------------------------------- the displacement since the Verlet lists were built
// (the largest of them: an exact maximum, whatever the order in which the threads meet)
__global__ void k_vdmax(int n, const double *x, const double *y, const double *z, const double *xr, const double *yr, const double *zr, double *part)
{
  __shared__ double s[256];
  const int t = threadIdx.x;
  double m = 0.0;
  for (int i = blockIdx.x * blockDim.x + t; i < n; i += gridDim.x * blockDim.x)
  {
    const double dx = x[i] - xr[i], dy = y[i] - yr[i], dz = z[i] - zr[i];
    const double d2 = (dx * dx + dy * dy) + dz * dz;
    m = d2 > m ? d2 : m;
  }
  s[t] = m;
  __syncthreads();
  for (int h = 128; h > 0; h >>= 1)
  {
    if (t < h) s[t] = s[t + h] > s[t] ? s[t + h] : s[t];
    __syncthreads();
  }
  if (t == 0) part[blockIdx.x] = s[0];
}

__global__ void k_vdmax_final(int nb, const double *part, double *dmax)
{
  __shared__ double s[256];
  const int t = threadIdx.x;
  double m = 0.0;
  for (int b = t; b < nb; b += blockDim.x) m = part[b] > m ? part[b] : m;
  s[t] = m;
  __syncthreads();
  for (int h = 128; h > 0; h >>= 1)
  {
    if (t < h) s[t] = s[t + h] > s[t] ? s[t + h] : s[t];
    __syncthreads();
  }
  if (t == 0) dmax[0] = sqrt(s[0]);
}

// ---------------------------------------------------------------------------------------------------- the search of one particle
// The candidates of the particle at (xi,yi,zi) found from the root of the tree: every particle at a squared distance not above Rc2,
// with the squared distance of the CPU.  They go to ck (the index of the particle times 2, plus 1 for a gas particle) and cd, element
// k at [k*st], in the order of the index.  Returns their number, or -1 if there are more than CAP.
//
// The visit of the tree makes the tests in single precision only and queues one 64-bit entry for each place that passes, (key << 32) | place in the
// tree, where key is what goes to ck (2*index + gas).  The keys are all different, so any sort of the entries as numbers gives the order of the
// index.  The entries are sorted, and then one loop in that order makes the test in double precision, from the place of the entry, and keeps the
// candidates where they are: the same set, in the same order, with the same squared distance as the test inside the visit gave.
//
// Why: the lanes of a warp are in different leaves, so the test in double precision inside the visit was run by the warp for almost every place
// of a leaf (about 5 times the number of candidates), while in a loop of its own all the lanes run it together, and the double precision of these
// GPUs is 1/64 of the single.  The sort by insertion of two vectors (the keys and the squared distances) had a loop whose length depends on the
// data, and a warp waits for its lane that moves most.  Measured on an RTX PRO 6000 Blackwell (one device, ngpusph = 7, time of k_solve on the wall
// for a call without Verlet lists): 8.9 -> 5.9 ms at 500k particles, 64.7 -> 56.8 ms at 2M, and in the collision of two stars 5.7 -> 4.5 (start),
// 6.7 -> 5.3 (t = 8), 7.2 -> 6.3 (periastron), 8.3 -> 7.2 ms (late); the call that builds the Verlet lists, where the sort is half of the time,
// 17.9 -> 9.7 ms at 500k.  The dense test alone, with an insertion sort of the vector of entries, gives 6.5 ms at 500k; Batcher's sorting network
// instead of the merge sort 6.2 ms; leaves of 16 places and 8 nodes under a node did not help.  Tried before and not kept: the rank of each key by
// comparison with all the others, and a walk of the tree by the 32 threads of a warp together for a group of particles, the places sorted once.
typedef unsigned long long u64;
static_assert(CAP <= 128, "the merge sort of collect has four passes (runs of 8 to 128)");

// the entries are held in cd, as the bits of a double (no floating-point operation touches them)
__device__ __forceinline__ u64 ent_ld(const double *v, int a) { return (u64)__double_as_longlong(v[a]); }
__device__ __forceinline__ void ent_st(double *v, int a, u64 e) { v[a] = __longlong_as_double((long long)e); }

// compare-exchange of two entries (a minimum and a maximum, no branch)
__device__ __forceinline__ void cex(u64 &a, u64 &b)
{
  const u64 lo = a < b ? a : b, hi = a < b ? b : a;
  a = lo;
  b = hi;
}

// the sort of 8 entries by a network of 19 compare-exchanges
__device__ __forceinline__ void sort8(u64 *v)
{
  cex(v[0], v[2]); cex(v[1], v[3]); cex(v[4], v[6]); cex(v[5], v[7]);
  cex(v[0], v[4]); cex(v[1], v[5]); cex(v[2], v[6]); cex(v[3], v[7]);
  cex(v[0], v[1]); cex(v[2], v[3]); cex(v[4], v[5]); cex(v[6], v[7]);
  cex(v[2], v[4]); cex(v[3], v[5]);
  cex(v[1], v[4]); cex(v[3], v[6]);
  cex(v[1], v[2]); cex(v[3], v[4]); cex(v[5], v[6]);
}

// one pass of the merge sort: the runs of w entries of src (element k at [k*sst]), two by two, merged into dst (element k at [k*dst_st]).  The loop is over the
// index of the output and the choice of the run it is taken from is made without a branch (an entry past the end of a run reads as the largest value), so the
// lanes of a warp go together
__device__ __forceinline__ void merge_pass(const double *src, int sst, double *dst, int dst_st, int n, int w)
{
  for (int lo = 0; lo < n; lo += 2 * w)
  {
    const int mid = min(lo + w, n), hi = min(lo + 2 * w, n);
    int i = lo, j = mid;
    u64 a = ent_ld(src, i * sst);
    u64 b = j < hi ? ent_ld(src, j * sst) : ~0ull;
    for (int k = lo; k < hi; k++)
    {
      const bool ta = a <= b;
      ent_st(dst, k * dst_st, ta ? a : b);
      i += ta ? 1 : 0;
      j += ta ? 0 : 1;
      const int idx = ta ? i : j, lim = ta ? mid : hi;
      const u64 v = ent_ld(src, min(idx, n - 1) * sst);
      const u64 nv = idx < lim ? v : ~0ull;
      a = ta ? nv : a;
      b = ta ? b : nv;
    }
  }
}

// The np entries in cd (element k at [k*st]) sorted: runs of 8 (the last one filled with the largest value) by the network, from cd to tmp, then merge passes that
// double the runs, from tmp to cd and back.  tmp is a vector of its own, element k at [k]: the stride st is that of cd and ck, the vectors of the caller.
// Returns the number of passes: the sorted entries are in cd if it is odd, in tmp if it is even.
__device__ __forceinline__ int sort_entries(double *cd, double *tmp, int np, int st)
{
  const int n = (np + 7) & ~7;
  for (int r = 0; r < n; r += 8)
  {
    u64 v[8];
#pragma unroll
    for (int q = 0; q < 8; q++) v[q] = (r + q < np) ? ent_ld(cd, (r + q) * st) : ~0ull;
    sort8(v);
#pragma unroll
    for (int q = 0; q < 8; q++) ent_st(tmp, r + q, v[q]);
  }
  int passes = 0;
  if (8 < n) { merge_pass(tmp, 1, cd, st, n, 8); passes = 1; }
  if (16 < n) { merge_pass(cd, st, tmp, 1, n, 16); passes = 2; }
  if (32 < n) { merge_pass(tmp, 1, cd, st, n, 32); passes = 3; }
  if (64 < n) { merge_pass(cd, st, tmp, 1, n, 64); passes = 4; }
  return passes;
}

// a < b and a <= b for two doubles that are not negative and not NaN, on the integer bits (the order of the bits is the order of the numbers): the same answer as the double
// comparison, made by the integer unit (the double precision unit of most GPUs does a fraction of an operation of the others a cycle)
__device__ __forceinline__ bool dlt(double a, double b) { return __double_as_longlong(a) < __double_as_longlong(b); }
__device__ __forceinline__ bool dle(double a, double b) { return __double_as_longlong(a) <= __double_as_longlong(b); }

// the test in double precision of the entries in src (element k at [k*sst]), in order: the place read from each, the squared distance of the CPU, the candidates
// kept in the order they are in, in ck and cd (element k at [k*st]).  The index written never passes the index read, so src may be cd itself (with sst = st).
// Returns the number of candidates.
__device__ __forceinline__ int dense_exact(const View &V, double xi, double yi, double zi, double Rc2, const double *src, int sst, unsigned *ck, double *cd, int np, int st)
{
  int m = 0;
  for (int a = 0; a < np; a++)
  {
    const u64 e = ent_ld(src, a * sst);
    const D4 r = V.rec[(unsigned)(e & 0xffffffffull)];
    const double ex = xi - r.x, ey = yi - r.y, ez = zi - r.z;
    const double d2 = (ex * ex + ey * ey) + ez * ez;
    if (dle(d2, Rc2))
    {
      ck[m * st] = (unsigned)(e >> 32);
      cd[m * st] = d2;
      m++;
    }
  }
  return m;
}

// The search as it was before: the test in double precision inside the visit and an insertion sort of the candidates.  It is the way of a search with more
// than CAP places through the single precision test, which finds out whether the candidates are more than CAP (a place that passes is not always a
// candidate) and gives them when they are not.
__device__ int collect_by_insertion(const View &V, double xi, double yi, double zi, double Rc2, unsigned *ck, double *cd, int st)
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
  int sp = 0, np = 0;
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
          if (np == CAP) return collect_by_insertion(V, xi, yi, zi, Rc2, ck, cd, st);
          const unsigned key = ((unsigned)V.perm[p] << 1) | (f.w > 0.5f ? 1u : 0u);
          ent_st(cd, np * st, ((u64)key << 32) | (unsigned)p);
          np++;
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
  double tmp[CAP];
  const int passes = sort_entries(cd, tmp, np, st);
  // (the exact test is written twice so that each call reads a vector the compiler knows)
  if (passes & 1) return dense_exact(V, xi, yi, zi, Rc2, cd, st, ck, cd, np, st);
  return dense_exact(V, xi, yi, zi, Rc2, tmp, 1, ck, cd, np, st);
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
  int vslot;       // Verlet lists: the place of the particle in them (-1: not in use)
  bool vhas;       // it has a list
  double vRc2;     // the squared radius to which that list is complete now (negative: none)
};

// The candidates of the particle from its Verlet list: every entry at a squared distance up to Rc2 (the list is complete out to that radius), in the order of the
// index, with the squared distance of the CPU: the same candidates, in the same order, as collect() finds from the root
__device__ int collect_vl(const View &V, int slot, double xi, double yi, double zi, double Rc2, unsigned *ck, double *cd, int st)
{
  const int len = V.vl_len[slot];
  const int *li = V.vl_idx + (size_t)slot * CAP;
  int m = 0;
  for (int a = 0; a < len; a++)
  {
    const int j = li[a];
    const double ex = xi - V.px[j], ey = yi - V.py[j], ez = zi - V.pz[j];
    const double d2 = (ex * ex + ey * ey) + ez * ez;
    if (d2 <= Rc2)
    {
      if (m == CAP) return -1;
      ck[m * st] = ((unsigned)j << 1) | (V.pu[j] != 0.0 ? 1u : 0u);
      cd[m * st] = d2;
      m++;
    }
  }
  return m;
}

// the candidates again, out to a squared radius r2 (with the margin): from the Verlet list if it reaches that far, else from the root of the tree; false if there are too many
__device__ bool recollect(const View &V, double xi, double yi, double zi, double r2, Cand &C)
{
  const double Rn = r2 * (GROW * GROW);
  int m;
  if (C.vslot >= 0 && Rn <= C.vRc2)
    m = collect_vl(V, C.vslot, xi, yi, zi, Rn, C.ck, C.cd, C.st);
  else
  {
    if (C.vslot >= 0 && C.vhas && C.nrec == 0) atomicAdd(V.nfall, 1);   // its list does not reach: the first search goes to the tree
    m = collect(V, xi, yi, zi, Rn, C.ck, C.cd, C.st);
  }
  C.nrec++;
  if (m < 0) return false;
  C.m = m;
  C.Rc2 = Rn;
  return true;
}

// the first search of a particle: in the call that builds the Verlet lists it reaches vskin further and its candidates are kept (unless they are too many)
__device__ bool first_search(const View &V, int slot, double xi, double yi, double zi, double r2, Cand &C)
{
  if (V.vmode == 2)
  {
    const double f = 1.0 + V.vskin;
    const double Rb = r2 * (GROW * GROW) * (f * f);
    const int m = collect(V, xi, yi, zi, Rb, C.ck, C.cd, C.st);
    if (m >= 0)
    {
      int *li = V.vl_idx + (size_t)slot * CAP;
      for (int a = 0; a < m; a++) li[a] = (int)(C.ck[a * C.st] >> 1);
      V.vl_len[slot] = m;
      V.vl_R[slot] = sqrt(Rb);
      C.nrec++;
      C.m = m;
      C.Rc2 = Rb;
      return true;
    }
    V.vl_len[slot] = -1;
  }
  return recollect(V, xi, yi, zi, r2, C);
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
    if (dlt(d2, r2) && (C.ck[k * C.st] & 1u))
    {
      double xq = ctoh2 * d2;
      xq = dlt(xq, V.c4) ? xq : V.c4;
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
// by drifting (label 10), then Newton steps kept inside the bracket by bisection (label 12).  It is written as a machine that makes one
// evaluation of the function at a time (SolveSt: where it stands between two evaluations), with the same operations in the same order
// as the loops of solve_particle, so that it can be stopped after any evaluation and go on elsewhere with the same results (the warp that makes it together, below).
// bfun is called from one place only, so the threads of a warp that are at different stages of their solves evaluate together: with a call in each of the three
// loops a warp paid the first evaluation, the most trips of the loop of label 10 and the most trips of the loop of label 12 one after the other.
struct SolveSt
{
  double xl, xh, fl, dx, ezr, xacc, dxmax, h, hsol;
  int neval, j;
  int phase;     // 0: the next evaluation is the first one; 1: the loop of label 10; 2: the loop of label 12; 3: done (hsol is the smoothing length without the floor); 4: failed
};

__device__ __forceinline__ void solve_init(const View &V, double hin, double divvi, SolveSt &S)
{
  const double h0 = hin - V.hfloor;
  S.xacc = 2.0e-6 * h0;
  const double a1 = h0 * (1.0 + (divvi * V.dt) / 3.0);
  const double hpguess = a1 > S.xacc ? a1 : S.xacc;
  const double b1 = V.dxfrac * h0, b2 = fabs((h0 * divvi) * V.dt);
  const double b3 = b1 > b2 ? b1 : b2, b4 = 0.5 * hpguess;
  const double dxm = b3 < b4 ? b3 : b4;
  S.xl = hpguess;
  S.dxmax = fabs(dxm);
  S.h = S.xl;
  S.xh = 0.0; S.fl = 0.0; S.dx = 0.0; S.ezr = 0.0; S.hsol = 0.0;
  S.neval = 0; S.j = 0; S.phase = 0;
}

// the start of an iteration of label 12 (the bisection when the Newton step leaves the bracket, the end of the solve when the step is small)
__device__ __forceinline__ void solve_l12(SolveSt &S)
{
  if ((S.ezr - S.xh) * (S.ezr - S.xl) >= 0.0)
  {
    S.dx = 0.5 * (S.xh - S.xl);
    S.ezr = S.xl + S.dx;
  }
  if (fabs(S.dx) < S.xacc) { S.hsol = S.ezr; S.phase = 3; return; }
  S.h = S.ezr;
  S.phase = 2;
}

// what is done with the value f of the function and its derivative df at S.h; sets S.h to the next point, or the phase to 3 or 4
__device__ __forceinline__ void solve_post(SolveSt &S, double f, double df)
{
  if (S.phase == 0)
  {
    S.neval++;
    S.fl = f;
    S.dx = f / df;
    if (fabs(S.dx) < S.xacc) { S.hsol = S.xl - S.dx; S.phase = 3; return; }
    if (S.fl == 0.0) { S.hsol = S.xl; S.phase = 3; return; }
    {
      const double c1 = fabs(S.dx) < S.dxmax ? fabs(S.dx) : S.dxmax;
      const double c2 = fabs(0.5 * S.xl);
      S.dx = copysign(c1 < c2 ? c1 : c2, S.dx);
    }
    // label 10: the root is not yet bracketed, drift a little more
    S.xh = S.xl - S.dx;
    S.h = S.xh;
    S.phase = 1;
    return;
  }
  if (S.phase == 1)
  {
    if (++S.neval > MAXEVAL) { S.phase = 4; return; }
    S.dx = f / df;
    if (fabs(S.dx) < S.xacc) { S.hsol = S.xh - S.dx; S.phase = 3; return; }
    if (f == 0.0) { S.hsol = S.xh; S.phase = 3; return; }
    S.dx = copysign(fabs(S.dx) < S.dxmax ? fabs(S.dx) : S.dxmax, S.dx);
    S.ezr = S.xh - S.dx;
    if (S.fl * f > 0.0)
    {
      S.dx = copysign(S.dx, f);
      S.xl = S.xh;
      S.xh = S.xl - S.dx;
      S.h = S.xh;
      return;
    }
    // the root is bracketed: xl must give a negative value
    if (S.fl > 0.0) { const double xtmp = S.xh; S.xh = S.xl; S.xl = xtmp; }
    S.j = 1;
    solve_l12(S);
    return;
  }
  // label 12
  if (++S.neval > MAXEVAL) { S.phase = 4; return; }
  if (f < 0.0) S.xl = S.ezr; else S.xh = S.ezr;
  S.dx = f / df;
  S.ezr = S.ezr - S.dx;
  if (++S.j > 50)
  {
    // the maximum number of iterations: the CPU code writes a warning and starts again from the middle of the bracket
    S.dx = 0.5 * (S.xh - S.xl);
    S.ezr = S.xl + S.dx;
    S.j = 1;
  }
  solve_l12(S);
}

// one evaluation and what follows it; false if the candidates could not be made
__device__ __forceinline__ bool solve_step(const View &V, double xi, double yi, double zi, Cand &C, SolveSt &S)
{
  double f, df;
  const double r2 = 4.0 * (S.h * S.h);
  if (!bfun(V, xi, yi, zi, C, r2, S.h, f, df)) { S.phase = 4; return false; }
  solve_post(S, f, df);
  return true;
}

// the candidates of the Verlet lists of the thread (the part of the start of k_solve that does not search)
__device__ __forceinline__ void vl_setup(const View &V, int slot, Cand &C)
{
  if (V.vmode == 1)
  {
    C.vslot = slot;
    C.vhas = V.vl_len[slot] >= 0;
    if (C.vhas)
    {
      const double Rc = V.vl_R[slot] - 2.0 * (*V.dmaxp);
      if (Rc > 0.0)
      {
        const double Rs = Rc * (1.0 - 1.0e-9);
        C.vRc2 = Rs * Rs;
      }
    }
  }
}

// what a thread makes of a particle once its smoothing length is solved
struct Fin
{
  double hfin, r2f, rho, omg, om0, psi, wn;
  int nn, flag;
};

// The smoothing length and the sums over the candidates (bonetsumwithphi) of a particle whose solve ended (ok: the search and the solve went well; hs the solution
// without the floor).  A point particle (u=0, dynhco=0) keeps its smoothing length (hin) and gets the density.
__device__ __forceinline__ void finish_particle(const View &V, Cand &C, bool ok, bool gas, double hin, double hs, const D4 &r, Fin &F)
{
  unsigned *ck = C.ck;
  double *cd = C.cd;
  F.nn = 0; F.flag = 0; F.rho = 0.0; F.omg = 0.0; F.om0 = 0.0; F.psi = 0.0; F.wn = 0.0;
  if (gas)
  {
    if (ok)
    {
      F.hfin = hs + V.hfloor;
      F.r2f = 4.0 * (F.hfin * F.hfin);
      if (F.r2f > C.Rc2) ok = recollect(V, r.x, r.y, r.z, F.r2f, C);
    }
  }
  else
  {
    F.hfin = hin;
    F.r2f = 4.0 * (F.hfin * F.hfin);
  }
  if (!ok)
  {
    F.flag = 1;
    F.hfin = hin;
    return;
  }
  const double r2f = F.r2f;
  int nn = 0;
  for (int a = 0; a < C.m; a++)
    if (dlt(cd[(a)], r2f)) nn++;
  F.nn = nn;
  const double hpi = F.hfin;
  const double h2 = hpi * hpi;
  const double ctoh2 = V.ctab / h2;
  double rho = 0.0, omg = 0.0, om0 = 0.0, psi = 0.0, wn = 0.0;
  if (!gas)
  {
    // bonetsumwithphi without the sums of the smoothing length (point particle, dynhco=0)
    for (int a = 0; a < C.m; a++)
    {
      const double r2p = cd[(a)];
      if (!dlt(r2p, r2f)) continue;
      if (ck[(a)] & 1u)
      {
        double xq = ctoh2 * r2p;
        xq = dlt(xq, V.c4) ? xq : V.c4;
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
      if (!dlt(r2p, r2f)) continue;
      const unsigned jj = ck[(a)] >> 1;
      const bool gj = (ck[(a)] & 1u) != 0u;
      const double amj = V.am[jj];
      const double xqr = ctoh2 * r2p;
      const double xq = dlt(xqr, V.c4) ? xqr : V.c4;
      const int it = (int)xq;
      const double fq = xq - (double)it;
      const D4 w0 = V.wdp[it], w1 = V.wdp[it + 1];
      if (dlt(r2p, 4.0 * h2))
      {
        if ((V.nselfg == 1 || !gj) && true)
          psi = psi + amj * (w0.z + fq * (w1.z - w0.z));
      }
      if (gj)
      {
        if (dlt(r2p, 4.0 * h2t))
        {
          // without a floor h2t is h2, the index and the fraction of the smoothing length without the floor are the same ones, unless the
          // table was clamped
          int itt = it;
          double fqt = fq;
          if (!(samet && dle(xqr, V.c4)))
          {
            const double xqt = ctoh2t * r2p;
            itt = (int)xqt;
            fqt = xqt - (double)itt;
          }
          const double2 g0 = V.gd[itt], g1 = V.gd[itt + 1];
          wn = wn + (g0.x + fqt * (g1.x - g0.x));
          om0 = om0 - (g0.y + fqt * (g1.y - g0.y));
        }
        if (dlt(r2p, 4.0 * h2))
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
  F.rho = rho; F.omg = omg; F.om0 = om0; F.psi = psi; F.wn = wn;
}

// The places of the lists are found by the warp, not by the block: a warp then waits only for its own threads (the iterations of the solve differ a lot from one
// particle to the next).  Every thread of the warp calls it (valid: it has a particle to give).  The list of the particle is written at its place, and its record.
typedef cub::WarpScan<int> WarpScanI;
__device__ __forceinline__ void emit_particle(const View &V, typename WarpScanI::TempStorage &wsm, bool valid, const Fin &F, const Cand &C, int i, int outidx)
{
  const unsigned *ck = C.ck;
  const double *cd = C.cd;
  int off, total;
  const int mine = (valid && !F.flag) ? F.nn : 0;
  WarpScanI(wsm).ExclusiveSum(mine, off, total);
  int wbase = 0;
  if ((threadIdx.x & 31) == 0 && total > 0) wbase = atomicAdd(V.ltot, total);
  wbase = __shfl_sync(0xffffffffu, wbase, 0);
  if (valid)
  {
    const int base = wbase + off;
    if (!F.flag)
    {
      int c = 0;
      const double r2f = F.r2f;
      for (int a = 0; a < C.m; a++)
        if (dlt(cd[(a)], r2f)) V.list[base + c++] = (int)(ck[(a)] >> 1) + 1;   // (the lists of the CPU hold the index from 1, as Fortran numbers the particles)
    }
    SphOut o;
    o.hp = F.hfin; o.rho = F.rho; o.omega = F.omg; o.om0 = F.om0; o.psi = F.psi; o.wn = F.wn;
    o.idx = i; o.nn = F.flag ? 0 : F.nn; o.first = base; o.flag = F.flag;
    V.out[outidx] = o;
  }
}

// The evaluations of the function of the solve of one thread (lane s) made by the 32 lanes of the warp together.  Lane j has the candidates j, j+32, ... in cka and cda (their places and
// squared distances).  The terms of the sums are made by the lanes at once and then added in the order of the index through the shuffles (every lane adds the terms of the candidates 0, 1, 2, ...
// in turn), so the sums are added term by term in the order bfun adds them, and the terms are the same numbers: the same bits.  A thread alone makes an evaluation with some hundreds of
// double precision instructions (on the GPUs of the RTX line and on the L40S a warp does one of them in 64 cycles of its scheduler), the warp together with some tens, and the slowest solves
// of a call take hundreds of evaluations one after the other: they are the longest chain of the kernel.  Every lane of the warp calls it with the same T (the state of the solve, copied from lane s)
// and it ends with the solve done.  When a trial radius goes past the candidates they are searched again from the root by lane s alone (in its own arrays and Cand C), and the lanes take them from there.
// False if that failed.
__device__ __forceinline__ bool coop_run(const View &V, int lane, int s, const D4 &r, Cand &C, SolveSt &T, unsigned (&cka)[CAP / 32], double (&cda)[CAP / 32], int &m, double &Rc2)
{
  const unsigned full = 0xffffffffu;
  while (T.phase < 3)
  {
    const double h = T.h;
    const double r2 = 4.0 * (h * h);
    if (r2 > Rc2)
    {
      bool okc = true;
      if (lane == s) okc = recollect(V, r.x, r.y, r.z, r2, C);
      okc = __shfl_sync(full, (int)okc, s) != 0;
      if (!okc) { T.phase = 4; return false; }
      m = __shfl_sync(full, C.m, s);
      Rc2 = __shfl_sync(full, C.Rc2, s);
#pragma unroll
      for (int q = 0; q < CAP / 32; q++) { cka[q] = 0u; cda[q] = 0.0; }
      for (int k = 0; k < m; k++)
      {
        const unsigned w = __shfl_sync(full, C.ck[k], s);
        const double v = __shfl_sync(full, C.cd[k], s);
        if ((k & 31) == lane)
        {
#pragma unroll
          for (int q = 0; q < CAP / 32; q++) if ((k >> 5) == q) { cka[q] = w; cda[q] = v; }
        }
      }
    }
    const double h2 = h * h;
    const double ctoh2 = V.ctab / h2;
    double tw[CAP / 32], to[CAP / 32];
    unsigned pm[CAP / 32];
#pragma unroll
    for (int q = 0; q < CAP / 32; q++)
    {
      const int k = lane + 32 * q;
      const double d2 = cda[q];
      const bool pass = (k < m) && dlt(d2, r2) && (cka[q] & 1u);
      double a = 0.0, b = 0.0;
      if (pass)
      {
        double xq = ctoh2 * d2;
        xq = dlt(xq, V.c4) ? xq : V.c4;
        const int it = (int)xq;
        const double fq = xq - (double)it;
        const double2 g0 = V.gd[it], g1 = V.gd[it + 1];
        a = g0.x + fq * (g1.x - g0.x);
        b = g0.y + fq * (g1.y - g0.y);
      }
      tw[q] = a; to[q] = b;
      pm[q] = __ballot_sync(full, pass);
    }
    double wn = 0.0, om0 = 0.0;
#pragma unroll
    for (int q = 0; q < CAP / 32; q++)
    {
      if (32 * q < m)
      {
        const int nk = min(32, m - 32 * q);
#pragma unroll 1   // (not unrolled: with the 32 shuffles of doubles written out four times ptxas of CUDA 13.1 refuses the kernel for sm_120)
        for (int kk = 0; kk < nk; kk++)
        {
          const double ta = __shfl_sync(full, tw[q], kk), tb = __shfl_sync(full, to[q], kk);
          if ((pm[q] >> kk) & 1u) { wn = wn + ta; om0 = om0 - tb; }
        }
      }
    }
    om0 = om0 / h;
    solve_post(T, wn - V.nnopt, om0);
  }
  return T.phase == 3;
}

// The threads of a warp that are not done after their first RESCUEK evaluations: when there are at most RESCUEMAX of them the warp solves them one after the other, all its lanes working on one
// (coop_run); with more, each goes on by itself, one evaluation at a time, until there are few.  Every lane of the warp calls it (want: this lane has a solve that is not done).
// Measured on dumps of the code (the kernel alone, bits equal to the version before): with the rescue the time of k_solve falls by 2 to 52% on an L40S and by 3 to 48% on an H100.
// Tried and not kept: the threads not done after RESCUEK evaluations handed over to a second kernel after k_solve, with their state and candidates in a buffer: one thread each was
// slower than without it at the collision (the longest solves then run after the main kernel instead of during it), one warp each took 47% off the kernel instead of 58%.
__device__ __forceinline__ void rescue_solve(const View &V, int lane, bool want, const D4 &r, Cand &C, SolveSt &S, bool &ok)
{
  const unsigned full = 0xffffffffu;
  for (;;)
  {
    const unsigned um = __ballot_sync(full, want && S.phase < 3);
    if (!um) break;
    if (__popc(um) > RESCUEMAX)
    {
      if (want && S.phase < 3) { solve_step(V, r.x, r.y, r.z, C, S); if (S.phase == 4) ok = false; }
      continue;
    }
    const int s = __ffs(um) - 1;
    SolveSt T;
    T.xl = __shfl_sync(full, S.xl, s); T.xh = __shfl_sync(full, S.xh, s); T.fl = __shfl_sync(full, S.fl, s); T.dx = __shfl_sync(full, S.dx, s);
    T.ezr = __shfl_sync(full, S.ezr, s); T.xacc = __shfl_sync(full, S.xacc, s); T.dxmax = __shfl_sync(full, S.dxmax, s); T.h = __shfl_sync(full, S.h, s);
    T.hsol = 0.0;
    T.neval = __shfl_sync(full, S.neval, s); T.j = __shfl_sync(full, S.j, s); T.phase = __shfl_sync(full, S.phase, s);
    int m = __shfl_sync(full, C.m, s);
    double Rc2 = __shfl_sync(full, C.Rc2, s);
    unsigned cka[CAP / 32];
    double cda[CAP / 32];
#pragma unroll
    for (int q = 0; q < CAP / 32; q++) { cka[q] = 0u; cda[q] = 0.0; }
    for (int k = 0; k < m; k++)
    {
      const unsigned w = __shfl_sync(full, C.ck[k], s);
      const double v = __shfl_sync(full, C.cd[k], s);
      if ((k & 31) == lane)
      {
#pragma unroll
        for (int q = 0; q < CAP / 32; q++) if ((k >> 5) == q) { cka[q] = w; cda[q] = v; }
      }
    }
    const bool good = coop_run(V, lane, s, r, C, T, cka, cda, m, Rc2);
    if (lane == s)
    {
      S = T;
      C.m = m; C.Rc2 = Rc2;
      if (!good) ok = false;
    }
  }
}

// The smoothing length of every particle of the rank, its density and the other sums over the neighbours (bonetsumwithphi), and its
// neighbour list, as solve_particle makes them.  A point particle (u=0, dynhco=0) keeps its smoothing
// length and gets the list and the density.
__global__ void __launch_bounds__(NTHR, NBLK) k_solve(View V)
{
  __shared__ typename WarpScanI::TempStorage ws[NTHR / 32];
  const int nt = *V.ntgt;
  const int t0 = (int)(((long long)nt * V.dev) / V.ndev);
  const int tl = blockIdx.x * NTHR + threadIdx.x;      // the particle in the chunk
  const int tid = t0 + V.tbase + tl;
  const bool valid = tl < V.cn;
  const int lane = threadIdx.x & 31;
  int i = -1;
  // the candidates of the thread (local memory: the compiler lays it out so that the threads of a warp read neighbouring words)
  unsigned ck[CAP];
  double cd[CAP];
  Cand C;
  C.ck = ck; C.cd = cd; C.st = 1; C.m = 0; C.Rc2 = -1.0; C.nrec = 0; C.nev = 0; C.vslot = -1; C.vhas = false; C.vRc2 = -1.0;
  SolveSt S;
  S.phase = 3; S.hsol = 0.0;
  D4 r = make_d4(0.0, 0.0, 0.0, 0.0);
  bool gas = false, ok = true;
  double hin = 0.0;
  if (valid)
  {
    const int k = V.tgt[tid];
    i = V.perm[k];
    r = V.rec[k];
    gas = V.frec[k].w > 0.5f;
    hin = V.hp[i - V.nlow];
    const int slot = i - V.nlow;
    vl_setup(V, slot, C);
    ok = !(V.flagmod > 0 && (i + 1) % V.flagmod == 0);
    if (!gas)
    {
      // the list at the smoothing length it has
      if (ok) ok = first_search(V, slot, r.x, r.y, r.z, 4.0 * (hin * hin), C);
    }
    else
    {
      // the first search reaches past the guess by the margin, the later ones past the radius asked for
      const double divvi = V.divv[i - V.nlow];
      const double hg = hin - V.hfloor;
      const double hguess0 = hg * (1.0 + (divvi * V.dt) / 3.0);
      const double hstart = hguess0 > 2.0e-6 * hg ? hguess0 : 2.0e-6 * hg;
      if (ok) ok = first_search(V, slot, r.x, r.y, r.z, 4.0 * (hstart * hstart), C);
      if (ok)
      {
        solve_init(V, hin, divvi, S);
        for (int t = 0; t < RESCUEK && S.phase < 3; t++) solve_step(V, r.x, r.y, r.z, C, S);
      }
    }
  }
  // the threads that are not done: when there are few of them in the warp the whole warp solves them, one after the other (every thread of the warp is here)
  rescue_solve(V, lane, valid && gas && ok, r, C, S, ok);
  if (valid && gas && ok) ok = S.phase == 3;
  Fin F;
  F.hfin = hin; F.r2f = 0.0; F.rho = 0.0; F.omg = 0.0; F.om0 = 0.0; F.psi = 0.0; F.wn = 0.0; F.nn = 0; F.flag = 0;
  if (valid) finish_particle(V, C, ok, gas, hin, S.hsol, r, F);
  emit_particle(V, ws[threadIdx.x >> 5], valid, F, C, i, V.tbase + tl);
}


// ---------------------------------------------------------------------------------------------------- one device
// What the CPU needs of the results of a call of the kernel, for the particles in the order of the particles (not in the order of the curve): hp, rho, nn and whether the GPU left
// the particle to the CPU.  The rest of the records (omega, 0mega, psi, wn, the places of the lists) is read only by the pair sums of the CPU, which are not made in a run that keeps
// the lists on the device; fetch_records brings it when the CPU needs it.
struct HRec
{
  double hp, rho;
  int nn, flag;
};
constexpr int FLCAP = 4096;   // the records left to the CPU that are listed on the device (the rest are found by a scan of the records)

// the records of the call in the order of the particles; the records left to the CPU are listed (their place in out and their particle)
__global__ void dev_unperm(int nown, int nlow, const SphOut *out, HRec *hrec, int *fcnt, int2 *flo)
{
  const int o = blockIdx.x * blockDim.x + threadIdx.x;
  if (o >= nown) return;
  const SphOut s = out[o];
  HRec r;
  r.hp = s.hp; r.rho = s.rho; r.nn = s.nn; r.flag = s.flag;
  hrec[s.idx - nlow] = r;
  if (s.flag)
  {
    const int q = atomicAdd(fcnt, 1);
    if (q < FLCAP) flo[q] = make_int2(o, s.idx);
  }
}

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
  HRec *hrecd = 0, *hrech = 0; // the results the CPU needs, in the order of the particles: on the device, and pinned
  int *fcntd = 0, *fcnth = 0;  // how many records the kernel left to the CPU: on the device, and pinned
  int2 *flod = 0, *floh = 0;   // the first FLCAP of them (place in out, particle): on the device, and pinned
  // Verlet lists (see View): the lists of the particles of the rank, where the positions were when they were built, the displacement since
  int *vl_len = 0, *vl_idx = 0, *vnfall = 0;
  double *vl_R = 0, *vxr = 0, *vyr = 0, *vzr = 0, *vdmax = 0, *vpart = 0;
  size_t vnown = 0, vn = 0;      // room for this many particles of the rank and in all
  bool vvalid = false;           // lists have been built, for the particles vbn of which vbnown are the rank's from vbnlow
  int vbn = 0, vbnlow = 0, vbnown = 0;
  bool vforce = false;           // the next call builds them again (too many lists did not reach)
  int vused = 0;                 // calls that have used the lists since they were built
  int voff = 0;                  // calls left in which the lists are not used (they were built again too soon)
  int vofflen = 24;              // how many the next time
  std::vector<cudaEvent_t> evS, evC;   // the kernels of chunk c done / the copies of chunk c done
};

std::vector<Dev> g_dev;
bool g_ready = false;

// what the last call of the smoothing-length kernel left on the devices (the forces can use it if it is for the same particles): valid after a call
// that made every particle on the GPU
struct BState
{
  bool valid;
  int n, nlow, nown, ch;
  // the particles that the kernel left to the CPU (more than CAP candidates): for each device the numbers of their
  // records; the forces put what the CPU made of them into those records and into the lists of the device
  int nflag;
  std::vector<std::vector<int> > frec;
};
BState g_bs = {false, 0, 0, 0, 0, 0, std::vector<std::vector<int> >()};

// where the lists made by the last call are on the devices and where they go in list() of the caller, for gpu_sph_getlists_
struct LState
{
  bool valid;
  int ch;
  std::vector<std::vector<long long> > base, tot;   // for each device and chunk: where in list() it goes (from 0) and how many entries
};
LState g_ls = {false, 0, std::vector<std::vector<long long> >(), std::vector<std::vector<long long> >()};

template <class T> void dfree(T *&p) { if (p) { cudaFree(p); p = 0; } }

void freedev(Dev &d)
{
  cudaSetDevice(d.id);
  dfree(d.x); dfree(d.y); dfree(d.z); dfree(d.am); dfree(d.u); dfree(d.hp); dfree(d.divv); dfree(d.bbpart); dfree(d.bb);
  dfree(d.keya); dfree(d.keyb); dfree(d.ida); dfree(d.idb); dfree(d.rec); dfree(d.frec); dfree(d.box);
  dfree(d.oflag); dfree(d.oscan); dfree(d.tgt); dfree(d.cnt); dfree(d.out); dfree(d.list); dfree(d.ccnt); dfree(d.tmp);
  if (d.hout) { cudaFreeHost(d.hout); d.hout = 0; }
  if (d.hcc) { cudaFreeHost(d.hcc); d.hcc = 0; }
  dfree(d.hrecd); dfree(d.fcntd); dfree(d.flod);
  if (d.hrech) { cudaFreeHost(d.hrech); d.hrech = 0; }
  if (d.fcnth) { cudaFreeHost(d.fcnth); d.fcnth = 0; }
  if (d.floh) { cudaFreeHost(d.floh); d.floh = 0; }
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
    dmalloc(d.hrecd, ocap); dmalloc(d.fcntd, 4); dmalloc(d.flod, FLCAP);
    dmalloc(d.list, ((size_t)ocap + chcap) * CAP);
    dmalloc(d.ccnt, 256);
    size_t b1 = 0, b2 = 0;
    CK(cub::DeviceRadixSort::SortPairs(NULL, b1, (unsigned long long *)0, (unsigned long long *)0, (int *)0, (int *)0, ncap, 0, 63));
    CK(cub::DeviceScan::ExclusiveSum(NULL, b2, (int *)0, (int *)0, ncap));
    d.tmpbytes = std::max(b1, b2) + 256;
    CK(cudaMalloc(&d.tmp, d.tmpbytes));
    CK(cudaMallocHost((void **)&d.hout, (size_t)ocap * sizeof(SphOut)));
    CK(cudaMallocHost((void **)&d.hcc, 256 * sizeof(int)));
    CK(cudaMallocHost((void **)&d.hrech, (size_t)ocap * sizeof(HRec)));
    CK(cudaMallocHost((void **)&d.fcnth, 4 * sizeof(int)));
    CK(cudaMallocHost((void **)&d.floh, FLCAP * sizeof(int2)));
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

// room for the Verlet lists of nown particles of the rank, n in all; false (and no lists) if the memory is not there
bool vl_reserve(Dev &d, int n, int nown)
{
  CK(cudaSetDevice(d.id));
  if ((size_t)nown > d.vnown || (size_t)n > d.vn || !d.vl_len)
  {
    cudaStreamSynchronize(d.st);
    dfree(d.vl_len); dfree(d.vl_idx); dfree(d.vnfall); dfree(d.vl_R); dfree(d.vxr); dfree(d.vyr); dfree(d.vzr); dfree(d.vdmax); dfree(d.vpart);
    d.vvalid = false;
    const size_t ocap = (size_t)nown + nown / 16 + 1024, ncap = (size_t)n + n / 16 + 1024;
    try
    {
      dmalloc(d.vl_len, ocap); dmalloc(d.vl_idx, ocap * CAP); dmalloc(d.vl_R, ocap);
      dmalloc(d.vxr, ncap); dmalloc(d.vyr, ncap); dmalloc(d.vzr, ncap);
      dmalloc(d.vdmax, 4); dmalloc(d.vpart, 256); dmalloc(d.vnfall, 4);
    }
    catch (const SphError &)
    {
      dfree(d.vl_len); dfree(d.vl_idx); dfree(d.vnfall); dfree(d.vl_R); dfree(d.vxr); dfree(d.vyr); dfree(d.vzr); dfree(d.vdmax); dfree(d.vpart);
      cudaGetLastError();
      return false;
    }
    d.vnown = ocap;
    d.vn = ncap;
  }
  return true;
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
  const bool mk = (&d == &g_dev[0]);
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
  bool skiplist;             // the lists stay on the device (they are copied by gpu_sph_getlists_ if the CPU needs them)
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
void drain(int nown, int nlow, int ch, const Host &H, int *ltot, int *ierr, long *nflagged, std::vector<std::vector<int> > *flagrec = 0)
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
  std::vector<std::vector<long long> > base(ndev, std::vector<long long>(nchmax + 1, 0)), tots(ndev, std::vector<long long>(nchmax + 1, 0));
  g_ls.valid = false;
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
        if (s.flag)
        {
          if (flagrec)
          {
#pragma omp critical(sphgpu_flagrec)
            (*flagrec)[k].push_back(o);
          }
          continue;
        }
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
        if (!H.skiplist) CK(cudaMemcpyAsync(H.list + hb, d.list + (size_t)c * ch * CAP, (size_t)tot * 4, cudaMemcpyDeviceToHost, d.st2));
      }
      CK(cudaEventRecord(d.evC[c], d.st2));
      tots[k][c] = tot;
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
  g_ls.valid = true; g_ls.ch = ch; g_ls.base = base; g_ls.tot = tots;
  if (nchmax > 0) spread(nchmax - 1);
  *ltot = (int)hb;
  if (nflagged) *nflagged = nfl;
}


// ---------------------------------------------------------------------------------------------------- the results, in the order of the particles (one device)
// With one device and the lists kept on it (skiplist), the copies of the results are queued behind the kernels before the CPUs start their own work, with no call of the CPU in
// between: a kernel puts the records in the order of the particles (dev_unperm), and what comes back is what the CPU reads, 24 bytes for each particle instead of 64, which the CPU
// stores in order instead of one place for each record in a random order.  The rest of the records is what only the pair sums of the CPU read (omega, 0mega, psi, wn and the places
// of the lists: nothing else on the CPU reads them, no dump writes them): the arrays of the CPU are not written with it, and fetch_records makes them as they were if the CPU needs
// them (it needs them where it needs the lists: getlists, and the forces when they are not made from what the device holds).  The values the arrays get are the same as the ones drain gives.
struct LazyRec
{
  bool pending = false;
  Host H;
  int nown = 0, nlow = 0, ch = 0;
  std::vector<long long> base;   // where the lists of each chunk go in list()
};
LazyRec g_lazy;

bool fast_drain_enabled()
{
  static const bool e = getenv("SPHGPU_SLOWDRAIN") == NULL;
  return e;
}

void queue_results(Dev &d, int nown, int nlow, int nch)
{
  CK(cudaSetDevice(d.id));
  for (int c = 0; c < nch; c++) CK(cudaStreamWaitEvent(d.st2, d.evS[c], 0));
  CK(cudaMemcpyAsync(d.hcc, d.ccnt, (size_t)nch * sizeof(int), cudaMemcpyDeviceToHost, d.st2));
  CK(cudaMemsetAsync(d.fcntd, 0, sizeof(int), d.st2));
  dev_unperm<<<nblocks(nown, 256), 256, 0, d.st2>>>(nown, nlow, d.out, d.hrecd, d.fcntd, d.flod);
  CK(cudaGetLastError());
  CK(cudaMemcpyAsync(d.hrech, d.hrecd, (size_t)nown * sizeof(HRec), cudaMemcpyDeviceToHost, d.st2));
  CK(cudaMemcpyAsync(d.fcnth, d.fcntd, sizeof(int), cudaMemcpyDeviceToHost, d.st2));
  CK(cudaMemcpyAsync(d.floh, d.flod, (size_t)FLCAP * sizeof(int2), cudaMemcpyDeviceToHost, d.st2));
  CK(cudaEventRecord(d.evC[0], d.st2));
}

// the same as drain for one device with the results queued by queue_results
void drain_fast(Dev &d, int nown, int nlow, int ch, int nch, const Host &H, int *ltot, int *ierr, long *nflagged, std::vector<std::vector<int> > *flagrec)
{
  CK(cudaSetDevice(d.id));
  CK(cudaEventSynchronize(d.evC[0]));
  std::vector<std::vector<long long> > base(1, std::vector<long long>(nch + 1, 0)), tots(1, std::vector<long long>(nch + 1, 0));
  long long hb = 0;
  bool full = false;
  for (int c = 0; c < nch; c++)
  {
    const long long tot = d.hcc[c];
    base[0][c] = hb;
    if ((size_t)(hb + tot) > H.listcap) full = true;
    tots[0][c] = tot;
    hb += tot;
  }
  CK(cudaStreamSynchronize(d.st2));
  CK(cudaStreamSynchronize(d.st));
  for (int q = 0; q < 4; q++) CK(cudaStreamSynchronize(d.sc[q]));
  g_ls.valid = false;
  g_lazy.pending = false;
  if (full) { *ierr = 2; return; }
  g_ls.valid = true; g_ls.ch = ch; g_ls.base = base; g_ls.tot = tots;
  long nfl = 0;
  const HRec *r = d.hrech;
#pragma omp parallel for schedule(static) reduction(+ : nfl)
  for (int q = 0; q < nown; q++)
  {
    const HRec s = r[q];
    H.sflag[q] = s.flag;
    nfl += s.flag;
    if (s.flag) continue;
    const int i = nlow + q;
    H.hp[i] = s.hp;
    H.rho[i] = s.rho;
    H.zeta[i] = 0.0;
    H.nn[i] = s.nn;
  }
  // the records left to the CPU: where they are in the order of the curve, and the particle of each (the forces read it from the pinned copy of the records)
  const int fc = d.fcnth[0];
  if (fc != nfl) throw SphError(3, "the number of particles left to the CPU is not what the results say");
  if (fc <= FLCAP)
  {
    for (int j = 0; j < fc; j++)
    {
      if (flagrec) (*flagrec)[0].push_back(d.floh[j].x);
      d.hout[d.floh[j].x].idx = d.floh[j].y;
    }
  }
  else
  {
    CK(cudaMemcpy(d.hout, d.out, (size_t)nown * sizeof(SphOut), cudaMemcpyDeviceToHost));
    for (int o = 0; o < nown; o++)
      if (d.hout[o].flag && flagrec) (*flagrec)[0].push_back(o);
  }
  *ltot = (int)hb;
  if (nflagged) *nflagged = nfl;
  g_lazy.pending = true; g_lazy.H = H; g_lazy.nown = nown; g_lazy.nlow = nlow; g_lazy.ch = ch; g_lazy.base = base[0];
}

// The rest of the records of the last call, to the arrays of the CPU, as drain gives them (the particles left to the CPU are not touched).
void fetch_records()
{
  if (!g_lazy.pending) return;
  g_lazy.pending = false;
  Dev &d = g_dev[0];
  CK(cudaSetDevice(d.id));
  const int nown = g_lazy.nown, ch = g_lazy.ch;
  const Host &H = g_lazy.H;
  CK(cudaMemcpy(d.hout, d.out, (size_t)nown * sizeof(SphOut), cudaMemcpyDeviceToHost));
#pragma omp parallel for schedule(static)
  for (int o = 0; o < nown; o++)
  {
    const SphOut s = d.hout[o];
    if (s.flag) continue;
    const int i = s.idx;
    H.omega[i] = s.omega;
    H.om0[i] = s.om0;
    H.psi[i] = s.psi;
    if (H.u[i] != 0.0) H.wn[i] = s.wn;
    H.first[i] = (int)(g_lazy.base[o / ch] + s.first);
  }
}


// ====================================================================================================================== the forces
// Bit 2 of ngpusph: the pair sums of uvdots (balAV3.f: getderivs, uvchunk and uvpairs) on the GPU, for nav=0..3 (the Riemann solver of
// nav=4 stays on the CPU), without periodic boxes and with dynhco=0.  Every term is formed here with the operations, in the order and
// with the tables of the CPU code (one multiplication or addition at a time, --fmad=false), so that each term of a pair is the number the
// CPU computes; the one difference is the root of the viscosity term, r2**0.5d0, which gfortran compiles as a call to pow and which is
// the correctly rounded sqrt here (it differs from pow's result in the last bit now and then).
//
// How the sums are made.  The CPU makes the term of the pair (i,j) for each j in the list of i and adds it to i (the "gather" part, in
// the order of the list) and to j (the "scatter" part, in a buffer of the chunk of pb_chunk=512 consecutive particles i, the buffers being
// added to j one chunk after the other: pairbuf.f90).  Here the same terms are made once, by one thread for each particle i (k_fpair:
// div v and curl v, then the pair loop, which keeps the terms of the pairs, csij and csijgc, in a place of its own for each entry of the
// lists).  The entries are then sorted by the particle they act on (a stable radix sort by j of the entries listed by i: for each j its
// entries come out in the order of i) and one thread for each j (k_fscat) adds them in that order, chunk by chunk, to what the gather
// part gave, exactly as the CPU does.  There is no atomic operation on floating-point numbers and the order of every sum is fixed by the
// index of the particles, whatever the order in which the threads run: the same bits in every run, with any number of devices (a device
// takes whole chunks).  With the viscosity off, or where no term goes through the root of the viscosity, the sums are the bits of the CPU.
constexpr int PBCHUNK = 512;   // pb_chunk of pairbuf.f90
#ifndef SPH_FTHR
#define SPH_FTHR 128
#endif
#ifndef SPH_FMINB
#define SPH_FMINB 1
#endif
constexpr int FTHR = SPH_FTHR;   // threads of a block of the kernels of the forces (-DSPH_FTHR, -DSPH_FMINB: tests of the launch)

struct FPar
{
  double alpha1, beta1, hfloor, gam, ctab, c4;
  int nav3;          // nav=3 (else the other viscosities of uvchunk)
  int nintvar2;      // nintvar=2
  int relaxnoheat;   // nrelax=1 and t<treloff: no heating by the viscosity
  int udotother;     // nrelax=0 or t>=treloff
  int ngr;           // the gravity correction is on
  int blocks;        // block steps: only the active particles take part
};

struct FOut
{
  double udot, uij, divv;
  int p;      // the particle (from 0); -1: no particle at this place of the work (the device does not take it)
  int pad;
};

// what the pair loop of a particle of the work needs: its smoothing length and sums, where its list starts in the list of the device and
// how long it is (valid=0: this device does not do this place of the work)
struct TPar
{
  double hp, rho, om, om0, psi;
  int p, nn, off, valid;
};

struct FView
{
  int n, nts, t0;                  // particles; targets of this device and the first one of them (in the arrays of all targets)
  const D4 *pos;                   // x y z mass
  const D4 *vel;                   // vx vy vz u
  const double *por2;
  const TPar *tpar;                // for each place of the work of this device (or, with inv, nothing)
  const int *inv;                  // with the records of the smoothing-length kernel: the record of each place of the work (-1: not done by this device)
  const SphOut *sout;              // those records
  int sch;                         // and the places of a chunk of its lists
  const int *list;
  const int *act;
  const int *tp, *tiw;             // the particle (from 1) and its place in the work of the rank (from 1), for all targets
  const int *eoff;                 // where the entries of the targets of this device start (nts+1 of them)
  const double *dw, *dg;
  const float *dwf, *dgf;          // the same tables in single precision (ngpusphprec=1)
  int prec;                        // 1: the terms of the pairs in single precision
  FPar P;
  double2 *c;                      // csij and csijgc of every entry
  int *key, *et, *val;             // for every entry: the particle it acts on (n: none), its target, its number
  const int *skey, *sval;          // the same sorted by key
  int nent;
  double *acc;                     // 6 arrays of n: the gather sums, then the sums with the scatter parts
  FOut *out;
  int *nlarge;
};

// the data of place tl of the work of this device: false if the device does not do it
__device__ __forceinline__ bool load_place(const FView &V, int tl, TPar &q)
{
  if (!V.inv)
  {
    q = V.tpar[tl];
    return q.valid != 0;
  }
  const int o = V.inv[tl];
  if (o < 0) return false;
  const SphOut s = V.sout[o];
  q.hp = s.hp; q.rho = s.rho; q.om = s.omega; q.om0 = s.om0; q.psi = s.psi;
  q.p = s.idx; q.nn = s.nn; q.off = (o / V.sch) * V.sch * CAP + s.first; q.valid = 1;
  return true;
}

// the term of the pairs of each particle p of the work: div v and curl v of p, then csij and csijgc of each entry of its list
__global__ void __launch_bounds__(FTHR, SPH_FMINB) k_fpair(FView V)
{
  const int tl = blockIdx.x * FTHR + threadIdx.x;
  if (tl >= V.nts) return;
  const int t = V.t0 + tl;
  TPar tq;
  if (!load_place(V, tl, tq))
  {
    FOut o;
    o.udot = 0.0; o.uij = 0.0; o.divv = 0.0; o.p = -1; o.pad = 0;
    V.out[tl] = o;
    return;
  }
  const int p = tq.p;
  const FPar &P = V.P;
  const int nn = tq.nn, fo = tq.off;
  const int e0 = V.eoff[tl];
  const D4 pp = V.pos[p], pv = V.vel[p];
  const double hpi = tq.hp, rhoi = tq.rho, por2i = V.por2[p], omi = tq.om, om0i = tq.om0, psii = tq.psi;
  const bool actp = P.blocks ? (V.act[p] != 0) : true;
  const double h2 = hpi * hpi;
  const double h5 = hpi * (h2 * h2);
  // getderivs.f
  double cx = 0.0, cy = 0.0, cz = 0.0, dvv = 0.0;
  if (rhoi > 0.0)
  {
    const double ctoh2 = P.ctab / h2;
    for (int k = 0; k < nn; k++)
    {
      const int j = V.list[fo + k] - 1;
      const D4 vj = V.vel[j];
      if (vj.w != 0.0)
      {
        const D4 pj = V.pos[j];
        const double dx = pp.x - pj.x, dy = pp.y - pj.y, dz = pp.z - pj.z;
        const double r2 = (dx * dx + dy * dy) + dz * dz;
        double xq = ctoh2 * r2;
        xq = xq < P.c4 ? xq : P.c4;
        const int it = (int)xq;
        const double fq = xq - (double)it;
        const double dwi = V.dw[it] + fq * (V.dw[it + 1] - V.dw[it]);
        const double dvx = pv.x - vj.x, dvy = pv.y - vj.y, dvz = pv.z - vj.z;
        const double cxi = dz * dvy - dy * dvz;
        const double cyi = dx * dvz - dz * dvx;
        const double czi = dy * dvx - dx * dvy;
        cx = cx + (pj.w * dwi) * cxi;
        cy = cy + (pj.w * dwi) * cyi;
        cz = cz + (pj.w * dwi) * czi;
        const double divin = (dx * dvx + dy * dvy) + dz * dvz;
        const double tij = dwi * divin;
        dvv = dvv - pj.w * tij;
      }
    }
    const double den = rhoi * h5;
    cx = cx / den;
    cy = cy / den;
    cz = cz / den;
    dvv = dvv / den;
  }
  const double abscurl = sqrt((cx * cx + cy * cy) + cz * cz);
  const double ci2 = (P.gam * por2i) * rhoi;
  const double ci = sqrt(ci2);
  // uvchunk
  const double hpt = hpi - P.hfloor;
  const double h2t = hpt * hpt;
  const double invh = 1.0 / hpi, invh2 = 1.0 / h2, invh2t = 1.0 / h2t, invh5 = 1.0 / h5;
  const bool gasp = pv.w != 0.0;
  double gx = 0.0, gy = 0.0, gz = 0.0, ggx = 0.0, ggy = 0.0, ggz = 0.0, ud = 0.0, um = 0.0;
  for (int k = 0; k < nn; k++)
  {
    const int j = V.list[fo + k] - 1;
    const D4 vj = V.vel[j];
    int keyv = j;
    double csij = 0.0, csijgc = 0.0;
    bool term = gasp && vj.w != 0.0;
    if (P.blocks)
    {
      const bool actj = V.act[j] != 0;
      if (!(actp || actj)) term = false;   // not a term of this substep (the CPU skips it)
      if (!actj) keyv = V.n;               // the particle it acts on does not keep what it gets
    }
    if (term)
    {
      const D4 pj = V.pos[j];
      const double dx = pp.x - pj.x, dy = pp.y - pj.y, dz = pp.z - pj.z;
      const double dvx = pv.x - vj.x, dvy = pv.y - vj.y, dvz = pv.z - vj.z;
      const double r2 = (dx * dx + dy * dy) + dz * dz;
      const double q4 = r2 * invh2;
      const double xq = P.ctab * (q4 < 4.0 ? q4 : 4.0);
      const int it = (int)xq;
      const double fq = xq - (double)it;
      double diwkin = 0.0;
      if (r2 < 4.0 * h2t)
      {
        const double xqt = (P.ctab * r2) * invh2t;
        const int itt = (int)xqt;
        const double fqt = xqt - (double)itt;
        diwkin = (V.dg[itt] + fqt * (V.dg[itt + 1] - V.dg[itt])) * invh2t;
      }
      const double dwijin = (V.dw[it] + fq * (V.dw[it + 1] - V.dw[it])) * invh5;
      const double divin = (dx * dvx + dy * dvy) + dz * dvz;
      double pijin = 0.0;
      if (divin < 0.0)
      {
        const double fi = fabs(dvv) / ((fabs(dvv) + abscurl) + (0.00001 * ci) * invh);
        double udbij = 0.0;
        if (r2 > 0.0)
        {
          if (P.nav3) udbij = divin / (ci * sqrt(r2)) * fi;
          else udbij = (hpi * divin) / (ci * r2) * fi;
        }
        if (fabs(udbij) > 1.0e8) atomicAdd(V.nlarge, 1);
        pijin = (por2i * (P.alpha1 + P.beta1 * udbij)) * udbij;
      }
      const double invbonnet = 1.0 / (om0i * pj.w);
      csij = por2i * (-(omi * invbonnet * diwkin) + dwijin);
      csijgc = (0.5 * pijin) * dwijin;
      const double csiju = P.relaxnoheat ? csij : csij + csijgc;
      if (P.nintvar2) ud = ud + (pj.w * csiju) * divin;
      else if (P.udotother) ud = ud + ((pj.w * pijin) * dwijin) * divin;
      if (actp)
      {
        const double uv = (por2i + 0.5 * pijin) * rhoi;
        um = um > uv ? um : uv;
      }
      if (P.ngr) csijgc = csijgc - ((0.5 * diwkin) * psii) * invbonnet;
      const double sx = csij * dx, sy = csij * dy, sz = csij * dz;
      const double sgx = csijgc * dx, sgy = csijgc * dy, sgz = csijgc * dz;
      gx = gx - pj.w * sx;
      gy = gy - pj.w * sy;
      gz = gz - pj.w * sz;
      ggx = ggx - pj.w * sgx;
      ggy = ggy - pj.w * sgy;
      ggz = ggz - pj.w * sgz;
    }
    else
    {
      // no term (a point particle, or a pair that the substep leaves out): the entry does not go to the sums on the neighbours
      keyv = V.n;
    }
    V.c[e0 + k] = make_double2(csij, csijgc);
    V.key[e0 + k] = keyv;
    V.et[e0 + k] = t;
    V.val[e0 + k] = e0 + k;
  }
  const size_t n = (size_t)V.n;
  V.acc[0 * n + p] = gx;
  V.acc[1 * n + p] = gy;
  V.acc[2 * n + p] = gz;
  V.acc[3 * n + p] = ggx;
  V.acc[4 * n + p] = ggy;
  V.acc[5 * n + p] = ggz;
  FOut o;
  o.udot = ud; o.uij = um; o.divv = dvv; o.p = p; o.pad = 0;
  V.out[tl] = o;
}

// The same as k_fpair with the terms of the pairs in single precision (ngpusphprec=1): div v and curl v stay in double precision
__global__ void __launch_bounds__(FTHR, SPH_FMINB) k_fpair_mix(FView V)
{
  const int tl = blockIdx.x * FTHR + threadIdx.x;
  if (tl >= V.nts) return;
  const int t = V.t0 + tl;
  TPar tq;
  if (!load_place(V, tl, tq))
  {
    FOut o;
    o.udot = 0.0; o.uij = 0.0; o.divv = 0.0; o.p = -1; o.pad = 0;
    V.out[tl] = o;
    return;
  }
  const int p = tq.p;
  const FPar &P = V.P;
  const int nn = tq.nn, fo = tq.off;
  const int e0 = V.eoff[tl];
  const D4 pp = V.pos[p], pv = V.vel[p];
  const double hpi = tq.hp, rhoi = tq.rho, por2i = V.por2[p], omi = tq.om, om0i = tq.om0, psii = tq.psi;
  const bool actp = P.blocks ? (V.act[p] != 0) : true;
  const double h2 = hpi * hpi;
  const double h5 = hpi * (h2 * h2);
  // getderivs.f
  double cx = 0.0, cy = 0.0, cz = 0.0, dvv = 0.0;
  if (rhoi > 0.0)
  {
    const double ctoh2 = P.ctab / h2;
    for (int k = 0; k < nn; k++)
    {
      const int j = V.list[fo + k] - 1;
      const D4 vj = V.vel[j];
      if (vj.w != 0.0)
      {
        const D4 pj = V.pos[j];
        const double dx = pp.x - pj.x, dy = pp.y - pj.y, dz = pp.z - pj.z;
        const double r2 = (dx * dx + dy * dy) + dz * dz;
        double xq = ctoh2 * r2;
        xq = xq < P.c4 ? xq : P.c4;
        const int it = (int)xq;
        const double fq = xq - (double)it;
        const double dwi = V.dw[it] + fq * (V.dw[it + 1] - V.dw[it]);
        const double dvx = pv.x - vj.x, dvy = pv.y - vj.y, dvz = pv.z - vj.z;
        const double cxi = dz * dvy - dy * dvz;
        const double cyi = dx * dvz - dz * dvx;
        const double czi = dy * dvx - dx * dvy;
        cx = cx + (pj.w * dwi) * cxi;
        cy = cy + (pj.w * dwi) * cyi;
        cz = cz + (pj.w * dwi) * czi;
        const double divin = (dx * dvx + dy * dvy) + dz * dvz;
        const double tij = dwi * divin;
        dvv = dvv - pj.w * tij;
      }
    }
    const double den = rhoi * h5;
    cx = cx / den;
    cy = cy / den;
    cz = cz / den;
    dvv = dvv / den;
  }
  const double abscurl = sqrt((cx * cx + cy * cy) + cz * cz);
  const double ci2 = (P.gam * por2i) * rhoi;
  const double ci = sqrt(ci2);
  // uvchunk with the terms of the pairs in single precision: the differences of positions and velocities are made in double and rounded to float, everything else of the
  // term is float (the tables too), and what is added up (the accelerations, du/dt, the maximum) is double
  const double hpt = hpi - P.hfloor;
  const double h2t = hpt * hpt;
  const double invh = 1.0 / hpi, invh2 = 1.0 / h2, invh2t = 1.0 / h2t, invh5 = 1.0 / h5;
  const float invh2f = (float)invh2, invh5f = (float)invh5, invh2tf = (float)invh2t, invhf = (float)invh, h2t4f = (float)(4.0 * h2t), hpif = (float)hpi;
  const float por2f = (float)por2i, cif = (float)ci, advf = (float)fabs(dvv), abscf = (float)abscurl, omf = (float)omi, om0f = (float)om0i, psif = (float)psii;
  const float a1f = (float)P.alpha1, b1f = (float)P.beta1, rhof = (float)rhoi;
  const bool gasp = pv.w != 0.0;
  double gx = 0.0, gy = 0.0, gz = 0.0, ggx = 0.0, ggy = 0.0, ggz = 0.0, ud = 0.0, um = 0.0;
  for (int k = 0; k < nn; k++)
  {
    const int j = V.list[fo + k] - 1;
    const D4 vj = V.vel[j];
    int keyv = j;
    double csij_d = 0.0, csijgc_d = 0.0;
    bool term = gasp && vj.w != 0.0;
    if (P.blocks)
    {
      const bool actj = V.act[j] != 0;
      if (!(actp || actj)) term = false;
      if (!actj) keyv = V.n;
    }
    if (term)
    {
      const D4 pj = V.pos[j];
      const double dx = pp.x - pj.x, dy = pp.y - pj.y, dz = pp.z - pj.z;
      const float dxf = (float)dx, dyf = (float)dy, dzf = (float)dz;
      const float dvxf = (float)(pv.x - vj.x), dvyf = (float)(pv.y - vj.y), dvzf = (float)(pv.z - vj.z);
      const float amj = (float)pj.w;
      const float r2 = (dxf * dxf + dyf * dyf) + dzf * dzf;
      const float q4 = r2 * invh2f;
      const double xq = P.ctab * (double)(q4 < 4.0f ? q4 : 4.0f);
      const int it = (int)xq;
      const float fq = (float)(xq - (double)it);
      float diwkin = 0.0f;
      if (r2 < h2t4f)
      {
        const double xqt = P.ctab * (double)(r2 * invh2tf);
        const int itt = (int)xqt;
        const float fqt = (float)(xqt - (double)itt);
        diwkin = (V.dgf[itt] + fqt * (V.dgf[itt + 1] - V.dgf[itt])) * invh2tf;
      }
      const float dwijin = (V.dwf[it] + fq * (V.dwf[it + 1] - V.dwf[it])) * invh5f;
      const float divin = (dxf * dvxf + dyf * dvyf) + dzf * dvzf;
      float pijin = 0.0f;
      if (divin < 0.0f)
      {
        const float fi = advf / ((advf + abscf) + (0.00001f * cif) * invhf);
        float udbij = 0.0f;
        if (r2 > 0.0f)
        {
          if (P.nav3) udbij = divin / (cif * sqrtf(r2)) * fi;
          else udbij = (hpif * divin) / (cif * r2) * fi;
        }
        if (fabsf(udbij) > 1.0e8f) atomicAdd(V.nlarge, 1);
        pijin = (por2f * (a1f + b1f * udbij)) * udbij;
      }
      const float invbonnet = 1.0f / (om0f * amj);
      float csij = por2f * (-(omf * invbonnet * diwkin) + dwijin);
      float csijgc = (0.5f * pijin) * dwijin;
      const float csiju = P.relaxnoheat ? csij : csij + csijgc;
      if (P.nintvar2) ud = ud + (double)((amj * csiju) * divin);
      else if (P.udotother) ud = ud + (double)(((amj * pijin) * dwijin) * divin);
      if (actp)
      {
        const float uv = (por2f + 0.5f * pijin) * rhof;
        um = um > (double)uv ? um : (double)uv;
      }
      if (P.ngr) csijgc = csijgc - ((0.5f * diwkin) * psif) * invbonnet;
      csij_d = (double)csij;
      csijgc_d = (double)csijgc;
      const double sx = csij_d * dx, sy = csij_d * dy, sz = csij_d * dz;
      const double sgx = csijgc_d * dx, sgy = csijgc_d * dy, sgz = csijgc_d * dz;
      gx = gx - pj.w * sx;
      gy = gy - pj.w * sy;
      gz = gz - pj.w * sz;
      ggx = ggx - pj.w * sgx;
      ggy = ggy - pj.w * sgy;
      ggz = ggz - pj.w * sgz;
    }
    else
    {
      keyv = V.n;
    }
    V.c[e0 + k] = make_double2(csij_d, csijgc_d);
    V.key[e0 + k] = keyv;
    V.et[e0 + k] = t;
    V.val[e0 + k] = e0 + k;
  }
  const size_t n = (size_t)V.n;
  V.acc[0 * n + p] = gx;
  V.acc[1 * n + p] = gy;
  V.acc[2 * n + p] = gz;
  V.acc[3 * n + p] = ggx;
  V.acc[4 * n + p] = ggy;
  V.acc[5 * n + p] = ggz;
  FOut o;
  o.udot = ud; o.uij = um; o.divv = dvv; o.p = p; o.pad = 0;
  V.out[tl] = o;
}

// the first entry with key >= q in the sorted keys
__device__ __forceinline__ int lower_key(const int *key, int m, int q)
{
  int lo = 0, hi = m;
  while (lo < hi)
  {
    const int mid = (lo + hi) >> 1;
    if (key[mid] < q) lo = mid + 1; else hi = mid;
  }
  return lo;
}

// the terms that act on each particle q as the scatter part: for each q the entries that name it, in the order of the work (the
// particle i, then the place in its list), added chunk by chunk to the gather sum that is in acc
__global__ void __launch_bounds__(FTHR, SPH_FMINB) k_fscat(FView V)
{
  const int q = blockIdx.x * FTHR + threadIdx.x;
  if (q >= V.n) return;
  const int lo = lower_key(V.skey, V.nent, q);
  if (lo >= V.nent || V.skey[lo] != q) return;
  const int hi = lower_key(V.skey, V.nent, q + 1);
  const size_t n = (size_t)V.n;
  double tot[6];
  for (int a = 0; a < 6; a++) tot[a] = V.acc[a * n + q];
  double b[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  int chunk = -1;
  const D4 pq = V.pos[q];
  for (int s = lo; s < hi; s++)
  {
    const int e = V.sval[s];
    const int t = V.et[e];
    const int ch = (V.tiw[t] - 1) / PBCHUNK;
    if (ch != chunk)
    {
      if (chunk >= 0) for (int a = 0; a < 6; a++) tot[a] = tot[a] + b[a];
      for (int a = 0; a < 6; a++) b[a] = 0.0;
      chunk = ch;
    }
    const int p = V.tp[t] - 1;
    const D4 pp = V.pos[p];
    const double dx = pp.x - pq.x, dy = pp.y - pq.y, dz = pp.z - pq.z;
    const double2 c = V.c[e];
    const double sx = c.x * dx, sy = c.x * dy, sz = c.x * dz;
    const double sgx = c.y * dx, sgy = c.y * dy, sgz = c.y * dz;
    b[0] = b[0] + pp.w * sx;
    b[1] = b[1] + pp.w * sy;
    b[2] = b[2] + pp.w * sz;
    b[3] = b[3] + pp.w * sgx;
    b[4] = b[4] + pp.w * sgy;
    b[5] = b[5] + pp.w * sgz;
  }
  for (int a = 0; a < 6; a++) tot[a] = tot[a] + b[a];
  for (int a = 0; a < 6; a++) V.acc[a * n + q] = tot[a];
}

// The particles a device touches in the forces with the data of the smoothing-length kernel (useb): its own targets (the gather sums are written for them) and the keys of
// the dispersion (the particles that get a scatter term).  Every other element of its sums is the +0.0 of the clearing, so the host adds only these.  The mask is
// written by many threads with the same value; the list of the touched particles comes out in the order of the index.
__global__ void k_tmark_own(int sl, const SphOut *out, int *mask)
{
  const int o = blockIdx.x * blockDim.x + threadIdx.x;
  if (o >= sl) return;
  mask[out[o].idx] = 1;
}

__global__ void k_tmark_keys(int nent, int n, const int *key, int *mask)
{
  const int e = blockIdx.x * blockDim.x + threadIdx.x;
  if (e >= nent) return;
  const int k = key[e];
  if (k < n) mask[k] = 1;
}

// the touched particles in the order of the index (pos: the exclusive sum of the mask), and how many
__global__ void k_tcompact(int n, const int *mask, const int *pos, int *idx, int *cnt)
{
  const int q = blockIdx.x * blockDim.x + threadIdx.x;
  if (q >= n) return;
  if (mask[q]) idx[pos[q]] = q;
  if (q == n - 1) cnt[0] = pos[q] + mask[q];
}

// the sums of the touched particles, packed (six numbers for each, in the order of the list), and the results of the targets in the order of the records
__global__ void k_tgather(int m, int n, const int *idx, const double *acc, double *v)
{
  const int j = blockIdx.x * blockDim.x + threadIdx.x;
  if (j >= m) return;
  const int q = idx[j];
  for (int a = 0; a < 6; a++) v[(size_t)j * 6 + a] = acc[(size_t)a * n + q];
}

__global__ void k_foutc(int sl, int nlow, const SphOut *sout, const FOut *out, FOut *oc)
{
  const int o = blockIdx.x * blockDim.x + threadIdx.x;
  if (o >= sl) return;
  oc[o] = out[sout[o].idx - nlow];
}

// the records of the particles that the CPU solved put into the records of the device (the places of the records in pidx)
__global__ void k_patch(int np, const int *pidx, const SphOut *rec, SphOut *out)
{
  const int j = blockIdx.x * blockDim.x + threadIdx.x;
  if (j >= np) return;
  out[pidx[j]] = rec[j];
}

__global__ void k_fpack(int n, const double *x, const double *y, const double *z, const double *am, const double *vx, const double *vy,
                        const double *vz, const double *u, D4 *pos, D4 *vel)
{
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= n) return;
  pos[i] = make_d4(x[i], y[i], z[i], am[i]);
  vel[i] = make_d4(vx[i], vy[i], vz[i], u[i]);
}

__global__ void k_fcontig(int nt, int p0, int *tp, int *tiw)
{
  const int t = blockIdx.x * blockDim.x + threadIdx.x;
  if (t >= nt) return;
  tp[t] = p0 + t;
  tiw[t] = t + 1;
}

__global__ void k_fcount(FView V, int *cnt)
{
  const int tl = blockIdx.x * blockDim.x + threadIdx.x;
  if (tl > V.nts) return;
  TPar q;
  cnt[tl] = (tl < V.nts && load_place(V, tl, q)) ? q.nn : 0;
}

#ifndef SPH_FTPAR_SCATTER
// the inverse of the order of the records: the record of each place of the work
__global__ void k_finv(int sl, int nlow, const SphOut *out, int *inv)
{
  const int o = blockIdx.x * blockDim.x + threadIdx.x;
  if (o >= sl) return;
  inv[out[o].idx - nlow] = o;
}
#endif

// the data of the places of the work from arrays indexed by the particle (the CPU's), for the targets t0.. of the arrays of all targets
__global__ void k_ftpar_u(int nts, int t0, const int *tp, const double *hp, const double *rho, const double *om, const double *om0, const double *psi,
                          const int *nn, const int *first, TPar *tpar)
{
  const int tl = blockIdx.x * blockDim.x + threadIdx.x;
  if (tl >= nts) return;
  const int p = tp[t0 + tl] - 1;
  TPar q;
  q.hp = hp[p]; q.rho = rho[p]; q.om = om[p]; q.om0 = om0[p]; q.psi = psi[p];
  q.p = p; q.nn = nn[p]; q.off = first[p]; q.valid = 1;
  tpar[tl] = q;
}

#ifdef SPH_FTPAR_SCATTER
// the same from what the smoothing-length kernel left on this device (the records of its slice, whose lists are in chunks of ch places of
// CAP entries): the place of the work of a particle is its index minus the first of the rank
__global__ void k_ftpar_r(int sl, int nlow, int ch, const SphOut *out, TPar *tpar)
{
  const int o = blockIdx.x * blockDim.x + threadIdx.x;
  if (o >= sl) return;
  const SphOut s = out[o];
  TPar q;
  q.hp = s.hp; q.rho = s.rho; q.om = s.omega; q.om0 = s.om0; q.psi = s.psi;
  q.p = s.idx; q.nn = s.nn; q.off = (o / ch) * ch * CAP + s.first; q.valid = 1;
  tpar[s.idx - nlow] = q;
}
#endif

// the buffers of one device for the forces (the particles of the call are uploaded whole: the arrays of the CPU are the source)
struct FDev
{
  int ncap = 0, ntcap = 0;
  size_t ecap = 0, lcap = 0;
  double *x = 0, *y = 0, *z = 0, *vx = 0, *vy = 0, *vz = 0, *am = 0, *u = 0, *hp = 0, *rho = 0, *por2 = 0, *om = 0, *om0 = 0, *psi = 0;
  D4 *pos = 0, *vel = 0;
  int *nn = 0, *first = 0, *act = 0;
  TPar *tpar = 0;
  int *inv = 0;
  int *hE = 0;           // pinned: the entries of the lists of this device
  double *acc = 0;
  int *tp = 0, *tiw = 0, *eoff = 0, *cnt = 0;
  FOut *out = 0;
  int *list = 0;
  double2 *c = 0;
  int *ka = 0, *kb = 0, *va = 0, *vb = 0, *et = 0;
  double *dw = 0, *dg = 0;
  float *dwf = 0, *dgf = 0;
  void *tmp = 0;
  size_t tmpbytes = 0;
  double *hacc = 0;      // pinned: the 6 arrays of n that come back
  size_t haccn = 0;
  FOut *hout = 0;        // pinned
  size_t houtn = 0;
  int *nlarge = 0, *hnlarge = 0;
  int *tmask = 0, *tpos = 0, *tidx = 0, *tcnt = 0;   // the particles this device touches (see k_tmark_own): mask, its exclusive sum, the list, how many
  int *hidx = 0, *hcnt = 0;                           // pinned: the list and the count
  double *tv = 0;                                     // the sums of the touched particles, six for each (see k_tgather)
  FOut *oc = 0;                                       // the results of the targets in the order of the records
  SphOut *prec = 0;                                   // the records of the particles the CPU solved (see k_patch) and the places of the records
  int *pidx = 0;
  size_t pcap = 0;
};
std::vector<FDev> g_fdev;

void ffree(FDev &f)
{
  double **dd[] = {&f.x, &f.y, &f.z, &f.vx, &f.vy, &f.vz, &f.am, &f.u, &f.hp, &f.rho, &f.por2, &f.om, &f.om0, &f.psi, &f.acc};
  for (size_t k = 0; k < sizeof dd / sizeof dd[0]; k++) dfree(*dd[k]);
  dfree(f.tpar); dfree(f.inv); if (f.hE) { cudaFreeHost(f.hE); f.hE = 0; }
  dfree(f.pos); dfree(f.vel); dfree(f.nn); dfree(f.first); dfree(f.act); dfree(f.tp); dfree(f.tiw); dfree(f.eoff); dfree(f.cnt); dfree(f.out);
  dfree(f.list); dfree(f.c); dfree(f.ka); dfree(f.kb); dfree(f.va); dfree(f.vb); dfree(f.et); dfree(f.tmp); dfree(f.nlarge);
  if (f.hacc) { cudaFreeHost(f.hacc); f.hacc = 0; }
  if (f.hout) { cudaFreeHost(f.hout); f.hout = 0; }
  if (f.hnlarge) { cudaFreeHost(f.hnlarge); f.hnlarge = 0; }
  dfree(f.tmask); dfree(f.tpos); dfree(f.tidx); dfree(f.tcnt); dfree(f.tv); dfree(f.oc); dfree(f.prec); dfree(f.pidx); f.pcap = 0;
  if (f.hidx) { cudaFreeHost(f.hidx); f.hidx = 0; }
  if (f.hcnt) { cudaFreeHost(f.hcnt); f.hcnt = 0; }
  f.ncap = f.ntcap = 0; f.ecap = f.lcap = 0; f.tmpbytes = 0; f.haccn = f.houtn = 0;
}

// The arrays of the CPU that the forces send are not registered with the driver (cudaHostRegister makes whole pages of the memory pinned, and the pages
// at the ends of an array are shared with the arrays next to it in the common blocks: a later copy of one of those, by the gravity for instance, then
// starts in a registered page and ends in a page that is not, which the driver refuses): they are copied by the threads of the rank into a pinned buffer
// of ours, and sent from it.
char *g_stage = 0;
size_t g_stagecap = 0;
struct StageItem
{
  const void *src;
  size_t bytes, off;
};
struct Stager
{
  std::vector<StageItem> items;
  size_t total = 0;
  size_t put(const void *src, size_t bytes)
  {
    const size_t o = total;
    total += (bytes + 255) / 256 * 256;
    StageItem it = {src, bytes, o};
    items.push_back(it);
    return o;
  }
  // allocates the buffer if it is too small and copies everything into it
  void fill()
  {
    if (total > g_stagecap || !g_stage)
    {
      if (g_stage) cudaFreeHost(g_stage);
      g_stage = 0;
      g_stagecap = total + total / 8 + 4096;
      CK(cudaHostAlloc((void **)&g_stage, g_stagecap, cudaHostAllocPortable));
    }
    const size_t blk = 1 << 20;
    std::vector<std::pair<int, size_t> > tasks;
    for (size_t i = 0; i < items.size(); i++)
      for (size_t b = 0; b < items[i].bytes; b += blk) tasks.push_back(std::make_pair((int)i, b));
#pragma omp parallel for schedule(static)
    for (long long t = 0; t < (long long)tasks.size(); t++)
    {
      const StageItem &it = items[tasks[t].first];
      const size_t b = tasks[t].second;
      memcpy(g_stage + it.off + b, (const char *)it.src + b, std::min(blk, it.bytes - b));
    }
  }
};
// the copy of one staged array to a device (what: its name, for the message of an error)
inline void h2d_staged(void *dst, size_t off, size_t bytes, cudaStream_t st, const char *what)
{
  if (!bytes) return;
  const cudaError_t e = cudaMemcpyAsync(dst, g_stage + off, bytes, cudaMemcpyHostToDevice, st);
  if (e != cudaSuccess)
  {
    char b_[256];
    snprintf(b_, sizeof b_, "CUDA error in the copy of %s (%zu bytes) to the device: %s", what, bytes, cudaGetErrorString(e));
    throw SphError(e == cudaErrorMemoryAllocation ? 4 : 3, b_);
  }
}
#define H2S(dst, off, bytes) h2d_staged(dst, off, bytes, st, #dst)

// What the forces need that is known before the smoothing lengths are solved: the velocities and u, which the predictor made (nothing between it and the forces writes them).  They are
// sent, and packed with the positions (k_fpack), while the devices solve (gpu_sph_fearly_, called by the work of gpu_sph_rhoh_), so that the forces send only the pressure term.  g_fe
// says the devices hold them for n particles; the forces use them if they are there and say so by clearing it.  g_estage is the pinned buffer of the copies (not g_stage: the forces fill
// that one while these copies may be on their way).
struct FEarly
{
  bool valid = false;
  int n = 0;
};
FEarly g_fe;
bool g_rhoh_queued = false;     // the call of gpu_sph_rhoh_ that is running has queued its kernels and its results (the work may send)
char *g_estage = 0;
size_t g_estagecap = 0;
bool fe_enabled()
{
  static const bool e = getenv("SPHGPU_NOEARLY") == NULL;
  return e;
}
// SPHGPU_FECHECK=1 (a test): the forces compare the arrays of the CPU with what was sent early, and stop at a difference
bool fe_check()
{
  static const bool e = getenv("SPHGPU_FECHECK") != NULL;
  return e;
}

// the sorted-by-key bits: enough for the keys 0..n (n = no key)
inline int key_bits(int n)
{
  int b = 1;
  while ((1LL << b) <= (long long)n) b++;
  return b;
}

// Copies the lists the last call left on the devices to list() (the places drain would have given them). ierr: 0 ok, 3 error of CUDA or no lists
void getlists(int *list, int *ierr)
{
  if (!g_ls.valid) { *ierr = 3; return; }
  try
  {
    fetch_records();      // (the CPU that needs the lists needs the sums of the same particles)
    const int ndev = (int)g_dev.size();
    for (int k = 0; k < ndev; k++)
    {
      Dev &d = g_dev[k];
      CK(cudaSetDevice(d.id));
      for (size_t c = 0; c < g_ls.tot[k].size(); c++)
      {
        const long long tot = g_ls.tot[k][c];
        if (tot > 0) CK(cudaMemcpyAsync(list + g_ls.base[k][c], d.list + (size_t)c * g_ls.ch * CAP, (size_t)tot * 4, cudaMemcpyDeviceToHost, d.st2));
      }
    }
    for (int k = 0; k < ndev; k++)
    {
      CK(cudaSetDevice(g_dev[k].id));
      CK(cudaStreamSynchronize(g_dev[k].st2));
    }
    CK(cudaSetDevice(g_dev[0].id));
  }
  catch (const SphError &e)
  {
    fprintf(stderr, "StarSmasher (SPH on the GPU): %s\n", e.what());
    *ierr = e.code;
  }
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
      // SPHGPU_NDEV=K (a test of the code, not for a run): the SPH uses K logical devices, dealt out round the nd physical ones of the process
      // (logical device k is on physical device first + k % nd), each with its own stream and buffers.  The gravity keeps the nd physical ones.
      // It gives, on one GPU, the bits of K devices (the times are not those of K GPUs).
      int nlog = nd;
      if (getenv("SPHGPU_NDEV") && atoi(getenv("SPHGPU_NDEV")) > 0) nlog = std::min(atoi(getenv("SPHGPU_NDEV")), 64);
      g_dev.resize(nlog);
      for (int k = 0; k < nlog; k++)
      {
        g_dev[k].id = first + k % nd;
        CK(cudaSetDevice(g_dev[k].id));
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

  // The lists of the last gpu_sph_rhoh_ call, which it left on the devices (skiplist=1), copied to list().  ierr: 0 ok, 3 error of CUDA or no lists.
  void gpu_sph_getlists_(int *list, int *ierr)
  {
    *ierr = 0;
    getlists(list, ierr);
  }

  // number of devices of this process
  void gpu_sph_ndev_(int *nd) { *nd = (int)g_dev.size(); }

  // What the forces will need and the CPU knows already, sent while the device solves the smoothing lengths: the velocities of the n particles (and u, which the device holds), packed
  // with the positions the devices hold into the arrays the pair kernels read.  Called from the work of gpu_sph_rhoh_ (one device, the lists kept on it); does nothing in any other
  // case.  The forces use it if the arrays of the CPU still hold these velocities and u (see FEarly); a call that cannot send leaves nothing changed that the forces read.
  void gpu_sph_fearly_(int *n_, double *vx, double *vy, double *vz, double *u, int *ierr)
  {
    *ierr = 0;
    g_fe.valid = false;
    try
    {
      if (!g_ready || !g_rhoh_queued || g_dev.size() != 1 || !fe_enabled()) return;
      const int n = *n_;
      if (g_fdev.size() != 1 || !g_fdev[0].x || n > g_fdev[0].ncap || !g_fdev[0].dw) return;
      Dev &d = g_dev[0];
      FDev &f = g_fdev[0];
      CK(cudaSetDevice(d.id));
      const size_t bytes = (size_t)n * 8;
      if (4 * bytes > g_estagecap || !g_estage)
      {
        if (g_estage) cudaFreeHost(g_estage);
        g_estage = 0;
        g_estagecap = 4 * bytes + 4 * bytes / 8 + 4096;
        CK(cudaHostAlloc((void **)&g_estage, g_estagecap, cudaHostAllocPortable));
      }
      const double *src[3] = {vx, vy, vz};
      const size_t blk = 1 << 20;
#pragma omp parallel for schedule(static) collapse(2)
      for (int a = 0; a < 3; a++)
        for (long long b = 0; b < (long long)bytes; b += blk)
          memcpy(g_estage + (size_t)a * bytes + b, (const char *)src[a] + b, std::min(blk, bytes - (size_t)b));
      if (fe_check()) memcpy(g_estage + 3 * bytes, u, bytes);     // (for the check of the forces: the u that the device holds is the one of this call)
      cudaStream_t st = d.st;
      CK(cudaMemcpyAsync(f.vx, g_estage, bytes, cudaMemcpyHostToDevice, st));
      CK(cudaMemcpyAsync(f.vy, g_estage + bytes, bytes, cudaMemcpyHostToDevice, st));
      CK(cudaMemcpyAsync(f.vz, g_estage + 2 * bytes, bytes, cudaMemcpyHostToDevice, st));
      CK(cudaMemcpyAsync(f.u, d.u, bytes, cudaMemcpyDeviceToDevice, st));
      k_fpack<<<nblocks(n, 256), 256, 0, st>>>(n, d.x, d.y, d.z, d.am, f.vx, f.vy, f.vz, f.u, f.pos, f.vel);
      CK(cudaGetLastError());
      g_fe.valid = true;
      g_fe.n = n;
    }
    catch (const SphError &e)
    {
      fprintf(stderr, "StarSmasher (SPH on the GPU): %s\n", e.what());
      g_fe.valid = false;
      *ierr = e.code;
    }
  }

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
  // skin > 0 (one device for the rank): Verlet lists, the candidates of the first search kept from the call that built them with a margin of skin of the radius
  // and used while they are complete (the same lists and sums, the same bits); skin <= 0: every call searches from the root.
  // work: a function (or NULL) that is called on the calling thread once the kernels are queued and before their results are collected.
  void gpu_sph_rhoh_(int *n_, int *nlow_, int *nup_, int *nmaxarr, double *x, double *y, double *z, double *am, double *u, double *hp,
                     double *divv, double *rho, double *omega, double *om0, double *psi, double *wn, double *zeta, int *nn, int *first,
                     int *sflag, int *list, int *listcap_, int *ltot, double *dt, double *hfloor, int *nnopt, int *nselfg, double *skin_, int *skiplist_,
                     void (**work)(void), int *ierr)
  {
    *ierr = 0;
    *ltot = 0;
    g_bs.valid = false;
    g_lazy.pending = false;
    g_fe.valid = false;
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
      // Verlet lists: built by this call (2), or used (1), or not at all (0)
      int vmode = 0;
      const double skin = *skin_;
      if (skin > 0.0 && ndev == 1 && nown > 0 && g_dev[0].voff > 0)
      {
        // the lists were built again too soon (see below): they are not used for a number of calls, counted, and then tried again
        g_dev[0].voff--;
        g_dev[0].vvalid = false;
      }
      else if (skin > 0.0 && ndev == 1 && nown > 0 && vl_reserve(g_dev[0], n, nown))
      {
        Dev &d = g_dev[0];
        CK(cudaSetDevice(d.id));
        const bool same = d.vvalid && d.vbn == n && d.vbnlow == nlow && d.vbnown == nown && !d.vforce;
        vmode = same ? 1 : 2;
        if (vmode == 2)
        {
          CK(cudaMemsetAsync(d.vl_len, 0xff, (size_t)nown * sizeof(int), d.st));
          CK(cudaMemcpyAsync(d.vxr, d.x, (size_t)n * 8, cudaMemcpyDeviceToDevice, d.st));
          CK(cudaMemcpyAsync(d.vyr, d.y, (size_t)n * 8, cudaMemcpyDeviceToDevice, d.st));
          CK(cudaMemcpyAsync(d.vzr, d.z, (size_t)n * 8, cudaMemcpyDeviceToDevice, d.st));
          CK(cudaMemsetAsync(d.vdmax, 0, sizeof(double), d.st));
        }
        else
        {
          const int nbm = std::min(256, nblocks(n, 256));
          k_vdmax<<<nbm, 256, 0, d.st>>>(n, d.x, d.y, d.z, d.vxr, d.vyr, d.vzr, d.vpart);
          k_vdmax_final<<<1, 256, 0, d.st>>>(nbm, d.vpart, d.vdmax);
        }
        CK(cudaMemsetAsync(d.vnfall, 0, sizeof(int), d.st));
        d.vforce = false;
      }
      double t1 = 0;
      if (plev) { for (int k = 0; k < ndev; k++) { CK(cudaSetDevice(g_dev[k].id)); CK(cudaStreamSynchronize(g_dev[k].st)); } t1 = now_ms(); }
      std::vector<cudaEvent_t> ev;
      for (int k = 0; k < ndev; k++)
      {
        Dev &d = g_dev[k];
        CK(cudaSetDevice(d.id));
        View V = makeview(d, ndev, k, n, nlow, nup, ph);
        if (vmode)
        {
          V.vmode = vmode; V.px = d.x; V.py = d.y; V.pz = d.z; V.pu = d.u; V.vl_len = d.vl_len; V.vl_idx = d.vl_idx; V.vl_R = d.vl_R;
          V.dmaxp = d.vdmax; V.vskin = skin; V.nfall = d.vnfall;
        }
        const int ns = slice_of(nown, ndev, k);
        if (k == 0) mark("(before the kernels)", d.st);
        if (ns > 0) launch_chunks(d, V, ns, ch, (k == 0 && plev >= 2) ? &ev : 0);
      }
      Host H;
      H.hp = hp; H.rho = rho; H.omega = omega; H.om0 = om0; H.psi = psi; H.wn = wn; H.zeta = zeta; H.u = u;
      H.nn = nn; H.first = first; H.sflag = sflag; H.list = list; H.listcap = listcap; H.skiplist = (*skiplist_ != 0);
      long nfl = 0;
      // The particles left to the CPU do not make the lists come back, and the forces do not leave the data on the devices: they take what the CPU made of those particles
      // to the device that owns them, so that a step with such particles is summed as the other steps are (with one device the bits are those of the plain path: the sums
      // of a step do not depend on how the work is divided when it is not divided)
      const bool keepres = true;
      std::vector<std::vector<int> > flagrec(ndev);
      // The work the caller has for its CPUs while the kernels run (the kernels are queued and the copies of the data are on their way: the
      // work must read the arrays of the particles only, and write nothing that the kernels or the copies after them read)
      // one device, the lists kept on it: the copies of the results are queued now (see queue_results)
      // (with more than 8 chunks the kernels end one after the other over many milliseconds, and drain copies and spreads each chunk while the next is computed: that is the better way)
      const bool fastdrain = fast_drain_enabled() && ndev == 1 && H.skiplist && nown > 0 && keepres && (nown + ch - 1) / ch <= 8;
      if (fastdrain) queue_results(g_dev[0], nown, nlow, (nown + ch - 1) / ch);
      g_rhoh_queued = fastdrain;
      if (work && *work) (*work)();
      g_rhoh_queued = false;
      if (fastdrain) drain_fast(g_dev[0], nown, nlow, ch, (nown + ch - 1) / ch, H, ltot, ierr, &nfl, &flagrec);
      else drain(nown, nlow, ch, H, ltot, ierr, &nfl, keepres ? &flagrec : 0);
      // lists left on the device that the CPU needs after all (particles it has to solve itself): copy them now
      if (*ierr == 0 && *skiplist_ != 0 && nfl > 0 && !keepres) getlists(list, ierr);
      if (vmode)
      {
        Dev &d = g_dev[0];
        CK(cudaSetDevice(d.id));
        if (*ierr == 0)
        {
          int nf = 0;
          double dm = 0.0;
          CK(cudaMemcpy(&nf, d.vnfall, sizeof nf, cudaMemcpyDeviceToHost));
          CK(cudaMemcpy(&dm, d.vdmax, sizeof dm, cudaMemcpyDeviceToHost));
          if (vmode == 2) { d.vvalid = true; d.vbn = n; d.vbnlow = nlow; d.vbnown = nown; d.vused = 0; }
          // the lists are built again when more than 2% of the particles have a list that no longer reaches the radius they need: a geometric criterion
          d.vforce = (vmode == 1) && ((long long)nf * 50 > (long long)nown);
          if (vmode == 1) d.vused++;
          // A list set that served fewer than VMINLIFE calls (counted, not timed) costs more than it saves, because building it is dearer than a plain search: the lists are
          // then left alone for vofflen calls, which doubles each time up to 384 and goes back to 24 when a set lasts; after them they are built again and judged again
          if (d.vforce)
          {
            if (d.vused + 1 < VMINLIFE) { d.voff = d.vofflen; d.vofflen = std::min(2 * d.vofflen, 384); d.vvalid = false; }
            else d.vofflen = 24;
          }
          if (plev) fprintf(stderr, "SPHGPU verlet: %s, largest displacement %.3e, %d of %d lists did not reach%s%s\n", vmode == 2 ? "built" : "used", dm, nf, nown, d.vforce ? " (built again next call)" : "", d.voff > 0 ? " (off)" : "");
        }
        else d.vvalid = false;
      }
      if (*ierr == 0 && (nfl == 0 || keepres))
      {
        g_bs.valid = true; g_bs.n = n; g_bs.nlow = nlow; g_bs.nown = nown; g_bs.ch = ch;
        g_bs.nflag = (int)nfl;
        g_bs.frec.clear();
        if (nfl > 0)
        {
          for (int k = 0; k < ndev; k++) std::sort(flagrec[k].begin(), flagrec[k].end());
          g_bs.frec = flagrec;
        }
      }
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

  // The tables dw/dh (dwtab) and dg (dgtab) of the kernel that the forces read (kernels.f), with the guard entry (ntab+1 each), copied to
  // every device of the process, once.  ierr: 0 ok, 3 error of CUDA, 4 no memory.
  void gpu_sph_ftables_(int *ntab, double *dwtab, double *dgtab, int *ierr)
  {
    *ierr = 0;
    try
    {
      if (!g_ready) throw SphError(3, "gpu_sph_ftables_ called before gpu_sph_setup_");
      if (g_fdev.size() != g_dev.size()) g_fdev.resize(g_dev.size());
      const size_t nt = (size_t)*ntab + 1;
      for (size_t k = 0; k < g_dev.size(); k++)
      {
        FDev &f = g_fdev[k];
        CK(cudaSetDevice(g_dev[k].id));
        dfree(f.dw); dfree(f.dg);
        dfree(f.dwf); dfree(f.dgf);
        dmalloc(f.dw, nt); dmalloc(f.dg, nt); dmalloc(f.dwf, nt); dmalloc(f.dgf, nt);
        {
          std::vector<float> a(nt), b(nt);
          for (size_t q = 0; q < nt; q++) { a[q] = (float)dwtab[q]; b[q] = (float)dgtab[q]; }
          CK(cudaMemcpy(f.dwf, a.data(), nt * sizeof(float), cudaMemcpyHostToDevice));
          CK(cudaMemcpy(f.dgf, b.data(), nt * sizeof(float), cudaMemcpyHostToDevice));
        }
        CK(cudaMemcpy(f.dw, dwtab, nt * sizeof(double), cudaMemcpyHostToDevice));
        CK(cudaMemcpy(f.dg, dgtab, nt * sizeof(double), cudaMemcpyHostToDevice));
      }
      CK(cudaSetDevice(g_dev[0].id));
    }
    catch (const SphError &e)
    {
      fprintf(stderr, "StarSmasher (SPH on the GPU): %s\n", e.what());
      *ierr = e.code;
    }
  }

  // The pair sums of uvdots for the targets tp(1:nt) (particles from 1, in the order of the work of the rank, tiw(t) being the place of
  // target t in that work, from 1; contig=1: the targets are the particles tp(1)..tp(1)+nt-1 and tiw(t)=t), with the lists, smoothing
  // lengths and sums that rho_and_h left in the arrays of the CPU.  On return the arrays of the CPU named in the last line hold, for every
  // particle of 1..n, the sums mvx..mgz of getderivs/uvchunk (the "gather" and "scatter" parts, the same terms in the same order as the
  // CPU), and for each target udot (before the factor of nintvar=1 and 3), uijmax (only where act, or if blocks=0) and divv.
  // work: a function (or NULL) that is called on the calling thread when the work of the devices is queued, before the call waits for them.
  // useb=1: the smoothing lengths and lists were just made by gpu_sph_rhoh_ for the same particles (shared steps, all of them on the
  // GPU): the positions, masses, smoothing lengths, sums and lists that are on the devices are used, and only the velocities, u and
  // the pressure term are sent.  ierr: 0 ok, 3 error of CUDA, 4 no memory; a return with ierr not 0 has changed none of the arrays of
  // the CPU.
  void gpu_sph_forces_(int *n_, int *nmaxarr, int *nt_, int *nloop_, int *contig_, int *tp, int *tiw, double *x, double *y, double *z,
                       double *vx, double *vy, double *vz, double *am, double *u, double *hp, double *rho, double *por2, double *om,
                       double *om0, double *psi, int *nn, int *first, int *list, int *listcap_, int *act, int *blocks_, double *alpha1,
                       double *beta1, double *hfloor, double *gam, int *nav, int *nintvar, int *relaxnoheat, int *udotother, int *ngr,
                       int *useb_, int *prec_, double *mvx, double *mvy, double *mvz, double *mgx, double *mgy, double *mgz, double *udot,
                       double *uijmax, double *divv, int *nlarge, void (**work)(void), int *ierr)
  {
    *ierr = 0;
    *nlarge = 0;
    const int n = *n_, nt = *nt_, nloop = *nloop_, contig = *contig_, blocks = *blocks_, prec = *prec_;
    try
    {
      if (!g_ready) throw SphError(3, "gpu_sph_forces_ called before gpu_sph_setup_");
      const int ndev = (int)g_dev.size();
      if ((int)g_fdev.size() != ndev) g_fdev.resize(ndev);
      for (int k = 0; k < ndev; k++)
        if (!g_fdev[k].dw) throw SphError(3, "the tables of the forces were not sent to the device");
      if (nt <= 0) return;
      const double tstart = now_ms();
      const int plev = prof_level();
      const int p0 = tp[0];
      // the data of the smoothing-length kernel on the devices can be used if it was just made for these particles
      bool useb = (*useb_ != 0) && contig && !blocks && g_bs.valid && g_bs.n == n && g_bs.nown == nt && g_bs.nlow == p0 - 1 && nloop == nt;
      if (!useb) fetch_records();      // (the forces made from the arrays of the CPU need all of the results of the call before)
      // The particles the CPU solved (g_bs.nflag > 0): their records and lists go to the device that owns them, as the GPU would have left
      // them, and the step is summed as a step without such particles.  The lists are placed after the lists of the slice of the device; if they do not fit, the lists
      // come back to the CPU and the forces are made from the arrays of the CPU as before.
      std::vector<std::vector<SphOut> > psrec(ndev);
      std::vector<std::vector<int> > pidx(ndev), plist(ndev);
      if (useb && g_bs.nflag > 0)
      {
        bool fits = true;
        for (int k = 0; k < ndev && fits; k++)
        {
          Dev &d = g_dev[k];
          const size_t pb = (size_t)slice_of(nt, ndev, k) * CAP;
          for (size_t a = 0; a < g_bs.frec[k].size(); a++)
          {
            const int o = g_bs.frec[k][a];
            const int i = d.hout[o].idx;
            SphOut r = d.hout[o];
            r.hp = hp[i]; r.rho = rho[i]; r.omega = om[i]; r.om0 = om0[i]; r.psi = psi[i]; r.wn = 0.0; r.nn = nn[i]; r.flag = 0;
            r.first = (int)(pb + plist[k].size() - (size_t)(o / g_bs.ch) * g_bs.ch * CAP);
            for (int e = 0; e < nn[i]; e++) plist[k].push_back(list[(size_t)first[i] + e]);
            psrec[k].push_back(r);
            pidx[k].push_back(o);
          }
          if (pb + plist[k].size() > ((size_t)d.nowncap + d.chcap) * CAP || pb + plist[k].size() > 2000000000ULL) fits = false;
        }
        if (!fits)
        {
          int ie = 0;
          getlists(list, &ie);
          if (ie != 0) throw SphError(ie, "the lists could not be brought back to the CPU");
          useb = false;
        }
      }
      // the longest list reached (the part of list() to send), and the work of each device: whole chunks of pb_chunk places (or, with the
      // data of the smoothing-length kernel, its slice of the particles)
      // three or more devices with the data of the smoothing-length kernel: a device brings back only what it touched.  With one device that is everything, with two it is more than a
      // half of the particles (the counters: 0.50-0.56 of n for each device in the collision and 0.52 in n500, 0.27-0.31 with four, 0.15 with eight) and the compaction does not pay
      // (measured: 2 devices, n500: 1% slower): they use the plain copy and the plain sum.  Both ways give the same bits.
      const bool cmp = useb && ndev > 2;
      long long listlen = 0, ebound = 0;
      std::vector<int> t0(ndev + 1, nt), nts(ndev, 0);
      std::vector<long long> ent(ndev, 0);
      if (!useb)
      {
#pragma omp parallel for schedule(static) reduction(max : listlen)
        for (int t = 0; t < nt; t++)
        {
          const int p = (contig ? p0 + t : tp[t]) - 1;
          const long long v = (long long)first[p] + nn[p];
          if (v > listlen) listlen = v;
        }
        const int nchk = (nloop + PBCHUNK - 1) / PBCHUNK;
        for (int k = 0; k < ndev; k++)
        {
          const long long cb = ((long long)nchk * k) / ndev;
          const long long iw0 = cb * PBCHUNK + 1;   // the first place of the work that the device takes
          if (k == 0) t0[k] = 0;
          else if (contig) t0[k] = (int)std::min<long long>(nt, iw0 - 1);
          else t0[k] = (int)(std::lower_bound(tiw, tiw + nt, (int)std::min<long long>(iw0, 2000000000LL)) - tiw);
        }
        for (int k = 0; k < ndev; k++)
        {
          long long e = 0;
#pragma omp parallel for schedule(static) reduction(+ : e)
          for (int t = t0[k]; t < t0[k + 1]; t++) e += nn[(contig ? p0 + t : tp[t]) - 1];
          ent[k] = e;
          nts[k] = t0[k + 1] - t0[k];
          ebound = std::max(ebound, e);
        }
      }
      else
      {
        long long e = 0;
#pragma omp parallel for schedule(static) reduction(+ : e)
        for (int t = 0; t < nt; t++) e += nn[p0 - 1 + t];
        ebound = e;      // (what a device does is a part of it)
        for (int k = 0; k < ndev; k++)
        {
          t0[k] = 0;
          nts[k] = slice_of(nt, ndev, k) > 0 ? nt : 0;
        }
      }
      if (ebound >= 2000000000LL) throw SphError(4, "too many entries in the lists for the forces on the GPU");
      // the velocities and u were sent and packed while the devices solved, if the devices still hold them (and no buffer has to be made larger: that frees them)
      bool fe = useb && ndev == 1 && g_fe.valid && g_fe.n == n;
      if (fe)
      {
        const FDev &f0 = g_fdev[0];
        if (n > f0.ncap || nt > f0.ntcap || (size_t)std::max(1LL, ebound) > f0.ecap || (size_t)listlen > f0.lcap || !f0.x) fe = false;
      }
      if (fe && fe_check())
      {
        // (a test) the arrays of the CPU are what the devices were sent
        const double *src[4] = {vx, vy, vz, u};
        for (int a = 0; a < 4; a++)
          if (memcmp(g_estage + (size_t)a * n * 8, src[a], (size_t)n * 8) != 0)
          {
            fprintf(stderr, "StarSmasher (SPH on the GPU): the velocities or u changed between the smoothing lengths and the forces (SPHGPU_FECHECK)\n");
            abort();
          }
        static int nchecked = 0;
        if (++nchecked == 1 || nchecked % 100 == 0) fprintf(stderr, "SPHGPU_FECHECK: the velocities and u sent early were the ones the CPU holds (%d calls)\n", nchecked);
      }
      g_fe.valid = false;
      // what is sent: the velocities, u and the pressure term always (unless they were sent early); the rest unless the devices hold it
      Stager stg;
      const size_t o_vx = fe ? 0 : stg.put(vx, (size_t)n * 8), o_vy = fe ? 0 : stg.put(vy, (size_t)n * 8), o_vz = fe ? 0 : stg.put(vz, (size_t)n * 8), o_u = fe ? 0 : stg.put(u, (size_t)n * 8),
                   o_por2 = stg.put(por2, (size_t)n * 8);
      size_t o_x = 0, o_y = 0, o_z = 0, o_am = 0, o_hp = 0, o_rho = 0, o_om = 0, o_om0 = 0, o_psi = 0, o_nn = 0, o_first = 0, o_list = 0, o_act = 0, o_tp = 0, o_tiw = 0;
      if (!useb)
      {
        o_x = stg.put(x, (size_t)n * 8); o_y = stg.put(y, (size_t)n * 8); o_z = stg.put(z, (size_t)n * 8); o_am = stg.put(am, (size_t)n * 8);
        o_hp = stg.put(hp, (size_t)n * 8); o_rho = stg.put(rho, (size_t)n * 8); o_om = stg.put(om, (size_t)n * 8); o_om0 = stg.put(om0, (size_t)n * 8);
        o_psi = stg.put(psi, (size_t)n * 8); o_nn = stg.put(nn, (size_t)n * 4); o_first = stg.put(first, (size_t)n * 4);
        o_list = stg.put(list, (size_t)listlen * 4);
        if (blocks) o_act = stg.put(act, (size_t)n * 4);
        if (!contig) { o_tp = stg.put(tp, (size_t)nt * 4); o_tiw = stg.put(tiw, (size_t)nt * 4); }
      }
      std::vector<size_t> o_prec(ndev, 0), o_pidx(ndev, 0), o_plist(ndev, 0);
      if (useb)
        for (int k = 0; k < ndev; k++)
        {
          if (psrec[k].empty()) continue;
          o_prec[k] = stg.put(psrec[k].data(), psrec[k].size() * sizeof(SphOut));
          o_pidx[k] = stg.put(pidx[k].data(), pidx[k].size() * sizeof(int));
          if (!plist[k].empty()) o_plist[k] = stg.put(plist[k].data(), plist[k].size() * sizeof(int));
        }
      stg.fill();
      FPar P;
      P.alpha1 = *alpha1; P.beta1 = *beta1; P.hfloor = *hfloor; P.gam = *gam;
      P.nav3 = (*nav == 3); P.nintvar2 = (*nintvar == 2); P.relaxnoheat = *relaxnoheat; P.udotother = *udotother; P.ngr = (*ngr != 0);
      P.blocks = blocks;
      const int kb = key_bits(n);
      std::vector<FView> views(ndev);
      std::vector<int> used(ndev, 0);
      // phase A on every device: the buffers, the data, the places of the work, the entries of the lists
      for (int k = 0; k < ndev; k++)
      {
        if (nts[k] <= 0) continue;
        used[k] = 1;
        Dev &d = g_dev[k];
        FDev &f = g_fdev[k];
        CK(cudaSetDevice(d.id));
        const size_t E = (size_t)std::max(1LL, ebound);
        const bool mk = (k == 0);
        if (n > f.ncap || nt > f.ntcap || E > f.ecap || (size_t)listlen > f.lcap || !f.x)
        {
          cudaStreamSynchronize(d.st);
          ffree(f);
          CK(cudaSetDevice(d.id));
          const size_t ncap = (size_t)n + n / 8 + 1024, ntcap = (size_t)nt + nt / 8 + 1024, ecap = E + E / 8 + 1024,
                       lcap = (size_t)listlen + listlen / 8 + 1024;
          double **dd[] = {&f.x, &f.y, &f.z, &f.vx, &f.vy, &f.vz, &f.am, &f.u, &f.hp, &f.rho, &f.por2, &f.om, &f.om0, &f.psi};
          for (size_t q = 0; q < sizeof dd / sizeof dd[0]; q++) dmalloc(*dd[q], ncap);
          dmalloc(f.acc, 6 * ncap);
          dmalloc(f.pos, ncap); dmalloc(f.vel, ncap);
          dmalloc(f.nn, ncap); dmalloc(f.first, ncap); dmalloc(f.act, ncap);
          dmalloc(f.tp, ntcap); dmalloc(f.tiw, ntcap); dmalloc(f.eoff, ntcap + 2); dmalloc(f.cnt, ntcap + 2); dmalloc(f.out, ntcap);
          dmalloc(f.tpar, ntcap);
          dmalloc(f.inv, ntcap);
          dmalloc(f.list, lcap);
          dmalloc(f.c, ecap);
          dmalloc(f.ka, ecap); dmalloc(f.kb, ecap); dmalloc(f.va, ecap); dmalloc(f.vb, ecap); dmalloc(f.et, ecap);
          dmalloc(f.nlarge, 4);
          dmalloc(f.tmask, ncap); dmalloc(f.tpos, ncap); dmalloc(f.tidx, ncap); dmalloc(f.tcnt, 4);
          dmalloc(f.tv, 6 * ncap); dmalloc(f.oc, ntcap);
          size_t b1 = 0, b2 = 0, b3 = 0;
          CK(cub::DeviceRadixSort::SortPairs(NULL, b1, (int *)0, (int *)0, (int *)0, (int *)0, (int)ecap, 0, 32));
          CK(cub::DeviceScan::ExclusiveSum(NULL, b2, (int *)0, (int *)0, (int)ntcap + 1));
          CK(cub::DeviceScan::ExclusiveSum(NULL, b3, (int *)0, (int *)0, (int)ncap));
          b2 = std::max(b2, b3);
          f.tmpbytes = std::max(b1, b2) + 256;
          CK(cudaMalloc(&f.tmp, f.tmpbytes));
          CK(cudaMallocHost((void **)&f.hacc, 6 * ncap * sizeof(double)));
          CK(cudaMallocHost((void **)&f.hout, ntcap * sizeof(FOut)));
          CK(cudaMallocHost((void **)&f.hnlarge, 4 * sizeof(int)));
          CK(cudaMallocHost((void **)&f.hE, 4 * sizeof(int)));
          CK(cudaMallocHost((void **)&f.hidx, ncap * sizeof(int)));
          CK(cudaMallocHost((void **)&f.hcnt, 4 * sizeof(int)));
          f.haccn = ncap; f.houtn = ntcap;
          f.ncap = (int)ncap; f.ntcap = (int)ntcap; f.ecap = ecap; f.lcap = lcap;
        }
        cudaStream_t st = d.st;
        if (mk) mark("f start", st);
        CK(cudaMemsetAsync(f.acc, 0, 6 * (size_t)n * sizeof(double), st));
        CK(cudaMemsetAsync(f.nlarge, 0, sizeof(int), st));
        if (cmp) CK(cudaMemsetAsync(f.tmask, 0, (size_t)n * sizeof(int), st));
        if (useb && !psrec[k].empty())
        {
          // the particles the CPU solved: their records and lists into the device's
          const size_t np = psrec[k].size();
          if (np > f.pcap)
          {
            cudaStreamSynchronize(st);
            dfree(f.prec); dfree(f.pidx);
            f.pcap = np + np / 2 + 64;
            dmalloc(f.prec, f.pcap); dmalloc(f.pidx, f.pcap);
          }
          H2S(f.prec, o_prec[k], np * sizeof(SphOut));
          H2S(f.pidx, o_pidx[k], np * sizeof(int));
          if (!plist[k].empty()) h2d_staged(d.list + (size_t)slice_of(nt, ndev, k) * CAP, o_plist[k], plist[k].size() * sizeof(int), st, "the lists of the particles solved by the CPU");
          k_patch<<<nblocks((int)np, 128), 128, 0, st>>>((int)np, f.pidx, f.prec, d.out);
        }
        if (!fe)
        {
          H2S(f.vx, o_vx, (size_t)n * 8); H2S(f.vy, o_vy, (size_t)n * 8); H2S(f.vz, o_vz, (size_t)n * 8);
          H2S(f.u, o_u, (size_t)n * 8);
        }
        H2S(f.por2, o_por2, (size_t)n * 8);
        const double *dx_ = f.x, *dy_ = f.y, *dz_ = f.z, *dam_ = f.am;
        if (!useb)
        {
          H2S(f.x, o_x, (size_t)n * 8); H2S(f.y, o_y, (size_t)n * 8); H2S(f.z, o_z, (size_t)n * 8);
          H2S(f.am, o_am, (size_t)n * 8);
          H2S(f.hp, o_hp, (size_t)n * 8); H2S(f.rho, o_rho, (size_t)n * 8);
          H2S(f.om, o_om, (size_t)n * 8); H2S(f.om0, o_om0, (size_t)n * 8); H2S(f.psi, o_psi, (size_t)n * 8);
          H2S(f.nn, o_nn, (size_t)n * 4); H2S(f.first, o_first, (size_t)n * 4);
          H2S(f.list, o_list, (size_t)listlen * 4);
          if (blocks) H2S(f.act, o_act, (size_t)n * 4);
          if (contig) k_fcontig<<<nblocks(nt, 256), 256, 0, st>>>(nt, p0, f.tp, f.tiw);
          else { H2S(f.tp, o_tp, (size_t)nt * 4); H2S(f.tiw, o_tiw, (size_t)nt * 4); }
        }
        else
        {
          dx_ = d.x; dy_ = d.y; dz_ = d.z; dam_ = d.am;
          k_fcontig<<<nblocks(nt, 256), 256, 0, st>>>(nt, p0, f.tp, f.tiw);
        }
        if (mk) mark("f h2d", st);
        if (!fe) k_fpack<<<nblocks(n, 256), 256, 0, st>>>(n, dx_, dy_, dz_, dam_, f.vx, f.vy, f.vz, f.u, f.pos, f.vel);
        if (mk) mark("f pack", st);
        FView &V = views[k];
        memset(&V, 0, sizeof V);
        V.n = n; V.nts = nts[k]; V.t0 = t0[k];
        V.pos = f.pos; V.vel = f.vel; V.por2 = f.por2; V.tpar = f.tpar;
        V.list = useb ? d.list : f.list; V.act = f.act; V.tp = f.tp; V.tiw = f.tiw; V.eoff = f.eoff;
        V.dw = f.dw; V.dg = f.dg; V.dwf = f.dwf; V.dgf = f.dgf;
        V.P = P; V.P.ctab = d.ctab; V.P.c4 = 4.0 * d.ctab;
        V.c = f.c; V.key = f.ka; V.et = f.et; V.val = f.va; V.skey = f.kb; V.sval = f.vb;
        V.acc = f.acc; V.out = f.out; V.nlarge = f.nlarge;
        if (!useb)
        {
          k_ftpar_u<<<nblocks(nts[k], 256), 256, 0, st>>>(nts[k], t0[k], f.tp, f.hp, f.rho, f.om, f.om0, f.psi, f.nn, f.first, f.tpar);
        }
        else
        {
          const int sl = slice_of(nt, ndev, k);
#ifdef SPH_FTPAR_SCATTER
          CK(cudaMemsetAsync(f.tpar, 0, (size_t)nt * sizeof(TPar), st));
          k_ftpar_r<<<nblocks(sl, 256), 256, 0, st>>>(sl, g_bs.nlow, g_bs.ch, d.out, f.tpar);
#else
          CK(cudaMemsetAsync(f.inv, 0xff, (size_t)nt * sizeof(int), st));
          k_finv<<<nblocks(sl, 256), 256, 0, st>>>(sl, g_bs.nlow, d.out, f.inv);
          V.inv = f.inv; V.sout = d.out; V.sch = g_bs.ch;
#endif
        }
        if (mk) mark("f places", st);
        k_fcount<<<nblocks(nts[k] + 1, 256), 256, 0, st>>>(V, f.cnt);
        size_t bytes = f.tmpbytes;
        CK(cub::DeviceScan::ExclusiveSum(f.tmp, bytes, f.cnt, f.eoff, nts[k] + 1, st));
        CK(cudaMemcpyAsync(f.hE, f.eoff + nts[k], sizeof(int), cudaMemcpyDeviceToHost, st));
        if (mk) mark("f count+scan", st);
      }
      // phase B: the pair terms, the sort, the sums on the neighbours, the results back
      for (int k = 0; k < ndev; k++)
      {
        if (!used[k]) continue;
        Dev &d = g_dev[k];
        FDev &f = g_fdev[k];
        CK(cudaSetDevice(d.id));
        cudaStream_t st = d.st;
        const bool mk = (k == 0);
        FView &V = views[k];
        CK(cudaStreamSynchronize(st));
        ent[k] = f.hE[0];
        if ((size_t)ent[k] > f.ecap) throw SphError(4, "the lists of the forces do not fit in the buffers");
        V.nent = (int)ent[k];
        if (prec) k_fpair_mix<<<nblocks(nts[k], FTHR), FTHR, 0, st>>>(V);
        else k_fpair<<<nblocks(nts[k], FTHR), FTHR, 0, st>>>(V);
        CK(cudaGetLastError());
        if (mk) mark("f pairs", st);
        if (cmp)
        {
          // the particles this device touches: its targets and the keys of its dispersion (everything else in its sums is +0.0)
          const int sl = slice_of(nt, ndev, k);
          if (sl > 0) k_tmark_own<<<nblocks(sl, 128), 128, 0, st>>>(sl, d.out, f.tmask);
          if (ent[k] > 0) k_tmark_keys<<<nblocks((int)ent[k], 256), 256, 0, st>>>((int)ent[k], n, f.ka, f.tmask);
          size_t tb = f.tmpbytes;
          CK(cub::DeviceScan::ExclusiveSum(f.tmp, tb, f.tmask, f.tpos, n, st));
          k_tcompact<<<nblocks(n, 256), 256, 0, st>>>(n, f.tmask, f.tpos, f.tidx, f.tcnt);
          CK(cudaMemcpyAsync(f.hcnt, f.tcnt, sizeof(int), cudaMemcpyDeviceToHost, st));
        }
        if (ent[k] > 0)
        {
          size_t bytes = f.tmpbytes;
          CK(cub::DeviceRadixSort::SortPairs(f.tmp, bytes, f.ka, f.kb, f.va, f.vb, (int)ent[k], 0, kb, st));
          if (mk) mark("f sort", st);
          k_fscat<<<nblocks(n, FTHR), FTHR, 0, st>>>(V);
          CK(cudaGetLastError());
          if (mk) mark("f scatter", st);
        }
        if (!cmp)
        {
          CK(cudaMemcpyAsync(f.hacc, f.acc, 6 * (size_t)n * sizeof(double), cudaMemcpyDeviceToHost, st));
          CK(cudaMemcpyAsync(f.hout, f.out, (size_t)nts[k] * sizeof(FOut), cudaMemcpyDeviceToHost, st));
          CK(cudaMemcpyAsync(f.hnlarge, f.nlarge, sizeof(int), cudaMemcpyDeviceToHost, st));
          if (mk) mark("f d2h", st);
        }
      }
      // The work the caller has for its CPUs while the devices are busy with phase B (work: a function or NULL).  It must read what the caller keeps
      // in its arrays and write nothing that this call writes (mvx.. udot, uijmax, divv) or reads from the CPU's arrays after this point.
      if (work && *work) (*work)();
      // with the data of the smoothing-length kernel only what a device touched comes back: the sums of those particles and the results of its targets
      if (cmp)
      {
        for (int k = 0; k < ndev; k++)
        {
          if (!used[k]) continue;
          Dev &d = g_dev[k];
          FDev &f = g_fdev[k];
          CK(cudaSetDevice(d.id));
          cudaStream_t st = d.st;
          CK(cudaStreamSynchronize(st));
          const int m = f.hcnt[0], sl = slice_of(nt, ndev, k);
          if (m > 0) k_tgather<<<nblocks(m, 128), 128, 0, st>>>(m, n, f.tidx, f.acc, f.tv);
          if (sl > 0) k_foutc<<<nblocks(sl, 128), 128, 0, st>>>(sl, g_bs.nlow, d.out, f.out, f.oc);
          CK(cudaMemcpyAsync(f.hacc, f.tv, (size_t)m * 6 * sizeof(double), cudaMemcpyDeviceToHost, st));
          CK(cudaMemcpyAsync(f.hidx, f.tidx, (size_t)m * sizeof(int), cudaMemcpyDeviceToHost, st));
          CK(cudaMemcpyAsync(f.hout, f.oc, (size_t)sl * sizeof(FOut), cudaMemcpyDeviceToHost, st));
          CK(cudaMemcpyAsync(f.hnlarge, f.nlarge, sizeof(int), cudaMemcpyDeviceToHost, st));
          if (k == 0) mark("f d2h", st);
        }
      }
      for (int k = 0; k < ndev; k++)
      {
        if (!used[k]) continue;
        CK(cudaSetDevice(g_dev[k].id));
        CK(cudaStreamSynchronize(g_dev[k].st));
      }
      const double tdev = now_ms();
      // nothing of the CPU has been changed up to here: the results now (the devices in order, so the sums do not depend on timing)
      double *outs[6] = {mvx, mvy, mvz, mgx, mgy, mgz};
      int first_used = -1;
      for (int k = 0; k < ndev; k++)
        if (used[k]) { first_used = k; break; }
      if (!cmp)
      {
        for (int a = 0; a < 6; a++)
        {
          double *dst = outs[a];
#pragma omp parallel for schedule(static)
          for (int q = 0; q < n; q++)
          {
            double s = g_fdev[first_used].hacc[(size_t)a * n + q];
            for (int k = first_used + 1; k < ndev; k++)
              if (used[k]) s = s + g_fdev[k].hacc[(size_t)a * n + q];
            dst[q] = s;
          }
        }
      }
      else
      {
        // With the data of the smoothing-length kernel a device holds +0.0 (the clearing) everywhere but at the particles it touched, so the sum of a particle is
        // that of the devices that touched it, added in the order of the devices, from the first one's value.  (The devices that did not touch it would add +0.0,
        // which changes a sum only if it is -0.0; no device's value is -0.0, see k_fpair and k_fscat: every sum there starts from +0.0 and a sum that comes out
        // as zero in round-to-nearest is +0.0, unless both terms are -0.0.)  The particles of one device are distinct, so the threads do not meet.
        for (int a = 0; a < 6; a++)
        {
          double *dst = outs[a];
#pragma omp parallel for schedule(static)
          for (int q = 0; q < n; q++) dst[q] = 0.0;
        }
        for (int k = first_used; k < ndev; k++)
        {
          if (!used[k]) continue;
          const FDev &f = g_fdev[k];
          const int m = f.hcnt[0];
          const double *v = f.hacc;      // six numbers for each touched particle, in the order of the list
          if (k == first_used)
          {
#pragma omp parallel for schedule(static)
            for (int j = 0; j < m; j++)
              for (int a = 0; a < 6; a++) outs[a][f.hidx[j]] = v[(size_t)j * 6 + a];
          }
          else
          {
#pragma omp parallel for schedule(static)
            for (int j = 0; j < m; j++)
              for (int a = 0; a < 6; a++) outs[a][f.hidx[j]] = outs[a][f.hidx[j]] + v[(size_t)j * 6 + a];
          }
        }
      }
      for (int k = 0; k < ndev; k++)
      {
        if (!used[k]) continue;
        const FDev &f = g_fdev[k];
        *nlarge += *f.hnlarge;
        const int nres = cmp ? slice_of(nt, ndev, k) : nts[k];   // (with the data of the smoothing-length kernel only the targets of the device, in the order of its records)
#pragma omp parallel for schedule(static)
        for (int o = 0; o < nres; o++)
        {
          const FOut s = f.hout[o];
          if (s.p < 0) continue;
          udot[s.p] = s.udot;
          divv[s.p] = s.divv;
          if (!blocks || act[s.p]) uijmax[s.p] = s.uij;
        }
      }
      CK(cudaSetDevice(g_dev[0].id));
      if (plev)
      {
        fprintf(stderr, "SPHGPU forces: %d targets, %lld entries, %d devices%s: devices %.2f ms, host %.2f ms, all %.2f ms\n", nt, ent[first_used], ndev,
                useb ? " (data of the smoothing-length kernel)" : "", tdev - tstart, now_ms() - tdev, now_ms() - tstart);
        if (plev >= 2) print_marks("forces phases (ms)");
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
