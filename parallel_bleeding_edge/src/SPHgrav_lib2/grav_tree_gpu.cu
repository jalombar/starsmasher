// Gravity with a tree on the GPUs (ngravtree=2).
//
// The tree, its moments and the list of the leaves are made on the CPUs by gravtree.f90 (the same code as ngravtree=1); only the visit
// of the tree is done here.  Physics of ngravtree=1, step by step (see gravtree.f90, routine walk_leaf):
//   - a node S is well separated from a leaf L if  theta^2 r^2 > (size_L + size_S)^2  and  r > size_L + size_S + max(hmax_L, hmax_S),
//     r being the distance of their centres of mass; then the quadrupole of S goes into the Taylor series of the force about the
//     centre of mass of L (m2l_loop);
//   - a leaf S that is not well separated is summed particle by particle with the pair formula of gravpair.f90 (grav_soft), which is
//     Newton's outside the softening of the two particles and the kernel's inside;
//   - any other node is opened.
// The force on a particle of L is the Taylor series at its position plus the sum over the leaves that are near.
//
// How it is done.  One warp for each leaf, one thread for each particle of the leaf.  The visit is cooperative: the warp takes up to 32
// nodes from a stack in shared memory, each lane tests one of them, and ballot/popc put the results in the lists in the order of the
// lanes: the nodes that are well separated (list far), the leaves that are near (list near), and the children of the nodes to open
// (back on the stack).  Nothing in that depends on anything but the tree, so the result of a leaf is the same on any GPU, with any
// number of blocks, and whatever the other leaves are.  The list far is done 32 nodes at a time, one node per lane (compute_m2l), the 20
// coefficients of the Taylor series are summed in the registers of each lane, and added over the lanes at the end of the visit by a
// butterfly of fixed order.  The list near is done one leaf at a time: its particles go to shared memory and every lane sums them for its
// particle (double-single sums).  A leaf with more than 32 particles (only points that coincide can make one) is done in passes of 32.
//
// Precision.  The positions go as pairs (hi, lo) of floats relative to the centre of mass of the system, as with ngravhilo=1 of the direct
// sum: separations a - b = (a.hi - b.hi) + (a.lo - b.lo) have almost double precision.  The same for the centres of mass of the nodes.
// Moments, Taylor series and pair formulas are in float; there is no double in the inner loops (the Blackwell cards have little FP64).
// The pairs of a source leaf are added in float and the sum of the leaf goes to double-single sums; a near leaf that is farther than every
// softening length is summed without the test of each pair (the same bits).  The result of a particle is added in double (4 double per
// particle: ax, ay, az, potential), in the order of the tree.
//
// Leaves with a very long list.  A leaf that contains a particle with a large softening length (the particles thrown out of a collision) cannot
// accept any node near it: its list of near leaves is the whole system, 40 times the average, and one warp doing it alone makes the kernel last as long as that
// leaf.  The first pass (tree_kernel<..., false>) notices it when the list of the near leaves has more than TREE_HEAVY entries at one of the times it is done: it
// does the far field and the entries done until then (the "head") and goes on with the visit only to count the entries.  The second pass (tree_kernel<..., true>)
// has a work item for each chunk of TREE_CHUNK consecutive entries of the list of a heavy leaf: the warp makes the visit of that leaf again (the list and its order
// depend only on the tree), sums only the entries of its chunk and leaves the sums in the pool; the warp that finishes the last chunk adds the sums of the chunks,
// in the order of the chunks, to the head.  A heavy leaf has the same result for any number of blocks, devices and threads.  The pool of the chunks has room for
// n/16 chunks (at least 8192); if it is not enough the call is done again with a larger pool (the result is the same).  The environment variable
// STARSMASHER_TREE_HEAVY replaces TREE_HEAVY, so that a test can make the ordinary models take the path of the heavy leaves.
//
// Sizes: the stack has TREE_STK entries per warp (gpu_tree_wait_ says if it was not enough, and Fortran stops with a message), the lists of far
// nodes and near leaves have room for 256 entries and are done when they have more than 224.
//
// Several devices for one process (ngpuperrank): each takes a run of consecutive leaves from the list of the leaves that have a
// particle to do, divided in proportion to the speed of its first pass measured in the calls before.  The leaves are independent, so the bits do not depend
// on the division.
//
// The interface is in C for the Fortran of gravtree.f90: gpu_tree_reserve_ gives pointers to the (pinned) host buffers that Fortran fills,
// gpu_tree_launch_ queues everything (copies, kernel, copy back) on the devices and returns, gpu_tree_wait_ waits and says if the stack or
// the pool of the chunks had no room.  Nothing is allocated again unless the problem grows.
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cuda_runtime.h>
#include "cutil.h"

// defined in grav_force_direct.cu: the devices of the process when ngpuperrank > 1
namespace SPHgrav
{
  extern int ndev, dev0;
}

#ifndef TREE_WARPS
#define TREE_WARPS 4        // warps (leaves at a time) in a block
#endif
#ifndef TREE_STK
#define TREE_STK 1024       // entries of the stack of a warp
#endif
#ifndef TREE_HEAVY
#define TREE_HEAVY 768      // a leaf whose list of near leaves is longer than this (at one of the times the list is done, see tree_kernel) is divided in chunks
#endif
#ifndef TREE_CHUNK
#define TREE_CHUNK 192      // entries of the list of near leaves in a chunk of a heavy leaf (at most TREE_NEARCAP - 32)
#endif
#define TREE_FARCAP 256     // entries of the list of the nodes that are well separated
#define TREE_NEARCAP 256    // entries of the list of the leaves that are near
#define FULLMASK 0xffffffffu

namespace
{

struct D4
{
  double x, y, z, w;
};

// chunks of heavy leaves found in a call (it goes on counting when the pool is full); err: bit 0 the stack did not have room, bit 1 the
// pool of the chunks was full
struct TreeStats
{
  unsigned long long nchunks;
  int err;
  int pad;
};

struct TreeArgs
{
  const float4 *pa;          // x y z (hi parts) and mass of the particles in the order of the tree
  const float4 *pb;          // x y z (lo parts) and 2h (negative for a point particle)
  const unsigned char *need; // 1 for the particles whose force is wanted
  const float4 *nf;          // 4 float4 for each node: com hi and size; com lo and hmax; mass and qxx qxy qxz; qyy qyz qzz and 0
  const int4 *ni;            // left child, right child (the left one is -1 for a leaf), first and last particle
  const int *work;           // the leaves to do
  int nwork;
  int *counter;
  float theta2;
  D4 *out;
  TreeStats *st;
  int hv, ch;                // heavy leaves: entries of the list near above which a leaf is divided, entries in a chunk
  int poolcap;               // chunks that fit in the pool
  int *hcount;               // [0] heavy leaves found so far, [1] chunks found so far (it goes on counting when the pool is full)
  int4 *hrec;                // for each heavy leaf: leaf, first chunk in the pool, number of chunks, entries done by the first pass
  int4 *pdesc;               // for each chunk of the pool: heavy leaf (its place in hrec), number of the chunk
  int *hdone;                // for each heavy leaf: chunks done
  D4 *part;                  // for each chunk of the pool: 32 lanes x 4 sums
};

// double-single sum, as the direct gravity of this library
struct tds
{
  float x, y;
  __device__ tds() : x(0.0f), y(0.0f) {}
  __device__ __forceinline__ void add(const float v)
  {
    const float vx = x + v;
    y = y - ((vx - x) - v);
    x = vx;
  }
  __device__ __forceinline__ double d() const { return (double)x + (double)y; }
};

// the softening kernels: force and potential factors for q = r/h in [0,1]  (the expressions of gravpair.f90)
template<int NK>
__device__ __forceinline__ void kern_ap(const float q, float &acc, float &pot)
{
  const float q2 = q * q;
  if (NK == 0)
  {
    // the cubic spline
    if (q < 0.5f)
    {
      acc = 10.666666666667f + q2 * (32.0f * q - 38.4f);
      pot = -2.8f + q2 * (5.333333333333f + q2 * (6.4f * q - 9.6f));
    }
    else
    {
      const float iq = 1.0f / q;
      const float q3 = q * q2;
      acc = 21.333333333333f - 48.0f * q + 38.4f * q2 - 10.666666666667f * q3 - 0.066666666667f * (iq * iq * iq);
      pot = -3.2f + 0.066666666667f * iq + q2 * (10.666666666667f + q * (-16.0f + q * (9.6f - 2.133333333333f * q)));
    }
  }
  else if (NK == 1)
  {
    // the Wendland W_{3,3} kernel
    acc = 28.4375f + q2 * (-187.6875f + q2 * (804.375f + q2 * (-4379.375f + q * (9009.0f + q * (-8957.8125f + q * (5005.0f + q * (-1515.9375f + 195.0f * q)))))));
    pot = -3.828125f + q2 * (14.21875f + q2 * (-46.921875f + q2 * (134.0625f + q2 * (-547.421875f + q * (1001.0f + q * (-895.78125f + q * (455.0f + q * (-126.328125f + 15.0f * q))))))));
  }
  else
  {
    // the Wendland C4 kernel
    acc = 20.625f + q2 * (-115.5f + q2 * (618.75f + q * (-1155.0f + q * (962.5f + q * (-396.0f + 65.625f * q)))));
    pot = -3.4375f + q2 * (10.3125f + q2 * (-28.875f + q2 * (103.125f + q * (-165.0f + q * (120.3125f + q * (-44.0f + 6.5625f * q))))));
  }
}

// The pair inside the softening of one of the two particles (grav_soft of gravpair.f90): the acceleration of i from j is gacc*(x_j - x_i)
// and the potential of i gets gpot.  hi, hj: 2h of the two; mj: mass of j; self: j is i itself; flagj: 1 for a gas particle, 0 for a
// point particle.  For a pair at zero distance: no force, and the potential of the particle itself (a point particle does not have it).
template<int NK>
__device__ __forceinline__ void soft_pair(const float r2, const float hi, const float hj, const float mj, const bool self,
                                          const float flagj, float &gacc, float &gpot)
{
  const float invhi = 1.0f / hi, invhj = 1.0f / hj;
  const float mj1i = mj * invhi, mj1j = mj * invhj;
  if (r2 <= 0.0f)
  {
    const float cself = NK == 0 ? 1.4f : (NK == 1 ? 1.9140625f : 1.71875f);
    gacc = 0.0f;
    gpot = -cself * (self ? 1.0f : 0.5f) * flagj * (mj1i + mj1j);
    return;
  }
  const float rinv = rsqrtf(r2);
  const float mrinv1 = rinv * mj;
  const float mrinv3 = rinv * rinv * mrinv1;
  const float r = r2 * rinv;
  const bool inI = r2 <= hi * hi, inJ = r2 <= hj * hj;
  float ai = 0.0f, pi = 0.0f, aj = 0.0f, pj = 0.0f;
  if (inI)
    kern_ap<NK>(r * invhi, ai, pi);
  if (inJ)
    kern_ap<NK>(r * invhj, aj, pj);
  gacc = 0.5f * ((inI ? mj1i * invhi * invhi * ai : mrinv3) + (inJ ? mj1j * invhj * invhj * aj : mrinv3));
  gpot = 0.5f * ((inI ? mj1i * pi : -mrinv1) + (inJ ? mj1j * pj : -mrinv1));
}

// Taylor coefficients of the force at the centre of the leaf from a node with mass q0 and second moments q.., at (dx,dy,dz) from it,
// added to f: f[0:3] the acceleration, [3:9] its derivatives xx xy xz yy yz zz, [9:19] its second derivatives xxx xxy xxz xyy xyz xzz yyy
// yyz yzz zzz, [19] the potential.  The expressions of m2l_loop of gravtree.f90 (after Phantom, kdtree.F90, Dehnen 2002; GPLv3).
__device__ __forceinline__ void m2l_add(float (&f)[20], const float dx, const float dy, const float dz, const float dr1, const float q0,
                                        const float qxx, const float qxy, const float qxz, const float qyy, const float qyz, const float qzz)
{
  const float dr12 = dr1 * dr1;
  const float dx2 = dx * dx, dx3 = dx * dx2, dy2 = dy * dy, dy3 = dy * dy2, dz2 = dz * dz, dz3 = dz * dz2;
  const float g0 = -dr1;
  const float g1 = dr12 * g0;
  const float g2 = -3.0f * dr12 * g1;
  const float g3 = -5.0f * dr12 * g2;
  const float g2dx = g2 * dx, g2dy = g2 * dy, g2dz = g2 * dz;
  const float d31 = 3.0f * g2dx + g3 * dx3;
  const float d32 = g2dy + g3 * dx2 * dy;
  const float d33 = g2dz + g3 * dx2 * dz;
  const float d34 = g2dx + g3 * dy2 * dx;
  const float d35 = g3 * dx * dy * dz;
  const float d36 = g2dx + g3 * dz2 * dx;
  const float d37 = 3.0f * g2dy + g3 * dy3;
  const float d38 = g2dz + g3 * dy2 * dz;
  const float d39 = g2dy + g3 * dz2 * dy;
  const float d310 = 3.0f * g2dz + g3 * dz3;
  const float d21 = g1 + g2 * dx2;
  const float d22 = g2dx * dy;
  const float d23 = g2dx * dz;
  const float d24 = g1 + g2 * dy2;
  const float d25 = g2dy * dz;
  const float d26 = g1 + g2 * dz2;
  const float d11 = g1 * dx, d12 = g1 * dy, d13 = g1 * dz;
  f[0] += d11 * q0 + 0.5f * (d31 * qxx + 2.0f * (d32 * qxy + d33 * qxz + d35 * qyz) + d34 * qyy + d36 * qzz);
  f[1] += d12 * q0 + 0.5f * (d32 * qxx + 2.0f * (d34 * qxy + d35 * qxz + d38 * qyz) + d37 * qyy + d39 * qzz);
  f[2] += d13 * q0 + 0.5f * (d33 * qxx + 2.0f * (d35 * qxy + d36 * qxz + d39 * qyz) + d38 * qyy + d310 * qzz);
  f[3] += d21 * q0;
  f[4] += d22 * q0;
  f[5] += d23 * q0;
  f[6] += d24 * q0;
  f[7] += d25 * q0;
  f[8] += d26 * q0;
  f[9] += d31 * q0;
  f[10] += d32 * q0;
  f[11] += d33 * q0;
  f[12] += d34 * q0;
  f[13] += d35 * q0;
  f[14] += d36 * q0;
  f[15] += d37 * q0;
  f[16] += d38 * q0;
  f[17] += d39 * q0;
  f[18] += d310 * q0;
  f[19] += g0 * q0 - 0.5f * (d21 * qxx + d24 * qyy + d26 * qzz + 2.0f * (d22 * qxy + d23 * qxz + d25 * qyz));
}


// The force and the potential at (dx,dy,dz) from the centre of the leaf from the Taylor series (expand_taylor of gravtree.f90).
__device__ __forceinline__ void expand_taylor(const float (&f)[20], const float dx, const float dy, const float dz, float &fxi, float &fyi,
                                              float &fzi, float &poti)
{
  const float ux = 0.5f * (dx * (dx * f[9] + dy * f[10] + dz * f[11]) + dy * (dx * f[10] + dy * f[12] + dz * f[13]) +
                           dz * (dx * f[11] + dy * f[13] + dz * f[14]));
  const float uy = 0.5f * (dx * (dx * f[10] + dy * f[12] + dz * f[13]) + dy * (dx * f[12] + dy * f[15] + dz * f[16]) +
                           dz * (dx * f[13] + dy * f[16] + dz * f[17]));
  const float uz = 0.5f * (dx * (dx * f[11] + dy * f[13] + dz * f[14]) + dy * (dx * f[13] + dy * f[16] + dz * f[17]) +
                           dz * (dx * f[14] + dy * f[17] + dz * f[18]));
  fxi = f[0] + dx * f[3] + dy * f[4] + dz * f[5] + ux;
  fyi = f[1] + dx * f[4] + dy * f[6] + dz * f[7] + uy;
  fzi = f[2] + dx * f[5] + dy * f[7] + dz * f[8] + uz;
  const float t3 = dx * ux + dy * uy + dz * uz;
  poti = f[19] - (f[0] * dx + f[1] * dy + f[2] * dz) -
         0.5f * (dx * (dx * f[3] + dy * f[4] + dz * f[5]) + dy * (dx * f[4] + dy * f[6] + dz * f[7]) + dz * (dx * f[5] + dy * f[7] + dz * f[8])) -
         t3 / 3.0f;
}

// The kernel.  CHUNK = false: the first pass, one warp for each leaf of the list work.  A leaf whose list of near leaves is very long (more than a.hv entries
// at one of the times the list is done) is not done to the end here: the near leaves that were already done (the first E0 of the list, in the order of the visit)
// are its "head", added to the far field and written to out; the others are left to the second pass, in chunks of a.ch consecutive entries of the list.  CHUNK =
// true: the second pass.  A work item is a chunk of a heavy leaf: the visit of the tree is done again for that leaf (it gives the same list in the same order, it
// does not depend on what is summed), only the entries of the chunk are summed, and the sums go to part.  The warp that finishes the last chunk of a leaf adds the
// sums of the chunks in the order of the chunks to the head, so the result does not depend on who did what, nor on the number of devices.
template<int NK, bool CHUNK>
__global__ void __launch_bounds__(TREE_WARPS * 32, 4) tree_kernel(const TreeArgs a)
{
  __shared__ int s_stk[TREE_WARPS][TREE_STK];
  __shared__ int s_far[TREE_WARPS][TREE_FARCAP];
  __shared__ int s_near[TREE_WARPS][TREE_NEARCAP];
  __shared__ float4 s_ta[TREE_WARPS][32];
  __shared__ float4 s_tb[TREE_WARPS][32];
  const int lane = threadIdx.x & 31;
  const int wid = threadIdx.x >> 5;
  int *const stk = s_stk[wid];
  int *const farl = s_far[wid];
  int *const nearl = s_near[wid];
  float4 *const ta = s_ta[wid];
  float4 *const tb = s_tb[wid];
  const unsigned lt = (1u << lane) - 1u;
  const float th2 = a.theta2;

  for (;;)
  {
    int w = 0;
    if (lane == 0)
      w = atomicAdd(a.counter, 1);
    w = __shfl_sync(FULLMASK, w, 0);
    int leaf, hslot = 0, hbase = 0, hnch = 0, elo = 0, ehi = 0;
    if (CHUNK)
    {
      // the second pass: the chunk w of the pool (the number of chunks is final: the first pass is over)
      if (w >= min(__ldcg(&a.hcount[1]), a.poolcap))
        break;
      const int4 pd = __ldcg(&a.pdesc[w]);
      const int4 hr = __ldcg(&a.hrec[pd.x]);
      leaf = hr.x;
      hslot = pd.x;
      hbase = hr.y;
      hnch = hr.z;
      elo = hr.w + pd.y * a.ch;
      ehi = elo + a.ch;
    }
    else
    {
      if (w >= a.nwork)
        break;
      leaf = __ldg(&a.work[w]);
    }
    const int4 li = __ldg(&a.ni[leaf]);
    const int l0 = li.z, u0 = li.w;
    const float4 La = __ldg(&a.nf[4 * leaf]), Lb = __ldg(&a.nf[4 * leaf + 1]);
    // a leaf of more than 32 particles (only points that coincide can make one) is done in passes of 32 and is never divided in chunks
    const bool onepass = (u0 - l0) < 32;

    // the particles of the leaf in passes of 32 (one pass but for points that coincide)
    for (int c0 = l0; c0 <= u0; c0 += 32)
    {
      const int pidx = c0 + lane;
      const bool has = pidx <= u0;
      const bool nd_i = has && __ldg(&a.need[pidx]) != 0;
      const unsigned needmask = __ballot_sync(FULLMASK, nd_i);
      if (needmask == 0)
        continue;
      float4 qa = make_float4(0.0f, 0.0f, 0.0f, 0.0f), qb = qa;
      if (has)
      {
        qa = __ldg(&a.pa[pidx]);
        qb = __ldg(&a.pb[pidx]);
      }
      const float hi = fabsf(qb.w);
      const float h2i = hi * hi;

      float f[20];
#pragma unroll
      for (int k = 0; k < 20; k++)
        f[k] = 0.0f;
      tds sax, say, saz, spt;
      int sp = 1;
      int nfar = 0, nnear = 0;
      // the first pass: the entries of the list near that were done (nflushed), the ones left out when the leaf is heavy (ndrop), and whether it is heavy.
      // the second pass: the number of entries of the list near so far (ecount)
      int nflushed = 0, ndrop = 0, ecount = 0;
      bool heavy = false;

      // the list far: 32 nodes at a time, a node for each lane.  fin: the last time (the lanes without node do nothing)
      auto flush_far = [&](const bool fin)
      {
        const int ng = fin ? ((nfar + 31) >> 5) : (nfar >> 5);
        for (int g = 0; g < ng; ++g)
        {
          const int idx = (g << 5) + lane;
          if (idx < nfar)
          {
            const int s = farl[idx];
            const float4 sa = __ldg(&a.nf[4 * s]), sb = __ldg(&a.nf[4 * s + 1]);
            const float dx = (La.x - sa.x) + (Lb.x - sb.x);
            const float dy = (La.y - sa.y) + (Lb.y - sb.y);
            const float dz = (La.z - sa.z) + (Lb.z - sb.z);
            const float dr1 = rsqrtf(dx * dx + dy * dy + dz * dz);
            const float4 q0 = __ldg(&a.nf[4 * s + 2]), q1 = __ldg(&a.nf[4 * s + 3]);
            m2l_add(f, dx, dy, dz, dr1, q0.x, q0.y, q0.z, q0.w, q1.x, q1.y, q1.z);
          }
        }
        if (fin)
        {
          nfar = 0;
        }
        else
        {
          const int done = ng << 5;
          const int rem = nfar - done;
          int tmp = 0;
          if (lane < rem)
            tmp = farl[done + lane];
          __syncwarp();
          if (lane < rem)
            farl[lane] = tmp;
          __syncwarp();
          nfar = rem;
        }
      };

      // the list near: a leaf at a time, its particles to shared memory, every lane sums them for its particle
      auto flush_near = [&]()
      {
        for (int e = 0; e < nnear; ++e)
        {
          const int ent = nearl[e];
          const int s = ent & 0x3fffffff;
          const bool lean = (ent & 0x40000000) != 0;
          const int4 si = __ldg(&a.ni[s]);
          for (int t0 = si.z; t0 <= si.w; t0 += 32)
          {
            const int cn = min(32, si.w - t0 + 1);
            if (lane < cn)
            {
              ta[lane] = __ldg(&a.pa[t0 + lane]);
              tb[lane] = __ldg(&a.pb[t0 + lane]);
            }
            __syncwarp();
            if (nd_i)
            {
              float px = 0.0f, py = 0.0f, pz = 0.0f, pp = 0.0f;
#define TACC_X(v) px += (v)
#define TACC_Y(v) py += (v)
#define TACC_Z(v) pz += (v)
#define TACC_P(v) pp += (v)
              if (lean)
              {
#pragma unroll 4
                for (int k = 0; k < cn; ++k)
                {
                  const float4 ja = ta[k], jb = tb[k];
                  const float dx = (ja.x - qa.x) + (jb.x - qb.x);
                  const float dy = (ja.y - qa.y) + (jb.y - qb.y);
                  const float dz = (ja.z - qa.z) + (jb.z - qb.z);
                  const float r2 = dx * dx + dy * dy + dz * dz;
                  const float rinv = rsqrtf(r2);
                  const float mr1 = rinv * ja.w;
                  const float mr3 = rinv * rinv * mr1;
                  TACC_X(mr3 * dx);
                  TACC_Y(mr3 * dy);
                  TACC_Z(mr3 * dz);
                  TACC_P(-mr1);
                }
              }
              else
              {
#pragma unroll 4
              for (int k = 0; k < cn; ++k)
              {
                const float4 ja = ta[k], jb = tb[k];
                const float dx = (ja.x - qa.x) + (jb.x - qb.x);
                const float dy = (ja.y - qa.y) + (jb.y - qb.y);
                const float dz = (ja.z - qa.z) + (jb.z - qb.z);
                const float r2 = dx * dx + dy * dy + dz * dz;
                const float hj = fabsf(jb.w);
                const float h2j = hj * hj;
                if (r2 >= fmaxf(h2i, h2j))
                {
                  const float rinv = rsqrtf(r2);
                  const float mr1 = rinv * ja.w;
                  const float mr3 = rinv * rinv * mr1;
                  TACC_X(mr3 * dx);
                  TACC_Y(mr3 * dy);
                  TACC_Z(mr3 * dz);
                  TACC_P(-mr1);
                }
                else
                {
                  float gacc, gpot;
                  soft_pair<NK>(r2, hi, hj, ja.w, t0 + k == pidx, jb.w >= 0.0f ? 1.0f : 0.0f, gacc, gpot);
                  TACC_X(gacc * dx);
                  TACC_Y(gacc * dy);
                  TACC_Z(gacc * dz);
                  TACC_P(gpot);
                }
              }
              }
              sax.add(px);
              say.add(py);
              saz.add(pz);
              spt.add(pp);
#undef TACC_X
#undef TACC_Y
#undef TACC_Z
#undef TACC_P
            }
            __syncwarp();
          }
        }
        nnear = 0;
      };

      // the visit
      if (lane == 0)
        stk[0] = 0;
      __syncwarp();
      while (sp > 0)
      {
        const int np = min(sp, 32);
        const bool act = lane < np;
        int s = 0;
        if (act)
          s = stk[sp - 1 - lane];
        __syncwarp();
        sp -= np;
        bool isfar = false, isleaf = false, isopen = false, safe = false;
        int lcn = 0, rcn = 0;
        if (act)
        {
          const float4 sa = __ldg(&a.nf[4 * s]), sb = __ldg(&a.nf[4 * s + 1]);
          const int4 si = __ldg(&a.ni[s]);
          const float dx = (La.x - sa.x) + (Lb.x - sb.x);
          const float dy = (La.y - sa.y) + (Lb.y - sb.y);
          const float dz = (La.z - sa.z) + (Lb.z - sb.z);
          const float r2 = dx * dx + dy * dy + dz * dz;
          const float sz = La.w + sa.w;
          const float rc = sz + fmaxf(Lb.w, sb.w);
          const bool ok = (th2 * r2 > sz * sz) && (r2 > rc * rc);
          isfar = ok;
          isleaf = !ok && si.x < 0;
          isopen = !ok && si.x >= 0;
          // a leaf that fails the angle but is farther than every softening length: all its pairs with the particles of this leaf are Newtonian
          safe = r2 > rc * rc * 1.00002f;
          lcn = si.x;
          rcn = si.y;
        }
        const unsigned mf = __ballot_sync(FULLMASK, isfar);
        const unsigned mn = __ballot_sync(FULLMASK, isleaf);
        const unsigned mo = __ballot_sync(FULLMASK, isopen);
        const int nopen = __popc(mo);
        const bool ovf = (sp + 2 * nopen > TREE_STK);
        if (ovf && lane == 0)
          atomicOr(&a.st->err, 1);
        if (CHUNK)
        {
          // only the entries ecount .. of the chunk elo .. ehi-1 of the list near are kept (they are in the order of the lanes)
          const int e = ecount + __popc(mn & lt);
          const bool keep = isleaf && e >= elo && e < ehi;
          const unsigned mk = __ballot_sync(FULLMASK, keep);
          if (keep)
            nearl[nnear + __popc(mk & lt)] = safe ? (s | 0x40000000) : s;
          nnear += __popc(mk);
          ecount += __popc(mn);
        }
        else
        {
          if (isfar)
            farl[nfar + __popc(mf & lt)] = s;
          if (isleaf)
            nearl[nnear + __popc(mn & lt)] = safe ? (s | 0x40000000) : s;
          nfar += __popc(mf);
          nnear += __popc(mn);
        }
        if (isopen && !ovf)
        {
          // the right child below the left one: the left one is taken first
          const int p = sp + 2 * __popc(mo & lt);
          stk[p] = rcn;
          stk[p + 1] = lcn;
        }
        sp = ovf ? 0 : sp + 2 * nopen;
        __syncwarp();
        if (CHUNK)
        {
          // the chunk is complete: the rest of the visit would only give entries after it
          if (ecount >= ehi)
            sp = 0;
        }
        else
        {
          if (nfar > TREE_FARCAP - 32)
            flush_far(false);
          if (nnear > TREE_NEARCAP - 32)
          {
            // a leaf with a very long list of near leaves is not done here: from the first time the list is that long on, the entries are only counted
            if (onepass && !heavy && nflushed + nnear > a.hv)
              heavy = true;
            if (heavy)
            {
              ndrop += nnear;
              nnear = 0;
            }
            else
            {
              nflushed += nnear;
              flush_near();
            }
          }
        }
      }
      if (CHUNK)
      {
        flush_near();
        // the sums of this chunk, and, for the warp that finishes the last chunk of the leaf, the total in the order of the chunks
        if (nd_i)
        {
          D4 pr;
          pr.x = sax.d();
          pr.y = say.d();
          pr.z = saz.d();
          pr.w = spt.d();
          a.part[(size_t)w * 32 + lane] = pr;
        }
        __threadfence();
        __syncwarp();
        int last = 0;
        if (lane == 0)
          last = (atomicAdd(&a.hdone[hslot], 1) == hnch - 1) ? 1 : 0;
        last = __shfl_sync(FULLMASK, last, 0);
        if (last)
        {
          __threadfence();
          if (nd_i)
          {
            double tx = 0.0, ty = 0.0, tz = 0.0, tp = 0.0;
            for (int c = 0; c < hnch; ++c)
            {
              const double *pp = reinterpret_cast<const double *>(a.part + ((size_t)(hbase + c) * 32 + lane));
              tx += __ldcg(pp);
              ty += __ldcg(pp + 1);
              tz += __ldcg(pp + 2);
              tp += __ldcg(pp + 3);
            }
            const double *ph = reinterpret_cast<const double *>(a.out + pidx);
            D4 o;
            o.x = __ldcg(ph) + tx;
            o.y = __ldcg(ph + 1) + ty;
            o.z = __ldcg(ph + 2) + tz;
            o.w = __ldcg(ph + 3) + tp;
            a.out[pidx] = o;
          }
        }
        continue;
      }
      flush_far(true);
      if (!heavy)
        flush_near();

      // the coefficients of the Taylor series added over the lanes, in an order that is always the same
#pragma unroll
      for (int k = 0; k < 20; k++)
      {
        float v = f[k];
#pragma unroll
        for (int o = 16; o > 0; o >>= 1)
          v += __shfl_xor_sync(FULLMASK, v, o);
        f[k] = v;
      }
      if (nd_i)
      {
        const float dx = (qa.x - La.x) + (qb.x - Lb.x);
        const float dy = (qa.y - La.y) + (qb.y - Lb.y);
        const float dz = (qa.z - La.z) + (qb.z - Lb.z);
        float tx, ty, tz, tp;
        expand_taylor(f, dx, dy, dz, tx, ty, tz, tp);
        D4 o;
        o.x = sax.d() + (double)tx;
        o.y = say.d() + (double)ty;
        o.z = saz.d() + (double)tz;
        o.w = spt.d() + (double)tp;
        a.out[pidx] = o;
      }
      if (heavy)
      {
        // the chunks of the entries E0 .. E0+nch*ch-1 of the list near, E0 = nflushed, are for the second pass
        const int nch = (nflushed + ndrop + nnear - nflushed + a.ch - 1) / a.ch;
        int slot = 0, base = 0;
        if (lane == 0)
        {
          slot = atomicAdd(&a.hcount[0], 1);
          base = atomicAdd(&a.hcount[1], nch);
        }
        slot = __shfl_sync(FULLMASK, slot, 0);
        base = __shfl_sync(FULLMASK, base, 0);
        if (base + nch <= a.poolcap)
        {
          if (lane == 0)
          {
            a.hrec[slot] = make_int4(leaf, base, nch, nflushed);
            a.hdone[slot] = 0;
          }
          for (int c = lane; c < nch; c += 32)
            a.pdesc[base + c] = make_int4(slot, c, 0, 0);
        }
        else if (lane == 0)
          atomicOr(&a.st->err, 2);
        if (lane == 0)
          atomicAdd(&a.st->nchunks, (unsigned long long)nch);
      }
    }
  }
}

// ---------------------------------------------------------------------------------------------------------------------------
// the host side
// ---------------------------------------------------------------------------------------------------------------------------

// host buffers (pinned, portable: every device of the process copies from and to them)
float4 *h_pa = NULL, *h_pb = NULL, *h_nf = NULL;
unsigned char *h_need = NULL;
int4 *h_ni = NULL;
int *h_work = NULL;
D4 *h_out = NULL;
size_t cap_n = 0, cap_node = 0;

// what a device holds
struct TreeDev
{
  int dev, nsm;
  float4 *pa, *pb, *nf;
  unsigned char *need;
  int4 *ni;
  int *work;
  D4 *out;
  int *counter;
  TreeStats *st;
  size_t cn, cnode, cwork;
  cudaEvent_t e1, e3, em;      // after the copies to the device, after the first pass, after the copy back
  bool events;
  // heavy leaves: the counter of the second pass, the counts, the records of the heavy leaves, the pool of the chunks
  int *counter2, *hcount, *hdone;
  int4 *hrec, *pdesc;
  D4 *part;
  size_t poolcap;
  TreeArgs a1, a2;        // the arguments of the two passes of the last call (to do it again with a larger pool)
  int nk;
  int off, cnt;      // the leaves of this call: work[off .. off+cnt)
  int plo, phi;      // the particles of those leaves: plo .. phi (0-based, in the order of the tree)
  int nlaunch;
  TreeDev() : dev(0), nsm(0), pa(NULL), pb(NULL), nf(NULL), need(NULL), ni(NULL), work(NULL), out(NULL), counter(NULL), st(NULL),
              cn(0), cnode(0), cwork(0), events(false), counter2(NULL), hcount(NULL), hdone(NULL), hrec(NULL), pdesc(NULL), part(NULL), poolcap(0),
              nk(0), off(0), cnt(0), plo(0), phi(-1), nlaunch(0) { memset(&a1, 0, sizeof(a1)); memset(&a2, 0, sizeof(a2)); }
};

// heavy leaves: the threshold above which a leaf is divided in chunks.  The environment variable STARSMASHER_TREE_HEAVY replaces it, so
// that a test can make the ordinary models take the path of the heavy leaves
int t_hv = TREE_HEAVY;
bool t_env_read = false;

void read_env()
{
  if (t_env_read)
    return;
  t_env_read = true;
  const char *e = getenv("STARSMASHER_TREE_HEAVY");
  if (e && atoi(e) > 0) t_hv = atoi(e);
}

std::vector<TreeDev> tdev;
std::vector<double> tspeed;      // leaves per ms of each device, relative to the fastest, smoothed over the calls
int t_first = 0, t_nd = 1;       // the devices of the call in flight
int t_pending = 0;
int t_nwork = 0;
TreeStats *t_host_st = NULL;      // pinned: one for each device, where the counts come back (the copy is asynchronous)

void free_host()
{
  if (h_pa) cudaFreeHost(h_pa);
  if (h_pb) cudaFreeHost(h_pb);
  if (h_nf) cudaFreeHost(h_nf);
  if (h_need) cudaFreeHost(h_need);
  if (h_ni) cudaFreeHost(h_ni);
  if (h_work) cudaFreeHost(h_work);
  if (h_out) cudaFreeHost(h_out);
  h_pa = h_pb = h_nf = NULL;
  h_need = NULL;
  h_ni = NULL;
  h_work = NULL;
  h_out = NULL;
}

void *hostalloc(const size_t bytes)
{
  void *p = NULL;
  CUDA_SAFE_CALL(cudaHostAlloc(&p, std::max<size_t>(bytes, 64), cudaHostAllocPortable));
  assert(p);
  return p;
}

// the pool of the chunks of the heavy leaves of a device: room for cap chunks
void set_pool(TreeDev &d, const size_t cap)
{
  if (d.part) { cudaFree(d.part); cudaFree(d.pdesc); }
  CUDA_SAFE_CALL(cudaMalloc((void **)&d.part, cap * 32 * sizeof(D4)));
  CUDA_SAFE_CALL(cudaMalloc((void **)&d.pdesc, cap * sizeof(int4)));
  d.poolcap = cap;
}

void ensure_dev(TreeDev &d, const size_t n, const size_t nnode, const size_t nwork)
{
  CUDA_SAFE_CALL(cudaSetDevice(d.dev));
  if (!d.events)
  {
    CUDA_SAFE_CALL(cudaEventCreate(&d.e1));
    CUDA_SAFE_CALL(cudaEventCreate(&d.e3));
    CUDA_SAFE_CALL(cudaEventCreate(&d.em));
    CUDA_SAFE_CALL(cudaMalloc((void **)&d.counter, sizeof(int)));
    CUDA_SAFE_CALL(cudaMalloc((void **)&d.counter2, sizeof(int)));
    CUDA_SAFE_CALL(cudaMalloc((void **)&d.hcount, 2 * sizeof(int)));
    CUDA_SAFE_CALL(cudaMalloc((void **)&d.st, sizeof(TreeStats)));
    int nsm = 0;
    CUDA_SAFE_CALL(cudaDeviceGetAttribute(&nsm, cudaDevAttrMultiProcessorCount, d.dev));
    d.nsm = nsm;
    d.events = true;
  }
  if (n > d.cn)
  {
    if (d.pa) { cudaFree(d.pa); cudaFree(d.pb); cudaFree(d.need); cudaFree(d.out); }
    CUDA_SAFE_CALL(cudaMalloc((void **)&d.pa, n * sizeof(float4)));
    CUDA_SAFE_CALL(cudaMalloc((void **)&d.pb, n * sizeof(float4)));
    CUDA_SAFE_CALL(cudaMalloc((void **)&d.need, n));
    CUDA_SAFE_CALL(cudaMalloc((void **)&d.out, n * sizeof(D4)));
    d.cn = n;
  }
  if (nnode > d.cnode)
  {
    if (d.nf) { cudaFree(d.nf); cudaFree(d.ni); }
    CUDA_SAFE_CALL(cudaMalloc((void **)&d.nf, nnode * 4 * sizeof(float4)));
    CUDA_SAFE_CALL(cudaMalloc((void **)&d.ni, nnode * sizeof(int4)));
    d.cnode = nnode;
  }
  if (nwork > d.cwork)
  {
    if (d.work) { cudaFree(d.work); cudaFree(d.hrec); cudaFree(d.hdone); }
    CUDA_SAFE_CALL(cudaMalloc((void **)&d.work, nwork * sizeof(int)));
    CUDA_SAFE_CALL(cudaMalloc((void **)&d.hrec, nwork * sizeof(int4)));
    CUDA_SAFE_CALL(cudaMalloc((void **)&d.hdone, nwork * sizeof(int)));
    d.cwork = nwork;
  }
  const size_t pool = std::max((size_t)8192, n / 16);
  if (pool > d.poolcap)
    set_pool(d, pool);
}

template<int NK, bool CHUNK>
void launch_variant(TreeDev &d, const TreeArgs &args)
{
  int occ = 0;
  CUDA_SAFE_CALL(cudaOccupancyMaxActiveBlocksPerMultiprocessor(&occ, tree_kernel<NK, CHUNK>, TREE_WARPS * 32, 0));
  occ = std::max(occ, 1);
  // the first pass has a block for four leaves of the call (not more than can be resident); the second pass does not know how many chunks there are, and the leaves of a device are few when it
  // holds heavy leaves (see the speed in gpu_tree_wait_), so it has as many blocks as can be resident
  const int blocks = CHUNK ? std::max(1, occ * d.nsm) : std::max(1, std::min(occ * d.nsm, (args.nwork + TREE_WARPS - 1) / TREE_WARPS));
  tree_kernel<NK, CHUNK><<<blocks, TREE_WARPS * 32>>>(args);
}

template<bool CHUNK>
void launch_nk(TreeDev &d, const TreeArgs &args, const int nk)
{
  if (nk == 0)      launch_variant<0, CHUNK>(d, args);
  else if (nk == 1) launch_variant<1, CHUNK>(d, args);
  else              launch_variant<2, CHUNK>(d, args);
}

// the two passes of a call on a device: the leaves, then the chunks of the heavy leaves
void launch_two(TreeDev &d)
{
  CUDA_SAFE_CALL(cudaMemsetAsync(d.counter2, 0, sizeof(int), 0));
  launch_nk<false>(d, d.a1, d.nk);
  CUDA_SAFE_CALL(cudaGetLastError());
  CUDA_SAFE_CALL(cudaEventRecord(d.em, 0));
  launch_nk<true>(d, d.a2, d.nk);
  CUDA_SAFE_CALL(cudaGetLastError());
}

// the division of nw leaves among nd devices in proportion to their speed (largest remainder; one leaf at least for each when there are enough)
void divide_work(const int nw, const int nd, std::vector<int> &cnt)
{
  cnt.assign(nd, 0);
  double sw = 0.0;
  for (int k = 0; k < nd; k++)
    sw += tspeed[k] > 0.0 ? tspeed[k] : 1.0;
  std::vector<double> rem(nd);
  int given = 0;
  for (int k = 0; k < nd; k++)
  {
    const double q = nw * (tspeed[k] > 0.0 ? tspeed[k] : 1.0) / sw;
    cnt[k] = (int)q;
    rem[k] = q - cnt[k];
    given += cnt[k];
  }
  while (given < nw)
  {
    int kb = 0;
    for (int k = 1; k < nd; k++)
      if (rem[k] > rem[kb]) kb = k;
    cnt[kb]++;
    rem[kb] = -1.0;
    given++;
  }
  if (nw >= nd)
    for (int k = 0; k < nd; k++)
      if (cnt[k] == 0)
      {
        int kb = 0;
        for (int q = 1; q < nd; q++)
          if (cnt[q] > cnt[kb]) kb = q;
        cnt[kb]--;
        cnt[k]++;
      }
}

} // namespace

extern "C"
{
  // 1 if this library can do the gravity with a tree (the Metal one and the stubs of the cpu build say 0)
  void gpu_tree_available_(int *ok)
  {
    *ok = 1;
  }

  // pointers to the host buffers for n particles and nnode nodes, which Fortran fills:
  //   pa(4,n) real*4: x y z (hi parts), mass;  pb(4,n) real*4: x y z (lo parts), 2h (negative for a point particle)
  //   need(n) integer*1;  nf(16,nnode) real*4;  ni(4,nnode) integer*4;  work(n) integer*4;
  //   out(4,n) real*8: the result, in the order of the tree
  // They are made again only if they are too small.
  void gpu_tree_reserve_(int *n_, int *nnode_, void **pa, void **pb, void **pneed, void **pnf, void **pni, void **pwork, void **pout)
  {
    const size_t n = *n_, nnode = *nnode_;
    if (n > cap_n || nnode > cap_node || h_pa == NULL)
    {
      free_host();
      const size_t nn = std::max(n + n / 8, (size_t)1024);
      const size_t nnd = std::max(nnode + nnode / 4, (size_t)1024);
      h_pa = (float4 *)hostalloc(nn * sizeof(float4));
      h_pb = (float4 *)hostalloc(nn * sizeof(float4));
      h_need = (unsigned char *)hostalloc(nn);
      h_nf = (float4 *)hostalloc(nnd * 4 * sizeof(float4));
      h_ni = (int4 *)hostalloc(nnd * sizeof(int4));
      h_work = (int *)hostalloc(nn * sizeof(int));
      h_out = (D4 *)hostalloc(nn * sizeof(D4));
      cap_n = nn;
      cap_node = nnd;
    }
    *pa = h_pa;
    *pb = h_pb;
    *pneed = h_need;
    *pnf = h_nf;
    *pni = h_ni;
    *pwork = h_work;
    *pout = h_out;
  }

  // Queue the copies, the kernel and the copy back of the result on the devices of the process, and return.
  void gpu_tree_launch_(int *n_, int *nnode_, int *nwork_, int *nkernel_, double *theta2_)
  {
    const int n = *n_, nnode = *nnode_, nwork = *nwork_, nk = *nkernel_;
    assert(!t_pending);
    read_env();
    int cur = 0;
    CUDA_SAFE_CALL(cudaGetDevice(&cur));
    int nd = SPHgrav::ndev > 1 ? SPHgrav::ndev : 1;
    const int first = nd > 1 ? SPHgrav::dev0 : cur;
    if ((int)tdev.size() != nd || t_first != first)
    {
      tdev.assign(nd, TreeDev());
      tspeed.assign(nd, 1.0);
      for (int k = 0; k < nd; k++)
        tdev[k].dev = first + k;
    }
    if (nd > 16)
    {
      fprintf(stderr, "\nStarSmasher: the tree gravity on the gpus is made for at most 16 devices per process, ngpuperrank gives %d.\n\n", nd);
      fflush(stderr);
      exit(EXIT_FAILURE);
    }
    if (t_host_st == NULL)
      t_host_st = (TreeStats *)hostalloc(16 * sizeof(TreeStats));
    t_first = first;
    t_nd = nd;
    t_nwork = nwork;
    std::vector<int> cnt;
    if (nd == 1)
      cnt.assign(1, nwork);
    else
      divide_work(nwork, nd, cnt);
    int off = 0;
    for (int k = 0; k < nd; k++)
    {
      TreeDev &d = tdev[k];
      d.off = off;
      d.cnt = cnt[k];
      off += cnt[k];
      d.nlaunch = 0;
      d.plo = 0;
      d.phi = -1;
      if (d.cnt > 0)
      {
        d.plo = h_ni[h_work[d.off]].z;
        d.phi = h_ni[h_work[d.off + d.cnt - 1]].w;
      }
    }
    for (int k = 0; k < nd; k++)
    {
      TreeDev &d = tdev[k];
      ensure_dev(d, (size_t)n, (size_t)nnode, (size_t)std::max(d.cnt, 1));
      if (d.cnt > 0)
      {
        CUDA_SAFE_CALL(cudaMemcpyAsync(d.pa, h_pa, (size_t)n * sizeof(float4), cudaMemcpyHostToDevice, 0));
        CUDA_SAFE_CALL(cudaMemcpyAsync(d.pb, h_pb, (size_t)n * sizeof(float4), cudaMemcpyHostToDevice, 0));
        CUDA_SAFE_CALL(cudaMemcpyAsync(d.need, h_need, (size_t)n, cudaMemcpyHostToDevice, 0));
        CUDA_SAFE_CALL(cudaMemcpyAsync(d.nf, h_nf, (size_t)nnode * 4 * sizeof(float4), cudaMemcpyHostToDevice, 0));
        CUDA_SAFE_CALL(cudaMemcpyAsync(d.ni, h_ni, (size_t)nnode * sizeof(int4), cudaMemcpyHostToDevice, 0));
        CUDA_SAFE_CALL(cudaMemcpyAsync(d.work, h_work + d.off, (size_t)d.cnt * sizeof(int), cudaMemcpyHostToDevice, 0));
        CUDA_SAFE_CALL(cudaMemsetAsync(d.counter, 0, sizeof(int), 0));
        CUDA_SAFE_CALL(cudaMemsetAsync(d.st, 0, sizeof(TreeStats), 0));
        CUDA_SAFE_CALL(cudaMemsetAsync(d.hcount, 0, 2 * sizeof(int), 0));
      }
      CUDA_SAFE_CALL(cudaEventRecord(d.e1, 0));
      if (d.cnt > 0)
      {
        TreeArgs args;
        args.pa = d.pa;
        args.pb = d.pb;
        args.need = d.need;
        args.nf = d.nf;
        args.ni = d.ni;
        args.work = d.work;
        args.nwork = d.cnt;
        args.counter = d.counter;
        args.theta2 = (float)(*theta2_);
        args.out = d.out;
        args.st = d.st;
        args.hv = t_hv;
        args.ch = TREE_CHUNK;
        args.poolcap = (int)std::min(d.poolcap, (size_t)0x3fffffff);
        args.hcount = d.hcount;
        args.hrec = d.hrec;
        args.pdesc = d.pdesc;
        args.hdone = d.hdone;
        args.part = d.part;
        d.a1 = args;
        d.a2 = args;
        d.a2.counter = d.counter2;
        d.nk = nk;
        launch_two(d);
        d.nlaunch = 1;
      }
      else
        CUDA_SAFE_CALL(cudaEventRecord(d.em, 0));
      if (d.cnt > 0)
      {
        CUDA_SAFE_CALL(cudaMemcpyAsync(h_out + d.plo, d.out + d.plo, (size_t)(d.phi - d.plo + 1) * sizeof(D4), cudaMemcpyDeviceToHost, 0));
        CUDA_SAFE_CALL(cudaMemcpyAsync(&t_host_st[k], d.st, sizeof(TreeStats), cudaMemcpyDeviceToHost, 0));
      }
      CUDA_SAFE_CALL(cudaEventRecord(d.e3, 0));
    }
    CUDA_SAFE_CALL(cudaSetDevice(cur));
    t_pending = 1;
  }

  // Wait for the call in flight.  ierr: 0, 1 if the stack had no room, 2 if the pool of the chunks had none.
  void gpu_tree_wait_(int *ierr)
  {
    int cur = 0;
    CUDA_SAFE_CALL(cudaGetDevice(&cur));
    *ierr = 0;
    if (!t_pending)
      return;
    std::vector<double> k1ms(t_nd, 0.0);
    for (int k = 0; k < t_nd; k++)
    {
      TreeDev &d = tdev[k];
      CUDA_SAFE_CALL(cudaSetDevice(d.dev));
      CUDA_SAFE_CALL(cudaEventSynchronize(d.e3));
      float b1 = 0.0f;
      CUDA_SAFE_CALL(cudaEventElapsedTime(&b1, d.e1, d.em));
      // the pool of the chunks was full: the call is done again on this device with a pool that has room for all (the result is the same)
      int tries = 0;
      while (d.cnt > 0 && (t_host_st[k].err & 2) && tries < 3)
      {
        const size_t need = (size_t)t_host_st[k].nchunks;
        fprintf(stderr, "gpu_tree: the pool of %lu chunks of the heavy leaves was full (%lu needed): the call is done again on device %d\n", (unsigned long)d.poolcap,
                (unsigned long)need, d.dev);
        set_pool(d, need + need / 4 + 64);
        d.a1.poolcap = (int)std::min(d.poolcap, (size_t)0x3fffffff);
        d.a1.part = d.part;
        d.a1.pdesc = d.pdesc;
        d.a2.poolcap = d.a1.poolcap;
        d.a2.part = d.part;
        d.a2.pdesc = d.pdesc;
        CUDA_SAFE_CALL(cudaMemsetAsync(d.counter, 0, sizeof(int), 0));
        CUDA_SAFE_CALL(cudaMemsetAsync(d.st, 0, sizeof(TreeStats), 0));
        CUDA_SAFE_CALL(cudaMemsetAsync(d.hcount, 0, 2 * sizeof(int), 0));
        launch_two(d);
        CUDA_SAFE_CALL(cudaMemcpyAsync(h_out + d.plo, d.out + d.plo, (size_t)(d.phi - d.plo + 1) * sizeof(D4), cudaMemcpyDeviceToHost, 0));
        CUDA_SAFE_CALL(cudaMemcpyAsync(&t_host_st[k], d.st, sizeof(TreeStats), cudaMemcpyDeviceToHost, 0));
        CUDA_SAFE_CALL(cudaDeviceSynchronize());
        tries++;
      }
      k1ms[k] = b1;
      if (d.cnt > 0 && t_host_st[k].err)
        *ierr = (t_host_st[k].err & 1) ? 1 : 2;
    }
    // the speed of the next call: leaves per ms of the first pass, relative to the fastest, half old and half new.  The second pass (the chunks of the heavy leaves) is not counted: its work is
    // not proportional to the leaves, and a device that held heavy leaves looked slow, was given fewer leaves, and its second pass ran on fewer warps (the step of the collision of two stars of
    // 500k particles with eight GPUs was 16% longer)
    if (t_nd > 1)
    {
      double smax = 0.0;
      std::vector<double> v(t_nd, 0.0);
      for (int k = 0; k < t_nd; k++)
      {
        v[k] = (tdev[k].cnt > 0 && k1ms[k] > 0.0) ? tdev[k].cnt / k1ms[k] : 0.0;
        smax = std::max(smax, v[k]);
      }
      for (int k = 0; k < t_nd; k++)
        if (v[k] > 0.0 && smax > 0.0)
          tspeed[k] = 0.5 * tspeed[k] + 0.5 * v[k] / smax;
    }
    CUDA_SAFE_CALL(cudaSetDevice(cur));
    t_pending = 0;
  }
}
