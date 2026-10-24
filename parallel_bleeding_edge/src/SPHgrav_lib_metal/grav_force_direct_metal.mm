#include <Metal/Metal.h>
#include <Foundation/Foundation.h>

#include <cassert>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

namespace {

constexpr NSUInteger kThreadsPerGroup = 256;
constexpr float kLargeCoordinate = 1.0e10f;

struct FloatPair {
  float hi;
  float lo;
};

struct MetalParticle {
  float px;
  float py;
  float pz;
  float mass;
  float h;
  float invh;
  float h2;
  float hflag;
};

// ngravhilo=1: the float remainder of each position, in a buffer of its own
// so that ngravhilo=0 reads and computes what it did before
struct MetalLo {
  float lx;
  float ly;
  float lz;
  float pad;
};

struct MetalForce {
  FloatPair ax;
  FloatPair ay;
  FloatPair az;
  FloatPair pot;
};

struct MetalParams {
  uint32_t nj;
  uint32_t nibeg;
  uint32_t ni;
  int32_t nkernel;
  uint32_t uselist;
};

static const char *kSPHgravMetalSource = R"METAL(
#include <metal_stdlib>
using namespace metal;

#define THREADS_PER_GROUP 256

struct FloatPair {
  float hi;
  float lo;
};

struct Particle {
  float px;
  float py;
  float pz;
  float mass;
  float h;
  float invh;
  float h2;
  float hflag;
};

struct Force {
  FloatPair ax;
  FloatPair ay;
  FloatPair az;
  FloatPair pot;
};

struct Params {
  uint nj;
  uint nibeg;
  uint ni;
  int nkernel;
  uint uselist;
};

inline void ds64_add(thread FloatPair &accum, const float x) {
  const float vx = accum.hi + x;
  const float vy = accum.lo - ((vx - accum.hi) - x);
  accum.hi = vx;
  accum.lo = vy;
}

// ngravhilo: set when the pipeline is built
constant bool kHilo [[function_constant(0)]];

inline float3 position(const Particle p) {
  return float3(p.px, p.py, p.pz);
}

// With kHilo each coordinate also carries the float remainder (ilo, jlo) of
// its double value.  Nearby floats subtract exactly, so the separation keeps
// float precision relative to itself however far the pair is from the
// reference point.
inline float3 separation(const Particle pi, const Particle pj, const float3 ilo, const float3 jlo) {
  return kHilo ? (position(pj) - position(pi)) + (jlo - ilo) : position(pj) - position(pi);
}

inline Force zero_force() {
  Force f;
  f.ax = {0.0f, 0.0f};
  f.ay = {0.0f, 0.0f};
  f.az = {0.0f, 0.0f};
  f.pot = {0.0f, 0.0f};
  return f;
}

inline Force force_ij_nkernel0(const Particle pi, const Particle pj, const float3 ilo, const float3 jlo,
                                thread Force iforce) {
  const float3 dr = separation(pi, pj, ilo, jlo);
  const float r2 = dot(dr, dr);

  const float rinv = r2 > 0.0f ? rsqrt(r2) : 0.0f;
  const float rinv2 = rinv * rinv;
  const float mrinv1 = rinv * pj.mass;
  const float mrinv3 = rinv2 * mrinv1;

  if (r2 >= max(pi.h2, pj.h2)) {
    ds64_add(iforce.ax, mrinv3 * dr.x);
    ds64_add(iforce.ay, mrinv3 * dr.y);
    ds64_add(iforce.az, mrinv3 * dr.z);
    ds64_add(iforce.pot, -mrinv1);
  } else if (r2 < min(pi.h2, pj.h2)) {
    const float2 invq = r2 > 0.0f ? float2(rinv * pi.h, rinv * pj.h) : float2(0.0f);
    const float2 q = r2 > 0.0f ? float2(1.0f / invq.x, 1.0f / invq.y) : float2(0.0f);
    const float2 q2 = float2(q.x * q.x, q.y * q.y);
    const float2 q3 = float2(q.x * q2.x, q.y * q2.y);
    const float2 invq2 = float2(invq.x * invq.x, invq.y * invq.y);
    const float2 invq3 = float2(invq.x * invq2.x, invq.y * invq2.y);

    const float2 f = float2(q.x < 0.5f ? 1.0f : 0.0f, q.y < 0.5f ? 1.0f : 0.0f);
    const float2 acc = float2(
        f.x * (10.666666666667f + q2.x * (32.0f * q.x - 38.4f)) +
            (1.0f - f.x) * (21.333333333333f - 48.0f * q.x + 38.4f * q2.x -
                             10.666666666667f * q3.x - 0.066666666667f * invq3.x),
        f.y * (10.666666666667f + q2.y * (32.0f * q.y - 38.4f)) +
            (1.0f - f.y) * (21.333333333333f - 48.0f * q.y + 38.4f * q2.y -
                             10.666666666667f * q3.y - 0.066666666667f * invq3.y));

    const float2 pot = float2(
        f.x * (-2.8f + q2.x * (5.333333333333f + q2.x * (6.4f * q.x - 9.6f))) +
            (1.0f - f.x) *
                (-3.2f + 0.066666666667f * invq.x +
                 q2.x * (10.666666666667f + q.x * (-16.0f + q.x * (9.6f - 2.133333333333f * q.x)))),
        f.y * (-2.8f + q2.y * (5.333333333333f + q2.y * (6.4f * q.y - 9.6f))) +
            (1.0f - f.y) *
                (-3.2f + 0.066666666667f * invq.y +
                 q2.y * (10.666666666667f + q.y * (-16.0f + q.y * (9.6f - 2.133333333333f * q.y)))));

    const float2 mj1 = float2(pj.mass * pi.invh, pj.mass * pj.invh);
    const float2 mj2 = float2(mj1.x * pi.invh * pi.invh, mj1.y * pj.invh * pj.invh);
    const float2 g = float2(r2 <= pi.h2 ? 1.0f : 0.0f, r2 <= pj.h2 ? 1.0f : 0.0f);
    const float gacc = r2 > 0.0f ? 0.5f * (g.x * mj2.x * acc.x + (1.0f - g.x) * mrinv3 +
                                            g.y * mj2.y * acc.y + (1.0f - g.y) * mrinv3)
                                  : 0.0f;
    const float gpot =
        (r2 > 0.0f || pi.h2 != pj.h2 || pi.mass != pj.mass)
            ? 0.5f * (g.x * mj1.x * pot.x + (g.x - 1.0f) * mrinv1 + g.y * mj1.y * pot.y +
                      (g.y - 1.0f) * mrinv1)
            : (-1.4f) * pj.hflag * (mj1.x + mj1.y);

    ds64_add(iforce.ax, gacc * dr.x);
    ds64_add(iforce.ay, gacc * dr.y);
    ds64_add(iforce.az, gacc * dr.z);
    ds64_add(iforce.pot, gpot);
  } else if (r2 < pi.h2) {
    const float invq = r2 > 0.0f ? rinv * pi.h : 0.0f;
    const float q = r2 > 0.0f ? 1.0f / invq : 0.0f;
    const float q2 = q * q;
    const float q3 = q * q2;
    const float invq2 = invq * invq;
    const float invq3 = invq * invq2;
    const float f = q < 0.5f ? 1.0f : 0.0f;
    const float acc = f * (10.666666666667f + q2 * (32.0f * q - 38.4f)) +
                      (1.0f - f) *
                          (21.333333333333f - 48.0f * q + 38.4f * q2 - 10.666666666667f * q3 -
                           0.066666666667f * invq3);
    const float pot = f * (-2.8f + q2 * (5.333333333333f + q2 * (6.4f * q - 9.6f))) +
                      (1.0f - f) *
                          (-3.2f + 0.066666666667f * invq +
                           q2 * (10.666666666667f + q * (-16.0f + q * (9.6f - 2.133333333333f * q))));
    const float2 mj1 = float2(pj.mass * pi.invh, pj.mass * pj.invh);
    const float mj2 = mj1.x * pi.invh * pi.invh;
    const float gacc = r2 > 0.0f ? 0.5f * (mj2 * acc + mrinv3) : 0.0f;
    const float gpot = r2 > 0.0f ? 0.5f * (mj1.x * pot - mrinv1) : (-1.4f) * pj.hflag * (mj1.x + mj1.y);

    ds64_add(iforce.ax, gacc * dr.x);
    ds64_add(iforce.ay, gacc * dr.y);
    ds64_add(iforce.az, gacc * dr.z);
    ds64_add(iforce.pot, gpot);
  } else {
    const float invq = r2 > 0.0f ? rinv * pj.h : 0.0f;
    const float q = r2 > 0.0f ? 1.0f / invq : 0.0f;
    const float q2 = q * q;
    const float q3 = q * q2;
    const float invq2 = invq * invq;
    const float invq3 = invq * invq2;
    const float f = q < 0.5f ? 1.0f : 0.0f;
    const float acc = f * (10.666666666667f + q2 * (32.0f * q - 38.4f)) +
                      (1.0f - f) *
                          (21.333333333333f - 48.0f * q + 38.4f * q2 - 10.666666666667f * q3 -
                           0.066666666667f * invq3);
    const float pot = f * (-2.8f + q2 * (5.333333333333f + q2 * (6.4f * q - 9.6f))) +
                      (1.0f - f) *
                          (-3.2f + 0.066666666667f * invq +
                           q2 * (10.666666666667f + q * (-16.0f + q * (9.6f - 2.133333333333f * q))));
    const float2 mj1 = float2(pj.mass * pi.invh, pj.mass * pj.invh);
    const float mj2 = mj1.y * pj.invh * pj.invh;
    const float gacc = r2 > 0.0f ? 0.5f * (mrinv3 + mj2 * acc) : 0.0f;
    const float gpot = r2 > 0.0f ? 0.5f * (-mrinv1 + mj1.y * pot) : (-1.4f) * pj.hflag * (mj1.x + mj1.y);

    ds64_add(iforce.ax, gacc * dr.x);
    ds64_add(iforce.ay, gacc * dr.y);
    ds64_add(iforce.az, gacc * dr.z);
    ds64_add(iforce.pot, gpot);
  }

  return iforce;
}

inline Force force_ij_nkernel1(const Particle pi, const Particle pj, const float3 ilo, const float3 jlo,
                                thread Force iforce) {
  const float3 dr = separation(pi, pj, ilo, jlo);
  const float r2 = dot(dr, dr);

  const float rinv = r2 > 0.0f ? rsqrt(r2) : 0.0f;
  const float rinv2 = rinv * rinv;
  const float mrinv1 = rinv * pj.mass;
  const float mrinv3 = rinv2 * mrinv1;

  if (r2 >= max(pi.h2, pj.h2)) {
    ds64_add(iforce.ax, mrinv3 * dr.x);
    ds64_add(iforce.ay, mrinv3 * dr.y);
    ds64_add(iforce.az, mrinv3 * dr.z);
    ds64_add(iforce.pot, -mrinv1);
  } else if (r2 < min(pi.h2, pj.h2)) {
    const float2 invq = r2 > 0.0f ? float2(rinv * pi.h, rinv * pj.h) : float2(0.0f);
    const float2 q = r2 > 0.0f ? float2(1.0f / invq.x, 1.0f / invq.y) : float2(0.0f);
    const float2 q2 = float2(q.x * q.x, q.y * q.y);

    const float2 acc = float2(
        28.4375f + q2.x * (-187.6875f + q2.x * (804.375f + q2.x * (-4379.375f + q.x * (9009.0f + q.x * (-8957.8125f + q.x * (5005.0f + q.x * (-1515.9375f + 195.0f * q.x))))))),
        28.4375f + q2.y * (-187.6875f + q2.y * (804.375f + q2.y * (-4379.375f + q.y * (9009.0f + q.y * (-8957.8125f + q.y * (5005.0f + q.y * (-1515.9375f + 195.0f * q.y))))))));

    const float2 pot = float2(
        r2 <= pi.h2 ? -3.828125f + q2.x * (14.21875f + q2.x * (-46.921875f + q2.x * (134.0625f + q2.x * (-547.421875f + q.x * (1001.0f + q.x * (-895.78125f + q.x * (455.0f + q.x * (-126.328125f + 15.0f * q.x)))))))) : 0.0f,
        r2 <= pj.h2 ? -3.828125f + q2.y * (14.21875f + q2.y * (-46.921875f + q2.y * (134.0625f + q2.y * (-547.421875f + q.y * (1001.0f + q.y * (-895.78125f + q.y * (455.0f + q.y * (-126.328125f + 15.0f * q.y)))))))) : 0.0f);

    const float2 mj1 = float2(pj.mass * pi.invh, pj.mass * pj.invh);
    const float2 mj2 = float2(mj1.x * pi.invh * pi.invh, mj1.y * pj.invh * pj.invh);
    const float2 g = float2(r2 <= pi.h2 ? 1.0f : 0.0f, r2 <= pj.h2 ? 1.0f : 0.0f);
    const float gacc = r2 > 0.0f ? 0.5f * (g.x * mj2.x * acc.x + (1.0f - g.x) * mrinv3 +
                                            g.y * mj2.y * acc.y + (1.0f - g.y) * mrinv3)
                                  : 0.0f;
    const float gpot =
        (r2 > 0.0f || pi.h2 != pj.h2 || pi.mass != pj.mass)
            ? 0.5f * (mj1.x * pot.x + (g.x - 1.0f) * mrinv1 + mj1.y * pot.y + (g.y - 1.0f) * mrinv1)
            : (-1.9140625f) * pj.hflag * (mj1.x + mj1.y);

    ds64_add(iforce.ax, gacc * dr.x);
    ds64_add(iforce.ay, gacc * dr.y);
    ds64_add(iforce.az, gacc * dr.z);
    ds64_add(iforce.pot, gpot);
  } else if (r2 < pi.h2) {
    const float invq = r2 > 0.0f ? rinv * pi.h : 0.0f;
    const float q = r2 > 0.0f ? 1.0f / invq : 0.0f;
    const float q2 = q * q;
    const float acc = 28.4375f + q2 * (-187.6875f + q2 * (804.375f + q2 * (-4379.375f + q * (9009.0f + q * (-8957.8125f + q * (5005.0f + q * (-1515.9375f + 195.0f * q)))))));
    const float pot = -3.828125f + q2 * (14.21875f + q2 * (-46.921875f + q2 * (134.0625f + q2 * (-547.421875f + q * (1001.0f + q * (-895.78125f + q * (455.0f + q * (-126.328125f + 15.0f * q))))))));
    const float2 mj1 = float2(pj.mass * pi.invh, pj.mass * pj.invh);
    const float mj2 = mj1.x * pi.invh * pi.invh;
    const float gacc = r2 > 0.0f ? 0.5f * (mj2 * acc + mrinv3) : 0.0f;
    const float gpot = r2 > 0.0f ? 0.5f * (mj1.x * pot - mrinv1) : (-1.9140625f) * pj.hflag * (mj1.x + mj1.y);

    ds64_add(iforce.ax, gacc * dr.x);
    ds64_add(iforce.ay, gacc * dr.y);
    ds64_add(iforce.az, gacc * dr.z);
    ds64_add(iforce.pot, gpot);
  } else {
    const float invq = r2 > 0.0f ? rinv * pj.h : 0.0f;
    const float q = r2 > 0.0f ? 1.0f / invq : 0.0f;
    const float q2 = q * q;
    const float acc = 28.4375f + q2 * (-187.6875f + q2 * (804.375f + q2 * (-4379.375f + q * (9009.0f + q * (-8957.8125f + q * (5005.0f + q * (-1515.9375f + 195.0f * q)))))));
    const float pot = -3.828125f + q2 * (14.21875f + q2 * (-46.921875f + q2 * (134.0625f + q2 * (-547.421875f + q * (1001.0f + q * (-895.78125f + q * (455.0f + q * (-126.328125f + 15.0f * q))))))));
    const float2 mj1 = float2(pj.mass * pi.invh, pj.mass * pj.invh);
    const float mj2 = mj1.y * pj.invh * pj.invh;
    const float gacc = r2 > 0.0f ? 0.5f * (mrinv3 + mj2 * acc) : 0.0f;
    const float gpot = r2 > 0.0f ? 0.5f * (-mrinv1 + mj1.y * pot) : (-1.9140625f) * pj.hflag * (mj1.x + mj1.y);

    ds64_add(iforce.ax, gacc * dr.x);
    ds64_add(iforce.ay, gacc * dr.y);
    ds64_add(iforce.az, gacc * dr.z);
    ds64_add(iforce.pot, gpot);
  }

  return iforce;
}

inline Force force_ij_nkernel2(const Particle pi, const Particle pj, const float3 ilo, const float3 jlo,
                                thread Force iforce) {
  const float3 dr = separation(pi, pj, ilo, jlo);
  const float r2 = dot(dr, dr);

  const float rinv = r2 > 0.0f ? rsqrt(r2) : 0.0f;
  const float rinv2 = rinv * rinv;
  const float mrinv1 = rinv * pj.mass;
  const float mrinv3 = rinv2 * mrinv1;

  if (r2 >= max(pi.h2, pj.h2)) {
    ds64_add(iforce.ax, mrinv3 * dr.x);
    ds64_add(iforce.ay, mrinv3 * dr.y);
    ds64_add(iforce.az, mrinv3 * dr.z);
    ds64_add(iforce.pot, -mrinv1);
  } else if (r2 < min(pi.h2, pj.h2)) {
    const float2 invq = r2 > 0.0f ? float2(rinv * pi.h, rinv * pj.h) : float2(0.0f);
    const float2 q = r2 > 0.0f ? float2(1.0f / invq.x, 1.0f / invq.y) : float2(0.0f);
    const float2 q2 = float2(q.x * q.x, q.y * q.y);

    const float2 acc = float2(
        20.625f + q2.x * (-115.5f + q2.x * (618.75f + q.x * (-1155.0f + q.x * (962.5f + q.x * (-396.0f + 65.625f * q.x))))),
        20.625f + q2.y * (-115.5f + q2.y * (618.75f + q.y * (-1155.0f + q.y * (962.5f + q.y * (-396.0f + 65.625f * q.y))))));

    const float2 pot = float2(
        -3.4375f + q2.x * (10.3125f + q2.x * (-28.875f + q2.x * (103.125f + q.x * (-165.0f + q.x * (120.3125f + q.x * (-44.0f + 6.5625f * q.x)))))),
        -3.4375f + q2.y * (10.3125f + q2.y * (-28.875f + q2.y * (103.125f + q.y * (-165.0f + q.y * (120.3125f + q.y * (-44.0f + 6.5625f * q.y)))))));

    const float2 mj1 = float2(pj.mass * pi.invh, pj.mass * pj.invh);
    const float2 mj2 = float2(mj1.x * pi.invh * pi.invh, mj1.y * pj.invh * pj.invh);
    const float gacc = r2 > 0.0f ? 0.5f * (mj2.x * acc.x + mj2.y * acc.y) : 0.0f;
    const float gpot =
        (r2 > 0.0f || pi.h2 != pj.h2 || pi.mass != pj.mass)
            ? 0.5f * (mj1.x * pot.x + mj1.y * pot.y)
            : (-1.71875f) * pj.hflag * (mj1.x + mj1.y);

    ds64_add(iforce.ax, gacc * dr.x);
    ds64_add(iforce.ay, gacc * dr.y);
    ds64_add(iforce.az, gacc * dr.z);
    ds64_add(iforce.pot, gpot);
  } else if (r2 < pi.h2) {
    const float invq = r2 > 0.0f ? rinv * pi.h : 0.0f;
    const float q = r2 > 0.0f ? 1.0f / invq : 0.0f;
    const float q2 = q * q;
    const float acc = 20.625f + q2 * (-115.5f + q2 * (618.75f + q * (-1155.0f + q * (962.5f + q * (-396.0f + 65.625f * q)))));
    const float pot = -3.4375f + q2 * (10.3125f + q2 * (-28.875f + q2 * (103.125f + q * (-165.0f + q * (120.3125f + q * (-44.0f + 6.5625f * q))))));
    const float2 mj1 = float2(pj.mass * pi.invh, pj.mass * pj.invh);
    const float mj2 = mj1.x * pi.invh * pi.invh;
    const float gacc = r2 > 0.0f ? 0.5f * (mj2 * acc + mrinv3) : 0.0f;
    const float gpot = r2 > 0.0f ? 0.5f * (mj1.x * pot - mrinv1) : (-1.71875f) * pj.hflag * (mj1.x + mj1.y);

    ds64_add(iforce.ax, gacc * dr.x);
    ds64_add(iforce.ay, gacc * dr.y);
    ds64_add(iforce.az, gacc * dr.z);
    ds64_add(iforce.pot, gpot);
  } else {
    const float invq = r2 > 0.0f ? rinv * pj.h : 0.0f;
    const float q = r2 > 0.0f ? 1.0f / invq : 0.0f;
    const float q2 = q * q;
    const float acc = 20.625f + q2 * (-115.5f + q2 * (618.75f + q * (-1155.0f + q * (962.5f + q * (-396.0f + 65.625f * q)))));
    const float pot = -3.4375f + q2 * (10.3125f + q2 * (-28.875f + q2 * (103.125f + q * (-165.0f + q * (120.3125f + q * (-44.0f + 6.5625f * q))))));
    const float2 mj1 = float2(pj.mass * pi.invh, pj.mass * pj.invh);
    const float mj2 = mj1.y * pj.invh * pj.invh;
    const float gacc = r2 > 0.0f ? 0.5f * (mrinv3 + mj2 * acc) : 0.0f;
    const float gpot = r2 > 0.0f ? 0.5f * (-mrinv1 + mj1.y * pot) : (-1.71875f) * pj.hflag * (mj1.x + mj1.y);

    ds64_add(iforce.ax, gacc * dr.x);
    ds64_add(iforce.ay, gacc * dr.y);
    ds64_add(iforce.az, gacc * dr.z);
    ds64_add(iforce.pot, gpot);
  }

  return iforce;
}

kernel void compute_forces(const device Particle *ptcl [[buffer(0)]],
                           device Force *force_out [[buffer(1)]],
                           constant Params &params [[buffer(2)]],
                           const device uint *ilist [[buffer(3)]],
                           const device float4 *plo [[buffer(4)]],
                           uint gid [[thread_position_in_grid]],
                           uint tid [[thread_index_in_threadgroup]]) {
  threadgroup Particle jpshared[THREADS_PER_GROUP];
  threadgroup float3 jloshared[THREADS_PER_GROUP];

  // uselist (block timesteps): force_out[gid] is for particle ilist[gid]
  // rather than for particle nibeg+gid
  Particle iptcl;
  float3 ilo = float3(0.0f);
  if (gid < params.ni) {
    const uint iidx = params.uselist != 0 ? ilist[gid] : params.nibeg + gid;
    iptcl = ptcl[iidx];
    if (kHilo) {
      ilo = plo[iidx].xyz;
    }
  } else {
    iptcl = {1.0e10f, 1.0e10f, 1.0e10f, 0.0f, 1.0f, 1.0f, 1.0f, 0.0f};
  }

  Force iforce = zero_force();

  for (uint j = 0; j < params.nj; j += THREADS_PER_GROUP) {
    jpshared[tid] = ptcl[j + tid];
    if (kHilo) {
      jloshared[tid] = plo[j + tid].xyz;
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);

    if (gid < params.ni) {
      for (uint jj = 0; jj < THREADS_PER_GROUP; ++jj) {
        if (params.nkernel == 0) {
          iforce = force_ij_nkernel0(iptcl, jpshared[jj], ilo, jloshared[jj], iforce);
        } else if (params.nkernel == 1) {
          iforce = force_ij_nkernel1(iptcl, jpshared[jj], ilo, jloshared[jj], iforce);
        } else {
          iforce = force_ij_nkernel2(iptcl, jpshared[jj], ilo, jloshared[jj], iforce);
        }
      }
    }

    threadgroup_barrier(mem_flags::mem_threadgroup);
  }

  if (gid < params.ni) {
    force_out[gid] = iforce;
  }
}
)METAL";

class SPHgravMetal {
 public:
  SPHgravMetal() = default;

  void setDevice(int requestedDevice) {
    @autoreleasepool {
      NSArray *devices = MTLCopyAllDevices();
      if (devices == nil || [devices count] == 0) {
        id<MTLDevice> fallback = MTLCreateSystemDefaultDevice();
        if (fallback == nil) {
          std::fprintf(stderr, "SPHgrav Metal could not find a Metal device\n");
          std::abort();
        }
        devices = @[ fallback ];
      }

      const NSUInteger deviceCount = [devices count];
      const NSUInteger deviceIndex =
          deviceCount > 0 ? static_cast<NSUInteger>(requestedDevice % static_cast<int>(deviceCount)) : 0;
      id<MTLDevice> nextDevice = [devices objectAtIndex:deviceIndex];

      if (device_ != nil && [device_ isEqual:nextDevice]) {
        return;
      }

      device_ = nextDevice;
      commandQueue_ = [device_ newCommandQueue];
      buildPipeline();

      std::fprintf(stderr, " SPHgrav Metal found %lu device(s)\n", static_cast<unsigned long>(deviceCount));
      for (NSUInteger i = 0; i < deviceCount; ++i) {
        id<MTLDevice> dev = [devices objectAtIndex:i];
        std::fprintf(stderr, "  Device= %lu: %s%s\n",
                     static_cast<unsigned long>(i),
                     [[dev name] UTF8String],
                     i == deviceIndex ? " [selected]" : "");
      }
    }
  }

  void setHilo(bool hilo) {
    hilo_ = hilo;
  }

  void first_half(int nj, int nibeg, int ni,
                  double *px, double *py, double *pz,
                  double *mass, double *h2, int nkernel) {
    assert(nibeg >= 0);
    const int paddedCount = loadParticles(nj, px, py, pz, mass, h2);
    pendingList_.clear();
    pendingNrange_ = 0;
    dispatch(paddedCount, nibeg, ni, nkernel, false);
  }

  // Block timesteps: forces only on the nlist particles in list (1-based,
  // all within nibeg+1..nibeg+nrange), out of all nj.  last_half zeroes the
  // others in that range.
  void first_half_list(int nj, int nibeg, int nrange, int nlist, const int *list,
                       double *px, double *py, double *pz,
                       double *mass, double *h2, int nkernel) {
    assert(nibeg >= 0);
    const int paddedCount = loadParticles(nj, px, py, pz, mass, h2);
    pendingList_.assign(list, list + nlist);
    for (int &k : pendingList_) {
      k -= 1;
    }
    pendingNrange_ = nrange;
    if (nlist > 0) {
      resizeListBuffer(nlist);
      auto *ilist = static_cast<uint32_t *>([listBuffer_ contents]);
      for (int i = 0; i < nlist; ++i) {
        ilist[i] = static_cast<uint32_t>(pendingList_[i]);
      }
    }
    dispatch(paddedCount, nibeg, nlist, nkernel, true);
  }

  void last_half(double *ax, double *ay, double *az, double *pot) {
    if (pendingUseList_) {
      for (int j = pendingNibeg_; j < pendingNibeg_ + pendingNrange_; ++j) {
        ax[j] = ay[j] = az[j] = pot[j] = 0.0;
      }
    }

    if (pendingNi_ > 0) {
      assert(pendingCommandBuffer_ != nil);
      [pendingCommandBuffer_ waitUntilCompleted];

      if ([pendingCommandBuffer_ status] != MTLCommandBufferStatusCompleted) {
        NSString *errorDescription = [[pendingCommandBuffer_ error] localizedDescription];
        std::fprintf(stderr, "SPHgrav Metal command buffer failed: %s\n",
                     errorDescription != nil ? [errorDescription UTF8String] : "unknown error");
        std::abort();
      }

      const auto *force = static_cast<const MetalForce *>([forceBuffer_ contents]);
      for (int i = 0; i < pendingNi_; ++i) {
        const int j = pendingUseList_ ? pendingList_[i] : i + pendingNibeg_;
        ax[j] = static_cast<double>(force[i].ax.hi) + static_cast<double>(force[i].ax.lo);
        ay[j] = static_cast<double>(force[i].ay.hi) + static_cast<double>(force[i].ay.lo);
        az[j] = static_cast<double>(force[i].az.hi) + static_cast<double>(force[i].az.lo);
        pot[j] = static_cast<double>(force[i].pot.hi) + static_cast<double>(force[i].pot.lo);
      }
    }

    const auto elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - lastStart_).count();
    std::fprintf(stderr, " >>> SPHgrav_lib took %g sec\n", elapsed);

    pendingCommandBuffer_ = nil;
    pendingNi_ = 0;
    pendingNibeg_ = 0;
    pendingNrange_ = 0;
    pendingUseList_ = false;
  }

 private:
  // Copies the particles into the shared buffer, padded to a whole number of
  // threadgroups, and returns the padded count.
  int loadParticles(int nj, double *px, double *py, double *pz,
                    double *mass, double *h2) {
    assert(device_ != nil);
    assert(commandQueue_ != nil);
    assert(pipelineState_ != nil);

    const int paddedCount = ((nj - 1) / static_cast<int>(kThreadsPerGroup) + 1) *
                            static_cast<int>(kThreadsPerGroup);

    resizeParticleBuffer(paddedCount);
    // the kernel declares the remainder buffer, so it is bound even when unused
    resizeLoBuffer(hilo_ ? paddedCount : 1);

    // Positions relative to the gas centre of mass, subtracted in double,
    // and with ngravhilo=1 their float remainders, as in
    // SPHgrav_lib2/grav_force_direct.cu (see the comment there).
    double wsum = 0.0, xref = 0.0, yref = 0.0, zref = 0.0;
    for (int i = 0; i < nj; ++i) {
      if (h2[i] > 0.0) {
        const double w = mass[i];
        wsum += w;
        xref += w * px[i];
        yref += w * py[i];
        zref += w * pz[i];
      }
    }
    if (wsum > 0.0) {
      xref /= wsum;
      yref /= wsum;
      zref /= wsum;
    }

    auto *particles = static_cast<MetalParticle *>([particleBuffer_ contents]);
    auto *lo = static_cast<MetalLo *>([loBuffer_ contents]);
    for (int i = 0; i < nj; ++i) {
      const float smoothing = std::sqrt(static_cast<float>(h2[i] > 0.0 ? h2[i] : -h2[i]));
      const double dx = px[i] - xref, dy = py[i] - yref, dz = pz[i] - zref;
      const float fx = static_cast<float>(dx), fy = static_cast<float>(dy), fz = static_cast<float>(dz);
      particles[i] = {
          fx,
          fy,
          fz,
          static_cast<float>(mass[i]),
          smoothing,
          smoothing > 0.0f ? 1.0f / smoothing : 0.0f,
          static_cast<float>(h2[i] > 0.0 ? h2[i] : -h2[i]),
          h2[i] > 0.0 ? 1.0f : 0.0f,
      };
      if (hilo_) {
        lo[i] = {static_cast<float>(dx - fx), static_cast<float>(dy - fy), static_cast<float>(dz - fz), 0.0f};
      }
    }

    for (int i = nj; i < paddedCount; ++i) {
      particles[i] = {kLargeCoordinate, kLargeCoordinate, kLargeCoordinate, 0.0f, 1.0f, 1.0f, 1.0f, 0.0f};
      if (hilo_) {
        lo[i] = {0.0f, 0.0f, 0.0f, 0.0f};
      }
    }
    return paddedCount;
  }

  // Forces on ni particles: nibeg+0..ni-1, or listBuffer_[0..ni-1] if uselist.
  // With ni == 0 nothing is dispatched (an empty grid is not allowed).
  void dispatch(int paddedCount, int nibeg, int ni, int nkernel, bool uselist) {
    pendingNi_ = ni;
    pendingNibeg_ = nibeg;
    pendingUseList_ = uselist;
    lastStart_ = std::chrono::steady_clock::now();
    if (ni == 0) {
      return;
    }

    resizeForceBuffer(ni);
    // the kernel declares the list buffer, so it is bound even when unused
    resizeListBuffer(1);

    MetalParams params = {
        static_cast<uint32_t>(paddedCount),
        static_cast<uint32_t>(nibeg),
        static_cast<uint32_t>(ni),
        static_cast<int32_t>(nkernel),
        uselist ? 1u : 0u,
    };
    std::memcpy([paramBuffer_ contents], &params, sizeof(params));

    @autoreleasepool {
      pendingCommandBuffer_ = [commandQueue_ commandBuffer];
      id<MTLComputeCommandEncoder> encoder = [pendingCommandBuffer_ computeCommandEncoder];
      [encoder setComputePipelineState:hilo_ ? pipelineStateHilo_ : pipelineState_];
      [encoder setBuffer:particleBuffer_ offset:0 atIndex:0];
      [encoder setBuffer:forceBuffer_ offset:0 atIndex:1];
      [encoder setBuffer:paramBuffer_ offset:0 atIndex:2];
      [encoder setBuffer:listBuffer_ offset:0 atIndex:3];
      [encoder setBuffer:loBuffer_ offset:0 atIndex:4];

      const NSUInteger nblocks =
          (static_cast<NSUInteger>(ni) + kThreadsPerGroup - 1) / kThreadsPerGroup;
      const MTLSize threadgroups = MTLSizeMake(nblocks, 1, 1);
      const MTLSize threadsPerThreadgroup = MTLSizeMake(kThreadsPerGroup, 1, 1);
      [encoder dispatchThreadgroups:threadgroups threadsPerThreadgroup:threadsPerThreadgroup];
      [encoder endEncoding];
      [pendingCommandBuffer_ commit];
    }
  }

  void buildPipeline() {
    @autoreleasepool {
      NSString *source = [NSString stringWithUTF8String:kSPHgravMetalSource];
      MTLCompileOptions *options = [[MTLCompileOptions alloc] init];
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
      [options setFastMathEnabled:YES];
#pragma clang diagnostic pop

      NSError *libraryError = nil;
      id<MTLLibrary> library = [device_ newLibraryWithSource:source options:options error:&libraryError];
      if (library == nil) {
        std::fprintf(stderr, "SPHgrav Metal failed to compile shader source: %s\n",
                     [[[libraryError localizedDescription] description] UTF8String]);
        std::abort();
      }

      // one pipeline for each value of kHilo, so either can be used
      pipelineState_ = makePipeline(library, false);
      pipelineStateHilo_ = makePipeline(library, true);
    }
  }

  id<MTLComputePipelineState> makePipeline(id<MTLLibrary> library, bool hilo) {
    MTLFunctionConstantValues *constants = [[MTLFunctionConstantValues alloc] init];
    [constants setConstantValue:&hilo type:MTLDataTypeBool atIndex:0];

    NSError *functionError = nil;
    id<MTLFunction> function = [library newFunctionWithName:@"compute_forces"
                                             constantValues:constants
                                                      error:&functionError];
    if (function == nil) {
      std::fprintf(stderr, "SPHgrav Metal could not find compute_forces in the shader library: %s\n",
                   functionError != nil ? [[functionError localizedDescription] UTF8String] : "unknown error");
      std::abort();
    }

    NSError *pipelineError = nil;
    id<MTLComputePipelineState> state = [device_ newComputePipelineStateWithFunction:function error:&pipelineError];
    if (state == nil) {
      std::fprintf(stderr, "SPHgrav Metal failed to create compute pipeline: %s\n",
                   [[[pipelineError localizedDescription] description] UTF8String]);
      std::abort();
    }
    return state;
  }

  void resizeParticleBuffer(int count) {
    const NSUInteger byteCount = static_cast<NSUInteger>(count) * sizeof(MetalParticle);
    if (particleBuffer_ != nil && [particleBuffer_ length] >= byteCount) {
      return;
    }
    particleBuffer_ = [device_ newBufferWithLength:byteCount options:MTLResourceStorageModeShared];
    if (particleBuffer_ == nil) {
      std::fprintf(stderr, "SPHgrav Metal could not allocate particle buffer\n");
      std::abort();
    }
  }

  void resizeLoBuffer(int count) {
    const NSUInteger byteCount = static_cast<NSUInteger>(count) * sizeof(MetalLo);
    if (loBuffer_ != nil && [loBuffer_ length] >= byteCount) {
      return;
    }
    loBuffer_ = [device_ newBufferWithLength:byteCount options:MTLResourceStorageModeShared];
    if (loBuffer_ == nil) {
      std::fprintf(stderr, "SPHgrav Metal could not allocate position remainder buffer\n");
      std::abort();
    }
  }

  void resizeListBuffer(int count) {
    const NSUInteger byteCount = static_cast<NSUInteger>(count) * sizeof(uint32_t);
    if (listBuffer_ != nil && [listBuffer_ length] >= byteCount) {
      return;
    }
    listBuffer_ = [device_ newBufferWithLength:byteCount options:MTLResourceStorageModeShared];
    if (listBuffer_ == nil) {
      std::fprintf(stderr, "SPHgrav Metal could not allocate active-particle list buffer\n");
      std::abort();
    }
  }

  void resizeForceBuffer(int count) {
    const NSUInteger forceBytes = static_cast<NSUInteger>(count) * sizeof(MetalForce);
    if (forceBuffer_ == nil || [forceBuffer_ length] < forceBytes) {
      forceBuffer_ = [device_ newBufferWithLength:forceBytes options:MTLResourceStorageModeShared];
      if (forceBuffer_ == nil) {
        std::fprintf(stderr, "SPHgrav Metal could not allocate force buffer\n");
        std::abort();
      }
    }
    if (paramBuffer_ == nil) {
      paramBuffer_ = [device_ newBufferWithLength:sizeof(MetalParams) options:MTLResourceStorageModeShared];
      if (paramBuffer_ == nil) {
        std::fprintf(stderr, "SPHgrav Metal could not allocate parameter buffer\n");
        std::abort();
      }
    }
  }

  id<MTLDevice> device_ = nil;
  id<MTLCommandQueue> commandQueue_ = nil;
  id<MTLComputePipelineState> pipelineState_ = nil;
  id<MTLComputePipelineState> pipelineStateHilo_ = nil;
  id<MTLBuffer> particleBuffer_ = nil;
  id<MTLBuffer> forceBuffer_ = nil;
  id<MTLBuffer> paramBuffer_ = nil;
  id<MTLBuffer> listBuffer_ = nil;
  id<MTLBuffer> loBuffer_ = nil;
  id<MTLCommandBuffer> pendingCommandBuffer_ = nil;
  int pendingNi_ = 0;
  int pendingNibeg_ = 0;
  int pendingNrange_ = 0;
  bool pendingUseList_ = false;
  bool hilo_ = false;
  std::vector<int> pendingList_;
  std::chrono::steady_clock::time_point lastStart_{};
};

SPHgravMetal gMetalGrav;

}  // namespace

extern "C" {

void firsthalf_grav_forces_(int *n, int *n_lower, int *my_length,
                            double *px, double *py, double *pz,
                            double *mass, double *range2, int *q_, int *nkernel) {
  (void)q_;
  gMetalGrav.first_half(*n, *n_lower - 1, *my_length, px, py, pz, mass, range2, *nkernel);
}

void firsthalf_grav_forces_list_(int *n, int *n_lower, int *my_length, int *nlist, int *list,
                                 double *px, double *py, double *pz,
                                 double *mass, double *range2, int *nkernel) {
  gMetalGrav.first_half_list(*n, *n_lower - 1, *my_length, *nlist, list,
                             px, py, pz, mass, range2, *nkernel);
}

void lasthalf_grav_forces_(int *n, double *ax, double *ay, double *az, double *pot, int *myrank) {
  (void)n;
  (void)myrank;
  gMetalGrav.last_half(ax, ay, az, pot);
}

void gpu_init_dev_(int *myrank, double *theta) {
  (void)theta;
  gMetalGrav.setDevice(*myrank);
}

// ngravhilo in sph.input: 1 sends each position as a float plus its float
// remainder, 0 as a float alone
void gpu_set_hilo_(int *ngravhilo) {
  gMetalGrav.setHilo(*ngravhilo != 0);
}

}  // extern "C"
