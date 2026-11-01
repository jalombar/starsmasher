#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <cub/cub.cuh>
#include <thrust/iterator/counting_iterator.h>

// particles per thread in the kernel of the new order (1, 2 or 4)
#ifndef GRAV_NI_DEFAULT
#define GRAV_NI_DEFAULT 1
#endif
#include "cuVector.h"
#include "cuVec3.h"

#if 1      /* use this if you want the potential to be in double precision */
#define _GPOTDP_
#endif

struct ds64
{
  float2 val;
  __host__ __device__ ds64() {}
  __host__ __device__ ds64(float x) : val(make_float2(x, x)) {}
  __host__ __device__ ds64 operator+=(const float x) 
  {
    const float vx = val.x + x;
    const float vy = val.y - ((vx - val.x) - x);
    val = make_float2(vx, vy);
    return *this;
  }
  __host__ __device__ double to_double() const { return (double)val.x + (double)val.y; }
  __host__ __device__ float to_float() const { return (float)((double)val.x + (double)val.y);}
};

typedef cuVec3<float> fcuVec3;
struct devParticle
{
  fcuVec3 pos;                 // 3
  float   mass;                // 4
  float   h, invh, h2, hflag;  // 8
  __host__ __device__ devParticle() {}
  __host__ __device__ devParticle(const fcuVec3 _pos, const float _mass, const float _h, const float _flag) :
    pos(_pos), mass(_mass), h(_h), invh(1.0f/_h), h2(_h*_h), hflag(_flag) {}
};
#define PTCL_LEN (sizeof(devParticle) / sizeof(float4))

// ngravhilo=1: the float remainder of each position (position - pos), kept in
// an array of its own, so that ngravhilo=0 reads and computes what it did before
struct devParticleLo
{
  fcuVec3 lo;                  // 3
  float   pad;                 // 4
  __host__ __device__ devParticleLo() {}
  __host__ __device__ devParticleLo(const fcuVec3 _lo) : lo(_lo), pad(0.0f) {}
};

struct devForce
{
  ds64 ax, ay, az;   // 6
#ifdef _GPOTDP_
  ds64 pot;          // 8
#else
  float pot;         // 7
  int  iPad;        // 8
#endif
  __host__ __device__ devForce() {}
  __device__ devForce(const float v) : ax(v), ay(v), az(v), pot(v) {}
};

#define __gravout

template<const bool HILO>
__forceinline__ __device__ devForce dev_force_ij_nkernel0(
    const devParticle pi,
    const devParticle pj,
    const fcuVec3 ilo,
    const fcuVec3 jlo,
    __gravout devForce    iforce)
{
  const fcuVec3 dr = HILO ? (pj.pos - pi.pos) + (jlo - ilo) : pj.pos - pi.pos;
  const float   r2 = dr*dr;
  
  const float  rinv  = rsqrtf(r2);
  const float  rinv2 = rinv   * rinv;
  const float mrinv1 = rinv   * pj.mass;
  const float mrinv3 = rinv2  * mrinv1;


  if (r2 >= fmaxf(pi.h2, pj.h2))
  {
    iforce.ax  +=   mrinv3 * dr.x;
    iforce.ay  +=   mrinv3 * dr.y;
    iforce.az  +=   mrinv3 * dr.z;
    iforce.pot += (-mrinv1);
  } 
  else if (r2 < fminf(pi.h2, pj.h2))
  {
    const float2 invq  = r2 > 0.0f ? make_float2(rinv * pi.h, rinv * pj.h) : make_float2(0.0f, 0.0f); 
    const float2 q     = r2 > 0.0f ? make_float2(1.0f/invq.x, 1.0f/invq.y) : make_float2(0.0f, 0.0f); 
    const float2 q2    = {q.x*q.x, q.y*q.y};

    const float2 q3    = {q.x*q2.x, q.y*q2.y};
    const float2 invq2 = {invq.x*invq. x, invq.y*invq .y};
    const float2 invq3 = {invq.x*invq2.x, invq.y*invq2.y};

    const float2 f   = {q.x < 0.5f, q.y < 0.5f};
    const float2 acc = 
      {
        (       f.x) * (10.666666666667f + q2.x * (32.0f * q.x + (-38.4f))) + 
        (1.0f - f.x) * (21.333333333333f + (-48.0f)*q.x + 38.4f*q2.x + (-10.666666666667f)*q3.x + (-0.066666666667f) * invq3.x),
        (       f.y) * (10.666666666667f + q2.y * (32.0f * q.y + (-38.4f))) + 
        (1.0f - f.y) * (21.333333333333f + (-48.0f)*q.y + 38.4f*q2.y + (-10.666666666667f)*q3.y + (-0.066666666667f) * invq3.y)
      };
    
    const float2 pot =
      {
        (       f.x) * ((-2.8f) + q2.x * (5.333333333333f  + q2.x * (6.4f * q.x + (-9.6f)))) + 
        (1.0f - f.x) * ((-3.2f) + 0.066666666667f * invq.x + q2.x * (10.666666666667f + q.x * ((-16.0f) + q.x * (9.6f + (-2.133333333333f) * q.x)))),
        (       f.y) * ((-2.8f) + q2.y * (5.333333333333f  + q2.y * (6.4f * q.y + (-9.6f)))) + 
        (1.0f - f.y) * ((-3.2f) + 0.066666666667f * invq.y + q2.y * (10.666666666667f + q.y * ((-16.0f) + q.y * (9.6f + (-2.133333333333f) * q.y))))
      };
     
    const float2 mj1 = {pj.mass * pi.invh,         pj.mass * pj.invh };
    const float2 mj2 = {mj1.x   * pi.invh*pi.invh, mj1.y   * pj.invh*pj.invh};

    const float2 g = {r2 <= pi.h2, r2 <= pj.h2};
    const float gacc = r2 > 0.0f ? 
        0.5f*(g.x*mj2.x*acc.x + (1.0f - g.x)*mrinv3 + 
              g.y*mj2.y*acc.y + (1.0f - g.y)*mrinv3) : 0.0f;
    const float gpot = r2 > 0.0f || pi.h2 != pj.h2 || pi.mass != pj.mass ? 
        0.5f*(g.x*mj1.x*pot.x + (g.x - 1.0f)*mrinv1 + 
              g.y*mj1.y*pot.y + (g.y - 1.0f)*mrinv1) : (-1.4f)*pj.hflag*(mj1.x + mj1.y);

    iforce.ax  += gacc * dr.x;
    iforce.ay  += gacc * dr.y;
    iforce.az  += gacc * dr.z;
    iforce.pot += gpot;
  }
  else if (r2 < pi.h2)
  {
    const float invq  = r2 > 0.0f ? rinv * pi.h : 0.0f; 
    const float q     = r2 > 0.0f ? 1.0f/invq : 0.0f; 
    const float q2    = q*q;

    const float q3    = q*q2;
    const float invq2 = invq*invq;
    const float invq3 = invq*invq2;

    const float f   = q < 0.5f;
    const float acc = 
        (       f) * (10.666666666667f + q2 * (32.0f * q + (-38.4f))) + 
        (1.0f - f) * (21.333333333333f + (-48.0f)*q + 38.4f*q2 + (-10.666666666667f)*q3 + (-0.066666666667f) * invq3);
    
    const float pot =
        (       f) * ((-2.8f) + q2 * (5.333333333333f  + q2 * (6.4f * q + (-9.6f)))) + 
        (1.0f - f) * ((-3.2f) + 0.066666666667f * invq + q2 * (10.666666666667f + q * ((-16.0f) + q * (9.6f + (-2.133333333333f) * q))));
     
    const float2 mj1 = {pj.mass * pi.invh,         pj.mass * pj.invh };
    const float mj2 = mj1.x   * pi.invh*pi.invh;

/*    const float2 g = {1.0f, 0.0f}; */
    const float gacc = r2 > 0.0f ? 0.5f*(mj2*acc + mrinv3) : 0.0f;
    const float gpot = r2 > 0.0f ? 0.5f*(mj1.x*pot -mrinv1) : (-1.4f)*pj.hflag*(mj1.x + mj1.y);

    iforce.ax  += gacc * dr.x;
    iforce.ay  += gacc * dr.y;
    iforce.az  += gacc * dr.z;
    iforce.pot += gpot;
  }
  else
  {
    const float invq  = r2 > 0.0f ? rinv * pj.h : 0.0f; 
    const float q     = r2 > 0.0f ? 1.0f/invq : 0.0f; 
    const float q2    = q*q;

    const float q3    = q*q2;
    const float invq2 = invq*invq;
    const float invq3 = invq*invq2;

    const float f   = q < 0.5f;
    const float acc = 
        (       f) * (10.666666666667f + q2 * (32.0f * q + (-38.4f))) + 
        (1.0f - f) * (21.333333333333f + (-48.0f)*q + 38.4f*q2 + (-10.666666666667f)*q3 + (-0.066666666667f) * invq3);
    
    const float pot =
        (       f) * ((-2.8f) + q2 * (5.333333333333f  + q2 * (6.4f * q + (-9.6f)))) + 
        (1.0f - f) * ((-3.2f) + 0.066666666667f * invq + q2 * (10.666666666667f + q * ((-16.0f) + q * (9.6f + (-2.133333333333f) * q))));
     
    const float2 mj1 = {pj.mass * pi.invh,         pj.mass * pj.invh };
    const float mj2 = mj1.y   * pj.invh*pj.invh;

/*    const float2 g = {0.0f, 1.0f}; */
    const float gacc = r2 > 0.0f ? 0.5f*(mrinv3 + mj2*acc) : 0.0f;
    const float gpot = r2 > 0.0f ? 0.5f*(-mrinv1 + mj1.y*pot) : (-1.4f)*pj.hflag*(mj1.x + mj1.y);

    iforce.ax  += gacc * dr.x;
    iforce.ay  += gacc * dr.y;
    iforce.az  += gacc * dr.z;
    iforce.pot += gpot;
  }
  return iforce;
}

template<const bool HILO>
__forceinline__ __device__ devForce dev_force_ij_nkernel1(
    const devParticle pi,
    const devParticle pj,
    const fcuVec3 ilo,
    const fcuVec3 jlo,
    __gravout devForce    iforce)
{
  const fcuVec3 dr = HILO ? (pj.pos - pi.pos) + (jlo - ilo) : pj.pos - pi.pos;
  const float   r2 = dr*dr;
  
  const float  rinv  = rsqrtf(r2);
  const float  rinv2 = rinv   * rinv;
  const float mrinv1 = rinv   * pj.mass;
  const float mrinv3 = rinv2  * mrinv1;


  if (r2 >= fmaxf(pi.h2, pj.h2))
  {
    iforce.ax  +=   mrinv3 * dr.x;
    iforce.ay  +=   mrinv3 * dr.y;
    iforce.az  +=   mrinv3 * dr.z;
    iforce.pot += (-mrinv1);
  } 
  else if (r2 < fminf(pi.h2, pj.h2))
  {
    const float2 invq  = r2 > 0.0f ? make_float2(rinv * pi.h, rinv * pj.h) : make_float2(0.0f, 0.0f); 
    const float2 q     = r2 > 0.0f ? make_float2(1.0f/invq.x, 1.0f/invq.y) : make_float2(0.0f, 0.0f); 
    const float2 q2    = {q.x*q.x, q.y*q.y};

    const float2 acc = 
        {
        28.4375f + q2.x * ( (-187.6875f) + q2.x * (804.375f + q2.x * ( (-4379.375f) + q.x * (9009.0f + q.x * ( (-8957.8125f) + q.x * (5005.0f + q.x * ( (-1515.9375f) + 195.0f*q.x))))))) ,
        28.4375f + q2.y * ( (-187.6875f) + q2.y * (804.375f + q2.y * ( (-4379.375f) + q.y * (9009.0f + q.y * ( (-8957.8125f) + q.y * (5005.0f + q.y * ( (-1515.9375f) + 195.0f*q.y)))))))
        };
    
    const float2 pot =
        {
        r2 <= pi.h2 ? (-3.828125f) + q2.x * (14.21875f + q2.x * ((-46.921875f) + q2.x * (134.0625f + q2.x * ( (-547.421875f) + q.x * (1001.0f + q.x * ( (-895.78125f) + q.x * (455.0f + q.x * ((-126.328125f) + 15.0f*q.x)))))))) : 0.0f , 
        r2 <= pj.h2 ? (-3.828125f) + q2.y * (14.21875f + q2.y * ((-46.921875f) + q2.y * (134.0625f + q2.y * ( (-547.421875f) + q.y * (1001.0f + q.y * ( (-895.78125f) + q.y * (455.0f + q.y * ((-126.328125f) + 15.0f*q.y)))))))) : 0.0f
        };
     
     const float2 mj1 = {pj.mass * pi.invh,         pj.mass * pj.invh };
     const float2 mj2 = {mj1.x   * pi.invh*pi.invh, mj1.y   * pj.invh*pj.invh};

     const float2 g = {r2 <= pi.h2, r2 <= pj.h2};
     const float gacc = r2 > 0.0f ? 
      	       0.5f*(g.x*mj2.x*acc.x + (1.0f - g.x)*mrinv3 + 
                     g.y*mj2.y*acc.y + (1.0f - g.y)*mrinv3) : 0.0f;
     const float gpot = r2 > 0.0f || pi.h2 != pj.h2 || pi.mass != pj.mass ? 
     	       0.5f*(mj1.x*pot.x + (g.x - 1.0f)*mrinv1 + 
                     mj1.y*pot.y + (g.y - 1.0f)*mrinv1) : (-1.9140625f)*pj.hflag*(mj1.x + mj1.y);

    iforce.ax  += gacc * dr.x;
    iforce.ay  += gacc * dr.y;
    iforce.az  += gacc * dr.z;
    iforce.pot += gpot;
  }
  else if (r2 < pi.h2)
  {
    const float invq  = r2 > 0.0f ? rinv * pi.h : 0.0f; 
    const float q     = r2 > 0.0f ? 1.0f/invq : 0.0f; 
    const float q2    = q*q;

    const float acc = 
        28.4375f + q2 * ( (-187.6875f) + q2 * (804.375f + q2 * ( (-4379.375f) + q * (9009.0f + q * ( (-8957.8125f) + q * (5005.0f + q * ( (-1515.9375f) + 195.0f*q)))))));
    
    const float pot =
        (-3.828125f) + q2 * (14.21875f + q2 * ((-46.921875f) + q2 * (134.0625f + q2 * ( (-547.421875f) + q * (1001.0f + q * ( (-895.78125f) + q * (455.0f + q * ((-126.328125f) + 15.0f*q))))))));
     
     const float2 mj1 = {pj.mass * pi.invh,         pj.mass * pj.invh };
     const float mj2 = mj1.x   * pi.invh*pi.invh;

/*     const float2 g = {1.0f, 0.0f}; */
     const float gacc = r2 > 0.0f ? 0.5f*(mj2*acc + mrinv3) : 0.0f;
     const float gpot = r2 > 0.0f ? 0.5f*(mj1.x*pot - mrinv1) : (-1.9140625f)*pj.hflag*(mj1.x + mj1.y);

    iforce.ax  += gacc * dr.x;
    iforce.ay  += gacc * dr.y;
    iforce.az  += gacc * dr.z;
    iforce.pot += gpot;
  }
  else
  {
    const float invq  = r2 > 0.0f ? rinv * pj.h : 0.0f; 
    const float q     = r2 > 0.0f ? 1.0f/invq : 0.0f; 
    const float q2    = q*q;

    const float acc = 
        28.4375f + q2 * ( (-187.6875f) + q2 * (804.375f + q2 * ( (-4379.375f) + q * (9009.0f + q * ( (-8957.8125f) + q * (5005.0f + q * ( (-1515.9375f) + 195.0f*q)))))));
    
    const float pot =
        (-3.828125f) + q2 * (14.21875f + q2 * ((-46.921875f) + q2 * (134.0625f + q2 * ( (-547.421875f) + q * (1001.0f + q * ( (-895.78125f) + q * (455.0f + q * ((-126.328125f) + 15.0f*q))))))));
     
     const float2 mj1 = {pj.mass * pi.invh,         pj.mass * pj.invh };
     const float mj2 = mj1.y   * pj.invh*pj.invh;

/*     const float2 g = {0.0f, 1.0f}; */
     const float gacc = r2 > 0.0f ? 0.5f*(mrinv3 + mj2*acc) : 0.0f;
     const float gpot = r2 > 0.0f ? 0.5f*(-mrinv1 + mj1.y*pot) : (-1.9140625f)*pj.hflag*(mj1.x + mj1.y);

    iforce.ax  += gacc * dr.x;
    iforce.ay  += gacc * dr.y;
    iforce.az  += gacc * dr.z;
    iforce.pot += gpot;
  }
  return iforce;
}

template<const bool HILO>
__forceinline__ __device__ devForce dev_force_ij_nkernel2(
    const devParticle pi,
    const devParticle pj,
    const fcuVec3 ilo,
    const fcuVec3 jlo,
    __gravout devForce    iforce)
{
  const fcuVec3 dr = HILO ? (pj.pos - pi.pos) + (jlo - ilo) : pj.pos - pi.pos;
  const float   r2 = dr*dr;
  
  const float  rinv  = rsqrtf(r2);
  const float  rinv2 = rinv   * rinv;
  const float mrinv1 = rinv   * pj.mass;
  const float mrinv3 = rinv2  * mrinv1;


  if (r2 >= fmaxf(pi.h2, pj.h2))
  {
    iforce.ax  +=   mrinv3 * dr.x;
    iforce.ay  +=   mrinv3 * dr.y;
    iforce.az  +=   mrinv3 * dr.z;
    iforce.pot += (-mrinv1);
  } 
  else if (r2 < fminf(pi.h2, pj.h2))
  {
    const float2 invq  = r2 > 0.0f ? make_float2(rinv * pi.h, rinv * pj.h) : make_float2(0.0f, 0.0f); 
    const float2 q     = r2 > 0.0f ? make_float2(1.0f/invq.x, 1.0f/invq.y) : make_float2(0.0f, 0.0f); 
    const float2 q2    = {q.x*q.x, q.y*q.y};

    const float2 acc = 
    	{
	20.625f + q2.x * ( (-115.5f) + q2.x * (618.75f + q.x * ( (-1155.0f) + q.x * (962.5f + q.x * ( (-396.0f) + 65.625f * q.x))))),
	20.625f + q2.y * ( (-115.5f) + q2.y * (618.75f + q.y * ( (-1155.0f) + q.y * (962.5f + q.y * ( (-396.0f) + 65.625f * q.y)))))
    	};
 
    const float2 pot =
    	{
	(-3.4375f) + q2.x * (10.3125f + q2.x * ( (-28.875f)  + q2.x * (103.125f + q.x *( (-165.0f) + q.x * (120.3125f + q.x * ( (-44.0f) + 6.5625f*q.x)))))), 
	(-3.4375f) + q2.y * (10.3125f + q2.y * ( (-28.875f)  + q2.y * (103.125f + q.y *( (-165.0f) + q.y * (120.3125f + q.y * ( (-44.0f) + 6.5625f*q.y))))))
    	};
     
    const float2 mj1 = {pj.mass * pi.invh,         pj.mass * pj.invh };
    const float2 mj2 = {mj1.x   * pi.invh*pi.invh, mj1.y   * pj.invh*pj.invh};

/*    const float2 g = {1.0f, 1.0f}; */
    const float gacc = r2 > 0.0f ? 
      	    0.5f*(mj2.x*acc.x + mj2.y*acc.y) : 0.0f;
    const float gpot = r2 > 0.0f || pi.h2 != pj.h2 || pi.mass != pj.mass ? 
      	      0.5f*(mj1.x*pot.x + mj1.y*pot.y) : (-1.71875f)*pj.hflag*(mj1.x + mj1.y);

    iforce.ax  += gacc * dr.x;
    iforce.ay  += gacc * dr.y;
    iforce.az  += gacc * dr.z;
    iforce.pot += gpot;
  }
  else if (r2 < pi.h2)
  {
    const float invq  = r2 > 0.0f ? rinv * pi.h : 0.0f; 
    const float q     = r2 > 0.0f ? 1.0f/invq : 0.0f; 
    const float q2    = q*q;

    const float acc = 
	20.625f + q2 * ( (-115.5f) + q2 * (618.75f + q * ( (-1155.0f) + q * (962.5f + q * ( (-396.0f) + 65.625f * q)))));
 
    const float pot =
	(-3.4375f) + q2 * (10.3125f + q2 * ( (-28.875f)  + q2 * (103.125f + q *( (-165.0f) + q * (120.3125f + q * ( (-44.0f) + 6.5625f*q))))));
     
    const float2 mj1 = {pj.mass * pi.invh,         pj.mass * pj.invh };
    const float mj2 = mj1.x   * pi.invh*pi.invh;

/*    const float2 g = {1.0f, 0.0f}; */
    const float gacc = r2 > 0.0f ? 
      	    0.5f*(mj2*acc + mrinv3) : 0.0f;
    const float gpot = r2 > 0.0f ? 
      	      0.5f*(mj1.x*pot - mrinv1) : (-1.71875f)*pj.hflag*(mj1.x + mj1.y);

    iforce.ax  += gacc * dr.x;
    iforce.ay  += gacc * dr.y;
    iforce.az  += gacc * dr.z;
    iforce.pot += gpot;
  }
  else
  {
    const float invq  = r2 > 0.0f ? rinv * pj.h : 0.0f; 
    const float q     = r2 > 0.0f ? 1.0f/invq : 0.0f; 
    const float q2    = q*q;

    const float acc = 
	20.625f + q2 * ( (-115.5f) + q2 * (618.75f + q * ( (-1155.0f) + q * (962.5f + q * ( (-396.0f) + 65.625f * q)))));
 
    const float pot =
	(-3.4375f) + q2 * (10.3125f + q2 * ( (-28.875f)  + q2 * (103.125f + q *( (-165.0f) + q * (120.3125f + q * ( (-44.0f) + 6.5625f*q))))));
     
    const float2 mj1 = {pj.mass * pi.invh,         pj.mass * pj.invh };
    const float mj2 = mj1.y   * pj.invh*pj.invh;

/*    const float2 g = {0.0f, 1.0f}; */
    const float gacc = r2 > 0.0f ? 
      	    0.5f*(mrinv3 + mj2*acc) : 0.0f;
    const float gpot = r2 > 0.0f ? 
      	      0.5f*(- mrinv1 + mj1.y*pot) : (-1.71875f)*pj.hflag*(mj1.x + mj1.y);

    iforce.ax  += gacc * dr.x;
    iforce.ay  += gacc * dr.y;
    iforce.az  += gacc * dr.z;
    iforce.pot += gpot;
  }
  return iforce;
}

// The kernel took 52 or 55 registers per thread, which lets 4 blocks of 256
// threads share a multiprocessor.  Capped at 48 registers (a few bytes of
// spill) 5 fit, and on an RTX 4070 Ti SUPER it runs 8-10% faster: more warps
// to hide the latencies with.  6 blocks (42 registers) is slower again.
// Turing holds 1024 threads per multiprocessor, so there 4 is the most.
#if defined(__CUDA_ARCH__) && __CUDA_ARCH__ == 750
#define MIN_BLOCKS_PER_SM 4
#else
#define MIN_BLOCKS_PER_SM 5
#endif

template<const int NTHREAD, const bool HILO>
__global__ void __launch_bounds__(NTHREAD, MIN_BLOCKS_PER_SM) dev_compute_forces(
    const int nj,
    const int nibeg,
    const int ni,
    const devParticle *ptcl_in,
    const devParticleLo *plo_in,
    __gravout devForce    *force_out,
    const int nkernel,
    const int *ilist,
    const int jchunk)
{
  const int tid = threadIdx.x;
  const int idx = blockIdx.x * blockDim.x + tid;

  // blockIdx.y picks the block of jchunk consecutive j that this block sums
  // over, and force_out then holds one partial force for each of them,
  // blockIdx.y * ni + idx, to be added up by dev_sum_forces
  const int jbeg = blockIdx.y * jchunk;
  const int jend = min(nj, jbeg + jchunk);

  // ilist (block timesteps): force_out[idx] is for particle ilist[idx]
  // rather than for particle nibeg+idx
  const int iidx = ilist == NULL ? idx + nibeg : ilist[idx < ni ? idx : 0];
  const devParticle iptcl = ptcl_in[iidx];
  const fcuVec3 ilo = HILO ? plo_in[iidx].lo : fcuVec3(0.0f);
  devForce iforce(0.0f);

  __shared__ devParticle jpshared[NTHREAD];
  __shared__ fcuVec3 jloshared[HILO ? NTHREAD : 1];

  if(nkernel==0)
  {
    for (int j = jbeg; j < jend; j += NTHREAD)
    {
#if 1
      float4 *src = (float4*)&ptcl_in[j];
      float4 *dst = (float4*)jpshared;
#pragma unroll
      for (int it = 0; it < PTCL_LEN; it++)
      {
        dst[tid] = src[tid];
        dst += NTHREAD;
        src += NTHREAD;
      }
#else
      jpshared[tid] = ptcl_in[j + tid];
#endif
      if (HILO)
        jloshared[tid] = plo_in[j + tid].lo;
      __syncthreads();

      if (idx < ni)
      {
#pragma unroll 8
        for (int jj = 0; jj < NTHREAD; jj++)
          iforce = dev_force_ij_nkernel0<HILO>(iptcl, jpshared[jj], ilo, jloshared[HILO ? jj : 0], iforce);
      }
      __syncthreads();
    }
  }
  else if (nkernel==1)
  {
    for (int j = jbeg; j < jend; j += NTHREAD)
    {
#if 1
      float4 *src = (float4*)&ptcl_in[j];
      float4 *dst = (float4*)jpshared;
#pragma unroll
      for (int it = 0; it < PTCL_LEN; it++)
      {
        dst[tid] = src[tid];
        dst += NTHREAD;
        src += NTHREAD;
      }
#else
      jpshared[tid] = ptcl_in[j + tid];
#endif
      if (HILO)
        jloshared[tid] = plo_in[j + tid].lo;
      __syncthreads();

      if (idx < ni)
      {
#pragma unroll 8
        for (int jj = 0; jj < NTHREAD; jj++)
          iforce = dev_force_ij_nkernel1<HILO>(iptcl, jpshared[jj], ilo, jloshared[HILO ? jj : 0], iforce);
      }
      __syncthreads();
    }
  }
  else if (nkernel==2)
  {
    for (int j = jbeg; j < jend; j += NTHREAD)
    {
#if 1
      float4 *src = (float4*)&ptcl_in[j];
      float4 *dst = (float4*)jpshared;
#pragma unroll
      for (int it = 0; it < PTCL_LEN; it++)
      {
        dst[tid] = src[tid];
        dst += NTHREAD;
        src += NTHREAD;
      }
#else
      jpshared[tid] = ptcl_in[j + tid];
#endif
      if (HILO)
        jloshared[tid] = plo_in[j + tid].lo;
      __syncthreads();

      if (idx < ni)
      {
#pragma unroll 8
        for (int jj = 0; jj < NTHREAD; jj++)
          iforce = dev_force_ij_nkernel2<HILO>(iptcl, jpshared[jj], ilo, jloshared[HILO ? jj : 0], iforce);
      }
      __syncthreads();
    }
  }

  if (idx < ni)
    force_out[blockIdx.y * ni + idx] = iforce;
}

// Adds the nsplit partial forces of each particle in the order of the splits,
// in double-single arithmetic: the same result every time.
__global__ void dev_sum_forces(
    const int ni,
    const int nsplit,
    const devForce *part,
    __gravout devForce *force_out,
    const int *omap = NULL)   // omap (the new order): where the force of entry idx goes in force_out
{
  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx >= ni)
    return;
  devForce f = part[idx];
  for (int s = 1; s < nsplit; s++)
  {
    const devForce g = part[s * ni + idx];
    f.ax += g.ax.val.x;  f.ax += g.ax.val.y;
    f.ay += g.ay.val.x;  f.ay += g.ay.val.y;
    f.az += g.az.val.x;  f.az += g.az.val.y;
#ifdef _GPOTDP_
    f.pot += g.pot.val.x;  f.pot += g.pot.val.y;
#else
    f.pot += g.pot;
#endif
  }
  force_out[omap != NULL ? omap[idx] : idx] = f;
}

// ------------------------------------------------------------------------------------------------------------------
// Spatial order and the fast path for tiles that are far from the block.
//
// The particles are put in the order of a Morton key (sort_by_key on the device), so that the 256 particles of a
// block of i and the 256 of a tile of j are close to each other.  For every tile of j the box of its particles and
// the largest h^2 in it are found; the block of i has its own box.  If the squared distance between the two boxes
// is at least the larger of the two h^2, every pair of the tile is farther apart than the softening of either
// particle, whatever the kernel (r^2 >= max(h_i^2, h_j^2) is the first test of dev_force_ij_nkernel*), and its
// force is Newtonian.  The decision is the same for all the threads of the block.  The Newtonian pair is computed
// with the positions of the j taken relative to the centre of the block, in the shared memory, so that a pair
// costs three subtractions instead of nine when the float remainders of ngravhilo=1 are used; the remainders are
// added in once per tile and particle, not once per pair.  The terms are added in float over 64 j and the sum is
// then added to the double-single force of the particle (before, every pair went through the double-single
// addition, 16 operations).
// Tiles that are near the block, and the one with the padding, go through the original code, on the particles in
// the new order.
// ------------------------------------------------------------------------------------------------------------------
struct TileBox
{
  float4 lo;   // min x, y, z of the tile, max h^2
  float4 hi;   // max x, y, z, and 1 if the tile holds padding (never taken for far)
};

__device__ __forceinline__ unsigned long long spread3(unsigned int v)
{
  unsigned long long x = v & 0x1fffffu;
  x = (x | x << 32) & 0x1f00000000ffffULL;
  x = (x | x << 16) & 0x1f0000ff0000ffULL;
  x = (x | x << 8)  & 0x100f00f00f00f00fULL;
  x = (x | x << 4)  & 0x10c30c30c30c30c3ULL;
  x = (x | x << 2)  & 0x1249249249249249ULL;
  return x;
}

__global__ void dev_morton_keys(const int nreal, const devParticle *ptcl, const float bx, const float by, const float bz,
    const float inv, unsigned long long *keys, int *vals)
{
  const int k = blockIdx.x * blockDim.x + threadIdx.x;
  if (k >= nreal)
    return;
  const fcuVec3 p = ptcl[k].pos;
  const unsigned int ix = (unsigned int)fminf(fmaxf((p.x - bx) * inv, 0.0f), 2097151.0f);
  const unsigned int iy = (unsigned int)fminf(fmaxf((p.y - by) * inv, 0.0f), 2097151.0f);
  const unsigned int iz = (unsigned int)fminf(fmaxf((p.z - bz) * inv, 0.0f), 2097151.0f);
  keys[k] = spread3(ix) | (spread3(iy) << 1) | (spread3(iz) << 2);
  vals[k] = k;
}

// sorted copy: particle k of the new order is particle perm[k] of the old one; the padding stays at the end
__global__ void dev_gather_sorted(const int n, const int nreal, const int *perm, const devParticle *ptcl, const devParticleLo *plo,
    devParticle *sptcl, devParticleLo *splo)
{
  const int k = blockIdx.x * blockDim.x + threadIdx.x;
  if (k >= n)
    return;
  const int o = k < nreal ? perm[k] : k;
  sptcl[k] = ptcl[o];
  if (plo != NULL)
    splo[k] = plo[o];
}

template<int NTHREAD>
__global__ void dev_tile_boxes(const int nreal, const devParticle *sp, TileBox *tb)
{
  const int j = blockIdx.x * NTHREAD + threadIdx.x;
  const bool real = j < nreal;
  float lx = 1.0e30f, ly = 1.0e30f, lz = 1.0e30f, hx = -1.0e30f, hy = -1.0e30f, hz = -1.0e30f, hh = 0.0f, pd = real ? 0.0f : 1.0f;
  if (real)
  {
    const devParticle p = sp[j];
    lx = hx = p.pos.x;  ly = hy = p.pos.y;  lz = hz = p.pos.z;  hh = p.h2;
  }
  for (int o = 16; o > 0; o >>= 1)
  {
    lx = fminf(lx, __shfl_xor_sync(0xffffffffu, lx, o));  hx = fmaxf(hx, __shfl_xor_sync(0xffffffffu, hx, o));
    ly = fminf(ly, __shfl_xor_sync(0xffffffffu, ly, o));  hy = fmaxf(hy, __shfl_xor_sync(0xffffffffu, hy, o));
    lz = fminf(lz, __shfl_xor_sync(0xffffffffu, lz, o));  hz = fmaxf(hz, __shfl_xor_sync(0xffffffffu, hz, o));
    hh = fmaxf(hh, __shfl_xor_sync(0xffffffffu, hh, o));  pd = fmaxf(pd, __shfl_xor_sync(0xffffffffu, pd, o));
  }
  __shared__ float red[NTHREAD / 32][8];
  const int w = threadIdx.x / 32;
  if ((threadIdx.x & 31) == 0)
  {
    red[w][0] = lx; red[w][1] = ly; red[w][2] = lz; red[w][3] = hx; red[w][4] = hy; red[w][5] = hz; red[w][6] = hh; red[w][7] = pd;
  }
  __syncthreads();
  if (threadIdx.x == 0)
  {
    for (int q = 1; q < NTHREAD / 32; q++)
    {
      lx = fminf(lx, red[q][0]); ly = fminf(ly, red[q][1]); lz = fminf(lz, red[q][2]);
      hx = fmaxf(hx, red[q][3]); hy = fmaxf(hy, red[q][4]); hz = fmaxf(hz, red[q][5]);
      hh = fmaxf(hh, red[q][6]); pd = fmaxf(pd, red[q][7]);
    }
    // (the values of warp 0 are the ones this thread already holds)
    lx = fminf(lx, red[0][0]); ly = fminf(ly, red[0][1]); lz = fminf(lz, red[0][2]);
    hx = fmaxf(hx, red[0][3]); hy = fmaxf(hy, red[0][4]); hz = fmaxf(hz, red[0][5]);
    hh = fmaxf(hh, red[0][6]); pd = fmaxf(pd, red[0][7]);
    TileBox b;
    b.lo = make_float4(lx, ly, lz, hh);
    b.hi = make_float4(hx, hy, hz, pd);
    tb[blockIdx.x] = b;
  }
}

// the particles (in the new order) whose original index lies in [lo, hi): the i of this process
struct InOriginalRange
{
  const int *perm;
  int lo, hi;
  __device__ bool operator()(const int k) const { const int o = perm[k]; return o >= lo && o < hi; }
};

__global__ void dev_make_omap(const int ni, const int nibeg, const int *isort, const int *perm, int *omap)
{
  const int e = blockIdx.x * blockDim.x + threadIdx.x;
  if (e < ni)
    omap[e] = perm[isort[e]] - nibeg;
}

// 1/sqrt(x) without the test for a denormal argument that rsqrtf makes (three instructions more per pair);
// in a far tile r^2 is at least h^2, which is not tiny
__device__ __forceinline__ float rsqrt_ftz(const float x)
{
  float y;
  asm("rsqrt.approx.ftz.f32 %0, %1;" : "=f"(y) : "f"(x));
  return y;
}

template<const int NK, const bool HILO>
__forceinline__ __device__ devForce dev_force_ij(const devParticle pi, const devParticle pj, const fcuVec3 ilo, const fcuVec3 jlo, devForce f)
{
  if (NK == 0)
    return dev_force_ij_nkernel0<HILO>(pi, pj, ilo, jlo, f);
  else if (NK == 1)
    return dev_force_ij_nkernel1<HILO>(pi, pj, ilo, jlo, f);
  return dev_force_ij_nkernel2<HILO>(pi, pj, ilo, jlo, f);
}

// blocks per multiprocessor asked of the compiler for the kernel of the new order, j between two additions to the double-single
// force, and the unrolling of the loop over j (tuned on Blackwell sm_120 at 5*10^5 particles)
#ifndef MINB_S
#define MINB_S 3
#endif
#ifndef FAR_FLUSH
#define FAR_FLUSH 128
#endif
#ifndef FAR_UNROLL
#define FAR_UNROLL 8
#endif
#define DO_PRAGMA(x) _Pragma(#x)
#define UNROLL_N(n) DO_PRAGMA(unroll n)

// Forces on the ni particles whose positions in the new order are isort[0..ni): the block handles NTHREAD*NI of them,
// the thread NI (the ones tid, tid+NTHREAD, ...).  force_out holds one partial force per j-chunk (blockIdx.y).
template<const int NTHREAD, const int NI, const bool HILO, const int NK>
__global__ void __launch_bounds__(NTHREAD, MINB_S) dev_compute_forces_s(
    const int nj,
    const int ni,
    const devParticle *sp,
    const devParticleLo *slo,
    const int *isort,
    const TileBox *tbox,
    const TileBox *gbox,
    const int usefast,
    const float pad,
    __gravout devForce *force_out,
    const int jchunk,
    unsigned char *mask,
    const int mode,
    const int ntiles)
{
  const int tid = threadIdx.x;
  const int ib = blockIdx.x * (NTHREAD * NI);
  const int jbeg = blockIdx.y * jchunk;
  const int jend = min(nj, jbeg + jchunk);

  __shared__ devParticle jpshared[NTHREAD];
  __shared__ fcuVec3 jloshared[HILO ? NTHREAD : 1];
  __shared__ float4 jrel[NTHREAD];
  __shared__ float bred[NTHREAD / 32][8];
  __shared__ float bbox[8];

  int si[NI];
  bool val[NI];
  devForce iforce[NI];
  fcuVec3 ipos[NI], ilo[NI];
  float lx = 1.0e30f, ly = 1.0e30f, lz = 1.0e30f, hx = -1.0e30f, hy = -1.0e30f, hz = -1.0e30f, hh = 0.0f;
#pragma unroll
  for (int k = 0; k < NI; k++)
  {
    const int e = ib + tid + k * NTHREAD;
    val[k] = e < ni;
    si[k] = isort[val[k] ? e : 0];
    iforce[k] = devForce(0.0f);
    const devParticle ip = sp[si[k]];
    ipos[k] = ip.pos;
    ilo[k] = HILO ? slo[si[k]].lo : fcuVec3(0.0f);
    if (val[k])
    {
      lx = fminf(lx, ip.pos.x);  hx = fmaxf(hx, ip.pos.x);
      ly = fminf(ly, ip.pos.y);  hy = fmaxf(hy, ip.pos.y);
      lz = fminf(lz, ip.pos.z);  hz = fmaxf(hz, ip.pos.z);
      hh = fmaxf(hh, ip.h2);
    }
  }
  // the box of the block of i, and (all lanes of a warp holding the same values after the exchange) that of the warp
  for (int o = 16; o > 0; o >>= 1)
  {
    lx = fminf(lx, __shfl_xor_sync(0xffffffffu, lx, o));  hx = fmaxf(hx, __shfl_xor_sync(0xffffffffu, hx, o));
    ly = fminf(ly, __shfl_xor_sync(0xffffffffu, ly, o));  hy = fmaxf(hy, __shfl_xor_sync(0xffffffffu, hy, o));
    lz = fminf(lz, __shfl_xor_sync(0xffffffffu, lz, o));  hz = fmaxf(hz, __shfl_xor_sync(0xffffffffu, hz, o));
    hh = fmaxf(hh, __shfl_xor_sync(0xffffffffu, hh, o));
  }
  const float wlx = lx, wly = ly, wlz = lz, whx = hx, why = hy, whz = hz, wh2 = hh;
  if ((tid & 31) == 0)
  {
    bred[tid / 32][0] = lx; bred[tid / 32][1] = ly; bred[tid / 32][2] = lz;
    bred[tid / 32][3] = hx; bred[tid / 32][4] = hy; bred[tid / 32][5] = hz; bred[tid / 32][6] = hh;
  }
  __syncthreads();
  if (tid == 0)
  {
    lx = bred[0][0]; ly = bred[0][1]; lz = bred[0][2]; hx = bred[0][3]; hy = bred[0][4]; hz = bred[0][5]; hh = bred[0][6];
    for (int q = 1; q < NTHREAD / 32; q++)
    {
      lx = fminf(lx, bred[q][0]); ly = fminf(ly, bred[q][1]); lz = fminf(lz, bred[q][2]);
      hx = fmaxf(hx, bred[q][3]); hy = fmaxf(hy, bred[q][4]); hz = fmaxf(hz, bred[q][5]);
      hh = fmaxf(hh, bred[q][6]);
    }
    bbox[0] = lx; bbox[1] = ly; bbox[2] = lz; bbox[3] = hx; bbox[4] = hy; bbox[5] = hz; bbox[6] = hh;
  }
  __syncthreads();
  const float blx = bbox[0], bly = bbox[1], blz = bbox[2], bhx = bbox[3], bhy = bbox[4], bhz = bbox[5], bh2 = bbox[6];
  // positions of the i relative to the centre of the block (the remainder added in)
  const float cx = 0.5f * (blx + bhx), cy = 0.5f * (bly + bhy), cz = 0.5f * (blz + bhz);
  float irx[NI], iry[NI], irz[NI];
#pragma unroll
  for (int k = 0; k < NI; k++)
  {
    irx[k] = (ipos[k].x - cx) + ilo[k].x;
    iry[k] = (ipos[k].y - cy) + ilo[k].y;
    irz[k] = (ipos[k].z - cz) + ilo[k].z;
  }

  for (int j0 = jbeg; j0 < jend; j0 += NTHREAD)
  {
    // mode 2: the tiles that the far pass (mode 1) has done are not done again
    if (mode == 2 && mask[(size_t)blockIdx.x * ntiles + j0 / NTHREAD] != 0)
      continue;
    bool far = false;
    if (usefast)
    {
      const TileBox tb = tbox[j0 / NTHREAD];
      const float gx = fmaxf(fmaxf(fmaxf(tb.lo.x - bhx, blx - tb.hi.x), 0.0f) - pad, 0.0f);
      const float gy = fmaxf(fmaxf(fmaxf(tb.lo.y - bhy, bly - tb.hi.y), 0.0f) - pad, 0.0f);
      const float gz = fmaxf(fmaxf(fmaxf(tb.lo.z - bhz, blz - tb.hi.z), 0.0f) - pad, 0.0f);
      const float sep2 = gx * gx + gy * gy + gz * gz;
      far = tb.hi.w == 0.0f && sep2 >= 1.00002f * fmaxf(bh2, tb.lo.w);
    }
    // mode 1: only the far tiles, and a note of which they are for the near pass
    if (mode == 1)
    {
      if (tid == 0)
        mask[(size_t)blockIdx.x * ntiles + j0 / NTHREAD] = far ? 1 : 0;
      if (!far)
        continue;
    }
    if (far)
    {
      const float4 q = *(const float4*)&sp[j0 + tid];
      const fcuVec3 lo = HILO ? slo[j0 + tid].lo : fcuVec3(0.0f);
      jrel[tid] = make_float4((q.x - cx) + lo.x, (q.y - cy) + lo.y, (q.z - cz) + lo.z, q.w);
      __syncthreads();
      for (int jb = 0; jb < NTHREAD; jb += FAR_FLUSH)
      {
        float pax[NI], pay[NI], paz[NI], pp[NI];
#pragma unroll
        for (int k = 0; k < NI; k++)
          pax[k] = pay[k] = paz[k] = pp[k] = 0.0f;
UNROLL_N(FAR_UNROLL)
        for (int jj = jb; jj < jb + FAR_FLUSH; jj++)
        {
          const float4 qj = jrel[jj];
#pragma unroll
          for (int k = 0; k < NI; k++)
          {
            const float dx = qj.x - irx[k];
            const float dy = qj.y - iry[k];
            const float dz = qj.z - irz[k];
            const float r2 = dx * dx + dy * dy + dz * dz;
            const float rinv = rsqrt_ftz(r2);
            const float mrinv1 = rinv * qj.w;
            const float mrinv3 = (rinv * rinv) * mrinv1;
            pax[k] += mrinv3 * dx;
            pay[k] += mrinv3 * dy;
            paz[k] += mrinv3 * dz;
            pp[k] -= mrinv1;
          }
        }
#pragma unroll
        for (int k = 0; k < NI; k++)
        {
          iforce[k].ax += pax[k];
          iforce[k].ay += pay[k];
          iforce[k].az += paz[k];
          iforce[k].pot += pp[k];
        }
      }
      __syncthreads();
    }
    else
    {
      // The tile is near the block somewhere.  Its particles are taken in groups of 32, each of which the warp compares with
      // its own 32 i: the groups that are far from them take the Newtonian loop, the others the full pair function.
      float4 *src = (float4*)&sp[j0];
      float4 *dst = (float4*)jpshared;
#pragma unroll
      for (int it = 0; it < PTCL_LEN; it++)
      {
        dst[tid] = src[tid];
        dst += NTHREAD;
        src += NTHREAD;
      }
      if (HILO)
        jloshared[tid] = slo[j0 + tid].lo;
      {
        const float4 q = *(const float4*)&sp[j0 + tid];
        const fcuVec3 lo = HILO ? slo[j0 + tid].lo : fcuVec3(0.0f);
        jrel[tid] = make_float4((q.x - cx) + lo.x, (q.y - cy) + lo.y, (q.z - cz) + lo.z, q.w);
      }
      __syncthreads();
      float pax[NI], pay[NI], paz[NI], pp[NI];
#pragma unroll
      for (int k = 0; k < NI; k++)
        pax[k] = pay[k] = paz[k] = pp[k] = 0.0f;
      for (int g = 0; g < NTHREAD / 32; g++)
      {
        bool gfar = false;
        if (usefast)
        {
          const TileBox gb = gbox[j0 / 32 + g];
          const float gx = fmaxf(fmaxf(fmaxf(gb.lo.x - whx, wlx - gb.hi.x), 0.0f) - pad, 0.0f);
          const float gy = fmaxf(fmaxf(fmaxf(gb.lo.y - why, wly - gb.hi.y), 0.0f) - pad, 0.0f);
          const float gz = fmaxf(fmaxf(fmaxf(gb.lo.z - whz, wlz - gb.hi.z), 0.0f) - pad, 0.0f);
          gfar = gb.hi.w == 0.0f && gx * gx + gy * gy + gz * gz >= 1.00002f * fmaxf(wh2, gb.lo.w);
        }
        if (gfar)
        {
#pragma unroll 8
          for (int jj = g * 32; jj < g * 32 + 32; jj++)
          {
            const float4 qj = jrel[jj];
#pragma unroll
            for (int k = 0; k < NI; k++)
            {
              const float dx = qj.x - irx[k];
              const float dy = qj.y - iry[k];
              const float dz = qj.z - irz[k];
              const float r2 = dx * dx + dy * dy + dz * dz;
              const float rinv = rsqrt_ftz(r2);
              const float mrinv1 = rinv * qj.w;
              const float mrinv3 = (rinv * rinv) * mrinv1;
              pax[k] += mrinv3 * dx;
              pay[k] += mrinv3 * dy;
              paz[k] += mrinv3 * dz;
              pp[k] -= mrinv1;
            }
          }
        }
        else
        {
#pragma unroll
          for (int k = 0; k < NI; k++)
          {
            if (val[k])
            {
              const devParticle iptcl = sp[si[k]];
              const fcuVec3 il = HILO ? slo[si[k]].lo : fcuVec3(0.0f);
              devForce f = iforce[k];
#pragma unroll 8
              for (int jj = g * 32; jj < g * 32 + 32; jj++)
                f = dev_force_ij<NK, HILO>(iptcl, jpshared[jj], il, jloshared[HILO ? jj : 0], f);
              iforce[k] = f;
            }
          }
        }
      }
#pragma unroll
      for (int k = 0; k < NI; k++)
      {
        iforce[k].ax += pax[k];
        iforce[k].ay += pay[k];
        iforce[k].az += paz[k];
        iforce[k].pot += pp[k];
      }
      __syncthreads();
    }
  }
#pragma unroll
  for (int k = 0; k < NI; k++)
    if (val[k])
      force_out[blockIdx.y * ni + ib + tid + k * NTHREAD] = iforce[k];
}

// as dev_sum_forces for the splits of two passes: the nsplit of the far pass (parta) and then those of the near pass (partb), in this order
__global__ void dev_sum_forces2(
    const int ni,
    const int nsplit,
    const devForce *parta,
    const devForce *partb,
    __gravout devForce *force_out,
    const int *omap)
{
  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx >= ni)
    return;
  devForce f = parta[idx];
  for (int s = 1; s < 2 * nsplit; s++)
  {
    const devForce g = s < nsplit ? parta[s * ni + idx] : partb[(s - nsplit) * ni + idx];
    f.ax += g.ax.val.x;  f.ax += g.ax.val.y;
    f.ay += g.ay.val.x;  f.ay += g.ay.val.y;
    f.az += g.az.val.x;  f.az += g.az.val.y;
#ifdef _GPOTDP_
    f.pot += g.pot.val.x;  f.pot += g.pot.val.y;
#else
    f.pot += g.pot;
#endif
  }
  force_out[omap != NULL ? omap[idx] : idx] = f;
}

template<int NTHREAD, bool PINNED>
struct SPHgrav_direct
{
  cuVector<devParticle, PINNED> ptcl;
  cuVector<devParticleLo, PINNED> plo;
  bool hilo;
  cuVector<devForce,    PINNED> force;
  cuVector<int,         PINNED> ilist;

  // partial forces of the splits in j (see launch), on the device only
  devForce *part;
  size_t partcap;

  // uselist: the last first_half computed forces only on the ni particles
  // in ilist (0-based), and last_half zeroes the rest of nibeg..nibeg+nrange-1
  bool uselist;
  int nrange;

  // Spatial order and fast path (see dev_compute_forces_s).  nreal is the number of particles (ptcl.size is rounded
  // up to a multiple of NTHREAD), bbmin/bbmax the box of their positions.
  int nreal;
  float bbmin[3], bbmax[3];
  devParticle *sptcl;
  devParticleLo *splo;
  unsigned long long *keys_a, *keys_b;
  int *perm_a, *perm_b, *isort, *omap, *dnum;
  TileBox *tbox, *gbox;
  // the far and near passes (ngravoverlap): the second sorted copy, boxes and partial forces of the near pass, the record of the far tiles
  // (one byte per block of i and tile of j), and what the far pass chose
  devParticle *sptcl2;
  TileBox *tbox2, *gbox2;
  devForce *part2;
  unsigned char *dmask;
  size_t part2cap, dmaskcap;
  int s2cap, t2cap, g2cap;
  int lmode, nsplit_far;
  void *ctemp;
  size_t ctempcap;
  int scap, slocap, tcap, gcap, icap, kcap, pcap, ocap, dcap;
  
  // the events that time one call; created by the first call and kept
  cudaEvent_t start, stop;
  bool events_made;

  int nibeg, ni;
  int ndevice;
  std::vector<bool> can_use_device;
  SPHgrav_direct(const int rank = 0) : hilo(false), part(NULL), partcap(0), nreal(0),
    sptcl(NULL), splo(NULL), keys_a(NULL), keys_b(NULL), perm_a(NULL), perm_b(NULL), isort(NULL), omap(NULL), dnum(NULL), tbox(NULL), gbox(NULL),
    sptcl2(NULL), tbox2(NULL), gbox2(NULL), part2(NULL), dmask(NULL), part2cap(0), dmaskcap(0), s2cap(0), t2cap(0), g2cap(0), lmode(0), nsplit_far(0),
    ctemp(NULL), ctempcap(0), scap(0), slocap(0), tcap(0), gcap(0), icap(0), kcap(0), pcap(0), ocap(0), dcap(0), events_made(false)
  {
    assert(cudaGetDeviceCount(&ndevice) == 0);
    if (rank == 0)
      fprintf(stderr, " SPHgrav found %d CUDA devices \n", ndevice);
    assert(ndevice > 0);
    can_use_device = std::vector<bool>(ndevice, false);
    int no_supported = 0;
    for (int dev = 0; dev < ndevice; dev++)
    {
      cudaDeviceProp p;
      assert(cudaGetDeviceProperties(&p, dev) == cudaSuccess);
      int computeMode = -1;
      cudaDeviceGetAttribute(&computeMode, cudaDevAttrComputeMode, dev);
      const bool supported = p.major > 1 || (p.major == 1 && p.minor >= 3);
      if (rank == 0)
        fprintf(stderr,"  Device= %d: %s computeMode= %d computeCapability= %d_%d  supported= %s\n", 
            dev, p.name, computeMode, p.major, p.minor, supported ? "YES" : "NO");
      no_supported += supported ? 1 : 0;
      can_use_device[dev] = true;
    }
    assert(no_supported > 0);
  }

  void setDevice(const int device)
  {
    /* Report this rather than assert on it: an assertion here produces a bare
       stack trace that gives no hint that the cause is a setting in sph.input,
       and assertions vanish entirely if NDEBUG is ever defined. */
    if (device < 0 || device >= ndevice)
    {
      fprintf(stderr,
          "\n"
          "StarSmasher: a process asked for GPU %d, but this node has only %d CUDA device%s.\n"
          "  This almost always means ngravprocs in sph.input is larger than the number\n"
          "  of GPUs available. Set ngravprocs to at most %d, or to -%d to request that\n"
          "  many GPUs per node, or leave it at 0 to have it detected automatically.\n"
          "\n",
          device, ndevice, ndevice == 1 ? "" : "s", ndevice, ndevice);
      fflush(stderr);
      exit(EXIT_FAILURE);
    }
    const cudaError_t err = cudaSetDevice(device);
    if (err != cudaSuccess)
    {
      fprintf(stderr, "\nStarSmasher: cudaSetDevice(%d) failed: %s\n\n",
              device, cudaGetErrorString(err));
      fflush(stderr);
      exit(EXIT_FAILURE);
    }
  }

  // The grid has one block of NTHREAD particles i along x, and the sum over
  // j is split in nsplit pieces along y.  With the whole sum in each block
  // the 391 blocks of a 10^5-particle run need 1.48 times the blocks that
  // fit on the device at once, so the second round leaves half of it idle.
  // Splitting makes the rounds many and short.  nsplit is the one for which
  // the rounds cost the least, taking the smallest of those within 1%.
  template<bool HI>
  int choose_split(const int nblocks, const int nj)
  {
    int dev = 0, nsm = 0, nb = 0;
    cudaGetDevice(&dev);
    cudaDeviceGetAttribute(&nsm, cudaDevAttrMultiProcessorCount, dev);
    cudaOccupancyMaxActiveBlocksPerMultiprocessor(&nb,
        dev_compute_forces<NTHREAD, HI>, NTHREAD, 0);
    const int resident = nsm * nb;
    // pieces of at least 8 tiles of NTHREAD, so that a block's loop stays long
    const int smax = std::max(1, std::min(32, nj / NTHREAD / 8));
    if (resident <= 0)
      return 1;
    double best = 1.0e30;
    for (int s = 1; s <= smax; s++)
      best = std::min(best, (double)((nblocks * s + resident - 1) / resident) / s);
    for (int s = 1; s <= smax; s++)
      if ((double)((nblocks * s + resident - 1) / resident) / s <= 1.01 * best)
        return s;
    return 1;
  }

  // the kernel of the new order for NI particles per thread
  template<int NI, bool HI, int NK>
  void launch_s_one(const int ni, const int nsplit, const int jchunk, devForce *out, const dim3 grid,
      const devParticle *spp, const TileBox *tb, const TileBox *gb, unsigned char *mask, const int mode, const int ntiles)
  {
    const int nj = ptcl.size;
    const float pad = 4.0e-7f * std::max(std::max(std::max(fabsf(bbmin[0]), fabsf(bbmax[0])), std::max(fabsf(bbmin[1]), fabsf(bbmax[1]))),
                                         std::max(fabsf(bbmin[2]), fabsf(bbmax[2])));
    dev_compute_forces_s<NTHREAD, NI, HI, NK><<<grid, dim3(NTHREAD, 1, 1)>>>(nj, ni, spp, HI ? splo : (const devParticleLo*)NULL,
        isort, tb, gb, 1, pad, out, jchunk, mask, mode, ntiles);
    CUDA_SAFE_CALL(cudaGetLastError());
  }

  template<int NI, bool HI>
  int choose_split_s(const int ni, const int nj)
  {
    int dev = 0, nsm = 0, nb = 0;
    cudaGetDevice(&dev);
    cudaDeviceGetAttribute(&nsm, cudaDevAttrMultiProcessorCount, dev);
    cudaOccupancyMaxActiveBlocksPerMultiprocessor(&nb, dev_compute_forces_s<NTHREAD, NI, HI, 2>, NTHREAD, 0);
    const int nblocks = (ni - 1) / (NTHREAD * NI) + 1;
    const int resident = nsm * nb;
    const int smax = std::max(1, std::min(32, nj / NTHREAD / 8));
    if (resident <= 0)
      return 1;
    double best = 1.0e30;
    for (int s = 1; s <= smax; s++)
      best = std::min(best, (double)((nblocks * s + resident - 1) / resident) / s);
    for (int s = 1; s <= smax; s++)
      if ((double)((nblocks * s + resident - 1) / resident) / s <= 1.01 * best)
        return s;
    return 1;
  }

  template<int NI, bool HI>
  void launch_s_ni(const int ni, const int nkernel)
  {
    const int nj = ptcl.size;
    const int nblocks = (ni - 1) / (NTHREAD * NI) + 1;
    int nsplit = choose_split_s<NI, HI>(ni, nj);
    const int ntiles = nj / NTHREAD;
    const int tiles = (ntiles + nsplit - 1) / nsplit;
    nsplit = (ntiles + tiles - 1) / tiles;
    const int jchunk = tiles * NTHREAD;
    // the near pass has the same split as the far pass (it depends on ni and nj only), and adds its partial forces to those of the far pass
    assert(lmode != 2 || nsplit == nsplit_far);
    const size_t need = (size_t)nsplit * ni;
    if (need > partcap)
    {
      if (part != NULL)
        CUDA_SAFE_CALL(cudaFree(part));
      CUDA_SAFE_CALL(cudaMalloc((void**)&part, need * sizeof(devForce)));
      partcap = need;
    }
    const dim3 grid(nblocks, nsplit, 1);
    const int nt = ptcl.size / NTHREAD;
    if (lmode == 0)
    {
      if (nkernel == 0)      launch_s_one<NI, HI, 0>(ni, nsplit, jchunk, part, grid, sptcl, tbox, gbox, NULL, 0, nt);
      else if (nkernel == 1) launch_s_one<NI, HI, 1>(ni, nsplit, jchunk, part, grid, sptcl, tbox, gbox, NULL, 0, nt);
      else                   launch_s_one<NI, HI, 2>(ni, nsplit, jchunk, part, grid, sptcl, tbox, gbox, NULL, 0, nt);
      dev_sum_forces<<<(ni - 1) / 256 + 1, 256>>>(ni, nsplit, part, force, omap);
      CUDA_SAFE_CALL(cudaGetLastError());
    }
    else if (lmode == 1)
    {
      nsplit_far = nsplit;
      const size_t mneed = (size_t)nblocks * nt;
      if (mneed > dmaskcap)
      {
        if (dmask != NULL)
          CUDA_SAFE_CALL(cudaFree(dmask));
        CUDA_SAFE_CALL(cudaMalloc((void**)&dmask, mneed));
        dmaskcap = mneed;
      }
      if (nkernel == 0)      launch_s_one<NI, HI, 0>(ni, nsplit, jchunk, part, grid, sptcl, tbox, gbox, dmask, 1, nt);
      else if (nkernel == 1) launch_s_one<NI, HI, 1>(ni, nsplit, jchunk, part, grid, sptcl, tbox, gbox, dmask, 1, nt);
      else                   launch_s_one<NI, HI, 2>(ni, nsplit, jchunk, part, grid, sptcl, tbox, gbox, dmask, 1, nt);
    }
    else
    {
      if (need > part2cap)
      {
        if (part2 != NULL)
          CUDA_SAFE_CALL(cudaFree(part2));
        CUDA_SAFE_CALL(cudaMalloc((void**)&part2, need * sizeof(devForce)));
        part2cap = need;
      }
      if (nkernel == 0)      launch_s_one<NI, HI, 0>(ni, nsplit, jchunk, part2, grid, sptcl2, tbox2, gbox2, dmask, 2, nt);
      else if (nkernel == 1) launch_s_one<NI, HI, 1>(ni, nsplit, jchunk, part2, grid, sptcl2, tbox2, gbox2, dmask, 2, nt);
      else                   launch_s_one<NI, HI, 2>(ni, nsplit, jchunk, part2, grid, sptcl2, tbox2, gbox2, dmask, 2, nt);
      dev_sum_forces2<<<(ni - 1) / 256 + 1, 256>>>(ni, nsplit, part, part2, force, omap);
      CUDA_SAFE_CALL(cudaGetLastError());
    }
  }
  template<typename T>
  void ensure(T **p, int &cap, const int need)
  {
    if (need <= cap)
      return;
    if (*p != NULL)
      CUDA_SAFE_CALL(cudaFree(*p));
    CUDA_SAFE_CALL(cudaMalloc((void**)p, (size_t)need * sizeof(T)));
    cap = need;
  }

  template<typename T>
  void ensure_pair(T **p, T **q, int &cap, const int need)
  {
    if (need <= cap)
      return;
    if (*p != NULL)
    {
      CUDA_SAFE_CALL(cudaFree(*p));
      CUDA_SAFE_CALL(cudaFree(*q));
    }
    CUDA_SAFE_CALL(cudaMalloc((void**)p, (size_t)need * sizeof(T)));
    CUDA_SAFE_CALL(cudaMalloc((void**)q, (size_t)need * sizeof(T)));
    cap = need;
  }

  // Put the particles in the order of a Morton key, find the box of every tile, and list the i of this process
  // (those whose original index is nibeg..nibeg+ni-1) in the new order.  All on the device.
  void prepare_sorted(const int ni)
  {
    const int n = ptcl.size;
    ensure(&sptcl, scap, n);
    if (hilo)
      ensure(&splo, slocap, n);
    // keys_a/keys_b and perm_a/perm_b are allocated together, so they share a capacity
    ensure_pair(&keys_a, &keys_b, kcap, nreal);
    ensure_pair(&perm_a, &perm_b, pcap, nreal);
    ensure(&isort, icap, std::max(ni, 1));
    ensure(&omap, ocap, std::max(ni, 1));
    ensure(&dnum, dcap, 1);
    const int ntiles = n / NTHREAD;
    ensure(&tbox, tcap, ntiles);
    ensure(&gbox, gcap, n / 32);
    float ext = std::max(std::max(bbmax[0] - bbmin[0], bbmax[1] - bbmin[1]), bbmax[2] - bbmin[2]);
    const float inv = 2097151.0f / std::max(ext, 1.0e-30f);
    dev_morton_keys<<<(nreal - 1) / 256 + 1, 256>>>(nreal, ptcl, bbmin[0], bbmin[1], bbmin[2], inv, keys_a, perm_a);
    size_t bytes = 0;
    cub::DeviceRadixSort::SortPairs(NULL, bytes, keys_a, keys_b, perm_a, perm_b, nreal, 0, 63);
    // the temporary space serves the selection too
    size_t bytes2 = 0;
    thrust::counting_iterator<int> counting(0);
    cub::DeviceSelect::If(NULL, bytes2, counting, isort, dnum, nreal, InOriginalRange{perm_b, nibeg, nibeg + ni});
    bytes = std::max(bytes, bytes2);
    if (bytes > ctempcap)
    {
      if (ctemp != NULL)
        CUDA_SAFE_CALL(cudaFree(ctemp));
      CUDA_SAFE_CALL(cudaMalloc(&ctemp, bytes));
      ctempcap = bytes;
    }
    cub::DeviceRadixSort::SortPairs(ctemp, bytes, keys_a, keys_b, perm_a, perm_b, nreal, 0, 63);
    dev_gather_sorted<<<(n - 1) / 256 + 1, 256>>>(n, nreal, perm_b, ptcl, hilo ? (const devParticleLo*)plo : (const devParticleLo*)NULL, sptcl, splo);
    dev_tile_boxes<NTHREAD><<<ntiles, NTHREAD>>>(nreal, sptcl, tbox);
    dev_tile_boxes<32><<<n / 32, 32>>>(nreal, sptcl, gbox);
    cub::DeviceSelect::If(ctemp, bytes, counting, isort, dnum, nreal, InOriginalRange{perm_b, nibeg, nibeg + ni});
    dev_make_omap<<<(ni - 1) / 256 + 1, 256>>>(ni, nibeg, isort, perm_b, omap);
    CUDA_SAFE_CALL(cudaGetLastError());
  }

  void launch_sorted(const int ni, const int nkernel)
  {
    prepare_sorted(ni);
    launch_sorted_only(ni, nkernel);
  }

  void launch_sorted_only(const int ni, const int nkernel)
  {
    if (hilo) launch_s_ni<GRAV_NI_DEFAULT, true>(ni, nkernel); else launch_s_ni<GRAV_NI_DEFAULT, false>(ni, nkernel);
  }

  void launch(const int ni, const int nkernel, const int *ilist_dev)
  {
    if (ilist_dev == NULL)
    {
      launch_sorted(ni, nkernel);
      return;
    }
    const int nj = ptcl.size;
    const int nblocks = (ni - 1) / NTHREAD + 1;
    int nsplit = hilo ? choose_split<true>(nblocks, nj) : choose_split<false>(nblocks, nj);
    const int ntiles = nj / NTHREAD;
    const int tiles = (ntiles + nsplit - 1) / nsplit;
    nsplit = (ntiles + tiles - 1) / tiles;
    const int jchunk = tiles * NTHREAD;

    devForce *out = force;
    if (nsplit > 1)
    {
      const size_t need = (size_t)nsplit * ni;
      if (need > partcap)
      {
        if (part != NULL)
          CUDA_SAFE_CALL(cudaFree(part));
        CUDA_SAFE_CALL(cudaMalloc((void**)&part, need * sizeof(devForce)));
        partcap = need;
      }
      out = part;
    }

    const dim3 block(NTHREAD, 1, 1);
    const dim3 grid (nblocks, nsplit, 1);
    if (hilo)
      dev_compute_forces<NTHREAD, true ><<<grid, block>>>(nj, nibeg, ni, ptcl, plo, out, nkernel, ilist_dev, jchunk);
    else
      dev_compute_forces<NTHREAD, false><<<grid, block>>>(nj, nibeg, ni, ptcl, (const devParticleLo*)NULL, out, nkernel, ilist_dev, jchunk);
    CUDA_SAFE_CALL(cudaGetLastError());
    if (nsplit > 1)
      dev_sum_forces<<<(ni - 1) / 256 + 1, 256>>>(ni, nsplit, part, force);
    CUDA_SAFE_CALL(cudaGetLastError());
  }

  void make_events()
  {
    if (events_made)
      return;
    CUDA_SAFE_CALL(cudaEventCreate(&start));
    CUDA_SAFE_CALL(cudaEventCreate(&stop));
    events_made = true;
  }

  void first_half(const int nibeg, const int ni, const int nkernel)
  {
    make_events();
    cudaEventRecord( start, 0 );

    this->nibeg = nibeg;
    this->ni    = ni;
    this->uselist = false;
    force.allocate(ni);
    ptcl.h2d();
    if (hilo) plo.h2d();


    assert(ptcl.size % NTHREAD == 0);
    launch(ni, nkernel, (const int*)NULL);
  }

  // ngravoverlap.  The far pass (called before the smoothing lengths are solved, with a limit for each h): the sorted order, the boxes of
  // the tiles with the limits, and the kernel for the tiles that are far whatever the h, which it notes in dmask.  The near pass (with the
  // h found): the second sorted copy, with the new h, the kernel for the rest of the tiles, and the sum of the partial forces of the two
  // passes (the far ones first).  Everything is queued on the default stream, so the host does not wait for the far pass.
  void first_half_far(const int nibeg, const int ni, const int nkernel)
  {
    make_events();
    cudaEventRecord( start, 0 );
    this->nibeg = nibeg;
    this->ni    = ni;
    this->uselist = false;
    force.allocate(ni);
    ptcl.h2d();
    if (hilo) plo.h2d();
    assert(ptcl.size % NTHREAD == 0);
    prepare_sorted(ni);
    lmode = 1;
    launch_sorted_only(ni, nkernel);
    lmode = 0;
  }

  void first_half_near(const int nkernel)
  {
    const int n = ptcl.size;
    ensure(&sptcl2, s2cap, n);
    ensure(&tbox2, t2cap, n / NTHREAD);
    ensure(&gbox2, g2cap, n / 32);
    ptcl.h2d_async();
    dev_gather_sorted<<<(n - 1) / 256 + 1, 256>>>(n, nreal, perm_b, ptcl, (const devParticleLo*)NULL, sptcl2, splo);
    dev_tile_boxes<NTHREAD><<<n / NTHREAD, NTHREAD>>>(nreal, sptcl2, tbox2);
    dev_tile_boxes<32><<<n / 32, 32>>>(nreal, sptcl2, gbox2);
    CUDA_SAFE_CALL(cudaGetLastError());
    lmode = 2;
    launch_sorted_only(ni, nkernel);
    lmode = 0;
  }

  // the far pass is not wanted after all (a smoothing length went past its limit): wait for it, and the usual first_half follows
  void abort_far()
  {
    CUDA_SAFE_CALL(cudaDeviceSynchronize());
  }

  // forces on the ni particles listed (0-based) in ilist.host_pointer only
  void first_half_list(const int nibeg, const int nrange, const int ni, const int nkernel)
  {
    make_events();
    cudaEventRecord( start, 0 );

    this->nibeg   = nibeg;
    this->nrange  = nrange;
    this->ni      = ni;
    this->uselist = true;
    if (ni == 0) return;
    force.allocate(ni);
    ptcl.h2d();
    if (hilo) plo.h2d();
    ilist.h2d();

    assert(ptcl.size % NTHREAD == 0);
    launch(ni, nkernel, (const int*)ilist);
  }

  float last_half()
  {
    if (ni > 0) force.d2h();
    cudaEventRecord( stop, 0 );
    CUDA_SAFE_CALL(cudaDeviceSynchronize());
    float elapsed_time_ms;
    cudaEventElapsedTime( &elapsed_time_ms, start, stop );

    return elapsed_time_ms/1000.0;
  }
};


#ifndef GRAV_NTHREAD
#define GRAV_NTHREAD 256
#endif
namespace SPHgrav
{
  const int NTHREAD = GRAV_NTHREAD;
  const bool PINNED = true;

  SPHgrav_direct<NTHREAD, PINNED> grav;

  void load_particles(
      const int nj,
      double *px, double *py, double *pz,
      double *mass, double *h2)
  {
    const float LARGE = 1.0e10;
    const int n = ((nj - 1)/NTHREAD + 1) * NTHREAD;
    grav.nreal = nj;
    grav.ptcl.allocate(n);
    if (grav.hilo) grav.plo.allocate(n);

    // The kernel works in single precision, so positions are sent relative
    // to a reference point, with the subtraction done in double.  Otherwise
    // gas far from the origin, as for a star near a supermassive black hole,
    // loses its close separations to rounding.  The reference is the gas
    // centre of mass (h2 > 0 marks gas, so a black hole does not pull it
    // away from the star).  A star's rounding grows with its distance from
    // the reference, so this favours the more massive of two stars, and sits
    // midway between two identical ones.  Only separations enter the forces
    // and potential.
    // With ngravhilo=1 each coordinate also carries its float
    // remainder (lo), and the kernel forms dr = (pj.pos - pi.pos) + (pj.lo -
    // pi.lo): nearby floats subtract exactly, so a separation keeps float
    // precision relative to itself wherever the pair is, and the reference
    // point no longer matters.  Without it, two separated stars cannot both
    // be near the reference point, and the farther one's close separations
    // lose digits.
    // The sums are made over blocks of 4096 particles, added in the order of the blocks, so that the
    // reference point does not depend on how many threads do the blocks (the code is built with
    // OpenMP when the rest of StarSmasher is).
    const int BLK = 4096;
    const int nblk = (nj + BLK - 1) / BLK;
    std::vector<double> part(10 * (size_t)nblk);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int b = 0; b < nblk; b++)
    {
      double bw = 0.0, bx = 0.0, by = 0.0, bz = 0.0;
      double mn[3] = {1.0e300, 1.0e300, 1.0e300}, mx[3] = {-1.0e300, -1.0e300, -1.0e300};
      const int i1 = std::min(nj, (b + 1) * BLK);
      for (int i = b * BLK; i < i1; i++)
      {
        if (h2[i] > 0.0)
        {
          const double w = mass[i];
          bw += w;
          bx += w*px[i];
          by += w*py[i];
          bz += w*pz[i];
        }
        mn[0] = std::min(mn[0], px[i]);  mx[0] = std::max(mx[0], px[i]);
        mn[1] = std::min(mn[1], py[i]);  mx[1] = std::max(mx[1], py[i]);
        mn[2] = std::min(mn[2], pz[i]);  mx[2] = std::max(mx[2], pz[i]);
      }
      part[10*b] = bw;
      part[10*b+1] = bx;
      part[10*b+2] = by;
      part[10*b+3] = bz;
      for (int c = 0; c < 3; c++)
      {
        part[10*b+4+c] = mn[c];
        part[10*b+7+c] = mx[c];
      }
    }
    double wsum = 0.0, xref = 0.0, yref = 0.0, zref = 0.0;
    double bmn[3] = {1.0e300, 1.0e300, 1.0e300}, bmx[3] = {-1.0e300, -1.0e300, -1.0e300};
    for (int b = 0; b < nblk; b++)
    {
      wsum += part[10*b];
      xref += part[10*b+1];
      yref += part[10*b+2];
      zref += part[10*b+3];
      for (int c = 0; c < 3; c++)
      {
        bmn[c] = std::min(bmn[c], part[10*b+4+c]);
        bmx[c] = std::max(bmx[c], part[10*b+7+c]);
      }
    }
    if (wsum > 0.0)
    {
      xref /= wsum;
      yref /= wsum;
      zref /= wsum;
    }
    // the box of the positions as the device sees them (relative to the reference point)
    grav.bbmin[0] = (float)(bmn[0] - xref);  grav.bbmax[0] = (float)(bmx[0] - xref);
    grav.bbmin[1] = (float)(bmn[1] - yref);  grav.bbmax[1] = (float)(bmx[1] - yref);
    grav.bbmin[2] = (float)(bmn[2] - zref);  grav.bbmax[2] = (float)(bmx[2] - zref);

    // each particle is converted on its own: no thread writes another's element
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int i = 0; i < nj; i++)
    {
      assert(h2  [i] != 0.0);
      assert(mass[i] >  0.0);
      const double dx = px[i] - xref, dy = py[i] - yref, dz = pz[i] - zref;
      grav.ptcl[i] = devParticle(fcuVec3(dx, dy, dz), mass[i],
          std::sqrt(h2[i] > 0.0 ? h2[i] : -h2[i]),
          h2[i] > 0.0 ? 1.0 : 0.0);
      if (grav.hilo)
        grav.plo[i] = devParticleLo(fcuVec3(dx - (double)grav.ptcl[i].pos.x,
                                            dy - (double)grav.ptcl[i].pos.y,
                                            dz - (double)grav.ptcl[i].pos.z));
    }

    for (int i = nj; i < n; i++)
    {
      grav.ptcl[i] = devParticle(fcuVec3(LARGE), 0.0f, 1.0f, 0.0f);
      if (grav.hilo) grav.plo[i] = devParticleLo(fcuVec3(0.0f));
    }
  }

  void firsthalf_grav_force(
      const int nj, const int nibeg, const int ni,
      double *px, double *py, double *pz,
      double *mass, double *h2, const int nkernel)
  {
    assert(nibeg >= 0);
    load_particles(nj, px, py, pz, mass, h2);
    grav.first_half(nibeg, ni, nkernel);
  }

  // ngravoverlap: the far pass, before the smoothing lengths are solved (h2 is a limit for each of them), and the near pass with the h found.
  // ok is 0 if the far pass could not be started.
  void firsthalf_grav_force_far(
      const int nj, const int nibeg, const int ni,
      double *px, double *py, double *pz,
      double *mass, double *h2, const int nkernel, int *ok)
  {
    assert(nibeg >= 0);
    load_particles(nj, px, py, pz, mass, h2);
    grav.first_half_far(nibeg, ni, nkernel);
    *ok = 1;
  }

  void firsthalf_grav_force_near(
      const int nj, double *px, double *py, double *pz,
      double *mass, double *h2, const int nkernel)
  {
    load_particles(nj, px, py, pz, mass, h2);
    grav.first_half_near(nkernel);
  }

  void gravity_far_abort()
  {
    grav.abort_far();
  }

  // block timesteps: forces only on the nlist particles in list (1-based,
  // all within nibeg+1..nibeg+nrange), out of all nj
  void firsthalf_grav_force_list(
      const int nj, const int nibeg, const int nrange,
      const int nlist, const int *list,
      double *px, double *py, double *pz,
      double *mass, double *h2, const int nkernel)
  {
    load_particles(nj, px, py, pz, mass, h2);
    if (nlist > 0)
    {
      grav.ilist.allocate(nlist);
      for (int i = 0; i < nlist; i++)
        grav.ilist[i] = list[i] - 1;
    }
    grav.first_half_list(nibeg, nrange, nlist, nkernel);
  }

  void lasthalf_grav_force(
      const int ni0,
      double *ax, double *ay, double *az, double *pot)
  {
    //fprintf(stderr, " -- 0: ni= %d  ni= %d  grav.ni= %d \n", ni, grav.ni);
    // assert(ni == grav.ni);
    const int ni = grav.ni;
    const float dt = grav.last_half();

    if (grav.uselist)
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
      for (int j = grav.nibeg; j < grav.nibeg + grav.nrange; j++)
        ax[j] = ay[j] = az[j] = pot[j] = 0.0;
    }

    // the ni results go to distinct elements (the entries of ilist are distinct)
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int i = 0; i < ni; i++)
    {
      const int j = grav.uselist ? grav.ilist[i] : i + grav.nibeg;
      ax [j] = grav.force[i].ax.to_double();
      ay [j] = grav.force[i].ay.to_double();
      az [j] = grav.force[i].az.to_double();
#ifdef _GPOTDP_
      pot[j] = grav.force[i].pot.to_double();
#else
      pot[j] = grav.force[i].pot;
#endif
    }


#if 1
    fprintf(stderr, " >>> SPHgrav_lib took %g sec\n", dt);
#endif
  }

  void setDevice(const int device)
  {
    grav.setDevice(device);
  }

  void setHilo(const bool hilo)
  {
    grav.hilo = hilo;
  }
}

extern "C"
{
  void firsthalf_grav_forces_(
      int *n, int *n_lower, int *my_length, 
      double *px, double *py, double *pz,
      double *mass, double *range2, int *q_, int *nkernel)
  {
    SPHgrav::firsthalf_grav_force(*n, *n_lower - 1, *my_length, px, py, pz, mass, range2, *nkernel);
  }

  void firsthalf_grav_forces_far_(
      int *n, int *n_lower, int *my_length,
      double *px, double *py, double *pz,
      double *mass, double *range2, int *q_, int *nkernel, int *ok)
  {
    SPHgrav::firsthalf_grav_force_far(*n, *n_lower - 1, *my_length, px, py, pz, mass, range2, *nkernel, ok);
  }

  void firsthalf_grav_forces_near_(
      int *n, double *px, double *py, double *pz,
      double *mass, double *range2, int *nkernel)
  {
    SPHgrav::firsthalf_grav_force_near(*n, px, py, pz, mass, range2, *nkernel);
  }

  void gravity_far_abort_()
  {
    SPHgrav::gravity_far_abort();
  }

  void firsthalf_grav_forces_list_(
      int *n, int *n_lower, int *my_length, int *nlist, int *list,
      double *px, double *py, double *pz,
      double *mass, double *range2, int *nkernel)
  {
    SPHgrav::firsthalf_grav_force_list(*n, *n_lower - 1, *my_length, *nlist, list,
        px, py, pz, mass, range2, *nkernel);
  }

  void lasthalf_grav_forces_(
      int *n, 
      double *ax, double *ay, double *az,
      double *pot,int *myrank) 
  {
    SPHgrav::lasthalf_grav_force(*n, ax, ay, az, pot);
  }

#if 0
  void gpu_init_dev_(int *myrank)
  {
    SPHgrav::setDevice(*myrank);
  }
#else 
  void gpu_init_dev_(int *myrank, double *theta)
  {
    SPHgrav::setDevice(*myrank);
  }
#endif

  // ngravhilo in sph.input: 1 sends each position as a float plus its float
  // remainder, 0 as a float alone
  void gpu_set_hilo_(int *ngravhilo)
  {
    SPHgrav::setHilo(*ngravhilo != 0);
  }
}

