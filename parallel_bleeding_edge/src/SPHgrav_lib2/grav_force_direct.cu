#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <vector>
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

template<const int NTHREAD, const bool HILO>
__global__ void dev_compute_forces(
    const int nj,
    const int nibeg,
    const int ni,
    const devParticle *ptcl_in,
    const devParticleLo *plo_in,
    __gravout devForce    *force_out,
    const int nkernel,
    const int *ilist)
{
  const int tid = threadIdx.x;
  const int idx = blockIdx.x * blockDim.x + tid;

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
    for (int j = 0; j < nj; j += NTHREAD)
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
    for (int j = 0; j < nj; j += NTHREAD)
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
    for (int j = 0; j < nj; j += NTHREAD)
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
    force_out[idx] = iforce;
}

template<int NTHREAD, bool PINNED>
struct SPHgrav_direct
{
  cuVector<devParticle, PINNED> ptcl;
  cuVector<devParticleLo, PINNED> plo;
  bool hilo;
  cuVector<devForce,    PINNED> force;
  cuVector<int,         PINNED> ilist;

  // uselist: the last first_half computed forces only on the ni particles
  // in ilist (0-based), and last_half zeroes the rest of nibeg..nibeg+nrange-1
  bool uselist;
  int nrange;
  
  cudaEvent_t start, stop;

  int nibeg, ni;
  int ndevice;
  std::vector<bool> can_use_device;
  SPHgrav_direct(const int rank = 0) : hilo(false)
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

  void first_half(const int nibeg, const int ni, const int nkernel)
  {
    cudaEventCreate( &start );
    cudaEventCreate( &stop  );
    cudaEventRecord( start, 0 );

    this->nibeg = nibeg;
    this->ni    = ni;
    this->uselist = false;
    force.allocate(ni);
    ptcl.h2d();
    if (hilo) plo.h2d();


    const int nthreads =        NTHREAD;
    const int nblocks  = (ni-1)/NTHREAD + 1;

    const dim3 block(nthreads, 1, 1);
    const dim3 grid (nblocks,  1, 1);

    const int nj = ptcl.size;
    assert(nj % NTHREAD == 0);

    if (hilo)
      dev_compute_forces<NTHREAD, true ><<<grid, block>>>(nj, nibeg, ni, ptcl, plo, force, nkernel, (const int*)NULL);
    else
      dev_compute_forces<NTHREAD, false><<<grid, block>>>(nj, nibeg, ni, ptcl, (const devParticleLo*)NULL, force, nkernel, (const int*)NULL);
  }

  // forces on the ni particles listed (0-based) in ilist.host_pointer only
  void first_half_list(const int nibeg, const int nrange, const int ni, const int nkernel)
  {
    cudaEventCreate( &start );
    cudaEventCreate( &stop  );
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

    const int nthreads =        NTHREAD;
    const int nblocks  = (ni-1)/NTHREAD + 1;

    const dim3 block(nthreads, 1, 1);
    const dim3 grid (nblocks,  1, 1);

    const int nj = ptcl.size;
    assert(nj % NTHREAD == 0);

    if (hilo)
      dev_compute_forces<NTHREAD, true ><<<grid, block>>>(nj, nibeg, ni, ptcl, plo, force, nkernel, (const int*)ilist);
    else
      dev_compute_forces<NTHREAD, false><<<grid, block>>>(nj, nibeg, ni, ptcl, (const devParticleLo*)NULL, force, nkernel, (const int*)ilist);
  }

  float last_half()
  {
    if (ni > 0) force.d2h();
    cudaEventRecord( stop, 0 );
    cudaDeviceSynchronize();
    float elapsed_time_ms;
    cudaEventElapsedTime( &elapsed_time_ms, start, stop );

    return elapsed_time_ms/1000.0;
  }
};


namespace SPHgrav
{
  const int NTHREAD = 256;
  const bool PINNED = false;

  SPHgrav_direct<NTHREAD, PINNED> grav;

  void load_particles(
      const int nj,
      double *px, double *py, double *pz,
      double *mass, double *h2)
  {
    const float LARGE = 1.0e10;
    const int n = ((nj - 1)/NTHREAD + 1) * NTHREAD;
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
    double wsum = 0.0, xref = 0.0, yref = 0.0, zref = 0.0;
    for (int i = 0; i < nj; i++)
      if (h2[i] > 0.0)
      {
        const double w = mass[i];
        wsum += w;
        xref += w*px[i];
        yref += w*py[i];
        zref += w*pz[i];
      }
    if (wsum > 0.0)
    {
      xref /= wsum;
      yref /= wsum;
      zref /= wsum;
    }

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
      for (int j = grav.nibeg; j < grav.nibeg + grav.nrange; j++)
        ax[j] = ay[j] = az[j] = pot[j] = 0.0;

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

