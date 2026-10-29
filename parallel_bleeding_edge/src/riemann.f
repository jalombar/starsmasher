c     Riemann-solver dissipation (nav=4).  Each interacting pair of SPH
c     particles has a 1D Riemann problem along the line joining them,
c     between the two particles' own states.  The
c     difference between its star pressure and that of the same problem
c     with the two particles at rest takes the place of the artificial
c     viscosity (see balAV3.f).  The pressure force and the energy
c     equation are otherwise those of Gaburov et al. (2010), so with no
c     relative motion along the pair the equations of motion are the
c     variational ones.
************************************************************************
      subroutine riemann_alloc
c     allocate and initialize the nav=4 arrays, once
      use riemann_data
      include 'starsmasher.h'
      if(allocated(gam1p)) return
      allocate(gam1p(nmax))
      gam1p=gam
      return
      end
************************************************************************
      subroutine riemann_solve(rl,vl,pl,gl,rr,vr,pr,gr,nrs,ps)
c     Star-region pressure ps of the exact Riemann problem between a left
c     state (rl,vl,pl) and a right state (rr,vr,pr), velocities measured
c     along the left-to-right direction.  Each side is a gamma-law gas with
c     its own exponent, gl and gr (Toro 2009, Riemann Solvers and Numerical
c     Methods for Fluid Dynamics, ch. 4, whose pressure functions hold side
c     by side for different gammas).  For a general equation of state each
c     particle's Gamma_1 stands in for gamma.
c     nrs=1 (nriemann, the default): the first guess is that of Toro's
c     adaptive non-iterative solver (Toro 2009, sec. 9.5.2), the
c     linearized (PVRS) pressure when the two pressures are close and it
c     lies between them, otherwise the two-rarefaction (TRRS) or the
c     two-shock (TSRS) approximation.  One Newton step from there matches
c     the converged solution closely at about two thirds of the cost.
c     The step is taken in q=p**z, z=(gamma-1)/(2 gamma), not in p.
c     f_L+f_R+vr-vl rises with p and is concave in p, so a Newton step in p
c     from a guess above the root overshoots below it, and for a fast
c     approach across a large density contrast it can land on pmin, which
c     makes p*-pstat negative and cools the pair.  In q the rarefaction
c     branch is linear and the shock branch convex, so a Newton step never
c     lands below the root, and it is exact when both waves are
c     rarefactions (equal gammas).  The step is Halley's, Newton's times
c     1/(1-L) with L=g g''/(2 g'^2) in q, which removes most of Newton's
c     overestimate (about 8% of p*-p*_0 in a Mach 20 shock, 0.3% with
c     Halley's) at the cost of a few multiplications.  L is capped at 1/2.
c     nrs=0: Newton iteration to convergence from the PVRS guess, usually
c     two to four steps.
c     Both calls of a pair (p* and p*_0) take the same path, so ps is the
c     same for both when the pair has no relative motion.
      implicit none
      real*8 rl,vl,pl,gl,rr,vr,pr,gr,ps
      integer nrs
      real*8 cl,cr,pmin,p,pn,fl,fr,dl,dr,d2l,d2r
      real*8 pmn,pmx,ppv,gm,z,p0,gsl,gsr,q,gf,dpdq,gq,gqq,hl
      integer it

      cl=sqrt(gl*pl/rl)
      cr=sqrt(gr*pr/rr)
      pmin=1.d-12*min(pl,pr)
      if(nrs.eq.0) then
         p=max(0.5d0*(pl+pr)-0.125d0*(vr-vl)*(rl+rr)*(cl+cr),pmin)
         do it=1,30
            call riemann_f(p,rl,pl,cl,gl,fl,dl,d2l)
            call riemann_f(p,rr,pr,cr,gr,fr,dr,d2r)
            pn=max(p-(fl+fr+vr-vl)/(dl+dr),pmin)
            if(abs(pn-p).lt.1.d-10*(pn+p)) then
               ps=pn
               return
            endif
            p=pn
         enddo
         ps=p
         return
      endif
      pmn=min(pl,pr)
      pmx=max(pl,pr)
      ppv=max(0.5d0*(pl+pr)-0.125d0*(vr-vl)*(rl+rr)*(cl+cr),0.d0)
      if(pmx.le.2.d0*pmn .and. ppv.ge.pmn .and. ppv.le.pmx) then
         p=ppv
      elseif(ppv.lt.pmn) then
         gm=0.5d0*(gl+gr)
         z=(gm-1.d0)/(2.d0*gm)
         p=max(cl+cr-0.5d0*(gm-1.d0)*(vr-vl),0.d0)
     $        /(cl/pl**z+cr/pr**z)
         p=p**(1.d0/z)
      else
         p0=max(ppv,pmin)
         gsl=sqrt(2.d0/((gl+1.d0)*rl)/(p0+(gl-1.d0)/(gl+1.d0)*pl))
         gsr=sqrt(2.d0/((gr+1.d0)*rr)/(p0+(gr-1.d0)/(gr+1.d0)*pr))
         p=(gsl*pl+gsr*pr-(vr-vl))/(gsl+gsr)
      endif
      p=max(p,pmin)
c     one Halley step in q=p**z.  g=f_L+f_R+vr-vl and its derivatives in
c     q follow from those in p, with dp/dq=p/(z q) and
c     d2p/dq2=(1/z)(1/z-1) p/q**2
      call riemann_f(p,rl,pl,cl,gl,fl,dl,d2l)
      call riemann_f(p,rr,pr,cr,gr,fr,dr,d2r)
      gm=0.5d0*(gl+gr)
      z=(gm-1.d0)/(2.d0*gm)
      q=p**z
      dpdq=p/(z*q)
      gf=fl+fr+vr-vl
      gq=(dl+dr)*dpdq
      gqq=(d2l+d2r)*dpdq**2+(dl+dr)*(1.d0/z)*(1.d0/z-1.d0)*p/q**2
      hl=min(gf*gqq/(2.d0*gq**2),0.5d0)
      q=q-gf/(gq*(1.d0-hl))
      ps=max(q,pmin**z)**(1.d0/z)
      return
      end
************************************************************************
      subroutine riemann_f(p,rk,pk,ck,g,f,df,d2f)
c     Toro's pressure function f_K(p) and its first two derivatives, for
c     the state K
      implicit none
      real*8 p,rk,pk,ck,g,f,df,d2f,a,b,s
      if(p.gt.pk) then
c     shock
         a=2.d0/((g+1.d0)*rk)
         b=(g-1.d0)/(g+1.d0)*pk
         s=sqrt(a/(p+b))
         f=(p-pk)*s
         df=s*(1.d0-0.5d0*(p-pk)/(b+p))
         d2f=-0.5d0*df/(p+b)-0.5d0*s*(b+pk)/(b+p)**2
      else
c     rarefaction
         f=2.d0*ck/(g-1.d0)*((p/pk)**((g-1.d0)/(2.d0*g))-1.d0)
         df=1.d0/(rk*ck)*(p/pk)**(-(g+1.d0)/(2.d0*g))
         d2f=-(g+1.d0)/(2.d0*g)*df/p
      endif
      return
      end
