c     Riemann-solver dissipation (nav=4).  Each interacting pair of SPH
c     particles has a 1D Riemann problem along the line joining them, with
c     the two particles' states reconstructed to the pair midpoint.  The
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
      allocate(gam1p(nmax),rgrad(15,nmax))
      gam1p=gam
      rgrad=0.d0
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
c     branch is linear and the shock branch convex, so the step never
c     lands below the root, and it is exact when both waves are
c     rarefactions (equal gammas).
c     nrs=0: Newton iteration to convergence from the PVRS guess, usually
c     two to four steps.
c     Both calls of a pair (p* and p*_0) take the same path, so ps is the
c     same for both when the pair has no relative motion.
      implicit none
      real*8 rl,vl,pl,gl,rr,vr,pr,gr,ps
      integer nrs
      real*8 cl,cr,pmin,p,pn,fl,fr,dl,dr
      real*8 pmn,pmx,ppv,gm,z,p0,gsl,gsr,q
      integer it

      cl=sqrt(gl*pl/rl)
      cr=sqrt(gr*pr/rr)
      pmin=1.d-12*min(pl,pr)
      if(nrs.eq.0) then
         p=max(0.5d0*(pl+pr)-0.125d0*(vr-vl)*(rl+rr)*(cl+cr),pmin)
         do it=1,30
            call riemann_f(p,rl,pl,cl,gl,fl,dl)
            call riemann_f(p,rr,pr,cr,gr,fr,dr)
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
c     one Newton step in q=p**z, with dq/dp=z q/p
      call riemann_f(p,rl,pl,cl,gl,fl,dl)
      call riemann_f(p,rr,pr,cr,gr,fr,dr)
      gm=0.5d0*(gl+gr)
      z=(gm-1.d0)/(2.d0*gm)
      q=p**z
      q=q-(fl+fr+vr-vl)*z*q/((dl+dr)*p)
      ps=max(q,pmin**z)**(1.d0/z)
      return
      end
************************************************************************
      subroutine riemann_f(p,rk,pk,ck,g,f,df)
c     Toro's pressure function f_K(p) and its derivative, for the state K
      implicit none
      real*8 p,rk,pk,ck,g,f,df,a,b,s
      if(p.gt.pk) then
c     shock
         a=2.d0/((g+1.d0)*rk)
         b=(g-1.d0)/(g+1.d0)*pk
         s=sqrt(a/(p+b))
         f=(p-pk)*s
         df=s*(1.d0-0.5d0*(p-pk)/(b+p))
      else
c     rarefaction
         f=2.d0*ck/(g-1.d0)*((p/pk)**((g-1.d0)/(2.d0*g))-1.d0)
         df=1.d0/(rk*ck)*(p/pk)**(-(g+1.d0)/(2.d0*g))
      endif
      return
      end
************************************************************************
      subroutine riemann_gradients
c     SPH estimates of grad rho, grad p and the velocity gradient tensor,
c         grad q_i = sum_j (m_j/rho_j) (q_j - q_i) grad_i W(r_ij, h_i),
c     over the neighbour lists of the force loop, for the particles this
c     rank solves.  They are then shared, since the pair loop needs the
c     gradients of every neighbour.  With block timesteps only the
c     particles solved in this substep are updated, and the others keep
c     the gradients of their last solve, as they keep their pressures.
      use riemann_data
      include 'starsmasher.h'
      integer i,j,in,k,itab,iw,nil
      real*8 dx,dy,dz,r2,h2,invh5,dw,f,pres_i,pres_j,d(3)
      real*8 buf(nmax)
      save buf

      if(blkdist) then
         nil=nwmine
      else
         nil=n_upper-n_lower+1
      endif
      do iw=1,nil
         if(blkdist) then
            i=wmine(iw)
         else
            i=n_lower+iw-1
            if(nblock.eq.1 .and. .not.(actblk(i).or.refblk(i))) cycle
         endif
         do k=1,15
            rgrad(k,i)=0.d0
         enddo
         if(u(i).eq.0.d0) cycle
         h2=hp(i)**2
         invh5=1.d0/(hp(i)*h2*h2)
         pres_i=por2(i)*rho(i)**2
         do in=1,nn(i)
            j=list(first(i)+in)
            if(j.eq.i .or. u(j).eq.0.d0) cycle
            dx=x(i)-x(j)
            dy=y(i)-y(j)
            dz=z(i)-z(j)
            if(periodic) call pbcwrap(dx,dy,dz)
            r2=dx*dx+dy*dy+dz*dz
            if(r2.ge.4.d0*h2) cycle
            itab=int(ctab*r2/h2)+1
c     grad_i W = dw*(r_i-r_j), as in the force loop
            dw=dwtab(itab)*invh5
            f=am(j)/rho(j)*dw
            pres_j=por2(j)*rho(j)**2
            d(1)=dx
            d(2)=dy
            d(3)=dz
            do k=1,3
               rgrad(k,i)=rgrad(k,i)+f*(rho(j)-rho(i))*d(k)
               rgrad(3+k,i)=rgrad(3+k,i)+f*(pres_j-pres_i)*d(k)
               rgrad(6+k,i)=rgrad(6+k,i)+f*(vx(j)-vx(i))*d(k)
               rgrad(9+k,i)=rgrad(9+k,i)+f*(vy(j)-vy(i))*d(k)
               rgrad(12+k,i)=rgrad(12+k,i)+f*(vz(j)-vz(i))*d(k)
            enddo
         enddo
      enddo

c     share, one component at a time through a contiguous buffer
      do k=1,15
         do i=1,ntot
            buf(i)=rgrad(k,i)
         enddo
         if(blkdist) then
            call sync_w(buf)
         else
            call sync_all(buf)
         endif
         do i=1,ntot
            rgrad(k,i)=buf(i)
         enddo
      enddo
      return
      end
************************************************************************
      subroutine riemann_limit(qi,qj,si,sj,ai,aj)
c     Midpoint states.  qi, qj are the values at particles i and j, and
c     si, sj the changes that the gradients at i and at j predict over the
c     separation from i to j, grad q.(r_j-r_i).  Where either slope
c     disagrees in sign with qj-qi the pair keeps the particle values
c     (first order, Inutsuka 2002, sec. 3.3).  Otherwise i is extrapolated
c     half way towards j and j half way back towards i, and both results
c     are kept between qi and qj.  If the two would cross, both take their
c     mean, so the reconstructed difference never has the opposite sign to
c     qj-qi.  For the velocity this keeps the dissipation of every pair
c     positive: p*-pstat has the sign of the approach speed, so the pair is
c     heated whether it approaches or recedes.
      implicit none
      real*8 qi,qj,si,sj,ai,aj,dq,lo,hi
      dq=qj-qi
      if(si*dq.le.0.d0 .or. sj*dq.le.0.d0) then
         ai=qi
         aj=qj
         return
      endif
      lo=min(qi,qj)
      hi=max(qi,qj)
      ai=min(max(qi+0.5d0*si,lo),hi)
      aj=min(max(qj-0.5d0*sj,lo),hi)
      if((aj-ai)*dq.lt.0.d0) then
         ai=0.5d0*(ai+aj)
         aj=ai
      endif
      return
      end
