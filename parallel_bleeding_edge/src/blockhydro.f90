!***********************************************************************
      subroutine sync_all(arr)
!     Every rank computes arr for its own particles (n_lower..n_upper).
!     Gather the full array onto every rank.
      include 'starsmasher.h'
      include 'mpif.h'
      real*8 arr(nmax)
      integer ierr
      if(nprocs.gt.1) call mpi_allgatherv(mpi_in_place,0,mpi_datatype_null,&
           arr,recvcounts,displs,mpi_double_precision,mpi_comm_world,ierr)
      return
      end
!***********************************************************************
      subroutine uvdots_active
!     Block timesteps (nblock=1): accelerations and du/dt for the active
!     particles only (actblk).  Equivalent to uvdots, but each pair (a,j) with
!     a active is handled from a's side in full: the a-side ("gather") term
!     with a's quantities when r<2h_a, and the j-side ("scatter") term with
!     j's quantities -- current if j is active, predicted/stored otherwise --
!     when r<2h_j.  Neighbours are found by a direct loop over all particles,
!     the same order of cost as the N_active x N gravity sum.
!     Pair terms: eqs. (A11)-(A12), (A14)-(A15) of Gaburov et al. (2010) and
!     the Balsara-switched artificial viscosity, exactly as in uvdots.
!     Requires nintvar=2, nrelax=0, ncooling=0, hfloor=0 (checked in advance_block).
      include 'starsmasher.h'
      include 'mpif.h'
      real*8 uijmax(nmax)
      common/uijmax/ uijmax
      double precision divv(nmax)
      common/commdivv/divv
      real*8 curlabs(nmax)
      common/blockcurl/curlabs
      real*8 myax(nmax),myay(nmax),myaz(nmax),myud(nmax),myuij(nmax)
      real*8 alpha1,beta1,dx,dy,dz,divin,r2,ca,cj,fa,fj,udb
      real*8 csa,cga,csj,cgj,pij,dwij,diwk,hd,rinvh2
      real*8 curlvxi,curlvyi,curlvzi,divvi,grpottot(nmax)
      integer a,j,itab,ierr

      alpha1=-2.d0*alpha
      beta1=2.d0*beta
      do a=1,ntot
         myax(a)=0.d0
         myay(a)=0.d0
         myaz(a)=0.d0
         myud(a)=0.d0
         myuij(a)=0.d0
      enddo

!     div v and |curl v| of active gas particles (from their own lists)
      do a=n_lower,n_upper
         if(.not.(actblk(a) .or. refblk(a)) .or. u(a).eq.0.d0) cycle
         call getderivs(a,curlvxi,curlvyi,curlvzi,divvi)
         divv(a)=divvi
         curlabs(a)=sqrt(curlvxi**2+curlvyi**2+curlvzi**2)
      enddo
      call sync_all(divv)
      call sync_all(curlabs)

      do a=n_lower,n_upper
         if(.not.actblk(a)) cycle
         if(u(a).ne.0.d0) then
            ca=sqrt(gam*por2(a)*rho(a))
            fa=abs(divv(a))/(abs(divv(a))+curlabs(a)+0.00001d0*ca/hp(a))
         endif
         do j=1,ntot
            if(j.eq.a) cycle
            dx=x(a)-x(j)
            dy=y(a)-y(j)
            dz=z(a)-z(j)
            r2=dx*dx+dy*dy+dz*dz
            if(r2.ge.4.d0*max(hp(a),hp(j))**2) cycle
            divin=dx*(vx(a)-vx(j))+dy*(vy(a)-vy(j))+dz*(vz(a)-vz(j))
            csa=0.d0
            cga=0.d0
            csj=0.d0
            cgj=0.d0
!     a-side term (j inside a's kernel)
            if(r2.lt.4.d0*hp(a)**2) then
               rinvh2=1.d0/hp(a)**2
               itab=int(ctab*r2*rinvh2)+1
               diwk=dgtab(itab)*rinvh2
               dwij=dwtab(itab)*rinvh2**2.5d0
               if(u(a).ne.0.d0 .and. u(j).ne.0.d0) then
                  pij=0.d0
                  if(divin.lt.0.d0) then
                     if(nav.eq.3) then
                        udb=divin/(ca*sqrt(r2))*fa
                     else
                        udb=hp(a)*divin/(ca*r2)*fa
                     endif
                     pij=por2(a)*(alpha1+beta1*udb)*udb
                  endif
                  csa=por2(a)*(-bonet_omega(a)/(bonet_0mega(a)*am(j))*diwk+dwij)
                  cga=0.5d0*pij*dwij
                  myud(a)=myud(a)+am(j)*(csa+cga)*divin
                  myuij(a)=max(myuij(a),(por2(a)+0.5d0*pij)*rho(a))
                  if(ngr.ne.0) cga=cga-0.5d0*diwk*bonet_psi(a)/(bonet_0mega(a)*am(j))
               else if(u(a).eq.0.d0 .and. u(j).ne.0.d0 .and. dynhco.ge.1&
                    .and. ngr.ne.0) then
                  if(hcolim) then
                     hd=hdynco(a)
                     if(r2.lt.4.d0*hd**2) then
                        diwk=dgtab(int(ctab*r2/hd**2)+1)/hd**2
                     else
                        diwk=0.d0
                     endif
                  endif
                  cga=-0.5d0*diwk*bonet_psi(a)/(bonet_0mega(a)*am(j))
               endif
            endif
!     j-side term (a inside j's kernel), with j's quantities
            if(r2.lt.4.d0*hp(j)**2) then
               rinvh2=1.d0/hp(j)**2
               itab=int(ctab*r2*rinvh2)+1
               diwk=dgtab(itab)*rinvh2
               dwij=dwtab(itab)*rinvh2**2.5d0
               if(u(a).ne.0.d0 .and. u(j).ne.0.d0) then
                  pij=0.d0
                  if(divin.lt.0.d0) then
                     cj=sqrt(gam*por2(j)*rho(j))
                     fj=abs(divv(j))/(abs(divv(j))+curlabs(j)+0.00001d0*cj/hp(j))
                     if(nav.eq.3) then
                        udb=divin/(cj*sqrt(r2))*fj
                     else
                        udb=hp(j)*divin/(cj*r2)*fj
                     endif
                     pij=por2(j)*(alpha1+beta1*udb)*udb
                  endif
                  csj=por2(j)*(-bonet_omega(j)/(bonet_0mega(j)*am(a))*diwk+dwij)
                  cgj=0.5d0*pij*dwij
                  if(ngr.ne.0) cgj=cgj-0.5d0*diwk*bonet_psi(j)/(bonet_0mega(j)*am(a))
               else if(u(j).eq.0.d0 .and. u(a).ne.0.d0 .and. dynhco.ge.1&
                    .and. ngr.ne.0) then
                  if(hcolim) then
                     hd=hdynco(j)
                     if(r2.lt.4.d0*hd**2) then
                        diwk=dgtab(int(ctab*r2/hd**2)+1)/hd**2
                     else
                        diwk=0.d0
                     endif
                  endif
                  cgj=-0.5d0*diwk*bonet_psi(j)/(bonet_0mega(j)*am(a))
               endif
            endif
!     acceleration of a:  -m_j [(a-side) + (j-side)] (x_a - x_j)
            myax(a)=myax(a)-am(j)*(csa+cga+csj+cgj)*dx
            myay(a)=myay(a)-am(j)*(csa+cga+csj+cgj)*dy
            myaz(a)=myaz(a)-am(j)*(csa+cga+csj+cgj)*dz
         enddo
      enddo

!     gravity on the active particles
      if(ngr.ne.0 .and. myrank.lt.ngravprocs) then
         call get_gravity_active
         if(ngravprocs.gt.1) then
            call mpi_reduce(grpot,grpottot,ntot,mpi_double_precision,&
                 mpi_sum,0,mpi_comm_world,ierr)
            if(myrank.eq.0) then
               do a=1,ntot
                  grpot(a)=grpottot(a)
               enddo
            endif
         endif
         do a=1,ntot
            myax(a)=myax(a)+gx(a)
            myay(a)=myay(a)+gy(a)
            myaz(a)=myaz(a)+gz(a)
         enddo
      endif

      call mpi_allreduce(myax,vxdot,ntot,mpi_double_precision,mpi_sum,&
           mpi_comm_world,ierr)
      call mpi_allreduce(myay,vydot,ntot,mpi_double_precision,mpi_sum,&
           mpi_comm_world,ierr)
      call mpi_allreduce(myaz,vzdot,ntot,mpi_double_precision,mpi_sum,&
           mpi_comm_world,ierr)
      call mpi_allreduce(myud,udot,ntot,mpi_double_precision,mpi_sum,&
           mpi_comm_world,ierr)
!     uijmax of inactive particles is kept (tstep uses it for their dtpart)
      call mpi_allreduce(mpi_in_place,myuij,ntot,mpi_double_precision,mpi_max,&
           mpi_comm_world,ierr)
      do a=1,ntot
         if(actblk(a)) uijmax(a)=myuij(a)
      enddo
      return
      end
