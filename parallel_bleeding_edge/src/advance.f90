! 
! pretty much where everything happens
!
!*********************************************************
      subroutine advance
!     implements second-order accurate leap-frog sph step
      include 'starsmasher.h'
      real*8 uo(nmax),vxo(nmax),vyo(nmax),vzo(nmax)
      common/oldarrays/uo,vxo,vyo,vzo
      integer i
      real*8 dthnew

!     variables used for radiative cooling portion of the code:
!         uorig=specific internal energy u particle would have achieved if no cooling
!         ttherm=thermal timescale
!         erad=energy lost to radiation
!         amnotte=mass not in thermal equilibrium
      real*8 uorig,ttherm,erad,amnotte
      common/lostenergy/ erad
      real*8 displacexdot, displaceydot,displacezdot,amdisplace
      real*8 displacex, displacey,displacez
      integer ndisplace
      common/displace/displacex,displacey,displacez,ndisplace

      ! coming into this routine, positions, udot, and accelerations are a half-timestep
      ! behind the specific internal energies u and the velocities
      
      dth=0.5d0*dt

      displacexdot=0d0
      displaceydot=0d0
      displacezdot=0d0

      if(ndisplace.ne.0) then
         ! Let's figure out how much the center of mass of the stars (but not the black
         ! hole!) will shift in this timestep

         amdisplace=0d0
         do i=1,n-1 ! ends at n-1 because don't want to consider the black hole
            displacexdot=displacexdot+am(i)*vx(i)
            displaceydot=displaceydot+am(i)*vy(i)
            displacezdot=displacezdot+am(i)*vz(i)
            amdisplace=amdisplace+am(i)
         enddo

         if(u(n).ne.0) then
            write(69,*) 'Was assuming last particle was a black hole.'
            write(69,*) 'If that is not the case, then adjust src code.'
            stop
         endif

         displacexdot=displacexdot/amdisplace
         displaceydot=displaceydot/amdisplace
         displacezdot=displacezdot/amdisplace
         displacex=displacex+displacexdot*dt
         displacey=displacey+displaceydot*dt
         displacez=displacez+displacezdot*dt

         if(myrank.eq.0) write(69,'(4g17.9)')'displace{x,y,z}=',displacex,displacey,displacez

      endif

      do i=1,n
         x(i)=x(i)+dt*(vx(i)-displacexdot)    ! update x to half-timestep
         y(i)=y(i)+dt*(vy(i)-displaceydot)    ! update y to half-timestep
         z(i)=z(i)+dt*(vz(i)-displacezdot)    ! update z to half-timestep
         vxo(i)=vx(i)           ! store old value of vx
         vyo(i)=vy(i)           ! store old value of vy
         vzo(i)=vz(i)           ! store old value of vz
         vx(i)=vx(i)+dth*vxdot(i) ! update vx to half-timestep
         vy(i)=vy(i)+dth*vydot(i) ! update vy to half-timestep
         vz(i)=vz(i)+dth*vzdot(i) ! update vz to half-timestep
         if(u(i).ne.0.d0) then
            uo(i)=u(i)             ! store old value of u
            if(ncooling.eq.0)then
               u(i)=u(i)+dth*udot(i)  ! update u to half-timestep
            else
!               ttherm=(ueq(i)-u(i))/uraddoti
               ttherm=abs(tthermal(i))
               u(i)=u(i)*exp(-dth/ttherm)+ueq(i)*(1-exp(-dth/ttherm))+udot(i)*dth
            endif
         else if(udot(i).ne.0.d0) then
            write(6,*)'warning: black hole particle has udot=',udot(i),i,myrank
            if(myrank.eq.0) write(69,*)'warning: black hole particle has udot=',udot(i),i
            stop
         endif

      enddo

      t=t+dth

!     if binary relaxation, adjust center of mass positions and velocities:
      if (nrelax.ge.2 .and. .not.gonedynamic) call cmadj

      if(nrelax.gt.0 .and. t.lt.tresplintmuoff) call resplintmu

      ! Note: we call rho_and_h before gravforce because the gravity calculation can
      ! need the most up-to-date smoothing lengths for the integrator to achieve
      ! second order accuracy
      call rho_and_h            ! evaluate rho and h at half-timestep

!     get potential energy at half-timestep:
!      if( ngr.ne.0 .and. mod(nit,nitpot).eq.0 ) call gravpot
      if(ngr.ne.0) call gravforce

      if(myrank.eq.nprocs-1) call cpu_time(time1)
      call uvdots                ! evaluate udot and accelerations at half-timestep
      if(myrank.eq.nprocs-1) then
         call cpu_time(time2)
         write (6,'(a,f6.3,a,i4)')&
              'uvdots: ',time2-time1
      endif

!     output stuff to standard output (for log file) each iteration:
!	  if(myrank.eq.0)
      call enout(.true.)

!      if(ngr.ne.0) call gravforce

      call tstep
      dthnew=0.5d0*dt

!     advance velocities:
      do i=1,n
         vx(i)=vxo(i)+(dth+dthnew)*vxdot(i) ! update vx next half-timestep
         vy(i)=vyo(i)+(dth+dthnew)*vydot(i) ! update vy next half-timestep
         vz(i)=vzo(i)+(dth+dthnew)*vzdot(i) ! update vz next half-timestep
      enddo
!     advance specific internal energies (nintvar=2) or entropies (nintvar=1):
      amnotte=0.d0
      do i=1,n
         if(u(i).ne.0.d0) then
            if(ncooling.eq.0) then
               u(i)=uo(i)+(dth+dthnew)*udot(i) ! update u to full-timestep
            else
               uorig=uo(i)+(dth+dthnew)*udot(i) ! update u to full-timestep
!               ttherm=(ueq(i)-uo(i))/uraddoti
               ttherm=abs(tthermal(i))
               
               if(myrank.eq.0 .and. ttherm.lt.1000) then
                  !               print *,'particle i=',i,'has a thermal timescale=',ttherm
                  amnotte=amnotte+am(i)
               endif

!           note if the timestep dt=dth+dthnew<<ttherm, then exp(-dt/ttherm)=1-dt/therm+...
!           and the next line of code is approximately equivalent to
!           u(i) = uo(i) + (ueq(i)-uo(i))*dt/ttherm + dt*udot(i)
!                = uo(i) + dt*(uraddoti+udot(i))
!           which is exactly what we expect for large thermal timescales
!
!           if instead the timestep dt=dth+dthnew>>ttherm, then exp(-dt/ttherm)~0
!           and the next line of code is approximately equivalent to
!           u(i) = ueq(i) + dt*udot(i)
!           which seems reasonable since if if ttherm is small then uo(i) ~ ueq(i)
!           and expansion, compressions, or shocks should still allow u(i) to change through
!           the dt*udot(i) term.
               u(i)=uo(i)*exp(-(dth+dthnew)/ttherm)+ueq(i)*(1-exp(-(dth+dthnew)/ttherm))+(dth+dthnew)*udot(i)
               erad=erad+am(i)*(uorig-u(i))

!               if(erad.ne.erad)then
!                  write(70,*)i,ttherm,ueq(i),udot(i),uo(i),u(i)
!                  write(70,*)'first term:',uo(i)*exp(-(dth+dthnew)/ttherm)
!                  write(70,*)'secon term:',ueq(i)*(1-exp(-(dth+dthnew)/ttherm))
!                  write(70,*)'third term:',(dth+dthnew)*udot(i)
!                  stop
!               endif

            endif
         endif
      enddo

      if(myrank.eq.0 .and. ncooling.ne.0)&
           write(69,*),'mass not in thermal equilibrium=',amnotte

      t=t+dthnew

!     if binary relaxation, adjust center of mass positions and velocities:
      if (nrelax.ge.2 .and. .not.gonedynamic) call cmadj

      return
      end
!***********************************************************************
      subroutine cmadj                                                 
!*****************************************************************
!     release 1.0
!     adjust center-of-mass for binary relaxation
!     called by advance
!******************************************************************
      include 'starsmasher.h'
      real*8 xcm1,ycm1,zcm1,xcm2,ycm2,zcm2,am1,am2
      real*8 vxcm1,vycm1,vzcm1,vxcm2,vycm2,vzcm2
      common/centersofmass/xcm1,ycm1,zcm1,xcm2,ycm2,zcm2,am1,am2
      common/vcentersofmass/vxcm1,vycm1,vzcm1,vxcm2,vycm2,vzcm2
      real*8 uo(nmax),vxo(nmax),vyo(nmax),vzo(nmax)
      common/oldarrays/uo,vxo,vyo,vzo

      real*8 delx1,dely1,delz1,delx2,dely2,delz2
      real*8 delvx1,delvy1,delvz1,delvx2,delvy2,delvz2
      real*8 sep1,tscanoff,scanconst
      integer i
      integer corepts1,corepts2
      common/core/corepts1,corepts2

      tscanoff=min(tf,treloff)

! set tscanon to a negative number if you don't want to scan
! or set sepfinal to the same value as sep0
      if(tscanon.ge.0.d0 .and. t.ge.tscanon .and. t.lt.tscanoff) then
!         sep1=sep0+(sepfinal-sep0)*(t-tscanon)/(tscanoff-tscanon)
         scanconst=log(sep0/sepfinal)/(tscanoff-tscanon)
         sep1=sep0*exp(-scanconst*(t-tscanon))
         
         if(myrank.eq.0)write(69,*)'binary separation ',sep1,' at time=',t
      else if(t.ge.tscanoff)then
         sep1=sepfinal
      else
         sep1=sep0
      endif

      call getcoms

      delx1=-xcm1-am2*sep1/(am1+am2)                                   
      dely1=-ycm1
      delz1=-zcm1
      delx2=-xcm2+am1*sep1/(am1+am2)
      dely2=-ycm2
      delz2=-zcm2
      delvx1=-vxcm1
      delvy1=-vycm1
      delvz1=-vzcm1
      delvx2=-vxcm2
      delvy2=-vycm2
      delvz2=-vzcm2

      do i=1,n1
         x(i)=x(i)+delx1
         y(i)=y(i)+dely1
         z(i)=z(i)+delz1
         vx(i)=vx(i)+delvx1
         vy(i)=vy(i)+delvy1
         vz(i)=vz(i)+delvz1
      enddo
      do i=n1+1,ntot
         x(i)=x(i)+delx2
         y(i)=y(i)+dely2
         z(i)=z(i)+delz2                                                
         vx(i)=vx(i)+delvx2
         vy(i)=vy(i)+delvy2
         vz(i)=vz(i)+delvz2
      enddo

      if(n1.eq.ntot/2 .and. am1.eq.am2) then
         print *,'force a precisely symmetric equal mass binary:'
         do i=1,n1
            x(i)=0.5d0*(x(i)-x(i+n1))
            vx(i)=0.5d0*(vx(i)-vx(i+n1))
            vxdot(i)=0.5d0*(vxdot(i)-vxdot(i+n1))
            gx(i)=0.5d0*(gx(i)-gx(i+n1))
            y(i)=0.5d0*(y(i)-y(i+n1))
            vy(i)=0.5d0*(vy(i)-vy(i+n1))
            vydot(i)=0.5d0*(vydot(i)-vydot(i+n1))
            gy(i)=0.5d0*(gy(i)-gy(i+n1))
            z(i)=0.5d0*(z(i)+z(i+n1))
            vz(i)=0.5d0*(vz(i)+vz(i+n1))
            vzdot(i)=0.5d0*(vzdot(i)+vzdot(i+n1))
            gz(i)=0.5d0*(gz(i)+gz(i+n1))
            hp(i)=0.5d0*(hp(i)+hp(i+n1))
            rho(i)=0.5d0*(rho(i)+rho(i+n1))
            por2(i)=0.5d0*(por2(i)+por2(i+n1))
            u(i)=0.5d0*(u(i)+u(i+n1))
            udot(i)=0.5d0*(udot(i)+udot(i+n1))
            uo(i)=0.5d0*(uo(i)+uo(i+n1))
            vxo(i)=0.5d0*(vxo(i)-vxo(i+n1))
            vyo(i)=0.5d0*(vyo(i)-vyo(i+n1))
            vzo(i)=0.5d0*(vzo(i)-vzo(i+n1))

            x(i+n1)=-x(i)
            vx(i+n1)=-vx(i)
            vxdot(i+n1)=-vxdot(i)
            gx(i+n1)=-gx(i)
            y(i+n1)=-y(i)
            vy(i+n1)=-vy(i)
            vydot(i+n1)=-vydot(i)
            gy(i+n1)=-gy(i)
            z(i+n1)=z(i)
            vz(i+n1)=vz(i)
            vzdot(i+n1)=vzdot(i)
            gz(i+n1)=gz(i)
            hp(i+n1)=hp(i)
            rho(i+n1)=rho(i)
            por2(i+n1)=por2(i)
            u(i+n1)=u(i)
            udot(i+n1)=udot(i)
            uo(i+n1)=uo(i)
            vxo(i+n1)=-vxo(i)
            vyo(i+n1)=-vyo(i)
            vzo(i+n1)=-vzo(i)
         enddo
      endif

! the delx1*x(i)/xcm1 and delx2*x(2)/xcm2 stuff is motivated by not wanting to open up a gap at x=0 in contact configurations...
!      do i=1,ntot
!         if(x(i).le.0.d0) then
!            x(i)=x(i)+delx1*x(i)/xcm1
!            y(i)=y(i)+dely1
!            z(i)=z(i)+delz1
!         endif
!         if(x(i).ge.0.d0) then
!            x(i)=x(i)+delx2*x(i)/xcm2
!            y(i)=y(i)+dely2
!            z(i)=z(i)+delz2
!         endif                                                            
!      enddo                               

!  		 write(69,'(a,6g15.7)')'star:',xcm1,xcm2,sep1,am1,am2,myrank
      if(myrank.eq.0)then
         write(69,'(a,4g15.7)')'cmadj: star 1:',delx1,dely1,delz1,am1
         write(69,'(a,4g15.7)')'cmadj: star 2:',delx2,dely2,delz2,am2
      endif
      return                                                           
      end
!***********************************************************************
      subroutine rho_and_h
!     calculate hp(i) and rho(i) according to the method of bonet
      use kdtree2_module
      include 'starsmasher.h'
      include 'mpif.h'
      double precision drhodhi,dphidhi
      real*8 divv(nmax),hpguess
      common/commdivv/divv
      integer maxit,i,j
      double precision ezrtsafe,xtmp,xacc,dxmax
      parameter (maxit=50)
      double precision df,dx,ff,fh,fl,xh,xl
      integer k
      real(kdkind), dimension(:,:), allocatable :: my_array
      real*8 r2
      integer cnt
      type(kdtree2_result), allocatable :: results(:)
      type(kdtree2), pointer :: tree2
      real*8 r2last
      integer in, jn
      real*8 dist2
      integer mylength,irank,ierr

      if(myrank.eq.nprocs-1) call cpu_time(time1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      call linkedlists
!     is being replaced with the creation of a kd tree...
!
      allocate(my_array(3,n))
      do k=1,n
         my_array(1,k)=x(k)
         my_array(2,k)=y(k)
         my_array(3,k)=z(k)
!         myhp(k)=0
      enddo
      tree2 => kdtree2_create(my_array,sort=.false.,rearrange=.true.) ! this is how you create a tree.
!     the above have replaced the call to linkedlists
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      allocate(results(n)) 
      do i=n_lower,n_upper

         if(i.eq.n_lower) then
            first(i)=0
         else
            first(i)=first(i-1)+nn(i-1)
         endif

         if(u(i).eq.0.d0) then
!            if(myrank.eq.0)write(69,*)'hp(core)=',hp(i),i
            r2last=-1.d30
            goto 60909
         endif

         ! Change hp(i) to mean hptilde(i) while solving eq.(A1) of GLPZ 2010.
         hp(i)=hp(i) - hfloor

         xacc=1.d-4*hp(i)
         hpguess=max(hp(i)*(1.d0+divv(i)*dt/3.d0),xacc)
         dxmax=min(max(xacc,abs(hp(i)*divv(i)*dt)),0.5d0*hpguess)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!         hp(i) = ezrtsafe(hpguess,xacc,i,dxmax)
!     will be replaced with the following spaghetti code...

         xl=hpguess
         dxmax=abs(dxmax)

         hp(i)=xl

         r2=4.d0*hp(i)**2.d0
         r2last=r2

         call kdtree2_r_nearest_around_point(tp=tree2,idxin=i,correltime=-1,&
              r2=r2,nfound=cnt,nalloc=n,results=results)
         nn(i)=cnt
         list(first(i)+1:first(i)+nn(i))=results(1:nn(i))%idx
         
         call bonetfuncd(i,fl,df)

         dx=fl/df
         if(abs(dx).lt.xacc) then
            ezrtsafe=xl-dx
            goto 60908
         elseif(fl.eq.0) then
            ezrtsafe=xl
            goto 60908
         endif

         dx=sign(min(min(abs(dx),dxmax),abs(0.5d0*xl)),dx)

10       xh=xl-dx

         hp(i)=xh

         r2=4.d0*hp(i)**2.d0
         if(r2.le.r2last) then
            ! there's no need to search the entire tree, because the
            ! h we are trying is smaller than the last h we tried!
            cnt=0
            do in=1,nn(i)
               jn=list(first(i)+in)
               dist2=(x(i)-x(jn))**2+(y(i)-y(jn))**2+(z(i)-z(jn))**2
               if(dist2.lt.r2) then
                  cnt=cnt+1
                  list(first(i)+cnt)=jn
               endif
            enddo
            nn(i)=cnt
         else
            call kdtree2_r_nearest_around_point(tp=tree2,idxin=i,&
                 correltime=-1,&
                 r2=r2,nfound=cnt,nalloc=n,results=results)
            nn(i)=cnt
            list(first(i)+1:first(i)+nn(i))=results(1:nn(i))%idx
         endif
         r2last=r2
      
         call bonetfuncd(i,fh,df)
         dx=fh/df
         
         if(abs(dx).lt.xacc) then
            ezrtsafe=xh-dx
            goto 60908
         elseif(fh.eq.0) then
            ezrtsafe=xh
            goto 60908
         endif
         
         dx=sign(min(abs(dx),dxmax),dx)
         
         ezrtsafe=xh-dx
         
         if(fl*fh.gt.0.d0) then
            dx=sign(dx,fh)
            xl=xh
            goto 10 ! root is not yet bracketed, drift a little more
         endif
         
         !     at this point, the root is bracketed.
         !     make sure that xl really does give a negative function value
         if(fl.gt.0.d0)then
            xtmp=xh
            xh=xl
            xl=xtmp
         endif
         
12       continue
         
         do j=1,maxit

            if( (ezrtsafe-xh)*(ezrtsafe-xl).ge.0.d0) then
               !     it turns out that newton-raphson is not
               !     such a good idea for this step
               !     let's change our next evaluation point
               dx=0.5d0*(xh-xl)
               ezrtsafe=xl+dx
            endif
            
            if(abs(dx).lt.xacc) goto 60908
            hp(i)=ezrtsafe

            r2=4.d0*hp(i)**2.d0
            if(r2.le.r2last) then
               ! there's no need to search the entire tree, because the
               ! h we are trying is smaller than the last h we tried!
               cnt=0
               do in=1,nn(i)
                  jn=list(first(i)+in)
                  dist2=(x(i)-x(jn))**2&
                       +(y(i)-y(jn))**2+(z(i)-z(jn))**2
                  if(dist2.lt.r2) then
                     cnt=cnt+1
                     list(first(i)+cnt)=jn
                  endif
               enddo
               nn(i)=cnt
            else
               call kdtree2_r_nearest_around_point(tp=tree2,idxin=i,&
                    correltime=-1,r2=r2,nfound=cnt,nalloc=n,&
                    results=results)
               nn(i)=cnt
               list(first(i)+1:first(i)+nn(i))=results(1:nn(i))%idx
            endif
            r2last=r2

            call bonetfuncd(i,ff,df)
            
            if(ff.lt.0.d0) then
               xl=ezrtsafe
            else
               xh=ezrtsafe
            endif
            
            !     let's assume for the moment that the next step will be a 
            !     newton-raphson one
            dx=ff/df
            ezrtsafe=ezrtsafe-dx
         enddo

         if(myrank.eq.0)write(69,*)&
              'warning: ezrtsafe exceeded maximum iterations...recovering',&
              i,ezrtsafe,ff,df,dx
         
         dx=0.5d0*(xh-xl)
         ezrtsafe=xl+dx
         goto 12
         
         
60908    hp(i)=ezrtsafe

         ! Change hp(i) back to the true smoothing length now that done solving eq. (A1) of GLPZ 2010.  
         hp(i) = hp(i) + hfloor
         
60909    continue  ! code is jumped to here for point particles (with u=0)

!     here ends the replacement for the line of code
!         hp(i) = ezrtsafe(hpguess,xacc,i,dxmax)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         r2=4.d0*hp(i)**2.d0
         if(r2.le.r2last) then
            ! there's no need to search the entire tree, because the
            ! h we are trying is smaller than the last h we tried!
            cnt=0
            do in=1,nn(i)
               jn=list(first(i)+in)
               dist2=(x(i)-x(jn))**2+(y(i)-y(jn))**2+(z(i)-z(jn))**2
               if(dist2.lt.r2) then
                  cnt=cnt+1
                  list(first(i)+cnt)=jn
               endif
            enddo
            nn(i)=cnt
         else
            call kdtree2_r_nearest_around_point(tp=tree2,idxin=i,&
                 correltime=-1,&
                 r2=r2,nfound=cnt,nalloc=n,results=results)
            nn(i)=cnt
            list(first(i)+1:first(i)+nn(i))=results(1:nn(i))%idx
         endif
         
         call bonetsumwithphi(i,rho(i), bonet_omega(i),&
              bonet_0mega(i),&
              bonet_psi(i), bonet_wn(i))

!         write(16,*) i,rho(i), bonet_omega(i),&
!             bonet_0mega(i),&
!             bonet_psi(i), bonet_wn(i)
!         write(17,*) i,x(i),y(i),z(i),hp(i)

        
!     compute f(i) array for springel and hernquist method.
!     the variable b is not exactly monaghan's b... 
!     instead it is his b divided by aa(i)
!     write(69,*) i, rho(i), bonet_wn(i)/hp(i)**3
        
!     f(i)=(1.d0+bb(i)*hp(i)**2*rho(i)**(-2.d0/3.d0)
!     $           *drhodhi/3.d0 )**(-1)
!     zeta(i)=-bb(i)*hp(i)**2*rho(i)**(-2.d0/3.d0)
!     $           *dphidhi/3.d0
         zeta(i) = 0;
         if (bonet_0mega(i) == 0) then
            bonet_0mega(i) = 1.0d-4
            bonet_omega(i) = 0
            bonet_psi(i)   = 0
         end if
!     bonet_psi(i) = 0;
!     bonet_omega(i) = 0;
!     write(69,*)i,hp(i),dxmax,hpguess,rho(i),nn(i),f(i)
!     hp(i)=hp(i)

!         if(i.eq.7855) write(78,'(999g11.4)') t,hp(i),rho(i),nn(i)!,results(1:nn(i))%idx
!         if(i.eq.13153) write(131,'(999g11.4)') t,hp(i),rho(i),nn(i)!,results(1:nn(i))%idx
!         if(i.eq.14562) write(145,'(999g11.4)') t,hp(i),rho(i),nn(i)!,results(1:nn(i))%idx
!         if(i.eq.24452) write(244,'(999g11.4)') t,hp(i),rho(i),nn(i)!,results(1:nn(i))%idx
!         if(i.eq.26249) write(262,'(999g11.4)') t,hp(i),rho(i),nn(i)!,results(1:nn(i))%idx
!         if(i.eq.28467) write(284,'(999g11.4)') t,hp(i),rho(i),nn(i)!,results(1:nn(i))%idx
!         if(i.eq.29660) write(296,'(999g11.4)') t,hp(i),rho(i),nn(i)!,results(1:nn(i))%idx
!         if(i.eq.39877) write(398,'(999g11.4)') t,hp(i),rho(i),nn(i)!,results(1:nn(i))%idx

!      myhp(i)=hp(i)

      enddo
      call kdtree2_destroy(tree2)
      deallocate(results)
      deallocate(my_array)

! make mpi call to insure all processors doing gravity have the same
! hp values (only nodes doing gravity need the hp values):
      mylength=n_upper-n_lower+1
      do irank=0,ngravprocs-1
         if(myrank.ne.irank) then
            call mpi_gatherv(hp(n_lower), mylength, mpi_double_precision,&
                 hp, recvcounts, displs, mpi_double_precision, irank,&
                 mpi_comm_world, ierr)
         else
            call mpi_gatherv(mpi_in_place, mylength, mpi_double_precision,&
                 hp, recvcounts, displs, mpi_double_precision, irank,&
                 mpi_comm_world, ierr)
         endif
      enddo

      call pressure
      if(myrank.eq.nprocs-1) then
         call cpu_time(time2)
         write (6,'(a,f6.3,a,i4)')&
              'rho\&h: ',time2-time1
      endif

      return
      end
!***********************************************************************
      subroutine bonetsum(i,bonet1_wn, bonet1_0mega)
!     compute density of particle i:
      include 'starsmasher.h'                                         
      real*8 r2,h3,h2,hpi, h4
      real*8 bonet1_wn, bonet1_0mega
      integer itab,j,in,i
      real*8 ctaboverh2
      real*8 dx, dy, dz
      real*8 fac
      integer offset, nni
!      if(u(i).eq.0.d0) then
!         write(69,*)'bonetsum should not be called for core points'
!         stop
!      endif

      bonet1_wn    = 0 
      bonet1_0mega = 0

      hpi=hp(i)
      h2=hpi**2
      ctaboverh2=ctab/h2
      h3=hpi*h2
      h4=h2*h2
!     accumulate contributions:
      offset = first(i)
      nni    = nn(i)
      do in=1,nni
        j=list(offset + in)
!       eg: to vectorize this loop there need to be a way to cast logical into real
!           will be working on it
        if(u(j).ne.0.d0) then
          fac = 1.d0
        else
           fac = 0.d0
        end if
        dx = x(i) - x(j)
        dy = y(i) - y(j)
        dz = z(i) - z(j)
        r2 = dx*dx + dy*dy + dz*dz;
        itab=int(ctaboverh2*r2) + 1
        bonet1_wn     = bonet1_wn    + fac*gtab(itab)
        bonet1_0mega  = bonet1_0mega - fac*dgdhtab(itab)
      enddo
      bonet1_0mega = bonet1_0mega/hpi

!       write(69,*) 'i= ', i, 'h= ', hp(i), 'nn= ', nn(i), 'wn= ', bonet1_wn
      
      return
      end
!***********************************************************************
      subroutine bonetsumwithphi(i,bonet1_rho,&
          bonet1_omega, bonet1_0mega, bonet1_psi, bonet1_wn)
!     compute density of particle i:
      include 'starsmasher.h'                                         
      real*8 r2,h3,h2,hpi, h4
      real*8 bonet1_rho, bonet1_omega, bonet1_0mega
      real*8 bonet1_psi, bonet1_wn
      integer itab,j,in,i
      real*8 ctaboverh2
      real*8 bonetdphidh
      real*8 hpitilde,h2tilde,ctaboverh2tilde
      integer itabtilde

      hpi = hp(i)
      h2=hpi**2
      ctaboverh2=ctab/h2
      if(u(i).ne.0.d0) then

         bonet1_rho   = 0.d0
         bonet1_omega = 0.d0
         bonet1_0mega = 0.d0
         bonet1_psi   = 0.d0
         bonet1_wn    = 0.d0
         
         hpitilde=hpi - hfloor
         h2tilde=hpitilde**2
         ctaboverh2tilde=ctab/h2tilde
         h3=hpi*h2
         h4=h2*h2
         !     accuumulate contributions:
         do in=1,nn(i)
            j=list(first(i)+in)
            r2=(x(i)-x(j))**2.d0+(y(i)-y(j))**2.d0+(z(i)-z(j))**2.d0
            itab=int(ctaboverh2*r2)+1

            ! the u(j).eq.0 is necessary in the next line because we always do gravity
            ! with point particles
            if(nselfgravity.eq.1 .or. u(j).eq.0.d0)&
                 bonet1_psi   = bonet1_psi   + am(j)*dphidhtab(itab)
            if(u(j).ne.0.d0) then
               if(r2.lt.4d0*h2tilde) then
                  itabtilde=int(ctaboverh2tilde*r2)+1
                  bonet1_wn    = bonet1_wn    + gtab(itabtilde)
                  bonet1_0mega = bonet1_0mega - dgdhtab(itabtilde)
               endif
               bonet1_rho   = bonet1_rho   + am(j)*wtab(itab)
               bonet1_omega = bonet1_omega + am(j)*dwdhtab(itab)
            endif
         enddo
         bonet1_rho   = bonet1_rho/h3
         bonet1_omega = bonet1_omega/h4
         bonet1_psi   = bonet1_psi/h2
         bonet1_0mega = bonet1_0mega/hpitilde
      else
         bonet1_rho   = 0.d0
         do in=1,nn(i)
            j=list(first(i)+in)
            r2=(x(i)-x(j))**2.d0+(y(i)-y(j))**2.d0+(z(i)-z(j))**2.d0
            itab=int(ctaboverh2*r2)+1
            if(u(j).ne.0.d0) then
               bonet1_rho   = bonet1_rho   + am(j)*wtab(itab)
            endif
         enddo
         bonet1_rho   = bonet1_rho/hpi**3
         bonet1_omega = 0.d0
         bonet1_psi   = 0.d0
         bonet1_0mega = 1.d30
      endif

      return
      end

!***********************************************************************
      subroutine bonetfuncd(i,bonetfunc,dbonetfuncdh)
      include 'starsmasher.h'
      real*8 bonetfunc,hpi,bonetrhos,bonetdrhodh,dbonetfuncdh
      integer i

      call bonetsum(i, bonetrhos, bonetdrhodh)
!      if(u(i).ne.0.d0) then
         bonetfunc    = bonetrhos - nnopt
!      else
!         write(69,*) 'bonetfuncd should not be called for core'
!         stop
!      endif
      dbonetfuncdh = bonetdrhodh

      return
      end
!***********************************************************************
      subroutine vdotsm
!     compute smoothed acceleration of particle i:
      include 'starsmasher.h'                                         
      real*8 r2,h3,h2,hpi
      integer itab,j,in,i
      real*8 invh2

      do i=n_lower,n_upper
         vxdotsm(i)=0.d0
         vydotsm(i)=0.d0
         vzdotsm(i)=0.d0
         hpi=hp(i)
         h2=hpi**2
         h3=hpi*h2
         invh2 = 1.d0/h2
         
!     accumulate contributions:
         do in=1,nn(i)
            j=list(first(i)+in)
            if(u(j).ne.0.d0)then
               r2=(x(i)-x(j))**2.d0+(y(i)-y(j))**2.d0+(z(i)-z(j))**2.d0
               itab=int(ctab*r2*invh2)+1
               vxdotsm(i)=vxdotsm(i)+am(j)*&
                    wtab(itab)*vxdot(j)
               vydotsm(i)=vydotsm(i)+am(j)*&
                    wtab(itab)*vydot(j)
               vzdotsm(i)=vzdotsm(i)+am(j)*&
                    wtab(itab)*vzdot(j)
            endif
         enddo
         vxdotsm(i)=vxdotsm(i)/(h3*rho(i))
         vydotsm(i)=vydotsm(i)/(h3*rho(i))
         vzdotsm(i)=vzdotsm(i)/(h3*rho(i))
      enddo

      return
      end
!***********************************************************************
      subroutine getcoms
!     get the centers of mass
      include 'starsmasher.h'
      real*8 xcm1,ycm1,zcm1,xcm2,ycm2,zcm2,am1,am2
      real*8 vxcm1,vycm1,vzcm1,vxcm2,vycm2,vzcm2
      common/centersofmass/xcm1,ycm1,zcm1,xcm2,ycm2,zcm2,am1,am2
      common/vcentersofmass/vxcm1,vycm1,vzcm1,vxcm2,vycm2,vzcm2
      integer i
      integer corepts1,corepts2
      common/core/corepts1,corepts2

      xcm1=0.d0
      ycm1=0.d0
      zcm1=0.d0
      am1=0.d0                                                        
      xcm2=0.d0
      ycm2=0.d0
      zcm2=0.d0
      vxcm1=0.d0
      vycm1=0.d0
      vzcm1=0.d0
      vxcm2=0.d0
      vycm2=0.d0
      vzcm2=0.d0
      am2=0.d0


      if(myrank.eq.0) write(69,*) 'n,ntot,corepts1',n,ntot,corepts1

      do i=1,n1
         am1=am1+am(i)
         xcm1=xcm1+am(i)*x(i)                                           
         ycm1=ycm1+am(i)*y(i)                                           
         zcm1=zcm1+am(i)*z(i)                                           
         vxcm1=vxcm1+am(i)*vx(i)
         vycm1=vycm1+am(i)*vy(i)
         vzcm1=vzcm1+am(i)*vz(i)        
      enddo
      do i=n1+1,ntot
         am2=am2+am(i)
         xcm2=xcm2+am(i)*x(i)
         ycm2=ycm2+am(i)*y(i)
         zcm2=zcm2+am(i)*z(i)
         vxcm2=vxcm2+am(i)*vx(i)
         vycm2=vycm2+am(i)*vy(i)
         vzcm2=vzcm2+am(i)*vz(i)
      enddo

      if(am1.ne.0.d0) then
         xcm1=xcm1/am1
         ycm1=ycm1/am1
         zcm1=zcm1/am1  
         vxcm1=vxcm1/am1
         vycm1=vycm1/am1
         vzcm1=vzcm1/am1  
      endif
      if(am2.ne.0.d0) then
         xcm2=xcm2/am2
         ycm2=ycm2/am2
         zcm2=zcm2/am2 
         vxcm2=vxcm2/am2
         vycm2=vycm2/am2
         vzcm2=vzcm2/am2 
      endif

      if(myrank.eq.0)write(69,'(a,7g11.3)')'getcoms:',xcm1,xcm2

      return
      end

