      subroutine compbest3(makeunboundstarcomp4)
c     pplot61: ri{1,2,3} can be zero  works well!
c     pplot62: ri(1,2,3} can not be less than smoothing length also works well!

c     analyze the results of a collision run and output important
c     numbers for a summary table
      include 'starsmasher.h'
      integer ncitmax,i,ic,j,icenter,icompfocus,nchng,ncit,in
      integer nicomp1,nicomp2,nicomp3,nicomp4
      real*8 spinz,spiny,spinx
      real*8 zrel,yrel,xrel
      real*8 rhocgs,spintotal,vzrel,vyrel,vxrel
      real*8 betaofi,prad,pgas,ucgs
      real*8 theta,dotproduct,rad,gamma,rhomax
      real*8 angularvelocity,rcyl,rhomaxfocus,temperature
      real*8 bi3,v2i3,bi2,v2i2,bi1,v2i1
      real*8 enth(nmax),am4
      real*8 x1,y1,z1,vx1,vy1,vz1,am1
      real*8 x2,y2,z2,vx2,vy2,vz2,am2
      real*8 x3,y3,z3,vx3,vy3,vz3,am3
      integer icomp(nmax)
      common/compbettercom3/am1,x1,y1,z1,vx1,vy1,vz1,
     $     am2,x2,y2,z2,vx2,vy2,vz2,
     $     am3,x3,y3,z3,vx3,vy3,vz3,am4,
     $     icomp
      parameter (ncitmax=100)
      logical db,pp
      character*18 fname
      integer ncomp1,ncomp2,ncomp3,ncomp4,corepts,corepts1
      real*8 amtiny
      integer irhomax4,irhomax1,irhomax2,irhomax3
      save irhomax4,irhomax1,irhomax2,irhomax3
      real*8 ri1,ri2,ri3

      logical firstt
      save firstt
      data firstt /.true./
      integer numstars,numstarslast
      common/nsl/ numstarslast
      real*8 am1last,am2last,am3last
      save am1last,am2last,am3last
      real*8 rhomax1,rhomax2,rhomax3,rhomax4
      integer swapcase,compa,compb,irhomaxtmp
      real*8 massa, massb
      logical checkifswapped
      common/cis/ checkifswapped
      real*8 amenclosed(nmax)
      logical makeunboundstarcomp4
      real*8 e1,e2,e3,mu
      integer unboundcomp
      real*8 x12,y12,z12,x13,y13,z13,x23,y23,z23
      real*8 vx12,vy12,vz12,vx13,vy13,vz13,vx23,vy23,vz23
      integer ccparents(4)
      integer numparents
      include 'mpif.h'
      integer mylength, ierr, irank

      corepts=0
      corepts1=corepts/2
      if(myrank.eq.0)write(69,*)'corepts=',corepts,' corepts1=',corepts1
      
      db=.true.                 ! if true, produce detailed output on convergence
      pp=.false.                ! if true, write pp files with component information    

c     Initializing ri1,ri2,ri3
      ri1=1.d30
      ri2=1.d30
      ri3=1.d30

 123  continue

      if(firstt) then
         ccparents(1)=cc(1)
         numparents=1
         do i=1,n
            do j=1,numparents
               if(cc(i).eq.ccparents(j)
     $              .or.(u(i).eq.0 .and. i.lt.n) ! this will put a point mass particle with star 1
     $              ) goto 1239
            enddo
            numparents=numparents+1
            ccparents(numparents)=cc(i)
c 1239       icomp(i)=mod(-j+1,3)+1! for component 1, then 3, then 2
c 1239       icomp(i)=mod(j+1,3)+1! for component 3, then 1, then 2
 1239       icomp(i)=j
c            if(am(i).gt.10.d0) then 
c               icomp(i)=2
c            else
c               icomp(i)=1
c            endif
         enddo

         if(numparents.lt.2) then
            numparents=numparents+1
            do i=n/2+1,n
               if(u(i).ne.0.d0) icomp(i)=numparents
            enddo
         endif
      endif

c     make sure all processors know the most recent rho values for all particles
      mylength=n_upper-n_lower+1
      do irank=0,nprocs-1
         if(myrank.ne.irank)then
            call mpi_gatherv(rho(n_lower), mylength, mpi_double_precision,
     $           rho, recvcounts, displs, mpi_double_precision, irank,
     $           mpi_comm_world, ierr) 
         else
            call mpi_gatherv(mpi_in_place, mylength, mpi_double_precision,
     $           rho, recvcounts, displs, mpi_double_precision, irank,
     $           mpi_comm_world, ierr) 
         endif
      enddo
         
      amtiny=1.d30
      do i=1,n
         amtiny=min(amtiny,am(i))
         if(u(i).eq.0.d0) then
            enth(i)=0
            corepts=corepts+1
            rho(i)=1.d30
         else
            if(nintvar.eq.1) then
c     p=(gam-1)*rho*u=a*rho^gam, so u=a*rho^(gam-1)/(gam-1)
               enth(i)=u(i)*rho(i)**(gam-1.d0)/(gam-1.d0)
            else
               enth(i)=u(i)     !use this if u(i) is actually specific energy
            endif
         endif
      enddo
      
      if(myrank.eq.0) then
         if(nintvar.eq.1) then
            write(69,*) '***assumed u(i) is truly p/rho^gamma***'
         else
            write(69,*) '**assumed u(i) is specific internal energy u**'
         endif
         write(69,*) 'smallest particle mass in system is',amtiny
      
c      write(69,*)'old irhomaxes=',irhomax1,irhomax2,irhomax3,irhomax4
c      write(69,*)icomp(irhomax1),icomp(irhomax2),icomp(irhomax3)
      endif

      nicomp1=0
      nicomp2=0
      nicomp3=0
      nicomp4=0

      rhomax1=0.d0
      rhomax2=0.d0
      rhomax3=0.d0
      rhomax4=0.d0
      do i=1,n
         if(icomp(i).eq.1) then
            if(rho(i).gt.rhomax1) then
               rhomax1=rho(i)
               irhomax1=i
            endif
            nicomp1=nicomp1+1
         else if(icomp(i).eq.2) then
            if(rho(i).gt.rhomax2) then
               rhomax2=rho(i)
               irhomax2=i
            endif
            nicomp2=nicomp2+1
         else if(icomp(i).eq.3) then
            if(rho(i).gt.rhomax3) then
               rhomax3=rho(i)
               irhomax3=i
            endif
            nicomp3=nicomp3+1
         else if(icomp(i).eq.4) then
            if(rho(i).gt.rhomax4) then
               rhomax4=rho(i)
               irhomax4=i
            endif
            nicomp4=nicomp4+1
         endif
      enddo

      if(myrank.eq.0) then
         write(69,*)'irhomaxes=',irhomax1,irhomax2,irhomax3,irhomax4

c      write(69,*)'rhomaxes=',rho(irhomax1),rho(irhomax2),
c     $     rho(irhomax3),rho(irhomax4)
c      write(69,*)'rhomaxes=',rhomax1,rhomax2,rhomax3,rhomax4
         write(69,*)'number of component 1 particles in initial guess=',
     $        nicomp1
         write(69,*)'number of component 2 particles in initial guess=',
     $        nicomp2
         write(69,*)'number of component 3 particles in initial guess=',
     $        nicomp3
         write(69,*)'number of component 4 particles in initial guess=',
     $        nicomp4
      endif
      
      call calccom

      if (db .and. myrank.eq.0) then
         write (69,'(6g13.6)') 'nit','nchng','m1','m2','m3'
         write (69,'(5g13.6)') 0,nicomp1+nicomp2+nicomp3,
     $        am1,am2,am3
      endif

c     iterate over component determination
      ncit=0
      nchng=1
      swapcase=0
      do while (nchng.ne.0 .and. ncit.lt.ncitmax)
         nchng=0
c         if(t.gt.122) write(69,*)ncit,t

         x1=x(irhomax1)
         y1=y(irhomax1)
         z1=z(irhomax1)
c         vx1=vx(irhomax1)
c         vy1=vy(irhomax1)
c         vz1=vz(irhomax1)
         x2=x(irhomax2)
         y2=y(irhomax2)
         z2=z(irhomax2)
c         vx2=vx(irhomax2)
c         vy2=vy(irhomax2)
c         vz2=vz(irhomax2)
         x3=x(irhomax3)
         y3=y(irhomax3)
         z3=z(irhomax3)
c         vx3=vx(irhomax3)
c         vy3=vy(irhomax3)
c         vz3=vz(irhomax3)

         do i=1,n
c     assign particles to components

c     calculate "total enthalpy" bi1 to component 1
            if(am1.gt.am(i)+0.9d0*amtiny) then
               ri1=sqrt((x(i)-x1)**2+(y(i)-y1)**2+(z(i)-z1)**2)+1d-30
               v2i1=(vx(i)-vx1)**2+(vy(i)-vy1)**2+(vz(i)-vz1)**2
               bi1=0.5d0*v2i1+enth(i)-am1/ri1
            else
               bi1=1.d30
            endif

c     calculate "total enthalpy" bi2 to component 2
            if(am2.gt.am(i)+0.9d0*amtiny) then
               ri2=sqrt((x(i)-x2)**2+(y(i)-y2)**2+(z(i)-z2)**2)+1d-30
               v2i2=(vx(i)-vx2)**2+(vy(i)-vy2)**2+(vz(i)-vz2)**2
               bi2=0.5d0*v2i2+enth(i)-am2/ri2
            else
               bi2=1.d30
            endif

c     calculate "total enthalpy" bi3 to component 3
            if(am3.gt.am(i)+0.9d0*amtiny) then
               ri3=sqrt((x(i)-x3)**2+(y(i)-y3)**2+(z(i)-z3)**2)+1d-30
               v2i3=(vx(i)-vx3)**2+(vy(i)-vy3)**2+(vz(i)-vz3)**2
               bi3=0.5d0*v2i3+enth(i)-am3/ri3
            else
               bi3=1.d30
            endif

c            if(i.eq.2514) then
c               write(69,*)'bs',i,bi1,bi2,bi3
c               write(69,*)'bs',i,ri1,ri2,ri3
c            endif

c     assign particle to one of 3 components
            if ( bi1.lt.0.d0 .and.
     $           (bi2.ge.0.d0 .or. ri1.lt.ri2) .and.
     $           (bi3.ge.0.d0 .or. ri1.lt.ri3)
     $           ) then
c     particle belongs to component 1
               ic=1
            elseif(bi2.lt.0.d0 .and.
     $              (bi1.ge.0.d0 .or. ri2.lt.ri1) .and.
     $              (bi3.ge.0.d0 .or. ri2.lt.ri3)
     $              )then
c     particle belongs to component 2
               ic=2
            elseif(bi3.lt.0.d0 .and.
     $              (bi1.ge.0.d0 .or. ri3.lt.ri1) .and.
     $              (bi2.ge.0.d0 .or. ri3.lt.ri2)
     $              )then
c     particle belongs to component 3
               ic=3
            else
c     particle is unbound (component 4)
               ic=4
c               if((icomp(i).eq.2 .or. icomp(i).eq.3) .and.
c     $              (ri2.lt.1 .or. ri3.lt.1))then
c                  write(69,*)'being put in 4 ... particle',i,icomp(i)
c                  write(69,*)bi1,bi2,bi3
c                  write(69,*)ri1,ri2,ri3
c                  stop
c                  endif
            endif

c 9/3/04: don't let point masses have their component change
            if(u(i).eq.0.d0) ic=icomp(i)

c            if(icomp(i).eq.2) stop
            
            if (ic.ne.icomp(i)) then
c     (particle assignment has changed)
c     adjust components
c     remove from old component

               if(am(i).gt.1d0 .and. myrank.eq.0) then
c               if(.true.) then
                  write(69,*)'i=',i,'m=',am(i),
     $                 'from',icomp(i),'to',ic
                  write(69,*)'bis=',bi1,bi2
                  write(69,*)'v2s=',v2i1, v2i2
                  write(69,*)'pots=',-(am1-am(i))/ri1,-(am2-am(i))/ri2
                  write(69,*)'rs=',ri1,ri2
               endif

               if (icomp(i).eq.1) then
c                  am1=am1-am(i)
                  if(i.eq.irhomax1) then
                     rhomax=0.d0
                     irhomax1=0
                     do j=1,n
                        if(j.ne.i .and. icomp(j).eq.1 .and.
     $                       rho(j).gt.rhomax) then
                           rhomax=rho(j)
                           irhomax1=j
                        endif
                     enddo
c                     write(69,*)'removing',i,'from comp 1; new=',irhomax1
                  endif
               else if (icomp(i).eq.2) then
c                  am2=am2-am(i)
                  if(i.eq.irhomax2) then
                     rhomax=0.d0
                     irhomax2=0
                     do j=1,n
                        if(j.ne.i .and. icomp(j).eq.2 .and.
     $                       rho(j).gt.rhomax) then
                           rhomax=rho(j)
                           irhomax2=j
                        endif
                     enddo
c                     write(69,*)'removing',i,'from comp 2; new=',irhomax2
                  endif
               else if (icomp(i).eq.3) then
c                  am3=am3-am(i)
                  if(i.eq.irhomax3) then
                     rhomax=0.d0
                     irhomax3=0
                     do j=1,n
                        if(j.ne.i .and. icomp(j).eq.3 .and.
     $                       rho(j).gt.rhomax) then
                           rhomax=rho(j)
                           irhomax3=j
                        endif
                     enddo
c                     write(69,*)'removing',i,'from comp 3; new=',irhomax3
                  endif
               endif

c     add to new component
               if (ic.eq.1) then
                  if(rho(i).gt.rho(irhomax1)) then
                     if(icomp(i).ne.4 .and. myrank.eq.0)
     $                    write(69,*)'in comp 1: replacing',irhomax1,
     $                    'with',i,'from comp',icomp(i)
                     irhomax1=i
                     if(icomp(i).ne.4) then
                        if(swapcase.ne.0 .and.
     $                       swapcase.ne.icomp(i)+ic .and. myrank.eq.0)then
                           write(69,*)'double swap???'
c                           stop
                        endif
                        swapcase=icomp(i)+ic
                        if(myrank.eq.0)write(69,*)'swapcase',i,icomp(i),ic
                     endif
                  endif
               else if (ic.eq.2) then
                  if(rho(i).gt.rho(irhomax2)) then
                     if(icomp(i).ne.4 .and. myrank.eq.0)
     $                    write(69,*)'in comp 2: replacing',irhomax2,
     $                    'with',i,'from comp',icomp(i)
                     irhomax2=i
                     if(icomp(i).ne.4) then
                        if(swapcase.ne.0 .and.
     $                       swapcase.ne.icomp(i)+ic .and. myrank.eq.0)then
                           write(69,*)'double swap???'
c                           stop
                        endif
                        swapcase=icomp(i)+ic
                        if(myrank.eq.0)write(69,*)'swapcase',i,icomp(i),ic
                     endif
                  endif
               else if (ic.eq.3) then
                  if(rho(i).gt.rho(irhomax3)) then
                     if(icomp(i).ne.4 .and. myrank.eq.0)
     $                    write(69,*)'in comp 3: replacing',irhomax3,
     $                    'with',i,'from comp',icomp(i)
                     irhomax3=i
                     if(icomp(i).ne.4) then
                        if(swapcase.ne.0 .and.
     $                       swapcase.ne.icomp(i)+ic .and. myrank.eq.0)then
                           write(69,*)'double swap???'
c                           stop
                        endif
                        swapcase=icomp(i)+ic
                        if(myrank.eq.0)write(69,*)'swapcase',i,icomp(i),ic
                     endif
                  endif
               endif

c               write(69,*)i,'is being moved from',icomp(i),'to',ic
c               write(69,*)x(i),y(i),z(i)

               icomp(i)=ic
               nchng=nchng+1
            endif
         enddo

         call calccom

         ncit=ncit+1
         
         if (db .and. myrank.eq.0) write (69,'(5g13.6)') ncit,nchng,am1,am2,am3

         if (ncit.ge.ncitmax .and. myrank.eq.0) then
            write(69,*) 'compbest: no convergence ???'
         endif

      enddo

      if(swapcase.ne.0) then
         if(swapcase.eq.3) then
            compa=1
            compb=2
            irhomaxtmp=irhomax1
            irhomax1=irhomax2
            irhomax2=irhomaxtmp
         elseif(swapcase.eq.4)then
            compa=1
            compb=3
            irhomaxtmp=irhomax1
            irhomax1=irhomax3
            irhomax3=irhomaxtmp
         elseif(swapcase.eq.5)then
            compa=2
            compb=3
            irhomaxtmp=irhomax2
            irhomax2=irhomax3
            irhomax3=irhomaxtmp
         else
            write(69,*)'swapcase problem'
            stop
         endif
         if(myrank.eq.0)write(69,*)'time to merge the stars at t=',t
         if(compa.eq.1) massa=am1last
         if(compa.eq.2) massa=am2last
         if(compa.eq.3) massa=am3last
         if(compb.eq.1) massb=am1last
         if(compb.eq.2) massb=am2last
         if(compb.eq.3) massb=am3last
         if(myrank.eq.0)write(69,*)'will merge components to recover...'
         if(myrank.eq.0)write(69,*)'last masses=',massa,massb
         if(massa.gt.massb)then
            if(myrank.eq.0)write(69,*)'switching all particles from',compb,'to',compa
            do i=1,n
               if(icomp(i).eq.compb)icomp(i)=compa
            enddo
         else
            if(myrank.eq.0)write(69,*)'switching all particles from',compa,'to',compb
            do i=1,n
               if(icomp(i).eq.compa)icomp(i)=compb
            enddo
         endif
         if(myrank.eq.0)write(69,*)'let us try again'
         checkifswapped=.true.
         firstt=.false.
         goto 123
      endif

      ncomp1=0
      ncomp2=0
      ncomp3=0
      ncomp4=0
      am4=0.d0
      do i=1,n
         if (icomp(i).eq.1) then
            ncomp1=ncomp1+1
         else if (icomp(i).eq.2) then
            ncomp2=ncomp2+1
         else if (icomp(i).eq.3) then
            ncomp3=ncomp3+1
         else if (icomp(i).eq.4) then
            am4=am4+am(i)
            ncomp4=ncomp4+1
         endif
      enddo

      numstars=min(ncomp1,1)+min(ncomp2,1)+min(ncomp3,1)

      if (ncomp1.eq.0 .and. ncomp2.eq.0 .and. ncomp3.eq.0) then
         if(myrank.eq.0)write(69,*)'warning: no stars found!'
         if(.not. firstt) then
            if(myrank.eq.0)then
               write(69,*)'trying a different initial guess:'
               write(69,*)'if stars are found, there is no guarantee that',
     $              'they will be the same component numbers as in the',
     $              'previous iteration'
            endif
            firstt=.true.
            goto 123
         else
            stop
         endif
      else
         firstt=.false.
      endif

      numstarslast=numstars
      am1last=am1
      am2last=am2
      am3last=am3

      if(myrank.eq.0)write(69,'(a,7a10)') '       ','mass   ','x   ','y   ','z   ',
     $     'vx   ','vy   ','vz   '

      if(myrank.eq.0) then
         if(am1.gt.0d0) then
            write(69,'(a,7g10.3)') 'star 1:',am1,x1,y1,z1,vx1,vy1,vz1
         endif
         if(am2.gt.0d0) then
            write(69,'(a,7g10.3)') 'star 2:',am2,x2,y2,z2,vx2,vy2,vz2
         endif
         if(am3.gt.0d0) then
            write(69,'(a,7g10.3)') 'star 3:',am3,x3,y3,z3,vx3,vy3,vz3
         endif
      endif

      if (pp) then

         if(am3.ge.am1 .and. am3.ge.am2) icompfocus=3
         if(am2.ge.am1 .and. am2.ge.am3) icompfocus=2
         if(am1.ge.am2 .and. am1.ge.am3) icompfocus=1

         rhomaxfocus=0.d0
         do i=1,n
            amenclosed(i)=0.d0
            if(icomp(i).eq.icompfocus) then
               if(rho(i).gt.rhomaxfocus) then
                  icenter=i
                  rhomaxfocus=rho(i)
               endif
            endif
         enddo

         spinx=0.d0
         spiny=0.d0
         spinz=0.d0
         do i=1,n
            if(icomp(i).eq.icompfocus) then
               xrel=x(i)-x(icenter)
               yrel=y(i)-y(icenter)
               zrel=z(i)-z(icenter)
               vxrel=vx(i)-vx(icenter)
               vyrel=vy(i)-vy(icenter)
               vzrel=vz(i)-vz(icenter)
               spinx=spinx+am(i)*(yrel*vzrel-zrel*vyrel)
               spiny=spiny+am(i)*(zrel*vxrel-xrel*vzrel)
               spinz=spinz+am(i)*(xrel*vyrel-yrel*vxrel)

               do j=1,n
                  if(icomp(j).eq.icompfocus .and. rho(j).lt.rho(i))then
                     amenclosed(j)=amenclosed(j)+am(i)
                  endif
               enddo
               amenclosed(i)=amenclosed(i)+0.5d0*am(i)

            endif
         enddo
         spintotal=(spinx**2+spiny**2+spinz**2)**0.5d0
         
         if(myrank.eq.0) then
            write(69,*)
            write(69,'(a,3g13.5)')'icompfocus=',icompfocus
            write(69,'(a,3g13.5)')'coordinates=',x(icenter),y(icenter),
     $           z(icenter)
            write(69,'(a,3g13.5)')'velocity vector=',vx(icenter),vy(icenter),
     $           vz(icenter)
            write(69,'(a,3g13.5)')'angular momentum vector=',spinx,spiny,spinz
            write(69,*)
            write (fname,102) nout
 102        format ('pp',i4.4,'.sph')
            open (13,file=fname)
         endif

         if(neos.ne.1) then
c     rachel, there is no need to do anything here because compbest3 is not used for our project
            print *,'this pp portion of combest3 currently works only'
            print *,'if the eos'
            print *,'is ideal gas + radiation pressure'
         endif

         do i=1,n
            rhocgs=rho(i)*munit/runit**3.d0
            ucgs=u(i)*gravconst*munit/runit
            call gettemperature(qconst*rhocgs/meanmolecular(i),
     $           -ucgs*rhocgs/arad,temperature)
            pgas=rhocgs*boltz*temperature/meanmolecular(i)
            prad=arad*temperature**4/3.d0
            betaofi=pgas/(pgas+prad)
            gamma=(32.d0-24.d0*betaofi-3.d0*betaofi**2) /
     $           (24.d0-21.d0*betaofi)
            rad=((x(i)-x(icenter))**2+
     $           (y(i)-y(icenter))**2+
     $           (z(i)-z(icenter))**2)**0.5d0
            dotproduct=spinx*(x(i)-x(icenter))+
     $           spiny*(y(i)-y(icenter))+
     $           spinz*(z(i)-z(icenter))
            if(rad.eq.0) then
               theta=3.14159265358979323846264d0/2.d0
               rcyl=rad*sin(theta)
               angularvelocity=0.d0
            else
               theta=acos(dotproduct/spintotal/rad)
               rcyl=rad*sin(theta)
               vxrel=vx(i)-vx(icenter)
               vyrel=vy(i)-vy(icenter)
               vzrel=vz(i)-vz(icenter)
               angularvelocity=( (vxrel**2+vyrel**2+vzrel**2)/
     $              rcyl**2)**0.5d0
            endif

            if(myrank.eq.0) then
               write (13,99) x(i),y(i),z(i),icomp(i),rho(i),temperature,
     $              rad,theta*180/3.14159265358979323846264d0,rcyl,
     $              amenclosed(i),pgas,prad,gamma,angularvelocity,grpot(i),
     $              u(i)
 99            format(3g13.5,i2,99g13.5)
            endif
         enddo
         if(myrank.eq.0) close (13)
      endif

      if(myrank.eq.0)write(69,*)'combest3:',numstars

      if(numstars.gt.2 .and. makeunboundstarcomp4) then
c     figure out which component can make it to infinity and move
c     that component to component 4, i.e. treat it like ejecta:
         mu=am1*(am2+am3)/(am1+am2+am3)
         x23=(am2*x2+am3*x3)/(am2+am3)
         y23=(am2*y2+am3*y3)/(am2+am3)
         z23=(am2*z2+am3*z3)/(am2+am3)
         vx23=(am2*vx2+am3*vx3)/(am2+am3)
         vy23=(am2*vy2+am3*vy3)/(am2+am3)
         vz23=(am2*vz2+am3*vz3)/(am2+am3)
         e1=0.5d0*mu*((vx1-vx23)**2+(vy1-vy23)**2+(vz1-vz23)**2)
     $        -am1*(am2+am3)/
     $        ((x1-x23)**2+(y1-y23)**2+(z1-z23)**2)**0.5d0

         mu=am2*(am1+am3)/(am2+am1+am3)
         x13=(am1*x1+am3*x3)/(am1+am3)
         y13=(am1*y1+am3*y3)/(am1+am3)
         z13=(am1*z1+am3*z3)/(am1+am3)
         vx13=(am1*vx1+am3*vx3)/(am1+am3)
         vy13=(am1*vy1+am3*vy3)/(am1+am3)
         vz13=(am1*vz1+am3*vz3)/(am1+am3)
         e2=0.5d0*mu*((vx2-vx13)**2+(vy2-vy13)**2+(vz2-vz13)**2)
     $        -am2*(am1+am3)/
     $        ((x2-x13)**2+(y2-y13)**2+(z2-z13)**2)**0.5d0

         mu=am3*(am1+am2)/(am3+am1+am2)
         x12=(am1*x1+am2*x2)/(am1+am2)
         y12=(am1*y1+am2*y2)/(am1+am2)
         z12=(am1*z1+am2*z2)/(am1+am2)
         vx12=(am1*vx1+am2*vx2)/(am1+am2)
         vy12=(am1*vy1+am2*vy2)/(am1+am2)
         vz12=(am1*vz1+am2*vz2)/(am1+am2)
         e3=0.5d0*mu*((vx3-vx12)**2+(vy3-vy12)**2+(vz3-vz12)**2)
     $        -am3*(am1+am2)/
     $        ((x3-x12)**2+(y3-y12)**2+(z3-z12)**2)**0.5d0

         if(e1.ge.max(e2,e3) .and. e1.gt.0.d0)then
            unboundcomp=1
            am1=0.d0
         else if(e2.ge.max(e1,e3) .and. e2.gt.0.d0)then
            unboundcomp=2
            am2=0.d0
         else if(e3.ge.max(e2,e1) .and. e3.gt.0.d0)then
            unboundcomp=3
            am3=0.d0
         else
            write(69,*)'trouble finding ejected star',e1,e2,e3
            stop
         endif

         if(myrank.eq.0)write(69,*)'ejected star is component',unboundcomp
         do i=1,ntot
            if(icomp(i).eq.unboundcomp) icomp(i)=4
         enddo

      endif
         
      return
      end
