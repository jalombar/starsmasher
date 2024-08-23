      subroutine multiequalmass
c     creates a star with varying equalmass values, spanning from +1 along the +x-axis to -1 along the -x-axis.
      include 'starsmasher.h'
      include 'mpif.h'
      integer numlines,i
      integer idumb,ip,ix,iy,iz
      real*8 anumden,rhotry,rhoex,rtry,rhomax,hc,xcm,ycm,zcm,amtot,
     $     ammin,ammax,xtry,ytry,ztry,ri,rhoi
      integer irtry
      real*8 amass,masscgs,radius
      real*8 tem(kdm),pres(kdm),
     &     rhoarray(kdm),uarray(kdm),rarray(kdm),
     $     rhoarray2(kdm),uarray2(kdm),
     $     muarray(kdm),muarray2(kdm)
      real*8 egsol
c     astronomical constants:
      parameter(egsol=1.9891d+33)
c     derived constants:
      real*8 integratednum
      common/splinestuff/rarray,uarray,muarray,rhoarray,
     $     uarray2,muarray2,rhoarray2,amass,radius,
     $     integratednum,maxmu,minmu,numlines
      real*8 utot2,wtota2
      integer ixmax,iymax,izmax,corepts
      double precision cellvolume,a1
      real*8 maxmu,minmu,drhodhi
      double precision utottest
      real*8 epot
      real*8 redge1,redge2
      real*8 hpguess,xacc,dxmax
      real*8 mci
      integer nmin
      real*8 amass1,amass2,hmax
      common/forcompbest/ amass1,amass2
      common/presarray/ pres,i
      common/hack/tem,redge1,masscgs,utot2,wtota2
      integer mylength,mygravlength, ierr
      integer comm_worker, irank
      common/gravworkers/comm_worker
      integer status(mpi_status_size)
      real*8 hmin
      integer numangle
      parameter (numangle=301)
      real*8 integral(numangle)
      real*8 rpreedge(numangle), rpost
      real*8 rprearray(kdm,numangle),rarray2(kdm,numangle)
      real*8 xcenter,ycenter,rcenter
      integer jangle0,ixmin,j,minix,maxix,miniy,maxiy,miniz,maxiz
      real*8 minequalmass,maxequalmass
      parameter (minequalmass=-0.75d0,maxequalmass=0d0)
      integer jangle
      real*8 cosangle,rightedgeoverradius
      parameter (rightedgeoverradius=1d0)
      
      call splinesetup

      idumb=-2391
      rhomax=rhoarray(1)
      if(myrank.eq.0)write(69,*)'maximum central density=',rhomax

      if(myrank.eq.0)write(69,*)'85 percent radius redge1=',redge1

      if(n.lt.0) then
            n=abs(n)*masscgs/egsol
      endif
      if(gflag.eq.0)then
         hc=0.5d0*radius*(1.8d0*nnopt/n)**(1.d0/3.d0) ! neighbor number is about 1.8*nnopt, roughly
      else
         hc=0.5d0*radius*(1.3d0*nnopt/n)**(1.d0/3.d0) !trying to better estimate the nearest neighbor number... about 1.3*nnopt
      endif

      redge2=radius-1.78d0*hc
      if(myrank.eq.0)write(69,*)'1.78h away from surface is redge2=',redge2

      redge=max(redge1,redge2)
      if(myrank.eq.0)write(69,*)'redge=max(redge1,redge2)=',redge

      if(myrank.eq.0)write(69,*)'hc=',hc
      
      nmin=max(156,int(12.d0*(0.5d0*redge/hc)**3.d0))
      if(myrank.eq.0)write(69,*) 'n should be at least',nmin

      n=max(n,nmin)

      if(myrank.eq.0) write(69,*) 'will try for n=',n

      jangle0=nint((numangle-1d0)/2)+1

      do jangle=1,numangle
!     in general, cos(angle)=xtry/rtry, where rtry is the spherical polar coordinate r
!     angle = 0 (x-axis) -> cos(angle)=1
!     angle = pi/2 -> cos(angle)=0
!     angle = pi (-x-axis) -> cos(angle)=-1
!         cos(angle)=2*(jangle-1d0)/(numangle-1d0)-1d0 runs from -1 (for j=1) to 1 (for j=numangle)
!     Solve for jangle:
!     jangle=nint((cos(angle)+1d0)*(numangle-1d0)/2)+1
!     For example, if cos(angle)=0, then jangle=jangle0=nint((numangle-1d0)/2)+1

         cosangle=2*(jangle-1d0)/(numangle-1d0)-1d0 !runs from -1 (for j=1) to 1 (for j=numangle)

         i=2
         integral(jangle)=0
         do while(rarray(i).lt.redge)
            xtry=cosangle*rarray(i)
            equalmass=(maxequalmass-minequalmass)*
     $           (min(xtry/radius,rightedgeoverradius) + 1)
     $           /(rightedgeoverradius+1)
     $           + minequalmass
            integral(jangle)=integral(jangle)+pi*(rarray(i)+rarray(i-1))**2*
     $           (rarray(i)-rarray(i-1))*
     $           (0.5d0*(rhoarray(i)+rhoarray(i-1))/rhomax)
     $           **equalmass
            i=i+1
         enddo
         if(myrank.eq.0)write(69,*)'integral(',jangle,')=',integral(jangle),4.d0/3.d0*pi*redge**3.d0
         
         rprearray(1,jangle)=rarray(1)
         rpreedge(jangle) = rprearray(1,jangle)
         do i = 2, numlines
            xtry=cosangle*rarray(i)
            equalmass=(maxequalmass-minequalmass)*
     $           (min(xtry/radius,rightedgeoverradius) + 1)
     $           /(rightedgeoverradius+1)
     $           + minequalmass
            rprearray(i,jangle) = rprearray(i-1,jangle) + 
     $           (rarray(i)-rarray(i-1))*
     $           (0.5d0*(rhoarray(i)+rhoarray(i-1))/rhomax)**equalmass*
     $           (rarray(i-1)/rprearray(i-1,jangle))**2
            if (rarray(i).le.redge) rpreedge(jangle) = rprearray(i,jangle)
         enddo
         call sph_spline(rprearray(1,jangle),rarray,numlines,1.d30,1.d30,rarray2(1,jangle))

      enddo

cccccccccccccccccccccccccccccccccccccccccccc
c     Call sph_spine in loop above here
cccccccccccccccccccccccccccccccccccccccccccc
      
c      mci=rhomax*4.d0/3.d0*pi*redge**3/n
c      mci=rhomax*integral(jangle0)/n








      if(myrank.eq.0) then
         write(69,*)'jangle0=',jangle0
c         write(69,*)'a central particle would have mass',mci   
      endif

      if(nnopt.le.0) then
         nnopt=max(nint(n*(hc/redge)**3.d0*8.d0),13)
      endif
      if(myrank.eq.0)write(69,*) 'using nnopt=',nnopt

c      if(4*mci.ge.amass) then
cc     if corepts=1 then the core point will be particle 1... first sph particle will be particle 2
c         corepts=1
c         if(myrank.eq.0)write(69,*)'we will use a core point'
c      else
c     if corepts=0 then there is no core particle and the first sph particle will be particle 1
         corepts=0
         if(myrank.eq.0)write(69,*)'we will not use a core point'
c      endif
      ip=corepts

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     make an hcp lattice
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      redge=radius
c      
c      if(amass.gt.50) then
c         maxi=1
c      else
c         maxi=20
c      endif
      
cc     do i=1,20              ! seems to work better when integrating u
cc     do i=1,1               ! seems to work better when integrating a
c      do i=1,maxi
c         hc=(nnopt/(8.d0*n))**(1.d0/3.d0)*redge
c         redge=radius-2.d0*hc
c         if(myrank.eq.0)write(69,*)'redge=',redge
c      enddo
      if(myrank.eq.0)write(69,*)'keeping particles up to a distance',radius-redge,
     $     'less than the full radius',radius

c     the fraction of particles at a region of density rhoex that
c     will be kept is (rhoex/rhomax)**equalmass.  so if the number
c     density of lattice points that will be tried is n, then we expect
c     the number of particles to be N=n*integrate[4*pi*r**2*
c     (rhoex(r)/rhomax)**equalmass,{r,0,redge}].  we can solve this for
c     n, and then use that the cell volume is 2/n (as there are two
c     particles per cell)

      if(myrank.eq.0) then
         write(69,*) '************STRETCHY HCP*************************'
         write(69,*) 'rarray(1)=',rarray(1)
      endif

c      if(equalmass.ne.0.d0) then
c         i=2
c         integral=0
c         do while(rarray(i).lt.redge)
c            integral=integral+pi*(rarray(i)+rarray(i-1))**2*
c     $           (rarray(i)-rarray(i-1))*
c     $           (0.5d0*(rhoarray(i)+rhoarray(i-1))/rhomax)
c     $           **equalmass
c            i=i+1
c         enddo
c         if(myrank.eq.0)write(69,*)'integral=',integral,4.d0/3.d0*pi*redge**3.d0
c      else
c         integral=4.d0/3.d0*pi*redge**3.d0
c         if(myrank.eq.0)write(69,*)'volume integral=',integral,4.d0/3.d0*pi*redge**3.d0
c      endif

      cellvolume=2.d0*integral(jangle0)/n
      a1=(cellvolume/2.d0**0.5d0)**(1.d0/3.d0)
      if(myrank.eq.0)write(69,*)'a1=',a1,cellvolume,integral(jangle0)
      
c     looking at figure 9(b) and page 18 of kittel
c     (a1 vector)=-0.5*a1*(x hat)-3^0.5/2*a1*(y hat)
c     (a2 vector)=a1*(x hat)
c     (a3 vector)=(8/3)^0.5*a1*(z hat)

      minix=9999999
      miniy=minix
      miniz=minix
      maxix=-minix
      maxiy=-minix
      maxiz=-minix
      
      ixmin=-int(rpreedge(1)/a1)-2 ! j=1 along the -x-axis
      ixmax=int(rpreedge(numangle)/a1)+2 ! j=numangle along the +x-axis
      iymax=-ixmin !int(rpreedge(j0)/(3.d0**0.5d0/2.d0*a1))+2
      izmax=-ixmin !int(rpreedge(j0)/(0.5d0*(8.d0/3.d0)**0.5d0*a1))+2

      if(myrank.eq.0) write(69,*)'ixmin,ixmax,iymin,iymax,izmin,izmax=',
     $     ixmin,ixmax,-iymax,iymax,-izmax,izmax 
      
      do ix=ixmin,ixmax
         do iy=-iymax,iymax
            do iz=-izmax,izmax
               xtry=(ix-0.5d0)*a1+mod(abs(iy),2)*0.5d0*a1
               ytry=iy*3.d0**0.5d0/2.d0*a1
     $              -(mod(abs(iz),2)-0.5d0)*1.d0/3.d0**0.5d0*a1
               ztry=(iz-0.5d0)*0.5d0*(8.d0/3.d0)**0.5d0*a1

               rtry=sqrt(xtry**2.d0+ytry**2.d0+ztry**2.d0)
!     in general, cos(angle)=xtry/rtry, where rtry is the spherical polar coordinate r
!     angle = 0 (x-axis) -> cos(angle)=1
!     angle = pi/2 -> cos(angle)=0
!     angle = pi (-x-axis) -> cos(angle)=-1
!         cos(angle)=2*(jangle-1d0)/(numangle-1d0)-1d0 runs from -1 (for j=1) to 1 (for j=numangle)
!     Solve for jangle:
!     jangle=nint((cos(angle)+1d0)*(numangle-1d0)/2)+1
!     For example, if cos(angle)=0, then jangle=jangle0=nint((numangle-1d0)/2)+1
               equalmass=(maxequalmass-minequalmass)*
     $              (min(xtry/radius,rightedgeoverradius) + 1)
     $              /(rightedgeoverradius+1)
     $              + minequalmass
               cosangle=xtry/rtry
               jangle=nint((cosangle+1d0)*(numangle-1d0)/2)+1
               
               if(rtry.lt.rpreedge(jangle)) then
                  call sph_splint(rprearray(1,jangle),rarray,rarray2(1,jangle),numlines,
     $                 rtry,rpost)
c               if(rpost.lt.redge) then
c                  call sph_splint(rarray,rhoarray,rhoarray2,numlines,
c     $                 rtry,rhoex)
c                  rhotry=rhomax**equalmass*ran1(idumb)      
c                  if (rhotry.le.rhoex**equalmass) then

c     The number densty n of particles is proportional to rho**equalmass,
c     so that means the cell spacing a1 should scale like rho**(-equalmass/3).
c     That is a1 = a1center*(rhoex/rhomax)**(-equalmass/3).

c     (particle is accepted)    
                  ip=ip+1                 
                  x(ip)=xtry*rpost/rtry
                  y(ip)=ytry*rpost/rtry 
                  z(ip)=ztry*rpost/rtry
                  if(rtry.le.a1 .and. myrank.eq.0) then
                     write(69,'(4i6,4e12.4)')ip,ix,iy,iz,
     $                    x(ip),y(ip),z(ip),rpost
                  endif
                  minix=min(minix,ix)
                  miniy=min(miniy,iy)
                  miniz=min(miniz,iz)
                  maxix=max(maxix,ix)
                  maxiy=max(maxiy,iy)
                  maxiz=max(maxiz,iz)
               endif
            enddo
         enddo
      enddo

      if(myrank.eq.0) then
         write(69,*)'minix,maxix,miniy,maxiy,miniz,maxiz=',
     $        minix,maxix,miniy,maxiy,miniz,maxiz
         write(69,*) '************DONE WITH STRETCHY HCP***************'
      endif


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     (note that n will be slightly changed!)
      n=ip-corepts
      if(myrank.eq.0) write (69,*) 'parent: n=',n
      amass1=n
      amass2=0
      if (n+corepts.gt.nmax) stop 'parent: n>nmax ???'
c     assign particle masses (to represent density):
      amtot=0.d0
      xcm=0.d0
      ycm=0.d0
      zcm=0.d0

c     start following loop at the first sph particle (with index 2, not 1, if there is a core particle)
      do i=1+corepts,n+corepts
         ri=sqrt(x(i)**2+y(i)**2+z(i)**2)
         call sph_splint(rarray,rhoarray,rhoarray2,numlines,ri,rhoi)
         if(rhoi.le.0.d0) then
            if(myrank.eq.0)write(69,*)'warning: rho(',i,')<=0 at r=',ri,'???'
         endif
c         equalmass=x(i)/ri
         equalmass=(maxequalmass-minequalmass)*
     $        (min(x(i)/radius,rightedgeoverradius) + 1)
     $        /(rightedgeoverradius+1)
     $        + minequalmass

         am(i)=integral(jangle0)/n*rhoi**(1.d0-equalmass)*rhomax**equalmass

c         write(300+myrank,*)i,am(i),integral(jangle0),n,rhoi,equalmass,rhomax,maxequalmass,minequalmass,x(i),radius

         xcm=xcm+am(i)*x(i)
         ycm=ycm+am(i)*y(i)
         zcm=zcm+am(i)*z(i)
         amtot=amtot+am(i)
         call sph_splint(rarray,uarray,uarray2,numlines,ri,u(i))
         call sph_splint(rarray,muarray,muarray2,numlines,ri,
     $        meanmolecular(i))
         anumden=rhoi/am(i)
c     the actual number of neighbors is closer to 1.9*nnopt
         hp(i)=(3.d0/32.d0/3.1415926535897932384626d0*
     $        1.9d0*nnopt/anumden)**(1.d0/3.d0) + hfloor

      enddo
      
      xcm=xcm/amtot
      ycm=ycm/amtot
      zcm=zcm/amtot
      if(myrank.eq.0) then
         write(69,'(a,3g13.4)') 'center of mass at',xcm,ycm,zcm
         write(69,'(a,g11.4,a,g11.4,a)')
     $        'parent: total mass sph fluid was m=',amtot
     $        ,
     $        ', total mass of star will equal', amass
      endif

      hmin=1.d30
      hmax=0.d0
      do i=1+corepts,n+corepts
         hmin=min(hmin,hp(i))
         hmax=max(hmax,hp(i))
      enddo
      if(myrank.eq.0)then
         write(69,*)'hmin=',hmin
         write(69,*)'hmax=',hmax
      endif

c     above here, n=number of sph particles
      n=n+corepts
c     below here, n=total number of particles

      ntot=n
      if(myrank.eq.0) write(69,*) 'parent: n=',n,' ntot=', ntot
      if (n.gt.nmax) stop 'parent: n>nmax ???'
      if(ntot.gt.nmax) then
         if(myrank.eq.0)write(69,*) 'n is too large'
         stop
      endif

      if(corepts.gt.0) then
         x(1)=0.d0
         y(1)=0.d0
         z(1)=0.d0
         vx(1)=0.d0
         vy(1)=0.d0
         vz(1)=0.d0
         am(1)=amass-amtot
         if(myrank.eq.0)write(69,*) 'mass of core = ',am(1),'msun'
         if(hco.gt.0.d0) then
            hp(1)=hco
         else
            hp(1)=hmin
         endif
         if(myrank.eq.0)write(69,*)'hp(core mass)',hp(1)
         cc(1)=int(2.d0*log(1.35d0*a1*
     $        (integral(jangle0)/(4.d0/3.d0*pi*redge**3))**(1.d0/3.d0))
     $        /log(2.d0)+15.d0)
         cc(1)=max(cc(1),0)
         cc(1)=min(cc(1),ntypes-1)
         u(1)=0.d0
         call sph_splint(rarray,muarray,muarray2,numlines,0.d0,
     $        meanmolecular(1))
         if (am(1).le.0) then
            if(myrank.eq.0)write(69,*)'core mass < 0'
            stop
         endif
      else
         ammin=1.d30
         ammax=0.d0
         do i=1,n
            am(i)=am(i)/amtot*amass
            if(am(i).gt.0.01d0*amass) then
               if(myrank.eq.0)write(69,*)'warning: particle',i,
     $              'has mass',am(i)
               if(myrank.eq.0)write(69,*)'x,y,z=',x(i),y(i),z(i)
               if(myrank.eq.0)write(69,*)'r,rho=',ri,rhoi
            endif
            ammin=min(ammin,am(i))
            ammax=max(ammax,am(i))
            rtry=sqrt(x(i)**2+y(i)**2+z(i)**2)
            call sph_splint(rarray,rhoarray,rhoarray2,numlines,rtry,
     $           rhoex)
            anumden=rhoex/am(i)
            hp(i)=(3.d0/32.d0/3.1415926535897932384626d0*
     $           1.9d0*nnopt/anumden)**(1.d0/3.d0) + hfloor

         enddo

     
         if(myrank.eq.0)then
            write (69,*) 'parent: min particle mass=',ammin,ammin*munit/egsol
            write (69,*) 'parent: max particle mass=',ammax,ammax*munit/egsol
         endif
      endif

c     because the particles were distributed uniformly in the sphere, the
c     the number density is 3*n/(4*pi*radius**3) and we need to choose
c     the smoothing length hc such that 
c     4*pi*(2*hc)**3/3 * (number density)=nnopt.
c     this gives 8*hc**3 * n/radius**3= nnopt, or:
c     hc=(nnopt/(8.d0*n))**(1.d0/3.d0)*redge
c     (should give a number of nearest neighbors close to nnopt)
         
      call stride_setup
      if(ntot.gt.nmax) then
         if(myrank.eq.0)write(69,*) 'n is too large'
         stop
      endif

      if(myrank.eq.0)write(69,*)'corepts=',corepts,ntot

      if(myrank.eq.0)write(69,*) masscgs/egsol,'solar masses has cc=',
     $     int(10000*masscgs/egsol)
      do i=1,ntot
         cc(i)=int(10000*masscgs/egsol)
      enddo

      if(corepts.eq.0 .and. treloff.le.0.d0) then
         if(myrank.eq.0) write(69,*)'will try to get correct u and w...'

         utottest=0.d0
         do i=1,n
            utottest=utottest+u(i)*am(i)
         enddo
         if(myrank.eq.0) write (69,'(a,g11.4,a,g11.4,a)')
     $        'parent: total internal energy of star was',utottest
     $        ,
     $        ', which should equal', utot2/(gravconst*munit**2/runit),
     $        '... renormalizing...'
         do i=1,n
            u(i)=u(i)*utot2/utottest/(gravconst*munit**2/runit)
         enddo

c     do loop to get total gravitational potential energy:
         if(ngr.ne.0)then
            call gravquant      ! must call for inital set up
            call rho_and_h
            call gravforce
            if(myrank.lt.ngravprocs)then
               if(nusegpus.eq.1)then
                  call lasthalf_grav_forces(ntot, gx, gy, gz, grpot)
               else
                  call get_gravity_using_cpus
               endif
               if(ngravprocs.gt.1) then
                  mygravlength=ngrav_upper-ngrav_lower+1
                  if(myrank.ne.0)then
                     call mpi_gatherv(grpot(ngrav_lower), mygravlength, mpi_double_precision,
     $                    grpot, gravrecvcounts, gravdispls, mpi_double_precision, 0,
     $                    comm_worker, ierr)
                  else
                     call mpi_gatherv(mpi_in_place, mygravlength, mpi_double_precision,
     $                    grpot, gravrecvcounts, gravdispls, mpi_double_precision, 0,
     $                    comm_worker, ierr)
                  endif
               endif
            endif

            if(myrank.eq.0) then
               epot=0.d0
               if(ngr.ne.0)then
                  do i=1,ntot
                     epot=epot+am(i)*grpot(i)
                  enddo
                  epot=0.5d0*epot
               endif
               do irank=1,nprocs-1
                  call mpi_send(epot, 1, mpi_double_precision,
     $                 irank,irank,mpi_comm_world,ierr)
               enddo
            else
               call mpi_recv(epot, 1, mpi_double_precision,
     $              0,myrank,mpi_comm_world,status,ierr)
            endif
            
            if(myrank.eq.0)write (69,'(a,g11.4,a,g11.4,a)')
     $           'parent: gravitational potential energy of star was',epot
     $           ,
     $           ', which should equal', wtota2/(gravconst*munit**2/runit),
     $           '... rescaling positions...'
            do i=1,n
               x(i)=x(i)*epot/(wtota2/(gravconst*munit**2/runit))
               y(i)=y(i)*epot/(wtota2/(gravconst*munit**2/runit))
               z(i)=z(i)*epot/(wtota2/(gravconst*munit**2/runit))
               hp(i)=hp(i)*epot/(wtota2/(gravconst*munit**2/runit))
               ri=sqrt(x(i)**2+y(i)**2+z(i)**2)
               call sph_splint(rarray,muarray,muarray2,numlines,ri,
     $              meanmolecular(i))
               meanmolecular(i)=min(maxmu,meanmolecular(i))
               meanmolecular(i)=max(minmu,meanmolecular(i))
               if(meanmolecular(i).lt.1d-25 .or.
     $              meanmolecular(i).gt.1.d-23) then
                  if(myrank.eq.0)write(69,*)'mean molecular problem...',i,meanmolecular(i)
                  stop
               endif
            enddo
         endif
      endif

      if(myrank.eq.0 .and. corepts.gt.0)
     $     write(69,*)'properties of core:',x(1),y(1),z(1),vx(1),
     $     vy(1),vz(1),am(1),hp(1),cc(1),u(1)

      if(myrank.eq.0)write(69,*) 'radius   density   u'
      do i=1,numlines,numlines/10
         write(*,*) rarray(i),rhoarray(i),uarray(i)
      enddo
      
      if(myrank.eq.0)then
         write(69,*)'rarray(1) and rarray(2)=',rarray(1),rarray(2)
         write(69,*)'rhoarray(1) and rhoarray(2)=',
     $     rhoarray(1),rhoarray(2)
         write(69,*)'rarray(numlines) =',rarray(numlines)
         write(69,*)'rhoarray(numlines) =',rhoarray(numlines)
         write(69,*)'near origin:'
         write(69,*)' r    rho'
      endif
      do irtry=0,nint(1000*rarray(2)),nint(100*rarray(2)+0.5)
         rtry=irtry/1000.d0
         call sph_splint(rarray,rhoarray,rhoarray2,numlines,rtry,rhoi)
         if(myrank.eq.0)write(69,*) rtry,rhoi
      enddo
      if(myrank.eq.0)then
         write(69,*)'near surface:'
         write(69,*)' r    rho'
      endif
      do irtry=nint(50*radius),nint(100*radius),nint(5*radius+0.5)
         rtry=min(irtry/100.d0,radius)
         call sph_splint(rarray,rhoarray,rhoarray2,numlines,rtry,rhoi)
         if(myrank.eq.0)write(69,*) rtry,rhoi
      enddo

c     assign velocities (all zero)
      do i=1,ntot
         vx(i)=0.d0
         vy(i)=0.d0
         vz(i)=0.d0
         vxdot(i)=0.d0
         vydot(i)=0.d0
         vzdot(i)=0.d0
         udot(i)=0.d0
      enddo

c     prepare leap-frog scheme for first iteration:
      call lfstart

      hmin=1.d30
      hmax=0.d0
      do i=1+corepts,ntot
         hmin=min(hmin,hp(i))
         hmax=max(hmax,hp(i))
      enddo

c      if(corepts.gt.0) then
c         hp(1)=0.5d0*hmin
c         if(myrank.eq.0)write(69,*)'hp(core mass)=',hp(1)
c      endif
      if(myrank.eq.0)then
         write(69,*)'hmin=',hmin
         write(69,*)'hmax=',hmax
         write(69,*)'parent: tf=',tf,dtout
         write(69,*)'parent: exiting parent'
      endif

      return

      end
