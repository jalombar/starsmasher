      subroutine parent
c     creates a star from the data file yrec output
      include 'starsmasher.h'
      include 'mpif.h'
      integer numlines,i
      integer idumb,ip,ix,iy,iz
      real*8 anumden,rhotry,rhoex,rtry,rhomax,hc,xcm,ycm,zcm,amtot,
     $     ammin,ammax,xtry,ytry,ztry,ri,rhoi,ran1
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
      real*8 integral
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
      integer mygravlength, ierr
      integer comm_worker, irank
      common/gravworkers/comm_worker
      integer status(mpi_status_size)
      real*8 hmin

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

      redge2=radius-3.d0*hc
      if(myrank.eq.0)write(69,*)'3h away from surface is redge2=',redge2

      redge=max(redge1,redge2)
      if(myrank.eq.0)write(69,*)'redge=max(redge1,redge2)=',redge

      if(myrank.eq.0)write(69,*)'hc=',hc
      
      nmin=max(156,int(12.d0*(0.5d0*redge/hc)**3.d0))
      if(myrank.eq.0)write(69,*) 'n should be at least',nmin

      n=max(n,nmin)

      if(myrank.eq.0)write(69,*) 'will try for n=',n

      mci=rhomax*4.d0/3.d0*pi*redge**3/n
      if(myrank.eq.0)write(69,*)'a central particle would have mass',mci   

      if(nnopt.le.0) then
         nnopt=max(nint(n*(hc/redge)**3.d0*8.d0),13)
      endif
      if(myrank.eq.0)write(69,*) 'using nnopt=',nnopt

      if(4*mci.ge.amass) then
c     if corepts=1 then the core point will be particle 1... first sph particle will be particle 2
         corepts=1
         if(myrank.eq.0)write(69,*)'we will use a core point'
      else
c     if corepts=0 then there is no core particle and the first sph particle will be particle 1
         corepts=0
         if(myrank.eq.0)write(69,*)'we will not use a core point'
      endif
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
c     the number of particles to be n=n*integrate[4*pi*r**2*
c     (rhoex(r)/rhomax)**equalmass,{r,0,redge}].  we can solve this for
c     n, and then use that the cell volume is 2/n (as there are two
c     particles per cell)
      
      if(equalmass.ne.0.d0) then
         i=2
         integral=0
         do while(rarray(i).lt.redge)
            integral=integral+pi*(rarray(i)+rarray(i-1))**2*
     $           (rarray(i)-rarray(i-1))*
     $           (0.5d0*(rhoarray(i)+rhoarray(i-1))/rhomax)
     $           **equalmass
            i=i+1
         enddo
         if(myrank.eq.0)write(69,*)'integral=',integral,4.d0/3.d0*pi*redge**3.d0
      else
         integral=4.d0/3.d0*pi*redge**3.d0
         if(myrank.eq.0)write(69,*)'volume integral=',integral,4.d0/3.d0*pi*redge**3.d0
      endif

      cellvolume=2.d0*integral/n
      a1=(cellvolume/2.d0**0.5d0)**(1.d0/3.d0)
      if(myrank.eq.0)write(69,*)'a1=',a1,cellvolume,integral
      
c     looking at figure 9(b) and page 18 of kittel
c     (a1 vector)=-0.5*a1*(x hat)-3^0.5/2*a1*(y hat)
c     (a2 vector)=a1*(x hat)
c     (a3 vector)=(8/3)^0.5*a1*(z hat)
      ixmax=int(redge/a1)+2
      iymax=int(redge/(3.d0**0.5d0/2.d0*a1))+2
      izmax=int(redge/(0.5d0*(8.d0/3.d0)**0.5d0*a1))+2
      do ix=-ixmax,ixmax
         do iy=-iymax,iymax
            do iz=-izmax,izmax
               xtry=(ix-0.5d0)*a1+mod(abs(iy),2)*0.5d0*a1
               ytry=iy*3.d0**0.5d0/2.d0*a1
     $              -(mod(abs(iz),2)-0.5d0)*1.d0/3.d0**0.5d0*a1
               ztry=(iz-0.5d0)*0.5d0*(8.d0/3.d0)**0.5d0*a1
               rtry=sqrt(xtry**2.d0+ytry**2.d0+ztry**2.d0)
               if(rtry.lt.redge) then
                  call sph_splint(rarray,rhoarray,rhoarray2,numlines,
     $                 rtry,rhoex)
                  rhotry=rhomax**equalmass*ran1(idumb)      
                  if (rhotry.le.rhoex**equalmass) then
c     (particle is accepted)    
                     ip=ip+1                 
                     if(rtry.le.a1 .and. myrank.eq.0) then
                        write(69,'(4i6,4e12.4)')ip,ix,iy,iz,
     $                       xtry,ytry,ztry,rtry
                     endif
                     x(ip)=xtry        
                     y(ip)=ytry            
                     z(ip)=ztry
                  endif
               endif
            enddo
         enddo
      enddo

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
         am(i)=amass/n*(integral*rhoi/amass)**(1.d0-equalmass)
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
         hp(1)=hmin
         if(myrank.eq.0)write(69,*)'hp(core mass)',hp(1)
         cc(1)=int(2.d0*log(1.35d0*a1*
     $        (integral/(4.d0/3.d0*pi*redge**3))**(1.d0/3.d0))
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
      subroutine splinesetup
c     Read in a stellar evolution code file
      include 'starsmasher.h'
      integer numlines,i,j
      real*8 amass,radiuscgs,masscgs,radius
      real*8 integratedmass2
      integer*4 io
      parameter(io=13)
      real*8 tem(kdm),pres(kdm),
     &     rhoarray(kdm),uarray(kdm),rarray(kdm),
     $     rhoarray2(kdm),uarray2(kdm),xm(kdm),
     $     muarray(kdm),muarray2(kdm)
      real*8 egsol,solrad
c     astronomical constants:
      parameter(egsol=1.9891d+33,solrad=6.9599d10)
c     derived constants:
      integer ndim
      parameter(ndim=9)
      real*8 xx(ndim,kdm)
      integer lz(ndim),ln(ndim)
      real*8 integratednum,integratednum2
      common/splinestuff/rarray,uarray,muarray,rhoarray,
     $     uarray2,muarray2,rhoarray2,amass,radius,
     $     integratednum,maxmu,minmu,numlines
      real*8 utot2,wtota2,sum
      real*8 maxmu,minmu
      real*8 zeroin
      real*8 redge1
      real*8 amass1,amass2
      common/forcompbest/ amass1,amass2
      real*8 temperaturefunction
      external temperaturefunction
      real*8 ufunction,useeostable,temupperlimit,uupperlimit
      external ufunction
      common/presarray/ pres,i
      common/hack/tem,redge1,masscgs,utot2,wtota2
      integer k
      integer model_number, num_zones, zone
      real*8 initial_mass,initial_z,              
     $     star_age,time_step,Teff,photosphere_L,photosphere_r,         
     $     center_eta,center_h1,center_he3,center_he4,center_c12,       
     $     center_n14,center_o16,center_ne20,star_mass,star_mdot,       
     $     star_mass_h1,star_mass_he3,star_mass_he4,star_mass_c12,      
     $     star_mass_n14,star_mass_o16,star_mass_ne20,he_core_mass,     
     $     c_core_mass,o_core_mass,si_core_mass,fe_core_mass,           
     $     neutron_rich_core_mass                 
      real*8 qq, logR, logRho, logT, logP, logPgas,
     $     luminosity, x_mass_fraction_H,y_mass_fraction_He,            
     $     z_mass_fraction_metals,gamma1,opacity,pp,cno,gradr,gradT,    
     $     grada,actual_gradT,total_energy,total_energy_integral,       
     $     scale_height,mu,dummy
      integer ii,one,maxnumcol,col
      parameter(maxnumcol=99)
      integer in(maxnumcol),numcol
      character*41 headers(maxnumcol)
      integer imass,ilogR,ilogT,
     $     ilogRho,ilogP,ix_mass_fraction_H,iy_mass_fraction_He,
     $     iz_mass_fraction_metals      
      real*8 mesadata(2:maxnumcol)

      if(myrank.eq.0) then
         write(69,*) 'splinesetup: stellarevolutioncodetype=',
     $        stellarevolutioncodetype
         write(69,*) '             about to read file ',trim(profilefile)
      endif
 10   open(io,file=profilefile,status='old')
      if(stellarevolutioncodetype.eq.0) then
c     This is a *.s2mm file made by TWIN on starsmasher.allegheny.edu
c     h, he4 , c12, n14, o16, ne20
         ln(1)=0
         ln(2)=2
         ln(3)=6
         ln(4)=7
         ln(5)=8
         ln(6)=10
         ln(7)=12
         ln(8)=14
         ln(9)=30
         
         lz(1)=1
         lz(2)=2
         lz(3)=6
         lz(4)=7
         lz(5)=8
         lz(6)=10
         lz(7)=12
         lz(8)=14
         lz(9)=26
         
c     get profiles:
         i=0
         do k=1,kdm
            i=i+1
            read(io,*, end=21) xm(i),rarray(i),pres(i),rhoarray(i),
     $           (xx(j,i),j=1,ndim)
            if(xm(1).eq.1)then
               if(myrank.eq.0)then
                  write(69,*)'Uuuh.. This seems like a MESA file..'
                  write(69,*)'Will reset the stellarevolutioncodetype=1'
               endif
               stellarevolutioncodetype=1
               close (io)
               goto 10
            endif
            if(i.gt.1 .and. rarray(i).le.rarray(i-1)) then
               if(myrank.eq.0) write(69,*)'ignoring shell',k
               i=i-1
            endif
         enddo
         if(myrank.eq.0)write(69,*)'need to increase kdm'
         stop
 21      close(io)
         numlines=i-1
         xm(numlines+1)=xm(numlines)
         rarray(numlines+1)=rarray(numlines)
         pres(numlines+1)=0.d0
         rhoarray(numlines+1)=0.d0
         do j=1,ndim
            xx(j,numlines+1)=xx(j,numlines)
         enddo
         if(myrank.eq.0)write(69,*)'number of shells=',numlines

         maxmu=0.d0
         minmu=1.d30
         do i=1,numlines
 123        continue
            muarray(i)=0.d0
            sum=0.d0
            do j=1,ndim  
               muarray(i)=muarray(i)+xx(j,i)
     $              *dble(1+lz(j))/dble(ln(j)+lz(j))
               sum=sum+xx(j,i)
            enddo
            if(dabs(sum-1.d0).gt.1.d-8) then
               if(myrank.eq.0)
     $              write(69,*) 'problem with abundances',sum-1.d0,i
               if(i.gt.0.9*numlines) then
                  if(myrank.eq.0)
     $                 write(69,*)
     $                 'solving by using previous shell values'
                  do j=1,ndim
                     xx(j,i)=xx(j,i-1)
                  enddo
                  goto 123
               else
                  stop
               endif
            endif
            muarray(i)=1.d0/muarray(i)*1.67262158d-24
            maxmu=max(maxmu,muarray(i))
            minmu=min(minmu,muarray(i))
         enddo
      else if(stellarevolutioncodetype.eq.1) then
c     profilefile comes from MESA
         read(io,*,err=22) one  ! It's safe to ignore a compiler warning here

         if (one.ne.1)then
 22         if(myrank.eq.0)then ! It's safe to ignore a compiler warning here 
               write(69,*)'Uuuh.. This doesn''t seem like a MESA file..'
               write(69,*)'We will reset the stellarevolutioncodetype=0'
            endif
            stellarevolutioncodetype=0
            close (io)
            goto 10
         endif

         read(io,*)             ! names of scalar variables associated with those integers                   
         read(io,*)             ! read past all the global data... we'll read this in later
         read(io,*)             ! blank line
         
         






         do i=1,maxnumcol
            in(i)=0
         enddo
         
         read(io,*,err=99) (in(i),i=1,maxnumcol)
         print *,'should never see this'
         
 99      continue
         numcol=0
         do while(in(numcol+1).ne.0)
            numcol=numcol+1
         enddo
         close(io)
         if(myrank.eq.0) then
            write(69,*) 'number of columns in MESA file=',numcol
         endif
         open(io,file=profilefile,status='old')
         read(io,*) one         ! list of integers                   
         if(one.ne.1) then
            write(69,*) 'wait... I thought this was a MESA file...'
            stop
         endif
         read(io,*)             ! names of scalar variables associated with those integers                   
         read(io,*) model_number,num_zones,initial_mass,initial_z,           
     $        star_age,time_step,Teff,photosphere_L,photosphere_r,           
     $        center_eta,center_h1,center_he3,center_he4,center_c12,         
     $        center_n14,center_o16,center_ne20,star_mass,star_mdot,         
     $        star_mass_h1,star_mass_he3,star_mass_he4,star_mass_c12,        
     $        star_mass_n14,star_mass_o16,star_mass_ne20,he_core_mass,       
     $        c_core_mass,o_core_mass,si_core_mass,fe_core_mass,             
     $        neutron_rich_core_mass
         read(io,*)             ! blank line
         read(io,*)             ! list of integers                   

         read(io,'(99a)')(headers(i), i=1,numcol)

         do i=1,numcol
            headers(i)=adjustl(headers(i))
            if(myrank.eq.0) write(69,*) i,headers(i)
            if(trim(headers(i)).eq.'mass') then
               imass=i
            elseif(trim(headers(i)).eq.'logR') then
               ilogR=i
            elseif(trim(headers(i)).eq.'logT') then
               ilogT=i
            elseif(trim(headers(i)).eq.'logRho') then
               ilogRho=i
            elseif(trim(headers(i)).eq.'logP') then
               ilogP=i
            elseif(trim(headers(i)).eq.'x_mass_fraction_H' .or.
     $             trim(headers(i)).eq.'x') then
               ix_mass_fraction_H=i
            elseif(trim(headers(i)).eq.'y_mass_fraction_He' .or.
     $             trim(headers(i)).eq.'y') then
               iy_mass_fraction_He=i
            elseif(trim(headers(i)).eq.'z_mass_fraction_metals' .or.
     $             trim(headers(i)).eq.'z') then
               iz_mass_fraction_metals=i
            endif
         enddo
         
         if(myrank.eq.0) then
            write(69,*)'Relevant data in these columns of MESA file:'
            write(69,*) imass,': mass'
            write(69,*) ilogR,': log_10 radius'
            write(69,*) ilogT,': log_10 temperature'
            write(69,*) ilogRho,': log_10 density'
            write(69,*) ilogP,' : log_10 pressure'
            write(69,*) ix_mass_fraction_H,' : X'
            write(69,*) iy_mass_fraction_He,' : Y'
            write(69,*) iz_mass_fraction_metals,' : Z'
         endif


         if(myrank.eq.0) then
            write(69,*)'model_number=',model_number            
            write(69,*)'num_zones=',num_zones                  
         endif
         maxmu=0.d0                  
         minmu=1.d30                 
         i=num_zones                 
         do k=1,num_zones            
c            read(io,*) zone, xm(i), logR, logT, logRho, logP,
c     $           x_mass_fraction_H,y_mass_fraction_He,             
c     $           z_mass_fraction_metals

c     This is for the format of files that Pablo sent
c            read(io,*) zone, logT, logRho, logP, logR, (dummy, ii=1,29),
c     $           x_mass_fraction_H,y_mass_fraction_He,             
c     $           z_mass_fraction_metals, (dummy, ii=1,26),xm(i)

            read(io,*) zone, (mesadata(col), col=2,numcol)
            xm(i)= mesadata(imass)
            logR=mesadata(ilogR)
            logT=mesadata(ilogT)
            logRho=mesadata(ilogRho)
            logP=mesadata(ilogP)
            x_mass_fraction_H=mesadata(ix_mass_fraction_H)
            y_mass_fraction_He=mesadata(iy_mass_fraction_He)
            z_mass_fraction_metals=mesadata(iz_mass_fraction_metals)

            if(  abs(x_mass_fraction_H  + y_mass_fraction_He +
     $           z_mass_fraction_metals - 1d0) .gt. 1.d-8) then
               write(69,*)' X + Y + Z - 1 should equal 0 but equals',
     $              x_mass_fraction_H  + y_mass_fraction_He +
     $              z_mass_fraction_metals - 1d0
               stop
            endif

            xm(i)=xm(i)*1.9892d33

c            read(io,*) zone, qq, logR, logRho, logT, logP, logPgas,             
c     $           luminosity,x_mass_fraction_H,y_mass_fraction_He,             
c     $           z_mass_fraction_metals,gamma1,opacity,pp,cno,gradr,
c     $           gradT,grada,actual_gradT,total_energy,
c     $           total_energy_integral,scale_height,mu           
c            xm(i)=star_mass*qq*1.9892d33     

            if(k.eq.1 .and. myrank.eq.0) then
               write(69,*) 'Surface x,y,z=',x_mass_fraction_H,
     $              y_mass_fraction_He,z_mass_fraction_metals
            endif

            rarray(i)=10d0**logR*6.9598d10
            if(i.lt.num_zones .and. rarray(i+1).le.rarray(i)) then             
               write(69,*)'ignoring shell',k
               stop      
            endif        
            rhoarray(i)=10d0**logRho       
            pres(i)=10d0**logP             
            muarray(i)=1/(2*x_mass_fraction_H+0.75d0*
     $              y_mass_fraction_He+0.5d0*z_mass_fraction_metals)*1.67262158d-24   
            maxmu=max(maxmu,muarray(i))    
            minmu=min(minmu,muarray(i))    
            if(abs(x_mass_fraction_H + y_mass_fraction_He +      
     $           z_mass_fraction_metals -1d0) .gt. 1d-15) then  
               print *, x_mass_fraction_H + y_mass_fraction_He +               
     $              z_mass_fraction_metals -1d0,'mass problem in zone',k          
               stop      
            endif
            i=i-1        
         enddo           
         close(io)       
         numlines=num_zones                
         xm(numlines+1)=xm(numlines)       
         rarray(numlines+1)=rarray(numlines)                 
         pres(numlines+1)=0.d0             
         rhoarray(numlines+1)=0.d0         
         if(myrank.eq.0) then
            write(69,*)'number of shells=',numlines              
            write(69,'(9a13)')'mass','radius','pressure','density','mu'
            do i=1,numlines 
               write(69,'(9g13.5)')xm(i),rarray(i),pres(i),rhoarray(i),
     $              muarray(i)    
            enddo      
         endif
         
      endif   
      if(myrank.eq.0)
     $     write(69,*)'done reading ',trim(profilefile)
      
c     this routine matches the pressure and
c     density profiles from the stellar evolution code.  this way, hydrostatic
c     equilbrium dp/dr=-g*rho is maintained even though the stellar evolution
c     code's equation of state is different than sph's equation of state. the
c     code below solves for the temperature profile that
c     is necessary to give the desired pressure and density profiles.  then it
c     uses this temperature to get the internal energy assuming an eos that is
c     ideal gas plus radiation pressure (if neos=1).  if the stellar model used a
c     different eos of state, the temperature and internal energy profiles in
c     the sph code will be a bit different than that from the stellar
c     evolution code.  it is only the pressure and density
c     profiles that should match between the sph and stellar codes.  because
c     the gravitational potential depends only on density, that profile should be
c     the same in the sph and stellar evolution models as well.
      
      i=1 ! this line is important because shell number i is shared in common block
      if(neos.eq.1) then
         tem(i)=zeroin(0.d0,pres(i)*muarray(i)/(rhoarray(i)*boltz),
     $        temperaturefunction,1.d-11)
         uarray(i)=1.5d0*boltz*tem(i)/muarray(i)
     $        +arad*tem(i)**4/rhoarray(i)
      else if(neos.eq.2) then
         uupperlimit=1.5d0*boltz*temupperlimit/muarray(i)
     $        +arad*temupperlimit**4/rhoarray(i)+1d13
 1233    continue

         uarray(i)=zeroin(0.d0,uupperlimit,ufunction,1.d-11)
         if(ufunction(uarray(i)).gt.1d-10) then
            uupperlimit=2*uupperlimit
            goto 1233
         endif

         tem(i)=useeostable(uarray(i),rhoarray(i),muarray(i),1)
c         write(101,*) i,uarray(i)

      else
         print *,'this initialization script currently assumes that',
     $        'the eos is ideal gas + radiation pressure (neos=1)',
     $        'or the eos is tabulated (neos=2)'
         stop
      endif
      if(myrank.eq.0)then
         write(69,*)'central mu=',muarray(1)
         write(69,*)'central u=',uarray(1)
         write(69,*)'central temperature=',tem(1)
      endif

      integratednum=0.d0
      integratedmass2=0.d0
      integratednum2=0.d0
      utot2=0.d0
      wtota2=0.d0
      do i=2,numlines
c     rhoarray(i) is the average density in the vicinity of rarray(i)
         integratedmass2=integratedmass2+rhoarray(i)*4.d0*pi*
     $        rarray(i)**2*0.5d0*(rarray(i+1)-rarray(i-1))
         integratednum2=integratednum2+rhoarray(i)*4.d0*pi*
     $        rarray(i)**2*0.5d0*(rarray(i+1)-rarray(i-1))/
     $        muarray(i)
         wtota2=wtota2-3.d0*pres(i)*4.d0*pi*
     $        rarray(i)**2*0.5d0*(rarray(i+1)-rarray(i-1))

         temupperlimit=pres(i)*muarray(i)/(rhoarray(i)*boltz)
         if(neos.eq.1) then
            tem(i)=zeroin(0.d0,temupperlimit,
     $           temperaturefunction,1.d-11)
            uarray(i)=1.5d0*boltz*tem(i)/muarray(i)
     $           +arad*tem(i)**4/rhoarray(i)
         else if(neos.eq.2) then
            uupperlimit=1.5d0*boltz*temupperlimit/muarray(i)
     $           +arad*temupperlimit**4/rhoarray(i)+1d13
 1234       continue
            uarray(i)=zeroin(0.d0,uupperlimit,ufunction,1.d-11)
            if(ufunction(uarray(i)).gt.1d-10) then
               uupperlimit=2*uupperlimit
               goto 1234
            endif
            tem(i)=useeostable(uarray(i),rhoarray(i),muarray(i),1)
c     make sure i=1 shell is done same way up above
c            write(102,*) i,uarray(i),tem(i)

         endif
         utot2=utot2+uarray(i)*rhoarray(i)*4.d0*pi*
     $        rarray(i)**2*0.5d0*(rarray(i+1)-rarray(i-1))
      enddo

      if(myrank.eq.0)then
         write(69,*)'mass from integrating rho profile=',
     $        integratedmass2/egsol,'msun'
         write(69,*)'number from integrating=',integratednum2
         write(69,*)
     $        'utot (in sph units)=',utot2/(gravconst*munit**2/runit)
         write(69,*)
     $        'wtot (in sph units)=',wtota2/(gravconst*munit**2/runit)
      endif

      masscgs=xm(numlines)
      if(myrank.eq.0)
     $     write(69,*)'total mass=',xm(numlines),xm(numlines)/egsol

      do i=1,numlines
         if(xm(i).gt.0.85d0*masscgs) then
            redge1=rarray(i)/runit
            goto 143
         endif
      enddo
 143  continue

      radiuscgs=rarray(numlines)
      radius=radiuscgs/runit
      amass=masscgs/munit
      if(myrank.eq.0)then
         write(69,'(3(a,g11.4),a)')
     $        '  mass =',masscgs,  'grams =',masscgs/egsol,'msun =',
     $        amass,'code units'
         write(69,'(3(a,g11.4),a)')
     $        'radius =',radiuscgs,'cm    =',radiuscgs/solrad,'rsun =',
     $        radius,'code units'
      endif

      if(myrank.eq.0)then
         write(69,*)'   i   rarray(i)  rhoarray(i) pres(i) ',
     $        '     muarray(i) tem(i)     uarray(i)'
         do i=1,numlines,numlines/10
            write(69,'(i5,9g12.5)') i,
     $           rarray(i), rhoarray(i), pres(i),
     $           muarray(i),tem(i), uarray(i)
         enddo
      endif

      do i=1,numlines
c     convert radius,density and pressure to code units
         rarray(i)=rarray(i)/radiuscgs*radius
         rhoarray(i)=rhoarray(i)*(radiuscgs**3/masscgs)
     $        /(radius**3/amass)
c     pres(i)=pres(i)*(radiuscgs**2/masscgs)**2/gravconst
c     $        /(radius**2/amass)**2
         uarray(i)=uarray(i)/(gravconst*munit)*runit
         if(i.gt.1) then
            if(rarray(i).lt.rarray(i-1) .and. myrank.eq.0) then
               write(69,*) 'radius warning'
               write(69,*) i-1,rarray(i-1),rhoarray(i-1),uarray(i-1),
     $              pres(i-1)
               write(69,*) i,rarray(i),rhoarray(i),uarray(i),pres(i)
            endif
            if(rhoarray(i).gt.rhoarray(i-1) .and. myrank.eq.0)
     $           write(69,'(a,3g11.4)')
     $           'rho increases outward warning near shell',
     $           i,rhoarray(i-1),rhoarray(i)
            if(pres(i).gt.pres(i-1) .and. myrank.eq.0)
     $           write(69,*) 'pressure warning'
         endif
      enddo

      
      open(21, file='parent.sph',status='unknown')
      do i=1,numlines
         write(21,'(9g15.7)') rarray(i),pres(i)*
     $        (radiuscgs**2/masscgs)**2/gravconst/
     $        (radius**2/amass)**2,
     $        rhoarray(i),tem(i),muarray(i),uarray(i)
      enddo
      close(21)
      
      call sph_spline(rarray,rhoarray,numlines,1.d30,1.d30,rhoarray2)
      call sph_spline(rarray,uarray,numlines,1.d30,1.d30,uarray2)
      call sph_spline(rarray,muarray,numlines,1.d30,1.d30,muarray2)
      
      end
