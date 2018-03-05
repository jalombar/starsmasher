      subroutine polyes
c     creates a star from the data file yrec output
      include 'starsmasher.h'
      include 'mpif.h'
      integer i
      integer ntry,idumb,ip,ix,iy,iz,maxtry
      parameter(maxtry=300000000)
      real*8 anumden,rhotry,rhoex,rtry,rhomax,hc,xcm,ycm,zcm,amtot,
     $     ammin,ammax,xtry,ytry,ztry,ri,rhoi,ran1
      integer irtry
      real*8 amass,radius
      integer ixmax,iymax,izmax,corepts
      double precision cellvolume,a1
      real*8 integral,deltar
      real*4 npoly
      integer nrgrid
      parameter(nrgrid=10000)                                         
      real rhopol(nrgrid)
      real ak
      real*8 rarrayi,rarrayim1
      double precision utottest
      real*8 epot
      real*8 hpguess,xacc,dxmax
      real*8 drhodhi
      real*8 utot2,wtota2
      integer mygravlength, ierr
      integer comm_worker, irank
      common/gravworkers/comm_worker
      integer status(mpi_status_size)

      integer maxtablesize
      parameter(maxtablesize=1000)
      integer numrho,numu,numx,iu,irho,iup
      real*8 eostable(maxtablesize,maxtablesize,maxnumx,3)
      real*8 zzz,steprho,stepu,stepx,rhotable1,utable1,xtable1,
     $     rhotablelast,utablelast,xtablelast
      common/eoscom/ zzz,rhotable1,utable1,xtable1,
     $     steprho,stepu,stepx,eostable,numrho,numu,numx

      real*8 rhocgs,log10rho,ucgsguess,pressurecgs,log10u
      real*8 rholow,rhohigh
      real*8 plow,phigh,stepp,steppp,plowp,phighp,f00,f10,f01,f11

      npoly=1.d0/(gam-1.d0)
      amass=starmass
      radius=starradius
      
      if(myrank.eq.0) then
         write(69,*)'making polytrope with n=',npoly
         write(69,*)'                    gam=',gam
         write(69,*)'                      m=',amass
         write(69,*)'                      r=',radius
      endif

c     get profiles:
      call poly(npoly,sngl(amass),sngl(radius),nrgrid,rhopol,ak)                 
      
      idumb=-2391
      rhomax=dble(rhopol(1))
      if(myrank.eq.0)then
         write(69,*)'central density=',rhomax
         write(69,*)'surface density=',rhopol(nrgrid)
      endif
      ip=0               
      corepts=0
      hc=(nnopt/(8.d0*n))**(1.d0/3.d0)*radius
      redge=radius-2.d0*hc
      if(myrank.eq.0)write(69,*)'keeping particles up to a distance',radius-redge,
     $     'less than the full radius',radius
      if(equalmass.lt.1.d0)then
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     make an hcp lattice
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         
c     the fraction of particles at a region of density rhoex that
c     will be kept is (rhoex/rhomax)**equalmass.  so if the number
c     density of lattice points that will be tried is n, then we expect
c     the number of particles to be n=n*integrate[4*pi*r**2*
c     (rhoex(r)/rhomax)**equalmass,{r,0,redge}].  we can solve this for
c     n, and then use that the cell volume is 2/n (as there are two
c     particles per cell)
         
         if(equalmass.gt.0) then
            i=2
            integral=0
            rarrayi=radius*(i-0.5d0)/dble(nrgrid-1)
            do while(rarrayi.lt.redge)
               rarrayim1=radius*(i-1.5d0)/dble(nrgrid-1)
               integral=integral+pi*(rarrayi+rarrayim1)**2*
     $              (rarrayi-rarrayim1)*
     $              (0.5d0*dble(rhopol(i)+rhopol(i-1))/rhomax)
     $              **equalmass
               i=i+1
               rarrayi=radius*(i-0.5d0)/dble(nrgrid-1)
            enddo
            if(myrank.eq.0)write(69,*)'integral=',integral
            integral=0.d0
            deltar=redge/1000.d0
            do irtry=nint(1000*deltar),nint(1000*redge),nint(1000*deltar+0.5)
               rtry=irtry/1000.d0
               rhoex=dble(rhopol(1+int(rtry/radius*dble(nrgrid-1))))
               integral=integral+4.d0*pi*rtry**2*deltar*
     $              (rhoex/rhomax)**equalmass
            enddo
         else
            if(myrank.eq.0)write(69,*)'maybe integral=',4.d0/3.d0*pi*radius**3
            integral=4.d0/3.d0*pi*redge**3
         endif
         if(myrank.eq.0)write(69,*)'using integral=',integral
         cellvolume=2.d0*integral/n
         a1=(cellvolume/2.d0**0.5d0)**(1.d0/3.d0)
         if(myrank.eq.0)write(69,*)'a1=',a1
         
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
     $                 -(mod(abs(iz),2)-0.5d0)*1.d0/3.d0**0.5d0*a1
                  ztry=(iz-0.5d0)*0.5d0*(8.d0/3.d0)**0.5d0*a1
                  rtry=sqrt(xtry**2.d0+ytry**2.d0+ztry**2.d0)
                  if(rtry.lt.redge) then
                     rhoex=dble(
     $                    rhopol(1+int(rtry/radius*dble(nrgrid-1))))
                     rhotry=rhomax**equalmass*ran1(idumb)      
                     if (rhotry.le.rhoex**equalmass) then
c     (particle is accepted)    
                        ip=ip+1                 
                        if(rtry.le.a1) then
                           if(myrank.eq.0)write(69,'(4i5,4e11.4)')ip,ix,iy,iz,
     $                          xtry,ytry,ztry,rtry
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
         n=ip
         if(myrank.eq.0)write (69,*) 'parent: n=',n
         if (n.gt.nmax) stop 'parent: n>nmax ???'
c     assign particle masses (to represent density):
c     cden=(4.d0*pi)/(3.d0*n)*redge**3.d0
         ammin=1.d30
         ammax=0.d0
         amtot=0.d0
         xcm=0.d0
         ycm=0.d0
         zcm=0.d0
         do i=1,n
            ri=sqrt(x(i)**2+y(i)**2+z(i)**2)
            rhoi=dble(rhopol(1+int(ri/radius*dble(nrgrid-1))))
            
            if(rhoi.le.0.d0 .and. myrank.eq.0) then
               write(69,*)'warning: rho(',i,')<=0 at r=',ri,'???'
            endif
            am(i)=amass/n*(integral*rhoi/amass)**(1.d0-equalmass)
            xcm=xcm+am(i)*x(i)
            ycm=ycm+am(i)*y(i)
            zcm=zcm+am(i)*z(i)
            amtot=amtot+am(i)
            if(neos.lt.2) then
               if(nintvar.eq.1) then
c     p=(ak)*rho^(1+1/n)=a*rho^gam, so a=(ak)*rho^(1+1/n-gam):
                  u(i)=dble(ak)*rhoi**(1.d0+1.d0/npoly-gam)
               else
c     p=(ak)*rho^(1+1/n)=(gam-1)*rho*u, so u=(ak)*rhoi**(1/n)/(gam-1)
                  u(i)=dble(ak)*rhoi**(1.d0/npoly)/(gam-1.d0)
               endif
            else
               rhocgs=rhoi*munit/runit**3.d0
               log10rho=log10(rhocgs)
               irho = int((log10rho-rhotable1)/steprho +1)
               if(irho.le.0 .or. irho.ge.numrho)then
                  write(69,*)'irho problem',i,irho,log10rho
                  stop
               endif
               irho=min(max(irho,1),numrho)
               ucgsguess=dble(ak)*rhoi**(1.d0/npoly)/(gam-1.d0)*
     $              gravconst*munit/runit
               log10u=log10(ucgsguess)
               iu =   int((log10u-utable1)/stepu +1)
               iu=min(max(iu,1),numu)

               pressurecgs=dble(ak)*rhoi**(1.d0+1.d0/npoly)*
     $              gravconst*munit**2/runit**4
c               write(35,*)'numu=',numu
c               numberu=numu
c               write(35,*)eostable(1,irho,3),eostable(2,irho,3),'...',
c     $              eostable(229,irho,3),eostable(230,irho,3),'...',
c     $              eostable(numu,irho,3)
c               write(35,*)i,iu,log10u,irho,log10rho,pressurecgs

               call hunt(eostable(1,irho,1,3),numu,pressurecgs,iu)
               call hunt(eostable(1,irho+1,1,3),numu,pressurecgs,iup)
               if(iu.le.0 .or. iu.ge.numu .or.
     $              iup.le.0 .or. iup.ge.numu)then
                  write(69,*)'iu problem',i,iu,iup,log10u,irho,log10rho
                  stop
               endif

               rholow=log10rho-(rhotable1+(irho-1)*steprho)
               rhohigh=rhotable1+irho*steprho-log10rho ! equals steprho-rholow
               plow=pressurecgs-eostable(iu,irho,1,3)
               phigh=eostable(iu+1,irho,1,3)-pressurecgs
               stepp=eostable(iu+1,irho,1,3)-eostable(iu,irho,1,3)
               plowp=pressurecgs-eostable(iup,irho+1,1,3)
               phighp=eostable(iup+1,irho+1,1,3)-pressurecgs
               steppp=eostable(iup+1,irho+1,1,3)-eostable(iup,irho+1,1,3)
c     use bi-linear interpolation among the four cartesian
c     grid points (irho,iu), (irho+1,iu), (irho,iu+1), and (irho+1,iu+1)
               f00=rholow*plowp/(steprho*steppp)
               f10=rhohigh*plow/(steprho*stepp)
               f01=rholow*phighp/(steprho*steppp)
               f11=rhohigh*phigh/(steprho*stepp)

c     if everything is correct, the following expression should equal pressurecgs
c                      f00*eostable(iup+1,irho+1,3)
c     $              + f10*eostable(iu+1, irho,  3)
c     $              + f01*eostable(iup,  irho+1,3)
c     $              + f11*eostable(iu,   irho,  3)

               if(abs(f00+f10+f01+f11-1d0).gt.0.0000000000001d0)then
                  write(69,*)'interpolation problem',i,f00,f10,f01,f11
                  stop
               endif

               u(i)=10d0**(utable1+
     $              stepu*(f00*iup+f10*iu+f01*(iup-1)+f11*(iu-1)))/
     $              (gravconst*munit/runit)

c               log10u=utable1+(iu-1)*stepu
c               log10rho=rhotable1+(irho-1)*steprho
c               write(35,*)i,iu,log10u,irho,log10rho,pressurecgs
c               log10u=utable1+(iup-1)*stepu
c               log10rho=rhotable1+(irho+1-1)*steprho
c               write(35,*)i,iup,log10u,irho+1,log10rho,pressurecgs
c               write(35,*)dble(ak)*rhoi**(1.d0/npoly)/(gam-1.d0),u(i)

c     write(35,*)    eostable(iu,irho,3), eostable(iu+1,irho,3)

            endif
            meanmolecular(i)=0.10329d-23
         enddo
         xcm=xcm/amtot
         ycm=ycm/amtot
         zcm=zcm/amtot
         if(myrank.eq.0) then
            write (69,'(a,g11.4,a,g11.4,a)')
     $           'parent: total mass of star was m=',amtot
     $           ,
     $           ', which should equal', amass,'... renormalizing...'
            write(69,*) 'center of mass at',xcm,ycm,zcm
         endif
         do i=1,n
            am(i)=am(i)/amtot*amass
            if(nintvar.eq.2) then
               u(i)=u(i)/amass*amtot
            endif
            if(am(i).gt.0.01d0*amass .and. myrank.eq.0) then
               write(69,*)'warning: particle',i,'has mass',am(i)
               write(69,*)'x,y,z=',x(i),y(i),z(i)
               write(69,*)'r,rho=',ri,rhoi
            endif
            ammin=min(ammin,am(i))
            ammax=max(ammax,am(i))
            rtry=sqrt(x(i)**2+y(i)**2+z(i)**2)
            rhoex=dble(rhopol(1+int(rtry/radius*dble(nrgrid-1))))
            anumden=rhoex/am(i)
            hp(i)=(3.d0/32.d0/3.1415926535897932384626d0*
     $           1.9d0*nnopt/anumden)**(1.d0/3.d0)
         enddo
         if(myrank.eq.0)then
            write (69,*) 'parent: min particle mass=',ammin
            write (69,*) 'parent: max particle mass=',ammax
         endif
         
c     because the particles were distributed uniformly in the sphere, the
c     the number density is 3*n/(4*pi*radius**3) and we need to choose
c     the smoothing length hc such that 
c     4*pi*(2*hc)**3/3 * (number density)=nnopt.
c     this gives 8*hc**3 * n/radius**3= nnopt, or:
c     hc=(nnopt/(8.d0*n))**(1.d0/3.d0)*redge
c     (should give a number of nearest neighbors close to nnopt)
         
      else
         if(myrank.eq.0)then
            if(equalmass.eq.1.d0) then
               write(69,*)'using equal mass particles'
            else
               write(69,*)'compromising with particle masses ...',
     $              equalmass
            endif
         endif
c     lay down particles to represent density (rejection method): 
         amtot=0.d0
         xcm=0.d0
         ycm=0.d0
         zcm=0.d0
         ammin=1.d30
         ammax=0.d0
         do ntry=1,maxtry   
            xtry=redge*(-1.d0+2.d0*ran1(idumb))  
            ytry=redge*(-1.d0+2.d0*ran1(idumb))
            ztry=redge*(-1.d0+2.d0*ran1(idumb))
            rtry=sqrt(xtry**2+ytry**2+ztry**2) 
            if (rtry.lt.redge) then          
               rhoex=dble(rhopol(1+int(rtry/radius*dble(nrgrid-1))))
               rhotry=rhomax**equalmass*ran1(idumb)      
               if(rhoex.gt.0.d0.and.rhotry.lt.rhoex**equalmass)then
c     (particle is accepted)    
                  ip=ip+1                 
                  x(ip)=xtry        
                  y(ip)=ytry            
                  z(ip)=ztry
c     if equalmass=0 ... am=volume/n*rho
c     if equalmass=1 ... am=amass/n
c     in general     ... am=amass*(rho/averagerho)^(1-equalmass)/n
                  am(ip)=amass/n*
     $                 (4.d0/3.d0*pi*redge**3*rhoex/amass)
     $                 **(1.d0-equalmass)
                  
                  if(nintvar.eq.1) then
                     u(ip)=dble(ak)*rhoex**(1.d0+1.d0/npoly-gam)
                  else
                     u(ip)=dble(ak)*rhoex**(1.d0/npoly)/(gam-1.d0)
                  endif
                  meanmolecular(ip)=0.10329d-23
                  ammin=min(ammin,am(ip))
                  ammax=max(ammax,am(ip))
                  xcm=xcm+am(ip)*x(ip)
                  ycm=ycm+am(ip)*y(ip)
                  zcm=zcm+am(ip)*z(ip)
                  amtot=amtot+am(ip)
               endif
            endif                     
            if (ip.eq.n) goto 1       
         end do
         if(myrank.eq.0)write(69,*) 'only got to ip=',ip
         stop 'init: not enough particles ???'
 1       continue                
         xcm=xcm/amtot
         ycm=ycm/amtot
         zcm=zcm/amtot
         if(myrank.eq.0)write (69,'(a,g11.4,a,g11.4,a)')
     $        'parent: total mass was m=',amtot,
     $        ', which should equal', amass,'... renormalizing...'
         do i=1,n
            am(i)=am(i)/amtot*amass
            rtry=sqrt(x(i)**2+y(i)**2+z(i)**2)
            rhoex=dble(rhopol(1+int(rtry/radius*dble(nrgrid-1))))
            anumden=rhoex/am(i)
            
c     if particles were distributed uniformly in the sphere 
c     (they're not here) then the number density would be 
c     3*n/(4*pi*radius**3) 
c     and we'd need to choose the smoothing length hc such that 
c     4*pi*(2*hc)**3/3 * (number density)=nnopt.
c     this gives 
c     8*hc**3 * n/radius**3= nnopt, or
c     hc=(nnopt/(8.d0*n))**(1.d0/3.d0)*radius
c     which is what we use when equalmass=0.
c     if there is not uniform distribution (as done here) then we 
c     should instead let the number density be 
c     anumden=rhoex/ami(i)
c     then 
c     4*pi*(2*hp(i))**3/3 * anumden=nnopt gives
            
            hp(i)=(3.d0/32.d0/3.1415926535897932384626d0*
     $           1.9d0*nnopt/anumden)**(1.d0/3.d0)
            
c     the following makes the hp smaller to attempt to avoid having 
c     too many neighbors... hp's are increased gradually with a 
c     do while loop below....
c     this line may be helpful only if using the grape board 
c     (to avoid potentially crashing it if a particle has too many 
c     neighbors), even then it may not really be useful, as maybe 
c     the code recovers anyway:
c     hp(i)=hp(i)*(rhoex/rhomax)**(1.d0/3.d0)
         enddo
         ammin=ammin/amtot*amass
         ammax=ammax/amtot*amass
      endif      
      ntot=n+corepts
      
      call stride_setup
      if(ntot.gt.nmax) then
         if(myrank.eq.0)write(69,*) 'n is too large'
         stop
      endif

      if(myrank.eq.0)write(69,*) amass,'solar masses has cc=',
     $     int(10000*amass)
      do i=1,n
         cc(i)=int(10000*amass)
      enddo

      if(corepts.eq.0 .and. treloff.le.0.d0) then
         if(myrank.eq.0) write(69,*)'will try to get correct u and w...'

         utottest=0.d0
         do i=1,n
            utottest=utottest+u(i)*am(i)
         enddo

         utot2=1.5d0/(5.d0-npoly)*amass**2/radius
         if(myrank.eq.0) write (69,'(a,g11.4,a,g11.4,a)')
     $        'parent: total internal energy of star was',utottest
     $        ,
     $        ', which should equal', utot2, 
     $        '... renormalizing...'
         do i=1,n
            u(i)=u(i)*utot2/utottest
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

            wtota2=-2.d0*utot2
            if(myrank.eq.0)write (69,'(a,g11.4,a,g11.4,a)')
     $           'parent: gravitational potential energy of star was',epot
     $           ,
     $           ', which should equal', wtota2,
     $           '... rescaling positions...'
            do i=1,n
               x(i)=x(i)*epot/wtota2
               y(i)=y(i)*epot/wtota2
               z(i)=z(i)*epot/wtota2
               hp(i)=hp(i)*epot/wtota2
            enddo
         endif
      endif

      if(myrank.eq.0)then
         write(69,'(5g11.4)') 'i', 'nn(i)','hp(i)'
         do i=1,n
            if(mod(i,1000).eq.1)
     $           write(69,'(5g11.4)') i, nn(i),hp(i)
         enddo
         write(69,*) 'i   density'
         do i=1,nrgrid,nrgrid/10
            write(*,*) i,rhopol(i)
         enddo
         
         write(69,*)'rhopol(1) and rhopol(2)=',
     $        rhopol(1),rhopol(2)
         write(69,*)'rhopol(nrgrid-1),rhopol(nrgrid) =',
     $        rhopol(nrgrid-1),
     $        rhopol(nrgrid)
         write(69,*)'near origin:'
         write(69,*)' r    rho'
         do irtry=0,nint(10*radius),nint(radius+0.5)
            rtry=irtry/100.d0
            rhoi=dble(rhopol(1+int(rtry/radius*dble(nrgrid-1))))
            write(69,*) rtry,rhoi
         enddo
         write(69,*)'near surface:'
         write(69,*)' r    rho'
         do irtry=nint(90*radius),nint(100*radius),nint(radius+0.5)
            rtry=min(irtry/100.d0,radius)
            rhoi=dble(rhopol(1+int(rtry/radius*dble(nrgrid-1))))
            write(69,*) rtry,rhoi
         enddo
      endif

c     assign velocities (all zero):
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
      if(myrank.eq.0)then
         write(69,*)'polyes: tf=',tf,dtout
         write(69,*) 'polyes: exiting polyes'
      endif
      return

      end
**************************************************************
c     this file is part of the starcrash code
c     version 1.0
c
c     copyright 2003, by joshua faber
c     this code is publicly available under the terms of the
c     gnu general public license.  see our documentation, or
c     http://www.gnu.org/licenses/gpl.html for more details.
***************************************************************

      subroutine poly(an,am,r,nr,rho,ak)
*****************************************************
c     release 1.0
c     calculates specific entropy a for polytrope with index n,
c     mass am, radius r, and gives density rho at nr points in radius
c     this routine does not "include spha.h", thus rho is unambiguous
c     called by setup1es,setup1em,setup2cs,setup2cm
************************************************** 
c      implicit double precision(a-h,o-z)
      implicit none
      double precision xmax,pi,andbl
      integer maxstp,nrm
      parameter (xmax=7.d0,maxstp=15,pi=3.14159265359d0,nrm=20000)
      integer ir,nr,nstp
      real*4 an,am,r,rho(nr),ak
      double precision xtab(nrm),ytab(nrm),yptab(nrm)
      common/polidx/andbl
      double precision xs,dr,y1,yp1,rhoc,yps,x1!,ys
      
      andbl=dble(an)

      xtab(1)=0.d0
      ytab(1)=1.d0
      yptab(1)=0.d0
      xs=xmax
      do 20 nstp=1,maxstp
c advance solution to r=dr using series expansion near origin:
       dr=xs/dble(nr-1)
       y1=1.d0-dr**2/6.d0+andbl*dr**4/120.d0
       yp1=-dr/3.d0+andbl*dr**3/30.d0
       x1=dr
       call rktab(y1,yp1,x1,xs,nr-1,xtab(2),ytab(2),yptab(2))
       do 10 ir=2,nr
10      if (ytab(ir).lt.0.d0) goto 11
11     xs=xtab(ir)
c       ys=ytab(ir)
       yps=yptab(ir)
       if (ir.eq.nr) goto 21
20    continue
      stop 'poly: no convergence ???'
21    continue
      ytab(nr)=0.d0

      rhoc=xs*dble(am)/(4.d0*pi*dabs(yps)*dble(r)**3)
      ak=sngl(4.d0*pi*dble(r)**2*rhoc**(1.d0-1.d0/andbl)/((andbl+1.d0)
     +                                                        *xs**2))
      do 30 ir=1,nr
30     rho(ir)=sngl(rhoc*ytab(ir)**andbl)

      return
      end
************************************************************************
      subroutine derivs(x,v,dv)
      implicit none
      double precision x,v(2),dv(2),an
      common/polidx/an
      dv(1)=v(2)
      dv(2)=-2.d0*v(2)/x-dabs(v(1))**an
      return
      end
************************************************************************
      subroutine rktab(y,yp,x1,x2,ntab,xtab,ytab,yptab)

c      implicit double precision(a-h,o-z)
      implicit none
      integer ntab
      double precision xtab(ntab),ytab(ntab),yptab(ntab)
      double precision v(2),dv(2)
      double precision x1, y, yp, x, h, x2
      integer k

      xtab(1)=x1
      ytab(1)=y
      yptab(1)=yp
      x=x1
      h=(x2-x1)/dble(ntab-1)
      v(1)=y
      v(2)=yp
      do 13 k=1,ntab-1
        call derivs(x,v,dv)
        call rk4(v,dv,2,x,h,v)
        if (x+h.eq.x) stop 'rktab: stepsize not significant ???'
        x=x+h
        xtab(k+1)=x
        ytab(k+1)=v(1)
        yptab(k+1)=v(2)
13    continue
      return
      end
************************************************************************
c this is the standard 4th order runge-kutta algorithm:
      subroutine rk4(y,dydx,n,x,h,yout)
c      implicit double precision(a-h,o-z)
      implicit none
      integer nmax
      integer i,n
      parameter (nmax=10)
      double precision y(n),dydx(n),yout(n),yt(nmax),dyt(nmax),dym(nmax)
      double precision h, h6, hh, x, xh

      if (n.gt.nmax) stop 'rk4: n>nmax ???'

      hh=h*0.5d0
      h6=h/6.d0
      xh=x+hh
      do 11 i=1,n
        yt(i)=y(i)+hh*dydx(i)
11    continue
      call derivs(xh,yt,dyt)
      do 12 i=1,n
        yt(i)=y(i)+hh*dyt(i)
12    continue
      call derivs(xh,yt,dym)
      do 13 i=1,n
        yt(i)=y(i)+h*dym(i)
        dym(i)=dyt(i)+dym(i)
13    continue
      call derivs(x+h,yt,dyt)
      do 14 i=1,n
        yout(i)=y(i)+h6*(dydx(i)+dyt(i)+2.d0*dym(i))
14    continue
      return
      end
