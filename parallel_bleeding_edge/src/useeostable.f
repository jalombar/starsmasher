      subroutine readineostable
      include 'starsmasher.h'
      integer numrho,numu,iu,irho
      integer maxtablesize
      parameter(maxtablesize=1000)
      real*8 rhotable(maxtablesize),
     $     utable(maxtablesize)

      real*8 eostable(maxtablesize,maxtablesize,3)
      real*8 steprho,stepu,rhotable1,utable1,
     $     rhotablelast,utablelast
      common/eoscom/ numrho,numu,rhotable1,utable1,
     $     steprho,stepu,eostable
      integer i
      real*8 steprhotest,steputest
      real*8 xxx,yyy,zzz,starmu
      
c      if(myrank.eq.0) write(69,*)'about to read sph.eos'
      open (43,file=eosfile)
      read(43,*) xxx
      read(43,*) yyy
      read(43,*) zzz
c     see equation (2-15) of clayton:
      starmu=1.67262158d-24/(2*xxx+0.75d0*yyy+0.5d0*zzz)
      if(myrank.eq.0) then
         write(69,*)'reading ',trim(eosfile)
         write(69,*)'x,y,z=',xxx,yyy,zzz
         write(69,*) 'starmu=',starmu
      endif
      if(abs(xxx+yyy+zzz-1.d0).gt.1d-16)then
         write(69,*) 'one minus 1=',xxx+yyy+zzz-1.d0
         stop
      endif
      do i=1,2
         read(43,*)
      enddo
      read(43,*) numrho,rhotable1,rhotablelast,steprho
      read(43,*) numu,utable1,utablelast,stepu
      do i=1,2
         read(43,*)
      enddo

      do irho=1,numrho
         do iu=1,numu
            read(43,*) rhotable(irho),utable(iu),
     $           (eostable(iu,irho,i),i=1,3)
            eostable(iu,irho,1)=10d0**eostable(iu,irho,1) ! record temperature instead of log temperature
            eostable(iu,irho,2)=eostable(iu,irho,2)*1.67262158d-24 ! record mean molecular mass in cgs units
         enddo
      enddo

      close(43)

      print *,'done reading eos table'

      print *, 'numrho=',numrho,'numu=',numu
      
      if(rhotable1.ne.rhotable(1))then
         print *,'rhotable(1) mismatch',rhotable1,rhotable(1),
     $        rhotable1-rhotable(1)
         stop
      endif
      if(utable1.ne.utable(1))then
         print *,'utable(1) mismatch'
         stop
      endif
      if(rhotablelast.ne.rhotable(numrho))then
         print *,'rhotable(numrho) mismatch'
         stop
      endif
      if(utablelast.ne.utable(numu))then
         print *,'utable(numu) mismatch'
         stop
      endif

      steprhotest=(rhotable(numrho)-rhotable(1))/(numrho-1)
      steputest=(utable(numu)-utable(1))/(numu-1)

      if(abs(steprhotest-steprho).gt.1d-9) then
         print *,'steprho problem'
         stop
      endif
      if(abs(steputest-stepu).gt.1d-9) then
         print *,'stepu problem'
         stop
      endif


      end
      
      function useeostable(ucgs,rhocgs,which)
c     which=1 gives temperature
c     which=2 gives mean molecular mass mu
c     which=3 gives pressure
      include 'starsmasher.h'
      real*8 pgas,prad
      real*8 rhocgs,ucgs,beta1,temperature,gam1
      real*8 useeostable
      integer iu,irho
      integer numrho,numu
      real*8 utable1,rhotable1
      integer maxtablesize
      parameter(maxtablesize=1000)
      real*8 eostable(maxtablesize,maxtablesize,3)
      real*8 meanmu

      real*8 steprho,stepu
      common/eoscom/ numrho,numu,rhotable1,utable1,
     $     steprho,stepu,eostable

      real*8 f00,f01,f10,f11,log10rho,log10u,
     $     rholow,rhohigh,ulow,uhigh
      integer which

c     utable(iu)=utable(1)+(iu-1)*stepu
c     so, iu= (utable(iu)-utable(1))/stepu + 1

      log10rho=log10(rhocgs)
      log10u=log10(ucgs)

c      irho = min(max(1,int((log10rho-rhotable1)/steprho +1)),numrho-1)
      irho = int((log10rho-rhotable1)/steprho +1)
      iu =   int((log10u-utable1)/stepu +1)

      if(irho.ge.1 .and. irho.le.numrho-1 .and. iu.ge.1 .and.
     $     iu.le.numu-1) then
c      if(irho.ge.1 .and. irho.le.numrho-1) then

         rholow=log10rho-(rhotable1+(irho-1)*steprho)
         rhohigh=rhotable1+irho*steprho-log10rho

         if(iu.ge.1 .and. iu.le.numu-1) then
            
c     print *,'  iu=',log10u,log10rho,iu,irho
            
            ulow=log10u-(utable1+(iu-1)*stepu)
            uhigh=utable1+iu*stepu-log10u
            
c     use bi-linear interpolation among the four cartesian
c     grid points (irho,iu), (irho+1,iu), (irho,iu+1), and (irho+1,iu+1)
            f00=rholow*ulow
            f10=rhohigh*ulow
            f01=rholow*uhigh
            f11=rhohigh*uhigh
            
            useeostable=(f00*eostable(iu+1,irho+1,which)
     $           + f10*eostable(iu+1,      irho,  which)
     $           + f01*eostable(iu,        irho+1,which)
     $           + f11*eostable(iu,        irho,which))/(steprho*stepu)
            
         else if(iu.le.0) then
            
c     print *,'l iu=',log10u,log10rho,iu,irho
            
            if(which.ne.2) then
c     we are at very low specific internal energy u, where the pressure
c     or temperature should be nearly proportional to rho*u (at fixed
c     composition)
               
c     use linear interpolation between the two cartesian
c     grid points (1,irho), (1,irho+1) and in addition extrapolate
c     to smaller u
               useeostable=ucgs/10d0**utable1*
     $              (rholow*eostable(1,irho+1,which)
     $              + rhohigh*eostable(1,irho,which))/steprho
            else
               useeostable=
     $              (rholow*eostable(1,irho+1,which)
     $              + rhohigh*eostable(1,irho,which))/steprho
            endif
            
         else if(iu.ge.numu) then
            
c     print *,'h iu=',log10u,log10rho,iu,irho
            
            if(which.eq.3) then
c     we are at very high specific internal energy u, where the pressure
c     nearly proportional to u (at fixed rho and composition)
               
c     use linear interpolation between the two cartesian
c     grid points (irho,numu), (irho+1,numu) and in addition extrapolate
c     to larger u
               useeostable=ucgs/10d0**(utable1+(numu-1)*stepu)*
     $              (rholow*eostable(numu,irho+1,which)
     $              + rhohigh*eostable(numu,irho,which))/steprho
            else if(which.eq.1) then
c     we are at very high specific internal energy u, where the temperature
c     is nearly proportional to u^(1/4) (at fixed rho and composition)
               
c     use linear interpolation between the two cartesian
c     grid points (irho,numu), (irho+1,numu) and in addition extrapolate
c     to larger u
               useeostable=
     $              (ucgs/10d0**(utable1+(numu-1)*stepu))**0.25d0*
     $              (rholow*eostable(numu,irho+1,which)
     $              + rhohigh*eostable(numu,irho,which))/steprho
            else
               useeostable=
     $              (rholow*eostable(numu,irho+1,which)
     $              + rhohigh*eostable(numu,irho,which))/steprho
            endif
            
         endif
      else
c     at extreme densities we will use ideal gas + radiation pressure

         if(irho.lt.1) irho=1
         if(irho.gt.numrho) irho=numrho

         if(iu.lt.1) then
            meanmu=eostable(1,irho,2)
         else if (iu.ge.numu) then
            meanmu=eostable(numu,irho,2)
         else
            ulow=log10u-(utable1+(iu-1)*stepu)
            uhigh=utable1+iu*stepu-log10u
            meanmu=(ulow*eostable(iu+1,irho,2)
     $           + uhigh*eostable(iu,  irho,2))/stepu
         endif

c     print *,'input=',qconst*rhocgs/meanmu,
c     $        -ucgs*rhocgs/arad

         if(which.eq.2) then
            useeostable=meanmu
            return
         endif
         call gettemperature(qconst*rhocgs/meanmu,
     $        -ucgs*rhocgs/arad,temperature)

c         print *,'irho, iu, meanmu, t=',irho,iu,meanmu,temperature

         if(which.eq.1) then
            useeostable=temperature
            return
         endif
         pgas=rhocgs*boltz*temperature/meanmu
         prad=arad*temperature**4/3.d0
         beta1=pgas/(pgas+prad)
         gam1=(32.d0-24.d0*beta1-3.d0*beta1**2) /
     $        (24.d0-21.d0*beta1)
         
         if(gam1.lt.0.999*4.d0/3.d0 .or. gam1.gt.1.001*5.d0/3.d0) then
            write(69,*)'warning gam1=',gam1
            write(69,*) beta1,pgas,prad,temperature,
     $           -ucgs*rhocgs/arad
            stop 'gam1 value does not make sense'
         endif
         useeostable=pgas+prad
         
      endif
         
      end
