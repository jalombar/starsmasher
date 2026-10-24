      subroutine readineostable
      include 'starsmasher.h'
      integer numrho,numu,numx,iu,irho,ix
      integer maxtablesize
      parameter(maxtablesize=1000)
      real*8 rhotable(maxtablesize),
     $     utable(maxtablesize)
      real*8 eostable(maxtablesize,maxtablesize,maxnumx,3)
      real*8 steprho,stepu,stepx,rhotable1,utable1,xtable1,
     $     rhotablelast,utablelast,xtablelast
      real*8 xxx,yyy,zzz,eosmu
      common/eoscom/ zzz,rhotable1,utable1,xtable1,
     $     steprho,stepu,stepx,eostable,numrho,numu,numx
      integer i
      real*8 steprhotest,steputest !, ugas, pgas

      if(myrank.eq.0) write(69,*)'about to read EOS file ',trim(eosfile)
      open (43,file=eosfile)
      read(43,*, err=98) numx,xtable1,xtablelast,stepx
 8    if(myrank.eq.0)
     $     write(69,*)'There are',numx,
     $     'hydrogen abundances, ranging from',xtable1,'to',
     $     xtablelast,'in steps of',stepx
      if(numx.gt.maxnumx) then
         write(69,*) 'Increase maxnumx to accomodate EOS file data'
         stop
      endif
      goto 102
 98   if(myrank.eq.0)
     $     write(69,*)'EOS file is for a single composition'
      numx=1
      stepx=0
      rewind(43)
 102  continue

      do ix=1,numx

         read(43,*) xxx
         read(43,*) yyy
         read(43,*) zzz
c     see equation (2-15) of clayton:
         eosmu=1.67262158d-24/(2*xxx+0.75d0*yyy+0.5d0*zzz)
         if(myrank.eq.0) then
c     write(69,*)'reading ',trim(eosfile)
            write(69,*)'x,y,z=',xxx,yyy,zzz
            write(69,*) 'EOS mu=',eosmu
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
     $              (eostable(iu,irho,ix,i),i=1,3)
               eostable(iu,irho,ix,1)=10d0**eostable(iu,irho,ix,1) ! record temperature instead of log temperature
               eostable(iu,irho,ix,2)=eostable(iu,irho,ix,2)*1.67262158d-24 ! record mean molecular mass in cgs units

c               if(myrank.eq.0) then
c                  ugas = 10d0**(utable1+(iu-1)*stepu)
c     $                 -arad*eostable(iu,irho,ix,1)**4
c     $                 /10d0**(rhotable1+(irho-1)*steprho)
c                  if(ugas.le.0) write (69,*)'WARNING: ugas=',ugas,irho,
c     $                 iu, 10d0**(utable1+(iu-1)*stepu), arad*eostable(iu,irho,ix,1)**4
c     $                 /10d0**(rhotable1+(irho-1)*steprho)
c                  pgas = eostable(iu,irho,ix,3)
c     $                 -arad*eostable(iu,irho,ix,1)**4/3
c                  if(pgas.le.0) write (69,*)'WARNING: pgas=',pgas,irho,
c     $                 iu, eostable(iu,irho,ix,3), arad*eostable(iu,irho,ix,1)**4/3
c               endif

            enddo
         enddo
         
         print *,'done reading eos table'
         
         print *, 'numrho=',numrho,'numu=',numu
         
         if(rhotable1.ne.rhotable(1))then
            print *,'rhotable(1) mismatch',rhotable1,rhotable(1),
     $           rhotable1-rhotable(1)
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

      enddo

      close(43)

      if(myrank.eq.0) write(69,*) 'readineostable: exiting'
         
      end
      
      function useeostable(ucgs,rhocgs,particlemu,which)
c     which=1 gives temperature
c     which=2 gives mean molecular mass mu
c     which=3 gives pressure
      include 'starsmasher.h'
      real*8 pgas,prad
      real*8 rhocgs,ucgs,beta1,temperature,gam1
      real*8 useeostable
      integer iu,irho,ix
      integer numrho,numu,numx
      real*8 utable1,rhotable1,xtable1
      integer maxtablesize
      parameter(maxtablesize=1000)
      real*8 eostable(maxtablesize,maxtablesize,maxnumx,3)

      real*8 steprho,stepu,stepx
      real*8 xxx,zzz,particlemu,meanmu
      common/eoscom/ zzz,rhotable1,utable1,xtable1,
     $     steprho,stepu,stepx,eostable,numrho,numu,numx

      real*8 f00,f01,f10,f11,log10rho,log10u,
     $     rholow,rhohigh,ulow,uhigh,xlow,xhigh
      real*8 useeostable0,useeostable1,meanmu0,meanmu1
      integer which, iuuse
      real*8 meanT4_0, meanT4_1, T4boundary, meanP0, meanP1, Pboundary,
     $     uboundary, rhoboundary, ugasboundary, newfirstargument,
     $     pgasboundary, Tboundary

c     utable(iu)=utable(1)+(iu-1)*stepu
c     so, iu= (utable(iu)-utable(1))/stepu + 1

      log10rho=log10(rhocgs)
      log10u=log10(ucgs)

c      irho = min(max(1,int((log10rho-rhotable1)/steprho +1)),numrho-1)
      irho = int((log10rho-rhotable1)/steprho +1)
      iu =   int((log10u-utable1)/stepu +1)

c     use particlemu and zzz to get hydrogen abundance xxx:
      xxx=(1.67262158d-24/particlemu+0.25d0*zzz-0.75d0)/1.25d0
      ix = min(int((xxx-xtable1)/stepx +1),maxnumx-1)

c     Note: xlow+xhigh = stepx
c           if xlow=0 then want to ignore the larger xxx value
c           if xhigh=0 then want to ignore the smaller xxx value
      xlow=xxx-(xtable1+(ix-1)*stepx)
      xhigh=xtable1+ix*stepx-xxx

c      if(myrank.eq.0 .and. ix.ne.1) then
c         write(69,*) 'Particle with mu=',particlemu,'has X=',xxx,
c     $        'and ix=',ix
c         stop
c      endif

      if((ix.lt.1 .or. ix.ge.numx) .and. numx.gt.1) then
         if(myrank.eq.0) then
            write(69,*)'xtable1=',xtable1
            write(69,*)'stepx=',stepx
            write(69,*)'numx=',numx
            write(69,*)'particle Z=',zzz
            write(69,*)'particle mu=',particlemu 
            write(69,*)'particle X=',xxx
            write(69,*)'ix=',ix
            write(69,*)'Expand X range covered by EOS file'
         endif

         if(ix.lt.1) then
            ix=1
         else
            stop
         endif
      endif

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
            
            if(numx.gt.1) then
               useeostable0=
     $                 f00*eostable(iu+1,irho+1,ix,which)
     $              +  f10*eostable(iu+1,irho,  ix,which)
     $              +  f01*eostable(iu,  irho+1,ix,which)
     $              +  f11*eostable(iu,  irho,  ix,which)
               useeostable1=
     $                 f00*eostable(iu+1,irho+1,ix+1,which)
     $              +  f10*eostable(iu+1,irho,  ix+1,which)
     $              +  f01*eostable(iu,  irho+1,ix+1,which)
     $              +  f11*eostable(iu,  irho,  ix+1,which)
               useeostable=(xlow*useeostable1+xhigh*useeostable0)
     $              /(steprho*stepu*stepx)
            else
               useeostable=
     $              (f00*eostable(iu+1,irho+1,1,which)
     $              +  f10*eostable(iu+1,irho,  1,which)
     $              +  f01*eostable(iu,  irho+1,1,which)
     $              +  f11*eostable(iu,  irho,  1,which))/(steprho*stepu)
            endif

         else if(iu.le.0) then
            
c     print *,'l iu=',log10u,log10rho,iu,irho
            
c     First we get the value at the edge of the table with smallest u 
            if(numx.gt.1) then
               useeostable0=rholow*eostable(1,irho+1,ix,  which)
     $              +      rhohigh*eostable(1,irho,  ix,  which)
               useeostable1=rholow*eostable(1,irho+1,ix+1,which)
     $              +      rhohigh*eostable(1,irho,  ix+1,which)
               useeostable=
     $              (xlow*useeostable1+xhigh*useeostable0)
     $              /(steprho*stepx)
            else
               useeostable=(rholow*eostable(1,irho+1,1,which)
     $              +     rhohigh*eostable(1,irho,  1,which))/steprho
            endif
            if(which.ne.2) then
c     Reminder: which=2 returns mean molecular weight from EOS table
               
c     we are at very low specific internal energy u, where the pressure
c     or temperature should be nearly proportional to rho*u (at fixed
c     composition).  We use this to extrapolate to smaller u:
               useeostable=ucgs/10d0**utable1*useeostable
            endif

         else if(iu.ge.numu) then
            
c     print *,'h iu=',log10u,log10rho,iu,irho
            
            if(which.eq.3) then
c     we are at very high specific internal energy u, where the pressure
c     nearly proportional to u (at fixed rho and composition)
               
c     use linear interpolation between the two cartesian
c     grid points (irho,numu), (irho+1,numu) and in addition extrapolate
c     to larger u
               if(numx.gt.1) then
                  useeostable0=rholow*eostable(numu,irho+1,ix,  which)
     $                 +      rhohigh*eostable(numu,irho,  ix,  which)
                  useeostable1=rholow*eostable(numu,irho+1,ix+1,which)
     $                 +      rhohigh*eostable(numu,irho,  ix+1,which)
                  useeostable=ucgs/10d0**(utable1+(numu-1)*stepu)*
     $                 (xlow*useeostable1+xhigh*useeostable0)
     $                 /(steprho*stepx)
               else
                  useeostable=(rholow*eostable(numu,irho+1,1,which)
     $                 +      rhohigh*eostable(numu,irho, 1,which))/steprho
               endif
            else if(which.eq.1) then
c     we are at very high specific internal energy u, where the temperature
c     is nearly proportional to u^(1/4) (at fixed rho and composition)
               
c     use linear interpolation between the two cartesian
c     grid points (irho,numu), (irho+1,numu) and in addition extrapolate
c     to larger u
               if(numx.gt.1) then
                  useeostable0=rholow*eostable(numu,irho+1,ix,  which)
     $                 +      rhohigh*eostable(numu,irho,  ix,  which)
                  useeostable1=rholow*eostable(numu,irho+1,ix+1,which)
     $                 +      rhohigh*eostable(numu,irho,  ix+1,which)
                  useeostable=
     $                 (ucgs/10d0**(utable1+(numu-1)*stepu))**0.25d0*
     $                 (xlow*useeostable1+xhigh*useeostable0)
     $                 /(steprho*stepx)
               else
                  useeostable=(rholow*eostable(numu,irho+1,1,which)
     $                 +      rhohigh*eostable(numu,irho, 1,which))/steprho
               endif
            else
               if(numx.gt.1) then
                  useeostable0=rholow*eostable(numu,irho+1,ix,  which)
     $                 +      rhohigh*eostable(numu,irho,  ix,  which)
                  useeostable1=rholow*eostable(numu,irho+1,ix+1,which)
     $                 +      rhohigh*eostable(numu,irho,  ix+1,which)
                  useeostable=
     $                 (xlow*useeostable1+xhigh*useeostable0)
     $                 /(steprho*stepx)
               else
                  useeostable=(rholow*eostable(numu,irho+1,1,which)
     $                 +      rhohigh*eostable(numu,irho, 1,which))/steprho
               endif
            endif
            
         endif
      else
c     at extreme densities we will use ideal gas + radiation pressure
c     The most common way to get into this block of code is for a particle to
c     have a lower density than the lowest density in the eosfile.  However,
c     this code is the backstop whenever the density or specific internal
c     energy of a particle is outside the bounds of the EOS table in eosfile.

         if(irho.lt.1) then
            irho=1
         else if(irho.gt.numrho) then
            irho=numrho
         endif

         if(which.eq.2) then
c     The SPH code doesn't actually ever call this routine with which=2, but this
c     block of code is retained in case someone ever does want to ask for mean
c     molecular weight.
            if(iu.lt.1 .or. iu.ge.numu) then
               iuuse = min(max(1,iu),numu)
               if(numx.gt.1) then
c                 This is for when the eosfile contains data for multiple compositions
                  meanmu0=eostable(iuuse,irho,ix,  2)
                  meanmu1=eostable(iuuse,irho,ix+1,2)
                  meanmu=(xlow*meanmu1+xhigh*meanmu0)/stepx
               else
c                 This is for when the eosfile contains data for only one composition
                  meanmu=eostable(iuuse,irho,1,2)
               endif
            else
               ulow=log10u-(utable1+(iu-1)*stepu)
               uhigh=utable1+iu*stepu-log10u
               if(numx.gt.1) then
                  meanmu0=ulow*eostable(iu+1,irho,ix,  2)
     $                 + uhigh*eostable(iu,  irho,ix,  2)
                  meanmu1=ulow*eostable(iu+1,irho,ix+1,2)
     $                 + uhigh*eostable(iu,  irho,ix+1,2)
                  meanmu=(xlow*meanmu1+xhigh*meanmu0)/(stepu*stepx)
               else
                  meanmu=(ulow*eostable(iu+1,irho, 1,2)
     $                 + uhigh*eostable(iu,  irho,1,2))/stepu
               endif
            endif
            useeostable=meanmu
            return
         endif

         if(iu.lt.1 .or. iu.ge.numu) then
            iuuse = min(max(1,iu),numu)
            if(numx.gt.1) then
               meanT4_0=eostable(iuuse,irho,ix,  1)**4
               meanT4_1=eostable(iuuse,irho,ix+1,1)**4
               T4boundary=(xlow*meanT4_1+xhigh*meanT4_0)/stepx
            else
               T4boundary=eostable(iuuse,irho,1,1)**4
            endif
            uboundary=10d0**(utable1+(iuuse-1)*stepu)
         else
            ulow=log10u-(utable1+(iu-1)*stepu)
            uhigh=utable1+iu*stepu-log10u
            if(numx.gt.1) then
               meanT4_0=ulow*eostable(iu+1,irho,ix,  1)**4
     $                +uhigh*eostable(iu,  irho,ix,  1)**4
               meanT4_1=ulow*eostable(iu+1,irho,ix+1,1)**4
     $                +uhigh*eostable(iu,  irho,ix+1,1)**4
               T4boundary=(xlow*meanT4_1+xhigh*meanT4_0)/(stepu*stepx)
            else
               T4boundary=(ulow*eostable(iu+1,irho,1,1)**4
     $                  + uhigh*eostable(iu,  irho,1,1)**4)/stepu
            endif
            uboundary=ucgs
         endif
         Tboundary = T4boundary**0.25d0

         rhoboundary = 10d0**(rhotable1+(irho-1)*steprho)
         ugasboundary = uboundary - arad*T4boundary/rhoboundary
!        oneovergammaminus1 = ugasboundary/(boltz*Tboundary/meanmu)
!        newqconst = oneovergammaminus1*boltz/arad
!        newqconst = ugasboundary*meanmu/(arad*Tboundary)
         newfirstargument = ugasboundary*rhocgs/(arad*Tboundary)

         If(newfirstargument.gt.0) then
!     In previous verion, we assumed monatomic gas, so
!     qconst=1.5d0*boltz/arad
!     call gettemperature(qconst*rhocgs/meanmu,
!     $        -ucgs*rhocgs/arad,temperature)
            
!     This next line finds the temperature T such that T = Tboundary*ugas/ugasboundary
!     We have to use gettemperature because ugas is itself dependent on T.  Specifically,
!     ugas = ucgs - arad*T^4/rhocgs
            call gettemperature(newfirstargument,
     $           -ucgs*rhocgs/arad,temperature)
         else
            ! ugasboundary=0 when radiation pressure dominates, so presumably radiation pressure still dominates here as well
            temperature = (ucgs*rhocgs/arad)**0.25d0
         endif

         if(which.eq.1) then
            useeostable=temperature
            return
         endif

c     If we get to this portion of the code, then which=3.  That is, the user is looking
c     for pressure.

c     In the old code we had the following, which assumes that meanmu is known.
c         pgas=rhocgs*boltz*temperature/meanmu
c     Here we want to pgas without knowning explicitly the mu value.
c     We will start by getting pgas boundary, and we will then get pgas assuming
c     that mu is the same at the two different positions in parameter space.
         if(iu.lt.1 .or. iu.ge.numu) then
            iuuse = min(max(1,iu),numu)
            if(numx.gt.1) then
               meanP0=eostable(iuuse,irho,ix,  3)
               meanP1=eostable(iuuse,irho,ix+1,3)
               Pboundary=(xlow*meanP1+xhigh*meanP0)/stepx
            else
               Pboundary=eostable(iuuse,irho,1,3)
            endif
         else
            ulow=log10u-(utable1+(iu-1)*stepu)
            uhigh=utable1+iu*stepu-log10u
            if(numx.gt.1) then
               meanP0=ulow*eostable(iu+1,irho,ix,  3)
     $              +uhigh*eostable(iu,  irho,ix,  3)
               meanP1=ulow*eostable(iu+1,irho,ix+1,3)
     $              +uhigh*eostable(iu,  irho,ix+1,3)
               Pboundary=(xlow*meanP1+xhigh*meanP0)/(stepu*stepx)
            else
               Pboundary=(ulow*eostable(iu+1,irho, 1,3)
     $              + uhigh*eostable(iu,  irho,1,3))/stepu
            endif
         endif
         pgasboundary = Pboundary - arad*Tboundary**4/3
         if(pgasboundary.gt.0) then
            pgas = pgasboundary*temperature/Tboundary*rhocgs/rhoboundary
         else
            ! the table's gas pressure is lost in its rounding where radiation dominates
            ! (P is given to 4 figures), so presumably radiation pressure still dominates here as well
            pgas = 0d0
         endif
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
