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
      real*8 steprhotest,steputest

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
      integer which

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
         write(69,*)'xtable1=',xtable1
         write(69,*)'stepx=',stepx
         write(69,*)'numx=',numx
         write(69,*)'particle X=',xxx
         write(69,*)'ix=',ix
         write(69,*)'Expand X range covered by EOS file'
         stop
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

         if(irho.lt.1) irho=1
         if(irho.gt.numrho) irho=numrho

         if(iu.lt.1) then
            if(numx.gt.1) then
               meanmu0=eostable(1,irho,ix,  2)
               meanmu1=eostable(1,irho,ix+1,2)
               meanmu=(xlow*meanmu1+xhigh*meanmu0)/stepx
            else
               meanmu=eostable(1,irho,1,2)
            endif
         else if (iu.ge.numu) then
            if(numx.gt.1) then
               meanmu0=eostable(numu,irho,ix  ,2)
               meanmu1=eostable(numu,irho,ix+1,2)
               meanmu=(xlow*meanmu1+xhigh*meanmu0)/stepx
            else
               meanmu=eostable(numu,irho,1,2)
            endif   
         else
            ulow=log10u-(utable1+(iu-1)*stepu)
            uhigh=utable1+iu*stepu-log10u
            if(numx.gt.1) then
               meanmu0=ulow*eostable(iu+1,irho,ix,  2)
     $              + uhigh*eostable(iu,  irho,ix,  2)
               meanmu1=ulow*eostable(iu+1,irho,ix+1,2)
     $              + uhigh*eostable(iu,  irho,ix+1,2)
               meanmu=(xlow*meanmu1+xhigh*meanmu0)/(stepu*stepx)
            else
               meanmu=(ulow*eostable(iu+1,irho, 1,2)
     $              + uhigh*eostable(iu,  irho,1,2))/stepu
            endif
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
