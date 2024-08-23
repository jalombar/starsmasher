      subroutine rescale
c     creates a star with varying equalmass values, spanning from +1 along the +x-axis to -1 along the -x-axis.
      include 'starsmasher.h'
      real*8 divv(nmax)
      common/commdivv/divv
      integer numlines,i
      integer idumb,ip,ix,iy,iz
      real*8 anumden,rhotry,rhoex,rtry,rhomax,hc,xcm,ycm,zcm,amtot,
     $     ammin,ammax,xtry,ytry,ztry,ri,rhoi
      integer irtry
      real*8 amass,masscgs,radius
      real*8 tem(kdm),pres(kdm),
     $     rhoarray(kdm),uarray(kdm),rarray(kdm),
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
      real*8 hmin
      integer numangle
      parameter (numangle=301)
      real*8 integral(numangle)
      real*8 rpreedge(numangle), rpost
      real*8 rprearray(kdm,numangle),rarray2(kdm,numangle)
      real*8 xcenter,ycenter,rcenter
      integer jangle0,ixmin,j,minix,maxix,miniy,maxiy,miniz,maxiz
      real*8 minequalmass,maxequalmass
      parameter (minequalmass=0d0,maxequalmass=1d0)
      integer jangle
      real*8 cosangle,rightedgeoverradius
      parameter (rightedgeoverradius=-0.5d0)
      integer nnoptold,noutold,nitold,
     $     navold,ngrold,nrelaxold,nchk
      real*8 hcoold,hfloorold,sep0old,
     $     tfold,dtoutold,told,
     $     alphaold,betaold,
     $     trelaxold,dtold,rescalefactor

      call splinesetup

      if(myrank.eq.0) write (69,*) 'split: reading startfile1 ', trim(startfile1)

      open(12,file=startfile1,form='unformatted')

      read(12) n,nnoptold,hcoold,hfloorold,sep0old,
     $     tfold,dtoutold,noutold,nitold,told,
     $     navold,alphaold,betaold,tjumpahead,
     $     ngrold,
     $     nrelaxold,trelaxold,dtold,omega2
      amtot=0.d0
      if(myrank.eq.0) write(69,*)'n=',n
      xcm=0
      ycm=0
      zcm=0
      rescalefactor=1.0175d0

      do i=1,n
         read (12) x(i),y(i),z(i),am(i),hp(i),rho(i),vx(i),vy(i),
     $        vz(i),vxdot(i),vydot(i),vzdot(i),u(i),udot(i),
     $        grpot(i),meanmolecular(i),
     $        cc(i),divv(i)
         
         x(i)=x(i)*rescalefactor
         y(i)=y(i)*rescalefactor
         z(i)=z(i)*rescalefactor
         hp(i)=hp(i)*rescalefactor
         rho(i)=rho(i)/rescalefactor**3
         grpot(i)=grpot(i)/rescalefactor
         u(i)=u(i)/rescalefactor

         xcm=xcm+am(i)*x(i)
         ycm=ycm+am(i)*y(i)
         zcm=zcm+am(i)*z(i)
         ! Assume star is supposed to be motionless
         vx(i)=0
         vy(i)=0
         vz(i)=0
         vxdot(i)=0
         vydot(i)=0
         vzdot(i)=0
         udot(i)=0 
         ueq(i)=0
         tthermal(i)=1d30
         amtot=amtot+am(i)
      enddo
      read(12) nchk
      close(12)

      ntot=n

      xcm=xcm/amtot
      ycm=ycm/amtot
      zcm=zcm/amtot
      do i=1,n
         x(i)=x(i)-xcm
         y(i)=y(i)-ycm
         z(i)=z(i)-zcm
      enddo

      if (n+corepts.gt.nmax) stop 'parent: n>nmax ???'

      call stride_setup

c     prepare leap-frog scheme for first iteration:
      call lfstart

      hmin=1.d30
      hmax=0.d0
      do i=1+corepts,ntot
         hmin=min(hmin,hp(i))
         hmax=max(hmax,hp(i))
      enddo

      if(myrank.eq.0)then
         write(69,*)'hmin=',hmin
         write(69,*)'hmax=',hmax
         write(69,*)'rescale: tf=',tf,dtout
         write(69,*)'rescale: exiting'
      endif

      return

      end
