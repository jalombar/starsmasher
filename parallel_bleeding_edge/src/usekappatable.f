      subroutine readinkappatable
      include 'starsmasher.h'
      integer numrho,numtem,item,irho
      integer maxtablesize
      parameter(maxtablesize=1000)
      real*8 rhotable(maxtablesize),temtable(maxtablesize)
      real kappatable(maxtablesize,maxtablesize)
      real localkappatable(maxtablesize,maxtablesize)
      real*8 steprho,steptem,rhotable1,temtable1
      common/kappacom/ numrho,numtem,rhotable1,temtable1,
     $     steprho,steptem,kappatable,localkappatable
      
      open (43,file=opacityfile)
      read(43,*)
      read(43,*)
      read(43,*) numrho,numtem

      if(myrank.eq.0) then
         write(69,*)'reading ',trim(opacityfile),
     $        ': numrho,numtem=',
     $        numrho,numtem
      endif
      
      do irho=1,numrho
         do item=1,numtem
            read(43,*) rhotable(irho),temtable(item),
     $           localkappatable(irho,item),
     $           kappatable(irho,item)
         enddo
      enddo

      rhotable1=rhotable(1)
      temtable1=temtable(1)

      steprho=(rhotable(numrho)-rhotable(1))/(numrho-1)
      steptem=(temtable(numtem)-temtable(1))/(numtem-1)
      print *,'steprho,steptem=',steprho,steptem

      end
      
      subroutine usetable(rho,tem,localkappa,pseudokappa)
      implicit none
      real*8 rho,tem
      real*8 pseudokappa,localkappa
      integer item,irho
      integer numrho,numtem
      real*8 temtable1,rhotable1
      integer maxtablesize
      parameter(maxtablesize=1000)
      real kappatable(maxtablesize,maxtablesize)
      real localkappatable(maxtablesize,maxtablesize)

      real*8 steprho,steptem
      common/kappacom/ numrho,numtem,rhotable1,temtable1,
     $     steprho,steptem,kappatable,localkappatable

      real*8 f00,f01,f10,f11,log10rho,log10tem,
     $     rholow,rhohigh,temlow,temhigh

c     temtable(item)=temtable(1)+(item-1)*steptem
c     so, item= (temtable(item)-temtable(1))/steptem + 1

      log10rho=log10(rho)
      log10tem=log10(tem)

      irho = min(max(1,int((log10rho-rhotable1)/steprho) + 1),numrho-1)
      item = min(int((log10tem-temtable1)/steptem) + 1 ,numtem-1)

      rholow=log10rho-(rhotable1+(irho-1)*steprho)
      rhohigh=rhotable1+irho*steprho-log10rho
      
ccc      write(90,*)rho,irho,rholow,rhohigh
ccc      write(90,*)tem,item

      if(item.ge.1) then

         temlow=log10tem-(temtable1+(item-1)*steptem)
         temhigh=temtable1+item*steptem-log10tem
         
c     use bi-linear interpolation among the four cartesian
c     grid point (irho,item), (irho+1,item), (irho,item+1), and (irho+1,item+1)
         f00=rholow*temlow
         f10=rhohigh*temlow
         f01=rholow*temhigh
         f11=rhohigh*temhigh
         
         pseudokappa=(f00*kappatable(irho+1,item+1)
     $        + f10*kappatable(irho,  item+1)
     $        + f01*kappatable(irho+1,item)
     $        + f11*kappatable(irho,  item))/(steprho*steptem)
         localkappa=(f00*localkappatable(irho+1,item+1)
     $        + f10*localkappatable(irho,  item+1)
     $        + f01*localkappatable(irho+1,item)
     $        + f11*localkappatable(irho,  item))/(steprho*steptem)

      else
         pseudokappa=(rholow*kappatable(irho+1,1)
     $        + rhohigh*kappatable(irho,  1))/steprho
     $        *(tem/10d0**temtable1)**2
         localkappa=(rholow*localkappatable(irho+1,1)
     $        + rhohigh*localkappatable(irho,  1))/steprho
     $        *(tem/10d0**temtable1)**2

      endif

ccc      write(90,*) pseudokappa,localkappa

      end
