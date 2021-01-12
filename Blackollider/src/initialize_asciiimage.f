      subroutine asciiimage
      include 'starsmasher.h'
      integer i,j,imax,jmax,jmaxmax,ip,ix,iy
c      parameter(imax=900,jmax=229)
      parameter(imax=900,jmaxmax=999)
      character *1 line(imax)
      real*8 density(imax,jmaxmax)
      real*8 filledcell,xtry,ytry

      filledcell=0
      open(12,file='sph.image')
      do j=1,jmaxmax
         read(12,'(900a1)',end=20)(line(i), i=1,imax)
         do i=1,imax
            if(line(i).eq.' ') then
               density(i,j)=0d0
            else
               density(i,j)=1d0
               filledcell=filledcell+density(i,j)
            endif
         enddo
      enddo            
 20   close(12)

      jmax=j-1
      if(myrank.eq.0) then
         write(69,*) 'jmax=',jmax
         write(69,*) 'filledcell=',filledcell
      endif

      ip=0
      do ix=1,imax
         do iy=1,jmax
            xtry=ix
            ytry=jmax-iy
            if(density(ix,iy).gt.0d0)then
               ip=ip+1
               x(ip)=xtry
               y(ip)=ytry
               z(ip)=0d0
               vx(ip)=(xtry/imax-0.5d0)/(0.5d0*imax*jmax)**0.25
               vy(ip)=(ytry/jmax-0.5d0)/(0.5d0*imax*jmax)**0.25
               vz(ip)=0d0
               u(ip)=1d0/sqrt(0.5d0*imax*jmax)
               vxdot(ip)=0d0
               vydot(ip)=0d0
               vzdot(ip)=0d0
               udot(ip)=0d0
               hp(ip)=0.5d0*
     $              (3d0*1.9d0*nnopt/4d0/pi*imax*jmax/filledcell)**(1d0/3d0)
               meanmolecular(ip)=
     $              1.67262158d-24/(2*0.7d0+0.75d0*0.28d0+0.5d0*0.02d0)
            endif
         enddo
      enddo
      n=ip
      ntot=n
      do ip=1,n
         am(ip)=1d0/n
      enddo

      if(myrank.eq.0) write(69,*) 'number of particles n=',n

      call lfstart

      end
