      subroutine asciiimage
      include 'starsmasher.h'
      integer i,j,imax,jmax,jmaxmax,ip,ix,iy
c      parameter(imax=900,jmax=229)
      parameter(imax=900,jmaxmax=999)
      character *1 line(imax)
      real*8 density(imax,jmaxmax)

      ip=0
      open(12,file=imagefile)
      do j=1,jmaxmax
         read(12,'(900a1)',end=20)(line(i), i=1,imax)
         do i=1,imax
            if (line(i) == ' ') then
               density(i,j) = 0d0
            else if (line(i) == '.') then
               density(i,j) = 0.15d0
            else if (line(i) == ':') then
               density(i,j) = 0.25d0
            else if (line(i) == '-') then
               density(i,j) = 0.40d0
            else if (line(i) == '=') then
               density(i,j) = 0.55d0
            else if (line(i) == '+') then
               density(i,j) = 0.60d0
            else if (line(i) == '*') then
               density(i,j) = 0.75d0
            else if (line(i) == '#') then
               density(i,j) = 0.92d0
            else if (line(i) == '%') then
               density(i,j) = 0.98d0
            else if (line(i) == '@') then
               density(i,j) = 1.00d0
            else if (line(i) == '$') then
               density(i,j) = 0.97d0
            else if (line(i) == '&') then
               density(i,j) = 0.90d0
            else
               density(i,j) = 0.5d0
            endif
            if(density(i,j).gt.0) then
               ip = ip + 1
            endif

         enddo
      enddo            
 20   close(12)
      n = ip
      ntot = n


      jmax=j-1
      if(myrank.eq.0) then
         write(69,*) 'number of particles n=',n
         write(69,*) 'jmax=',jmax
         call flush(69)
      endif

      ip=0
      do ix=1,imax
         do iy=1,jmax
            if(density(ix,iy).gt.0d0)then
               ip=ip+1
               x(ip)=ix * 0.55d0
               y(ip)=jmax-iy
               z(ip)=0d0
               vx(ip)=(dble(ix)/imax-0.5d0)/(0.5d0*imax*jmax)**0.25
               vy(ip)=(y(ip)/jmax-0.5d0)/(0.5d0*imax*jmax)**0.25
               vz(ip)=0d0
               u(ip)=1d0/sqrt(0.5d0*imax*jmax)
               vxdot(ip)=0d0
               vydot(ip)=0d0
               vzdot(ip)=0d0
               udot(ip)=0d0
               hp(ip)=0.55d0 * 0.5d0*
     $              (3d0*1.9d0*nnopt/4d0/pi*imax*jmax/n)**(1d0/3d0)
               meanmolecular(ip)=
     $              1.67262158d-24/(2*0.7d0+0.75d0*0.28d0+0.5d0*0.02d0)
               am(ip)=density(ix,iy)/n
            endif
         enddo
      enddo

      call lfstart

      end
