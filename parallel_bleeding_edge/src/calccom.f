      subroutine calccom
      include'starsmasher.h'
      real*8 am1,am2,am3,am4
      real*8 x1,y1,z1,vx1,vy1,vz1
      real*8 x2,y2,z2,vx2,vy2,vz2
      real*8 x3,y3,z3,vx3,vy3,vz3
      integer icomp(nmax),i
      common/compbettercom3/am1,x1,y1,z1,vx1,vy1,vz1,
     $     am2,x2,y2,z2,vx2,vy2,vz2,
     $     am3,x3,y3,z3,vx3,vy3,vz3,am4,
     $     icomp

c     determine components cm positions and velocities
      am1=0.d0
      x1=0.d0
      y1=0.d0
      z1=0.d0
      vx1=0.d0
      vy1=0.d0
      vz1=0.d0
      am2=0.d0
      x2=0.d0
      y2=0.d0
      z2=0.d0
      vx2=0.d0
      vy2=0.d0
      vz2=0.d0
      am3=0.d0
      x3=0.d0
      y3=0.d0
      z3=0.d0
      vx3=0.d0
      vy3=0.d0
      vz3=0.d0
      am4=0.d0

c      write(69,*)'calccom: n=',n
      do i=1,n
         if (icomp(i).eq.1) then
            am1=am1+am(i)
            x1=x1+am(i)*x(i)
            y1=y1+am(i)*y(i)
            z1=z1+am(i)*z(i)
            vx1=vx1+am(i)*vx(i)
            vy1=vy1+am(i)*vy(i)
            vz1=vz1+am(i)*vz(i)
         else if (icomp(i).eq.2) then
            am2=am2+am(i)
            x2=x2+am(i)*x(i)
            y2=y2+am(i)*y(i)
            z2=z2+am(i)*z(i)
            vx2=vx2+am(i)*vx(i)
            vy2=vy2+am(i)*vy(i)
            vz2=vz2+am(i)*vz(i)
         else if (icomp(i).eq.3) then
            am3=am3+am(i)
            x3=x3+am(i)*x(i)
            y3=y3+am(i)*y(i)
            z3=z3+am(i)*z(i)
            vx3=vx3+am(i)*vx(i)
            vy3=vy3+am(i)*vy(i)
            vz3=vz3+am(i)*vz(i)
         else if (icomp(i).eq.4) then
            am4=am4+am(i)
         endif
      enddo
      if(am1.gt.0.d0) then
         x1=x1/am1
         y1=y1/am1
         z1=z1/am1
         vx1=vx1/am1
         vy1=vy1/am1
         vz1=vz1/am1
      endif
      if(am2.gt.0.d0) then
         x2=x2/am2
         y2=y2/am2
         z2=z2/am2
         vx2=vx2/am2
         vy2=vy2/am2
         vz2=vz2/am2
      endif
      if(am3.gt.0.d0) then
         x3=x3/am3
         y3=y3/am3
         z3=z3/am3
         vx3=vx3/am3
         vy3=vy3/am3
         vz3=vz3/am3
      endif

      return
      end
