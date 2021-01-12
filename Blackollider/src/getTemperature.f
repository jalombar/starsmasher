      subroutine gettemperature(q,r,x3)
c     subroutine to solve 4th order equations to determine the temperature x3
c     for an equation of state with both ideal gas and radiation pressure.
c     written by scott fleming 10/04/02 and james lombardi 2002-2003

c     the fourth order equation comes from u_gas+ u_rad = u, with
c     u_gas proportional to t and u_rad proportional to t^4

c     in general, we can transform a 4th order equation to x^4+px^2+qx+r=0
c     (see pages 57-58 of stillwell's "mathematics and its history" text)
c     but we fortunately don't even have an x^2 term (that is, p=0).

c     follow stillwell, we can transform this into a cubic equation:
c     first solve for y by using the fact that b^2-4ac=0
c     equation is then:  y^3=ry+q^2/8
c     using the solution of cubic equations found in stillwell page 55:

      implicit none
      real*8 q,r,k,b,piece1,piece2
      real*8 y1,y2,yy,aa,b2,c2,x3,kh
c      real*8 piece2old,x3old

      k = 0.125d0*q**2
      kh=0.5d0*k
      if(kh**2-(r/3.d0)**3.le.0.d0)then
         write(69,*) k,r,kh**2-(r/3.d0)**3
         stop 'bad input: imaginary results?'
      endif

      piece1 = kh+(kh**2-(r/3.d0)**3)**0.5d0
c      piece2old = kh-(kh**2-(r/3.d0)**3)**0.5d0
      piece2 = (r/3.d0)**3.d0/piece1

c      write(69,*)piece2old,piece2

      y1 = piece1**(1.d0/3.d0)
      
c     fortran can't handle cube roots of neg. #'s
      y2 = -dabs(piece2)**(1.d0/3.d0)
      yy = y1+y2

c     equation to solve: (x^2+p+y)^2=ax^2+bx+c
c     now take square root of both sides with:

      aa = 2.d0*yy
      b = -q
c      c = -r+y**2

c     re-writing ax^2+bx+c as a square then solving the equation we
c     obtain 2 results:
c     x^2 + (-(a^(1/2)))x + (-b/(2(a)^(1/2))+p+y) = 0 (1)
c     or
c     x^2 + (a^(1/2))x + (b/(2(a)^(1/2))+p+y) = 0     (2)

c     our solution we're interested in:
      b2 = aa**0.5d0
      c2 = 0.5d0*b/b2 + yy

c     therefore, we once again have x^2+bx+c=0, and our answer we want is:
c      x3old = 0.5d0*(-b2 + (b2**2-4.d0*c2)**0.5d0)
      x3 = -2.d0*c2/(b2 + (b2**2-4.d0*c2)**0.5d0)

c      write(69,*) 'temperature components', x3old, x3

      if(piece1.lt.0.d0) write(69,*)
     $     'piece 1 lt 0',k,r,piece1,piece2
c      if(piece1.eq.-piece2)then
      if(b2.eq.0.d0)then
c         write(69,*)
c     $        'piece 1 eq -piece 2 (rad pressure dominates)',
c     $        k,r,piece1,piece2,b2,c2,x3,(-r)**0.25d0
         x3=(-r -q*( -r-q*(-r)**0.25d0 )**0.25d0)**0.25d0
c         write(69,*) x3
c         write(69,*) piece1,piece2
      endif
      if(piece2.ge.0.d0) x3=-(r+(r/q)**4)/q
      
      end
      
