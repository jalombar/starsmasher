      subroutine sph_splint(xa,ya,y2a,n,x,y)
      implicit none
      integer n
      real*8 x,y,xa(n),y2a(n),ya(n)
      integer k,khi,klo
      real*8 a,b,h
      klo=1
      khi=n
 1    if (khi-klo.gt.1) then
         k=(khi+klo)/2
         if(xa(k).gt.x)then
            khi=k
         else
            klo=k
         endif
         goto 1
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.d0) then
         write(69,*)'bad xa input in sph_splint'
         stop
      endif
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a**3.d0-a)*y2a(klo)+(b**3.d0-b)*y2a(khi))
     $     *(h**2)/6.d0
      return
      end
