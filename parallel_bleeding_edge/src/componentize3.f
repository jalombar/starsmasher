      recursive subroutine componentize3(i,ic)
      include 'starsmasher.h'
      integer, intent(in) :: i
      integer, intent(in) :: ic
      integer in, j
      real*8 am4
      real*8 xcom1,ycom1,zcom1,vxcom1,vycom1,vzcom1,am1
      real*8 xcom2,ycom2,zcom2,vxcom2,vycom2,vzcom2,am2
      real*8 xcom3,ycom3,zcom3,vxcom3,vycom3,vzcom3,am3
      integer icomp(nmax)
      common/compbettercom3/am1,xcom1,ycom1,zcom1,vxcom1,vycom1,vzcom1,
     $     am2,xcom2,ycom2,zcom2,vxcom2,vycom2,vzcom2,
     $     am3,xcom3,ycom3,zcom3,vxcom3,vycom3,vzcom3,am4,
     $     icomp
c      real*8 rhoj(nmax)

      if(ic.eq.4) return

      do in=1,nn(i)
         j=list(first(i)+in)
c         write(69,*)j,'is a neighbor of',i
         if(rho(j).le.rho(i) .and. icomp(j).eq.4) then
c         if(grpot(j).ge.grpot(i) .and. icomp(j).eq.4) then
            icomp(j)=ic
            call componentize3(j,ic)
         endif
      enddo
      end
