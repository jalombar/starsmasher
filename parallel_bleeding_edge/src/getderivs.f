      subroutine getderivs(i,curlvxi,curlvyi,curlvzi,divvi)
c calculate curl and divergence of velocity
      include 'starsmasher.h' 
      double precision dwijin,curlxin,curlyin,curlzin

      double precision hpi,h2,h5,rhoi,r2
      double precision curlvxi,curlvyi,curlvzi,divvi
      integer i,in,itab,j
      double precision divin,tijin,ctaboverh2
      double precision dx,dy,dz

      hpi=hp(i)
      h2=hpi**2
      h5=hpi*h2**2
      rhoi=rho(i)
      ctaboverh2=ctab/h2
c     compute the 3 components of the curl of the velocity:
c     curlvx(i), curlvy(i) and curlvz(i)
c     "gather" part of the sum (i-j pair contributes to curl v_i):
      curlvxi=0.d0
      curlvyi=0.d0
      curlvzi=0.d0
      divvi=0.d0
c     a point particle with no gas inside its kernel has rho=0, and the
c     divisions below would give 0/0.  Its div v and curl v are zero.
      if(rhoi.le.0.d0) return
      do in=1,nn(i)
         j=list(first(i)+in)
         if(u(j).ne.0.d0) then
            dx=x(i)-x(j)
            dy=y(i)-y(j)
            dz=z(i)-z(j)
            if(periodic) call pbcwrap(dx,dy,dz)
            r2=dx**2+dy**2+dz**2
            itab=int(min(ctaboverh2*r2,4.d0*ctab))+1
            dwijin=dwtab(itab)
            curlxin=dz*(vy(i)-vy(j))
     $           -dy*(vz(i)-vz(j))
            curlyin=dx*(vz(i)-vz(j))
     $           -dz*(vx(i)-vx(j))
            curlzin=dy*(vx(i)-vx(j))
     $           -dx*(vy(i)-vy(j))
            curlvxi=curlvxi+am(j)*dwijin*curlxin
            curlvyi=curlvyi+am(j)*dwijin*curlyin
            curlvzi=curlvzi+am(j)*dwijin*curlzin
            
            divin=dx*(vx(i)-vx(j))+
     $           dy*(vy(i)-vy(j))+
     $           dz*(vz(i)-vz(j))
            tijin=dwijin*divin
            divvi=divvi-am(j)*tijin
         endif
      enddo
      curlvxi=curlvxi/(rhoi*h5)
      curlvyi=curlvyi/(rhoi*h5)
      curlvzi=curlvzi/(rhoi*h5)
      divvi=divvi/(rhoi*h5)

      return
      end








