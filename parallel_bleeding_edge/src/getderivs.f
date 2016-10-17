      subroutine getderivs(i,curlvxi,curlvyi,curlvzi,divvi)
c calculate curl and divergence of velocity
      include 'starsmasher.h' 
      double precision dwijin,curlxin,curlyin,curlzin

      double precision hpi,h2,h5,rhoi,r2
      double precision curlvxi,curlvyi,curlvzi,divvi
      integer i,in,itab,j
      double precision divin,tijin,ctaboverh2

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
      do in=1,nn(i)
         j=list(first(i)+in)
         if(u(j).ne.0.d0) then
            r2=(x(i)-x(j))**2.d0+(y(i)-y(j))**2.d0+(z(i)-z(j))**2.d0
            itab=int(ctaboverh2*r2)+1
            dwijin=dwtab(itab)
            curlxin=(z(i)-z(j))*(vy(i)-vy(j))
     $           -(y(i)-y(j))*(vz(i)-vz(j))
            curlyin=(x(i)-x(j))*(vz(i)-vz(j))
     $           -(z(i)-z(j))*(vx(i)-vx(j))
            curlzin=(y(i)-y(j))*(vx(i)-vx(j))
     $           -(x(i)-x(j))*(vy(i)-vy(j))
            curlvxi=curlvxi+am(j)*dwijin*curlxin
            curlvyi=curlvyi+am(j)*dwijin*curlyin
            curlvzi=curlvzi+am(j)*dwijin*curlzin
            
            divin=(x(i)-x(j))*(vx(i)-vx(j))+
     $           (y(i)-y(j))*(vy(i)-vy(j))+
     $           (z(i)-z(j))*(vz(i)-vz(j))
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








