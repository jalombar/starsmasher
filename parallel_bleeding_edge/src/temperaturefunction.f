      double precision function temperaturefunction(temi)

      include 'starsmasher.h'
      real*8 rhoarray(kdm),uarray(kdm),rarray(kdm),
     $     muarray(kdm),pres(kdm)
      integer i
      real*8 temi
      common/splinestuff/rarray,uarray,muarray,rhoarray
      common/presarray/ pres,i

      temperaturefunction=1.d0-(rhoarray(i)*boltz*temi/muarray(i)
     $     +arad*temi**4/3.d0)/pres(i)

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      double precision function ufunction(ucgsi)

      include 'starsmasher.h'
      real*8 rhoarray(kdm),uarray(kdm),rarray(kdm),
     $     muarray(kdm),pres(kdm)
      integer i
      real*8 ucgsi,useeostable
      common/splinestuff/rarray,uarray,muarray,rhoarray
      common/presarray/ pres,i

      ufunction=1.d0-useeostable(ucgsi,rhoarray(i),muarray(i),3)
     $     /pres(i)

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      double precision function ueqfunction(ucgsi)

      real*8 ucgsi,useeostable
      real*8 rhocgs,teq,mucgs
      common/ueqstuff/rhocgs,teq,mucgs

      ueqfunction=1.d0-useeostable(ucgsi,rhocgs,mucgs,1)/teq

      return
      end
