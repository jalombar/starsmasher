      subroutine resplintmu
      include 'starsmasher.h'
      real*8 rtry,amass,radius
      integer numlines,i
      real*8 rhoarray(kdm),uarray(kdm),rarray(kdm),
     $     rhoarray2(kdm),uarray2(kdm),muarray(kdm),muarray2(kdm)
      real*8 integratednum,meanmolecularexi,maxmu,minmu
      common/splinestuff/rarray,uarray,muarray,rhoarray,
     $     uarray2,muarray2,rhoarray2,amass,radius,
     $     integratednum,maxmu,minmu,numlines

      do i=1,n
         rtry=sqrt(x(i)**2.d0+y(i)**2.d0+z(i)**2.d0)
         call sph_splint(rarray,muarray,muarray2,numlines,rtry,
     $        meanmolecularexi)
         if(rtry.ge.rarray(numlines)) meanmolecularexi=muarray(numlines)
c         meanmolecular(i)=(1.d0-2.d0*dt/trelax)*meanmolecular(i)+
c     $        2.d0*dt/trelax*meanmolecularexi         
         meanmolecular(i)=meanmolecularexi         
         meanmolecular(i)=min(maxmu,meanmolecular(i))
         meanmolecular(i)=max(minmu,meanmolecular(i))
      enddo

      return
      end
