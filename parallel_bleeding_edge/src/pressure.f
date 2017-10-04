      subroutine pressure
      include 'starsmasher.h'
      real*8 pgas,prad
      real*8 rhocgs,ucgs,beta1,temperature,gam1,useeostable
      integer i

      if(neos.eq.0) then
         if(nintvar.eq.1) then
c     do i=1,ntot
            do i=n_lower,n_upper
c     p=a*rho^gam, so p/rho^2=a*rho^(gam-2.d0)
               por2(i)=u(i)*rho(i)**(gam-2.d0)
            enddo
         else
            do i=n_lower,n_upper
c     p=(gam-1)*rho*u, so p/rho^2=(gam-1)*u/rho
               if(u(i).ne.0.d0) por2(i)=(gam-1)*u(i)/rho(i)
            enddo
         endif
      else if(neos.eq.1) then
c         do i=1,ntot
         do i=n_lower,n_upper
            if(u(i).ne.0.d0) then
               rhocgs=rho(i)*munit/runit**3.d0
               if(nintvar.eq.1) then
                  ucgs=u(i)*rho(i)**(gam-1.d0)/(gam-1.d0)
     $                 *gravconst*munit/runit
               else
                  ucgs=u(i)*gravconst*munit/runit
               endif
               call gettemperature(qconst*rhocgs/meanmolecular(i),
     $              -ucgs*rhocgs/arad,temperature)
               pgas=rhocgs*boltz*temperature/meanmolecular(i)
               prad=arad*temperature**4/3.d0
               beta1=pgas/(pgas+prad)
               gam1=(32.d0-24.d0*beta1-3.d0*beta1**2) /
     $              (24.d0-21.d0*beta1)
               
               if(gam1.lt.0.999*4.d0/3.d0 .or. gam1.gt.1.001*5.d0/3.d0) then
                  write(69,*)'warning gam1=',gam1,'at i=',i
                  write(69,*) beta1,pgas,prad,temperature,rho(i),
     $                 meanmolecular(i),rhocgs,ucgs,x(i),y(i),z(i),hp(i),
     $                 qconst*rhocgs/meanmolecular(i),
     $                 -ucgs*rhocgs/arad
                  stop 'gam1 value does not make sense'
               endif               
               por2(i)=(pgas+prad)/rho(i)**2/punit
            else
               por2(i)=0.d0
            endif
         enddo
      else if(neos.eq.2) then
c     use tabulated eos here!

c     we need to use the tabulated
c     eos to get por2 (pressure over rho squared) from a table that uses specific
c     internal energy u and density rho as input variables.  my recommendation is to
c     make such a table for a gam=2 eos.  then make sure the code behaves the same
c     when neos=0 and when neos=2.

         do i=n_lower,n_upper
            if(u(i).ne.0.d0) then
               rhocgs=rho(i)*munit/runit**3.d0
               if(nintvar.eq.1) then
                  ucgs=u(i)*rho(i)**(gam-1.d0)/(gam-1.d0)
     $                 *gravconst*munit/runit
               else
                  ucgs=u(i)*gravconst*munit/runit
               endif
               por2(i)=useeostable(ucgs,rhocgs,meanmolecular(i),3)
     $              /rho(i)**2/punit
            else
               por2(i)=0.d0
            endif
         enddo

      endif

      return
      end
