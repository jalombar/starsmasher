      subroutine polymces
c     creates a star from the data file yrec output
      include 'starsmasher.h'
      integer i
      integer idumb,ip,ix,iy,iz!,maxtry
c      parameter(maxtry=300000000)
      real*8 anumden,rhotry,rhoex,rtry,rhomax,hc,xcm,ycm,zcm,amtot,
     $     ammin,ammax,xtry,ytry,ztry,ri,rhoi,ran1
      integer irtry
      real*8 amass,radius
      integer ixmax,iymax,izmax,corepts
      double precision cellvolume,a1
      real*8 integral
!      real*8 npoly
      integer nrgrid
      parameter(nrgrid=10000)                                         
      real*8 rarray(nrgrid)
      real*8 rhoarray(nrgrid)
      real*8 rhoarray2(nrgrid)
      real*8 ak
c      real*8 rarrayi,rarrayim1
      integer numlines
      real*8 coremass,amasstot,hmax,hmin
      
!      npoly=1.d0/(gam-1.d0)
      amasstot=starmass
      radius=starradius
      
      if(equalmass.ne.0.d0) then
         write(69,*)'right now equalmass must be zero'
         stop
      endif

c     get profiles:
      call polymc(npoly,amasstot,radius,nrgrid,rhoarray,rarray,ak,
     $     numlines,corepts,coremass)
      call sph_spline(rarray,rhoarray,numlines,1.d30,1.d30,rhoarray2)
c     call polymc(npoly,sngl(amass),sngl(radius),nrgrid,rhopol,ak)
      amass=amasstot-coremass
      idumb=-2391
      rhomax=rhoarray(1)
      if(myrank.eq.0) then
         write(69,*)'innermost density value=',rhomax,'at r=',rarray(1)
         write(69,*)'outermost density=',rhoarray(numlines),'at r=',
     $        rarray(numlines)
      endif
      ip=0               
c     corepts=0
      if(equalmass.lt.1.d0)then
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     make an hcp lattice
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         redge=radius
c     do i=1,20         ! seems to work better when integrating u
         do i=1,1               ! seems to work better when integrating a
            hc=((nnopt*39.d0/22.d0)/(8.d0*n))**(1.d0/3.d0)*redge
            redge=radius-3.d0*hc
            if(myrank.eq.0) write(69,*)'redge=',redge
         enddo
         if(myrank.eq.0) then
            write(69,*)'keeping particles up to a distance',radius-redge,
     $           'less than the full radius',radius,'(subrtracted 3*hc)'
            
c     the fraction of particles at a region of density rhoex that
c     will be kept is (rhoex/rhomax)**equalmass.  so if the number
c     density of lattice points that will be tried is n, then we expect
c     the number of particles to be n=n*integrate[4*pi*r**2*
c     (rhoex(r)/rhomax)**equalmass,{r,0,redge}].  we can solve this for
c     n, and then use that the cell volume is 2/n (as there are two
c     particles per cell)
            
            write(69,*)'maybe integral=',4.d0/3.d0*pi*radius**3
         endif
         integral=4.d0/3.d0*pi*redge**3
         
         if(myrank.eq.0)write(69,*)'using integral=',integral
         cellvolume=2.d0*integral/n
         a1=(cellvolume/2.d0**0.5d0)**(1.d0/3.d0)
         if(myrank.eq.0)write(69,*)'a1=',a1
         
c     looking at figure 9(b) and page 18 of kittel
c     (a1 vector)=-0.5*a1*(x hat)-3^0.5/2*a1*(y hat)
c     (a2 vector)=a1*(x hat)
c     (a3 vector)=(8/3)^0.5*a1*(z hat)
         ixmax=int(redge/a1)+2
         iymax=int(redge/(3.d0**0.5d0/2.d0*a1))+2
         izmax=int(redge/(0.5d0*(8.d0/3.d0)**0.5d0*a1))+2
         do ix=-ixmax,ixmax
            do iy=-iymax,iymax
               do iz=-izmax,izmax
                  xtry=(ix-0.5d0)*a1+mod(abs(iy),2)*0.5d0*a1
                  ytry=iy*3.d0**0.5d0/2.d0*a1
     $                 -(mod(abs(iz),2)-0.5d0)*1.d0/3.d0**0.5d0*a1
                  ztry=(iz-0.5d0)*0.5d0*(8.d0/3.d0)**0.5d0*a1
                  rtry=sqrt(xtry**2.d0+ytry**2.d0+ztry**2.d0)
                  if(rtry.lt.redge) then
c     rhoex=dble(
c     $                       rhopol(1+int(rtry/radius*dble(nrgrid-1))))
                     if(rtry.lt.rarray(1) .or.
     $                    rtry.gt.rarray(numlines)) then
                        write(69,*)'worry about extrapolating at radii',
     $                       rtry,rarray(1),rarray(numlines)
                        stop
                     endif
                     
                     call sph_splint(rarray,rhoarray,rhoarray2,numlines,
     $                    rtry,rhoex)
                     rhotry=ran1(idumb)      
                     if (rhotry.le.rhoex**equalmass) then
c     (particle is accepted)    
                        ip=ip+1                 
                        if(rtry.le.a1) then
                           if(myrank.eq.0) write(69,'(4i5,4e12.4)')ip,ix,iy,iz,
     $                          xtry,ytry,ztry,rtry
                        endif
                        x(ip)=xtry        
                        y(ip)=ytry            
                        z(ip)=ztry
                     endif
                  endif
               enddo
            enddo
         enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     (note that n will be slightly changed!)
         n=ip
         if(myrank.eq.0) write (69,*) 'parent: n=',n
         if (n.gt.nmax) stop 'parent: n>nmax ???'
c     assign particle masses (to represent density):
c     cden=(4.d0*pi)/(3.d0*n)*redge**3.d0
         ammin=1.d30
         ammax=0.d0
         amtot=0.d0
         xcm=0.d0
         ycm=0.d0
         zcm=0.d0
         do i=1,n
            ri=sqrt(x(i)**2+y(i)**2+z(i)**2)
c     rhoi=dble(rhopol(1+int(ri/radius*dble(nrgrid-1))))
            call sph_splint(rarray,rhoarray,rhoarray2,numlines,
     $           ri,rhoi)
            
            if(rhoi.le.0.d0) then
               if(myrank.eq.0) write(69,*)'warning: rho(',i,')<=0 at r=',ri,'???'
            endif
            am(i)=amass/n*(integral*rhoi/amass)**(1.d0-equalmass)
            xcm=xcm+am(i)*x(i)
            ycm=ycm+am(i)*y(i)
            zcm=zcm+am(i)*z(i)
            amtot=amtot+am(i)
            if(nintvar.eq.1) then
               u(i)=dble(ak)*rhoi**(1.d0+1.d0/npoly-gam)
            else
c     p=(gam-1)*rho*u=a*rho^gam, so u=a*rho^(gam-1)/(gam-1)
               u(i)=dble(ak)*rhoi**(gam-1.d0)/(gam-1.d0)
            endif
            meanmolecular(i)=meanmolecularweight ! meanmolecular weight can be set in sph.input
         enddo
         xcm=xcm/amtot
         ycm=ycm/amtot
         zcm=zcm/amtot
         if(myrank.eq.0) then
            write (69,'(a,g11.4,a,g11.4,a)')
     $           'parent: total mass of star was m=',amtot
     $           ,
     $           ', which should equal', amass,'... renormalizing...'
            write(69,'(a,3g15.7)') 'center of mass at',xcm,ycm,zcm
         endif
         open(100,file='sph.passivelyAdvected')                      
         do i=1,n
            am(i)=am(i)/amtot*amass
            if(nintvar.eq.2) then
               u(i)=u(i)/amass*amtot
            endif
            if(am(i).gt.0.01d0*amass .and. myrank.eq.0) then
               write(69,*)'warning: particle',i,'has mass',am(i)
               write(69,*)'x,y,z=',x(i),y(i),z(i)
               write(69,*)'r,rho=',ri,rhoi
            endif
            ammin=min(ammin,am(i))
            ammax=max(ammax,am(i))
            rtry=sqrt(x(i)**2+y(i)**2+z(i)**2)
c     rhoex=dble(rhopol(1+int(rtry/radius*dble(nrgrid-1))))
            call sph_splint(rarray,rhoarray,rhoarray2,numlines,
     $           rtry,rhoex)
            
            
            anumden=rhoex/am(i)
            hp(i)=(3.d0/32.d0/3.1415926535897932384626d0*
     $           (nnopt*39.d0/22.d0)/anumden)**(1.d0/3.d0)

            aa(i)=1d30                                            
            dd(i)=hfloor                                          
            bb(i)=(1/(hp(i)-dd(i)) - 1/aa(i))/rhoex**(1d0/3d0)    
            if(myrank.eq.0) write(100+myrank,*) aa(i),bb(i),dd(i) 

         enddo
         close(100)
         if(myrank.eq.0) then
            write (69,*) 'parent: min particle mass=',ammin
            write (69,*) 'parent: max particle mass=',ammax
         endif

c     because the particles were distributed uniformly in the sphere, the
c     the number density is 3*n/(4*pi*radius**3) and we need to choose
c     the smoothing length hc such that 
c     4*pi*(2*hc)**3/3 * (number density)=nnopt.
c     this gives 8*hc**3 * n/radius**3= nnopt, or:
c     hc=(nnopt/(8.d0*n))**(1.d0/3.d0)*redge
c     (should give a number of nearest neighbors close to nnopt)
         
      endif

      hmin=1.d30
      hmax=0.d0
      do i=1,n
         hmin=min(hmin,hp(i))
         hmax=max(hmax,hp(i))
      enddo
      if(myrank.eq.0) then
         write(69,*)'hmin=',hmin
         write(69,*)'hmax=',hmax
      endif

c     ntot=n+corepts      
      n=n+corepts
      ntot=n

      if(corepts.gt.0)then
c         ntot=n+corepts
         am(ntot)=coremass
         x(ntot)=0.d0
         y(ntot)=0.d0
         z(ntot)=0.d0
         hp(ntot)=hmin
         u(ntot)=0.d0
         meanmolecular(ntot)=0.d0
c         nn(ntot)=0
      endif
      
      call stride_setup
      if(ntot.gt.nmax) then
         write(69,*) 'n is too large'
         stop
      endif
c      call reset_h

      do i=1,n
         cc(i)=int(10000*amasstot)
      enddo
      if(myrank.eq.0)write(69,*) 'using cc=',cc(1),
     $     'for all sph particles'
      cc(ntot)=1
      if(myrank.eq.0) then
         write(69,*)'core point',ntot,'has cc=',cc(ntot)
         write(69,'(5g11.4)') 'i', 'nn(i)','hp(i)'
         do i=1,n
            if(mod(i,1000).eq.1)
     $           write(69,'(5g11.4)') i, nn(i),hp(i)
         enddo
      endif
      call rho_and_h

      if(myrank.eq.0) then
         write(69,*) 'i   density'
         do i=1,numlines,numlines/10
            write(*,*) i,rhoarray(i)
         enddo
         
         write(69,*)'rhoarray(1) and rhoarray(2)=',
     $        rhoarray(1),rhoarray(2)
         write(69,*)'rhoarray(numlines-1),rhoarray(numlines) =',
     $        rhoarray(numlines-1),
     $        rhoarray(numlines)
         write(69,*)'near origin:'
         write(69,*)' r    rho'
         do irtry=0,nint(100*radius),nint(10*radius+.5)
            rtry=irtry/1000.d0
c     rhoi=dble(rhopol(1+int(rtry/radius*dble(nrgrid-1))))
            call sph_splint(rarray,rhoarray,rhoarray2,numlines,
     $           rtry,rhoi)
            write(69,*) rtry,rhoi
         enddo
         write(69,*)'near surface:'
         write(69,*)' r    rho'
         do irtry=nint(90*radius),nint(100*radius),nint(radius+.5)
            rtry=min(irtry/100.d0,radius)
c     rhoi=dble(rhopol(1+int(rtry/radius*dble(nrgrid-1))))
            call sph_splint(rarray,rhoarray,rhoarray2,numlines,
     $           rtry,rhoi)
            write(69,*) rtry,rhoi
         enddo
      endif
      
c     assign velocities (all zero):
      do i=1,ntot
         vx(i)=0.d0
         vy(i)=0.d0
         vz(i)=0.d0
         vxdot(i)=0.d0
         vydot(i)=0.d0
         vzdot(i)=0.d0
         udot(i)=0.d0
      enddo
c     prepare leap-frog scheme for first iteration:

      if(myrank.eq.0) write(69,*)'polymces.f: about to call lfstart'

      call lfstart

c      write(69,*)'polymces.f: about to call reset_h'
c      call reset_h

      hmin=1.d30
      hmax=0.d0
      do i=1,n
         hmin=min(hmin,hp(i))
         hmax=max(hmax,hp(i))
      enddo
c      write(69,*)'hp(core mass)=-hmin=',-hmin
c      hp(ntot)=-hmin
      hp(ntot)=hmin
      if(corepts.gt.0) write(69,*)'hp(core mass)=',hp(ntot)
      if(myrank.eq.0) then
         write(69,*)'hmin=',hmin
         write(69,*)'hmax=',hmax
         
         write(69,*)'polymces: tf=',tf,dtout
         write(69,*) 'polymces: exiting'
      endif
      return

      end
c     edited feb 2006
c     the y variable in this code used to stand for theta.
c     now, y = alpha = xi^2 * theta
**************************************************************
c     this file is part of the starcrash code
c     version 1.0
c
c     copyright 2003, by joshua faber
c     this code is publicly available under the terms of the
c     gnu general public license.  see our documentation, or
c     http://www.gnu.org/licenses/gpl.html for more details.
***************************************************************

c      function poly(an,am,r,nr,rho,ak,ytab1,yptab1)
      subroutine polymc(an,am,r,nr,rhoarray,rarray,ak,numlines,
     $     corepts,coremass)
c     calculates specific entropy a for polytrope with index n,
c     mass am, radius r, and gives density rho at nr points in radius
      implicit none
      real*8 pi,an,a
      integer nrm,nr
      parameter (pi=3.14159265359d0,nrm=200000)
      real*8 andbl,am,r,rhoarray(nr),ak,rarray(nr)
      real*8 xtab(nrm),ybothtab(2,nrm)
      common/polidx/andbl
      real*8 dr,yp1,rhoc,x1,ypc
      integer ir
c     rkqs stuff:
      real*8 htry,eps
      external derivsmc,rkqs
      integer nbad,nok
      real*8 hmin
      integer kmax,kount,nvar
      real*8 dxsav,ystart(2)
      common/path/ kmax,kount,dxsav,xtab,ybothtab
      integer numlines,corepts
      real*8 coremass,ee

      nvar=2
      dxsav=0.d0
      kmax=nrm
      eps=1.d-13
      andbl=an
      hmin=0.0d0
c      do i=1,nr
c         yscal(i)=1.d-30
c      enddo
      
c     set core mass:
c      yp1=0.d0

c      x1=3.65375d0
      x1=1.d0

c     enter core mass during execution:
c      write(69,*)'enter e:'
c      read*,yp1

c      ee=0.19841618098796d0       !yields 0.99000000
c      ee=2.109167603556806d0    !gives a coremass of 0.9
c      ee=3.2463087530355d0      !yields 0.85000000
c      ee=4.4443139716630d0        !yields 0.80000000
c      ee= 5.7095988116149d0     !gives a coremass of 0.75
c      ee=7.0496227929576d0      !yields 0.70000000
c      ee=8.4731249154836d0      !yields 0.65000000
c      ee=9.9904284872328d0      !yields 0.60000000
c      ee=11.613840854792d0        !yields 0.55000000
c      ee=13.3581851250238d0     !gives a coremass of 0.5
c      ee=17.286120424968d0      !yields 0.40000000
c      ee=21.97824502596553d0    !gives a coremass of 0.3
c      ee=24.707161583d0         !gives a coremass of 0.25
c      ee=27.7673143116708d0     !gives a coremass of 0.2
c      ee=29.446145221851407d0   !gives a coremass of 0.175
c      ee=31.2406649999890d0     !gives a coremass of 0.15  
c      ee=33.1660659703d0        !gives a coremass of 0.125
c      ee=35.2403485466562d0     !gives a coremass of 0.1
c      ee=39.9250152980455d0     !gives a coremass of 0.05
      ee=45.480816d0            !gives a coremass of +0


      if(ee.ge.45.4808d0)then
         corepts=0
      else 
         corepts=1
      endif

      yp1=-(ee/(an+1)**an/x1**(an+1))**(1/(an-1))*x1**2

      write(6,*)'slope at surface=',yp1

      dr=x1/dble(nr-1)
      htry=-dr
      
c     y1=yptab(1)*dr+dr**2.d0 !-dr**4.d0/6.d0+an*dr**6.d0/120.d0
c     yp1=yptab(1)+2.d0*dr   !-2.d0/3.d0*dr**3.d0+an*dr**5.d0/20.d0
c     ytab(2)=y1
c     yptab(2)=yp1
c     xtab(2)=dr

      ystart(1)=0.d0
      ystart(2)=yp1
c      ystart(1)=yp1*dr+dr**2.d0 !-dr**4.d0/6.d0+an*dr**6.d0/120.d0
c      ystart(2)=yp1+2.d0*dr     !-2.d0/3.d0*dr**3.d0+an*dr**5.d0/20.d0

c      write(69,*) 'starting with', ystart(1),ystart(2)

      call odeint(ystart,nvar,x1,dr,eps,htry,hmin,nok,nbad,derivsmc,
     $     rkqs)

      write(6,*)'polymc: kount=',kount
      write(6,35)'ir    ','xtab    ','ytab    ','yptab    '
      do ir=1,kount,kount/10
 35      format(5g15.7)
         write(6,*)ir,xtab(ir),
     $        ybothtab(1,ir),ybothtab(2,ir)
      enddo

c 11   continue
cc      xc=xtab(kount)
cc      yc=ybothtab(1,kount)
      ypc=ybothtab(2,kount)
cc     if (ir.eq.nr) goto 21
c      goto 21
c      stop 'poly: no convergence ???'
c 21   continue

      a=r/x1
      rhoc=am/4.d0/pi/a**3/dabs(yp1)
      ak=4*pi*a**2/rhoc**(1.d0/an-1.d0)/(an+1.d0)

      open(26,file='polymc.col', status='unknown')
      do ir=kount,1,-1
         rarray(kount-ir+1)=xtab(ir)*a
         rhoarray(kount-ir+1)=rhoc*(ybothtab(1,ir)/xtab(ir)**2)**an
         write(26,*) rarray(kount-ir+1),
     $        ak*rhoarray(kount-ir+1)**(1.d0+1.d0/an),
     $        rhoarray(kount-ir+1)
      enddo
      close(26)
      numlines=kount

c      rhoc=(4.d0*pi*a**2.d0/ak/(an+1.d0))**(an/(1.d0-an))
c      ak=0.424188107d0
      if(yp1.eq.0.d0)then
         write(6,*)'should be:'
         write(6,*)'xi_1 = 3.65375'
         write(6,*)'y_1 = 0'
         write(6,*)'xi_1*yprime = 2.71406'
         write(6,*)
         write(6,*)'actual:'
         write(6,*)'xi_1 = ',x1
         write(6,*)'alpha_1 = ',ybothtab(1,1)
         write(6,*)'xi_1*alphaprime = ',yp1
         write(6,*)
      
         write(6,*)'k = ',ak
         write(6,*)'e=',(an+1)**an*x1**(an+1)*(-yp1/x1**2)**(an-1)
      endif

      coremass=4.d0*pi*a**3.d0*rhoc*ypc !,am*yp1/dabs(yps)
      write(6,*)'core mass: ',coremass

      if(corepts.eq.0) then
         write(6,*)'changing to coremass=0'
         coremass=0
      endif

      write(6,*)'total mass: ',4.d0*pi*a**3.d0*rhoc*dabs(yp1) !,am
      write(6,*)'total radius: ',r

c     poly=4.d0*pi*a**3.d0*rhoc*(xs*dabs(yps))-1.d0

c      ytab(nr)=0.d0
c      do ir=1,nr
c         rho(ir)=rhoc*ytab(ir)**an
c      enddo

      return
      end
************************************************************************
      subroutine derivsmc(x,v,dv)
      implicit none
      real*8 x,v(2),dv(2),an
      common/polidx/an
      dv(1)=v(2)
      dv(2)=-dabs(v(1))**an*x**(2.d0-2.d0*an)-2.d0*v(1)/x**2+2.d0/x*v(2)
      return
      end
************************************************************************
      subroutine odeint(ystart,nvar,x1,x2,eps,h1,hmin,nok,nbad,derivsmc,
     $     rkqs)
      integer nbad,nok,nvar,nmax,kmaxx,maxstp
      real*8 eps,h1,hmin,x1,x2,ystart(nvar),tiny
      parameter (maxstp=100000,nmax=2,kmaxx=200000,tiny=1.d-30)
      integer i,kmax,kount,nstp
      real*8 dxsav,h,hdid,hnext,x,xsav,dydx(nmax),xp(kmaxx),
     $     y(nmax),yp(nmax,kmaxx),yscal(nmax)
      common /path/ kmax,kount,dxsav,xp,yp
      external derivsmc,rkqs
      x=x1
      h=sign(h1,x2-x1)
      nok=0
      nbad=0
      kount=0
      do i=1,nvar
         y(i)=ystart(i)
      enddo
      if (kmax.gt.0) xsav=x-2.d0*dxsav
      do nstp=1,maxstp
         call derivsmc(x,y,dydx)
         do i=1,nvar
            yscal(i)=abs(y(i))+abs(h*dydx(i))+tiny
         enddo
         if(kmax.gt.0)then
            if(abs(x-xsav).gt.abs(dxsav)) then
               if(kount.lt.kmax-1)then
                  kount=kount+1
                  xp(kount)=x
                  do i=1,nvar
c                     if(y(i).lt.0.d0)then
c                        write(69,*)'stopping first time in'
c                        write(69,*)i,y(i),kount
c                        return
c                     endif
                     yp(i,kount)=y(i)
                  enddo
                  xsav=x
               endif
            endif
         endif
         if((x+h-x2)*(x+h-x1).gt.0.d0) h=x2-x
         call rkqs(y,dydx,nvar,x,h,eps,yscal,hdid,hnext,derivsmc)
         if(hdid.eq.h)then
            nok=nok+1
         else
            nbad=nbad+1
         endif
         if((x-x2)*(x2-x1).ge.0.d0)then
            do i=1,nvar
               ystart(i)=y(i)
            enddo
            if(kmax.ne.0)then
               kount=kount+1
               xp(kount)=x
               do i=1,nvar
c                  if(y(i).lt.0.d0)then
c                     write(69,*)'stopping second time in'
c                     write(69,*)i,y(i),kount
c                     return
c                  endif
                  yp(i,kount)=y(i)
               enddo
            endif
            return
         endif
         if(abs(hnext).lt.hmin) then
            write(69,*)'stepsize smaller than minimum in odeint'
            stop
         endif
         h=hnext
      enddo
      write(69,*) 'too many steps in odeint'
      stop
      return
      end
************************************************************************
c     5th order runge-kutta with quality control
      subroutine rkqs(y,dydx,n,x,htry,eps,yscal,hdid,hnext,derivsmc)
      integer i,n,nmax
      real*8 eps,hdid,hnext,htry,x,dydx(n),y(n),yscal(n)
      external derivsmc
      parameter(nmax=50)
c     htry is the stepsize to be tried
c     eps is the required accuracy
c     yscal is the error scaling
c     y and x are set to their new values
c     hdid is the stepsize actually used
c     hnext is an estimation for the next stepsize
      real*8 errmax,h,htemp,xnew,yerr(nmax),ytemp(nmax),safety,pgrow,
     $     pshrink,errcon
      parameter(safety=0.9d0,pgrow=-0.2d0,pshrink=-0.25d0,
     $     errcon=1.89d-4)

      h=htry
 10   call rkck(y,dydx,n,x,h,ytemp,yerr,derivsmc)
      errmax=0.d0
      do i=1,n
         errmax=max(errmax,dabs(yerr(i)/yscal(i)))
      enddo
      errmax=errmax/eps
      if(errmax.gt.1)then
         htemp=safety*h*errmax**pshrink
         h=sign(max(dabs(htemp),0.1d0*dabs(h)),h)
         xnew=x+h
         if(xnew.eq.x) then
            write(6,*) 'rkqs: stepsize underflow',x,y
         endif
         goto 10
      else
         if(errmax.gt.errcon)then
            hnext=safety*h*errmax**pgrow
         else
            hnext=5.d0*h
         endif
         hdid=h
         x=x+h
         do i=1,n
            y(i)=ytemp(i)
         enddo
         return
      endif
      end
************************************************************************
c     cash-karp runge-kutta step
      subroutine rkck(y,dydx,n,x,h,yout,yerr,derivsmc)
      implicit none
      integer n,nmax
      real*8 h,x,dydx(n),y(n),yerr(n),yout(n)
      external derivsmc
      parameter(nmax=50)
c     given n variables y and their derivatives dydx
c     advance over h, return values in yout
c     need derivsmc(x,y,dydx) that gives dydx at x
      integer i
      real*8 ak2(nmax),ak3(nmax),ak4(nmax),ak5(nmax),ak6(nmax),
     $     ytemp(nmax),a2,a3,a4,a5,a6,b21,b31,b32,b41,b42,b43,b51,
     $     b52,b53,b54,b61,b62,b63,b64,b65,c1,c3,c4,c6,dc1,dc3,
     $     dc4,dc5,dc6
      parameter(a2=.2d0,a3=.3d0,a4=.6d0,a5=1.d0,a6=.875d0,b21=.2d0,
     $     b31=3.d0/40.d0,b32=9.d0/40.d0,b41=.3d0,b42=-.9d0,b43=1.2d0,
     $     b51=-11.d0/54.d0,b52=2.5d0,b53=-70.d0/27.d0,b54=35.d0/27.d0,
     $     b61=1631.d0/55296.d0,b62=175.d0/512.d0,b63=575.d0/13824.d0,
     $     b64=44275.d0/110592.d0,b65=253.d0/4096.d0,c1=37.d0/378.d0,
     $     c3=250.d0/621.d0,c4=125.d0/594.d0,c6=512.d0/1771.d0,
     $     dc1=c1-2825.d0/27648.d0,dc3=c3-18575.d0/48384.d0,
     $     dc4=c4-13525.d0/55296.d0,dc5=-277.d0/14336.d0,dc6=c6-.25d0)
      do i=1,n
         ytemp(i)=y(i)+b21*h*dydx(i)
      enddo
      call derivsmc(x+a2*h,ytemp,ak2)
      do i=1,n
         ytemp(i)=y(i)+h*(b31*dydx(i)+b32*ak2(i))
      enddo
      call derivsmc(x+a3*h,ytemp,ak3)
      do i=1,n
         ytemp(i)=y(i)+h*(b41*dydx(i)+b42*ak2(i)+b43*ak3(i))
      enddo
      call derivsmc(x+a4*h,ytemp,ak4)
      do i=1,n
         ytemp(i)=y(i)+h*(b51*dydx(i)+b52*ak2(i)+b53*ak3(i)+b54*ak4(i))
      enddo
      call derivsmc(x+a5*h,ytemp,ak5)
      do i=1,n
         ytemp(i)=y(i)+h*(b61*dydx(i)+b62*ak2(i)+b63*ak3(i)+
     $        b64*ak4(i)+b65*ak5(i))
      enddo
      call derivsmc(x+a6*h,ytemp,ak6)
      do i=1,n
         yout(i)=y(i)+h*(c1*dydx(i)+c3*ak3(i)+c4*ak4(i)+c6*ak6(i))
      enddo
      do i=1,n
         yerr(i)=h*(dc1*dydx(i)+dc3*ak3(i)+dc4*ak4(i)+dc5*ak5(i)+
     $        dc6*ak6(i))
      enddo

      return
      end
