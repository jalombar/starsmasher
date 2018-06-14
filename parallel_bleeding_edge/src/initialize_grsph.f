      subroutine grsph
c     creates a star from the data file yrec output
      include 'starsmasher.h'

      integer i,j,k, l,p

                                !! sph-particles parameters
      real*8 anumden,amtot,hc

                                !! grid parameters
      real*8 xi,xf, yi,yf, zi,zf
      integer nx,ny,nz,ndata,di,dj,dk
      integer nxmax,nymax,nzmax
      parameter(nxmax=120,nymax=120,nzmax=120)
      integer gridmax
      parameter(gridmax=nxmax*nymax*nzmax)
      real*8 press(gridmax), dens(gridmax),vxgrid(gridmax),
     $     vygrid(gridmax),vzgrid(gridmax)
c      parameter(nx=80,ny=80,nz=40)
      real*8 dx,dy,dz
      real*8 presscartgrid(0:nxmax-1,0:nymax-1,0:nzmax-1)
      real*8 denscartgrid(0:nxmax-1,0:nymax-1,0:nzmax-1)
      integer hcpmax
      parameter(hcpmax=200)     !(hcpmax=135)  !parameter(hcpmax=128)
      real*8 dxhcp(0:hcpmax,0:hcpmax)
      real*8 dyhcp(0:hcpmax,0:hcpmax)
      real*8 dzhcp(0:hcpmax)
      real*8 dphcpgrid(0:hcpmax,0:hcpmax,0:hcpmax)
      real*8 ucartgrid(0:nxmax-1,0:nymax-1,0:nzmax-1)
      real*8 dpcartgrid(0:nxmax-1,0:nymax-1,0:nzmax-1)
      real*8 vxcartgrid(0:nxmax-1,0:nymax-1,0:nzmax-1)
      real*8 vycartgrid(0:nxmax-1,0:nymax-1,0:nzmax-1)
      real*8 vzcartgrid(0:nxmax-1,0:nymax-1,0:nzmax-1)
      real*8 xcartgrid(0:nxmax-1),ycartgrid(0:nymax-1),zcartgrid(0:nzmax-1)
      integer imin,imax,jmin,jmax,kmin,kmax
      integer in,jn,kn
      real*8 fullamtot,fullekin,fulleint,fullepot
      real*8 cellvolume,a1
      integer ixmax,iymax,izmax,ix,iy,iz
      real*8 xtry,ytry,ztry
      real*8 deltax,deltay,deltaz
      real*8 f000,f100,f010,f110,f001,f101,f011,f111,amax
      real*8 t2
      double precision hpmin,hpmax
      integer nnemax,nnemin,nnesig,nneavr
      real*8 xharmonic, yharmonic, zharmonic, rharmonic, ratio, rdotv
      real*8 xbh, ybh, zbh, ambh,tooclose

      t=0.d0
      nselfgravity=0

cc get parameters from input file:
c      open(12,file='sph.input',err=100)
c      read(12,input)
c      close(12)

      if(n.gt.nmax) then
         if(myrank.eq.0)write(69,*)'error: n>nmax'
         stop
      endif

      open(11,file='boundaries.dat', status='old')
      read(11,*) xi
      read(11,*) xf
      read(11,*) yi
      read(11,*) yf
      read(11,*) zi
      read(11,*) zf
      read(11,*) nx
      read(11,*) ny
      read(11,*) nz
      read(11,*) xbh
      read(11,*) ybh
      read(11,*) zbh
      read(11,*) ambh
      close(11)
c      xi = -100.d0      !-140.d0      !-100.d0
c      xf = 150.d0       !100.d0       !150.d0
c      yi = -180.d0      !-100.d0      !-180.d0
c      yf = 110.d0       !180.d0       !110.d0
c      zi = 0.d0
c      zf = 35.d0        !40.d0        !35.d0
c      xi = xmingrid
c      xf = xmaxgrid
c      yi = ymingrid
c      yf = ymaxgrid
c      zi = zmingrid
c      zf = zmaxgrid
      dx = (xf-xi)/(nx-1)
      dy = (yf-yi)/(ny-1)
      dz = (zf-zi)/(nz-1)
      print*,'grid: ',xi,xf,yi,yf,zi,zf

 
      ndata = nx*ny*nz

c      munit=1.9891d33
c      runit=gravconst*munit/crad**2

      if(myrank.eq.0) then
c         write(69,*)'mass and radius unit have been forced to...'
c         write(69,*)'munit=',munit
c         write(69,*)'runit=',runit
c     dump paramteres:
         write(69,*)'{xi,xf,dx} =', xi,xf,dx,'; nx=',nx
         write(69,*)'{yi,yf,dy} =', yi,yf,dy,'; ny=',ny
         write(69,*)'{zi,zf,dz} =', zi,zf,dz,'; nz=',nz      
      endif

c      if(nintvar.ne.1) then
c         nintvar=1
c         if(myrank.eq.0) write(69,*) 'forced nintvar=1'
c      endif

      if(gam.ne.2) then
         if(myrank.eq.0)
     $        write(69,*)'set gam=2 in sph.input, and try again',
     $        gam
         stop
      endif

c     reading data...
      open(11,file='density.dat', status='old')
      open(12,file='pressure.dat', status='old')
      open(13,file='vx.dat', status='old')
      open(14,file='vy.dat', status='old')
      open(15,file='vz.dat', status='old')
      
c     reading the data...
      read(11,*) t, (dens(i),i=1,ndata)
      read(12,*) t2, (press(i),i=1,ndata)
      if(t.ne.t2 .and. myrank.eq.0)
     $     write(69,*)'times do not match (density and pressure)'
      read(13,*) t, (vxgrid(i),i=1,ndata)
      if(t.ne.t2 .and. myrank.eq.0)
     $     write(69,*)'times do not match (pressure and vx)'
      read(14,*) t2, (vygrid(i),i=1,ndata)
      if(t.ne.t2 .and. myrank.eq.0)
     $     write(69,*)'times do not match (vx and vy)'
      read(15,*) t, (vzgrid(i),i=1,ndata)
      if(t.ne.t2 .and. myrank.eq.0)
     $     write(69,*)'times do not match (vy and vz)'

      close(11)
      close(12)
      close(13)
      close(14)
      close(15)

                                !! get the positions x,y,z  --- assuming an uniform cartesian grid... [*]

      do i = 0,nx-1
         xcartgrid(i)=xi+i*dx
      enddo
      do j = 0,ny-1
         ycartgrid(j)=yi+j*dy
      enddo
      do k = 0,nz-1
         zcartgrid(k)=zi+k*dz
      enddo

      amax=0.d0
      fullamtot=0.d0
      fullekin=0.d0
      fulleint=0.d0
      fullepot=0.d0
      l=0
      do k = 0,nz-1
         do j = 0,ny-1
            do i = 0,nx-1
               l=l+1
               if(k.eq.0) then
                  fullamtot=fullamtot+0.5d0*dens(l)
                  fullekin=fullekin+0.5d0*dens(l)*
     $                 (vx(l)**2+vy(l)**2+vz(l)**2)
                  fulleint=fulleint+0.5d0*press(l)/(gam-1)
                  fullepot=fullepot-0.5d0*dens(l)/sqrt(
     $                 (xcartgrid(i)-xbh)**2+
     $                 (ycartgrid(j)-ybh)**2+
     $                 (zcartgrid(k)-zbh)**2)
               else
                  fullamtot=fullamtot+dens(l)
                  fullekin=fullekin+dens(l)*
     $                 (vxgrid(l)**2+vygrid(l)**2+vzgrid(l)**2)
                  fulleint=fulleint+press(l)/(gam-1)
                  fullepot=fullepot-dens(l)/sqrt(
     $                 (xcartgrid(i)-xbh)**2+
     $                 (ycartgrid(j)-ybh)**2+
     $                 (zcartgrid(k)-zbh)**2)
               endif

               if(dens(l).lt.1.d-12) then
                  dens(l)=0.d0
                  ucartgrid(i,j,k)=0.d0
               else
c                  if(nintvar.eq.1) then
c                     ucartgrid(i,j,k)=press(l)/dens(l)**gam
c                  else
cc     p=(gam-1)*rho*u, so u=p/((gam-1)*rho)
                     ucartgrid(i,j,k)=press(l)/((gam-1)*dens(l))
c                  endif

                  if(ucartgrid(i,j,k).gt.amax) then
                     amax=ucartgrid(i,j,k)
                     if(myrank.eq.0)write(69,*) 'amax tracker:',
     $                    l,press(l),dens(l),amax
                  endif

                  if(ucartgrid(i,j,k).gt.amax) then
                     write(69,*) 'trouble',i,j,k,press(l),dens(l),
     $                    ucartgrid(i,j,k),l,press(l)/dens(l)**gam
                     stop
                  endif
               endif

               presscartgrid(i,j,k)=press(l)
               denscartgrid(i,j,k)=dens(l)
               dpcartgrid(i,j,k)=dens(l)*press(l)
               vxcartgrid(i,j,k)=vxgrid(l)
               vycartgrid(i,j,k)=vygrid(l)
               vzcartgrid(i,j,k)=vzgrid(l)
            enddo
         enddo
      enddo
      fullamtot=2*fullamtot*dx*dy*dz
      fullekin=fullekin*dx*dy*dz
      fulleint=2*fulleint*dx*dy*dz
      fullepot=ambh*2*fullepot*dx*dy*dz
      if(myrank.eq.0) then
         write(69,*) 'after throwing away ultra-low densities, amax=',
     $        amax
      endif

      print *, 'will try for n=',n

      cellvolume=2*(xf-xi)*(yf-yi)*(zf-zi)/n
      anumden = 1.d0/cellvolume
      hc=(3.d0/32.d0/3.1415926535897932384626d0*
     $              1.9d0*nnopt/anumden)**(1.d0/3.d0)
      if(myrank.eq.0) write(69,*) 'hc=',hc
 
      a1=(cellvolume/2.d0**0.5d0)**(1.d0/3.d0)
      ixmax=int((xf-xi)/a1)-1
      iymax=int((yf-yi)/(3.d0**0.5d0/2.d0*a1))-1
      izmax=int((zf-zi)/(0.5d0*(8.d0/3.d0)**0.5d0*a1))-1

      if(ixmax.gt.hcpmax .or. iymax.gt.hcpmax .or.izmax.gt.hcpmax)then
         print *,'ixmax,iymax,izmax=',ixmax,iymax,izmax
         stop
      endif
      
      if(myrank.eq.0) then
         write(69,*)'best estimate of total mass=',fullamtot
         write(69,*)'best estimate of total kinetic energy',fullekin
         write(69,*)'best estimate of total internal energy',fulleint
         write(69,*)'best estimate of total potential energy',fullepot

         write(69,*)'cellvolume=',cellvolume
         write(69,*)'a1=',a1
         write(69,*) 'ixmax, iymax, izmax=',ixmax,iymax,izmax
         write(69,*) 'min x=',xi
         write(69,*) 'max x=',xi + ixmax*a1 + 0.5d0*a1
         write(69,*) 'min y=',yi 
         write(69,*) 'max y=',yi + iymax*3.d0**0.5d0/2.d0*a1
     $        +1.d0/3.d0**0.5d0*a1
         write(69,*) 'min z=',zi
         write(69,*) 'max z=',zi + izmax*0.5d0*(8.d0/3.d0)**0.5d0*a1
      endif

      do iz=0,izmax
         deltaz=iz*0.5d0*(8.d0/3.d0)**0.5d0*a1
         dzhcp(iz)=deltaz
         k=int(deltaz/dz)

         do iy=0,iymax
            deltay=iy*3.d0**0.5d0/2.d0*a1
     $           -(mod(abs(iz),2)-1d0)*1.d0/3.d0**0.5d0*a1
            dyhcp(iy,iz)=deltay
            j=int(deltay/dy)
            do ix=0,ixmax
               deltax=ix*a1+mod(abs(iy),2)*0.5d0*a1
               dxhcp(ix,iy)=deltax
               i=int(deltax/dx)

c     for now, use tri-linear interpolation among the eight cartesian
c     grid point (i,j,k), (i+1,j,k), (i,j+1,k), (i+1,j+1,k), (i,j,k+1),
c     (i+1,j,k+1), (i,j+1,k+1), and (i+1,j+1,k+1).
               f000=(deltax-i*dx)*(deltay-j*dy)*(deltaz-k*dz)
               f100=((i+1)*dx-deltax)*(deltay-j*dy)*(deltaz-k*dz)
               f010=(deltax-i*dx)*((j+1)*dy-deltay)*(deltaz-k*dz)
               f110=((i+1)*dx-deltax)*((j+1)*dy-deltay)*(deltaz-k*dz)
               f001=(deltax-i*dx)*(deltay-j*dy)*((k+1)*dz-deltaz)
               f101=((i+1)*dx-deltax)*(deltay-j*dy)*((k+1)*dz-deltaz)
               f011=(deltax-i*dx)*((j+1)*dy-deltay)*((k+1)*dz-deltaz)
               f111=((i+1)*dx-deltax)*((j+1)*dy-deltay)*
     $              ((k+1)*dz-deltaz)
               
c     the following array will help determine what hcp lattice points
c     will and will not have a particle
               dphcpgrid(ix,iy,iz)=(f000*dpcartgrid(i+1,j+1,k+1)
     $              + f100*dpcartgrid(i,  j+1,k+1)
     $              + f010*dpcartgrid(i+1,j,  k+1)
     $              + f110*dpcartgrid(i,  j,  k+1)
     $              + f001*dpcartgrid(i+1,j+1,k)
     $              + f101*dpcartgrid(i,  j+1,k)
     $              + f011*dpcartgrid(i+1,j,  k)
     $              + f111*dpcartgrid(i,  j,  k))/(dx*dy*dz)

            enddo
         enddo
      enddo


      p = 0
      l = 0
      tooclose=3*hc ! particles that are < tooclose from vacuum will not be kept
                                ! note: hc is less than the h near the surface, and
                                ! 3*hc is approximately the 2*h_surface
      do iz=0,izmax
         deltaz=iz*0.5d0*(8.d0/3.d0)**0.5d0*a1
         ztry= zi + deltaz
         k=int(deltaz/dz)
         dk=int(tooclose/(0.5d0*(8.d0/3.d0)**0.5d0*a1))
         kmin=max(iz-dk,0)
         kmax=min(iz+dk,izmax)
         do iy=0,iymax
            deltay=iy*3.d0**0.5d0/2.d0*a1
     $           -(mod(abs(iz),2)-1d0)*1.d0/3.d0**0.5d0*a1
            ytry= yi + deltay
            j=int(deltay/dy)
            dj=int(tooclose/(3.d0**0.5d0/2.d0*a1))+1
            jmin=max(iy-dj,0)
            jmax=min(iy+dj,iymax)
            do ix=0,ixmax
               l=l+1

               deltax=ix*a1+mod(abs(iy),2)*0.5d0*a1
               xtry= xi + deltax
               i=int(deltax/dx)
               di=int(tooclose/a1)+1
               imin=max(ix-di,0)
               imax=min(ix+di,ixmax)

               if(i+1.ge.nx .or. j+1.ge.ny .or. k+1.ge.nz) then
                  print *, 'off the cartgrid',i,j,k
                  stop
               endif

c     reject particle locations that will have smoothing kernels extending into
c     vacuum:
               do kn=kmin,kmax
                  do jn=jmin,jmax
                     do in=imin,imax
                        if(  (dxhcp(ix,iy)-dxhcp(in,jn))**2+
     $                       (dyhcp(iy,iz)-dyhcp(jn,kn))**2+
     $                       (dzhcp(iz)-dzhcp(kn))**2
     $                       .le.tooclose**2
     $                       .and.
     $                       dphcpgrid(in,jn,kn)
c     $                       .eq.0.d0) goto 42
     $                       .le.1d-28) goto 42
                     enddo
                  enddo
               enddo

c     put a particle at this location:
               p = p+1
               x(p) = xtry
               y(p) = ytry
               z(p) = ztry

c     for now, use tri-linear interpolation among the eight cartesian
c     grid point (i,j,k), (i+1,j,k), (i,j+1,k), (i+1,j+1,k), (i,j,k+1),
c     (i+1,j,k+1), (i,j+1,k+1), and (i+1,j+1,k+1).
               f000=(deltax-i*dx)*(deltay-j*dy)*(deltaz-k*dz)
               f100=((i+1)*dx-deltax)*(deltay-j*dy)*(deltaz-k*dz)
               f010=(deltax-i*dx)*((j+1)*dy-deltay)*(deltaz-k*dz)
               f110=((i+1)*dx-deltax)*((j+1)*dy-deltay)*(deltaz-k*dz)
               f001=(deltax-i*dx)*(deltay-j*dy)*((k+1)*dz-deltaz)
               f101=((i+1)*dx-deltax)*(deltay-j*dy)*((k+1)*dz-deltaz)
               f011=(deltax-i*dx)*((j+1)*dy-deltay)*((k+1)*dz-deltaz)
               f111=((i+1)*dx-deltax)*((j+1)*dy-deltay)*
     $              ((k+1)*dz-deltaz)

               vx(p)=(f000*vxcartgrid(i+1,j+1,k+1)
     $              + f100*vxcartgrid(i,  j+1,k+1)
     $              + f010*vxcartgrid(i+1,j,  k+1)
     $              + f110*vxcartgrid(i,  j,  k+1)
     $              + f001*vxcartgrid(i+1,j+1,k)
     $              + f101*vxcartgrid(i,  j+1,k)
     $              + f011*vxcartgrid(i+1,j,  k)
     $              + f111*vxcartgrid(i,  j,  k))/(dx*dy*dz)

               vy(p)=(f000*vycartgrid(i+1,j+1,k+1)
     $              + f100*vycartgrid(i,  j+1,k+1)
     $              + f010*vycartgrid(i+1,j,  k+1)
     $              + f110*vycartgrid(i,  j,  k+1)
     $              + f001*vycartgrid(i+1,j+1,k)
     $              + f101*vycartgrid(i,  j+1,k)
     $              + f011*vycartgrid(i+1,j,  k)
     $              + f111*vycartgrid(i,  j,  k))/(dx*dy*dz)

               vz(p)=(f000*vzcartgrid(i+1,j+1,k+1)
     $              + f100*vzcartgrid(i,  j+1,k+1)
     $              + f010*vzcartgrid(i+1,j,  k+1)
     $              + f110*vzcartgrid(i,  j,  k+1)
     $              + f001*vzcartgrid(i+1,j+1,k)
     $              + f101*vzcartgrid(i,  j+1,k)
     $              + f011*vzcartgrid(i+1,j,  k)
     $              + f111*vzcartgrid(i,  j,  k))/(dx*dy*dz)

               am(p)= f000*denscartgrid(i+1,j+1,k+1)
     $              + f100*denscartgrid(i,  j+1,k+1)
     $              + f010*denscartgrid(i+1,j,  k+1)
     $              + f110*denscartgrid(i,  j,  k+1)
     $              + f001*denscartgrid(i+1,j+1,k)
     $              + f101*denscartgrid(i,  j+1,k)
     $              + f011*denscartgrid(i+1,j,  k)
     $              + f111*denscartgrid(i,  j,  k)

               u(p) =(f000*ucartgrid(i+1,j+1,k+1)
     $              + f100*ucartgrid(i,  j+1,k+1)
     $              + f010*ucartgrid(i+1,j,  k+1)
     $              + f110*ucartgrid(i,  j,  k+1)
     $              + f001*ucartgrid(i+1,j+1,k)
     $              + f101*ucartgrid(i,  j+1,k)
     $              + f011*ucartgrid(i+1,j,  k)
     $              + f111*ucartgrid(i,  j,  k))/(dx*dy*dz)

               if(u(p).gt.amax) then
c               if(p.eq.298580 .or. p+1.eq.298580) then
                  write(69,*) 'p=',p,'u(p)=',u(p)
                  write(69,*) i,j,k
                  write(69,*) ix,iy,iz
                  write(69,*) vx(p),vy(p),vz(p)
                  write(69,*) dphcpgrid(ix,iy,iz)
                  write(69,*)'last number is dhcpgrid'
                  write(69,*)'u:'
                  write(69,*)ucartgrid(i+1,j+1,k+1),
     $              ucartgrid(i,  j+1,k+1),
     $              ucartgrid(i+1,j,  k+1),
     $              ucartgrid(i,  j,  k+1),
     $              ucartgrid(i+1,j+1,k),
     $              ucartgrid(i,  j+1,k),
     $              ucartgrid(i+1,j,  k),
     $              ucartgrid(i,  j,  k)
                  write(69,*)'dens:'
                  write(69,*)denscartgrid(i+1,j+1,k+1),
     $              denscartgrid(i,  j+1,k+1),
     $              denscartgrid(i+1,j,  k+1),
     $              denscartgrid(i,  j,  k+1),
     $              denscartgrid(i+1,j+1,k),
     $              denscartgrid(i,  j+1,k),
     $              denscartgrid(i+1,j,  k),
     $              denscartgrid(i,  j,  k)
                  write(69,*)'press:'
                  write(69,*)presscartgrid(i+1,j+1,k+1),
     $              presscartgrid(i,  j+1,k+1),
     $              presscartgrid(i+1,j,  k+1),
     $              presscartgrid(i,  j,  k+1),
     $              presscartgrid(i+1,j+1,k),
     $              presscartgrid(i,  j+1,k),
     $              presscartgrid(i+1,j,  k),
     $              presscartgrid(i,  j,  k)
                  write(69,*)'dp:'
                  write(69,*)dpcartgrid(i+1,j+1,k+1),
     $              dpcartgrid(i,  j+1,k+1),
     $              dpcartgrid(i+1,j,  k+1),
     $              dpcartgrid(i,  j,  k+1),
     $              dpcartgrid(i+1,j+1,k),
     $              dpcartgrid(i,  j+1,k),
     $              dpcartgrid(i+1,j,  k),
     $              dpcartgrid(i,  j,  k)
                  stop
               endif

               if(nintvar.eq.1) then
c     at this point, u(p) is a specific internal energy,
c     but because nintvar==1, we want u(p) to be changed to a instead.
cc     p=(gam-1)*rho*u=a*rho^gam, so a=(gam-1)*u*rho^(1-gam)
                  u(p)=u(p)*(gam-1)*(am(p)/(dx*dy*dz))**(1-gam)
               endif


               hp(p) = hc

c               aa(p)=hmax
c               bb(p)=(1/hp(p)-1/hmax)/am(p)**(1.d0/3.d0)
               
               if(iz.gt.0) then
                  p=p+1
                  x(p) = x(p-1)
                  y(p) = y(p-1)
                  z(p) = -z(p-1)
                  u(p) = u(p-1)
                  am(p) = am(p-1)
                  hp(p) = hp(p-1)
                  vx(p) = vx(p-1)
                  vy(p) = vy(p-1)
                  vz(p) = -vz(p-1)
c                  aa(p) = aa(p-1)
c                  bb(p) = bb(p-1)
               endif

 42            continue

c               write(69,*) ix,iy,iz,xtry,ytry,ztry
            enddo
         enddo
      enddo


c      ntot=p
c      n=ntot
      n=p

      if(.false.) then
         do i=1,n
            xharmonic=x(i)-xbh
            yharmonic=y(i)-ybh
            zharmonic=z(i)-zbh
            rharmonic=(xharmonic**2+yharmonic**2+zharmonic**2)**0.5d0
            ratio=1+ambh/rharmonic
            x(i)=xharmonic*ratio+xbh
            y(i)=yharmonic*ratio+ybh
            z(i)=zharmonic*ratio+zbh
            rdotv=xharmonic*vx(i)+yharmonic*vy(i)+zharmonic*vz(i)
            vx(i)=vx(i)*ratio-xharmonic*ambh*rdotv/rharmonic**3
            vy(i)=vy(i)*ratio-yharmonic*ambh*rdotv/rharmonic**3
            vz(i)=vz(i)*ratio-zharmonic*ambh*rdotv/rharmonic**3
         enddo
         if(myrank.eq.0) write(69,*)
     $        '*done* with harmonic coordinate transformation'         
      else
         if(myrank.eq.0) write(69,*)
     $        '*not* doing harmonic coordinate transformation'         
      endif

      if(myrank.eq.0)
     $     write(69,*)'number of particles=',n

      n=n+1
      ntot=n
      i=ntot
      am(i)=ambh
      x(i)=xbh
      y(i)=ybh
      z(i)=zbh
      vx(i)=0.d0
      vy(i)=0.d0
      vz(i)=0.d0
      gx(i)=0.d0
      gy(i)=0.d0
      gz(i)=0.d0
      vxdot(i)=0.d0
      vydot(i)=0.d0
      vzdot(i)=0.d0
      u(i)=0.d0
      udot(i)=0.d0
      if(hco.gt.0.d0) then
         hp(i)=hco
      else
         hp(i)=hp(ntot-1)
      endif
      cc(i)=0
      if(myrank.eq.0) write(69,*)'black_hole_mass',am(i)
      if(myrank.eq.0) write(69,*)'black_hole_smoothing_length',hp(i)

c     set black hole to origin and stationary initially
c      xbh=-0.3266d0      !-0.43d0
c      ybh=-0.7796d0      !-0.73d0
c      zbh=0.d0
c      ambh=11.0d0        !8.8d0__wrong mass!!!!
c      xbh=xposbh
c      ybh=yposbh
c      zbh=zposbh
c      ambh=bhmass
c      vxbh=0.d0
c      vybh=0.d0
c      vzbh=0.d0
c      print*,'bh:',xbh,ybh,zbh,ambh
c      print*,'grid: ',xi,xf,yi,yf,zi,zf
 
c      call gravquant
c      if(myrank.eq.0) write(69,*) 'back from gravquant'

      t=0.d0


c      do i=1,n
c         nn(i)=nnopt
c      enddo
c      call rho_and_h
c
c      if(myrank.eq.0) write(69,*) 'back from initial rho_and_h'

c      call linkedlists
c      do i=1,n
c         call newnene(i)
c      enddo
c      do j=1,1
c      do j=1,6
cc         print *,myrank,'j=',j,myrank
c         nnemin=10000
c         nnemax=0
c         nneavr=0
c         nnesig=0
c         do i=1,n                                
c            nnemin=min(nnemin,nn(i))
c            if(nnemin.eq.0) then
c               print *,myrank,'nn=0',i,myrank
c               stop
c            endif
c            nnemax=max(nnemax,nn(i))
c            nneavr=nneavr+nn(i)
c            nnesig=nnesig+nn(i)**2
c         enddo                                  
c         nneavr=int(float(nneavr)/float(n))
c         nnesig=int(sqrt(float(nnesig)/float(n)-float(nneavr)**2))
cc         if(myrank.eq.0)
cc     $        write (6,*)'nnmin:',nnemin,' nnmax:',nnemax
cc         if(myrank.eq.0)
cc     $        write(6,*)' avg:',nneavr,' sig:',nnesig
c         
cc         if(j.gt.1 .and. nnemax-nnemin.le.0.2d0*nneavr) goto 422
c
cc         call adjust            !note: this includes calls to newnene(i) at end of routine
c         
c         call rho_and_h
c         hpmin=1.d30                             
c         hpmax=0.d0
c         do i=1,n
c            hpmin=min(hpmin,hp(i))                
c            hpmax=max(hpmax,hp(i))
c         enddo
c         if(myrank.eq.0)write (69,*) 'opthp2: hpmin=',hpmin       
c         if(myrank.eq.0)write (69,*) '       hpmax=',hpmax       
c
c      enddo
c
c 422  continue

c      if(nintvar.eq.2) then
cc     the u(i) array will actually be specific internal energy u
cc     p=(gam-1)*rho*u=a*rho^gam, so u=a*rho^(gam-1)/(gam-1)
c         do i=1,ntot-1
cc            write(80+myrank,*)i,u(i),rho(i),hp(i)
c            u(i)=u(i)*(am(i)/(dx*dy*dz))**(gam-1)/(gam-1)
c
cc            write(80+myrank,*)i,u(i)
cc            stop
c
c         enddo
c      endif

      amtot=0.d0
      do i=1,ntot-1
         amtot=amtot+am(i)
      enddo
c      write(69,*)'before rescaling total mass in gas=',amtot
      do i=1,ntot-1
         am(i)=am(i)/amtot*fullamtot
      enddo
      amtot=0.d0
      do i=1,ntot-1
         amtot=amtot+am(i)
      enddo
      if(myrank.eq.0)then
         write(69,*) 'total rest mass in gas=',amtot
         i=1
         write(69,*)x(i),y(i),z(i),dens(i),press(i),i
         i=ntot-1
         write(69,*)x(i),y(i),z(i),dens(i),press(i),i
      endif

      if(myrank.eq.0) write(69,*) 'about to call lfstart'      

c     prepare leap-frog scheme for first iteration:
      call lfstart

      if(myrank.eq.0) write(69,*) 'back from lfstart'


      return

c     error condition:
 100  stop 'init:  error reading input file ???'

      end
      
