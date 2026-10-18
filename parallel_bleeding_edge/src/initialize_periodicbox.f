      subroutine periodicbox
c     iname='pbx': start from particles in a periodic box, read from
c     startfile1 in the out*.sph format.  Positions, masses, smoothing
c     lengths, velocities and u are taken as they are, so the file can be a
c     lattice laid down by a script or a dump from an earlier run (for
c     example two relaxed halves of a shock tube put side by side).
      include 'starsmasher.h'
      real*8 divv(nmax)
      common/commdivv/divv
      integer i,nchk
      integer nnoptold,noutold,nitold,navold,ngrold,nrelaxold
      real*8 hcoold,hfloorold,sep0old,tfold,dtoutold,told,
     $     alphaold,betaold,trelaxold,dtold,tjumpaheadold,omega2old

      if(.not.periodic) stop 'periodicbox: set boxx, boxy or boxz'

      if(myrank.eq.0) write (69,*) 'periodicbox: reading startfile1 ',
     $     trim(startfile1)

      open(12,file=startfile1,form='unformatted',status='old')
      read(12) n,nnoptold,hcoold,hfloorold,sep0old,
     $     tfold,dtoutold,noutold,nitold,told,
     $     navold,alphaold,betaold,tjumpaheadold,
     $     ngrold,nrelaxold,trelaxold,dtold,omega2old
      if (n.gt.nmax) stop 'periodicbox: n>nmax ???'
      do i=1,n
         read (12) x(i),y(i),z(i),am(i),hp(i),rho(i),vx(i),vy(i),
     $        vz(i),vxdot(i),vydot(i),vzdot(i),u(i),udot(i),
     $        grpot(i),meanmolecular(i),
     $        cc(i),divv(i)
         if(u(i).eq.0.d0) stop 'periodicbox: no point particles'
         vxdot(i)=0
         vydot(i)=0
         vzdot(i)=0
         udot(i)=0
         ueq(i)=0
         tthermal(i)=1d30
      enddo
      read(12) nchk
      close(12)
      if (nchk.ne.n) stop 'periodicbox: problem with startfile1 ???'

      ntot=n
      if(myrank.eq.0) then
         write(69,*)'periodicbox: n=',n
         write(69,*)'periodicbox: boxx,boxy,boxz=',boxx,boxy,boxz
      endif

      call pbcpositions
      call stride_setup

c     prepare leap-frog scheme for first iteration:
      call lfstart

      return
      end
************************************************************************
      subroutine pbcwrap(dx,dy,dz)
c     replace a separation by that of the nearest periodic image
      include 'starsmasher.h'
      real*8 dx,dy,dz
      if(boxx.gt.0.d0) dx=dx-boxx*anint(dx/boxx)
      if(boxy.gt.0.d0) dy=dy-boxy*anint(dy/boxy)
      if(boxz.gt.0.d0) dz=dz-boxz*anint(dz/boxz)
      return
      end
************************************************************************
      subroutine pbcpositions
c     put every particle back inside the box, -box/2 <= x < box/2
      include 'starsmasher.h'
      integer i
      do i=1,n
         if(boxx.gt.0.d0) x(i)=x(i)-boxx*floor(x(i)/boxx+0.5d0)
         if(boxy.gt.0.d0) y(i)=y(i)-boxy*floor(y(i)/boxy+0.5d0)
         if(boxz.gt.0.d0) z(i)=z(i)-boxz*floor(z(i)/boxz+0.5d0)
      enddo
      return
      end
************************************************************************
      subroutine pbcimages(pad,store,nimg,ximg,imap)
c     Periodic images of the particles that lie within pad of a face, for
c     the neighbour tree.  Image k sits at ximg(:,k) and is a copy of
c     particle imap(k).  pad(d) is at most half the period in direction d,
c     so each image is the only one of its particle within pad of the box.
c     With store=.false. the images are only counted, so that the caller
c     can allocate ximg and imap to size nimg.
      include 'starsmasher.h'
      real*8 pad(3),ximg(3,*),box(3),xi(3),xs(3)
      integer nimg,imap(*),i,sx,sy,sz,lo(3),hi(3),s(3),d
      logical keep,store

      box(1)=boxx
      box(2)=boxy
      box(3)=boxz
      do d=1,3
         lo(d)=0
         hi(d)=0
         if(box(d).gt.0.d0) then
            lo(d)=-1
            hi(d)=1
         endif
      enddo
      nimg=0
      do i=1,n
         xi(1)=x(i)
         xi(2)=y(i)
         xi(3)=z(i)
         do sx=lo(1),hi(1)
         do sy=lo(2),hi(2)
         do sz=lo(3),hi(3)
            if(sx.eq.0 .and. sy.eq.0 .and. sz.eq.0) cycle
            s(1)=sx
            s(2)=sy
            s(3)=sz
            keep=.true.
            do d=1,3
               xs(d)=xi(d)+s(d)*box(d)
               if(s(d).ne.0 .and. abs(xs(d)).ge.0.5d0*box(d)+pad(d))
     $              keep=.false.
            enddo
            if(keep) then
               nimg=nimg+1
               if(.not.store) cycle
               ximg(1,nimg)=xs(1)
               ximg(2,nimg)=xs(2)
               ximg(3,nimg)=xs(3)
               imap(nimg)=i
            endif
         enddo
         enddo
         enddo
      enddo
      return
      end
************************************************************************
      subroutine pbcmap(lst,cnt,imap,r2,padmin)
c     Turn the tree indices of a neighbour list into particle indices, an
c     image k (tree index n+k) becoming its particle imap(k).  A search
c     reaching past the padding could have missed images, so stop then.
      include 'starsmasher.h'
      integer cnt,lst(cnt),imap(*),k
      real*8 r2,padmin
      if(r2.gt.padmin**2) then
         write(6,*)'pbcmap: search radius',sqrt(r2),'exceeds padding',
     $        padmin,'; 2h must stay below half of each period'
         stop 'pbcmap: search radius exceeds periodic padding'
      endif
      do k=1,cnt
         if(lst(k).gt.n) lst(k)=imap(lst(k)-n)
      enddo
      return
      end
************************************************************************
      real*8 function pbcdist2(i,j)
c     squared separation of particles i and j via the nearest image
      include 'starsmasher.h'
      integer i,j
      real*8 dx,dy,dz
      dx=x(i)-x(j)
      dy=y(i)-y(j)
      dz=z(i)-z(j)
      call pbcwrap(dx,dy,dz)
      pbcdist2=dx*dx+dy*dy+dz*dz
      return
      end
