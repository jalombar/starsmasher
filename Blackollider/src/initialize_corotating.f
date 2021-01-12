      subroutine corotating                                         
********************************************************
c     release 1.0
c     sets up a synchronized binary with equal-mass ns,
c     using previously calculated relaxed single-star models,
c     processed by b2pos into a file called "pos.sph".
c     called by init
c     calls gravquant,lfstart
**********************************************************
      include 'starsmasher.h'
      real*8 xcm1,ycm1,zcm1,xcm2,ycm2,zcm2,am1,am2
      common/centersofmass/xcm1,ycm1,zcm1,xcm2,ycm2,zcm2,am1,am2
      integer n2,i,nchk,corepts
      integer nnoptold,noutold,nitold,navold,ngrold,nrelaxold
      real*8 hcoold,hfloorold,sep0old,tfold,dtoutold,told,
     $     alphaold,betaold,trelaxold,dtold
      integer istart,corepts1,corepts2
      common/core/corepts1,corepts2
      logical resetsep0,twofiles
      real*8 am1chk,am2chk,realdummy1,realdummy2,realdummy3
      real*8 deltax1,deltay1,deltaz1,deltax2,deltay2,deltaz2
      real*8 egsol,solrad
      parameter(egsol=1.9891d+33,solrad=6.9599d10)
      real*8 deltaxbin,deltaybin,deltazbin
      character*7 dummy

      inquire(file='input.bs',exist=resetsep0)
      if(resetsep0)then
         write(69,*) 'will read file input.bs'
         open(30,file='input.bs')
         read(30,*) dummy
         if(dummy.ne.'binary:')then
            write(69,*) 'input.bs should have binary first'
            stop
         endif
         read(30,*) am1chk
         read(30,*) deltax1,deltay1,deltaz1
         read(30,*) realdummy1,realdummy2,realdummy3
         read(30,*) am2chk
         read(30,*) deltax2,deltay2,deltaz2
         close(30)

         deltax1=deltax1*solrad/runit
         deltay1=deltay1*solrad/runit
         deltaz1=deltaz1*solrad/runit
         deltax2=deltax2*solrad/runit
         deltay2=deltay2*solrad/runit
         deltaz2=deltaz2*solrad/runit

         if(myrank.eq.0) then
            write(69,*)'star 1 mass [msun]=',am1chk
            write(69,*)'star 2 mass [msun]=',am2chk
         endif

         deltaxbin=deltax2-deltax1
         deltaybin=deltay2-deltay1
         deltazbin=deltaz2-deltaz1

         if((deltaxbin**2+deltaybin**2+deltazbin**2)**0.5d0.gt.
     $        am1chk**0.85d0+am2chk**0.85d0) then
            if(myrank.eq.0) then
               write(69,*)'n.b.: using input.bs to set sep0 ...'
               write(69,*)' ****(ignoring value in sph.input) ****'
               write(69,*)'disabling any scanning'
            endif
            sep0=(deltaxbin**2+deltaybin**2+deltazbin**2)**0.5d0
            sepfinal=sep0
            tscanon=-1
         else
            if(myrank.eq.0)write(69,*)
     $           'n.b.: this seems to be a contact binary'

            if(sepfinal.lt.
     $           (deltaxbin**2+deltaybin**2+deltazbin**2)**0.5d0) then
               if(myrank.eq.0) then
                  write(69,*)'using input.bs to set sepfinal ...'
                  write(69,*)' ****(ignoring value in sph.input) ****'
               endif
               sepfinal=(deltaxbin**2+deltaybin**2+deltazbin**2)**0.5d0
               if(myrank.eq.0)write(69,*)'reassigning sepfinal=',sepfinal
            endif
            if(sep0.lt.am1chk+am2chk) then
               if(myrank.eq.0) then
                  write(69,*)'sep0 value in sph.input seems too small'
                  write(69,*)' ****(ignoring value in sph.input) ****'
               endif
               sep0=min(max(am1chk+am2chk,1.5d0*sepfinal),1.5d0*sepfinal)
               if(myrank.eq.0) write(69,*)'reassigning sep0=',sep0
            endif
            if(tscanon.lt.0.d0) then
               if(myrank.eq.0) write(69,*)'scanning should be enabled'
               tscanon=20.d0
               if(myrank.eq.0) write(69,*)'reassigning tscanon=',20.d0
            endif
         endif
      endif

      if(myrank.eq.0) write(69,*)'corotating: sep0=',sep0

      if(sep0.lt.sepfinal) then
         write(69,*)'did you mean to scan outward?'
         stop
      endif

      if(myrank.eq.0) write(69,*)'corotating: sep0/rsun=',sep0*runit/solrad
      corepts1=0
      corepts2=0
      if(myrank.eq.0) write (69,*) 'corotating: reading start files ...'
      
      open(12,file=startfile1,form='unformatted')
c     (the following read sequence must match exactly the write sequence
c     used in subroutine dump)
      read(12) n1,nnoptold,hcoold,hfloorold,sep0old,
     $     tfold,dtoutold,noutold,nitold,told,
     $     navold,alphaold,betaold,tjumpahead,
     $     ngrold,
     $     nrelaxold,trelaxold,dtold,omega2
      am1=0.d0
      do i=1,n1
         read (12) x(i),y(i),z(i),am(i),hp(i),rho(i),vx(i),vy(i),
     $        vz(i),vxdot(i),vydot(i),vzdot(i),u(i),udot(i),
     $        grpot(i),meanmolecular(i),
     $        cc(i), realdummy1
         ueq(i)=0
         tthermal(i)=1d30

         am1=am1+am(i)
         if(hp(i).le.0.d0 .or. u(i).eq.0.d0) then
            if(myrank.eq.0) write(69,*)'particle',i,
     $           'is a corepoint of mass',am(i)
         endif
      enddo
      read(12) nchk
      close(12)

      if(myrank.eq.0) write(69,*) 'mass1= ', am1

      if (nchk.ne.n1) stop 'corotating: problem with file sph.start1u'

      if(myrank.eq.0) write(69,*)'n1=',n1

      inquire(file=startfile2,exist=twofiles)
      if(twofiles) then
         open(12,file=startfile2,form='unformatted')
c     (the following read sequence must match exactly the write sequence
c     used in subroutine dump)
         read(12) n2,nnoptold,hcoold,hfloorold,sep0old,
     $        tfold,dtoutold,noutold,nitold,told,navold,
     $        alphaold,betaold,tjumpahead,ngrold,
     $        nrelaxold,trelaxold,dtold
         if(myrank.eq.0) write(69,*)'n2=',n2
         am2=0.d0
         istart=n1+1-corepts1
         do i=istart,istart+n2-1
            read (12) x(i),y(i),z(i),am(i),hp(i),rho(i),vx(i),vy(i),
     $           vz(i),vxdot(i),vydot(i),vzdot(i),u(i),udot(i),
c     $           gx(i),gy(i),gz(i),grpot(i),meanmolecular(i),
c     $           cc(i)
     $           grpot(i),meanmolecular(i),
     $           cc(i)
            ueq(i)=0
            tthermal(i)=1d30

            x(i)=-x(i)
            vx(i)=-vx(i)
            vxdot(i)=-vxdot(i)
            gx(i)=-gx(i)
            y(i)=-y(i)
            vy(i)=-vy(i)
            vydot(i)=-vydot(i)
            gy(i)=-gy(i)

            am2=am2+am(i)
            if(hp(i).le.0.d0 .or. u(i).eq.0.d0) then
               if(myrank.eq.0) write(69,*)'particle',i,
     $              'is a corepoint of mass',am(i)
            endif
         enddo
         read(12) nchk
         close(12)
         if (nchk.ne.n2) then
            write(69,*) 'corotating: problem with file sph.start2u'
            stop
         endif
      else
         if(myrank.eq.0) write(69,*)'secondary point mass mass',mbh
         i=n1+1
         n2=1
         am2=mbh
         am(i)=am2
         x(i)=0.d0
         y(i)=0.d0
         z(i)=0.d0
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
            hp(i)=hp(n1)
         endif
         aa(i)=hp(i)
         bb(i)=0
         cc(i)=0
         dd(i)=0
         if(myrank.eq.0) write(69,*)'secondary pt mass hp',hp(i)
      endif
      corepts=corepts1+corepts2
      ntot=n1+n2
      if (ntot.gt.nmax) then
         write(69,*) 'must increase nmax...'
         stop
      endif
      n=ntot-corepts

      if(myrank.eq.0) write(69,*) 'mass2= ', am2

      if(resetsep0) then
         if(dabs(am1*munit/egsol/am1chk-1.d0).gt.1.d-8 .or.
     $        dabs(am2*munit/egsol/am2chk-1.d0).gt.1.d-8)then
            write(69,*)'mass(es) in sph.start?u does not match with'
            write(69,*)'mass(es) in input.bs'
            stop
         endif
      endif

      if(ntot.gt.nmax) then
         write(69,*)'error: n1+n2>nmax'
         stop
      endif

      if(myrank.eq.0) write (69,*)'corotating: n=',n,'ntot=',ntot
c     shift stars along x-axis so that cm is at the origin
c     and separation = sep0
      if(myrank.eq.0) write(69,*)'sep0=',sep0
      do i=1,n1-corepts1
         x(i)=x(i)-sep0*am2/(am1+am2)
      enddo
      do i=n1+1-corepts1,ntot-corepts1
         x(i)=x(i)+sep0*am1/(am1+am2)
      enddo
      do i=ntot-corepts1+1,ntot
         x(i)=x(i)-sep0*am2/(am1+am2)
      enddo
      if(myrank.eq.0) then
         write(69,*)'star 1 has mass',am1,n1-corepts1,ntot-corepts1+1,ntot
         write(69,*)'star 2 has mass',am2,n1+1-corepts1,ntot-corepts1
         write(69,*)'star 1 shifted by',-sep0*am2/(am1+am2)
         write(69,*)'star 2 shifted by',sep0*am1/(am1+am2)
      endif

      call stride_setup

c     prepare leap-frog scheme for first iteration:
      call lfstart
      if(myrank.eq.0) write (69,*) 'corotating:          ... done'

      return
      end
