      block data gravfar_data
c     the counters of the far passes that were kept and redone (gravforce), which are not given a value anywhere else
      logical farstarted
      integer nfallback,nfar
      common/gravfari/farstarted,nfallback,nfar
      data farstarted,nfallback,nfar /.false.,0,0/
      end

      subroutine gravfar_start
c     ngravoverlap=1 (direct gravity on the gpus): start the part of the gpu's work that does not depend on the smoothing lengths now,
c     before they are solved, with a limit for each h (1.25 times the old one).  gravforce finishes it, or does it all again if an h
c     has grown past its limit.
      include 'starsmasher.h'
      integer i,mygravlength,ok,nfallback,nfar
      real*8 rbound(nmax)
      logical farstarted
      common/gravfarr/rbound
      common/gravfari/farstarted,nfallback,nfar
      farstarted=.false.
      if(ngravoverlap.eq.0 .or. ngravtree.ne.0 .or. nusegpus.eq.0) return
      if(ngr.eq.0 .or. nblock.eq.1) return
      if(myrank.ge.ngravprocs) return
      mygravlength=ngrav_upper-ngrav_lower+1
!$omp parallel do simd schedule(static) private(i)
      do i=1,ntot
         rbound(i) = 1.5625d0*4*hp(i)**2
         if(u(i).eq.0.d0) rbound(i)=-rbound(i)
      enddo
      call firsthalf_grav_forces_far(ntot, ngrav_lower, mygravlength,
     $     x, y, z, am, rbound, q, nkernel, ok)
      farstarted=(ok.eq.1)
      return
      end

      subroutine gravforce
      include 'starsmasher.h'
      include 'mpif.h'
      integer i, ierr
      real*8 range(ntot)
      real*8 amasstot,fcmx,fcmy,fcmz
      integer mylength, mygravlength
      integer irank      
      integer nlist,ilist(ntot)
      integer nfallback,nfar
      real*8 rbound(nmax)
      logical farstarted,donegrav
      integer nbad
      common/gravfarr/rbound
      common/gravfari/farstarted,nfallback,nfar

      if(nusegpus.eq.0)return
      if(ngravtree.ne.0)return ! the tree finds the gravity when it is collected (grav_collect)

      if(myrank.lt.ngravprocs) then
c     each processor will compute gravity for its chunk of particles
         call wall_time(time0)
!$omp parallel do simd schedule(static)
         do i=1,ntot
            range(i) = 4*hp(i)**2;
            if(u(i).eq.0.d0) range(i)=-range(i)
         enddo
         
         if(myrank.eq.ngravprocs-1) call wall_time(time1) 
ccc   ********* build the tree **********
         
         mygravlength=ngrav_upper-ngrav_lower+1

c     the far pass of ngravoverlap=1 is on the gpu already; it is good if no h is above its limit, and if one is the whole of the
c     gravity is done again below
         donegrav=.false.
         if(farstarted) then
            farstarted=.false.
            nbad=0
!$omp parallel do simd schedule(static) private(i) reduction(+:nbad)
            do i=1,ntot
               if(abs(range(i)).gt.abs(rbound(i))) nbad=nbad+1
            enddo
            if(nbad.eq.0) then
               call firsthalf_grav_forces_near(ntot, x, y, z, am, range,
     $              nkernel)
               donegrav=.true.
               nfar=nfar+1
            else
               call gravity_far_abort
               nfallback=nfallback+1
               write(6,'(a,i8,a,i8)') 'gravity far pass redone: ',
     $              nfallback,' times in ',nfar+nfallback
            endif
         endif

c     block steps: when any particle is inactive, the gpu computes gravity
c     on the active particles only, and the others' gx,gy,gz,grpot come
c     back as zero (they keep the accelerations of the start of their
c     step).  Each particle's force is the same sum over all particles
c     either way, so unlike the symmetric pair loop on cpus (see balAV3.f)
c     there is no fraction active above which doing them all is cheaper
         if(donegrav) then
c     nothing more to launch
         else if(nblock.eq.1 .and. count(actblk(1:ntot)).lt.ntot) then
            nlist=0
            do i=ngrav_lower,ngrav_upper
               if(actblk(i)) then
                  nlist=nlist+1
                  ilist(nlist)=i
               endif
            enddo
            call firsthalf_grav_forces_list(ntot, ngrav_lower,
     $           mygravlength, nlist, ilist, x, y, z, am, range, nkernel)
         else
            call firsthalf_grav_forces(ntot, ngrav_lower, mygravlength,
     $           x, y, z, am, range,q,nkernel)
         endif
         if(myrank.eq.ngravprocs-1) then
            call wall_time(time2)         
            write (6,'( a, f6.3, a,i4, a, i4)')'1sthalf:',time2-time1
         endif
      endif

c      amasstot=0.d0
c      fcmx=0.d0
c      fcmy=0.d0
c      fcmz=0.d0
c      do i=1,ntot
c         fcmx=fcmx-mass(i)*gx(i)
c         fcmy=fcmy-mass(i)*gy(i)
c         fcmz=fcmz-mass(i)*gz(i)
c         amasstot=amasstot-mass(i)
c      enddo
c      fcmx=fcmx/amasstot
c      fcmy=fcmy/amasstot
c      fcmz=fcmz/amasstot
      
c      write(69,'(1x,a,3g12.4)')
c     $     'grape gave c.o.m. acceleration (to be subtracted) =',
c     $     fcmx,fcmy,fcmz
      
c      do i=1,ntot
c         gx(i)=gx(i)-fcmx
c         gy(i)=gy(i)-fcmy
c         gz(i)=gz(i)-fcmz
c      enddo
      
c      call wall_time(time3)
      
c      write (6,'(3(a,g10.3),g10.3,a)')
c     $     'gravfo: ',time3-time0,' s' !,
cc     $     time2-time1,'+',time1-time0,time3-time2,'overhead'
      
cc     despite the name force(j,i) is actually the *acceleration* of particle i
cc     in the j direction...
c      amasstot=0.d0
c      fcmx=0.d0
c      fcmy=0.d0
c      fcmz=0.d0
c      do i=1,ntot
c         fcmx=fcmx-mass(i)*gx(i)
c         fcmy=fcmy-mass(i)*gy(i)
c         fcmz=fcmz-mass(i)*gz(i)
c         amasstot=amasstot-mass(i)
c      enddo
c      fcmx=fcmx/amasstot
c      fcmy=fcmy/amasstot
c      fcmz=fcmz/amasstot
c      
c      write(69,'(1x,a,3g12.4)')
c     $     'gravforce gave c.o.m. acceleration=',
c     $     fcmx,fcmy,fcmz
      
      return
      end
