      subroutine gravforce
      include 'starsmasher.h'
      include 'mpif.h'
      integer i, ierr
      real*8 range(ntot)
      real*8 amasstot,fcmx,fcmy,fcmz
      integer mylength, mygravlength
      integer irank      
      integer nlist,ilist(ntot)

      if(nusegpus.eq.0)return

      if(myrank.lt.ngravprocs) then
c     each processor will compute gravity for its chunk of particles
         call cpu_time(time0)
         do i=1,ntot
            range(i) = 4*hp(i)**2;
            if(u(i).eq.0.d0) range(i)=-range(i)
         enddo
         
         if(myrank.eq.ngravprocs-1) call cpu_time(time1) 
ccc   ********* build the tree **********
         
         mygravlength=ngrav_upper-ngrav_lower+1

c     block steps: when any particle is inactive, the gpu computes gravity
c     on the active particles only, and the others' gx,gy,gz,grpot come
c     back as zero (they keep the accelerations of the start of their
c     step).  Each particle's force is the same sum over all particles
c     either way, so unlike the symmetric pair loop on cpus (see balAV3.f)
c     there is no fraction active above which doing them all is cheaper
         if(nblock.eq.1 .and. count(actblk(1:ntot)).lt.ntot) then
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
            call cpu_time(time2)         
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
      
c      call cpu_time(time3)
      
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
