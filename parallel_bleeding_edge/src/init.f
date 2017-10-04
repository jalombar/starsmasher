      subroutine init
c     initialization of run (new or restart)
      include 'starsmasher.h'
      include 'mpif.h'
      real*8 xcm1,ycm1,zcm1,xcm2,ycm2,zcm2,am1,am2
      common/centersofmass/xcm1,ycm1,zcm1,xcm2,ycm2,zcm2,am1,am2

      integer i,nchk,nrelaxold,nnoptold,navold,ngrold
      integer corepts
      real*8 hcoold,hfloorold,tfold,dtoutold,
     $     alphaold,betaold,trelaxold
      logical restart
      character*3 iname
      namelist /initt/ iname
      common /inittcom/ iname
      integer corepts1
      common/core/corepts1
      real*8 divv(nmax)
      common/commdivv/divv
      real*8 amass1,amass2
      common/forcompbest/ amass1,amass2
      common /jumpcomm/ tjumpahead
      real*8 xcm,ycm,zcm,vxcm,vycm,vzcm,amtot
      real*8 erad,tjumpaheadold
      common/lostenergy/ erad
      integer iostatus,ncoolingold
      real*8 displacex,displacey,displacez
      integer ndisplace
      common/displace/displacex,displacey,displacez,ndisplace

      erad=0d0

      call energy_setup

c     read in hosts, find gravhost:
c      call rank_setup
c     set up hydro communicator:
c      call hydrocomm_setup
c     initialize "my" (parallel) arrays:
c      call myinit

c     get parameters from input file:
      call get_input
c     compute look-up tables:
      call tabulinit

      if(myrank.eq.0) write(69,*)'nav, nintvar, neos, ncooling=',
     $     nav,nintvar,neos,ncooling

      punit=gravconst*(munit/runit**2)**2

      if(myrank.eq.0) then
         write(69,*) 'These are the code units being used:'
         write(69,*) '   mass unit=',munit,'g'
         write(69,*) '   radius unit=',runit,'cm'
         write(69,*) '   pressure unit=',punit,'dyne/cm^2'
      endif  

c     check for existence of restart dump file:
      inquire(file='restartrad.sph',exist=restart)

c     if dump file exists, read it and restart from it:
      if(restart)then
         if(myrank.eq.0) then
            write(69,*)'init: continuing run'
            write(69,*)'init: reading restart dump file ...'
         endif
         corepts=0

         open(12,file='restartrad.sph',form='unformatted')
c     (the following read sequence must match exactly the write sequence
c     used in subroutine dump)
         read(12,iostat=iostatus) ntot,nnoptold,hcoold,hfloorold,sep0,tfold,
     $        dtoutold,nout,nit,t,navold,alphaold,betaold,tjumpaheadold,
     $        ngrold,nrelaxold,trelaxold,dt,omega2,ncoolingold,erad,
     $        ndisplace,displacex,displacey,displacez
         if(iostatus.ne.0) then
            ndisplace=0
            displacex=0d0
            displacey=0d0
            displacez=0d0
            if(myrank.eq.0) write(69,*)
     $           'Will try to read restartrad.sph again: trelaxold=',
     $           trelaxold
            close(12)
            open(12,file='restartrad.sph',form='unformatted')
            read(12,iostat=iostatus) ntot,nnoptold,hcoold,hfloorold,sep0,tfold,
     $           dtoutold,nout,nit,t,navold,alphaold,betaold,tjumpaheadold,
     $           ngrold,nrelaxold,trelaxold,dt,omega2,ncoolingold,erad
         endif
         if(iostatus.ne.0) then
            ncoolingold=0
            if(myrank.eq.0) write(69,*)
     $           'Will try to read restartrad.sph again: trelaxold=',
     $           trelaxold
            close(12)
            open(12,file='restartrad.sph',form='unformatted')
            read(12,iostat=iostatus) ntot,nnoptold,hcoold,hfloorold,sep0,tfold,
     $           dtoutold,nout,nit,t,navold,alphaold,betaold,tjumpaheadold,
     $           ngrold,nrelaxold,trelaxold,dt,omega2
         endif

         if(dtoutold.ne.dtout) then
            if(myrank.eq.0)write(69,*)'nout from restartrad.sph is',nout
            nout = t/dtout-1
            if(myrank.eq.0)write(69,*)'nout changed to',nout
         endif

         if(tjumpahead.gt.0) tjumpahead=-tjumpahead ! negative tjumpahead is so that tjumpahead will not be changed by changetf just because old eint value is not known in changetf

c         tf=sign(max(abs(tfold),abs(tf)),tf)

         if(myrank.eq.0) then
            write(69,*) myrank,'tjumpahead=',tjumpahead,'tf=',tf
         endif

c         if(nrelaxold.eq.0) then
c            iname='hyp'
c            open(16,file='m1m2rp.sph')
c            read(16,*) amass1,amass2
c            close(16)
cc            write(69,*) 'amass1,amass2=',amass1,amass2
c         endif
         omeg=sqrt(omega2)
         if(myrank.eq.0) write(69,*)'restartrad: omeg=',omeg

         if(ncoolingold.eq.0) then
            do i=1,ntot
               read(12) x(i),y(i),z(i),am(i),hp(i),rho(i),
     $              vx(i),vy(i),vz(i),vxdot(i),vydot(i),vzdot(i),
     $              u(i),udot(i), !gx(i),gy(i),gz(i),
     $              grpot(i),meanmolecular(i),cc(i),divv(i)
               ueq(i)=0d0
               tthermal(i)=1d30
            enddo
         else
            do i=1,ntot
               read(12) x(i),y(i),z(i),am(i),hp(i),rho(i),
     $              vx(i),vy(i),vz(i),vxdot(i),vydot(i),vzdot(i),
     $              u(i),udot(i), !gx(i),gy(i),gz(i),
     $              grpot(i),meanmolecular(i),cc(i),divv(i),
     $              ueq(i),tthermal(i)
            enddo
         endif
         read(12) nchk
         close(12)

         corepts1=(corepts+1)/2
         n=ntot-corepts
         if(nrelax.ge.2) then

c     use these lines of code if doing binaries for triple collisions:
            n1=1
            do while(cc(n1).eq.cc(1))
               n1=n1+1
               if(u(n1).eq.0.d0) goto 168
            enddo
 168        continue
            n1=n1-1
            if(myrank.eq.0) write(69,*)
     $         'number of particles in star 1 is n1=',n1
         else
            n1=ntot
         endif
         
c         if(nrelax.eq.1 .and. tresplintmuoff.gt.t) then
c            if(myrank.eq.0) write(69,*)'calling splinesetup...'
c            call splinesetup
c         endif

c         if(trelax.ne.trelaxold) then
c            if(myrank.eq.0) write(69,*)'***using trelax=',trelaxold
c            trelax=trelaxold
c         endif
c         if(nrelax.ne.nrelaxold) then
c            if(myrank.eq.0) write(69,*)'***using nrelax=',nrelaxold
c            nrelax=nrelaxold                  
c         endif
         if(nrelax.ne.2 .and. nrelaxold.eq.2) then
            nout=0
            nit=0
            t=0.d0
            if(treloff.le.0) then
               if(myrank.eq.0) then
                  write(69,*)
     $                 'assuming a dynamical calculation is starting.'
                  write(69,*)'resetting nout=nit=t=0'
                  write(69,*)'setting trelax=1.d30'
               endif
               trelax=1.d30
               if(alpha.le.0.d0 .and. beta.le.0.d0) then
                  alpha=1.d0
                  beta=2.d0
                  if(myrank.eq.0) then
                     write(69,*)'resetting alpha,beta=',alpha,beta
                  endif
               endif
            endif
         endif

         if(omega2.eq.0.d0 .or. (nrelax.eq.3 .and. treloff.le.0))
     $        gonedynamic=.true.
         if(myrank.eq.0) write(69,*)'gonedynamic=',gonedynamic
         if(myrank.eq.0) write(69,*)'n=',n,'ntot=',ntot
         
         if(nchk.ne.ntot) stop 'init: problem with dump file ???'
         
         if(myrank.eq.0) write(69,*) 'about to call stride_setup'
         
         call stride_setup
         
         dth=0.5d0*dt

         call gravquant         ! must call for inital set up
         
c         call cpu_time(time1)
c         call linkedlists
c         call cpu_time(time2)
c         if(myrank.eq.0) write(69,*) 'linkedlists lasted',time2-time1
c         do i=n_lower,n_upper
c            call newnene(i)
c         enddo
c         call cpu_time(time1)
c         if(myrank.eq.0) write(69,*) 'newnene calls lasted',time1-time2

         call pressure
      
         if(nrelaxold.ge.2) then
            if(nrelax.lt.2) then
               if(myrank.eq.0) write(69,*) 
     $               'did you mean to have nrelax>=2?'
c               stop
            endif
            call getcoms
            if(.not.gonedynamic) then
               if(myrank.eq.0) write(69,*) 'sep0=',sep0
               sepfinal=sep0
               if(myrank.eq.0) write(69,*) 'sepfinal=',sepfinal
c     use the following lines if want to restart a scanning run into a fixed
c     separation run that relaxes until t=treloff before going
c     dynamic:
ccccccccccccccc               treloff=dble(nint(t))
               sep0=dabs(xcm2-xcm1)
               sepfinal=sep0
               if(myrank.eq.0) write(69,*) 'sep0 & sepfinal set to',sep0
            endif

            hfloor=hfloorold
            if(myrank.eq.0) write(69,*) 'set hfloor=hfloorold=',hfloor

         endif
   
         if(myrank.eq.0) write(69,*)'init:         ... done(ish)'
         
c     otherwise create new data set:
      else
         nout=0
         nit=0
c         write(69,*) "setting t=0", myrank
         t=0.d0
c     initialize output parameters:
c     get 3-letter code for type of initial condition from init file
         open(12,file='sph.init',err=100)
         read(12,initt)
         close(12)
         if(myrank.eq.0) write(69,*)'init: new run, iname=',iname
         if(iname.eq.'1es') then
            call polyes
         else if (iname.eq.'1mc') then
            call polymces
         else if (iname.eq.'2cr') then
            call corotating
         elseif(iname.eq.'bps')then
            call bps            !binary plus single
         elseif(iname.eq.'bph')then
            call bpbh           !binary plus black hole
         elseif(iname.eq.'hyp')then
            call hyperbolic
         elseif(iname.eq.'hbs')then
            call hyperbolic_binary_single
         elseif(iname.eq.'erg')then
            call parent
            if(myrank.eq.0) write(69,*)'init: tf=',tf,dtout
         elseif(iname.eq.'tri')then
            call triple
         elseif(iname.eq.'bhe')then
            call smbh
         elseif(iname.eq.'rin')then
            call readin
         elseif(iname.eq.'grs')then
            call grsph
         elseif(iname.eq.'txt')then
            call asciiimage
         else
            if(myrank.eq.0) write(69,*)'init: unknown name',iname
            stop
         endif
      endif

      if(n.gt.nmax)then
         if(myrank.eq.0) write(69,*)'init: increase nmax to',n
         stop
      endif

      if(myrank.eq.0)then
         xcm=0d0
         ycm=0d0
         zcm=0d0
         vxcm=0d0
         vycm=0d0
         vzcm=0d0
         amtot=0d0
         do i=1,ntot
            xcm=xcm+am(i)*x(i)
            ycm=ycm+am(i)*y(i)
            zcm=zcm+am(i)*z(i)
            vxcm=vxcm+am(i)*vx(i)
            vycm=vycm+am(i)*vy(i)
            vzcm=vzcm+am(i)*vz(i)
            amtot=amtot+am(i)
         enddo
         xcm=xcm/amtot
         ycm=ycm/amtot
         zcm=zcm/amtot
         vxcm=vxcm/amtot
         vycm=vycm/amtot
         vzcm=vzcm/amtot
         write(69,*)'total mass=',amtot
         write(69,'(a,9g15.7)')'center of mass position=',xcm,ycm,zcm
         write(69,'(a,9g15.7)')'center of mass velocity=',vxcm,vycm,vzcm
      endif
         
      do i=n+1,ntot
         nn(i)=0
      enddo

c     write run parameters:
c here !!!!!
      if(myrank.eq.0) then
         write(69,*)'init: t=',t,' nit=',nit
         write(69,101) n,nnopt,hco,hfloor,dtout,nout,tf,sep0,nav,
     $     alpha,beta,ngr,nrelax,trelax,dt,ntot-n
 101     format (' init: parameters for this ideal + radiation run:',/,
     $        ' n=',i7,' nnopt=',i4,' hco=',g12.4,' hfloor=',g12.4,/,
     $        ' dtout=',g10.3,' nout=',i4,' tf=',g10.3,/,
     $        ' sep0=',g12.4,/,
     $        ' nav=',i2,' alpha=',f6.2,' beta=',f6.2,/,
     $        ' ngr=',i3,/,
     $        ' nrelax=',i2,' trelax=',g12.4,
     $        ' dt=',g12.4,
     $        ' corepts=',i2,/)
         if(nintvar.eq.1) then
            write(69,*)'integrating entropic variable a'
         elseif(nintvar.eq.2)then
            write(69,*)'integrating energy density u'
         else
            write(69,*)'must integrate either a or u'
            stop
         endif
         if(neos.eq.0) then
            write(69,*)'polytropic equation of state'
         elseif(neos.eq.1)then
            write(69,*)'ideal gas + radiation pressure eos'
         elseif(neos.eq.2)then
            write(69,*)'tabulated equation of state'
            if(nintvar.ne.2) then
               write(69,*)'must integrate u'
               stop
            endif
         else
            write(69,*)'invalid neos=',neos
            stop
         endif
         if(nusegpus.eq.0)then
            write(69,*)'cpus will be used for any gravity'
         else
            write(69,*)'gpus will be used for gravity'
         endif
         if(nselfgravity.eq.1)then
            write(69,*)'the gas is self-gravitating'
         else
            write(69,*)'the gas is *not* self-gravitating'
         endif
         write(69,*)'courant numbers=',cn1,cn2,cn3,cn4,cn5,cn6,cn7
         if(reat.ge.0.d0) write(69,*)'eating radius reat=',reat
         write(69,*)'init:            ... done'
         if(ncooling.eq.0)then
            write(69,*) 'no radiative cooling'
         else
            write(69,*) 'radiative cooling will be implemented'
         endif
         if(ndisplace.eq.0)then
            write(69,*)'no displacing origin to be near SPH particles'
         else
            write(69,*) 'displacing origin to be near SPH particles'
            write(69,'(4g17.9)') 'displace{x,y,z}=',
     $           displacex,displacey,displacez
         endif
      endif

      return

c     error condition:
 100  stop 'init:  error reading input file ???'
      end
************************************************************************
c      subroutine rank_setup
c      include 'starsmasher.h'
c      integer i,j
c      logical strideexists,notfound,notinstride
c      character*15 allhosts(nhmax),hosts(nhmax)
c      real*8 allstride(nhmax)
c
c      gravhost=-1
c      inquire(file='stride.in',exist=strideexists)
c      if(strideexists)then
c         open(23,file='stride.in',status='old')
c         do i=1,nhmax
c            read(23,*,end=33) allhosts(i),allstride(i)
c         enddo
c 33      continue
c         close(23)
c         if(myrank.eq.0) write(69,*) 'read stride.in file',nprocs
c         scalesum=0.d0
c         open(24,file='hosts',status='old')
c         do i=1,nprocs
c            read(24,*) hosts(i)
c            notfound=.true.
c            notinstride=.true.
c            j=1
c            do while(notfound.and.j.le.nhmax)
c               if(hosts(i).eq.allhosts(j))then
c                  stride(i)=allstride(j)
c                  scalesum=scalesum+stride(i)
c                  if(stride(i).gt.0) lasthydroproc=i-1
c                  notfound=.false.
c                  notinstride=.false.
c               else
c                  j=j+1
c               endif
c            enddo
c            if(notinstride)then
c               write(69,*)'need to put ',hosts(i),' in stride.in'
c               stop
c            endif
c            if(hosts(i).eq.'toe') gravhost=i-1
c         enddo
c         close(24)
c      else
c         scalesum=0.d0
c         open(24,file='hosts',status='old')
c         do i=1,nprocs
c            read(24,*) hosts(i)
c            if(hosts(i).eq.'toe') then
c               stride(i)=0
c               gravhost=i-1
c            else
c               stride(i)=1
c               scalesum=scalesum+stride(i)
c               lasthydroproc=i-1
c            endif
c         enddo
c         close(24)
c      endif
c      
c      return
c      end
************************************************************************
c      subroutine hydrocomm_setup
cc     make new mpi_comm for machines doing hydro calculations
cc     (exclude toe if it's stride value is 0)
c      include 'starsmasher.h'
c      integer color
c      if(gravhost.gt.0 .and. stride(gravhost+1).eq.0)then
c         color=1
c      endif
c      
c      end
************************************************************************
      subroutine get_input
      include 'starsmasher.h'
      logical autotf
      common/autotfblock/autotf
      common/orbitalelements/e0,semimajoraxis
      real*8 rhocgs,mucgs
      common/ueqstuff/rhocgs,teq,mucgs
      common /jumpcomm/ tjumpahead
      logical fileexists
      integer ierr
      real*8 displacex,displacey,displacez
      integer ndisplace
      common/displace/displacex,displacey,displacez,ndisplace

      ndisplace=0
      displacex=0d0
      displacey=0d0
      displacez=0d0

      semimajoraxis=0.d0
      bimpact=-1.d30
      e0=-1.d30
      vinf2=1.d30

c     set some default values, so that they don't necessarily have to be set in the sph.input file:
      tf=50000                 ! desired final time to stop simulation
      dtout=100                ! how often an out*.sph files should be dumped
      n=100000                 ! desired number of particles.  if n<0 then |n|=number of particles *per solar mass*.  used only if making a new star.
      gflag=1                   ! set to 0 for g function from appendix of gaburov et al. (2010); set to 1 for a g function that works better when there are black holes
      nnopt=22+gflag                 ! controls neighbor number.  leave it at 22 to get almost 40 neighbors.
      nav=3                    ! artificial viscosity (av) flag.  leave it at 3 to get a hybrid balsara-monaghan av.
      alpha=1                  ! av coefficient for term linear in mu
      beta=2                   ! av coefficient for mu^2 term
      ngr=3                    ! gravity flag.  leave it at 3.  if your want no gravity, ngr=0 might still work.
      hco=1d0                  ! softening/smoothing length for compact object or core particle
      hfloor=0d0               ! hp(i) = hptilde(i) + hfloor, where hp(i)=smoothing length and hptilde(i) is used in eq.(A1) of GLPZ 2010.
      nrelax=1                 ! relaxation flag.  0=dynamical calculation, 1=relaxation of single star, 2=relaxation of binary in corotating frame with centrifugal force, 3=calculation rotating frame with centrifugal and coriolis forces
      trelax=1.d30             ! timescale for artificial drag force.  keep it very large to turn off the drag force, which seems best even in relaxation runs (as the av can do the relaxation).
      sep0=200                 ! initial separation of two stars in a binary or collision calculation
      equalmass=0              ! particle mass is proportional to rho^(1-equalmass), so equalmass=1 has equal mass particles and equalmass=0 is for constant number density.
      treloff=0                ! time to end a relaxation and switch to a dynamical calculation
      tresplintmuoff=0.        ! time to stop resplinting the mean molecular weight.  leave this at 0.
      nitpot=1                 ! number of iterations between evaluation of the gravitational potential energy.
      tscanon=0                ! time that the scan of a binary starts
      sepfinal=1.d30           ! final separation for the scan of a binary
      nintvar=2                ! 1=integrate entropic variable a, 2=integrate internal energy u
      ngravprocs=-2            ! the number of gravity processors (must be <= min(nprocs,ngravprocsmax))
      qthreads=0               ! number of gpu threads per particle. typically set to 1, 2, 4, or 8.  set to a negative value to optimize the number of threads by timing.  set to 0 to guess the best number of threads without timing.
      mbh=20d0                 ! mass of black hole
      runit=6.9599d10          ! number of cm in the unit of length.  use 6.9599d10 if want solar radius.
      munit=1.9891d33          ! number of g in unit of mass.  use 1.9891d33 if want solar mass.
!     the courant numbers cn1, cn2, cn3, and cn4 are for sph particles:
!     dt_sph=1/(1/dt1 + 1/dt2 + 1/dt3 + 1/dt4)
      cn1=.6d0                 ! dt1=cn1*h/v_signal
      cn2=0.06d0               ! dt2=cn2*(h/|a-a_smoothed|)^0.5
      cn3=0.06d0               ! dt3=cn3*u/|du/dt|
      cn4=1.d30                ! dt4=cn4*v_signal/|a-a_smoothed|
!     the courant numbers cn5, cn6, and cn7 are for a particle i that is a compact object (co):
!     dt_co=1/(1/dt5 + 1/dt6)
      cn5=0.02d0               ! dt5=cn5*r_ij/v_ij (minimized over all other particles j)
      cn6=0.02d0               ! dt6=cn6*(r_ij/a_ij)^.5 (minimized over all other particles j)
      cn7=4.d0                 ! r_ij=(x_ij^2+y_ij^2+z_ij^2+cn7*h_i^2)^.5
!     the final timestep dt is the minimum of dt_sph and dt_co for all particles i
      computeexclusivemode=0   ! set this to 1 if on machine like grapefree with gpus in compute exclusive mode; set this to 0 on supercomputers like lincoln
      omega_spin=0.d0 ! angular rotation rate of star, used in nrelax=1 relaxations to give a rigidly rotating model
      ppn=12
      neos=1 ! 0 for polytropic equation of state (eos), 1 for ideal gas + radiation pressure, 2 for tabulated eos
      nselfgravity=1 ! 0 if just do gravity to point particles, 1 if self-gravitating
      gam=5.d0/3.d0 ! leave this set at a reasonable value even if using neos=1 or 2 (because the value of gam is used in estimating the local sound speed in balav3.f)
      reat=-1.d0
      starmass=1d0
      starradius=1d0
      ncooling=0 ! 0 if no cooling, otherwise radiative cooling
      nkernel=0
      teq=100d0
      tjumpahead=1d30
      startfile1='sph.start1u'
      startfile2='sph.start2u'
      eosfile='sph.eos'
      opacityfile='sph.opacity'
      profilefile='eg.last1.muse_s2mm'
      throwaway=.false.
      stellarevolutioncodetype=1

      open(12,file='sph.input',err=100)
      read(12,input)
      close(12)

      if(startfile2.ne.'sph.start2u') then
         inquire(file=startfile2,exist=fileexists)
         if(.not. fileexists) then
            write(69,*) 'startfile2 given in sph.input does not exist'
            call mpi_finalize(ierr)            
            stop
         endif
      endif

      if(myrank.eq.0 .and. ncooling.gt.0) write(69,*)
     $     'background temperature teq=',teq

      call set_nusegpus         ! if using gpus, this sets nusegpus=1 *and* nselfgravity=1

      if(cn1.lt.0.d0 .or. cn2.lt.0.d0 .or. cn3.lt.0.d0 .or.
     $     cn4.lt.0.d0 .or. cn5.lt.0.d0 .or. cn6.lt.0.d0 .or.
     $     cn7.lt.0.d0) then
         write(69,*) 'strange cn:',cn1,cn2,cn3,cn4,cn5,cn6,cn7
         stop
      endif

      if(ngr.ne.0 .and. nusegpus.eq.0) then
         ngravprocs=nprocs
         if(myrank.eq.0) then
            write(69,*)'using cpus to calculate gravity w/ ngravprocs=',
     $           ngravprocs
         endif
      endif

c      if(myrank.eq.0 .and. ngr.gt.0 .and. ngravprocs.eq.0) then
c         write(69,*) 'need to have at least one gravity process'
c         write(69,*) 'make sure ngravprocs is set in sph.input'
c         stop
c      endif

      if(ngr.eq.0 .and. ngravprocs.ne.0) then
         ngravprocs=0
         if(myrank.eq.0)
     $        write(69,*)'because ngr=0, we are setting ngravprocs=0'
      endif

      autotf=.false.
      if(tf.lt.0.d0) then
         autotf=.true.
         tf=abs(tf)
         open(23,file='ecc.sph')
         open(34,file='jumpahead.sph')
      endif

      return

 100  stop 'init: error reading input file'

      end
************************************************************************
      subroutine stride_setup
c     assign n_lower and n_upper to each processor
      include 'starsmasher.h'
c      integer i,stridesum,intstride(nprocs)
      integer i

c      write(69,*)'nprocs=',nprocs
c      write(69,*)'lasthydroproc=',lasthydroproc

c      if(nprocs.eq.1)then
c         n_lower=1
c         n_upper=n
c      else
c         stridesum=0
c         do i=1,myrank+1
c            n_lower=stridesum+1
c            if(i.eq.lasthydroproc+1.or.nprocs.eq.1)then
c               n_upper=n
c               intstride(i)=n_upper-n_lower+1
c            else
c               intstride(i)=int(stride(i)/scalesum*n+0.5d0)
c               n_upper=min(n_lower+intstride(i)-1,n)
c            endif
c            stridesum=stridesum+intstride(i)
c         enddo
c      endif

      if(nrelax.ne.1)then
         n=ntot
      else
         n=n
      endif

c      if(myrank.eq.0) write(69,*)'n_lower,n_upper,n',n_lower,n_upper,n

      do i=1,n
         vxdotsm(i)=0.d0
         vydotsm(i)=0.d0
         vzdotsm(i)=0.d0
      enddo

      return
      end
************************************************************************
      subroutine energy_setup
      include 'starsmasher.h'
      logical energyfilealreadyexists!,strideexists
      integer filenum
      character*11 energyfile
      character*8 logfile

      if(myrank.eq.0) then 
         energyfilealreadyexists=.true.
         filenum=-1
         do while(energyfilealreadyexists)
            filenum=filenum+1
            if(myrank.eq.0) write(energyfile,101)filenum
 101        format('energy',i1.1,'.sph')
            inquire(file=energyfile,exist=energyfilealreadyexists)
         enddo
      
         write(logfile,103)filenum
 103     format('log',i1.1,'.sph')
         open(22,file=energyfile,status='unknown')
         open(69,file=logfile,status='unknown')
         write(69,*)'writing energy data to ',energyfile
         write(69,*)'writing log data to ',logfile
      endif
      end
************************************************************************
      subroutine lfstart
c     prepare for first leap-frog iteration
      include 'starsmasher.h'                                          
      integer i!,status(mpi_status_size)

c     omega_spin>0 gives counterclockwise spin:
      if(omega_spin.ne.0.d0 .and. nrelax.eq.1) then
         if(myrank.eq.0) write(69,*) 'omega_spin=',omega_spin
         do i=1,ntot
            vx(i)=vx(i)-omega_spin*y(i)
            vy(i)=vy(i)+omega_spin*x(i)
         enddo
      endif

      call gravquant   ! must call for inital set up
      call rho_and_h
      if(ngr.ne.0) call gravforce
      call uvdots                ! evaluate udot and accelerations at half-timestep

      call enout(.true.)

c     advance velocities to half-timestep:
      call tstep

      dth=0.5d0*dt
c      do i=n_lower,n
c      do i=n_lower,n_upper
      do i=1,n
         vx(i)=vx(i)+vxdot(i)*dth
         vy(i)=vy(i)+vydot(i)*dth
         vz(i)=vz(i)+vzdot(i)*dth
      enddo
c      do i=n_lower,n_upper
      do i=1,n
         u(i)=u(i)+dth*udot(i)  ! update u to half-timestep
      enddo
      t=t+dth

      return
      end

************************************************************************
      subroutine readin
c     6/22/94 - all initial parameters except for n come from sph.input
      include 'starsmasher.h' 
      integer i,nchk,nrelaxold,nnoptold,noutold,nitold,navold,
     $     ngrold,corepts
      real*8 hcoold,hfloorold,sep0old,tfold,dtoutold,told,
     $     alphaold,betaold,trelaxold

c     read in data from previous run:      
      corepts=0
      if(myrank.eq.0) write(69,*)'readin: reading start file ...' 
      open(12,file='startu.sph',form='unformatted')
c     (the following read sequence must match exactly the write sequence
c     used in subroutine dump)
      read(12) ntot,nnoptold,hcoold,hfloorold,sep0old,
     $     tfold,dtoutold,noutold,nitold,told,navold,
     $     alphaold,betaold,tjumpahead,ngrold,nrelaxold,
     $     trelaxold,dt
      do i=1,ntot
         read (12) x(i),y(i),z(i),am(i),hp(i),rho(i),vx(i),vy(i),
     $        vz(i),vxdot(i),vydot(i),vzdot(i),u(i),udot(i),
     $        gx(i),gy(i),gz(i),grpot(i),meanmolecular(i),
     $        cc(i)
         if(u(i).eq.0.d0) then
            corepts=corepts+1
            if(myrank.eq.0) write(69,*)
     $         'should corepts really be set properly?'
            stop
         endif
      enddo
      read(12) nchk
      close(12)
      n=ntot-corepts
      
      if (nchk.ne.ntot) stop 'readin: problem with start file ???'
      call stride_setup
      if(myrank.eq.0) write(69,*)'readin:          ... done'
      return                                                           
      end
