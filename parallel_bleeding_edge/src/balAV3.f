      subroutine uvdots
      use riemann_data
c     evaluate right-hand sides of equations of motion:
c=======================================================================
c  balvdots creates arrays of the accelerations of each particle, using
c  balsara's form of the artificial viscosity pi_{ij}.
c  this incarnation of balvdots implements the form in which only the pi_{ij}
c     term is multiplied by the gradient of a symmetrized w_{ij}.
c     the p/rho^2 terms are each multiplied by the gradient of the kernel
c     evaluated with the smoothing length of that particle only.
c     see springel and hernquist's mnras article (astro-ph/0111016)
c     and monaghan's compressible turbulence article (astro-ph/0204118).
c=======================================================================
      include 'starsmasher.h'
      include 'mpif.h'
      double precision curlvxi,curlvyi,curlvzi
      real*8 csij,udbij,fi,ci,r2,
     $     abscurli,divvi,vxi,vyi,vzi,ci2,ami,por2i,h5,
     $     h2,hpi,csiju,csijgc
      integer itab,j,in,i,ierr,iw,nloop
      real*8 uijmax(nmax)
      common/uijmax/ uijmax
      real*8 dwijin,pijin,sxijin,syijin,szijin,divin, diwkin
      real*8 sxijingc,syijingc,szijingc
      logical switchedav
      double precision divv(nmax)
      common/commdivv/divv
      integer mylength, mygravlength
      integer comm_worker
      common/gravworkers/comm_worker
      real*8 myvxdot(ntot),myvydot(ntot),myvzdot(ntot)
      real*8 myvxdotgc(ntot),myvydotgc(ntot),myvzdotgc(ntot)
      integer irank
      real*8 invh, invh2, invh5
      real*8 dx, dy, dz
      real*8 dvx, dvy, dvz
      real*8 invbonnet
      real*8 alpha1, beta1
      integer offset, nni
c     variables used for radiative cooling portion of the code:
      real*8 uunit,udotunit,rhocgs,ucgs,temperature,kappa
      real*8 opacity,pgas,prad
      real*8 gradponrho,scaleheight,columnden
      real*8 kappar
      real*8 ueqfunction,useeostable,uupperlimit!,ulowerlimit
      external ueqfunction
      real*8 zeroin,mucgs
      common/ueqstuff/rhocgs,teq,mucgs
      real*8 grpottot(nmax),uraddoti
      save grpottot
      real*8 hpitilde,h2tilde,invh2tilde,htil
      integer itabtilde
c     table position in units of the table spacing, and its fraction
      real*8 xq,fq,xqtilde,fqtilde
c     nav=4: the pair's Riemann problem along r_ij (see riemann.f)
      real*8 rij,vir,vjr,pstar,pstat,vsig,cj
      real*8 rhl,rhr,vl,vr,pl,pr

      if(nav.eq.0) then
         if(nintvar.eq.1 .or. nintvar.eq.3)then
!$omp parallel do simd schedule(static)
            do i=1,n
               udot(i)=0.d0
            enddo
         else
            print *,'turn off av by setting alpha=beta=0'
            stop
         endif
      endif

c     nav=4: dissipation from a Riemann solver instead of alpha and beta.
c     The pair problem needs the density, pressure, Gamma_1 and gradients
c     of every neighbour.  With block timesteps these are shared in
c     rho_and_h, otherwise only this rank's own particles hold them.
      if(nav.eq.4) then
         call riemann_alloc
         if(.not.blkdist .and. nblock.eq.0) then
            call sync_all(rho)
            call sync_all(por2)
            if(neos.ne.0) call sync_all(gam1p)
         endif
      endif

      switchedav=.false.
      if(nrelax.eq.1
     $     .and. t.lt.treloff .and. alpha.eq.0.d0
     $     .and. beta.eq.0.d0) then
         switchedav=.true.
         alpha=1.d0
         beta=2.d0
         if(myrank.eq.0)
     $        write(69,*)'vdots: during relaxation alpha,beta=',alpha,beta
      endif

c  initialize dv/dt before accumulating:
!$omp parallel do schedule(static) private(i)
      do i=1,n
         myvxdot(i)=0.d0
         myvydot(i)=0.d0
         myvzdot(i)=0.d0
         myVXDOTgc(I)=0.d0
         myVYDOTgc(I)=0.d0
         myVZDOTgc(I)=0.d0
c     block steps: an inactive particle keeps its uijmax (tstep uses it)
         if(actblk(i)) uijmax(i)=0.d0
         udot(i)=0.d0
      end do

      if (nrelax.eq.1 .and.t.lt.treloff) then
         if(myrank.eq.0)write(69,*)
     $        'baludots is not using pi_{ij} right now'
      endif

c     the factor of 2d0 in the following lines is because our pi_{ij} below
c     is going to involve only p_i/rho_i^2 and not, as is often done,
c     (p_i/rho_i^2 + p_j/rho_j^2).  to keep alpha and beta values at the same
c     strength as in other implementations of sph, we multiply by 2 here.
      alpha1 = -2.d0 * alpha
      beta1  =  2.d0 * beta

c     for each particle (block steps with blkdist: this rank's round-robin
c     share of the particles to be solved, otherwise its own index range):
      if(blkdist) then
         nloop=nwmine
      else
         nloop=n_upper-n_lower+1
      endif
c     the pair loop, by chunks of particles on the threads (uvpairs, below)
      call uvpairs(nloop,alpha1,beta1,myvxdot,myvydot,myvzdot,
     $     myvxdotgc,myvydotgc,myvzdotgc)
c      write(6,'(a)')'hydrompi'

      if(nblock.eq.1) then
c     block steps: a point particle is never refreshed, so when it is
c     inactive it has no list this substep.  With dynamic softening, its
c     kernel term on active gas inside it is found directly.
         if(dynhco.ge.1 .and. ngr.ne.0) then
            do i=n_lower,n_upper
               if(u(i).ne.0.d0 .or. actblk(i) .or. refblk(i)) cycle
c     grad G at htilde, the length that solves eq. (a1)
               htil=hp(i)-hfloor
               if(dynhco.eq.3) htil=hdynco(i)
               invh2=1.d0/htil**2
               do j=1,ntot
                  if(.not.actblk(j) .or. u(j).eq.0.d0) cycle
                  dx=x(i)-x(j)
                  dy=y(i)-y(j)
                  dz=z(i)-z(j)
                  r2=dx*dx+dy*dy+dz*dz
                  if(r2.ge.4.d0*htil**2) cycle
                  xqtilde=ctab*r2*invh2
                  itabtilde=int(xqtilde)+1
                  fqtilde=xqtilde-dble(itabtilde-1)
                  diwkin=(dgtab(itabtilde)
     $                 +fqtilde*(dgtab(itabtilde+1)-dgtab(itabtilde)))
     $                 *invh2
                  if(hcolim) then
                     if(r2.lt.4.d0*hdynco(i)**2) then
                        xqtilde=ctab*r2/hdynco(i)**2
                        itabtilde=int(xqtilde)+1
                        fqtilde=xqtilde-dble(itabtilde-1)
                        diwkin=(dgtab(itabtilde)
     $                       +fqtilde*(dgtab(itabtilde+1)-dgtab(itabtilde)))
     $                       /hdynco(i)**2
                     else
                        diwkin=0.d0
                     endif
                  endif
                  csijgc=-0.5d0*diwkin*bonet_psi(i)/(bonet_0mega(i)*am(j))
                  myvxdotgc(j)=myvxdotgc(j)+am(i)*csijgc*dx
                  myvydotgc(j)=myvydotgc(j)+am(i)*csijgc*dy
                  myvzdotgc(j)=myvzdotgc(j)+am(i)*csijgc*dz
               enddo
            enddo
         endif
c     div v of the particles solved here is needed on every rank (the h and
c     rho of inactive particles are predicted from it)
         if(blkdist) then
            call sync_w(divv)
            call sync_w(uijmax)
         else
            call sync_all(divv)
         endif
      endif

      mylength=n_upper-n_lower+1
c     one collective gives every rank the whole array (formerly one gather per rank)
      if(blkdist) then
         call sync_w(udot)
      else if(nprocs.gt.1) then
         call mpi_allgatherv(mpi_in_place,0,mpi_datatype_null,udot,
     $        recvcounts,displs,mpi_double_precision,mpi_comm_world,ierr)
      endif

      if(ncooling.ne.0) then
c     to get the pressure scale height, we will need the v{x,y,z}dot arrays *without* the
c     gravitational acceleration included.  note: the myv{x,y,z}dot arrays are not
c     affected by the following mpi calls.  instead, they are just being summed up to give
c     the v{x,y,z}dot arrays.
         call sumvdot(myvxdot,myvydot,myvzdot)
         
         punit=gravconst*(munit/runit**2)**2
         uunit=gravconst*munit/runit
         udotunit=(gravconst*munit)**1.5d0/runit**2.5d0
c     if(myrank.eq.0) write(69,*) 'uunit,udotunit=',uunit,udotunit

         if(neos.eq.1)then
            do i=n_lower,n_upper
               if(u(i).ne.0.d0) then
                  rhocgs=rho(i)*munit/runit**3.d0
                  if(nintvar.eq.1) then
                     ucgs=u(i)*rho(i)**(gam-1.d0)/(gam-1.d0)
     $                    *uunit
                  else
                     ucgs=u(i)*uunit
                  endif
                  call gettemperature(qconst*rhocgs/meanmolecular(i),
     $                 -ucgs*rhocgs/arad,temperature)
                  pgas=rhocgs*boltz*temperature/meanmolecular(i)/punit
                  prad=arad*temperature**4/3.d0/punit
                  
                  gradponrho=(vxdot(i)**2+vydot(i)**2+vzdot(i)**2)**0.5d0
                  scaleheight=(pgas+prad)/(gradponrho*rho(i))
                  columnden=rho(i)*scaleheight*1.06d0*munit/runit**2
c                 kappa is local opacity; kappar is pseudo-mean opacity
                  call usetable(rhocgs,temperature,kappa,kappar)
                  uraddoti=4*sigma*(teq**4-temperature**4)/
     $                 (columnden**2*kappar + 1/kappa)
     $                 /udotunit
                  ueq(i)=(arad*teq**4/rhocgs
     $                 +1.5d0*boltz*teq/meanmolecular(i))
     $                 /uunit
                  
                  tthermal(i)=(ueq(i)-u(i))/uraddoti

c     approximate the equilibrium temperature... not really correct...
c     teq=(10.d0**4+udot(i)*(columnden**2*kappa+ 1/kappa)
c     $           *udotunit/(4*sigma))**0.25d0
                  
               else
                  tthermal(i)=1.d30
                  ueq(i)=0.d0
               endif
            enddo
         else if(neos.eq.2)then
            do i=n_lower,n_upper
               if(u(i).ne.0.d0) then
                  rhocgs=rho(i)*munit/runit**3.d0
                  mucgs=meanmolecular(i)
                  if(nintvar.eq.1) then
                     ucgs=u(i)*rho(i)**(gam-1.d0)/(gam-1.d0)
     $                    *uunit
                  else
                     ucgs=u(i)*uunit
                  endif
                  temperature = useeostable(ucgs,rhocgs,
     $                 meanmolecular(i),1)
                  gradponrho=(vxdot(i)**2+vydot(i)**2+vzdot(i)**2)**0.5d0
                  scaleheight=useeostable(ucgs,rhocgs,
     $                 meanmolecular(i),3)/(punit*gradponrho*rho(i))
                  columnden=rho(i)*scaleheight*1.06d0*munit/runit**2
c                 kappa is local opacity; kappar is pseudo-mean opacity
                  call usetable(rhocgs,temperature,kappa,kappar)
                  uraddoti=4*sigma*(teq**4-temperature**4)/
     $                 (columnden**2*kappar + 1/kappa)
     $                 /udotunit
                  
c     using eos table, solve for u_eq such that temperature t=teq 
                  uupperlimit=1.1d0*(1.5d0*boltz*teq/meanmolecular(i)
     $                 +arad*teq**4/rhocgs)
 1234             continue
                  ucgs=zeroin(0d0,uupperlimit,ueqfunction,1.d-11)
                  if(ueqfunction(ucgs).gt.1d-6) then
                     uupperlimit=2*uupperlimit
                     write(69,*)'Having trouble with u_eq for particle'
     $                    ,i,'at time',t,ueqfunction(ucgs)
                     ucgs=zeroin(0d0,uupperlimit,ueqfunction,1.d-11)
                     write(69,*)'Found',ucgs,'for'
     $                    ,i,'at time',t,ueqfunction(ucgs)
                  endif
                  ueq(i)=ucgs/uunit

                  tthermal(i)=(ueq(i)-u(i))/uraddoti

               else
                  tthermal(i)=1.d30
                  ueq(i)=0
               endif
            enddo
         else
            print *,'radiative cooling requires neos=1 or 2'
            stop
         endif

         mylength=n_upper-n_lower+1
c     one collective gives every rank the whole array (formerly one gather per rank)
         call mpi_allgatherv(mpi_in_place,0,mpi_datatype_null,tthermal,
     $     recvcounts,displs,mpi_double_precision,mpi_comm_world,ierr)
         call mpi_allgatherv(mpi_in_place,0,mpi_datatype_null,ueq,
     $     recvcounts,displs,mpi_double_precision,mpi_comm_world,ierr)
      endif

c      write(6,'(a)')'hydrompi_complete'
      if(ngr.ne.0 .and. myrank.lt.ngravprocs)then
         if(myrank.eq.ngravprocs-1) call wall_time(time0)
         if(nusegpus.eq.1)then
            call grav_collect(ntot, gx, gy, gz, grpot)
         else
c     block steps: gravity only on the active particles, unless at least
c     90% are active, when the symmetric pair loop over all pairs is cheaper
            if(nblock.eq.1 .and. count(actblk(1:ntot)).lt.0.9d0*ntot) then
               call get_gravity_active
            else
               call get_gravity_using_cpus
            endif
         endif   
         if(myrank.eq.ngravprocs-1) call wall_time(time1)
c     (with block steps the potential is needed only when every particle is active, for the energies)
         if(ngravprocs.gt.1 .and. (nblock.ne.1 .or. count(actblk(1:ntot)).eq.ntot)) then
            if(nusegpus.eq.1)then
               mygravlength=ngrav_upper-ngrav_lower+1
               if(myrank.ne.0) then
                  call mpi_gatherv(grpot(ngrav_lower), mygravlength, mpi_double_precision,
     $                 grpot, gravrecvcounts, gravdispls, mpi_double_precision, 0,
     $                 comm_worker, ierr)
               else
                  call mpi_gatherv(mpi_in_place, mygravlength, mpi_double_precision,
     $                 grpot, gravrecvcounts, gravdispls, mpi_double_precision, 0,
     $                 comm_worker, ierr)
               endif
            else
               call mpi_reduce(grpot, grpottot, n, mpi_double_precision,
     $              mpi_sum, 0, mpi_comm_world, ierr)
               if(myrank.eq.0) then
                  do i=1,n
                     grpot(i)=grpottot(i)
                  enddo
               endif
            endif
         endif
         if(myrank.eq.ngravprocs-1)then
            write(6,'(a,f6.3,a,g10.3,a,i4)')'2ndhalf:',time1-time0
         endif
         if(nusegpus.eq.1)then
!$omp parallel do simd schedule(static)
            do i=ngrav_lower,ngrav_upper
               myvxdot(i)=myvxdot(i)+gx(i)
               myvydot(i)=myvydot(i)+gy(i)
               myvzdot(i)=myvzdot(i)+gz(i)
            enddo
         else
!$omp parallel do simd schedule(static)
            do i=1,n
               myvxdot(i)=myvxdot(i)+gx(i)
               myvydot(i)=myvydot(i)+gy(i)
               myvzdot(i)=myvzdot(i)+gz(i)
            enddo
         endif
      endif

ccccccccccccccccccccccccccccccccccccccccc
c     reassemble acceleration arrays here
ccccccccccccccccccccccccccccccccccccccccc

c     use mpi_allreduce to collect subtotals, sum them, and redistribute
c     to all processes

c     Now that scaleheight has been calculated, we can include all
c     relevant terms in the calculation of the acceleration:
!$omp parallel do simd schedule(static)
      do i=1,n
         myvxdot(i)=myvxdot(i)+myvxdotgc(i)
         myvydot(i)=myvydot(i)+myvydotgc(i)
         myvzdot(i)=myvzdot(i)+myvzdotgc(i)
      enddo

      call sumvdot(myvxdot,myvydot,myvzdot)

cc     if this is a relaxation calculation, add centrifugal and drag forces:
c      if (nrelax.ne.0 .and. .not. gonedynamic ) then
          call relax        ! add drag and fictitious forces if necessary
c     endif

      if(switchedav)then
         alpha=0.d0
         beta=0.d0
      endif

c      write(69,*)'balav3: finished' 

      return
      end
c***********************************************************************
      subroutine sumvdot(ax,ay,az)
c     add the accelerations of all the ranks into vxdot, vydot and vzdot,
c     with one collective instead of one for each component
      include 'starsmasher.h'
      include 'mpif.h'
      real*8 ax(n),ay(n),az(n)
      real*8, allocatable, save :: buf(:)
      real*8, allocatable :: wbuf(:)
      integer i,k,ierr

c     with one rank there is nothing to add: the accelerations are the ones given
      if(nprocs.eq.1) then
!$omp parallel do simd schedule(static)
         do i=1,n
            vxdot(i)=ax(i)
            vydot(i)=ay(i)
            vzdot(i)=az(i)
         enddo
         return
      endif
      if(blkdist .and. nblock.eq.1) then
c     block steps: only the accelerations of the listed particles (active, or refreshed: wall) are used afterwards, so only those
c     are added up; the others keep what this rank has for them, which is replaced before anyone reads it
!$omp parallel do simd schedule(static) private(i)
         do i=1,n
            vxdot(i)=ax(i)
            vydot(i)=ay(i)
            vzdot(i)=az(i)
         enddo
         if(nwall.eq.0) return
         allocate(wbuf(3*nwall))
!$omp parallel do schedule(static) private(i,k)
         do k=1,nwall
            i=wall(k)
            wbuf(k)=ax(i)
            wbuf(nwall+k)=ay(i)
            wbuf(2*nwall+k)=az(i)
         enddo
         call mpi_allreduce(mpi_in_place,wbuf,3*nwall,mpi_double_precision,
     $        mpi_sum,mpi_comm_world,ierr)
!$omp parallel do schedule(static) private(i,k)
         do k=1,nwall
            i=wall(k)
            vxdot(i)=wbuf(k)
            vydot(i)=wbuf(nwall+k)
            vzdot(i)=wbuf(2*nwall+k)
         enddo
         deallocate(wbuf)
         return
      endif
      if(allocated(buf)) then
         if(size(buf).lt.3*n) deallocate(buf)
      endif
      if(.not.allocated(buf)) allocate(buf(3*n))
!$omp parallel do simd schedule(static)
      do i=1,n
         buf(i)=ax(i)
         buf(n+i)=ay(i)
         buf(2*n+i)=az(i)
      enddo
      call mpi_allreduce(mpi_in_place,buf,3*n,mpi_double_precision,
     $     mpi_sum,mpi_comm_world,ierr)
!$omp parallel do simd schedule(static)
      do i=1,n
         vxdot(i)=buf(i)
         vydot(i)=buf(n+i)
         vzdot(i)=buf(2*n+i)
      enddo

      return
      end
c***********************************************************************
      subroutine uvchunk(ic,nloop,alpha1,beta1,mvx,mvy,mvz,mgx,mgy,mgz,
     $     smap,bm,acc)
c     The pair loop of uvdots for chunk ic of the particles of this rank
c     (pb_chunk consecutive ones): the terms of each particle i, which go to i
c     into mvx.. and mgx.. (the "gather" part) and to its neighbours j into the
c     buffer of the chunk, smap/bm/acc being the working space of the thread
c     (see pairbuf.f90).  Called from a parallel region, so every variable
c     below is local to the call.
      use riemann_data
      use pairbuf_mod
      include 'starsmasher.h'
      include 'mpif.h'
      integer ic,nloop
      real*8 alpha1,beta1
      real*8 mvx(*),mvy(*),mvz(*),mgx(*),mgy(*),mgz(*)
      integer smap(*)
      integer*8 bm(*)
      real*8 acc(6,*)
      integer iw0,iw1,cnt,jmin,jmax,jslot,jw
      double precision curlvxi,curlvyi,curlvzi
      real*8 csij,udbij,fi,ci,r2,
     $     abscurli,divvi,vxi,vyi,vzi,ci2,ami,por2i,h5,
     $     h2,hpi,csiju,csijgc
      integer itab,j,in,i,iw
      real*8 uijmax(nmax)
      common/uijmax/ uijmax
      real*8 dwijin,pijin,sxijin,syijin,szijin,divin, diwkin
      real*8 sxijingc,syijingc,szijingc
      double precision divv(nmax)
      common/commdivv/divv
      real*8 invh, invh2, invh5
      real*8 dx, dy, dz
      real*8 dvx, dvy, dvz
      real*8 invbonnet
      integer offset, nni
      real*8 hpitilde,h2tilde,invh2tilde
      integer itabtilde
      real*8 xq,fq,xqtilde,fqtilde
      real*8 rij,vir,vjr,pstar,pstat,vsig,cj
      real*8 rhl,rhr,vl,vr,pl,pr

      iw0=(ic-1)*pb_chunk+1
      iw1=min(ic*pb_chunk,nloop)
      cnt=0
      jmin=nmax
      jmax=0
      do iw=iw0,iw1
         if(blkdist) then
            i=wmine(iw)
         else
            i=n_lower+iw-1
         endif
c     block steps (nblock=1): only particles with a neighbour list from this
c     substep take part, i.e. the active ones and the inactive ones that
c     were refreshed because an active particle lies inside their kernel.
c     Each list gives the terms of particle i's kernel, which act on i and
c     on its neighbour j, and only active particles keep their results.
c     This finds every term on an active particle.  With shared steps every
c     particle is active (actblk), so nothing is skipped.
         if(.not.(actblk(i) .or. refblk(i))) cycle
         hpi=hp(i)
         hpitilde=hpi - hfloor
         if(dynhco.eq.3) hpitilde=hdynco(i)
         h2=hpi**2
         h2tilde=hpitilde**2
         h5=hpi*h2**2
         invh  = 1.0/hpi
         invh2 = 1.0/h2
         invh2tilde = 1.0/h2tilde
         invh5 = 1.0/h5
         ami=am(i)
         por2i=por2(i)
         if(nav.eq.4) then
            ci2=gam1p(i)*por2i*rho(i)
         else
            ci2=gam*por2i*rho(i)
         endif
         ci=sqrt(ci2)
         vxi=vx(i)
         vyi=vy(i)
         vzi=vz(i)
         call getderivs(i,curlvxi,curlvyi,curlvzi,divvi)
         divv(i)=divvi
         abscurli=sqrt(curlvxi**2+curlvyi**2+curlvzi**2)

c     calculate gradwij's to all neighbors:
         offset = first(i)
         nni    = nn(i)
         do in=1,nni
            j=list(offset+in)
            if(.not.(actblk(i) .or. actblk(j))) cycle
       
            dx = x(i) - x(j)
            dy = y(i) - y(j)
            dz = z(i) - z(j)
            if(periodic) call pbcwrap(dx,dy,dz)
            dvx = vxi - vx(j)
            dvy = vyi - vy(j)
            dvz = vzi - vz(j)
            r2 = dx*dx + dy*dy + dz*dz 
c     a capped point particle's list reaches beyond 2h (to 2*htilde)
            xq=ctab*min(r2*invh2,4.d0)
            itab=int(xq)+1
            fq=xq-dble(itab-1)

            if(r2.lt.4d0*h2tilde) then
               xqtilde=ctab*r2*invh2tilde
               itabtilde=int(xqtilde)+1
               fqtilde=xqtilde-dble(itabtilde-1)
               diwkin = (dgtab(itabtilde)
     $              +fqtilde*(dgtab(itabtilde+1)-dgtab(itabtilde)))
     $              * invh2tilde
            else
               diwkin = 0d0
            endif

            dwijin = (dwtab(itab)+fq*(dwtab(itab+1)-dwtab(itab))) * invh5

            if(u(j).ne.0.d0 .and. u(i).ne.0) then
              divin = dx*dvx + dy*dvy + dz*dvz
               if (divin.lt.0.d0 .and. nav.ne.4) then
c                  r2=(x(i)-x(j))**2.d0+(y(i)-y(j))**2.d0
c     $                 +(z(i)-z(j))**2.d0
c     cj=sqrt(5.d0/3.d0*por2(j)*rho(j))
                  
c     calculate dinshaw balsara's (db) \mu_{ij}: see equation (a12) of http://adsabs.harvard.edu/abs/2012apj...745...71p
                  fi=abs(divvi)/(abs(divvi)+abscurli+
     $                 0.00001d0*ci*invh)
c                 fj=abs(divv(j))/(abs(divv(j))+
c     $                    sqrt(curlvx(j)**2+curlvy(j)**2+curlvz(j)**2)+
c     $                    0.00001d0*cj/hp(j))
c     modified from very old version:
                  if(r2.gt.0) then
                     if(nav.eq.3) then
                        udbij=divin/(ci*r2**0.5d0)*fi
                     else
                        udbij=hpi*divin/(ci*r2)*fi
                     endif
                  else
                     udbij=0d0
                  endif

c     old version:
c     udbij=div(in)/((ci+cj)*r2**0.5d0)*(fi+fj)
c     udbij=divin/(ci*r2**0.5d0)*fi
                  if(abs(udbij).gt.1.d8) then
                     write(69,*)'c warning large udbij...',i,j, udbij,divin,ci,
     $                    am(i),am(j),u(i),u(j)
                     if(dt.lt.1d-17)then
                        write(69,*)'unreasonably small dt. abort run.'
!$omp atomic write
                        pb_fail=1
                     endif
                  endif
c                               a     +    b * c   (multiply-add, cpu likes these)
                  pijin=por2i*(alpha1 + beta1*udbij)*udbij
               else
                  pijin=0.d0
               endif

c     nav=4: the Riemann problem along r_hat=(r_i-r_j)/|r_i-r_j|, on which
c     j is the left state and i the right, each with its particle's own
c     density, pressure and velocity along r_hat.  pstat is the star
c     pressure of the same problem with the two particles at rest, so
c     p*-pstat is the part that comes from their relative motion, and
c     (p*-pstat)(1/rho_i^2+1/rho_j^2) takes the place of the AV.  It is
c     the same from both ends of the pair, so the dissipation conserves
c     energy pair by pair, and it is zero without relative motion, where
c     the equations reduce to the variational ones.  (Second-order states
c     reconstructed to the pair midpoint from limited gradients were
c     tried and dropped: the limiter switching them on and off broke
c     energy conservation by percents in black hole-star encounters.  On
c     shock tubes with navsw=2 they were within a few per cent of these
c     first-order states, and with navsw=1 up to 1.5 times more accurate.)
               if(nav.eq.4 .and. r2.eq.0.d0) then
                  vir=0.d0
                  vjr=0.d0
               elseif(nav.eq.4) then
                  rij=sqrt(r2)
                  vir=(vxi*dx+vyi*dy+vzi*dz)/rij
                  vjr=(vx(j)*dx+vy(j)*dy+vz(j)*dz)/rij
                  pl=por2(j)*rho(j)**2
                  pr=por2i*rho(i)**2
                  rhl=rho(j)
                  rhr=rho(i)
                  vl=vjr
                  vr=vir
c     navsw selects a switch: 2 or 3 drops receding pairs, 1 or 3
c     multiplies by balsara's f_i as in nav=3
                  if((navsw.eq.2 .or. navsw.eq.3) .and. divin.ge.0.d0) then
                     pijin=0.d0
                  else
                     call riemann_solve(rhl,vl,pl,gam1p(j),rhr,vr,pr,
     $                    gam1p(i),nriemann,pstar)
                     call riemann_solve(rhl,0.d0,pl,gam1p(j),rhr,0.d0,pr,
     $                    gam1p(i),nriemann,pstat)
c     the exact p* falls as vr-vl rises, so p*-pstat must have the sign of
c     vl-vr.  If the one-step solution does not, solve this pair to
c     convergence (rare: never in the pairs of the collision dumps tested)
                     if(nriemann.ne.0 .and.
     $                    (pstar-pstat)*(vr-vl).gt.0.d0) then
                        call riemann_solve(rhl,vl,pl,gam1p(j),rhr,vr,pr,
     $                       gam1p(i),0,pstar)
                        call riemann_solve(rhl,0.d0,pl,gam1p(j),rhr,0.d0,
     $                       pr,gam1p(i),0,pstat)
                     endif
                     pijin=(pstar-pstat)*(1.d0/rho(i)**2+1.d0/rho(j)**2)
                     if(navsw.eq.1 .or. navsw.eq.3) pijin=pijin*abs(divvi)/
     $                    (abs(divvi)+abscurli+0.00001d0*ci*invh)
                  endif
               endif
c     hydro correction part.  here we are implementing corrections as shown
c     in equations (a11) and (a12) of gaburov et al. (2010):
               invbonnet = 1.0/(bonet_0mega(i)*am(j))

c     csij will not contain AV (pijin) contributions so that we can
c     later use the acceleration derived from it to calculate scale height.
c     All AV contributions to the acceleration will be accounted for with
c     the help of csijgc below. (gc=gravity correction... this is for terms
c     associated with the gravity correction *and* with AV)
               csij = por2i * (-bonet_omega(i) * invbonnet * diwkin
     $              + dwijin)
               csijgc= 0.5d0 * pijin * dwijin

c     hydro part
               if (nrelax.eq.1 .and.t.lt.treloff) then
                  csiju= csij ! No heating from AV during relaxation
               else
                  csiju= csij + csijgc
               endif

               if(nintvar.eq.2) then
                  udot(i)=udot(i)+am(j)*csiju*divin
               elseif (nrelax.eq.0 .or. t.ge.treloff) then
                  udot(i)=udot(i)+am(j)*pijin*dwijin*divin
               endif

               if(nav.eq.4) then
c     on the scale of the AV version's p/rho (c^2/Gamma_1): half the
c     pair's signal speed, squared, over Gamma_1
                  cj=sqrt(gam1p(j)*por2(j)*rho(j))
                  vsig=ci+cj+max(vjr-vir,0.d0)
                  if(actblk(i)) uijmax(i)=max(uijmax(i),
     $                 (0.5d0*vsig)**2/gam1p(i))
               elseif(actblk(i)) then
                  uijmax(i)=max(uijmax(i),
     $              (por2i + 0.5d0 * pijin)*rho(i))
               endif
c               uijmax(j)=max(uijmax(j),
c     $              (por2i + 0.5d0 * pijin)*rho(i))

               if(ngr.ne.0) then
c     gravity correction part
c                  csij = csij - 0.5d0*diwkin * bonet_psi(i) * invbonnet
                  csijgc = csijgc - 0.5d0*diwkin * bonet_Psi(i) * invBonnet
               endif
            else if(u(i).eq.0.d0 .and. u(j).ne.0.d0 .and. dynhco.ge.1
     $              .and. ngr.ne.0) then
c     point particle i with a dynamic softening length: no hydro or AV, but
c     the same gravity correction (a14),(a15) as for an sph particle, since
c     its h_i obeys eq. (a1) summed over its sph neighbours j
               csij=0.d0
               if(hcolim) then
c     grad G at the unlimited solution h_dyn (see rho_and_h)
                  if(r2.lt.4.d0*hdynco(i)**2) then
                     xqtilde=ctab*r2/hdynco(i)**2
                     itabtilde=int(xqtilde)+1
                     fqtilde=xqtilde-dble(itabtilde-1)
                     diwkin=(dgtab(itabtilde)
     $                    +fqtilde*(dgtab(itabtilde+1)-dgtab(itabtilde)))
     $                    /hdynco(i)**2
                  else
                     diwkin=0.d0
                  endif
               endif
               csijgc = -0.5d0*diwkin*bonet_psi(i)/(bonet_0mega(i)*am(j))
            else
               csij=0.d0
               csijgc=0.d0
            endif
            
            sxijin = csij * dx !(x(i) - x(j))
            syijin = csij * dy !(y(i) - y(j))
            szijin = csij * dz !(z(i) - z(j))
            sxijingc = csijgc * dx !(x(i) - x(j))
            syijingc = csijgc * dy !(y(i) - y(j))
            szijingc = csijgc * dz !(z(i) - z(j))

c     finally compute each component of dv/dt:

c     note: sxijin has been calculated using pi/rhoi^2.  so in the "gather"
c     part of the sum, everything is as you would expect.  but the gather part
c     of the sum only gives half of the contributions for any given particle
c     i.  for the sake of concreteness, let's say i=5.  the other half of the
c     contributions to the acceleration of particle 5 comes from the scatter
c     portion of the code with j=5 (note: this is *j*, not i).  in the scatter
c     portion, the index i will be a number other than 5, and so the p/rho^2
c     being used inside of sxijin will come from that particle with index
c     different than 5.  the gather contributions to particle 5 all happen
c     together in one big group.  the scatter contributions to particle 5 will
c     happen haphazardly as the loop over particle index i gradually covers
c     all the particles that have 5 as a neighbor.

c     "gather" part of the sum:
            mvx(i)=mvx(i)-am(j)*sxijin
            mvy(i)=mvy(i)-am(j)*syijin
            mvz(i)=mvz(i)-am(j)*szijin
c     "scatter" part of the sum, which goes to the j of this chunk (see pairbuf.f90):
c     the slot of j in the buffer of the chunk, made when j is first touched
            jslot=smap(j)
            if(jslot.eq.0) then
               cnt=cnt+1
               smap(j)=cnt
               jw=ishft(j-1,-6)+1
               bm(jw)=ior(bm(jw),ishft(1_8,iand(j-1,63)))
               jslot=cnt
               acc(1,jslot)=0.d0
               acc(2,jslot)=0.d0
               acc(3,jslot)=0.d0
               acc(4,jslot)=0.d0
               acc(5,jslot)=0.d0
               acc(6,jslot)=0.d0
               jmin=min(jmin,j)
               jmax=max(jmax,j)
            endif
            acc(1,jslot)=acc(1,jslot)+ami*sxijin
            acc(2,jslot)=acc(2,jslot)+ami*syijin
            acc(3,jslot)=acc(3,jslot)+ami*szijin
c     "Gather" part of the sum:
            mgx(i)=mgx(i)-am(j)*sxijingc
            mgy(i)=mgy(i)-am(j)*syijingc
            mgz(i)=mgz(i)-am(j)*szijingc
c     "Scatter" part of the sum:
            acc(4,jslot)=acc(4,jslot)+ami*sxijingc
            acc(5,jslot)=acc(5,jslot)+ami*syijingc
            acc(6,jslot)=acc(6,jslot)+ami*szijingc
         enddo

c     p=(gam-1)*rho*u=a*rho^gam, so a=(gam-1)*rho^(1-gam)*u
         if(nintvar.eq.1 .and. u(i).ne.0.d0)
     $        udot(i)=udot(i)*0.5d0*(gam-1.d0)*rho(i)**(1.d0-gam)

c     ln A is the specific entropy up to the factor 3k/2mu and a constant, so
c     d(lnA)/dt = (2 mu/3k) (du/dt)/T, with du/dt the dissipative part alone --
c     which is exactly what the sum above holds, since the branch that fills it
c     leaves out the adiabatic pdV work that does not change the entropy.  The
c     temperature comes from pressure, which runs earlier in the same step over
c     this same n_lower:n_upper range; the fallback covers the first call,
c     before pressure has ever filled it.
         if(nintvar.eq.3 .and. u(i).ne.0.d0) then
            if(tempcache(i).le.0.d0)
     $           call getT_from_lna(u(i),rho(i)*munit/runit**3.d0,
     $           meanmolecular(i),-1.d0,tempcache(i))
            udot(i)=udot(i)*0.5d0*gravconst*munit/runit
     $           *2.d0*meanmolecular(i)/(3.d0*boltz*tempcache(i))
         endif

      enddo
      call pb_store(ic,cnt,jmin,jmax,smap,bm,acc)

      return
      end
c***********************************************************************
      subroutine uvpairs(nloop,alpha1,beta1,mvx,mvy,mvz,mgx,mgy,mgz)
c     the pair loop of uvdots, by chunks of particles shared among the threads, and then
c     the terms on the neighbours added to the accelerations in the order of the chunks
      use omp_utils
      use pairbuf_mod
      include 'starsmasher.h'
      include 'mpif.h'
      integer nloop,ic,it,nthr,ierr
      real*8 alpha1,beta1
      real*8 mvx(*),mvy(*),mvz(*),mgx(*),mgy(*),mgz(*)
      pb_fail=0
      nthr=max_threads()
      call pb_begin(nloop,n,nthr)
!$omp parallel default(shared) private(ic,it) num_threads(nthr)
      it=thread_id()+1
!$omp do schedule(dynamic,1)
      do ic=1,pb_nchunk
         call uvchunk(ic,nloop,alpha1,beta1,mvx,mvy,mvz,mgx,mgy,mgz,
     $        pb_scr(it)%smap,pb_scr(it)%bm,pb_scr(it)%acc)
      enddo
!$omp end do
!$omp end parallel
!     only the master thread may call MPI (MPI_THREAD_FUNNELED), so a thread that cannot go on says so
      if(pb_fail.ne.0) call mpi_abort(mpi_comm_world,1,ierr)
      call pb_merge(mvx,mvy,mvz,mgx,mgy,mgz)
      return
      end
