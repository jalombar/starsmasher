!***********************************************************************
      subroutine advance_block
!     Block (hierarchical, power-of-2) timesteps: McMillan (1986), Makino
!     (1991), Springel (2005, Sec. 5).  Selected with nblock=1 in sph.input.
!
!     One call advances the system by one substep, from the current tick to
!     the next tick at which some particle's step ends.  Particle i has step
!     dtmaxblk/2**kbin(i).  A step of bin k spans 2**(nbinmax-k) ticks of an
!     integer timeline, so deciding which particles are active is exact.
!
!     Kick-drift-kick per particle: vh/uh hold each particle's half-step
!     velocity and u (constant through its step).  Every substep drifts all
!     positions.  Forces are evaluated with velocities and u predicted to the
!     current time.  Active particles then receive the closing half-kick of
!     their old step and the opening half-kick of their new one.
!
!     Stage 1: forces are computed for all particles every substep (inactive
!     particles keep the accelerations from the start of their step), so this
!     checks the bookkeeping, not the speed.
      include 'starsmasher.h'
      include 'mpif.h'
      integer kbin(nmax)
      integer(kind=8) tstart(nmax),slen(nmax),tick,nticks,tnext,span,tend,snew
      real*8 vhx(nmax),vhy(nmax),vhz(nmax),uh(nmax),t0b,dtick
      logical started
      double precision divv(nmax)
      common/commdivv/divv
      real*8 displacex,displacey,displacez
      integer ndisplace
      common/displace/displacex,displacey,displacez,ndisplace
      save kbin,tstart,slen,tick,nticks,vhx,vhy,vhz,uh,t0b,dtick,started
      data started/.false./
      real*8 axs(nmax),ays(nmax),azs(nmax),uds(nmax)
      integer knew(nmax)
      logical active(nmax)
      integer i,j,in,k,knat,kmx,ierr,nact,nwake,nalist
      integer alist(nmax)
!     refresh marking: below this many active particles, the direct loop
!     over active particles is cheaper than building and searching a tree
!     (and the tree is never built for fewer than 8, which kdtree2 needs)
      integer nreftree
      parameter(nreftree=64)
      integer(kind=8) nwaketot,nsub,nacttot,nreftot
      save nwaketot,nsub,nacttot,nreftot
      data nwaketot,nsub,nacttot,nreftot/0,0,0,0/
      real*8 dtt,tmid,hstep_old,hstep_new,dtmin

      if(nrelax.ne.0 .or. ncooling.ne.0 .or. ndisplace.ne.0) then
         if(myrank.eq.0) write(69,*)&
              'nblock=1 needs nrelax=0, ncooling=0 and ndisplace=0'
         stop 'nblock=1 needs nrelax=0, ncooling=0, ndisplace=0'
      endif

      if(.not.started) then
         if(dtmaxblk.le.0.d0) dtmaxblk=dtout
         nticks=2_8**nbinmax
         dtick=dtmaxblk/dble(nticks)
!     shared-step convention: t runs half a step ahead of the positions
         t0b=t-0.5d0*dt
         t=t0b
         tick=0
         if(nintvar.ne.2) stop 'nblock=1 needs nintvar=2'
!     The velocities and u in hand are half-step values of the shared-step
!     scheme for the current global dt.  Take them back to time t.
         do i=1,ntot
            vx(i)=vx(i)-0.5d0*dt*vxdot(i)
            vy(i)=vy(i)-0.5d0*dt*vydot(i)
            vz(i)=vz(i)-0.5d0*dt*vzdot(i)
            if(u(i).ne.0.d0) u(i)=u(i)-0.5d0*dt*udot(i)
         enddo
         do i=1,ntot
            actblk(i)=.true.
            refblk(i)=.false.
         enddo
!     point particles: zero div v before the first h solve (a value read
!     from a dump written before div v was guarded may be NaN)
         do i=1,ntot
            if(u(i).eq.0.d0) divv(i)=0.d0
         enddo
         call build_worklists
         call rho_and_h
         if(ngr.ne.0) call gravforce
         call uvdots
         call tstep
         do i=1,ntot
            active(i)=.true.
            kbin(i)=0
         enddo
         call assign_bins(active,kbin,knew,tick,nticks)
         do i=1,ntot
            kbin(i)=knew(i)
            tstart(i)=0
            slen(i)=2_8**(nbinmax-kbin(i))
            hstep_new=0.5d0*dtmaxblk/dble(2_8**kbin(i))
            vhx(i)=vx(i)+hstep_new*vxdot(i)
            vhy(i)=vy(i)+hstep_new*vydot(i)
            vhz(i)=vz(i)+hstep_new*vzdot(i)
            if(u(i).ne.0.d0) uh(i)=u(i)+hstep_new*udot(i)
         enddo
         started=.true.
         if(myrank.eq.0) write(69,*)'block timesteps: dtmax=',dtmaxblk,&
              ' nbinmax=',nbinmax,' t0=',t0b
      endif

!     Wake-up (Saitoh & Makino 2009): a particle in the middle of a step that
!     is now more than twice its natural step (dtpart from the last tstep,
!     which includes the BH-gas pairwise limit) has its step cut short, to end
!     at the next tick aligned with the shorter step.  With step length L,
!     vh = v(t_start) + a L/2, so changing L to L' changes vh by a (L'-L)/2,
!     and the positions already drifted with the old vh are corrected to match.
      nwake=0
      do i=1,ntot
         if(dtpart(i).ge.0.5d0*dble(slen(i))*dtick) cycle
         k=0
         do while(dtmaxblk/dble(2_8**k).gt.dtpart(i) .and. k.lt.nbinmax)
            k=k+1
         enddo
         snew=2_8**(nbinmax-k)
         tend=(tick/snew+1)*snew
         if(tend-tstart(i).ge.slen(i)) cycle
         hstep_new=0.5d0*dble(tend-tstart(i)-slen(i))*dtick
         x(i)=x(i)+hstep_new*vxdot(i)*dble(tick-tstart(i))*dtick
         y(i)=y(i)+hstep_new*vydot(i)*dble(tick-tstart(i))*dtick
         z(i)=z(i)+hstep_new*vzdot(i)*dble(tick-tstart(i))*dtick
         vhx(i)=vhx(i)+hstep_new*vxdot(i)
         vhy(i)=vhy(i)+hstep_new*vydot(i)
         vhz(i)=vhz(i)+hstep_new*vzdot(i)
         if(u(i).ne.0.d0) uh(i)=uh(i)+hstep_new*udot(i)
         slen(i)=tend-tstart(i)
         kbin(i)=max(kbin(i),k)
         nwake=nwake+1
      enddo
      nwaketot=nwaketot+nwake

!     next tick at which some particle's step ends
      tnext=huge(tnext)
      do i=1,ntot
         span=slen(i)
         tnext=min(tnext,tstart(i)+span)
      enddo
      dtt=dble(tnext-tick)*dtick

!     drift all positions with their half-step velocities
      do i=1,ntot
         x(i)=x(i)+dtt*vhx(i)
         y(i)=y(i)+dtt*vhy(i)
         z(i)=z(i)+dtt*vhz(i)
      enddo
      tick=tnext
      t=t0b+dble(tick)*dtick
      dt=dtt                 ! rho_and_h predicts h with divv*dt/3

!     velocities and u predicted to time t for the force evaluation
      nact=0
      do i=1,ntot
         span=slen(i)
         active(i)=(tstart(i)+span.eq.tick)
         if(active(i)) nact=nact+1
         tmid=t0b+(dble(tstart(i))+0.5d0*dble(span))*dtick
         vx(i)=vhx(i)+(t-tmid)*vxdot(i)
         vy(i)=vhy(i)+(t-tmid)*vydot(i)
         vz(i)=vhz(i)+(t-tmid)*vzdot(i)
         if(u(i).ne.0.d0) u(i)=uh(i)+(t-tmid)*udot(i)
         axs(i)=vxdot(i)
         ays(i)=vydot(i)
         azs(i)=vzdot(i)
         uds(i)=udot(i)
      enddo

      do i=1,ntot
         actblk(i)=active(i)
      enddo
      if(nblockfull.eq.1) then
!     testing: h, rho and hydro recomputed for all particles (every
!     inactive particle counts as refreshed, so uvdots uses its list)
         do i=1,ntot
            refblk(i)=.not.active(i)
         enddo
         call build_worklists
         call rho_and_h
         if(ngr.ne.0) call gravforce
         call uvdots
      else
!     inactive gas particles: rho and h predicted from div v over the substep
         do i=1,ntot
            if(.not.active(i) .and. u(i).ne.0.d0) then
               rho(i)=rho(i)*exp(-divv(i)*dtt)
!     (h=htilde+hfloor, and it is htilde that follows the density)
               if(dynhco.eq.3) then
                  hdynco(i)=hdynco(i)*exp(divv(i)*dtt/3.d0)
                  hp(i)=(hdynco(i)**hcopnorm+hfloor**hcopnorm)**(1.d0/hcopnorm)
               else
                  hp(i)=hfloor+(hp(i)-hfloor)*exp(divv(i)*dtt/3.d0)
               endif
            endif
         enddo
!     inactive particles j with an active particle inside their kernel
!     (r < 2 h_j) re-solve h, rho, omega, chi, psi and div v / curl v at the
!     current positions: the active particles' j-side pair terms (and the
!     h_j half of the softened gravity) use exactly these.  Their
!     accelerations and kicks are unchanged
         do i=1,ntot
            refblk(i)=.false.
         enddo
         nalist=0
         do i=1,ntot
            if(active(i)) then
               nalist=nalist+1
               alist(nalist)=i
            endif
         enddo
         if(nalist.ge.max(nreftree,8)) then
!     many active particles: search a tree of the active ones instead
            call refresh_by_tree(nalist,alist)
         else
            do j=1,ntot
               if(active(j) .or. u(j).eq.0.d0) cycle
               do in=1,nalist
                  i=alist(in)
                  if((x(i)-x(j))**2+(y(i)-y(j))**2+(z(i)-z(j))**2&
                       .lt.4.d0*hp(j)**2) then
                     refblk(j)=.true.
                     exit
                  endif
               enddo
            enddo
         endif
         nacttot=nacttot+nact
         nreftot=nreftot+count(refblk(1:ntot))
         call build_worklists
         call rho_and_h
!     forces on the active particles (uvdots handles block steps through
!     actblk/refblk, and inactive particles' results are discarded below)
         if(ngr.ne.0) call gravforce
         call uvdots
      endif

!     inactive particles keep the accelerations of the start of their step,
!     and active ones get the closing half-kick with the new accelerations
      do i=1,ntot
         if(active(i)) then
            span=slen(i)
            hstep_old=0.5d0*dble(span)*dtick
            vx(i)=vhx(i)+hstep_old*vxdot(i)
            vy(i)=vhy(i)+hstep_old*vydot(i)
            vz(i)=vhz(i)+hstep_old*vzdot(i)
            if(u(i).ne.0.d0) u(i)=uh(i)+hstep_old*udot(i)
         else
            vxdot(i)=axs(i)
            vydot(i)=ays(i)
            vzdot(i)=azs(i)
            udot(i)=uds(i)
         endif
      enddo

!     energies need every particle's potential: only at full synchronization
      blksync=(mod(tick,nticks).eq.0)
      if(blksync) call enout(.true.)
      call tstep
      call assign_bins(active,kbin,knew,tick,nticks)

!     opening half-kick of the new step for active particles
      dtmin=1.d30
      do i=1,ntot
         if(active(i)) then
            kbin(i)=knew(i)
            tstart(i)=tick
            slen(i)=2_8**(nbinmax-kbin(i))
            hstep_new=0.5d0*dtmaxblk/dble(2_8**kbin(i))
            vhx(i)=vx(i)+hstep_new*vxdot(i)
            vhy(i)=vy(i)+hstep_new*vydot(i)
            vhz(i)=vz(i)+hstep_new*vzdot(i)
            if(u(i).ne.0.d0) uh(i)=u(i)+hstep_new*udot(i)
         endif
         dtmin=min(dtmin,dtmaxblk/dble(2_8**kbin(i)))
      enddo
!     Between calls, expose the state in the shared-step convention: v and u
!     predicted to t_positions+dt/2 for one dt, with t half a step ahead of
!     the positions.  A dump written now therefore restarts correctly in
!     either mode (the start-up above inverts exactly this).  The internal
!     half-step values vh,uh are untouched.
      dt=dtmin
      do i=1,ntot
         span=slen(i)
         tmid=t0b+(dble(tstart(i))+0.5d0*dble(span))*dtick
         vx(i)=vhx(i)+(t+0.5d0*dt-tmid)*vxdot(i)
         vy(i)=vhy(i)+(t+0.5d0*dt-tmid)*vydot(i)
         vz(i)=vhz(i)+(t+0.5d0*dt-tmid)*vzdot(i)
         if(u(i).ne.0.d0) u(i)=uh(i)+(t+0.5d0*dt-tmid)*udot(i)
      enddo
      t=t0b+dble(tick)*dtick+0.5d0*dt
      nsub=nsub+1
      if(myrank.eq.0 .and. mod(tick,nticks).eq.0) then
         write(69,*)'block sync t=',t,' finest bin=',maxval(kbin(1:ntot)),&
              ' substeps=',nsub,' wake-ups=',nwaketot,&
              ' mean active=',dble(nacttot)/max(nsub,1_8),&
              ' mean refreshed=',dble(nreftot)/max(nsub,1_8)
         write(69,'(a,40i6)')' bin histogram (bins 0..finest):',&
              (count(kbin(1:ntot).eq.k),k=0,maxval(kbin(1:ntot)))
         nsub=0
         nwaketot=0
         nacttot=0
         nreftot=0
      endif
      return
      end
!***********************************************************************
      subroutine assign_bins(active,kbin,knew,tick,nticks)
!     New bins for active particles from their natural steps dtpart (tstep).
!     Larger steps only when the current tick is a multiple of the new step.
!     No step more than 2**nblimit (default 4) times a neighbour's (Saitoh &
!     Makino 2009).
      include 'starsmasher.h'
      include 'mpif.h'
      logical active(nmax)
      integer kbin(nmax),knew(nmax)
      integer(kind=8) tick,nticks
      integer i,j,in,k,knat,ierr,iw,nloop
      do i=1,ntot
         knew(i)=-1
      enddo
!     block steps (blkdist): each active particle is handled by the rank that
!     holds its neighbour list (its round-robin share), otherwise by its owner
      if(blkdist) then
         nloop=nwmine
      else
         nloop=n_upper-n_lower+1
      endif
      do iw=1,nloop
         if(blkdist) then
            i=wmine(iw)
         else
            i=n_lower+iw-1
         endif
         if(.not.active(i)) cycle
         knat=0
         do while(dtmaxblk/dble(2_8**knat).gt.dtpart(i) .and. knat.lt.nbinmax)
            knat=knat+1
         enddo
         if(knat.lt.kbin(i)) then
            k=knat
            do while(mod(tick,2_8**(nbinmax-k)).ne.0)
               k=k+1
            enddo
            knat=k
         endif
         do in=1,nn(i)
            j=list(first(i)+in)
            knat=max(knat,kbin(j)-nblimit)
         enddo
         if(dtmaxblk/dble(2_8**nbinmax).gt.dtpart(i)) write(69,*)&
              'warning: particle',i,' wants dt=',dtpart(i),&
              ' below the finest block step, so raise nbinmax'
         knew(i)=min(knat,nbinmax)
      enddo
      call mpi_allreduce(mpi_in_place,knew,ntot,mpi_integer,mpi_max,&
           mpi_comm_world,ierr)
      do i=1,ntot
         if(.not.active(i)) knew(i)=kbin(i)
      enddo
      return
      end
!***********************************************************************
      subroutine get_gravity_active
!     Gravitational acceleration and potential on the active particles only
!     (actblk), from all particles: N_active x N pairs.  Same Wendland C4
!     softened pair formulas as get_gravity_using_cpus (nkernel=2), but each
!     pair contributes only to the active particle i (no reaction on j), so
!     inactive particles' gx,gy,gz,grpot are left at zero here.
      include 'starsmasher.h'
      integer i,j,kact
      real*8 dr1,dr2,dr3,r2,rinv,rinv2,amj,twohpi,twohpj,fourhpi2,fourhpj2
      real*8 mrinv1,mrinv3,invq1,invq2,qq1,qq2,q21,q22
      real*8 acc1,acc2,pot1,pot2,mj11,mj12,mj21,mj22,g1,g2,gacc,gpot
      if(nkernel.ne.2 .or. nselfgravity.ne.1) &
           stop 'nblock=1 gravity needs nkernel=2 and nselfgravity=1'
      do i=1,ntot
         gx(i)=0.d0
         gy(i)=0.d0
         gz(i)=0.d0
         grpot(i)=0.d0
      enddo
!     the active particles are dealt out round-robin over the gravity ranks,
!     so the work is even however they are numbered
      kact=0
      do i=1,ntot
         if(.not.actblk(i)) cycle
         kact=kact+1
         if(mod(kact-1,ngravprocs).ne.myrank) cycle
         twohpi=2*hp(i)
         fourhpi2=twohpi*twohpi
         do j=1,ntot
            dr1=x(j)-x(i)
            dr2=y(j)-y(i)
            dr3=z(j)-z(i)
            r2=dr1**2+dr2**2+dr3**2
            amj=am(j)
            twohpj=2*hp(j)
            fourhpj2=twohpj*twohpj
            if(j.ne.i .and. r2.ge.max(fourhpi2,fourhpj2)) then
               rinv=1/sqrt(r2)
               mrinv1=rinv*amj
               mrinv3=rinv*rinv*mrinv1
               gx(i)=gx(i)+mrinv3*dr1
               gy(i)=gy(i)+mrinv3*dr2
               gz(i)=gz(i)+mrinv3*dr3
               grpot(i)=grpot(i)-mrinv1
               cycle
            endif
            if(r2.gt.0) then
               rinv=1/sqrt(r2)
               invq1=rinv*twohpi
               invq2=rinv*twohpj
               qq1=1/invq1
               qq2=1/invq2
            else
               rinv=0.d0
               qq1=0.d0
               qq2=0.d0
            endif
            mrinv1=rinv*amj
            mrinv3=rinv*rinv*mrinv1
            q21=qq1*qq1
            q22=qq2*qq2
            acc1=20.625d0+q21*((-115.5d0)+q21*(618.75d0+qq1*((-1155.d0)&
                 +qq1*(962.5d0+qq1*((-396.d0)+65.625d0*qq1)))))
            acc2=20.625d0+q22*((-115.5d0)+q22*(618.75d0+qq2*((-1155.d0)&
                 +qq2*(962.5d0+qq2*((-396.d0)+65.625d0*qq2)))))
            mj11=amj/twohpi
            mj12=amj/twohpj
            mj21=mj11/fourhpi2
            mj22=mj12/fourhpj2
            if(r2.le.fourhpi2) then
               g1=1
               pot1=(-3.4375d0)+q21*(10.3125d0+q21*((-28.875d0)+q21*(103.125d0&
                    +qq1*((-165.d0)+qq1*(120.3125d0+qq1*((-44.d0)+6.5625d0*qq1))))))
            else
               g1=0
               pot1=0
            endif
            if(r2.le.fourhpj2) then
               g2=1
               pot2=(-3.4375d0)+q22*(10.3125d0+q22*((-28.875d0)+q22*(103.125d0&
                    +qq2*((-165.d0)+qq2*(120.3125d0+qq2*((-44.d0)+6.5625d0*qq2))))))
            else
               g2=0
               pot2=0
            endif
            if(r2.gt.0) then
               gacc=0.5d0*(g1*mj21*acc1+(1.0d0-g1)*mrinv3+&
                    g2*mj22*acc2+(1.0d0-g2)*mrinv3)
               gpot=0.5d0*(mj11*pot1+(g1-1.0d0)*mrinv1+&
                    mj12*pot2+(g2-1.0d0)*mrinv1)
            else
               gacc=0.d0
               gpot=0.d0
               if(j.eq.i .and. u(i).ne.0) gpot=-1.71875d0*(mj11+mj12) ! self term, counted once
            endif
            gx(i)=gx(i)+gacc*dr1
            gy(i)=gy(i)+gacc*dr2
            gz(i)=gz(i)+gacc*dr3
            grpot(i)=grpot(i)+gpot
         enddo
      enddo
      return
      end
!***********************************************************************
      subroutine sync_all(arr)
!     Every rank computes arr for its own particles (n_lower..n_upper).
!     Gather the full array onto every rank.
      include 'starsmasher.h'
      include 'mpif.h'
      real*8 arr(nmax)
      integer ierr
      if(nprocs.gt.1) call mpi_allgatherv(mpi_in_place,0,mpi_datatype_null,&
           arr,recvcounts,displs,mpi_double_precision,mpi_comm_world,ierr)
      return
      end
!***********************************************************************
      subroutine build_worklists
!     This substep's particles to be solved (active or refreshed), in index
!     order, identical on every rank, and this rank's round-robin share.
      include 'starsmasher.h'
      integer i
      nwall=0
      nwmine=0
      do i=1,ntot
         mywork(i)=.false.
         if(actblk(i) .or. refblk(i)) then
            nwall=nwall+1
            wall(nwall)=i
            if(mod(nwall-1,nprocs).eq.myrank) then
               nwmine=nwmine+1
               wmine(nwmine)=i
               mywork(i)=.true.
            endif
         endif
      enddo
      blkdist=.true.
      return
      end
!***********************************************************************
      subroutine sync_w(arr)
!     arr of the listed particles (wall) as computed by the rank each was
!     dealt to: the others contribute exact zeros to the sum.
      include 'starsmasher.h'
      include 'mpif.h'
      real*8 arr(nmax),buf(nwall)
      integer k,ierr
      if(nprocs.eq.1 .or. nwall.eq.0) return
      do k=1,nwall
         if(mywork(wall(k))) then
            buf(k)=arr(wall(k))
         else
            buf(k)=0.d0
         endif
      enddo
      call mpi_allreduce(mpi_in_place,buf,nwall,mpi_double_precision,&
           mpi_sum,mpi_comm_world,ierr)
      do k=1,nwall
         arr(wall(k))=buf(k)
      enddo
      return
      end
!***********************************************************************
      subroutine sync_wi(iarr)
!     integer version of sync_w
      include 'starsmasher.h'
      include 'mpif.h'
      integer iarr(nmax),ibuf(nwall)
      integer k,ierr
      if(nprocs.eq.1 .or. nwall.eq.0) return
      do k=1,nwall
         if(mywork(wall(k))) then
            ibuf(k)=iarr(wall(k))
         else
            ibuf(k)=0
         endif
      enddo
      call mpi_allreduce(mpi_in_place,ibuf,nwall,mpi_integer,&
           mpi_sum,mpi_comm_world,ierr)
      do k=1,nwall
         iarr(wall(k))=ibuf(k)
      enddo
      return
      end

!***********************************************************************
      subroutine refresh_by_tree(nalist,alist)
!     Mark for refreshing every inactive gas particle j with an active
!     particle inside its kernel (r < 2 h_j), as the direct loop in
!     advance_block does, but with a kd tree of the nalist active particles.
!     The tree is searched with a radius larger by a part in 1e12, and each
!     candidate is then tested with the direct loop's own expression, so
!     exactly the same particles are marked.
      use kdtree2_module
      include 'starsmasher.h'
      integer nalist,alist(nmax)
      real(kdkind), allocatable :: ad(:,:)
      real(kdkind), target :: qv(3)
      type(kdtree2), pointer :: atree
      type(kdtree2_result), allocatable :: res(:)
      integer j,k,i,nf
      allocate(ad(3,nalist),res(nalist))
      do k=1,nalist
         ad(1,k)=x(alist(k))
         ad(2,k)=y(alist(k))
         ad(3,k)=z(alist(k))
      enddo
      atree => kdtree2_create(ad,sort=.false.,rearrange=.true.)
      do j=1,ntot
         if(actblk(j) .or. u(j).eq.0.d0) cycle
         qv(1)=x(j)
         qv(2)=y(j)
         qv(3)=z(j)
         call kdtree2_r_nearest(tp=atree,qv=qv,r2=4.d0*hp(j)**2*(1.d0+1.d-12),&
              nfound=nf,nalloc=nalist,results=res)
         do k=1,nf
            i=alist(res(k)%idx)
            if((x(i)-x(j))**2+(y(i)-y(j))**2+(z(i)-z(j))**2&
                 .lt.4.d0*hp(j)**2) then
               refblk(j)=.true.
               exit
            endif
         enddo
      enddo
      call kdtree2_destroy(atree)
      deallocate(ad,res)
      return
      end
