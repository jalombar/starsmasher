      subroutine checkpt(itype)
c     writes checkpointing file.
c     if itype=0, write to externally defined unit 88.
c     if itype=1, write to local 'restart.sph' file if nit is
c     a multiple of nitch.
c     With block steps (nblock=1) nit counts substeps, and between full
c     synchronizations most particles hold velocities and u predicted
c     for the current substep.  So write instead at the first full
c     synchronization at least nitch substeps after the last checkpoint,
c     where the state is exactly what an out*.sph dump holds.
      include 'starsmasher.h'                                          
      include 'mpif.h'
      integer nitch,itype
      parameter (nitch=1000) ! dump a restartrad.sph checkpoint file every nitch iterations
      real*8 divv(nmax)
      common/commdivv/divv
      integer mylength, ierr
      integer nitlast
      save nitlast
      data nitlast/0/
      logical writeit

c      if (mod(nit,nitch).eq.0 .or. t.ge.32661.425d0) then
      if(nblock.eq.1) then
         writeit=blkfull .and. nit-nitlast.ge.nitch
      else
         writeit=mod(nit,nitch).eq.0
      endif
      if (writeit) then
         nitlast=nit

c     myrank=0 needs the rho, divv values to make the output file
         mylength=n_upper-n_lower+1
         if(myrank.ne.0)then
            call mpi_gatherv(rho(n_lower), mylength, mpi_double_precision,
     $           rho, recvcounts, displs, mpi_double_precision, 0,
     $           mpi_comm_world, ierr) 
            call mpi_gatherv(divv(n_lower), mylength, mpi_double_precision,
     $           divv, recvcounts, displs, mpi_double_precision, 0,
     $           mpi_comm_world, ierr) 
         else
            call mpi_gatherv(mpi_in_place, mylength, mpi_double_precision,
     $           rho, recvcounts, displs, mpi_double_precision, 0,
     $           mpi_comm_world, ierr) 
            call mpi_gatherv(mpi_in_place, mylength, mpi_double_precision,
     $           divv, recvcounts, displs, mpi_double_precision, 0,
     $           mpi_comm_world, ierr) 
         endif

         if(myrank.ne.0) return

         call wall_time(time1)
         write (69,*) 'checkpt: writing local checkpt file at nit=',
     $        nit
         open(12,file='restartrad.sph',form='unformatted',err=100)
         call dump(12)
         close (12)
         call wall_time(time2)
         seconds = time2-time1
         write (6,*) 'restartrad:',seconds,'s'
      endif

      return
c     error condition:
 100  stop 'checkpt: error opening unit ???'
      end
********************************************************************
      subroutine dump(iu)
c     writes binary dump file to unit iu.
c     this routine contains the standard format for all
c     binary dump files.
      include 'starsmasher.h'
      integer iu,i
      real*8 divv(nmax)
      common/commdivv/divv
      common /jumpcomm/ tjumpahead,rjump,rjumpfactor,rjumppercent,
     $     reentryfactor
      real*8 erad
      common/lostenergy/ erad
      real*8 epoteat,ekineat,einteat
      common /eatenenergy/ epoteat,ekineat,einteat
      integer lastnumstars
      common /starcount/ lastnumstars

c     The last four header entries were ndisplace and displace{x,y,z} until
c     the origin shift was removed.  They now hold what a restart needs and
c     did not have: changetf's count of stars, and the energy carried off by
c     particles the black hole swallowed (reat>0), which energy*.sph adds up.
      write(iu) ntot,nnopt,hco,hfloor,sep0,tf,dtout,nout,nit,t,
     $     nav,alpha,beta,tjumpahead,ngr,nrelax,trelax,dt,omega2,
     $     ncooling,erad,lastnumstars,epoteat,ekineat,einteat
      if(ncooling.eq.0) then
         do i=1,ntot
            write(iu) x(i),y(i),z(i),am(i),hp(i),rho(i),
     $           vx(i),vy(i),vz(i),vxdot(i),vydot(i),vzdot(i),
     $           u(i),udot(i),  !gx(i),gy(i),gz(i),
     $           grpot(i),meanmolecular(i),cc(i),divv(i)!,
!     $           por2(i)*rho(i)**2
         enddo
      else
         do i=1,ntot
            write(iu) x(i),y(i),z(i),am(i),hp(i),rho(i),
     $           vx(i),vy(i),vz(i),vxdot(i),vydot(i),vzdot(i),
     $           u(i),udot(i),  !gx(i),gy(i),gz(i),
     $           grpot(i),meanmolecular(i),cc(i),divv(i),
     $           ueq(i),tthermal(i)!,
!     $           por2(i)*rho(i)**2
         enddo
      endif
      write(iu) ntot            ! to check that the file was not corrupted

      return
      end
********************************************************************
      subroutine output
c     output at every iteration
      include 'starsmasher.h'
      
c     output every dtout:
      if (dtout.gt.0.d0) then
         if ((t.ge.dble(nout)*dtout).or.(t.ge.tf)) then
c            if(ngr.ne.0)then
c               call gravpot
c            endif
            call duout
            nout=nout+1
         endif
      else if (dtout.lt.0.d0) then
c     (if dtout<0, interpret as a number of iterations)
         if ((nit.ge.nout*int(abs(dtout))).or.(t.ge.tf)) then
            call duout
            nout=nout+1
         endif         
      endif
      
      return
      end
********************************************************************
      subroutine enout(printingenergytofile)
      use omp_utils
      include 'starsmasher.h'
      include 'mpif.h'
      real*8 xcm1,ycm1,zcm1,xcm2,ycm2,zcm2,am1,am2
      common/centersofmass/xcm1,ycm1,zcm1,xcm2,ycm2,zcm2,am1,am2
      real*8 hpi,ui,rhoi,
     $     etot,epot,ekin,eint,v2i,stot,ajtot,vcm,ajx,
     $     ajy,ajz,vcmx,vcmy,vcmz,xcm,ycm,zcm,amtot,xmax,
     $     ymax,zmax,xmin,ymin,zmin,pressi
      integer ihpmin,ihpmax,iumin,iumax,irhomin,irhomax,
     $     ipmin,ipmax
      logical printingenergytofile,logit,firstenout
      save firstenout
      data firstenout/.true./
      real*8 etotlast
      common/energyblock/ etotlast,etot
      real*8, save :: einit = +128.0;
      integer mynnemin(2),nnemin(2),mynnemax(2),nnemax(2)
      real*8 myrhomin(2),rhomin(2),myrhomax(2),rhomax(2)
      real*8 myumin(2),umin(2),myumax(2),umax(2)
      real*8 myhpmin(2),hpmin(2),myhpmax(2),hpmax(2)
      real*8 mypmin(2),pmin(2),mypmax(2),pmax(2)
      integer ierr,i
c     the sums of enout are made of partial sums over fixed blocks of
c     particles, added in block order (see omp_utils), so that their bits
c     do not depend on the number of threads
      integer ib,nb,nspecial,nown
      real*8 bs(10),bm(6)
      real*8 bsm(5),bmm(20)
      double precision mynneavr,nneavr,mynnesig,nnesig
      double precision myeint,mystot,myekin
      real*8 rhocgs, ucgs, temperature, ugas

c     variables used for radiative cooling portion of the code:
      real*8 erad
      common/lostenergy/ erad

      real*8 xi,yi,zi
      real*8 epoteat,ekineat,einteat
      common /eatenenergy/ epoteat,ekineat,einteat


      if(myrank.eq.0) then
c     find system box and center of mass:
         xmin=1.d30
         ymin=1.d30
         zmin=1.d30
         xmax=-1.d30
         ymax=-1.d30
         zmax=-1.d30
         amtot=0.d0
         xcm=0.d0
         ycm=0.d0
         zcm=0.d0
         vcmx=0.d0
         vcmy=0.d0
         vcmz=0.d0
         ajx=0.d0
         ajy=0.d0
         ajz=0.d0
c     omeg=sqrt(omega2)
         if(nrelax.eq.2) write(69,*)'enout: omeg=',omeg
         nb=nredblocks(ntot)
         call reserve_rpart(16,nb)
!$omp parallel do schedule(static) private(ib,i,xi,yi,zi,bs,bm)
         do ib=1,nb
            bm(1)=1.d30
            bm(2)=1.d30
            bm(3)=1.d30
            bm(4)=-1.d30
            bm(5)=-1.d30
            bm(6)=-1.d30
            bs=0.d0
            do i=(ib-1)*redblock+1,min(ib*redblock,ntot)
               bm(1)=min(x(i),bm(1))
               bm(2)=min(y(i),bm(2))
               bm(3)=min(z(i),bm(3))
               bm(4)=max(x(i),bm(4))
               bm(5)=max(y(i),bm(5))
               bm(6)=max(z(i),bm(6))
               bs(1)=bs(1)+am(i)
               bs(2)=bs(2)+am(i)*x(i)
               bs(3)=bs(3)+am(i)*y(i)
               bs(4)=bs(4)+am(i)*z(i)
               bs(5)=bs(5)+am(i)*vx(i)
               bs(6)=bs(6)+am(i)*vy(i)
               bs(7)=bs(7)+am(i)*vz(i)
               if((nrelax.eq.2 .and. .not. gonedynamic) .or. nrelax.eq.3) then 
                  bs(8)=bs(8)+am(i)*(y(i)*vz(i)-z(i)*(vy(i)+omeg*x(i)))
                  bs(9)=bs(9)+am(i)*(z(i)*(vx(i)-omeg*y(i))-x(i)*vz(i))
                  bs(10)=bs(10)+
     $                 am(i)*(x(i)*(vy(i)+omeg*x(i))-y(i)*(vx(i)-omeg*y(i)))
               else
                  xi=x(i)
                  yi=y(i)
                  zi=z(i)
                  bs(8)=bs(8)+am(i)*(yi*vz(i)-zi*vy(i))
                  bs(9)=bs(9)+am(i)*(zi*vx(i)-xi*vz(i))
                  bs(10)=bs(10)+am(i)*(xi*vy(i)-yi*vx(i))
               endif
            enddo
            rpart(1:10,ib)=bs
            rpart(11:16,ib)=bm
         enddo
         do ib=1,nb
            amtot=amtot+rpart(1,ib)
            xcm=xcm+rpart(2,ib)
            ycm=ycm+rpart(3,ib)
            zcm=zcm+rpart(4,ib)
            vcmx=vcmx+rpart(5,ib)
            vcmy=vcmy+rpart(6,ib)
            vcmz=vcmz+rpart(7,ib)
            ajx=ajx+rpart(8,ib)
            ajy=ajy+rpart(9,ib)
            ajz=ajz+rpart(10,ib)
            xmin=min(rpart(11,ib),xmin)
            ymin=min(rpart(12,ib),ymin)
            zmin=min(rpart(13,ib),zmin)
            xmax=max(rpart(14,ib),xmax)
            ymax=max(rpart(15,ib),ymax)
            zmax=max(rpart(16,ib),zmax)
         enddo
         xcm=xcm/amtot
         ycm=ycm/amtot
         zcm=zcm/amtot
         vcmx=vcmx/amtot
         vcmy=vcmy/amtot
         vcmz=vcmz/amtot
         vcm=sqrt(vcmx**2+vcmy**2+vcmz**2)
c     ajx=ajx/amtot
c     ajy=ajy/amtot
c     ajz=ajz/amtot
         ajtot=sqrt(ajx**2+ajy**2+ajz**2)
!         if(mod(nit,nitpot).eq.0)then
c     With block steps enout writes an energy*.sph row at every full
c     synchronization, which can be many per dtmaxblk.  The summary in the
c     log stays at one per dtmaxblk, plus the first call, at the start.
         logit=nblock.ne.1 .or. blkfull .or. .not.printingenergytofile
     $        .or. firstenout
         firstenout=.false.
         if(logit) then
            write (69,90) nit,t,
     $           xmin,xmax,ymin,ymax,zmin,zmax
!         else
!            write (69,91) nit,t,
!     $           xmin,xmax,ymin,ymax,zmin,zmax
!         endif
            write(6,*)'end it.',nit,'t=',t
         endif
 90      format(/,
     $        ' output: end of iteration ',i8,'       time=',f12.5,/,
     $        '   system box= ',g10.3,'< x <',g10.3,/,
     $        '               ',g10.3,'< y <',g10.3,/,
     $        '               ',g10.3,'< z <',g10.3)      
! 91      format(/,
!     $        ' output: end of iteration ',i8,'       time=',f12.5,/,
!     $        '   system box= ',g10.3,'< x <',g10.3,/,
!     $        '               ',g10.3,'< y <',g10.3,/,
!     $        '               ',g10.3,'< z <',g10.3)

         epot=0.d0
         if(ngr.ne.0)then
c     do loop to get total gravitational potential energy:
c     the sum is made by blocks.  The loop that follows it is the one that
c     reports point masses and NaNs, and it adds up epot itself in the
c     old way, so it runs only when there is something to report
            nb=nredblocks(ntot)
            call reserve_rpart(2,nb)
!$omp parallel do schedule(static) private(ib,i)
            do ib=1,nb
               rpart(1,ib)=0.d0
               rpart(2,ib)=0.d0
               do i=(ib-1)*redblock+1,min(ib*redblock,ntot)
                  rpart(1,ib)=rpart(1,ib)+am(i)*grpot(i)
                  if(grpot(i).ne.grpot(i)) rpart(2,ib)=rpart(2,ib)+1.d0
                  if(logit .and. (hp(i).le.0.d0 .or. u(i).eq.0.d0))
     $                 rpart(2,ib)=rpart(2,ib)+1.d0
               enddo
            enddo
            nspecial=0
            do ib=1,nb
               epot=epot+rpart(1,ib)
               nspecial=nspecial+nint(rpart(2,ib))
            enddo
            if(nspecial.gt.0) then
               epot=0.d0
               do i=1,ntot
                  if(logit .and. (hp(i).le.0.d0 .or. u(i).eq.0.d0)) then
                     write(69,'(a,i6,3g17.9)')
     $                    'position of point mass',i,x(i),y(i),z(i)
                  endif
                  epot=epot+am(i)*grpot(i)
                  if(grpot(i).ne.grpot(i)) then
                     write(129,*) epot,i,am(i),grpot(i),hp(i),x(i),y(i),z(i),u(i)
                  endif
               enddo
            endif
            if (epot.eq.0) then
               write(69,*) 'Problem calculating gravity: epot=0'
               stop
            endif
            epot=0.5d0*epot
         endif
      endif

c     get min/max values of various quantities:
      mynnemin(1)=nmax
      mynnemax(1)=0
      myrhomin(1)=1d30
      myrhomax(1)=0
      myumin(1)=1d30
      myumax(1)=0
      myhpmin(1)=1d30
      myhpmax(1)=0
      mypmin(1)=1d30
      mypmax(1)=0
      mynnemin(2)=0
      mynnemax(2)=0
      myrhomin(2)=0
      myrhomax(2)=0
      myumin(2)=0
      myumax(2)=0
      myhpmin(2)=0
      myhpmax(2)=0
      mypmin(2)=0
      mypmax(2)=0
      mynneavr=0
      mynnesig=0
      myeint=0.d0
      mystot=0.d0
      myekin=0.d0
      nown=n_upper-n_lower+1
      nb=nblk(nown,redblock_fine)
      call reserve_rpart(25,nb)
c     the equation of state inversion costs more in some blocks than in others (the core, the envelope): blocks are handed out
c     dynamically; each block has its own partial sum, so the bits do not depend on who does it
!$omp parallel do schedule(dynamic,4) private(ib,i,bsm,bmm,pressi,rhocgs,ucgs,temperature,ugas,v2i)
      do ib=1,nb
         bmm(1)=nmax
         bmm(2)=0.d0
         bmm(3)=0.d0
         bmm(4)=0.d0
         bmm(5)=1d30
         bmm(6)=0.d0
         bmm(7)=0.d0
         bmm(8)=0.d0
         bmm(9)=1d30
         bmm(10)=0.d0
         bmm(11)=0.d0
         bmm(12)=0.d0
         bmm(13)=1d30
         bmm(14)=0.d0
         bmm(15)=0.d0
         bmm(16)=0.d0
         bmm(17)=1d30
         bmm(18)=0.d0
         bmm(19)=0.d0
         bmm(20)=0.d0
         bsm=0.d0
         do i=n_lower+(ib-1)*redblock_fine,
     $        min(n_lower+ib*redblock_fine-1,n_upper)
         if(nn(i).lt.bmm(1))then
            bmm(1)=nn(i)
            bmm(2)=i
         endif
         if(nn(i).gt.bmm(3))then
            bmm(3)=nn(i)
            bmm(4)=i
         endif
         if(rho(i).lt.bmm(5) .and. u(i).ne.0.d0)then
            bmm(5)=rho(i)
            bmm(6)=i
         endif
         if(rho(i).gt.bmm(7))then
            bmm(7)=rho(i)
            bmm(8)=i
         endif
         if(u(i).lt.bmm(9) .and. u(i).ne.0.d0)then
            bmm(9)=u(i)
            bmm(10)=i
         endif
         if(u(i).gt.bmm(11))then
            bmm(11)=u(i)
            bmm(12)=i
         endif
         if(hp(i).lt.bmm(13) .and. u(i).ne.0)then
            bmm(13)=hp(i)
            bmm(14)=i
         endif
         if(hp(i).gt.bmm(15) .and. u(i).ne.0)then
            bmm(15)=hp(i)
            bmm(16)=i
         endif
         pressi = por2(i)*rho(i)**2
         if(pressi.lt.bmm(17) .and. u(i).ne.0.d0)then
            bmm(17)=pressi
            bmm(18)=i
         endif
         if(pressi.gt.bmm(19))then
            bmm(19)=pressi
            bmm(20)=i
         endif
         bsm(1)=bsm(1)+nn(i)
         bsm(2)=bsm(2)+nn(i)**2

         if(u(i).gt.0.d0) then
            if(nintvar.eq.1) then
c     p=(gam-1)*rho*u=a*rho^gam, so u=a*rho^(gam-1)/(gam-1)
               bsm(3)=bsm(3)+am(i)*u(i)*rho(i)**(gam-1.d0)/(gam-1.d0)
               ucgs=u(i)*rho(i)**(gam-1.d0)/(gam-1.d0)
     $              *gravconst*munit/runit
            else if(nintvar.eq.3) then
c     ln A is stored, so the internal energy has to come back through the
c     temperature before it can be added to the energy budget.
               rhocgs=rho(i)*munit/runit**3.d0
               call getT_from_lna(u(i),rhocgs,meanmolecular(i),-1.d0,
     $              temperature)
               ucgs=1.5d0*boltz*temperature/meanmolecular(i)
     $              +arad*temperature**4/rhocgs
               bsm(3)=bsm(3)+am(i)*ucgs/(gravconst*munit/runit)
            else
               bsm(3)=bsm(3)+am(i)*u(i)
               ucgs=u(i)*gravconst*munit/runit
            endif
            if(neos.eq.1) then
               rhocgs=rho(i)*munit/runit**3.d0
c     pressure_one found this temperature from these u, rho and mean molecular weight when pressure was called
c     (same arguments of gettemperature): then it is not found again
               if(nintvar.eq.tmemo_nv .and. tmemo_u(i).eq.u(i) .and.
     $              tmemo_rho(i).eq.rho(i) .and.
     $              tmemo_mu(i).eq.meanmolecular(i)) then
                  temperature=tempcache(i)
               else
                  call gettemperature(qconst*rhocgs/meanmolecular(i),
     $                 -ucgs*rhocgs/arad,
     $                 temperature )
               endif
               ugas=1.5d0*boltz*temperature/meanmolecular(i)
               bsm(4)=bsm(4)+am(i)*munit*(
     $              1.5d0*boltz/meanmolecular(i)*
     $              log(ugas*rhocgs**(-2.d0/3.d0)) +
     $              4.d0/3.d0*arad*temperature**3.d0/rhocgs)
            elseif(neos.eq.0) then
               if(nintvar.eq.1) then
                  bsm(4)=bsm(4)+am(i)*
     $                 log(u(i))/(gam-1.d0)
               else
                  bsm(4)=bsm(4)+am(i)*
     $                 log(u(i)*rho(i)**(1.d0-gam))/(gam-1.d0)
               endif
            elseif(neos.eq.2)then
c     would need some way of getting entropy for tabulated eos
c     if we can calculate entropy from a table, then fine.  but the entropy
c     calculation is not really needed.  it doesn't affect how particles are moved.
c     so updating the entropy calculation is a low priority.
            endif
         endif
         if((nrelax.eq.2 .and. .not. gonedynamic) .or. nrelax.eq.3)then 
            v2i=(vx(i)-omeg*y(i))**2+(vy(i)+omeg*x(i))**2+vz(i)**2
         else 
            v2i=vx(i)**2+vy(i)**2+vz(i)**2
         endif
         bsm(5)=bsm(5)+am(i)*v2i

         enddo
         rpart(1:20,ib)=bmm
         rpart(21:25,ib)=bsm
      enddo
c     the partial results of the blocks, combined in block order (the first
c     particle to reach a minimum or maximum, as the plain loop gave)
      do ib=1,nb
         if(rpart(1,ib).lt.mynnemin(1)) then
            mynnemin(1)=nint(rpart(1,ib))
            mynnemin(2)=nint(rpart(2,ib))
         endif
         if(rpart(3,ib).gt.mynnemax(1)) then
            mynnemax(1)=nint(rpart(3,ib))
            mynnemax(2)=nint(rpart(4,ib))
         endif
         if(rpart(5,ib).lt.myrhomin(1)) then
            myrhomin(1)=rpart(5,ib)
            myrhomin(2)=rpart(6,ib)
         endif
         if(rpart(7,ib).gt.myrhomax(1)) then
            myrhomax(1)=rpart(7,ib)
            myrhomax(2)=rpart(8,ib)
         endif
         if(rpart(9,ib).lt.myumin(1)) then
            myumin(1)=rpart(9,ib)
            myumin(2)=rpart(10,ib)
         endif
         if(rpart(11,ib).gt.myumax(1)) then
            myumax(1)=rpart(11,ib)
            myumax(2)=rpart(12,ib)
         endif
         if(rpart(13,ib).lt.myhpmin(1)) then
            myhpmin(1)=rpart(13,ib)
            myhpmin(2)=rpart(14,ib)
         endif
         if(rpart(15,ib).gt.myhpmax(1)) then
            myhpmax(1)=rpart(15,ib)
            myhpmax(2)=rpart(16,ib)
         endif
         if(rpart(17,ib).lt.mypmin(1)) then
            mypmin(1)=rpart(17,ib)
            mypmin(2)=rpart(18,ib)
         endif
         if(rpart(19,ib).gt.mypmax(1)) then
            mypmax(1)=rpart(19,ib)
            mypmax(2)=rpart(20,ib)
         endif
         mynneavr=mynneavr+rpart(21,ib)
         mynnesig=mynnesig+rpart(22,ib)
         myeint=myeint+rpart(23,ib)
         mystot=mystot+rpart(24,ib)
         myekin=myekin+rpart(25,ib)
      enddo

      call mpi_reduce(mynnemin,nnemin,1,mpi_2integer,mpi_minloc,0,
     $     mpi_comm_world,ierr)
      call mpi_reduce(mynnemax,nnemax,1,mpi_2integer,mpi_maxloc,0,
     $     mpi_comm_world,ierr)
      call mpi_reduce(myrhomin,rhomin,1,mpi_2double_precision,
     $     mpi_minloc,0,mpi_comm_world,ierr)
      call mpi_reduce(myrhomax,rhomax,1,mpi_2double_precision,
     $     mpi_maxloc,0,mpi_comm_world,ierr)
      call mpi_reduce(myumin,umin,1,mpi_2double_precision,
     $     mpi_minloc,0,mpi_comm_world,ierr)
      call mpi_reduce(myumax,umax,1,mpi_2double_precision,
     $     mpi_maxloc,0,mpi_comm_world,ierr)
      call mpi_reduce(myhpmin,hpmin,1,mpi_2double_precision,
     $     mpi_minloc,0,mpi_comm_world,ierr)
      call mpi_reduce(myhpmax,hpmax,1,mpi_2double_precision,
     $     mpi_maxloc,0,mpi_comm_world,ierr)
      call mpi_reduce(mypmin,pmin,1,mpi_2double_precision,
     $     mpi_minloc,0,mpi_comm_world,ierr)
      call mpi_reduce(mypmax,pmax,1,mpi_2double_precision,
     $     mpi_maxloc,0,mpi_comm_world,ierr)
      call mpi_reduce(mynneavr,nneavr,1,mpi_double_precision,mpi_sum,0,
     $     mpi_comm_world,ierr)
      call mpi_reduce(mynnesig,nnesig,1,mpi_double_precision,mpi_sum,0,
     $     mpi_comm_world,ierr)
      call mpi_reduce(myeint,eint,1,mpi_double_precision,mpi_sum,0,
     $     mpi_comm_world,ierr)
      call mpi_reduce(myekin,ekin,1,mpi_double_precision,mpi_sum,0,
     $     mpi_comm_world,ierr)
      ekin=0.5d0*ekin
      call mpi_reduce(mystot,stot,1,mpi_double_precision,mpi_sum,0,
     $     mpi_comm_world,ierr)

      if(myrank.eq.0) then

         etot=epot+ekin+eint
         
!         if(mod(nit,nitpot).eq.0)then
            if(logit) write (69,904) epot,ekin,eint,etot,stot,vcm,ajtot
 904        format('   energies: W=',g10.3,' T=',f10.4,' U=',f10.4,/,
     $           '     Etot=',f10.4,' Stot=',g10.3,' vcm=',g10.3,
     $           ' Jtot=',g10.3)
!         else
!            write (69,905) epot,ekin,eint,etot,stot,vcm,ajtot
! 905        format('   energies:W''=',g10.3,' T=',f10.4,' U=',f10.4,/,
!     $           '    Etot''=',f10.4,' Stot=',g10.3,' vcm=',g10.3,
!     $           ' Jtot=',g10.3)
!         endif

         nneavr=nneavr/n
         nnesig=sqrt(nnesig/n-nneavr**2)
         if(logit) write (69,92) nneavr,nnesig,nnemin(1),nnemax(1),
     $        nnemin(2),nnemax(2)
 92      format('   neighbors: avr=',f4.0,' sig=',f4.0,'  min=',i4,
     $        ' max=',i5,i6,i6)
         
 906     format('   density: rhomin=',g10.3,' at (',g10.3,',',g10.3,',',
     $        g10.3,')',i6,/,
     $        '            rhomax=',g10.3,' at (',g10.3,',',g10.3,',',
     $        g10.3,')',i6,/,
     $        '    pressure: pmin=',g10.3,' at (',g10.3,',',g10.3,',',
     $        g10.3,')',i6,/,   
     $        '              pmax=',g10.3,' at (',g10.3,',',g10.3,',',
     $        g10.3,')',i6)      
         
         irhomin=nint(rhomin(2))
         irhomax=nint(rhomax(2))
         iumin=nint(umin(2))
         iumax=nint(umax(2))
         ihpmin=nint(hpmin(2))
         ihpmax=nint(hpmax(2))
         ipmin=nint(pmin(2))
         ipmax=nint(pmax(2))
         
         if(logit) write(69,906) rhomin(1),x(irhomin),y(irhomin),
     $        z(irhomin),irhomin,
     $        rhomax(1),x(irhomax),y(irhomax),z(irhomax),irhomax,
     $        pmin(1),x(ipmin),y(ipmin),z(ipmin),ipmin,
     $        pmax(1),x(ipmax),y(ipmax),z(ipmax),ipmax
         

         if(logit .and. nintvar.eq.1)then
            write(69,908) umin(1),x(iumin),y(iumin),z(iumin),iumin,
     $           umax(1),x(iumax),y(iumax),z(iumax),iumax,
     $           hpmin(1),x(ihpmin),y(ihpmin),z(ihpmin),ihpmin,
     $           hpmax(1),x(ihpmax),y(ihpmax),z(ihpmax),ihpmax !,
         else if(logit) then
            write(69,909) umin(1),x(iumin),y(iumin),z(iumin),iumin,
     $           umax(1),x(iumax),y(iumax),z(iumax),iumax,
     $           hpmin(1),x(ihpmin),y(ihpmin),z(ihpmin),ihpmin,
     $           hpmax(1),x(ihpmax),y(ihpmax),z(ihpmax),ihpmax !,
         endif
 908     format(
     $        '  in entropy: amin=',g10.3,' at (',g10.3,',',g10.3,',',
     $        g10.3,')',i6,/,
     $        '              amax=',g10.3,' at (',g10.3,',',g10.3,',',
     $        g10.3,')',i6,/,
     $        '   smoothing: hmin=',g10.3,' at (',g10.3,',',g10.3,',',
     $        g10.3,')',i6,/,
     $        '              hmax=',g10.3,' at (',g10.3,',',g10.3,',',
     $        g10.3,')',i6)
         
 909     format(
     $        '   in energy: umin=',g10.3,' at (',g10.3,',',g10.3,',',
     $        g10.3,')',i6,/,
     $        '              umax=',g10.3,' at (',g10.3,',',g10.3,',',
     $        g10.3,')',i6,/,
     $        '   smoothing: hmin=',g10.3,' at (',g10.3,',',g10.3,',',
     $        g10.3,')',i6,/,
     $        '              hmax=',g10.3,' at (',g10.3,',',g10.3,',',
     $        g10.3,')',i6)
         
c     append results of iteration to file:
         if(printingenergytofile)then
c            if (einit == 128.0) einit = etot
c            if (einit .eq. 0.0) then
c               print *, etot
c               stop 'error: einit = 0'
c            end if
            if(nrelax.ge.2)then
               call getcoms
               write(69,*)'separation=',
     $              sqrt((xcm1-xcm2)**2+(ycm1-ycm2)**2+(zcm1-zcm2)**2)
c               write(22,'(15g15.7)')t,epot,ekin,eint,etot,
c     $              (etot-einit)/abs(einit), stot,ajtot,
               if(ncooling.eq.0) then
                  write(22,'(14g15.7)')t,epot,ekin,eint,etot,stot,ajtot,
     $                 omeg,
     $                 sqrt((xcm1-xcm2)**2+(ycm1-ycm2)**2+(zcm1-zcm2)**2)
               else
                  write(22,'(14g15.7)')t,epot,ekin,eint,etot,stot,ajtot,
     $                 erad,omeg,
     $                 sqrt((xcm1-xcm2)**2+(ycm1-ycm2)**2+(zcm1-zcm2)**2)
               endif
            else
c               write(22,'(15g15.7)')t,epot,ekin,eint,etot,
c     $              (etot-einit)/abs(einit), stot,ajtot
               if(ncooling.eq.0) then
                  if(reat.le.0d0) then
                     write(22,'(14g15.7)')t,epot,ekin,eint,etot,stot,ajtot
                  else
                     write(22,'(14g15.7)')t,epot,ekin,eint,etot,stot,ajtot,
     $                    epoteat,ekineat,einteat
                  endif
               else
                  if(reat.le.0) then
                     write(22,'(14g15.7)')t,epot,ekin,eint,etot,stot,ajtot,
     $                    erad
                  else
                     write(22,'(14g15.7)')t,epot,ekin,eint,etot,stot,ajtot,
     $                    erad,epoteat,ekineat,einteat
                  endif
               endif
            endif
            call flush(22)
         else
            write(69,'(14g15.7)')'t','epot','ekin','eint','etot','stot',
     $           'ajtot'
            write(69,'(14g15.7)')t,epot,ekin,eint,etot,stot,ajtot
            write(34,'(14g15.7)')t,epot,ekin,eint,etot,stot,ajtot
         endif
         
      endif

      return
      end
********************************************************************
      subroutine duout
c     write binary dump file containing complete current results
      include 'starsmasher.h'
      include 'mpif.h' 
      character*16 outfn
      real*8 rhocgs,ucgs,temperature,rtemp
      integer i
      character*3 iname
      common /inittcom/ iname
      logical autotf,autojump
      common/autotfblock/autotf,autojump
      real*8 divv(nmax)
      common/commdivv/divv
      integer mylength, ierr, mygravlength
      integer comm_worker
      common/gravworkers/comm_worker
      real*8 gtot(nmax)
      real*8 tfsave
      real*8 tfcap,tfrestart
      integer ndtoutafter
      common/destroyedcom/tfcap,tfrestart,ndtoutafter

c     myrank=0 needs the rho,divv values to make the output file
      mylength=n_upper-n_lower+1
      if(myrank.ne.0)then
         call mpi_gatherv(rho(n_lower), mylength, mpi_double_precision,
     $        rho, recvcounts, displs, mpi_double_precision, 0,
     $        mpi_comm_world, ierr)
         call mpi_gatherv(divv(n_lower), mylength, mpi_double_precision,
     $        divv, recvcounts, displs, mpi_double_precision, 0,
     $        mpi_comm_world, ierr)
      else
         call mpi_gatherv(mpi_in_place, mylength, mpi_double_precision,
     $        rho, recvcounts, displs, mpi_double_precision, 0,
     $        mpi_comm_world, ierr)
         call mpi_gatherv(mpi_in_place, mylength, mpi_double_precision,
     $        divv, recvcounts, displs, mpi_double_precision, 0,
     $        mpi_comm_world, ierr)
      endif

      if(nrelax.eq.1 .and. mod(nout,10).eq.0) then
c     this is here for the generation of col*.sph files
         if(myrank.ne.0)then
            call mpi_gatherv(por2(n_lower), mylength, mpi_double_precision,
     $           por2, recvcounts, displs, mpi_double_precision, 0,
     $           mpi_comm_world, ierr)
         else
            call mpi_gatherv(mpi_in_place, mylength, mpi_double_precision,
     $           por2, recvcounts, displs, mpi_double_precision, 0,
     $           mpi_comm_world, ierr)
         endif
      endif

      if(myrank.eq.0) then

         if(nout.lt.10000) then
            write(outfn,1101) nout
 1101       format('out',i4.4,'.sph')
         else if(nout.lt.100000) then
            write(outfn,2101) nout
 2101       format('out',i5.5,'.sph')
         else if(nout.lt.1000000) then
            write(outfn,3101) nout
 3101       format('out',i6.6,'.sph')
         else
            write(outfn,4101) nout
 4101       format('out',i7.7,'.sph')
         endif
         
         write (69,*) 'duout: writing file ',outfn,'at t=',t
         
         open(12,file=outfn,form='unformatted')
         call dump(12)
         close (12)
      endif
   
      if(nrelax.eq.1 .and. mod(nout,10).eq.0) then
         mylength=n_upper-n_lower+1
         if(myrank.ne.0)then
            call mpi_gatherv(nn(n_lower), mylength, mpi_integer,
     $           nn, recvcounts, displs, mpi_integer, 0,
     $           mpi_comm_world, ierr)
         else
            call mpi_gatherv(mpi_in_place, mylength, mpi_integer,
     $           nn, recvcounts, displs, mpi_integer, 0,
     $           mpi_comm_world, ierr)
         endif
         if(ngr.ne.0 .and. myrank.lt.ngravprocs)then
c     myrank=0 needs the nn,gx,gy,gz values to make the col file
            if(ngravprocs.gt.1) then
               if(nusegpus.eq.1) then
c     gpu gravity: each rank holds the full g for its own slice
                  mygravlength=ngrav_upper-ngrav_lower+1
                  if(myrank.ne.0)then
                     call mpi_gatherv(gx(ngrav_lower), mygravlength, mpi_double_precision,
     $                    gx, gravrecvcounts, gravdispls, mpi_double_precision, 0,
     $                    comm_worker, ierr)
                     call mpi_gatherv(gy(ngrav_lower), mygravlength, mpi_double_precision,
     $                    gy, gravrecvcounts, gravdispls, mpi_double_precision, 0,
     $                    comm_worker, ierr)
                     call mpi_gatherv(gz(ngrav_lower), mygravlength, mpi_double_precision,
     $                    gz, gravrecvcounts, gravdispls, mpi_double_precision, 0,
     $                    comm_worker, ierr)
                  else
                     call mpi_gatherv(mpi_in_place, mygravlength, mpi_double_precision,
     $                    gx, gravrecvcounts, gravdispls, mpi_double_precision, 0,
     $                    comm_worker, ierr)
                     call mpi_gatherv(mpi_in_place, mygravlength, mpi_double_precision,
     $                    gy, gravrecvcounts, gravdispls, mpi_double_precision, 0,
     $                    comm_worker, ierr)
                     call mpi_gatherv(mpi_in_place, mygravlength, mpi_double_precision,
     $                    gz, gravrecvcounts, gravdispls, mpi_double_precision, 0,
     $                    comm_worker, ierr)
                  endif
               else
c     cpu_grav leaves each rank with a partial sum over all particles,
c     so the partial sums must be added, not gathered slice by slice
                  call mpi_reduce(gx, gtot, n, mpi_double_precision,
     $                 mpi_sum, 0, comm_worker, ierr)
                  if(myrank.eq.0) gx(1:n)=gtot(1:n)
                  call mpi_reduce(gy, gtot, n, mpi_double_precision,
     $                 mpi_sum, 0, comm_worker, ierr)
                  if(myrank.eq.0) gy(1:n)=gtot(1:n)
                  call mpi_reduce(gz, gtot, n, mpi_double_precision,
     $                 mpi_sum, 0, comm_worker, ierr)
                  if(myrank.eq.0) gz(1:n)=gtot(1:n)
               endif
            endif
         endif

         if(myrank.eq.0) then
            write(outfn,102) nout
 102        format('col',i4.4,'.sph')
            open(13,file=outfn,status='unknown')
            do i=1,n
               rhocgs=rho(i)*munit/runit**3.d0
               if(nintvar.eq.1) then
                  ucgs=u(i)*rho(i)**(gam-1.d0)/(gam-1.d0)
     $                 *gravconst*munit/runit
               else if(nintvar.eq.3 .and. u(i).ne.0.d0) then
                  call getT_from_lna(u(i),rhocgs,meanmolecular(i),-1.d0,
     $                 temperature)
                  ucgs=1.5d0*boltz*temperature/meanmolecular(i)
     $                 +arad*temperature**4/rhocgs
               else
                  ucgs=u(i)*gravconst*munit/runit
               endif
               if(nintvar.eq.3 .and. u(i).ne.0.d0) then
c     the temperature came out of the ln A inversion just above
               else if(arad.gt.0.d0.and.nintvar.eq.2
     $                 .and.u(i).ne.0.d0)then
                  call gettemperature(qconst*rhocgs/meanmolecular(i),
     $                 -ucgs*rhocgs/arad,
     $                 temperature )
               else
                  temperature=ucgs/(1.5d0*boltz)*meanmolecular(i)
               endif
               rtemp = (x(i)**2+y(i)**2+z(i)**2)**0.5d0
               write(13,'(16g15.7)')rtemp,por2(i)*rho(i)**2,rho(i),
     $              temperature,meanmolecular(i),am(i),hp(i),nn(i),
     $              (gx(i)*x(i)+gy(i)*y(i)+gz(i)*z(i))/rtemp,
     $              ((vxdot(i)-gx(i))*x(i)+(vydot(i)-gy(i))*y(i)+
     $              (vzdot(i)-gz(i))*z(i))/rtemp,x(i),y(i),z(i),grpot(i),
     $              u(i),vx(i)**2+vy(i)**2+vz(i)**2
            enddo
            close(13)
         endif
      endif

c     changetf is where rjump is tested, so it runs whenever jumping is on.  It
c     may only change tf, though, when the user asked for that with a negative
c     tf; otherwise the tf from sph.input stands.
      if(autotf .or. autojump) then
         tfsave=tf
         if(myrank.eq.0)
     $        write(69,*)'*************calling changetf***************',t
         call changetf
         if(myrank.eq.0)
     $        write(69,*)'*************done w/ changetf***************',t
         if(.not.autotf) tf=tfsave
c     the end changetf sets once a jump-ahead run's star is destroyed holds
c     even when tf is otherwise the user's
         tf=min(tf,tfcap)
      endif
      
      return
      end
