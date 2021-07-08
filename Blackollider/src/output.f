      subroutine checkpt(itype)
c     writes checkpointing file.
c     if itype=0, write to externally defined unit 88.
c     if itype=1, write to local 'restart.sph' file if nit is
c     a multiple of nitch.
      include 'starsmasher.h'                                          
      include 'mpif.h'
      integer nitch,itype
      parameter (nitch=1000) ! dump a restartrad.sph checkpoint file every nitch iterations
      real*8 divv(nmax)
      common/commdivv/divv
      integer mylength, ierr

c      if (mod(nit,nitch).eq.0 .or. t.ge.32661.425d0) then
      if (mod(nit,nitch).eq.0) then

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

         call cpu_time(time1)
         write (69,*) 'checkpt: writing local checkpt file at nit=',
     $        nit
         open(12,file='restartrad.sph',form='unformatted',err=100)
         call dump(12)
         close (12)
         call cpu_time(time2)
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
      common /jumpcomm/ tjumpahead
      real*8 erad
      common/lostenergy/ erad
      real*8 displacex, displacey,displacez
      integer ndisplace
      common/displace/displacex,displacey,displacez,ndisplace

      write(iu) ntot,nnopt,hco,hfloor,sep0,tf,dtout,nout,nit,t,
     $     nav,alpha,beta,tjumpahead,ngr,nrelax,trelax,dt,omega2,
     $     ncooling,erad,ndisplace,displacex,displacey,displacez
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
      logical printingenergytofile
      real*8 etotlast
      common/energyblock/ etotlast,etot
      real*8, save :: einit = +128.0;
      integer mynnemin(2),nnemin(2),mynnemax(2),nnemax(2)
      real*8 myrhomin(2),rhomin(2),myrhomax(2),rhomax(2)
      real*8 myumin(2),umin(2),myumax(2),umax(2)
      real*8 myhpmin(2),hpmin(2),myhpmax(2),hpmax(2)
      real*8 mypmin(2),pmin(2),mypmax(2),pmax(2)
      integer ierr,i
      double precision mynneavr,nneavr,mynnesig,nnesig
      double precision myeint,mystot,myekin
      real*8 rhocgs, ucgs, temperature, ugas

c     variables used for radiative cooling portion of the code:
      real*8 erad
      common/lostenergy/ erad
      real*8 xi,yi,zi,displacex,displacey,displacez
      integer ndisplace
      common/displace/displacex,displacey,displacez,ndisplace

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
         do i=1,ntot
            xmin=min(x(i),xmin)
            ymin=min(y(i),ymin)
            zmin=min(z(i),zmin)
            xmax=max(x(i),xmax)
            ymax=max(y(i),ymax)
            zmax=max(z(i),zmax)
            amtot=amtot+am(i)
            xcm=xcm+am(i)*x(i)
            ycm=ycm+am(i)*y(i)
            zcm=zcm+am(i)*z(i)
            vcmx=vcmx+am(i)*vx(i)
            vcmy=vcmy+am(i)*vy(i)
            vcmz=vcmz+am(i)*vz(i)
            if((nrelax.eq.2 .and. .not. gonedynamic) .or. nrelax.eq.3) then 
               ajx=ajx+am(i)*(y(i)*vz(i)-z(i)*(vy(i)+omeg*x(i)))
               ajy=ajy+am(i)*(z(i)*(vx(i)-omeg*y(i))-x(i)*vz(i))
               ajz=ajz+
     $              am(i)*(x(i)*(vy(i)+omeg*x(i))-y(i)*(vx(i)-omeg*y(i)))
            else
               xi=x(i)+displacex
               yi=y(i)+displacey
               zi=z(i)+displacez
               ajx=ajx+am(i)*(yi*vz(i)-zi*vy(i))
               ajy=ajy+am(i)*(zi*vx(i)-xi*vz(i))
               ajz=ajz+am(i)*(xi*vy(i)-yi*vx(i))
            endif
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
            write (69,90) nit,t,
     $           xmin,xmax,ymin,ymax,zmin,zmax
!         else
!            write (69,91) nit,t,
!     $           xmin,xmax,ymin,ymax,zmin,zmax
!         endif
         write(6,*)'end it.',nit,'t=',t
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
            do i=1,ntot
               if(hp(i).le.0.d0 .or. u(i).eq.0.d0) then
                  write(69,'(a,i6,3g17.9)')
     $                 'position of point mass',i,x(i),y(i),z(i)
               endif
c             if(u(i).ne.0.d0)then
               epot=epot+am(i)*grpot(i)

               if(grpot(i).ne.grpot(i)) then
                  write(129,*) epot,i,am(i),grpot(i),hp(i),x(i),y(i),z(i),u(i)
               endif

c             else
c              epot=epot+am(i)**2*42.d0/30.d0/hp(i) ! undo self-energy of point mass
c             endif
            enddo
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
      mynneavr=0
      mynnesig=0
      myeint=0.d0
      mystot=0.d0
      myekin=0.d0
      do i=n_lower,n_upper
         if(nn(i).lt.mynnemin(1))then
            mynnemin(1)=nn(i)
            mynnemin(2)=i
         endif
         if(nn(i).gt.mynnemax(1))then
            mynnemax(1)=nn(i)
            mynnemax(2)=i
         endif
         if(rho(i).lt.myrhomin(1) .and. u(i).ne.0.d0)then
            myrhomin(1)=rho(i)
            myrhomin(2)=i
         endif
         if(rho(i).gt.myrhomax(1))then
            myrhomax(1)=rho(i)
            myrhomax(2)=i
         endif
         if(u(i).lt.myumin(1) .and. u(i).ne.0.d0)then
            myumin(1)=u(i)
            myumin(2)=i
         endif
         if(u(i).gt.myumax(1))then
            myumax(1)=u(i)
            myumax(2)=i
         endif
         if(hp(i).lt.myhpmin(1) .and. u(i).ne.0)then
            myhpmin(1)=hp(i)
            myhpmin(2)=i
         endif
         if(hp(i).gt.myhpmax(1) .and. u(i).ne.0)then
            myhpmax(1)=hp(i)
            myhpmax(2)=i
         endif
         pressi = por2(i)*rho(i)**2
         if(pressi.lt.mypmin(1) .and. u(i).ne.0.d0)then
            mypmin(1)=pressi
            mypmin(2)=i
         endif
         if(pressi.gt.mypmax(1))then
            mypmax(1)=pressi
            mypmax(2)=i
         endif
         mynneavr=mynneavr+nn(i)
         mynnesig=mynnesig+nn(i)**2

         if(u(i).gt.0.d0) then
            if(nintvar.eq.1) then
c     p=(gam-1)*rho*u=a*rho^gam, so u=a*rho^(gam-1)/(gam-1)
               myeint=myeint+am(i)*u(i)*rho(i)**(gam-1.d0)/(gam-1.d0)
               ucgs=u(i)*rho(i)**(gam-1.d0)/(gam-1.d0)
     $              *gravconst*munit/runit
            else
               myeint=myeint+am(i)*u(i)
               ucgs=u(i)*gravconst*munit/runit
            endif
            if(neos.eq.1) then
               rhocgs=rho(i)*munit/runit**3.d0
               call gettemperature(qconst*rhocgs/meanmolecular(i),
     $              -ucgs*rhocgs/arad,
     $              temperature )
               ugas=1.5d0*boltz*temperature/meanmolecular(i)
               mystot=mystot+am(i)*munit*(
     $              1.5d0*boltz/meanmolecular(i)*
     $              log(ugas*rhocgs**(-2.d0/3.d0)) +
     $              4.d0/3.d0*arad*temperature**3.d0/rhocgs)
            elseif(neos.eq.0) then
               if(nintvar.eq.1) then
                  mystot=mystot+am(i)*
     $                 log(u(i))/(gam-1.d0)
               else
                  mystot=mystot+am(i)*
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
         myekin=myekin+am(i)*v2i

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
            write (69,904) epot,ekin,eint,etot,stot,vcm,ajtot
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
         write (69,92) nneavr,nnesig,nnemin(1),nnemax(1),nnemin(2),
     $        nnemax(2)
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
         
         write(69,906) rhomin(1),x(irhomin),y(irhomin),z(irhomin),irhomin,
     $        rhomax(1),x(irhomax),y(irhomax),z(irhomax),irhomax,
     $        pmin(1),x(ipmin),y(ipmin),z(ipmin),ipmin,
     $        pmax(1),x(ipmax),y(ipmax),z(ipmax),ipmax
         

         if(nintvar.eq.1)then
            write(69,908) umin(1),x(iumin),y(iumin),z(iumin),iumin,
     $           umax(1),x(iumax),y(iumax),z(iumax),iumax,
     $           hpmin(1),x(ihpmin),y(ihpmin),z(ihpmin),ihpmin,
     $           hpmax(1),x(ihpmax),y(ihpmax),z(ihpmax),ihpmax !,
         else
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
            if (einit == 128.0) einit = etot
            if (einit .eq. 0.0) then
               print *, etot
               stop 'error: einit = 0'
            end if
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
                  write(22,'(14g15.7)')t,epot,ekin,eint,etot,stot,ajtot
               else
                  write(22,'(14g15.7)')t,epot,ekin,eint,etot,stot,ajtot,
     $              erad
               endif
            endif
            call flush(22)
         else
            write(69,'(14g15.7)')'t','epot','ekin','eint','etot','stot',
     $           'ajtot'
            write(69,'(14g15.7)')t,epot,ekin,eint,etot,stot,ajtot
            write(34,'(14g15.7)')t,epot,ekin,eint,etot,stot,ajtot
            call flush(34)
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
      logical autotf
      common/autotfblock/autotf
      real*8 divv(nmax)
      common/commdivv/divv
      integer mylength, ierr, mygravlength
      integer comm_worker
      common/gravworkers/comm_worker

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
         else
            write(outfn,2101) nout
 2101       format('out',i5.5,'.sph')
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
               mygravlength=ngrav_upper-ngrav_lower+1
               if(myrank.ne.0)then
                  call mpi_gatherv(gx(ngrav_lower), mygravlength, mpi_double_precision,
     $                 gx, gravrecvcounts, gravdispls, mpi_double_precision, 0,
     $                 comm_worker, ierr)
                  call mpi_gatherv(gy(ngrav_lower), mygravlength, mpi_double_precision,
     $                 gy, gravrecvcounts, gravdispls, mpi_double_precision, 0,
     $                 comm_worker, ierr)
                  call mpi_gatherv(gz(ngrav_lower), mygravlength, mpi_double_precision,
     $                 gz, gravrecvcounts, gravdispls, mpi_double_precision, 0,
     $                 comm_worker, ierr)
               else
                  call mpi_gatherv(mpi_in_place, mygravlength, mpi_double_precision,
     $                 gx, gravrecvcounts, gravdispls, mpi_double_precision, 0,
     $                 comm_worker, ierr)
                  call mpi_gatherv(mpi_in_place, mygravlength, mpi_double_precision,
     $                 gy, gravrecvcounts, gravdispls, mpi_double_precision, 0,
     $                 comm_worker, ierr)
                  call mpi_gatherv(mpi_in_place, mygravlength, mpi_double_precision,
     $                 gz, gravrecvcounts, gravdispls, mpi_double_precision, 0,
     $                 comm_worker, ierr)
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
               else
                  ucgs=u(i)*gravconst*munit/runit
               endif
               if(arad.gt.0.d0.and.nintvar.eq.2.and.u(i).ne.0.d0)then
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

      if(autotf) then
         if(myrank.eq.0)
     $        write(69,*)'*************calling changetf***************',t
         call changetf
         if(myrank.eq.0)
     $        write(69,*)'*************done w/ changetf***************',t
      endif
      
      return
      end
