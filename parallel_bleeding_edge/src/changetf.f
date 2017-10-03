      subroutine changetf
c     pplot27: try to detect collisions
c     pplot28: if two stars are in process of merging, don't call it a new binary
c     pplot46: center of mass position and velocity used to calculate orbital elements, but position of lowest grpot particle in component is used when determining if particles are bound and which component is which color
      include 'starsmasher.h'
      include 'mpif.h'
      real*8 ax1,ay1,az1,ax2,ay2,az2,ax3,ay3,az3
      real*8 axr,ayr,azr,axg,ayg,azg,axb,ayb,azb
      integer icomp(nmax)
      real*8 am1,am2,am3,am4,amtot
      real*8 speedalpha
      real*8 speedbeta,sep
      real*8 actualspeedalphainf,actualspeedbetainf,thetaindegrees
      real*8 actualspeed3rdinf,actualspeedcominf,speed3rdinf,speedcominf
      real*8 speed3rd,speedcom,speedbetainf,speedalphainf
      real*8 vzdrft,vydrft,vxdrft
      real*8 x1,y1,z1,vx1,vy1,vz1
      real*8 x2,y2,z2,vx2,vy2,vz2
      real*8 x3,y3,z3,vx3,vy3,vz3
      common/compbettercom3/am1,x1,y1,z1,vx1,vy1,vz1,
     $     am2,x2,y2,z2,vx2,vy2,vz2,
     $     am3,x3,y3,z3,vx3,vy3,vz3,am4,
     $     icomp
      real*8 amr,amg,amb,xr,yr,zr,xg,yg,zg,xb,yb,zb
      save xr,yr,zr,xg,yg,zg,xb,yb,zb
      real*8 vxr,vyr,vzr,vxg,vyg,vzg,vxb,vyb,vzb
      save vxr,vyr,vzr,vxg,vyg,vzg,vxb,vyb,vzb
      save amr,amg,amb
      real*8 timelast
      save timelast
      logical ft
      save ft
      data ft /.true./
      real*8 agr,eccgr,lxgr,lygr,lzgr,dgr
      real*8 abr,eccbr,lxbr,lybr,lzbr,dbr
      real*8 abg,eccbg,lxbg,lybg,lzbg,dbg
      real*8 abin3,eccbin3,lxbin3,lybin3,lzbin3,lbin3
      real*8 a12,ecc12,lx12,ly12,lz12,l12
      real*8 xcm,ycm,zcm
      real*8 vxcm,vycm,vzcm
      real*8 xother,yother,zother
      real*8 vxother,vyother,vzother
      real*8 amother,ambin,ebin,ebin3ebin
      logical abinary
      real*8 mass1,mass2,mass3
      logical compinabinary(4)
      real*8 mcom,xcom,ycom,zcom,vxcom,vycom,vzcom
      real*8 m3rd,x3rd,y3rd,z3rd,vx3rd,vy3rd,vz3rd
      real*8 mce,mejecta
      real*8 ricom,v2icom,bicom
      integer blue,green,red
      save blue,green,red
c      character*1 stara,starb,starc
      integer istara,istarb,istarc
      character*9 state,laststate
      save laststate
      data laststate /'initial'/
      integer starnum
      character*128 allstates
      save allstates
      data allstates /''/
      real*8 rother
      character*9 starname(3)
      save starname
      integer numstars,lastnumstars
      data lastnumstars /3/
      save lastnumstars
      real*8 lastamb,lastamg,lastamr
      save lastamb,lastamg,lastamr
      data lastamb,lastamg,lastamr/0.d0,0.d0,0.d0/
      integer absorber,absorbee
      data absorbee/3/
      real*8 deltaamb,deltaamg,deltaamr
      real*8 rother1,rother2
      integer numstarslast
      common/nsl/ numstarslast
      real*8 xalpha,yalpha,zalpha
      real*8 vxalpha,vyalpha,vzalpha
      real*8 xbeta,ybeta,zbeta
      real*8 vxbeta,vybeta,vzbeta
      real*8 amalpha,ambeta
      real*8 atwo,ecctwo,rtwo
      logical bps               ! bps=binary plus single
      save bps
      logical checkifswapped
      common/cis/ checkifswapped
      real*8 r12
      real*8 tprintmassloss,eejecta,ekinejecta,eintejecta,epotejecta
      data tprintmassloss /1.d30/
      save tprintmassloss
      real*8 etotcheck
      integer numspeeds
      real*8 speed1, speed2, speedunit, eunit

      real*8 tfold
      integer i
      real*8 eintsave,eint,tjumpaheadold,dotproduct
      save eintsave
      data eintsave /-1.d30/
      common /jumpcomm/ tjumpahead
      real*8 dotproduct12
      integer mylength,ierr,irank
      real*8 ran1
      integer idumb
      save idumb
      data idumb /-2391/

c synchronize all variables:
      do i=1,n
         vx(i)=vx(i)-0.5d0*dt*vxdot(i)
         vy(i)=vy(i)-0.5d0*dt*vydot(i)
         vz(i)=vz(i)-0.5d0*dt*vzdot(i)
         if(u(i).ne.0.d0) then
            u(i)=u(i)-0.5d0*dt*udot(i)
         endif
      enddo
      t=t-0.5d0*dt

      inquire(file='input.3s',exist=bps)
c if the file input.3s exists in this directory, then make bps=.false.
c if the file input.3s does not exist, then assume binary+single
      bps= .not. bps

      if(myrank.eq.0)write(69,*)'************************* firsttime=',ft

c     make checkifswapped=.false. if the particle indices may change from
c     one outfile to the next.... (this happens in some of evghenii's gpu runs)
      checkifswapped=.true.
 88   continue

c     the rho array is used by compbest3 to figure out components, so let's assemble that array among processors
      mylength=n_upper-n_lower+1
      do irank=0,nprocs-1
         if(myrank.ne.irank)then
            call mpi_gatherv(rho(n_lower), mylength, mpi_double_precision,
     $           rho, recvcounts,displs, mpi_double_precision, irank,
     $           mpi_comm_world, ierr)
         else
            call mpi_gatherv(mpi_in_place, mylength, mpi_double_precision,
     $           rho, recvcounts,displs, mpi_double_precision, irank,
     $           mpi_comm_world, ierr)
         endif
      enddo

      call compbest3(.false.)

      ax1=0.d0
      ay1=0.d0
      az1=0.d0
      ax2=0.d0
      ay2=0.d0
      az2=0.d0
      ax3=0.d0
      ay3=0.d0
      az3=0.d0
      do i=1,n
         if(icomp(i).eq.1) then
            ax1=ax1+am(i)*vxdot(i)
            ay1=ay1+am(i)*vydot(i)
            az1=az1+am(i)*vzdot(i)
         elseif(icomp(i).eq.2) then
            ax2=ax2+am(i)*vxdot(i)
            ay2=ay2+am(i)*vydot(i)
            az2=az2+am(i)*vzdot(i)
         elseif(icomp(i).eq.3) then
            ax3=ax3+am(i)*vxdot(i)
            ay3=ay3+am(i)*vydot(i)
            az3=az3+am(i)*vzdot(i)
         endif
      enddo
      ax1=ax1/am1
      ay1=ay1/am1
      az1=az1/am1
      ax2=ax2/am2
      ay2=ay2/am2
      az2=az2/am2
      ax3=ax3/am3
      ay3=ay3/am3
      az3=az3/am3

c      write(69,'(a,3g17.9)')'mass of star 1=',am1
c      write(69,'(a,3g17.9)')'position of star 1=',x1,y1,z1
c      write(69,'(a,3g17.9)')'velocity of star 1=',vx1,vy1,vz1

c      write(69,'(a,3g17.9)')'mass of star 2=',am2
c      write(69,'(a,3g17.9)')'position of star 2=',x2,y2,z2
c      write(69,'(a,3g17.9)')'velocity of star 2=',vx2,vy2,vz2

c      write(69,'(a,3g17.9)')'mass of star 3=',am3
c      write(69,'(a,3g17.9)')'position of star 3=',x3,y3,z3
c      write(69,'(a,3g17.9)')'velocity of star 3=',vx3,vy3,vz3

c      do i=1,n
c         if(hp(i).le.0) then
c            write(69,'(a,3g17.9)')
c     $           'position of point mass=',x(i),y(i),z(i)
c            write(69,*)'        located in icomp',icomp(i)
c            write(69,'(a,3g17.9)')'     mass of particle=',am(i)
c         endif
c      enddo

      numstars=0
      if(am1.gt.0.d0) numstars=numstars+1
      if(am2.gt.0.d0) numstars=numstars+1
      if(am3.gt.0.d0) numstars=numstars+1
      if(myrank.eq.0)write(69,*)'number of stars=',numstars

      if(ft) then
         starname(1)='1'
         starname(2)='2'
         starname(3)='3'
c     first time through just make most massive component red(r), 2nd most massive green (g), and third most massive blue (b)

         if(myrank.eq.0)write(69,*)'this is the first time in this subroutine'
         if(myrank.eq.0)write(69,*)'making most massive star blue, next green, next red'
         if(am1.le.am2 .and. am1.le.am3) then
            red=1
            if(am2.le.am3) then
               green=2
               blue=3
            else
               green=3
               blue=2
            endif
         elseif(am2.le.am1 .and. am2.le.am3) then
            red=2
            if(am1.le.am3) then
               green=1
               blue=3
            else
               green=3
               blue=1
            endif
         else
            red=3
            if(am1.le.am2) then
               green=1
               blue=2
            else
               green=2
               blue=1
            endif
         endif
         timelast=t
c         if(myrank.eq.0)write(69,*)'timelast being set to',timelast
      endif

      if(red.eq.1)then
         if(myrank.eq.0)write(69,*)'component 1 is red'
         amr=am1
         xr=x1
         yr=y1
         zr=z1
         vxr=vx1
         vyr=vy1
         vzr=vz1
         axr=ax1
         ayr=ay1
         azr=az1
      elseif(red.eq.2) then
         if(myrank.eq.0)write(69,*)'component 2 is red'
         amr=am2
         xr=x2
         yr=y2
         zr=z2
         vxr=vx2
         vyr=vy2
         vzr=vz2
         axr=ax2
         ayr=ay2
         azr=az2
      else
         if(myrank.eq.0)write(69,*)'component 3 is red'
         amr=am3
         xr=x3
         yr=y3
         zr=z3
         vxr=vx3
         vyr=vy3
         vzr=vz3
         axr=ax3
         ayr=ay3
         azr=az3
      endif
      if(green.eq.1)then
         if(myrank.eq.0)write(69,*)'component 1 is green'
         amg=am1
         xg=x1
         yg=y1
         zg=z1
         vxg=vx1
         vyg=vy1
         vzg=vz1
         axg=ax1
         ayg=ay1
         azg=az1
      elseif(green.eq.2) then
         if(myrank.eq.0)write(69,*)'component 2 is green'
         amg=am2
         xg=x2
         yg=y2
         zg=z2
         vxg=vx2
         vyg=vy2
         vzg=vz2
         axg=ax2
         ayg=ay2
         azg=az2
      else
         if(myrank.eq.0)write(69,*)'component 3 is green'
         amg=am3
         xg=x3
         yg=y3
         zg=z3
         vxg=vx3
         vyg=vy3
         vzg=vz3
         axg=ax3
         ayg=ay3
         azg=az3
      endif
      if(blue.eq.1)then
         if(myrank.eq.0)write(69,*)'component 1 is blue'
         amb=am1
         xb=x1
         yb=y1
         zb=z1
         vxb=vx1
         vyb=vy1
         vzb=vz1
         axb=ax1
         ayb=ay1
         azb=az1
      elseif(blue.eq.2) then
         if(myrank.eq.0)write(69,*)'component 2 is blue'
         amb=am2
         xb=x2
         yb=y2
         zb=z2
         vxb=vx2
         vyb=vy2
         vzb=vz2
         axb=ax2
         ayb=ay2
         azb=az2
      else
         if(myrank.eq.0)write(69,*)'component 3 is blue'
         amb=am3
         xb=x3
         yb=y3
         zb=z3
         vxb=vx3
         vyb=vy3
         vzb=vz3
         axb=ax3
         ayb=ay3
         azb=az3
      endif

      if(red.le.0 .or. green.le.0 .or. blue.le.0 .or.
     $     red.eq.green .or. red.eq.blue .or. green.eq.blue
     $     .or.red+green+blue.ne.6) then
         if(myrank.eq.0)write(69,*)'star not assigned properly',red,green,blue
         stop
      endif

      dgr=((xg-xr)**2+(yg-yr)**2+(zg-zr)**2)**0.5d0
      dbr=((xr-xb)**2+(yr-yb)**2+(zr-zb)**2)**0.5d0
      dbg=((xg-xb)**2+(yg-yb)**2+(zg-zb)**2)**0.5d0

      amtot=0.d0
      do i=1,n
         amtot=amtot+am(i)
      enddo

      abinary=.false.
      compinabinary(1)=.false.
      compinabinary(2)=.false.
      compinabinary(3)=.false.
      compinabinary(4)=.false.
      a12=1.d30

      tfold=tf
      tjumpaheadold=tjumpahead

      if(amr.gt.0.d0 .and. amg.gt.0.d0) then
         call elements(amr,amg,xr-xg,yr-yg,zr-zg,
     $        vxr-vxg,vyr-vyg,vzr-vzg,
     $        agr,eccgr,lxgr,lygr,lzgr)
         dotproduct=(xg-xr)*(vxg-vxr)+
     $        (yg-yr)*(vyg-vyr)+
     $        (zg-zr)*(vzg-vzr)

         if(eccgr.lt.1.d0) then            
            if(myrank.eq.0)write(69,'(a,3g13.6)')
     $           'green and red are a binary with a,e,dotproduct=',
     $           agr,eccgr,dotproduct

c            if(min(amg,amr).lt.1.d0 .and. agr.lt.(amg+amr)/25.d0
c     $           .and. eccgr.gt.0.999d0 .and. eccgr.lt.1.d0
c     $           .and. abs(dotproduct).lt.0.1d0) then
            if((agr.lt.min((amg+amr)*0.04d0,10.d0)
     $           .and. eccgr.lt.1.d0) .or.
     $           (abs(dotproduct).lt.0.1d0 .and. eccgr.gt.0.98d0 .and.
     $           eccgr.lt.1.d0 .and. agr.lt.min(amg+amr,10.d0) ))then
               if(myrank.eq.0)write(69,*)'green and red seem to have merged'
               do i=1,n
                  if(amr.le.amg .and. icomp(i).eq.red) icomp(i)=green
                  if(amr.gt.amg .and. icomp(i).eq.green) icomp(i)=red
               enddo
               numstarslast=numstars-1
               checkifswapped=.true.
               goto 88
            endif

            compinabinary(green)=.true.
            compinabinary(red)=.true.
c            primary=green
c            secondary=red
c            intruder=blue
            abinary=.true.
            ebin=-0.5d0*amg*amr/agr
            ambin=amg+amr
            xcm=(amr*xr+amg*xg)/(amr+amg)
            ycm=(amr*yr+amg*yg)/(amr+amg)
            zcm=(amr*zr+amg*zg)/(amr+amg)
            vxcm=(amr*vxr+amg*vxg)/(amr+amg)
            vycm=(amr*vyr+amg*vyg)/(amr+amg)
            vzcm=(amr*vzr+amg*vzg)/(amr+amg)
            amother=amb
            xother=xb
            yother=yb
            zother=zb
            vxother=vxb
            vyother=vyb
            vzother=vzb
            rother1=((xother-xg)**2
     $           +(yother-yg)**2
     $           +(zother-zg)**2)**0.5d0
            rother2=((xother-xr)**2
     $           +(yother-yr)**2
     $           +(zother-zr)**2)**0.5d0
            mass1=amg
            mass2=amr
            mass3=amb
            ecc12=eccgr
            dotproduct12=dotproduct
            a12=agr
            r12=dgr
            lx12=lxgr
            ly12=lygr
            lz12=lzgr
c            if(ft .and. bps)then
c               starname(1)='i'  ! the numbers 1,2,3 correspond to a ranking of
c               starname(2)='p'  ! the initial masses, not to component numbers
c               starname(3)='s'
c               if(myrank.eq.0)write(69,*)'starnames=',(trim(starname(starnum)),starnum=1,3)
c            endif
         else
            if(a12.ge.1.d30 .and. ecc12.eq.0.d0) then
               mass1=amb
               mass2=amg
               mass3=amr
               ecc12=eccgr
               dotproduct12=dotproduct
               a12=agr
               r12=dgr
               lx12=lxgr
               ly12=lygr
               lz12=lzgr
            endif
         endif
      else
         agr=0.d0
         eccgr=-1.d0
      endif
      if(amr.gt.0.d0 .and. amb.gt.0.d0) then
         call elements(amr,amb,xr-xb,yr-yb,zr-zb,
     $        vxr-vxb,vyr-vyb,vzr-vzb,
     $        abr,eccbr,lxbr,lybr,lzbr)
         dotproduct=(xb-xr)*(vxb-vxr)+
     $        (yb-yr)*(vyb-vyr)+
     $        (zb-zr)*(vzb-vzr)
         if(eccbr.lt.1.d0 .and. (.not. abinary .or.
     $        abr*(1.d0+eccbr).le.a12*(1.d0+ecc12)))then 
            if(myrank.eq.0)write(69,'(a,3g13.6)')
     $           'blue and red are a binary with a,e,dotproduct=',
     $           abr,eccbr,dotproduct

c            if(min(amb,amr).lt.1.d0 .and. abr.lt.(amb+amr)/25.d0
c     $           .and. eccbr.gt.0.999d0 .and. eccbr.lt.1.d0
c     $           .and. abs(dotproduct).lt.0.1d0) then
            if((abr.lt.min((amb+amr)*0.04d0,10.d0)
     $           .and. eccbr.lt.1.d0) .or.
     $           (abs(dotproduct).lt.0.1d0 .and. eccbr.gt.0.98d0 .and.
     $           eccbr.lt.1.d0 .and. abr.lt.min(amb+amr,10.d0) ))then
               if(myrank.eq.0)write(69,*)'blue and red seem to have merged'
               do i=1,n
                  if(amr.le.amb .and. icomp(i).eq.red) icomp(i)=blue
                  if(amr.gt.amb .and. icomp(i).eq.blue) icomp(i)=red
               enddo
               numstarslast=numstars-1
               checkifswapped=.true.
               goto 88
            endif

            compinabinary(blue)=.true.
            compinabinary(red)=.true.
c            primary=blue
c            secondary=red
c            intruder=green
            if(abinary .and. abr*(1.d0+eccbr).le.a12*(1.d0+ecc12))then
               if(myrank.eq.0)write(69,*)'so green star not a binary member after all'
               compinabinary(green)=.false.
            endif
            if(.not. abinary .or. abr.le.a12)then
               abinary=.true.
               ebin=-0.5d0*amr*amb/abr
               ambin=amr+amb
               xcm=(amr*xr+amb*xb)/(amr+amb)
               ycm=(amr*yr+amb*yb)/(amr+amb)
               zcm=(amr*zr+amb*zb)/(amr+amb)
               vxcm=(amr*vxr+amb*vxb)/(amr+amb)
               vycm=(amr*vyr+amb*vyb)/(amr+amb)
               vzcm=(amr*vzr+amb*vzb)/(amr+amb)
               amother=amg
               xother=xg
               yother=yg
               zother=zg
               vxother=vxg
               vyother=vyg
               vzother=vzg
               rother1=((xother-xb)**2
     $              +(yother-yb)**2
     $              +(zother-zb)**2)**0.5d0
               rother2=((xother-xr)**2
     $              +(yother-yr)**2
     $              +(zother-zr)**2)**0.5d0
               mass1=amb
               mass2=amr
               mass3=amg
               ecc12=eccbr
               dotproduct12=dotproduct
               a12=abr
               r12=dbr
               lx12=lxbr
               ly12=lybr
               lz12=lzbr
c               if(ft .and. bps)then
c                  starname(1)='p' ! the numbers 1,2,3 correspond to a ranking of
c                  starname(2)='i' ! of initial masses, not the component numbers
c                  starname(3)='s'
c                  if(myrank.eq.0)write(69,*)'starnames=',(trim(starname(starnum)),
c     $                 starnum=1,3)
c               endif
            endif
         else
            if(a12.ge.1.d30 .and. ecc12.eq.0.d0) then
               mass1=amb
               mass2=amg
               mass3=amr
               ecc12=eccbr
               dotproduct12=dotproduct
               a12=abr
               r12=dbr
               lx12=lxbr
               ly12=lybr
               lz12=lzbr
            endif
         endif
      else
         abr=0.d0
         eccbr=0.d0
      endif
      if(amg.gt.0.d0 .and. amb.gt.0.d0) then
         call elements(amg,amb,xg-xb,yg-yb,zg-zb,
     $        vxg-vxb,vyg-vyb,vzg-vzb,
     $        abg,eccbg,lxbg,lybg,lzbg)
         dotproduct=(xb-xg)*(vxb-vxg)+
     $        (yb-yg)*(vyb-vyg)+
     $        (zb-zg)*(vzb-vzg)

         if(eccbg.lt.1.d0 .and. (.not. abinary .or.
     $        abg*(1.d0+eccbg).le.a12*(1.d0+ecc12)))then 
            if(myrank.eq.0)write(69,'(a,3g13.6)')
     $           'blue and green are a binary with a,e,dotproduct=',
     $           abg,eccbg,dotproduct

c            if(min(amb,amg).lt.1.d0 .and. abg.lt.(amb+amg)/25.d0
c     $           .and. eccbg.gt.0.999d0 .and. eccbg.lt.1.d0
c     $           .and. abs(dotproduct).lt.0.1d0) then
            if((abg.lt.min((amb+amg)*0.04d0,10.d0)
     $           .and. eccbg.lt.1.d0) .or.
     $           (abs(dotproduct).lt.0.1d0 .and. eccbg.gt.0.98d0 .and.
     $           eccbg.lt.1.d0 .and. abg.lt.min(amb+amg,10.d0) ))then
               if(myrank.eq.0)write(69,*)'blue and green seem to have merged'
               do i=1,n
                  if(amg.le.amb .and. icomp(i).eq.green) icomp(i)=blue
                  if(amg.gt.amb .and. icomp(i).eq.blue) icomp(i)=green
               enddo
               numstarslast=numstars-1
               checkifswapped=.true.
               goto 88
            endif

            compinabinary(blue)=.true.
            compinabinary(green)=.true.
c            primary=blue
c            secondary=green
c            intruder=red
            if(abinary .and. abg*(1.d0+eccbg).le.a12*(1.d0+ecc12))then
               if(myrank.eq.0)write(69,*)'so red star not a binary member after all'
               compinabinary(red)=.false.
            endif
            if(.not. abinary .or. abg.le.a12)then
               abinary=.true.
               ebin=-0.5d0*amg*amb/abg
               ambin=amg+amb
               xcm=(amg*xg+amb*xb)/(amg+amb)
               ycm=(amg*yg+amb*yb)/(amg+amb)
               zcm=(amg*zg+amb*zb)/(amg+amb)
               vxcm=(amg*vxg+amb*vxb)/(amg+amb)
               vycm=(amg*vyg+amb*vyb)/(amg+amb)
               vzcm=(amg*vzg+amb*vzb)/(amg+amb)
               amother=amr
               xother=xr
               yother=yr
               zother=zr
               vxother=vxr
               vyother=vyr
               vzother=vzr
               rother1=((xother-xb)**2
     $              +(yother-yb)**2
     $              +(zother-zb)**2)**0.5d0
               rother2=((xother-xg)**2
     $              +(yother-yg)**2
     $              +(zother-zg)**2)**0.5d0
               mass1=amb
               mass2=amg
               mass3=amr
               ecc12=eccbg
               dotproduct12=dotproduct
               a12=abg
               r12=dbg
               lx12=lxbg
               ly12=lybg
               lz12=lzbg
c               if(ft .and. bps)then
c                  starname(1)='p' ! the numbers 1,2,3 correspond to a ranking of
c                  starname(2)='s' ! the initial masses, not the component numbers
c                  starname(3)='i'
c                  if(myrank.eq.0)write(69,*)'starnames=',(trim(starname(starnum)),
c     $                 starnum=1,3)
c               endif
            endif
         else
            if(a12.ge.1.d30 .and. ecc12.eq.0.d0) then
               mass1=amb
               mass2=amg
               mass3=amr
               ecc12=eccbg
               dotproduct12=dotproduct
               a12=abg
               r12=dbg
               lx12=lxbg
               ly12=lybg
               lz12=lzbg
            endif
         endif
      else
         abg=0.d0
         eccbg=-1.d0
      endif

      timelast=t
      if(numstars.ne.lastnumstars) then
         if(myrank.eq.0)write(69,*)lastnumstars-numstars+1,'stars have merged'
         if(t.gt.1 .and. tjumpahead.ge.0.d0) then
            tf=min(tfold,dble(nint(t+4000.d0)))
            tjumpahead=1.d30
         endif
         if(lastnumstars-numstars.gt.1) then
            if(myrank.eq.0)write(69,*)'two stars disappeared at once?'
            tf=max(tfold,dble(nint(t+4000.d0)))
            tjumpahead=1.d30
c            stop
         endif
            
         if(amb.eq.0.d0 .and. lastamb.gt.0.d0) absorbee=1
         if(amg.eq.0.d0 .and. lastamg.gt.0.d0) absorbee=2
         if(amr.eq.0.d0 .and. lastamr.gt.0.d0) absorbee=3
         if(myrank.eq.0)write(69,*) amb,lastamb,amg,lastamg,amr,lastamr,absorbee
         deltaamb=amb-lastamb
         deltaamg=amg-lastamg
         deltaamr=amr-lastamr
         if(myrank.eq.0)write(69,*)'deltas=',
     $        deltaamb,deltaamg,deltaamr
         if(deltaamb.gt.deltaamg .and. deltaamb.gt.deltaamr) then
            absorber=1
         else if(deltaamg.gt.deltaamb .and. deltaamg.gt.deltaamr) then
            absorber=2
         else if(deltaamr.gt.deltaamg .and. deltaamr.gt.deltaamb) then
            absorber=3
         else
            if(myrank.eq.0)write(69,*)'cannot determine absorber',
     $           deltaamb,deltaamg,deltaamr
            stop
         endif
         lastnumstars=numstars
         if(myrank.eq.0)write(69,*)'star ',trim(starname(absorber)),
     $        ' has absorbed star ',trim(starname(absorbee)),
     $        absorber,absorbee
         starname(absorber)=
     $        '{'//trim(starname(absorber))//','
     $        //trim(starname(absorbee))//'}'
      else if(numstars.gt.1 .and. ecc12.ge.1.d0 .and.
     $        dotproduct12.gt.0) then
         if(myrank.eq.0)write(69,*)
     $        'the stars have not merged, and they never will: ecc=',ecc12
         tf=min(tfold,dble(nint(t+4000.d0)))
         tjumpahead=1.d30
      else if(numstars.gt.1) then
         if(ecc12.gt.1.d0) then
            if(myrank.eq.0)write(69,*)
     $           'the stars have not merged, but they may still'
            tjumpahead=1.d30
         else
            if(myrank.eq.0)write(69,*)
     $           'the stars have not merged, but they will'
            if(dotproduct12.gt.0 .and. a12*(1.d0+ecc12).gt.
     $           1000.d0) then
c               tjumpahead=1d30
               tjumpahead=min(tjumpahead,dble(nint(t+4000.d0)))
c               tjumpahead=min(tjumpahead,dble(nint(t+75.d0+25*ran1(idumb))))
               tf=max(dble(nint(tjumpahead+4000.d0)),tf)
            else
               if(myrank.eq.0)write(69,*)'not going to jump ahead (',dotproduct12,
     $              a12*(1.d0+ecc12),')!'
               tjumpahead=1.d30
            endif
         endif
         tf=dble(nint(t+4000.d0))
      endif
      tjumpahead=abs(tjumpahead)
      lastamb=amb
      lastamg=amg
      lastamr=amr

      if(eintsave.gt.0.d0) then
         eint=0.d0
         do i=1,n
            if(nintvar.eq.1) then
c     p=(gam-1)*rho*u=a*rho^gam, so u=a*rho^(gam-1)/(gam-1)
               eint=eint+am(i)*u(i)*rho(i)**(gam-1.d0)/(gam-1.d0)
            else
               eint=eint+am(i)*u(i)
            endif
         enddo
         if(abs(1.d0-eint/eintsave).gt.0.01d0) then
            if(myrank.eq.0)write(69,*)'u is changing a lot: will integrate longer'
            tf=max(tfold,dble(nint(t+1000.d0)))
c            tjumpahead=1.d30
         endif
      endif

      if(abinary) then
c     there is a binary:
         rother=((xother-xcm)**2
     $        +(yother-ycm)**2
     $        +(zother-zcm)**2)**0.5d0
         
         starnum=0
         if(compinabinary(blue))then
c            stara='b'
            istara=1
            starnum=starnum+1
         else
c            starc='b'
            istarc=1
         endif
         if(compinabinary(green))then
            if(starnum.eq.0) then
c               stara='g'
               istara=2
            else
c               starb='g'
               istarb=2
            endif
            starnum=starnum+1
         else
c            starc='g'
            istarc=2
         endif
         if(compinabinary(red))then
            if(starnum.eq.0) then
c               stara='r'
               istara=3
               if(myrank.eq.0)write(69,*)'stara should not be red'
               stop
            else
c               starb='r'
               istarb=3
            endif
            starnum=starnum+1
         else
c            starc='r'
            istarc=3
         endif

      else
c     there are just individual stars:
         if(myrank.eq.0)write(69,'(a,3g13.6)')
     $        'there are just individual stars with masses',amb,amg,amr
         starnum=0
         if(amb.gt.0.d0)then
c            stara='b'
            istara=1
            amalpha=amb
            xalpha=xb
            yalpha=yb
            zalpha=zb
            vxalpha=vxb
            vyalpha=vyb
            vzalpha=vzb
            starnum=starnum+1
         endif
         if(amg.gt.0.d0)then
            if(starnum.eq.0) then
c               stara='g'
               istara=2
               amalpha=amg
               xalpha=xg
               yalpha=yg
               zalpha=zg
               vxalpha=vxg
               vyalpha=vyg
               vzalpha=vzg
            else
c               starb='g'
               istarb=2
               ambeta=amg
               xbeta=xg
               ybeta=yg
               zbeta=zg
               vxbeta=vxg
               vybeta=vyg
               vzbeta=vzg
            endif
            starnum=starnum+1
         endif
         if(amr.gt.0.d0)then
            if(starnum.eq.0) then
c               stara='r'
               istara=3
               amalpha=amr
               xalpha=xr
               yalpha=yr
               zalpha=zr
               vxalpha=vxr
               vyalpha=vyr
               vzalpha=vzr
            else if(starnum.eq.1) then
c               starb='r'
               istarb=3
               ambeta=amr
               xbeta=xr
               ybeta=yr
               zbeta=zr
               vxbeta=vxr
               vybeta=vyr
               vzbeta=vzr
            else
c               starc='r'
               istarc=3
            endif
            starnum=starnum+1
         else
c            starc='r'
            istarc=3
         endif
 
         if(myrank.eq.0)write(69,*)'number of individual stars=',starnum
         
         if(starnum.eq.1) then
c     there is a single object
c            write(state,201) stara
            speedalpha=(vxalpha**2+vyalpha**2+vzalpha**2)**0.5d0
            write(state,201) trim(starname(istara))
 201        format(a)
            if(myrank.eq.0)write(69,*)'speed of ',trim(starname(istara)),
     $           speedalpha
            if(t.gt.tprintmassloss) then
               if(myrank.eq.0)write(69,'(a)')
     $              '> speed info (a single object); fl; e_ej'
               if(myrank.eq.0)write(69,'(a,g9.2,a,f7.1)')
     $              '  speed info (a single object): &',speedalpha,
     $              '  % t=',t
               numspeeds=1
               speed1=speedalpha
            endif
         else if(starnum.eq.2) then
c     there are two single objects
c            write(state,202) stara,starb
            vxdrft=(amalpha*vxalpha+ambeta*vxbeta)/(amalpha+ambeta)
            vydrft=(amalpha*vyalpha+ambeta*vybeta)/(amalpha+ambeta)
            vzdrft=(amalpha*vzalpha+ambeta*vzbeta)/(amalpha+ambeta)
            speedalpha=((vxalpha-vxdrft)**2+(vyalpha-vydrft)**2
     $           +(vzalpha-vzdrft)**2)**0.5d0
            speedbeta=((vxbeta-vxdrft)**2+(vybeta-vydrft)**2
     $           +(vzbeta-vzdrft)**2)**0.5d0
            if(myrank.eq.0)write(69,*)'speed of ',trim(starname(istara)),
     $           speedalpha     ! this is in c.o.m. frame of these two stars
            if(myrank.eq.0)write(69,*)'speed of ',trim(starname(istarb)),
     $           speedbeta      ! this is in c.o.m. frame of these two stars
            sep=((xalpha-xbeta)**2+
     $           (yalpha-ybeta)**2+
     $           (zalpha-zbeta)**2)**0.5d0
            speedalphainf=(
     $           speedalpha**2
     $           -2.d0*ambeta**2/(sep*(amalpha+ambeta))
     $           )**0.5d0
            speedbetainf=(
     $           speedbeta**2
     $           -2.d0*amalpha**2/(sep*(amalpha+ambeta))
     $           )**0.5d0
            if(myrank.eq.0)write(69,*)'at infinity, speed of ',trim(starname(istara)),
     $           speedalphainf  ! this is in c.o.m. frame of these two stars
            if(myrank.eq.0)write(69,*)'at infinity, speed of ',trim(starname(istarb)),
     $           speedbetainf   ! this is in c.o.m. frame of these two stars

            if(myrank.eq.0)write(69,'(a,3g13.6)')
     $           'drift velocity components:',vxdrft,vydrft,vzdrft
c     shift back to original inertial reference frame:
            actualspeedalphainf=(
     $           (vxdrft+speedalphainf*(vxalpha-vxdrft)/speedalpha)**2+
     $           (vydrft+speedalphainf*(vyalpha-vydrft)/speedalpha)**2+
     $           (vzdrft+speedalphainf*(vzalpha-vzdrft)/speedalpha)**2
     $           )**0.5d0
            actualspeedbetainf=(
     $           (vxdrft+speedbetainf*(vxbeta-vxdrft)/speedbeta)**2+
     $           (vydrft+speedbetainf*(vybeta-vydrft)/speedbeta)**2+
     $           (vzdrft+speedbetainf*(vzbeta-vzdrft)/speedbeta)**2
     $           )**0.5d0
            if(myrank.eq.0)write(69,'(2a,3g13.6)')
     $           'at infinity, actual speed of ',
     $           trim(starname(istara)),
     $           actualspeedalphainf,actualspeedalphainf*amalpha
            if(myrank.eq.0)write(69,'(2a,3g13.6)')
     $           'at infinity, actual speed of ',
     $           trim(starname(istarb)),
     $           actualspeedbetainf,actualspeedbetainf*ambeta

            if(t.gt.tprintmassloss) then
               if(myrank.eq.0)write(69,'(2(a,g9.2),a,f7.1)')
     $       '> @ inf, actual speed info (2 single objects);fl;e_ej'
               if(myrank.eq.0)write(69,'(2(a,g9.2),a,f7.1)')
     $       '  at infinity, actual speed info (two single objects): &',
     $              actualspeedalphainf,',',
     $              actualspeedbetainf,'  % t=',t
               numspeeds=2
               speed1=actualspeedalphainf
               speed2=actualspeedbetainf
            endif
            dotproduct=(xalpha-xbeta)*(vxalpha-vxbeta)+
     $           (yalpha-ybeta)*(vyalpha-vybeta)+
     $           (zalpha-zbeta)*(vzalpha-vzbeta)
            if(istara+istarb.eq.3) then
               atwo=abg
               ecctwo=eccbg
            else if(istara+istarb.eq.4) then
               atwo=abr
               ecctwo=eccbr
            else if(istara+istarb.eq.5) then
               atwo=agr
               ecctwo=eccgr
            endif
            rtwo=((xalpha-xbeta)**2
     $           +(yalpha-ybeta)**2
     $           +(zalpha-zbeta)**2)**0.5d0

            if(myrank.eq.0)write(69,*)'dp=',atwo,ecctwo,dotproduct,rtwo
            if(dotproduct.gt.0.d0 .and.
     $           rtwo.ge.amalpha+ambeta .and.
     $           rtwo.ge.abs(atwo))then
               write(state,202) trim(starname(istara)),
     $              trim(starname(istarb))
 202           format(a,',',a)
            else
               state=laststate
            endif
         else if(starnum.eq.3) then
c     there are three single objects
c            write(state,202) stara,starb,starc
            write(state,203) trim(starname(istara)),
     $           trim(starname(istarb)),
     $           trim(starname(istarc))
 203        format(a,',',a,',',a)
         endif

      endif
         
      if(abinary .and. amother.gt.0.d0) then
c     there is a binary and a third star:
         call elements(ambin,amother,xcm-xother,ycm-yother,zcm-zother,
     $        vxcm-vxother,vycm-vyother,vzcm-vzother,
     $        abin3,eccbin3,lxbin3,lybin3,lzbin3)
         dotproduct=(xcm-xother)*(vxcm-vxother)+
     $        (ycm-yother)*(vycm-vyother)+
     $        (zcm-zother)*(vzcm-vzother)
         
         if(myrank.eq.0)write(69,*)'dotproduct for other star and binary=',dotproduct

         if(eccbin3.lt.1.d0) then
            compinabinary(1)=.true.
            compinabinary(2)=.true.
            compinabinary(3)=.true.
         endif
         if(ft .and. bps) then
            write(state,102) trim(starname(istara)),
     $           trim(starname(istarb)),
     $           trim(starname(istarc))
         else if(eccbin3.lt.1.d0 .or.
     $           (dotproduct.lt.0.d0 .and.
     $           abin3*(1.d0-eccbin3).le.a12*(1.d0+ecc12)).or.
     $           min(min(rother,rother1),rother2).lt.
     $                 abin3*(1.d0-eccbin3) .or.
     $           min(min(rother,rother1),rother2).lt.
     $                 a12*(1.d0-ecc12)
     $           )then
c     there is a binary and a third star that is close or will eventually come in close to the binary:
            if(myrank.eq.0)write(69,'(a,9g13.6)')
     $           'this is a triple or resonance, star 3s a,e=',
     $           abin3,eccbin3

c     keep in same state unless third star is both sufficiently far away and moving away:
c            if(myrank.eq.0)write(69,*)'debug r12=',r12,ambin,a12,a12*(1.d0-ecc12)
            if(rother.ge.abin3*(1.d0-eccbin3)
     $           .and. rother.ge.a12*(1.d0+ecc12) .and.
     $           rother1.ge.abin3*(1.d0-eccbin3)
     $           .and. rother1.ge.a12*(1.d0+ecc12) .and.
     $           rother2.ge.abin3*(1.d0-eccbin3)
     $           .and. rother2.ge.a12*(1.d0+ecc12) .and.
     $           dotproduct.gt.0.d0 .and.
     $           r12.ge.ambin/10.d0 )then ! this last condition is just to prevent cases where the binary is about to merge as counting as a resonance... this condition could be improved if necessary
               if(myrank.eq.0)write(69,*)'rother,abin3,a12*(1+ecc12)=',
     $              rother,abin3,a12*(1.d0+ecc12)
               write(state,301) trim(starname(1)),
     $              trim(starname(2)),
     $              trim(starname(3))
 301           format('(',a,',',a,',',a,')')
            else if(bps) then
               state=laststate
            else
               state='(1,2,3)'
            endif

            tf=max(tfold,dble(nint(t+4000.d0)))
            tjumpahead=1.d30

         else
c     there is a binary and a third star that is both not bound and not headed close to the binary:
            write(state,102) trim(starname(istara)),
     $           trim(starname(istarb)),
     $           trim(starname(istarc))
 102        format('(',a,',',a,'),',a)

            tjumpahead=min(tjumpahead,dble(nint(t+100.d0)))
            tf=max(dble(nint(tjumpahead+4000.d0)),tf)

         endif
         
         ebin3ebin=-0.5d0*ambin*amother/abin3/ebin
         if(myrank.eq.0)write(69,*)'e_{binary,other}/e_{binary}=',ebin3ebin

         l12=(lx12**2+ly12**2+lz12**2)**0.5d0
         lbin3=(lxbin3**2+lybin3**2+lzbin3**2)**0.5d0
         thetaindegrees=acos((lx12*lxbin3+ly12*lybin3+lz12*lzbin3)
     $        /(l12*lbin3))*
     $        180/pi
         if(myrank.eq.0)write(69,'(3(a,g10.3),1(a,g10.2),3(a,g10.3),a,g10.3,3a)')
     $        '&','mass1','&','mass2','&','mass3',
     $        '&','e_12','&','a_12','&','periastron','&',
     $        'e_b3','&',
     $        'e_b3/e_12','&','theta','\\'
         if(myrank.eq.0)write(69,'(3(a,g10.3),1(a,f10.2),3(a,g10.3),a,g10.3,a,i10,a)')
     $        '&',mass1,'&',mass2,'&',mass3,
     $        '&',ecc12,'&',a12,'&',
     $        abin3*(1.d0-eccbin3),'&',eccbin3,'&',
     $        mass3*(mass1+mass2)/(mass1*mass2)*a12/abin3,'&',
     $        nint(thetaindegrees),
     $        '$^\circ$ \\'
         if(myrank.eq.0)write(69,'(3(a,g10.3),1(a,g10.2),3(a,g10.3),a,g10.3,3a)')
     $        '&','mass1','&','mass2','&','mass3',
     $        '&','x_1','&','y_1','&','z_1','&',
     $        '&','vx_1','&','vy_1','&','vz_1','&',
     $        '&','x_2','&','y_2','&','z_2','&',
     $        '&','vx_2','&','vy_2','&','vz_2','&',
     $        '\\'
         if(myrank.eq.0)write(69,'(99(a,g10.3))')
     $        '&',amb,'&',amg,'&',amr,
     $           '&',xb-vxb*dt,'&',yb-vyb*dt,'&',zb-vzb*dt,
     $        '&',vxb-axb*dt,'&',vyb-ayb*dt,'&',vzb-azb*dt,
     $           '&',xg-vxg*dt,'&',yg-vyg*dt,'&',zg-vzg*dt,
     $        '&',vxg-axg*dt,'&',vyg-ayg*dt,'&',vzg-azg*dt,
     $        '\\ %',xr-vxr*dt,'&',yr-vyr*dt,'&',zr-vzr*dt,
     $        '&',vxr-axr*dt,'&',vyr-ayr*dt,'&',vzr-azr*dt
      else
         abin3=0.d0
         eccbin3=-1.d0
         ebin3ebin=0.d0
         if(abinary) then
c     there is a binary but no third star:
c            write(state,103) stara,starb
            write(state,103) trim(starname(istara)),
     $           trim(starname(istarb))
 103        format('(',a,',',a,')')
         endif
      endif






      if(tf.ne.tfold) then
         eintsave=0.d0
         do i=1,n
            if(nintvar.eq.1) then
               eintsave=eintsave
     $              +am(i)*u(i)*rho(i)**(gam-1.d0)/(gam-1.d0)
            else
               eintsave=eintsave+am(i)*u(i)
            endif
         enddo
         if(myrank.eq.0)
     $        write(69,'(3(a,g13.5))')'tf changed from',tfold,'to',tf,
     $        'when eint=',eintsave
      endif         

      if(tjumpahead.ne.tjumpaheadold) then
         if(myrank.eq.0)write(69,'(3(a,g13.5))')'tjumpahead changed from',tjumpaheadold,
     $        'to',tjumpahead
      endif         





      ft=.false.

      if(myrank.eq.0)write(69,*)'state=',state

      if(state.ne.laststate) then
         laststate=state
         if(allstates.eq.'')then
            allstates=state
         else
            allstates=trim(allstates)//' -> '//state
            tprintmassloss=t+200
            if(myrank.eq.0)write(69,*)'time to print mass loss data set to',tprintmassloss
         endif
         if(myrank.eq.0)write(69,*)trim(allstates)
      endif

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c here we'll determine how much of the mass in component 4 is
c in a bound common envelope and how much is actually ejecta

      xcom=0.d0
      ycom=0.d0
      zcom=0.d0
      vxcom=0.d0
      vycom=0.d0
      vzcom=0.d0
      x3rd=0.d0
      y3rd=0.d0
      z3rd=0.d0
      vx3rd=0.d0
      vy3rd=0.d0
      vz3rd=0.d0
      mcom=0.d0
      m3rd=0.d0
      do i=1,n
         if(compinabinary(icomp(i))) then
            mcom=mcom+am(i)
            xcom=xcom+am(i)*x(i)
            ycom=ycom+am(i)*y(i)
            zcom=zcom+am(i)*z(i)
            vxcom=vxcom+am(i)*vx(i)
            vycom=vycom+am(i)*vy(i)
            vzcom=vzcom+am(i)*vz(i)
         else if(icomp(i).ne.4) then
            m3rd=m3rd+am(i)
            x3rd=x3rd+am(i)*x(i)
            y3rd=y3rd+am(i)*y(i)
            z3rd=z3rd+am(i)*z(i)
            vx3rd=vx3rd+am(i)*vx(i)
            vy3rd=vy3rd+am(i)*vy(i)
            vz3rd=vz3rd+am(i)*vz(i)
         endif
      enddo
      
      mce=0.d0                  !(mce=mass in common envelope)
      mejecta=0.d0
      if(mcom.gt.0.d0) then
         if(myrank.eq.0)write(69,*)'total mass in binary or multiple=',mcom
         xcom=xcom/mcom
         ycom=ycom/mcom
         zcom=zcom/mcom
         vxcom=vxcom/mcom
         vycom=vycom/mcom
         vzcom=vzcom/mcom
         if(myrank.eq.0)write(69,'(a,3g13.6)')
     $        'center of binary or multiple=',xcom,ycom,zcom
         if(myrank.eq.0)write(69,'(a,3g13.6)')
     $        'velocity of binary or multiple=',vxcom,vycom,vzcom

         ekinejecta=0.d0
         eintejecta=0.d0
         epotejecta=0.d0
         etotcheck=0.d0
         do i=1,n
            if(icomp(i).eq.4) then
               ricom=sqrt((x(i)-xcom)**2+(y(i)-ycom)**2+(z(i)-zcom)**2)
               v2icom=(vx(i)-vxcom)**2+(vy(i)-vycom)**2+(vz(i)-vzcom)**2
               bicom=0.5d0*v2icom+u(i)-mcom/ricom
c               bicom=0.5d0*v2icom-mcom/ricom
               if(bicom.le.0.d0) then
                  mce=mce+am(i)
               else
                  mejecta=mejecta+am(i)
                  ekinejecta=ekinejecta+0.5d0*am(i)*
     $                 (vx(i)**2+vy(i)**2+vz(i)**2)
                  eintejecta=eintejecta+am(i)*u(i)
                  epotejecta=epotejecta+am(i)*grpot(i)
               endif
            endif
            etotcheck=etotcheck+am(i)*u(i)+0.5d0*am(i)*
     $           (vx(i)**2+vy(i)**2+vz(i)**2)+0.5d0*am(i)*grpot(i)
         enddo
         eejecta=ekinejecta+eintejecta+epotejecta

         speedcom=(vxcom**2+vycom**2+vzcom**2)**0.5d0
         if(numstars.eq.2) then
            speedcom=(vxcom**2+vycom**2+vzcom**2)**0.5d0
            if(myrank.eq.0)write(69,*)'speed of binary:',
     $           speedcom,speedcom*mcom
            if(t.gt.tprintmassloss) then
               if(myrank.eq.0)write(69,'(a,g9.2,a,f7.1)')
     $              '> speed info (c.o.m. of binary); fl; e_ej'
               if(myrank.eq.0)write(69,'(a,g9.2,a,f7.1)')
     $              '  speed info (c.o.m. of binary): &',speedcom,
     $              '  % t=',t
               numspeeds=1
               speed1=speedcom
            endif
         elseif(numstars.eq.3 .and.m3rd.eq.0.d0) then
            speedcom=(vxcom**2+vycom**2+vzcom**2)**0.5d0
            if(myrank.eq.0)write(69,*)'speed of triple:',
     $           speedcom,speedcom*mcom
            if(t.gt.tprintmassloss) then
               if(myrank.eq.0)write(69,'(a,g9.2,a,f7.1)')
     $              '> speed info (c.o.m. of triple); fl; e_ej'
               if(myrank.eq.0)write(69,'(a,g9.2,a,f7.1)')
     $              '  speed info (c.o.m. of triple): &',speedcom,
     $              '  % t=',t
               numspeeds=1
               speed1=speedcom
            endif
         elseif(numstars.eq.3) then
c     there must be a binary and a third star
            x3rd=x3rd/m3rd
            y3rd=y3rd/m3rd
            z3rd=z3rd/m3rd
            vx3rd=vx3rd/m3rd
            vy3rd=vy3rd/m3rd
            vz3rd=vz3rd/m3rd
            vxdrft=(mcom*vxcom+m3rd*vx3rd)/(mcom+m3rd)
            vydrft=(mcom*vycom+m3rd*vy3rd)/(mcom+m3rd)
            vzdrft=(mcom*vzcom+m3rd*vz3rd)/(mcom+m3rd)
            speedcom=((vxcom-vxdrft)**2+(vycom-vydrft)**2
     $           +(vzcom-vzdrft)**2)**0.5d0
            speed3rd=((vx3rd-vxdrft)**2+(vy3rd-vydrft)**2
     $           +(vz3rd-vzdrft)**2)**0.5d0
            if(myrank.eq.0)write(69,*)'speed of binary:',
     $           speedcom,speedcom*mcom ! this is in c.o.m. frame of the binary + single
            if(myrank.eq.0)write(69,*)'speed of 3rd star:',
     $           speed3rd,speed3rd*m3rd ! this is in c.o.m. frame of the binary + single

            sep=((xcom-x3rd)**2+
     $           (ycom-y3rd)**2+
     $           (zcom-z3rd)**2)**0.5d0
            speedcominf=(
     $           speedcom**2
     $           -2.d0*m3rd**2/(sep*(mcom+m3rd))
     $           )**0.5d0
            speed3rdinf=(
     $           speed3rd**2
     $           -2.d0*mcom**2/(sep*(mcom+m3rd))
     $           )**0.5d0
            if(myrank.eq.0)write(69,*)'at infinity, speed of binary:',
     $           speedcominf,speedcominf*mcom ! this is in c.o.m. frame of the binary + single
            if(myrank.eq.0)write(69,*)'at infinity, speed of 3rd star:',
     $           speed3rdinf,speed3rdinf*m3rd ! this is in c.o.m. frame of the binary + single

            if(myrank.eq.0)write(69,'(a,3g13.6)')
     $           'drift velocity components:',vxdrft,vydrft,vzdrft
c     shift back to original inertial reference frame:
            actualspeedcominf=(
     $           (vxdrft+speedcominf*(vxcom-vxdrft)/speedcom)**2+
     $           (vydrft+speedcominf*(vycom-vydrft)/speedcom)**2+
     $           (vzdrft+speedcominf*(vzcom-vzdrft)/speedcom)**2
     $           )**0.5d0
            actualspeed3rdinf=(
     $           (vxdrft+speed3rdinf*(vx3rd-vxdrft)/speed3rd)**2+
     $           (vydrft+speed3rdinf*(vy3rd-vydrft)/speed3rd)**2+
     $           (vzdrft+speed3rdinf*(vz3rd-vzdrft)/speed3rd)**2
     $           )**0.5d0
            if(myrank.eq.0)write(69,'(a,3g13.6)')
     $           'at infinity, actual speed of binary:',
     $           actualspeedcominf,actualspeedcominf*mcom
            if(myrank.eq.0)write(69,'(a,3g13.6)')
     $           'at infinity, actual speed of 3rd star:',
     $           actualspeed3rdinf,actualspeed3rdinf*m3rd

            if(t.gt.tprintmassloss) then
               if(myrank.eq.0)write(69,'(2(a,g9.2),a,f7.1)')
     $    '> @ inf, actual speed info (com of binary, a single);fl;e_ej'
               if(myrank.eq.0)write(69,'(2(a,g9.2),a,f7.1)')
     $    '  at inf, actual speed info (c.o.m. of binary, a single): &',
     $              actualspeedcominf,',',
     $              actualspeed3rdinf,'  % t=',t
               numspeeds=2
               speed1=actualspeedcominf
               speed2=actualspeed3rdinf
            endif
         endif
      else
         if(myrank.eq.0)write(69,*) 'there are no binaries or other multiples'
         mejecta=am4
         ekinejecta=0.d0
         eintejecta=0.d0
         epotejecta=0.d0
         etotcheck=0.d0
         do i=1,n
            if(icomp(i).eq.4) then
               ekinejecta=ekinejecta+0.5d0*am(i)*
     $              (vx(i)**2+vy(i)**2+vz(i)**2)
               eintejecta=eintejecta+am(i)*u(i)
               epotejecta=epotejecta+am(i)*grpot(i)
            endif
            etotcheck=etotcheck+am(i)*u(i)+0.5d0*am(i)*
     $           (vx(i)**2+vy(i)**2+vz(i)**2)+0.5d0*am(i)*grpot(i)
         enddo
         eejecta=ekinejecta+eintejecta+epotejecta

      endif

      if(myrank.eq.0)write(69,*)'mcom=',mcom,'m3rd=',m3rd,'mcom+m3rd=',mcom+m3rd,
     $     'amtot=',amtot,'mejecta=',mejecta,'am4=',am4

      if(t.gt.tprintmassloss) then
         if(myrank.eq.0)write(69,'(5(a,g9.2),a,f7.1)')
     $      '  ejecta info (mass loss fraction and ejecta energies): &',
     $        mejecta/amtot,'&',
     $        eejecta,'&',
     $        ekinejecta,'&',
     $        eintejecta,'&',
     $        epotejecta,'  % t=',t         

         speedunit=(gravconst*munit/runit)**0.5d0*1.d-5
         eunit=gravconst*munit**2/runit*1.d-48
         if(myrank.eq.0)write(69,*)'speed and e units=',speedunit,eunit

         if(numspeeds.eq.1) then
            if(mejecta/amtot.ge.0.0095d0) then
               if(myrank.eq.0)write(69,'(a,g9.2,a,f9.3,a,g9.2,a,f7.1)')
     $              '> &',
     $              speed1*speedunit,'&',
     $              mejecta/amtot,'&',
     $              eejecta*eunit,'  \\% t=',t
            else
               if(myrank.eq.0)write(69,'(a,g9.1,a,g9.1,a,g9.1,a,f7.1)')
     $              '> &',
     $              speed1*speedunit,'&',
     $              mejecta/amtot,'&',
     $              eejecta*eunit,'  \\% t=',t               
            endif
         elseif(numspeeds.eq.2) then
            if(mejecta/amtot.ge.0.0095d0) then
               if(myrank.eq.0)write(69,'(2(a,i4),a,f9.3,a,g9.2,a,f7.1)')
     $              '> &',
     $              nint(speed1*speedunit),',',
     $              nint(speed2*speedunit),'&',
     $              mejecta/amtot,'&',
     $              eejecta*eunit,'  \\% t=',t
            else
               if(myrank.eq.0)write(69,'(2(a,i4),a,g9.1,a,g9.1,a,f7.1)')
     $              '> &',
     $              nint(speed1*speedunit),',',
     $              nint(speed2*speedunit),'&',
     $              mejecta/amtot,'&',
     $              eejecta*eunit,'  \\% t=',t
            endif
         else
            if(myrank.eq.0)write(69,*)'??? numspeeds=',numspeeds
            stop
         endif


         tprintmassloss=1.d30
      endif
      if(myrank.eq.0)write(69,*)'fractional mass loss:',mejecta/amtot
      if(myrank.eq.0)write(69,'(a,4g13.5)')
     $     'energies (e,t,u,w) in ejecta:',eejecta,
     $     ekinejecta,eintejecta,epotejecta
      
      if(myrank.eq.0)write(69,*)'etotcheck=',etotcheck,t


      if(abs(am4-mce-mejecta).gt.1.d-5) then
         if(myrank.eq.0)write(69,*)'conservation of mass problem:',am4,mce,mejecta
         stop
      endif
      if(myrank.eq.0)write(69,'(3(a,g13.6))')'m4=mce+mejecta',
     $     am4,'=',mce,'+',mejecta

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      write(23,'(d12.5,3d16.8,33d12.4)') t,amb,amg,amr,
     $     xb,yb,zb,xg,yg,zg,xr,yr,zr,
     $     dbg,dbr,dgr,am4/amtot,am4,
     $     abg,eccbg,abr,eccbr,agr,eccgr,abin3,eccbin3,
     $     ebin3ebin,mejecta

c asynchronize all variables:
      do i=1,n
         vx(i)=vx(i)+0.5d0*dt*vxdot(i)
         vy(i)=vy(i)+0.5d0*dt*vydot(i)
         vz(i)=vz(i)+0.5d0*dt*vzdot(i)
         if(u(i).ne.0.d0) then
            u(i)=u(i)+0.5d0*dt*udot(i)
         endif
      enddo
      t=t+0.5d0*dt

      return

      end
