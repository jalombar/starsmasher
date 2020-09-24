      SUBROUTINE robusttrajectories
      INCLUDE 'spha.h'
      real*8 ax1,ay1,az1,ax2,ay2,az2,ax3,ay3,az3
      real*8 X1,Y1,Z1,VX1,VY1,VZ1,am1
      real*8 X2,Y2,Z2,VX2,VY2,VZ2,am2
      real*8 X3,Y3,Z3,VX3,VY3,VZ3,am3
      real*8 am4
      real*8 axr,ayr,azr
      real*8 axg,ayg,azg
      real*8 axb,ayb,azb
      integer icomp(nmax)
      COMMON/COMPBETTERCOM3/AM1,X1,Y1,Z1,VX1,VY1,VZ1,
     $     AM2,X2,Y2,Z2,VX2,VY2,VZ2,
     $     AM3,X3,Y3,Z3,VX3,VY3,VZ3,am4,
     $     ICOMP
      real*8 amr,amg,amb,xr,yr,zr,xg,yg,zg,xb,yb,zb
      save xr,yr,zr,xg,yg,zg,xb,yb,zb
      real*8 vxr,vyr,vzr,vxg,vyg,vzg,vxb,vyb,vzb
      save vxr,vyr,vzr,vxg,vyg,vzg,vxb,vyb,vzb
      save amr,amg,amb
      logical ft
      SAVE FT
      data FT /.TRUE./
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
      real*8 RIcom,V2Icom,BIcom
      integer blue,green,red
      save blue,green,red
c      character*1 stara,starb,starc
      integer istara,istarb,istarc
      integer starnum
      real*8 rother
      character*9 starname(4)
      save starname
      integer numstars,lastnumstars
      data lastnumstars /3/
      save lastnumstars
      real*8 lastamb,lastamg,lastamr
      save lastamb,lastamg,lastamr
      integer absorber,absorbee
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
      real*8 eejecta,ekinejecta,eintejecta,epotejecta
      real*8 etotcheck
      integer numspeeds
      real*8 speed1, speed2, speedunit, eunit
      integer ndumpflag
      common/dumpflag/ ndumpflag
      integer i
      real*8 dotproduct,amtot,vxdrft,vydrft,vzdrft,speedbeta,sep,
     $     speedalphainf,speedbetainf,
     $     actualspeedalphainf,actualspeedbetainf
      real*8 speedcom,speed3rd,speedcominf,speed3rdinf
      real*8 speedalpha,actualspeedcominf,actualspeed3rdinf
      real*8 thetaindegrees
      real*8 xrel,yrel,zrel
      real*8 vxrel,vyrel,vzrel
      real*8 spinx,spiny,spinz
      real*8 rhomaxic
      real*8 spintotal(3),momentofinertia(3)
      integer ic,icenter(3)
      integer ii
      real*8 energy(3),partpot,dist2,userfunc9b
      real*8 energykin(3),energypot(3),energyint(3)
      real*8 xcenter(3),ycenter(3),zcenter(3)
      real*8 vxcenter(3),vycenter(3),vzcenter(3)

      real*8 rmin
      CHARACTER*11 FNAME,FNAME2
      real*8 vejavg,v2ejavg
      real*8 vunit
      real*8 vinkms,binwidth,vinkmsmax
      parameter(vinkmsmax=4000.d0)
      integer numbins
      parameter(numbins=100)
      integer binnumber
      real*8 binmass(numbins)      
      parameter(binwidth=vinkmsmax/numbins)
      logical binvelocitydata
      real*8 ambsome,amgsome
      real*8 amblast,amglast
      real*8 xbsome,ybsome,zbsome,vxbsome,vybsome,vzbsome
      real*8 xgsome,ygsome,zgsome,vxgsome,vygsome,vzgsome
      real*8 ambcircumbinary,amgcircumbinary
      integer counter

      inquire(file='input.3s',exist=bps)
c If the file input.3s exists in this directory, then make bps=.false.
c If the file input.3s does not exist, then assume binary+single
      bps= .not. bps

      print *,'************************* FIRSTTIME=',FT

c     Make checkifswapped=.false. if the particle indices may change from
c     one outfile to the next.... (this happens in some of evghenii's gpu runs)
      checkifswapped=.true.
 88   continue

      ndumpflag=-1
      call COMPBEST3

      xcenter(1)=x1
      ycenter(1)=y1
      zcenter(1)=z1
      vxcenter(1)=vx1
      vycenter(1)=vy1
      vzcenter(1)=vz1
      xcenter(2)=x2
      ycenter(2)=y2
      zcenter(2)=z2
      vxcenter(2)=vx2
      vycenter(2)=vy2
      vzcenter(2)=vz2
      xcenter(3)=x3
      ycenter(3)=y3
      zcenter(3)=z3
      vxcenter(3)=vx3
      vycenter(3)=vy3
      vzcenter(3)=vz3

      do ic=1,3
         rhomaxic=0.d0
         do i=1,n
            if(icomp(i).eq.ic) then
               if(rho(i).gt.rhomaxic) then
                  icenter(ic)=i
                  rhomaxic=rho(i)
               endif
               if(a(i).eq.0) then
                  icenter(ic)=i
                  rhomaxic=1.d30
               endif
            endif
         enddo
         print *,'particle at center of comp',ic,'=',icenter(ic)
         spinx=0.d0
         spiny=0.d0
         spinz=0.d0
         do i=1,n
            if(icomp(i).eq.ic) then
               xrel=x(i)-xcenter(ic)
               yrel=y(i)-ycenter(ic)
               zrel=z(i)-zcenter(ic)
               vxrel=vx(i)-vx(icenter(ic))
               vyrel=vy(i)-vy(icenter(ic))
               vzrel=vz(i)-vz(icenter(ic))
               spinx=spinx+am(i)*(yrel*vzrel-zrel*vyrel)
               spiny=spiny+am(i)*(zrel*vxrel-xrel*vzrel)
               spinz=spinz+am(i)*(xrel*vyrel-yrel*vxrel)
            endif
         enddo
         spintotal(ic)=(spinx**2+spiny**2+spinz**2)**0.5d0
      enddo

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

      numstars=0
      if(am1.gt.0.d0) numstars=numstars+1
      if(am2.gt.0.d0) numstars=numstars+1
      if(am3.gt.0.d0) numstars=numstars+1
      print *,'number of stars=',numstars

      runit=6.9599d10           ! number of cm in the unit of length.  Use 6.9599d10 if want solar radius.
      munit=1.9891d33           ! number of g in unit of mass.  Use 1.9891d33 if want solar mass.

      if(ft) then
         speedunit=(gravconst*munit/runit)**0.5d0*1.d-5
         eunit=gravconst*munit**2/runit*1.d-48
         print *,'ASSUMING unit of length runit (in cm)=',runit
         print *,'ASSUMING unit of mass runit (in g)=',munit

         print *,'fyi, speed and energy units (in km/s and 10^48 erg)=',
     $        speedunit,eunit

         starname(1)='1'
         starname(2)='2'
         starname(3)='3'
c     first time through just make most massive component red(r), 2nd most massive green (g), and third most massive blue (b)

         print *,'This is the first time in this subroutine'
         print *,'making most massive star blue, next green, next red'
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
      endif

      if(red.eq.1)then
         print *,'component 1 is red'
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
         print *,'component 2 is red'
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
         print *,'component 3 is red'
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
         print *,'component 1 is green'
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
         print *,'component 2 is green'
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
         print *,'component 3 is green'
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
         print *,'component 1 is blue'
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
         print *,'component 2 is blue'
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
         print *,'component 3 is blue'
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
         print *,'star not assigned properly',red,green,blue
         stop
      endif

      dgr=((xg-xr)**2+(yg-yr)**2+(zg-zr)**2)**0.5d0
      dbr=((xr-xb)**2+(yr-yb)**2+(zr-zb)**2)**0.5d0
      dbg=((xg-xb)**2+(yg-yb)**2+(zg-zb)**2)**0.5d0

      ambcircumbinary=amb
      amgcircumbinary=amg
      counter=0
 198  continue
      counter=counter+1
      amblast=amb
      amglast=amg
      amtot=0.d0
      ambsome=0d0
      amgsome=0d0
      xbsome=0d0
      ybsome=0d0
      zbsome=0d0
      vxbsome=0d0
      vybsome=0d0
      vzbsome=0d0
      xgsome=0d0
      ygsome=0d0
      zgsome=0d0
      vxgsome=0d0
      vygsome=0d0
      vzgsome=0d0
      do i=1,N
         amtot=amtot+am(i)
         if(icomp(i).eq.blue .and.
     $        (x(i)-x(icenter(blue)))**2+
     $        (y(i)-y(icenter(blue)))**2+
     $        (z(i)-z(icenter(blue)))**2.le.dbg**2) then
            ambsome=ambsome+am(i)
            xbsome=xbsome+am(i)*x(i)
            ybsome=ybsome+am(i)*y(i)
            zbsome=zbsome+am(i)*z(i)
            vxbsome=vxbsome+am(i)*vx(i)
            vybsome=vybsome+am(i)*vy(i)
            vzbsome=vzbsome+am(i)*vz(i)
         else if (icomp(i).eq.green .and.
     $           (x(i)-x(icenter(green)))**2+
     $           (y(i)-y(icenter(green)))**2+
     $           (z(i)-z(icenter(green)))**2.le.dbg**2) then
            amgsome=amgsome+am(i)
            xgsome=xgsome+am(i)*x(i)
            ygsome=ygsome+am(i)*y(i)
            zgsome=zgsome+am(i)*z(i)
            vxgsome=vxgsome+am(i)*vx(i)
            vygsome=vygsome+am(i)*vy(i)
            vzgsome=vzgsome+am(i)*vz(i)
         endif
      enddo

c      xb=xbsome/ambsome
c      yb=ybsome/ambsome
c      zb=zbsome/ambsome
c      vxb=vxbsome/ambsome
c      vyb=vybsome/ambsome
c      vzb=vzbsome/ambsome
c      xg=xgsome/amgsome
c      yg=ygsome/amgsome
c      zg=zgsome/amgsome
c      vxg=vxgsome/amgsome
c      vyg=vygsome/amgsome
c      vzg=vzgsome/amgsome
      xb=x(icenter(blue))
      yb=y(icenter(blue))
      zb=z(icenter(blue))
      vxb=vx(icenter(blue))
      vyb=vy(icenter(blue))
      vzb=vz(icenter(blue))
      xg=x(icenter(green))
      yg=y(icenter(green))
      zg=z(icenter(green))
      vxg=vx(icenter(green))
      vyg=vy(icenter(green))
      vzg=vz(icenter(green))
      amb=ambsome
      amg=amgsome
      dbg=((xg-xb)**2+(yg-yb)**2+(zg-zb)**2)**0.5d0
      if((ambsome.ne.amblast .or. amgsome.ne.amglast) .and.
     $     counter.le.100) then
         goto 198
      endif
      ambcircumbinary=ambcircumbinary-ambsome
      amgcircumbinary=amgcircumbinary-amgsome

c      print *,amb,ambsome,xb,xbsome,yb,ybsome,zb,zbsome,
c     $     vxb,vxbsome,vyb,vybsome,vzb,vzbsome
c      print *
c      print *,amg,amgsome,xg,xgsome,yg,ygsome,zg,zgsome,
c     $     vxg,vxgsome,vyg,vygsome,vzg,vzgsome
c      sToP 
c      print *,'counter=',counter

      abinary=.false.
      compinabinary(1)=.false.
      compinabinary(2)=.false.
      compinabinary(3)=.false.
      compinabinary(4)=.false.
      a12=1.d30
      if(amr.gt.0.d0 .and. amg.gt.0.d0) then
         call elements(amr,amg,xr-xg,yr-yg,zr-zg,
     $        vxr-vxg,vyr-vyg,vzr-vzg,
     $        agr,eccgr,lxgr,lygr,lzgr)

         if(eccgr.lt.1.d0) then
            dotproduct=(xg-xr)*(vxg-vxr)+
     $           (yg-yr)*(vyg-vyr)+
     $           (zg-zr)*(vzg-vzr)
            
            write(6,'(a,3g13.6)')
     $           'green and red are a binary with a,e,dotproduct=',
     $           agr,eccgr,dotproduct

            if((agr.lt.min((amg+amr)*0.04d0,10.d0)
     $           .and. eccgr.lt.1.d0) .or.
     $           (abs(dotproduct).lt.0.1d0 .and. eccgr.gt.0.98d0 .and.
     $           eccgr.lt.1.d0 .and. agr.lt.min(amg+amr,10.d0) ))then
               print *,'green and red seem to have merged'
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
            a12=agr
            r12=dgr
            lx12=lxgr
            ly12=lygr
            lz12=lzgr
c            if(ft .and. bps)then
c               starname(1)='i'  ! The numbers 1,2,3 correspond to a ranking of
c               starname(2)='p'  ! the initial masses, not to component numbers
c               starname(3)='s'
c               print*,'starnames=',(trim(starname(starnum)),starnum=1,3)
c            endif
         endif
      else
         agr=0.d0
         eccgr=-1.d0
      endif
      if(amr.gt.0.d0 .and. amb.gt.0.d0) then
         call elements(amr,amb,xr-xb,yr-yb,zr-zb,
     $        vxr-vxb,vyr-vyb,vzr-vzb,
     $        abr,eccbr,lxbr,lybr,lzbr)
         if(eccbr.lt.1.d0 .and. (.not. abinary .or.
     $        abr*(1.d0+eccbr).le.a12*(1.d0+ecc12)))then 
            dotproduct=(xb-xr)*(vxb-vxr)+
     $           (yb-yr)*(vyb-vyr)+
     $           (zb-zr)*(vzb-vzr)
            write(6,'(a,3g13.6)')
     $           'blue and red are a binary with a,e,dotproduct=',
     $           abr,eccbr,dotproduct

            if((abr.lt.min((amb+amr)*0.04d0,10.d0)
     $           .and. eccbr.lt.1.d0) .or.
     $           (abs(dotproduct).lt.0.1d0 .and. eccbr.gt.0.98d0 .and.
     $           eccbr.lt.1.d0 .and. abr.lt.min(amb+amr,10.d0) ))then
               print *,'blue and red seem to have merged'
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
            if(abinary .and. abr*(1.d0+eccbr).le.a12*(1.d0+ecc12))then
               print *,'so green star not a binary member after all'
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

         if(eccbg.lt.1.d0 .and. (.not. abinary .or.
     $        abg*(1.d0+eccbg).le.a12*(1.d0+ecc12)))then 
            dotproduct=(xb-xg)*(vxb-vxg)+
     $           (yb-yg)*(vyb-vyg)+
     $           (zb-zg)*(vzb-vzg)
            
            write(6,'(a,3g13.6)')
     $           'blue and green are a binary with a,e,dotproduct=',
     $           abg,eccbg,dotproduct

            if((abg.lt.min((amb+amg)*0.04d0,10.d0)
     $           .and. eccbg.lt.1.d0) .or.
     $           (abs(dotproduct).lt.0.1d0 .and. eccbg.gt.0.98d0 .and.
     $           eccbg.lt.1.d0 .and. abg.lt.min(amb+amg,10.d0) ))then
               print *,'blue and green seem to have merged'
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
            if(abinary .and. abg*(1.d0+eccbg).le.a12*(1.d0+ecc12))then
               print *,'so red star not a binary member after all'
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

      if(numstars.ne.lastnumstars) then
         print *,lastnumstars-numstars+1,'stars have merged'
         if(lastnumstars-numstars.gt.1) stop
         absorbee=3
         if(amb.eq.0.d0 .and. lastamb.gt.0.d0) absorbee=1
         if(amg.eq.0.d0 .and. lastamg.gt.0.d0) absorbee=2
         if(amr.eq.0.d0 .and. lastamr.gt.0.d0) absorbee=3
         deltaamb=amb-lastamb
         deltaamg=amg-lastamg
         deltaamr=amr-lastamr
         print *,'deltas=',
     $        deltaamb,deltaamg,deltaamr
         if(deltaamb.gt.deltaamg .and. deltaamb.gt.deltaamr) then
            absorber=1
         else if(deltaamg.gt.deltaamb .and. deltaamg.gt.deltaamr) then
            absorber=2
         else if(deltaamr.gt.deltaamg .and. deltaamr.gt.deltaamb) then
            absorber=3
         else
            absorber=4
            print *,'cannot determine absorber',
     $           deltaamb,deltaamg,deltaamr
         endif
         lastnumstars=numstars
         if(absorber.lt.4)then
            write(6,*)'star ',trim(starname(absorber)),
     $           ' has absorbed star ',trim(starname(absorbee)),
     $           absorber,absorbee
c            starname(absorber)=
c     $           '{'//trim(starname(absorber))//','
c     $           //trim(starname(absorbee))//'}'
         else
            write(6,*)
     $           'infinity has absorbed star ',trim(starname(absorbee)),
     $           absorber,absorbee
            starname(absorber)=
     $           '{infinity,'
     $           //trim(starname(absorbee))//'}'
         endif
      endif
      lastamb=amb
      lastamg=amg
      lastamr=amr

      if(abinary) then
c     There is a binary:
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
               print *,'stara should not be red'
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
c     There are just individual stars:
         write(6,'(a,3g13.6)')
     $        'There are just individual stars with masses',amb,amg,amr
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
 
         print *,'number of individual stars=',starnum
         
      endif
         
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c Here we'll determine how much of the mass in component 4 is
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
      do i=1,N
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
         print *,'total mass in binary or multiple=',mcom
         xcom=xcom/mcom
         ycom=ycom/mcom
         zcom=zcom/mcom
         vxcom=vxcom/mcom
         vycom=vycom/mcom
         vzcom=vzcom/mcom
         write(6,'(a,3g13.6)')
     $        'center of binary or multiple=',xcom,ycom,zcom
         write(6,'(a,3g13.6)')
     $        'velocity of binary or multiple=',vxcom,vycom,vzcom

         ekinejecta=0.d0
         eintejecta=0.d0
         epotejecta=0.d0
         etotcheck=0.d0
         do i=1,N
            if(icomp(i).eq.4) then
               RIcom=SQRT((X(I)-Xcom)**2+(Y(I)-Ycom)**2+(Z(I)-Zcom)**2)
               V2Icom=(VX(I)-VXcom)**2+(VY(I)-VYcom)**2+(VZ(I)-VZcom)**2
c               BIcom=0.5d0*V2Icom+A(I)-Mcom/RIcom
               BIcom=0.5d0*V2Icom-Mcom/RIcom
c               BIcom=0.5d0*V2Icom-Mcom/RIcom
               If(BIcom.le.0.d0) then
                  mce=mce+am(i)
               else
                  mejecta=mejecta+am(i)
                  ekinejecta=ekinejecta+0.5d0*am(i)*
     $                 (vx(i)**2+vy(i)**2+vz(i)**2)
                  eintejecta=eintejecta+am(i)*a(i)
                  epotejecta=epotejecta+am(i)*grpot(i)
               endif
            endif
            etotcheck=etotcheck+am(i)*a(i)+0.5d0*am(i)*
     $           (vx(i)**2+vy(i)**2+vz(i)**2)+0.5d0*am(i)*grpot(i)
         enddo
         eejecta=ekinejecta+eintejecta+epotejecta

         speedcom=(vxcom**2+vycom**2+vzcom**2)**0.5d0
         if(numstars.eq.2) then
            speedcom=(vxcom**2+vycom**2+vzcom**2)**0.5d0
            print *,'speed of binary:',
     $           speedcom,speedcom*mcom
         elseif(numstars.eq.3 .and.m3rd.eq.0.d0) then
            speedcom=(vxcom**2+vycom**2+vzcom**2)**0.5d0
            print *,'speed of triple:',
     $           speedcom,speedcom*mcom
         elseif(numstars.eq.3) then
c     There must be a binary AND a third star
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
            print *,'speed of binary:',
     $           speedcom,speedcom*mcom ! This is in c.o.m. frame OF THE BINARY + SINGLE
            print *,'speed of 3rd star:',
     $           speed3rd,speed3rd*m3rd ! This is in c.o.m. frame OF THE BINARY + SINGLE

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
            print *,'at infinity, speed of binary:',
     $           speedcominf,speedcominf*mcom ! This is in c.o.m. frame OF THE BINARY + SINGLE
            print *,'at infinity, speed of 3rd star:',
     $           speed3rdinf,speed3rdinf*m3rd ! This is in c.o.m. frame OF THE BINARY + SINGLE

            write(6,'(a,3g13.6)')
     $           'drift velocity components:',vxdrft,vydrft,vzdrft
c     Shift back to original inertial reference frame:
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
            write(6,'(a,3g13.6)')
     $           'at infinity, actual speed of binary:',
     $           actualspeedcominf,actualspeedcominf*mcom
            write(6,'(a,3g13.6)')
     $           'at infinity, actual speed of 3rd star:',
     $           actualspeed3rdinf,actualspeed3rdinf*m3rd

         endif
      else
         write(6,*) 'There are no binaries or other multiples'
         mejecta=am4
         ekinejecta=0.d0
         eintejecta=0.d0
         epotejecta=0.d0
         etotcheck=0.d0
         do i=1,N
            if(icomp(i).eq.4) then
               ekinejecta=ekinejecta+0.5d0*am(i)*
     $              (vx(i)**2+vy(i)**2+vz(i)**2)
               eintejecta=eintejecta+am(i)*a(i)
               epotejecta=epotejecta+am(i)*grpot(i)
            endif
            etotcheck=etotcheck+am(i)*a(i)+0.5d0*am(i)*
     $           (vx(i)**2+vy(i)**2+vz(i)**2)+0.5d0*am(i)*grpot(i)
         enddo
         eejecta=ekinejecta+eintejecta+epotejecta

      endif

      vunit=sqrt(gravconst*munit/runit)
      print *,'vunit=',vunit
         
      binvelocitydata=.false.
      if(binvelocitydata) then
         IF(NOUT.LE.9999) THEN
            WRITE (FNAME,1102) NOUT
         else
            WRITE (FNAME,1103) NOUT
         ENDIF
 1102    FORMAT ('vd',I4.4,'.sph')
 1103    FORMAT ('vd',I5.5,'.sph')      
         OPEN (33, FILE=FNAME)
         
         IF(NOUT.LE.9999) THEN
            WRITE (FNAME2,1104) NOUT
         else
            WRITE (FNAME2,1105) NOUT
         ENDIF
 1104    FORMAT ('bd',I4.4,'.sph')
 1105    FORMAT ('bd',I5.5,'.sph')      
         OPEN (42, FILE=FNAME2)

         mejecta=0.d0
         vejavg=0.d0
         v2ejavg=0.d0
         
         do binnumber=1,numbins
            binmass(binnumber)=0.d0
         enddo
         
         do i=1,N
            if(icomp(i).eq.4) then
               RIcom=SQRT((X(I)-Xcom)**2+(Y(I)-Ycom)**2+(Z(I)-Zcom)**2)
               V2Icom=(VX(I)-VXcom)**2+(VY(I)-VYcom)**2+(VZ(I)-VZcom)**2
c     BIcom=0.5d0*V2Icom+U(I)+grpot(i)
c               BIcom=0.5d0*V2Icom+A(I)-Mcom/RIcom
               BIcom=0.5d0*V2Icom-Mcom/RIcom
c     BIcom=0.5d0*V2Icom-Mcom/RIcom
               If(BIcom.gt.0.d0) then
                  mejecta=mejecta+am(i)
                  vejavg=vejavg+am(i)*V2Icom**0.5d0
                  v2ejavg=v2ejavg+am(i)*V2Icom
                  vinkms=V2Icom**0.5d0*vunit*1.d-5 ! converted to km/s
                  
                  WRITE (33, 98) i,vinkms, RIcom
 98               FORMAT(99g13.5)
                  
                  binnumber=int(vinkms/binwidth)+1
                  binmass(binnumber)=binmass(binnumber)+am(i)
                  
               endif
            endif
         enddo
         CLOSE (33)
         
         vejavg=vejavg/mejecta
         v2ejavg=v2ejavg/mejecta
         
         write(42,*) binwidth, vejavg*vunit*1d-5
         do binnumber=1,numbins
            write(42,*)binnumber,binmass(binnumber)
         enddo
         close(42)
         
         print *,'average speed=',vejavg*vunit,
     $        (v2ejavg-vejavg**2)**0.5d0*vunit,
     $        vinkmsmax
         
      endif
         
      print *,'mcom=',mcom,'m3rd=',m3rd,'mcom+m3rd=',mcom+m3rd,
     $     'amtot=',amtot,'mejecta=',mejecta,'am4=',am4
         
      print *,'fractional mass loss:',mejecta/amtot
      write(6,'(a,4g13.5)') 'energies (E,T,U,W) in ejecta:',eejecta,
     $     ekinejecta,eintejecta,epotejecta
      
      print *,'etotcheck=',etotcheck,t


      if(abs(am4-mce-mejecta).gt.1.d-5) then
         print *,'conservation of mass problem:',am4,mce,mejecta
         stop
      endif
      write(6,'(3(a,g13.6))')'m4=mce+mejecta',am4,'=',mce,'+',mejecta

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     assume a "bg" binary exists and calculate a period for it
c     from Kepler III.  Assume units are solar units and output
c     period in days:

      print '(a,9f10.2)','Orbital period (days)=',
     $     2*pi*(abg*runit)**1.5d0*
     $     (gravconst*(AMb+AMg)*munit)**(-0.5d0)/
     $     (60*60*24),2*pi*abg**1.5d0*(AMb+AMg)**(-0.5d0)*0.018445d0

      print '(a,9f10.2)','Orbital period (code units)=',
     $     2*pi*abg**1.5d0*(AMb+AMg)**(-0.5d0)

      write(24,'(d12.5,99g13.5)') T,AMb,AMg,
     $     Xb,Yb,Zb,Xg,Yg,Zg,
     $     dbg,am4,
     $     abg,eccbg,
     $     mejecta,spintotal(blue),spintotal(green),
     $     2*pi*abg**1.5d0*(AMb+AMg)**(-0.5d0)

      return

      END
