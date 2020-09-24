      SUBROUTINE COMPBEST3
C     Figure out what components each particle is in
      INCLUDE 'spha.h'
      REAL*8 ENTH(NMAX),am4,rhomax
      real*8 X1,Y1,Z1,VX1,VY1,VZ1,am1
      real*8 X2,Y2,Z2,VX2,VY2,VZ2,am2
      real*8 X3,Y3,Z3,VX3,VY3,VZ3,am3
      integer icomp(nmax),ic,ncitmax,nicomp1,nicomp2,nicomp3,nicomp4,
     $     ncit,nchng
      real*8 v2i1,bi1,v2i2,bi2,v2i3,bi3
      integer ncomp1,ncomp2,ncomp3,ncomp4
      COMMON/COMPBETTERCOM3/AM1,X1,Y1,Z1,VX1,VY1,VZ1,
     $     AM2,X2,Y2,Z2,VX2,VY2,VZ2,
     $     AM3,X3,Y3,Z3,VX3,VY3,VZ3,am4,
     $     ICOMP,ncomp1,ncomp2,ncomp3
      integer icompsave(nmax),nsave
      save icompsave,nsave
      PARAMETER (NCITMAX=100)
      LOGICAL DB!,PP
      real*8 amtiny
      integer irhomax4,irhomax1,irhomax2,irhomax3
      save irhomax4,irhomax1,irhomax2,irhomax3
      real*8 ri1,ri2,ri3

      logical FIRSTT
      SAVE FIRSTT
      data FIRSTT /.TRUE./
      integer numstars,numstarslast
      common/nsl/ numstarslast
      real*8 am1last,am2last,am3last
      save am1last,am2last,am3last
      real*8 rhomax1,rhomax2,rhomax3,rhomax4
      integer swapcase,compa,compb,irhomaxtmp
      real*8 massa, massb
      real*8 amsave(nmax),musave(nmax)!,asave(nmax)
      real*8 xsave(nmax),ysave(nmax),zsave(nmax)
c      real*8 vxsave(nmax),vysave(nmax),vzsave(nmax)
      save amsave,musave,xsave,ysave,zsave!,asave!,vxsave,vysave,vzsave
c      logical alreadyused(nmax)
      integer next
      integer bestnext
      real*8 r2next,r2nextbest
      integer ccparents(4)
      integer numparents
      common/parents/numparents, ccparents
      integer i,j
      real*8 sizesq1, sizesq2

      DB=.TRUE.
C     (if true, produce detailed output on convergence)
c      PP=.true.
C     (if true, write pp files with component information)    

c      firstt=.true.
 123  continue

      if(firstt) then
c         if(n.lt.39877) then
         if(.true.) then
            ccparents(1)=cc(1)
            numparents=1
            do i=1,n
               do j=1,numparents
                  if(cc(i).eq.ccparents(j)) goto 1239
               enddo
               numparents=j
               ccparents(numparents)=cc(i)
c     1239       icomp(i)=mod(j+1,3)+1! for component 3, then 1, then 2
 1239          icomp(i)=j
            enddo
         else
            numparents=3
            ccparents(1)=cc(1)
            do i=1,n/2
               icomp(i)=1
            enddo
            ccparents(2)=cc(n/2+1)
            do i=n/2+1,2*(n/2)
               icomp(i)=2
            enddo
            ccparents(3)=cc(n)
            icomp(n)=3
         endif
c         if(numparents.le.2) then
c            numparents=numparents+1
c            do i=n/2+1,n
c               if(a(i).ne.0.d0) icomp(i)=numparents
c            enddo
c         endif

      else
         if(nsave.ne.n) then
c     Particles were thrown away during jumpahead, so icomp array needs adjusted:            
            print *,'nsave,n=',nsave,n
            do i=1,n
               r2nextbest=1.d30
               do next=1,nsave
                  if(amsave(next).eq.am(i)
     $                 .and. musave(next).eq.wmeanmolecular(i))then
                     r2next=(x(i)-xsave(next))**2+
     $                    (y(i)-ysave(next))**2+
     $                    (z(i)-zsave(next))**2
c                     r2next=abs(a(i)-asave(next))
                     if(r2next.lt.r2nextbest) then
                        r2nextbest=r2next
                        bestnext=next
                     endif
                  endif
               enddo
               icomp(i)=icompsave(bestnext)
            enddo
         endif
      endif

 124  continue
      amtiny=1.d30
      do i=1,N
         amtiny=min(amtiny,am(i))
         if(HP(i).le.-998d0 .or. a(i).eq.0.d0) then
            ENTH(I)=0
            print *,'compact object rho(',i,')=',rho(i),'->1d30'
            rho(i)=1.d30
         else
c            ENTH(I)=A(I)*rho(i)**(2.d0/3.d0)/(2.d0/3.d0)
c            ENTH(I)=A(I)   !use this if A(I) is actually specific energy

c     To a good approximation: pressure/rho = (gamma-1)*u = (2/3)u, so enthalpy = u + P/rho = (5/3)*u
c            ENTH(I)= 5d0/3d0*A(I) !use this if A(I) is actually specific energy and want true enthalpy
            ENTH(I)=0.d0
         endif
      enddo
c      write(6,*) '***ASSUMED A(I) is truly P/rho^gamma***'
      write(6,*) '**ASSUMED A(I) is specific internal energy u**'
c      write(6,*)
c      write(6,*) '***********************************************'
c      write(6,*) '**USING ONLY KINETIC AND GRAVITATIONAL ENERGY**'
c      write(6,*) '**NOT ALSO INTERNAL..........................**'
c      write(6,*) '***********************************************'
c      write(6,*)
      write(6,*) 'Smallest particle mass in system is',amtiny
      
c      print*,'old irhomaxes=',irhomax1,irhomax2,irhomax3,irhomax4
c      print *,icomp(irhomax1),icomp(irhomax2),icomp(irhomax3)

      nicomp1=0
      nicomp2=0
      nicomp3=0
      nicomp4=0

      rhomax1=0.d0
      rhomax2=0.d0
      rhomax3=0.d0
      rhomax4=0.d0
      do i=1,N
         if(icomp(i).eq.1) then
            if(rho(i).gt.rhomax1) then
               rhomax1=rho(i)
               irhomax1=i
            endif
            nicomp1=nicomp1+1
         else if(icomp(i).eq.2) then
            if(rho(i).gt.rhomax2) then
               rhomax2=rho(i)
               irhomax2=i
            endif
            nicomp2=nicomp2+1
         else if(icomp(i).eq.3) then
            if(rho(i).gt.rhomax3) then
               rhomax3=rho(i)
               irhomax3=i
            endif
            nicomp3=nicomp3+1
         else if(icomp(i).eq.4) then
            if(rho(i).gt.rhomax4) then
               rhomax4=rho(i)
               irhomax4=i
            endif
            nicomp4=nicomp4+1
         endif
      enddo
      print*,'irhomaxes=',irhomax1,irhomax2,irhomax3,irhomax4

c      if(t.gt.3) sToP

c      print*,'rhomaxes=',rho(irhomax1),rho(irhomax2),
c     $     rho(irhomax3),rho(irhomax4)
c      print*,'rhomaxes=',rhomax1,rhomax2,rhomax3,rhomax4
      print *,'number of component 1 particles in initial guess=',
     $     nicomp1
      print *,'number of component 2 particles in initial guess=',
     $     nicomp2
      print *,'number of component 3 particles in initial guess=',
     $     nicomp3
      print *,'number of component 4 particles in initial guess=',
     $     nicomp4
      
      call calccom

      IF (DB) THEN
         WRITE (6,'(6g12.6)') 'nit','nchng','m1','m2','m3'
         IF (DB) WRITE (6,'(5g13.6)') 0,nicomp1+nicomp2+nicomp3,
     $        AM1,AM2,AM3
      ENDIF

C     Iterate over component determination
      NCIT=0
      NCHNG=1
      swapcase=0
      DO WHILE (NCHNG.NE.0 .and. NCIT.LT.NCITMAX)
         NCHNG=0
c         if(t.gt.122) print *,ncit,t

         x1=x(irhomax1)
         y1=y(irhomax1)
         z1=z(irhomax1)
         vx1=vx(irhomax1)
         vy1=vy(irhomax1)
         vz1=vz(irhomax1)
         x2=x(irhomax2)
         y2=y(irhomax2)
         z2=z(irhomax2)
         vx2=vx(irhomax2)
         vy2=vy(irhomax2)
         vz2=vz(irhomax2)
         x3=x(irhomax3)
         y3=y(irhomax3)
         z3=z(irhomax3)
         vx3=vx(irhomax3)
         vy3=vy(irhomax3)
         vz3=vz(irhomax3)

         DO I=1,N
C     Assign particles to components

C     Calculate "total enthalpy" BI1 to component 1
            if(am1.gt.am(i)+0.9d0*amtiny) then
               RI1=SQRT((X(I)-X1)**2+(Y(I)-Y1)**2+(Z(I)-Z1)**2)
c               RI1=max(SQRT((X(I)-X1)**2+(Y(I)-Y1)**2+(Z(I)-Z1)**2),
c     $              hp(i))
               V2I1=(VX(I)-VX1)**2+(VY(I)-VY1)**2+(VZ(I)-VZ1)**2
               BI1=0.5d0*V2I1+ENTH(I)-AM1/RI1
            else
               BI1=1.d30
            endif

C     Calculate "total enthalpy" BI2 to component 2
            if(am2.gt.am(i)+0.9d0*amtiny) then
               RI2=SQRT((X(I)-X2)**2+(Y(I)-Y2)**2+(Z(I)-Z2)**2)
c               RI2=max(SQRT((X(I)-X2)**2+(Y(I)-Y2)**2+(Z(I)-Z2)**2),
c     $              hp(i))
               V2I2=(VX(I)-VX2)**2+(VY(I)-VY2)**2+(VZ(I)-VZ2)**2
               BI2=0.5d0*V2I2+ENTH(I)-AM2/RI2
            else
               BI2=1.d30
            endif

C     Calculate "total enthalpy" BI3 to component 3
            if(am3.gt.am(i)+0.9d0*amtiny) then
               RI3=SQRT((X(I)-X3)**2+(Y(I)-Y3)**2+(Z(I)-Z3)**2)
c               RI3=max(SQRT((X(I)-X3)**2+(Y(I)-Y3)**2+(Z(I)-Z3)**2),
c     $              hp(i))
               V2I3=(VX(I)-VX3)**2+(VY(I)-VY3)**2+(VZ(I)-VZ3)**2
               BI3=0.5d0*V2I3+ENTH(I)-AM3/RI3
            else
               BI3=1.d30
            endif

c            if(i.eq.2514) then
c               print *,'bs',i,BI1,BI2,BI3
c               print *,'bs',i,RI1,RI2,RI3
c            endif

C     Assign particle to one of 3 components
            IF ( BI1.Lt.0.d0 .AND.
     $           (BI2.ge.0.d0 .or. RI1.lt.RI2) .and.
     $           (BI3.ge.0.d0 .or. RI1.lt.RI3)
     $           ) then
C     Particle belongs to component 1
               IC=1
            ELSEIF(BI2.Lt.0.d0 .and.
     $              (BI1.ge.0.d0 .or. RI2.lt.RI1) .and.
     $              (BI3.ge.0.d0 .or. RI2.lt.RI3)
     $              )then
C     Particle belongs to component 2
               IC=2
            ELSEIF(BI3.Lt.0.d0 .and.
     $              (BI1.ge.0.d0 .or. RI3.lt.RI1) .and.
     $              (BI2.ge.0.d0 .or. RI3.lt.RI2)
     $              )then
C     Particle belongs to component 3
               IC=3
            ELSE
C     Particle is unbound (component 4)
               IC=4
c               if((icomp(i).eq.2 .or. icomp(i).eq.3) .and.
c     $              (RI2.lt.1 .or. RI3.lt.1))then
c                  print *,'being put in 4 ... particle',i,icomp(i)
c                  print *,BI1,BI2,BI3
c                  print *,RI1,RI2,RI3
c                  stop
c                  endif
            ENDIF

c 9/3/04: Don't let point masses have their component change
            IF(HP(I).LT.-998.d0 .or. a(i).eq.0.d0) IC=ICOMP(I)

c            if(icomp(i).eq.2) stop
            
            IF (IC.NE.ICOMP(I)) THEN
C     (particle assignment has changed)
C     Adjust components
C     Remove from old component

               if(AM(I).gt.1d0) then
c               if(.true.) then
                  print *,'i=',I,'m=',AM(I),
     $                 'from',ICOMP(I),'to',IC
                  print *,'BIs=',BI1,BI2
                  print *,'V2s=',V2I1, V2I2
                  print *,'pots=',-(AM1-AM(I))/RI1,-(AM2-AM(I))/RI2
                  print *,'rs=',RI1,RI2
               endif

               IF (ICOMP(I).EQ.1) THEN
c                  AM1=AM1-AM(I)
                  if(i.eq.irhomax1) then
                     rhomax=0.d0
                     irhomax1=0
                     do j=1,N
                        if(j.ne.i .and. icomp(j).eq.1 .and.
     $                       rho(j).gt.rhomax) then
                           rhomax=rho(j)
                           irhomax1=j
                        endif
                     enddo
c                     print*,'Removing',i,'from comp 1; new=',irhomax1
                  endif
               ELSE IF (ICOMP(I).EQ.2) THEN
c                  AM2=AM2-AM(I)
                  if(i.eq.irhomax2) then
                     rhomax=0.d0
                     irhomax2=0
                     do j=1,N
                        if(j.ne.i .and. icomp(j).eq.2 .and.
     $                       rho(j).gt.rhomax) then
                           rhomax=rho(j)
                           irhomax2=j
                        endif
                     enddo
c                     print*,'Removing',i,'from comp 2; new=',irhomax2
                  endif
               ELSE IF (ICOMP(I).EQ.3) THEN
c                  AM3=AM3-AM(I)
                  if(i.eq.irhomax3) then
                     rhomax=0.d0
                     irhomax3=0
                     do j=1,N
                        if(j.ne.i .and. icomp(j).eq.3 .and.
     $                       rho(j).gt.rhomax) then
                           rhomax=rho(j)
                           irhomax3=j
                        endif
                     enddo
c                     print*,'Removing',i,'from comp 3; new=',irhomax3
                  endif
               ENDIF

C     Add to new component
               IF (IC.EQ.1) THEN
                  if(rho(i).gt.rho(irhomax1)) then
                     if(icomp(i).ne.4)
     $                    print *,'In comp 1: Replacing',irhomax1,
     $                    'with',i,'from comp',icomp(i)
                     irhomax1=i
                     if(icomp(i).ne.4) then
                        if(swapcase.ne.0 .and.
     $                       swapcase.ne.icomp(i)+ic)then
                           print *,'double swap???'
c                           stop
                        endif
                        swapcase=icomp(i)+ic
                        print *,'swapcase',i,icomp(i),ic
                     endif
                  endif
               ELSE IF (IC.EQ.2) THEN
                  if(rho(i).gt.rho(irhomax2)) then
                     if(icomp(i).ne.4)
     $                    print *,'In comp 2: Replacing',irhomax2,
     $                    'with',i,'from comp',icomp(i)
                     irhomax2=i
                     if(icomp(i).ne.4) then
                        if(swapcase.ne.0 .and.
     $                       swapcase.ne.icomp(i)+ic)then
                           print *,'double swap???'
c                           stop
                        endif
                        swapcase=icomp(i)+ic
                        print *,'swapcase',i,icomp(i),ic
                     endif
                  endif
               ELSE IF (IC.EQ.3) THEN
                  if(rho(i).gt.rho(irhomax3)) then
                     if(icomp(i).ne.4)
     $                    print *,'In comp 3: Replacing',irhomax3,
     $                    'with',i,'from comp',icomp(i)
                     irhomax3=i
                     if(icomp(i).ne.4) then
                        if(swapcase.ne.0 .and.
     $                       swapcase.ne.icomp(i)+ic)then
                           print *,'double swap???'
c                           stop
                        endif
                        swapcase=icomp(i)+ic
                        print *,'swapcase',i,icomp(i),ic
                     endif
                  endif
               ENDIF

c               print *,i,'is being moved from',icomp(i),'to',ic
c               print *,x(i),y(i),z(i)

               ICOMP(I)=IC
               NCHNG=NCHNG+1
            ENDIF
         ENDDO

         call calccom

         NCIT=NCIT+1
         
         IF (DB) WRITE (6,'(5g13.6)') NCIT,NCHNG,AM1,AM2,AM3

         IF (NCIT.GE.NCITMAX) then
            write(6,*) 'COMPBEST: NO CONVERGENCE ???'
         endif

      ENDDO

      if(swapcase.ne.0 .and. .false.) then

 1011    continue
         do i=1,n
            if(a(i).gt.0.d0) icomp(i)=1
         enddo
         
c         if(swapcase.eq.3) then
            compa=1
            compb=2
            irhomaxtmp=irhomax1
            irhomax1=irhomax2
            irhomax2=irhomaxtmp
c         elseif(swapcase.eq.4)then
c            compa=1
c            compb=3
c            irhomaxtmp=irhomax1
c            irhomax1=irhomax3
c            irhomax3=irhomaxtmp
c         elseif(swapcase.eq.5)then
c            compa=2
c            compb=3
c            irhomaxtmp=irhomax2
c            irhomax2=irhomax3
c            irhomax3=irhomaxtmp
c         else
c            print *,'swapcase problem'
c            stop
c         endif
         print *,'time to merge the stars at t=',t
         if(compa.eq.1) massa=am1last
         if(compa.eq.2) massa=am2last
         if(compa.eq.3) massa=am3last
         if(compb.eq.1) massb=am1last
         if(compb.eq.2) massb=am2last
         if(compb.eq.3) massb=am3last
         print *,'Will merge components to recover...'
         print *,'last masses=',massa,massb
         if(massa.gt.massb)then
            print *,'switching all particles from',compb,'to',compa
            do i=1,n
               if(icomp(i).eq.compb)icomp(i)=compa
            enddo
         else
            print *,'switching all particles from',compa,'to',compb
            do i=1,n
               if(icomp(i).eq.compa)icomp(i)=compb
            enddo
         endif
         print *,'Let us try again'
         firstt=.false.
         goto 124
      endif

      ncomp1=0
      ncomp2=0
      ncomp3=0
      ncomp4=0
      AM4=0.d0
      sizesq1=0
      sizesq2=0
      DO I=1,N
         IF (ICOMP(I).EQ.1) THEN
            sizesq1=sizesq1+am(i)*((x(i)-x1)**2+(y(i)-y1)**2
     $           +(z(i)-z1)**2)
            ncomp1=ncomp1+1
         ELSE IF (ICOMP(I).EQ.2) THEN
            sizesq2=sizesq2+am(i)*((x(i)-x2)**2+(y(i)-y2)**2
     $           +(z(i)-z2)**2)
            ncomp2=ncomp2+1
         ELSE IF (ICOMP(I).EQ.3) THEN
            ncomp3=ncomp3+1
         ELSE IF (ICOMP(I).EQ.4) THEN
            AM4=AM4+AM(I)
            ncomp4=ncomp4+1
         ENDIF
      ENDDO
      sizesq1=sizesq1/am1
      sizesq2=sizesq2/am2

      numstars=min(ncomp1,1)+min(ncomp2,1)+min(ncomp3,1)

      if (ncomp1.eq.0 .and. ncomp2.eq.0 .and. ncomp3.eq.0) then
         print *,'WARNING: no stars found!'
         if(firstt.eqv. .false.) then
            print *,'Trying a different initial guess:'
            print *,'IF STARS ARE FOUND, THERE IS NO GUARANTEE THAT',
     $           'THEY WILL BE THE SAME COMPONENT NUMBERS AS IN THE',
     $           'PREVIOUS ITERATION'
            firstt=.true.
            goto 123
         else
            stop
         endif
      endif

      numstarslast=numstars
      am1last=am1
      am2last=am2
      am3last=am3

      write(6,'(a,7a10)') '       ','mass   ','x   ','y   ','z   ',
     $     'vx   ','vy   ','vz   '

      if(am1.gt.0d0) then
         write(6,'(a,9g10.3)') 'star 1:',am1,x1,y1,z1,vx1,vy1,vz1,
     $        (x1**2+y1**2+z1**2)**0.5d0
      endif
      if(am2.gt.0d0) then
         write(6,'(a,9g10.3)') 'star 2:',am2,x2,y2,z2,vx2,vy2,vz2,
     $        (x2**2+y2**2+z2**2)**0.5d0
      endif
      if(am3.gt.0d0) then
         write(6,'(a,9g10.3)') 'star 3:',am3,x3,y3,z3,vx3,vy3,vz3,
     $        (x3**2+y3**2+z3**2)**0.5d0
      endif

c      print *,'merge?',(x1-x2)**2+(y1-y2)**2+(z1-z2)**2,
c     $     (am1+am2)**2,
c     $     (vx1-vx2)**2+(vy1-vy2)**2+(vz1-vz2)**2,
c     $     0.5d0*(am1+am2)/((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)**0.5d0,
c     $     sizesq1,sizesq2
c
cc      if((x1-x2)**2+(y1-y2)**2+(z1-z2)**2.lt.2*(am1+am2)**2 .and.
c      if((x1-x2)**2+(y1-y2)**2+(z1-z2)**2.lt.sizesq1+sizesq2 .and.
c     $     (vx1-vx2)**2+(vy1-vy2)**2+(vz1-vz2)**2.lt.
c     $     0.5d0*(am1+am2)/((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)**0.5d0
c     $     .and. am1.lt.1000 .and. am2.lt.1000)then
c         print *,'These stars are close enough to have merged.',
c     $        ((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)**0.5d0,
c     $        (am1+am2),
c     $        (vx1-vx2)**2+(vy1-vy2)**2+(vz1-vz2)**2
cc         sToP
c         goto 1011
c      endif



      firstt=.false.
      nsave=n
      do i=1,n
         amsave(i)=am(i)
         musave(i)=wmeanmolecular(i)
         icompsave(i)=icomp(i)
         xsave(i)=x(i)
         ysave(i)=y(i)
         zsave(i)=z(i)
      enddo

      RETURN
      END
