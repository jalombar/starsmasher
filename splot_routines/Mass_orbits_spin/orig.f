      SUBROUTINE READIT(NNIT,IFORM)
      INCLUDE 'spha.h'                                     
      CHARACTER*13 FNAME
      logical fileexists
      real*8 divv(nmax)
      common/commdivv/divv
      integer nnit,iform,i,nchk
      real*8 sep0,tjumpahead
      real*8 theta, xold,yold,vxold,vyold

      IF(NNIT.LE.9999) THEN
c         WRITE (FNAME,1101) NNIT
c      ELSE if(NNIT.le.9999) then
         WRITE (FNAME,1102) NNIT
      ELSE if(NNIT.le.99999) then
         WRITE (FNAME,2102) NNIT
      ELSE
         WRITE (FNAME,3102) NNIT
      ENDIF
 1101 FORMAT ('out',I3.3,'.sph')
 1102 FORMAT ('out',I4.4,'.sph')
 2102 FORMAT ('out',I5.5,'.sph')
 3102 FORMAT ('out',I6.6,'.sph')
 3111 continue
      inquire(FILE=FNAME,EXIST=fileexists)
      if(fileexists)then
         OPEN (12,FILE=FNAME,FORM='UNFORMATTED')
      else
         print*,'Ran out of output files'
         stop
c         iform=-1
c         return
cc         call sleep(10)
cc         goto 3111
      endif
      
      IF (IFORM.EQ.4) THEN
C    (The following READ sequence must match exactly the WRITE sequence
C     used in spha)
         READ(12) N,NNOPT,HMIN,HMAX,SEP0,
     $        TF,DTOUT,NOUT,NIT,T,
     $        NAV,ALPHA,BETA,tjumpahead,
     $        NGR,
     $        NRELAX,TRELAX,DT,OMEGA2
c         if(NNIT.eq.0) then

         if(N.gt.NMAX) then
            print *,'Increase NMAX in *.h files and recompile',N
            stop
         endif

         print *
            write(6,*) 'basic run parameters:',
     $           N,NNOPT,HMIN,HMAX,SEP0,
     $           TF,DTOUT,NOUT,NIT,T,
     $           NAV,ALPHA,BETA,ETA2,
     $           NGR,
     $           NRELAX,TRELAX,DT,OMEGA2
c         endif
c        stop ':)'
         DO I=1,N
            READ (12) X(I),Y(I),Z(I),AM(I),HP(I),RHO(I),VX(I),
     $           VY(I),VZ(I),VXDOT(I),VYDOT(I),VZDOT(I),A(I),
     $           ADOT(I),GRPOT(I),wmeanmolecular(i),
     $           CC(I),DIVV(I)
         ENDDO
c         READ (12) NCHK
         NCHK = N
         print *,'NCHK=',NCHK
         IF (NCHK.NE.N) STOP 'READIT: ERROR READING NCHK???'
         CLOSE (12)
      else
         STOP 'READIT: IFORM UNKNOWN ???'
      ENDIF

c Synchronize all variables:
      do i=1,n
         vx(i)=vx(i)-0.5d0*dt*vxdot(i)
         vy(i)=vy(i)-0.5d0*dt*vydot(i)
         vz(i)=vz(i)-0.5d0*dt*vzdot(i)
         if(a(i).ne.0.d0) then
            a(i)=a(i)-0.5d0*dt*adot(i)
         endif
      enddo
      t=t-0.5d0*dt

      if(omega2.ne.0.d0 .and. nrelax.eq.2) then
         theta=omega2**0.5d0*t
         print *,'omega2=',omega2,theta

         print *,
     $        'CALCULATION WAS IN ROTATING FRAME: UNDOING ROTATION!'

         do I=1,N
                                ! rotation is counterclockwise:
            xold=x(i)
            yold=y(i)
            X(i)=xold*cos(theta)-yold*sin(theta)
            y(i)=xold*sin(theta)+yold*cos(theta)
            vxold=vx(i)
            vyold=vy(i)
            vX(i)=vxold*cos(theta)-vyold*sin(theta)-omega2**0.5d0*y(i)
            vy(i)=vxold*sin(theta)+vyold*cos(theta)+omega2**0.5d0*x(i)
         enddo
      endif


      RETURN
      END
