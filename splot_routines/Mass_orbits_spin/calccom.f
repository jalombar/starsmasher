      subroutine calccom
      include'spha.h'
      real*8 am1,am2,am3,am4
      real*8 X1,Y1,Z1,VX1,VY1,VZ1
      real*8 X2,Y2,Z2,VX2,VY2,VZ2
      real*8 X3,Y3,Z3,VX3,VY3,VZ3
      integer icomp(nmax),i
      COMMON/COMPBETTERCOM3/AM1,X1,Y1,Z1,VX1,VY1,VZ1,
     $     AM2,X2,Y2,Z2,VX2,VY2,VZ2,
     $     AM3,X3,Y3,Z3,VX3,VY3,VZ3,am4,
     $     ICOMP

C     Determine components CM positions and velocities
      AM1=0.d0
      X1=0.d0
      Y1=0.d0
      Z1=0.d0
      VX1=0.d0
      VY1=0.d0
      VZ1=0.d0
      AM2=0.d0
      X2=0.d0
      Y2=0.d0
      Z2=0.d0
      VX2=0.d0
      VY2=0.d0
      VZ2=0.d0
      AM3=0.d0
      X3=0.d0
      Y3=0.d0
      Z3=0.d0
      VX3=0.d0
      VY3=0.d0
      VZ3=0.d0
      am4=0.d0

c      print *,'calccom: N=',N
      DO I=1,N
         IF (ICOMP(I).EQ.1) THEN
            AM1=AM1+AM(I)
            X1=X1+AM(I)*X(I)
            Y1=Y1+AM(I)*Y(I)
            Z1=Z1+AM(I)*Z(I)
            VX1=VX1+AM(I)*VX(I)
            VY1=VY1+AM(I)*VY(I)
            VZ1=VZ1+AM(I)*VZ(I)
         ELSE IF (ICOMP(I).EQ.2) THEN
            AM2=AM2+AM(I)
            X2=X2+AM(I)*X(I)
            Y2=Y2+AM(I)*Y(I)
            Z2=Z2+AM(I)*Z(I)
            VX2=VX2+AM(I)*VX(I)
            VY2=VY2+AM(I)*VY(I)
            VZ2=VZ2+AM(I)*VZ(I)
         ELSE IF (ICOMP(I).EQ.3) THEN
            AM3=AM3+AM(I)
            X3=X3+AM(I)*X(I)
            Y3=Y3+AM(I)*Y(I)
            Z3=Z3+AM(I)*Z(I)
            VX3=VX3+AM(I)*VX(I)
            VY3=VY3+AM(I)*VY(I)
            VZ3=VZ3+AM(I)*VZ(I)
         ELSE IF (ICOMP(I).EQ.4) THEN
            AM4=AM4+AM(I)
         ENDIF
      ENDDO
      if(AM1.gt.0.d0) then
         X1=X1/AM1
         Y1=Y1/AM1
         Z1=Z1/AM1
         VX1=VX1/AM1
         VY1=VY1/AM1
         VZ1=VZ1/AM1
      endif
      if(AM2.gt.0.d0) then
         X2=X2/AM2
         Y2=Y2/AM2
         Z2=Z2/AM2
         VX2=VX2/AM2
         VY2=VY2/AM2
         VZ2=VZ2/AM2
      endif
      if(AM3.gt.0.d0) then
         X3=X3/AM3
         Y3=Y3/AM3
         Z3=Z3/AM3
         VX3=VX3/AM3
         VY3=VY3/AM3
         VZ3=VZ3/AM3
      endif

      return
      end
