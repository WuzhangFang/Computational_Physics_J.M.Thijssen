C APW band structure calculation for copper, as described in 
C Sec. 6.5.2 of the book "Computational Physics" by Jos Thijssen,
C Cambridge University Press, 1998.
C Second edition: 2007
C Program written by Jos Thijssen, Spring 1999
C Two remarks concerning the potential
C The muffin-tin potential is taken from a self-consistent 
C layer-KKR calculation by Maziar Nekovee
C This potential is given in an accurate, parametrised form. 
C The integration grid is uniform. 
C The program "logapw.f" uses the potential on a logarithmic grid.
C 

      PROGRAM apw
C Main program 
      include  "apwglob"

      CALL Initialise

      CALL MakeLine (GammaPoint, KPoint)
c      CALL MakeLine (GammaPoint, XPoint)
c      CALL MakeLine (XPoint, WPoint)
c      CALL MakeLine (XPoint, UPoint)
c      CALL MakeLine (KPoint, GammaPoint) 
      CLOSE(7)
      END



      SUBROUTINE Initialise
C Initialisation of parameters; K-vectors are read in

      include "apwglob" 

      CALL InitValues
      CALL FillDist
      CALL ReadKVecs
      RMuffinTin = 2.41191D0
      OPEN(7, File='spectrum')

      END



      SUBROUTINE AtomInt(Energy)
C Integration of the radial Schrodinger equation  at energy E
C (actually performed inside the routine "atom").
C Value and derivative at MT boundary are Phi and PhiDer 
C respectively. These are stored for every L-value in the arrays 
C ChiR and ChiDotR

      include "apwglob"

      DOUBLE PRECISION Energy, PhiDer, Phi

      INTEGER L

      DO L=0, MaxL
         CALL Atom (Energy, L, Phi, PhiDer)
         ChiR(L) = Phi
         ChiDotR(L) = PhiDer
      ENDDO
      END



      SUBROUTINE FillFArr(First, Last, L, Ener)
C Fill the array FArr = F(R, L, E) (see below)
C for use in the Numerov routine

      include "apwglob"

      INTEGER I, L, First, Last
      DOUBLE PRECISION R, Ener, F, Vcu

      EXTERNAL F

      DO I=First, Last
        R = I*HStep
        FArr(I) = L*(L+1)/(R*R) - 2*(Vcu(R)/R + Ener)
      END DO
      END



      SUBROUTINE Atom (Ener, L, PhiMax, DerPhiMax)
C Atom initialises the starting values for the radial integration,
C it calls Numerov to perform this integration for angular momentum L
C and energy E, and returns the value and derivative of the solution
C in PhiMax and PhiDer Max respectively

      include "apwglob"

      DOUBLE PRECISION PhiNext, PhiStart, Ener, PhiMax, 
     .                 DerPhiMax, Solution(MaxSol)

      INTEGER L

      PotNum = 1000/(L+1)
      HStep=RMuffinTin/PotNum
      
      CALL FillFArr(1, PotNum, L, Ener)
      IF (L.EQ.0) THEN
        PhiStart = 2*HStep*HStep/12.D0*29.D0
C Analytic expression; 29 is the nuclear charge of copper
      ELSE
        PhiStart = 0.D0
      ENDIF
      PhiNext = HStep**(L+1)
      CALL Numerov(HStep, 0, PotNum, MaxSol, FArr, 
     .             .TRUE., PhiStart, PhiNext, Solution)
      PhiMax = Solution(PotNum) 
      PhiMax = PhiMax/RMuffinTin
C A somewhat complicated expression is used to calculate the derivative 
C at the last point, with accuracy of order h^3 
      DerPhiMax = (Solution(PotNum)-Solution(PotNum-1))/HStep
     .         +0.125D0*HStep*(3*FArr(PotNum)*Solution(PotNum)+
     .            FArr(PotNum-1)*Solution(PotNum-1))
      DerPhiMax = (DerPhiMax-PhiMax)/RMuffinTin
      END 







      SUBROUTINE InitValues

C ***  Coordinates of symmetrical points in the BZ are initialised ***
C ***  Potential-parameters are initialised.                       ***


      include "apwglob"

      DOUBLE PRECISION Point1(3), Point2(3), Point3(3), Point4(3), 
     .                 Point5(3)
      INTEGER I

      DATA  Point1    /0.0D0 , 0.0D0 , 0.0D0/,
     .      Point2    /0.5D0 , 0.5D0 , 0.5D0/,
     .      Point3    /0.75D0, 0.75D0, 0.0D0 /,
     .      Point4    /1.D0  , 0.D0  , 0.D0 /,
     .      Point5    /1.D0  , 0.5D0  , 0.0D0/

      DO I=1,3
        GammaPoint(I) = Point1(I)
        LPoint(I)     = Point2(I)
        KPoint(I)     = Point3(I)
        XPoint(I)     = Point4(I)
        WPoint(I)     = Point5(I)
      ENDDO

      PI = 4.D0*DATAN(1.D0)

      LatConst = 6.8219117D0
      RecLatConst = 2*PI/LatConst
      UnitVol = (LatConst**3)/4.D0

      END



      SUBROUTINE ReadKVecs

C *** The K-vectors of the reciprocal lattice which label
C *** the basis vectors in this calculation are read from 
C *** a file 'KVectors'. 


      include "apwglob"

      INTEGER I, J

      OPEN (UNIT=8, FILE='KVectors')
      DO I=1, NumberKA
        READ (8,*) (KVec(I, J), J=1, 3)
      ENDDO
      END



      SUBROUTINE FillDist

C *** The 'DistMatrix' is filled with the squares of the norms of the 
C *** relevant K-vectors.


      include "apwglob"

      INTEGER I, J, K

      DO I=-6,6
        DO J=-6,6
          DO K=-6,6
            DistMatrix(I,J,K) = 3*(I*I+J*J+K*K)-
     .                          2*I*J-2*J*K-2*I*K
          ENDDO
        ENDDO
      ENDDO
      END



      SUBROUTINE MakeLine (FirstPoint, SecPoint)

      include "apwglob"

C *** A line in the IBZ is scanned. On a series of points on this 
C *** line, the Hamiltonian is diagonalised.

      DOUBLE PRECISION FirstPoint(3), SecPoint(3), StepArr(3), Length, 
     .                 Step(3), StepLength, kBZ(3)
      INTEGER J, L, Num

      PARAMETER (StepLength = 0.02D0)

      DO L=1,3
        StepArr(L) = SecPoint(L)-FirstPoint(L)
      ENDDO
      Length = StepArr(1)*StepArr(1)+StepArr(2)*StepArr(2)+
     .         StepArr(3)*StepArr(3)
      Length = SQRT(Length)
      Num = NINT (Length/StepLength)
      WRITE (6,*) 'nr of steps', Num

      DO  L=1,3
        Step(L) = StepArr(L)/Num
      ENDDO

      DO L=1,3
        kBZ(L) =  FirstPoint(L)
      ENDDO
      DO J=0, Num
        DO L=1,3
          kBZ(L) = FirstPoint(L) + J*Step(L)
        ENDDO
        CALL Spectrum (kBZ)
      ENDDO
      END
        




      SUBROUTINE BuildABC(kBZ)
C The A, B and C-matrices are filled according to 
C Eq. (6.28)

      include "apwglob"

      INTEGER I, J, L, DeltaK(3)
      DOUBLE PRECISION Arg1, Arg2, Arg3, PreFac, kBZ(3),
     .       KTot(3), QTot(3), Dist, J1,
     .       JL1, JL2, InnerProd, KdotQ,
     .       KNorm, QNorm, X, SphBesJ, Legendre
      print *, 'enter kBZ'
      PreFac = 2*PI*RMuffinTin*RMuffinTin/UnitVol
      DO I=1, NumberKA
        DO J=1, I
          DO L=1, 3
            DeltaK(L) = KVec(I, L) - KVec(J,L)
          ENDDO
          MatA(I, J) = 0.D0
          MatB(I, J) = 0.D0
          Dist = SQRT(DBLE(DistMatrix(DeltaK(1),DeltaK(2),DeltaK(3))))*
     .           RecLatConst
          Ktot(1) = kBZ(1) - KVec(I,1) + KVec(I,2) + KVec(I,3)
          Ktot(2) = kBZ(2) + KVec(I,1) - KVec(I,2) + KVec(I,3)
          Ktot(3) = kBZ(3) + KVec(I,1) + KVec(I,2) - KVec(I,3)
          Qtot(1) = kBZ(1) - KVec(J,1) + KVec(J,2) + KVec(J,3)
          Qtot(2) = kBZ(2) + KVec(J,1) - KVec(J,2) + KVec(J,3)
          Qtot(3) = kBZ(3) + KVec(J,1) + KVec(J,2) - KVec(J,3)
          KdotQ = InnerProd(KTot, QTot)*RecLatConst*RecLatConst
          KNorm = SQRT(InnerProd(KTot, KTot))*RecLatConst
          QNorm = SQRT(InnerProd(QTot, QTot))*RecLatConst
          X = Dist*RMuffinTin
          IF (X.GT.1E-8) THEN
            J1 = (sin(x)/x-cos(x))/x
            J1 = J1/Dist
          ELSE
            J1 = RMuffinTin/3.d0
          ENDIF
          MatA(I, J) = -2.D0*PreFac*J1
          MatB(I, J) = MatA(I,J)*0.5D0*KdotQ
          Arg1 = QNorm*RMuffinTin
          Arg2 = KNorm*RMuffinTin
          IF (KNorm.GT.1E-8.AND.QNorm.GT.1E-8) THEN
            Arg3 = KdotQ/KNorm/QNorm
          ELSE
            Arg3 = 1.D0-1.D-10
          ENDIF
          DO L=0, MaxL
            IF (Arg1.GT.1E-8) THEN
              JL1=SphBesJ (L, Arg1)
            ELSEIF (L.EQ.0) THEN
              JL1 = 1.D0
            ELSE
              JL1 = 0.D0
            ENDIF
            IF (Arg2.GT.1E-8) THEN
              JL2=SphBesJ (L, Arg2)
            ELSEIF (L.EQ.0) THEN
              JL2 = 1.D0
            ELSE
              JL2 = 0.D0
            ENDIF
            MatC(I,J, L) = (2*L+1)*PreFac*Legendre(L,Arg3)*JL1*JL2
          ENDDO
          IF (I.EQ.J) THEN
            MatA(I,I) = MatA(I,I) + 1.D0
            MatB(I,I) = MatB(I,I) + 0.5D0*KdotQ
          ENDIF
        ENDDO
      ENDDO
      END



      SUBROUTINE FindDeter (Energy, Deter)
C Find the determinant. An LU decomposition is preformed on the
C Hamiltonian at energy "Energy", and the diagonal elements are 
C multiplied to obtain the determinant Deter. 

      include "apwglob"

      DOUBLE PRECISION Energy, HMatrix(NumberKA, NumberKA), D,
     .       Deter, Deter1, TMat(NumberKA, NumberKA), D1, 
     .       Work(NumberKA), PlayMate(NumberKA, NumberKA)

      INTEGER J, IPIVOT(NumberKA), INFO, I1, I2, Sign

      CALL FillH(HMatrix, Energy)
! The determinant did not come out right using DSYTRF!!
      CALL DGETRF(NumberKA, NumberKA, HMatrix, NumberKA, IPIVOT, 
     .                INFO)
      D = 1.D0
      DO J=1, NumberKA
        IF (IPivot(J).EQ.J) THEN
          D = D*HMatrix(J,J)
        ELSE
          D = -D*HMatrix(J,J)
        END IF
      ENDDO
      Deter = D
      END



      SUBROUTINE Spectrum (kBZ)
C Calculates the energy spectrum at the point kBZ in the
C Brillouin zone.
      
      include "apwglob"

      INTEGER I, J, StepNum

      DOUBLE PRECISION kBZ(3), Energy, Deter, 
     .       PrevDet, PPrevDet, EStep, 
     .       BestEner

!      PARAMETER (Step = 1.D-3)

      print '("start spectrum. KbZ = ", 3F10.5)', kBZ
      CALL BuildABC (kBZ)
      PPrevDet = 0.0D0
      PrevDet = 0.0D0
      I = 0
      StepNum = 100
      EStep = 0.38D0/StepNum
      DO J=1, StepNum
        Energy=-0.04D0+DBLE(J-1)/StepNum*0.38D0
        CALL FindDeter(Energy, Deter)
        IF (Deter*PrevDet.LT.0.D0) THEN
          I = I+1
C Linear interpolation in the energy to fine zero of the
C determinant
          WRITE (7,'(I8 F16.10)') I, 
     .           Energy+Deter*2*EStep/(PrevDet-Deter)
          WRITE (6,'(I8 F16.10)') I, 
     .           Energy+Deter*2*EStep/(PrevDet-Deter)
        ELSE
C Quadratic interpolation in energy to find the minimum of
C the determinant
     .    IF (((PrevDet-Deter)*(PrevDet-PPrevDet).GT.0.D0).AND.
     .      (ABS(PrevDet).LT.ABS(Deter)).AND.
     .      (ABS(PrevDet).LT.ABS(PPrevDet)).AND.
     .      (Energy.GT.-0.04D0+EStep)) THEN
          BestEner = Energy+0.5D0*EStep*(-PPrevDet-3*PrevDet+4*Deter)/
     .                                   (PPrevDet-2*PrevDet+Deter)
          I = I+1
          WRITE (7,'(I8 F16.10)') I, BestEner
          WRITE (6,'(I8 F16.10)') I, BestEner
        END IF
        PPrevDet = PrevDet
        PrevDet = Deter
      ENDDO

      END

      


               

      DOUBLE PRECISION FUNCTION InnerProd(Vec1, Vec2)
C Calculates Inner Product between the two three-dim vectors 
C Vec1 and Vec2

      DOUBLE PRECISION Vec1(3), Vec2(3)

      InnerProd = Vec1(1)*Vec2(1)+Vec1(2)*Vec2(2)+Vec1(3)*Vec2(3)
      
      END



      SUBROUTINE FillH (HMatrix, Energy)
C The Hamiltonian "HMatrix" is calculated at energy E,
C using Eq. (6.27)

      include "apwglob"

      DOUBLE PRECISION HMatrix(NumberKA, NumberKA), Energy

      INTEGER I, J, L

      CALL AtomInt (Energy)
      DO I=1, NumberKA
        DO J=1,I
          HMatrix(I, J) = 0.D0
          HMatrix(I, J) = -Energy*MatA(I, J) + MatB(I, J)
          DO L=0, MaxL
            HMatrix(I,J) = HMatrix(I,J) + MatC(I, J, L)*
     .                     ChiDotR(L)/ChiR(L)
          ENDDO
          HMatrix(J,I) = HMatrix(I,J)
        ENDDO
      ENDDO
      END



      DOUBLE PRECISION FUNCTION Vcu(R)
C Parametrisation of the potential

      IMPLICIT NONE

      DOUBLE PRECISION R

      Vcu = 29*EXP(-2.3151241717834D0*R**0.81266614122432D0+
     .              2.1984250222603D-02*R**4.2246376280056D0)
     .       -0.15595606773483D0*R-3.1350051440417D-03*R**2+
     .       5.1895222293006D-02*R**3 - 2.8027608685637D-02*R**4
      END
