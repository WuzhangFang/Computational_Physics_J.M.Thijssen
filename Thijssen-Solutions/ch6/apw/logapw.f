      PROGRAM apw

C APW band structure calculation for copper, as described in 
C Sec. 6.5.2 of the book "Computational Physics" by Jos Thijssen,
C Cambridge University Press, 1998.
C Program written by Jos Thijssen, Spring 1999
C
C This program calculates the band structure of a solid in the muffin-tin
C approximation using the APW method. 
C The program is more advanced than the 'standard' program apw.f.
C The difference is that a logarithmic grid is used for integrating the 
C Schrodinger equation. The potential is read from a file 'potential'.
C Comments for some subroutines can be found in the program "apw.f"

      include  "apwglob"

      CALL ReadPotential
      CALL Initialise

      CALL MakeLine (GammaPoint, KPoint)
c      CALL MakeLine (GammaPoint, XPoint)
c      CALL MakeLine (XPoint, WPoint)
c      CALL MakeLine (XPoint, UPoint)
c      CALL MakeLine (KPoint, GammaPoint)
 
      END



      SUBROUTINE Initialise

      include "apwglob" 

      CALL InitValues
      CALL FillDist
      CALL ReadKVecs
      RMuffinTin = 2.41191D0

      END



      SUBROUTINE AtomInt(Energy)

C The values of the radial solution and its derivative
C are found on the edge of the muffin-tin. 
C The radial integration is done in the routine "atom".

      include "apwglob"

      DOUBLE PRECISION Energy, PhiDer, Phi

      INTEGER L

      DO L=0, MaxL
         CALL Atom (Energy, L, Phi, PhiDer)
         ChiR(L) = Phi
         ChiDotR(L) = PhiDer
      ENDDO
      END


      SUBROUTINE ReadPotential

C The potential is read from the file 'potential' and written into
C the array 'potential'. The array 'IntPoints' contains the r-values
C of the logarithmic grid. The arrays 'MultArray' and 'DerivPoints'
C are filled for efficiency and for convenience. 

      include "apwglob"
      include "Pot.glob"

      INTEGER I

      DOUBLE PRECISION Expon

      OPEN (UNIT=4, File='potential')

      READ (4, *) RMuffinTin, Delta, PotNum
      R0 = RMuffinTin/(EXP(Delta*(Potnum-1))-1.D0)

      DO I=1, PotNum/5
         READ (4,*) Potential(I*5-4), Potential(I*5-3),
     .                      Potential(I*5-2), Potential(I*5-1),
     .                      Potential(I*5)
      ENDDO

      DO I=1, PotNum
         Expon = EXP(Delta*(I-1))
         IntPoints(I) = R0*(Expon-1.D0)
         DerivPoints(I) = (Delta*Expon*R0)**2
         MultArray(I) = SQRT(ExPon)
      ENDDO
      END


      SUBROUTINE FillFArr(First, Last, L, Ener)
C Fill the array FArr = F(R, L, E) 
C for use in the Numerov routine

      include "apwglob"
      include "Pot.glob"

      INTEGER I, L, First, Last
      DOUBLE PRECISION R, Ener, F, DeltaSq

      DeltaSq = Delta*Delta
      DO I=First, Last
        R = IntPoints(I)
        FArr(I) = DerivPoints(I)*(L*(L+1)/(R*R)-2*(Potential(I)/R+
     .            Ener))+0.25D0*DeltaSq
      END DO
      END


      SUBROUTINE Atom (Ener, L, PhiMax, DerPhiMax)
C First, the array FArr is filled. This array is then used in 
C the routine "Numerov" to find the radial solution for energy
C E and angular momentum number L.
C The value and the derivative of the radial solution
C are returned.

      include "apwglob"
      include "Pot.glob"

      DOUBLE PRECISION PhiNext, PhiStart, Ener, PhiMax, 
     .                 DerPhiMax, Solution(MaxSol), Step

      INTEGER L

      CALL FillFArr(2, PotNum, L, Ener)
      IF (L.EQ.0) THEN
        PhiStart = 2*HStep*HStep/12.D0*29.D0
      ELSE
        PhiStart = 0.D0
      ENDIF
      Step = IntPoints(2)
      PhiStart = Step**(L+1)/MultArray(2)
      Step = IntPoints(3)
      PhiNext = Step**(L+1)/MultArray(3)
      CALL Numerov(1.D0, 2, PotNum, MaxSol, FArr, 
     .             .FALSE., PhiStart, PhiNext, Solution)
      Solution(PotNum) = Solution(PotNum)*MultArray(PotNum)
      Solution(PotNum-1) = Solution(PotNum-1)*MultArray(PotNum-1)
      PhiMax = Solution(PotNum)
      PhiMax = PhiMax/RMuffinTin
      Step = IntPoints(PotNum)-IntPoints(PotNum-1)
      DerPhiMax = (Solution(PotNum)-Solution(PotNum-1))/Step
     .         +0.125D0*Step*(3*FArr(PotNum)*Solution(PotNum)+
     .            FArr(PotNum-1)*Solution(PotNum-1))
      DerPhiMax = (DerPhiMax-PhiMax)/RMuffinTin
      END 





      SUBROUTINE InitValues

C Points in the BZ are defined. The lattice constant, 
C 2*PI/LatConst, and the volume of the unit are calculated.

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

C The K-vectors of the reciprocal lattice that label
C the basis vectors in this calculation are read from 
C a file 'KVectors'. 


      include "apwglob"

      INTEGER I, J

      OPEN (UNIT=8, FILE='KVectors')
      DO I=1, NumberKA
        READ (8,*) (KVec(I, J), J=1, 3)
      ENDDO
      END



      SUBROUTINE FillDist

C The 'DistMatrix' is filled with the squares of the norms of the
C relevant K-vectors.                                            

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

C      include "apwglob"

C A line in the IBZ is scanned. On a series of points on this
C line, the Hamiltonian is diagonalised.                     

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
C The matrices A, B and C, which do not depend on the 
C energy, are calculated at the beginning of the program. 
C In the routine FillH, these matrices are used to fill 
C the actual Hamilton matrix. 

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
C The matrix H-EI is evaluated for Energy E; I is the unit matrix. 
C The determinant of this matrix is then calculated. 
 
      include "apwglob"

      DOUBLE PRECISION Energy, HMatrix(NumberKA, NumberKA), D,
     .       Deter

      INTEGER Indx(NumberKA), J, INFO, IPIVOT(NumberKA)

      CALL FillH(HMatrix, Energy)
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
C In this routine, the actual spectrum is calculated for 
C a particular vector k in the IBZ. This is done in a rather
C brute-force like way: for a fine energy grid, the determinant 
C is evaluated, and when there is a change of sign, the 
C energy where the determinant vanishes is approximated by
C a quadratic interpolation. If there is a minimum at which
C the value of the determinant is sufficiently small, the 
C minimum (actually guessed by quadratic interpolation) is
C also recorded as a spectral energy. 
C There might by a change of sign of the determinant on both sides
C of a singularity. This is recorded as a spectral value, 
C but it should be eliminated.
      
C      include "apwglob"

      INTEGER I, J, StepNum

      DOUBLE PRECISION kBZ(3), Energy, Deter, 
     .       PrevDet, PPrevDet, Step, 
     .       BestEner

      print '("start spectrum. KbZ = ", 3F10.5)', kBZ
      CALL BuildABC (kBZ)
      PPrevDet = 0.0D0
      PrevDet = 0.0D0
      I = 0
      StepNum = 100
      Step = 0.38D0/(StepNum+1)
      DO J=1, StepNum
        Energy=-0.04D0+DBLE(J-1)/StepNum*0.38D0
        CALL FindDeter(Energy, Deter)
        IF (Deter*PrevDet.LT.0.D0) THEN
          I = I+1
          WRITE (7,'(I8 F16.10)') I, 
     .           Energy+Deter*2*Step/(PrevDet-Deter)
          WRITE (6,'(I8 F16.10)') I, 
     .           Energy+Deter*2*Step/(PrevDet-Deter)
          if (energy.GT.0.27D0) then
            print *, deter, prevdet, pprevdet
          endif
        ELSE
     .    IF (((PrevDet-Deter)*(PrevDet-PPrevDet).GT.0.D0).AND.
     .      (ABS(PrevDet).LT.ABS(Deter)).AND.
     .      (ABS(PrevDet).LT.ABS(PPrevDet)).AND.
     .      (Energy.GT.-0.04D0+Step)) THEN
          BestEner = Energy+0.5D0*Step*(-PPrevDet-3*PrevDet+4*Deter)/
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

      DOUBLE PRECISION Vec1(3), Vec2(3)

      InnerProd = Vec1(1)*Vec2(1)+Vec1(2)*Vec2(2)+Vec1(3)*Vec2(3)
      
      END



      SUBROUTINE FillH (HMatrix, Energy)

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


