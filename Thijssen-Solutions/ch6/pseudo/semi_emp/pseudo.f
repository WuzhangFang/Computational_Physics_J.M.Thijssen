      PROGRAM Pseudo

C *********************************************************************
C * This program calculates the electronic structure of silicon.      *
C * The pseudo-potential method is used. Input for the calculation    *
C * are the pseudo-potential matrix elements.                         *
C * The states taken into account in the calculation are characterised*
C * by their K-values. Only states with |K|^2<= 20 are treated.       *
C * Part of these are treated in a direct way, the remaining ones     *
C * in a Lowdin  perturbation expansion.                              *
C *                                                                   *
C * The global variables are stored in the file "pseudo.glb".         *
C * They include :                                                    * 
C *     the K-values, stored in          Kvec;                        *
C *     the coordinates of the                                        *
C *     symmetrical points in the BZ:    XPoint etc.;                 *
C *     the matrix elements of                                        *
C *     the pseudo-potential             VSqrtN;                      *
C *     the reciprocal lattice constant: RecLatConst;                 *
C *     the matrix containing the                                     *
C *     norm of the various V-vectors    DistMatrix                   *
C *     the total number of K-vectors:   NumberKAB;                   *
C *     the number of K-vectors treated                               * 
C *                      exactly:        NumberKA.                    *
C *                                                                   *
C * Results are written into the file BandStr in the following form:  *
C * A line contains three real numbers; these are the x, y, and       *
C * z-components of the k-vector in the BZ.                           *
C * The next line then contains the energy spectrum for that k-vector.*
C * This program is described in section 6.7.1 of the book            *
C * "Computational Physics" by Jos Thijssen                           *
C * Cambridge University Press, 1999                                  *
C * Program written by J. Thijssen.                                   *
C * Final form: April 1998                                            * 
C *********************************************************************

      include  "pseudo.glb"

      CALL Initialise

      CALL InitSpectrum (LPoint)

      CALL MakeLine (LPoint, GammaPoint, 1)
c      CALL MakeLine (GammaPoint, XPoint, 2)
c      CALL MakeLine (XPoint, UPoint, 3)
c      CALL MakeLine (KPoint, GammaPoint, 4)

      END



      SUBROUTINE Initialise
C Initialise some parameters and read in K vectors (reciprocal lattice)

      include "pseudo.glb" 

      CALL InitValues
      CALL ReadKVecs
      CALL FillDist
      OPEN (UNIT=7, FILE='BandStr')
      WRITE (7,*) '# A line contains three real numbers; these are the'
      WRITE (7,*) '# x, y and C * z-components of the k-vector'
      WRITE (7,*) '# in the BZ.'
      WRITE (7,*) '# The next line then contains the energy spectrum'
      WRITE (7,*) '# for that k-vector. '
      END



      SUBROUTINE InitValues

C ***  Coordinates of symmetrical points in the BZ are initialised ***
C ***  Potential-parameters are initialised.                       ***


      include "pseudo.glb"

      DOUBLE PRECISION Point1(3), Point2(3), Point3(3), 
     .                 Point4(3), Point5(3)
      INTEGER I

      DATA  Point1    /0.0D0 , 0.0D0 , 0.0D0/,
     .      Point2    /0.5D0 , 0.5D0 , 0.5D0/,
     .      Point3    /0.75D0, 0.75D0, 0.0D0 /,
     .      Point4    /1.D0  , 0.D0  , 0.D0 /,
     .      Point5    /1.D0  , 0.25D0  , 0.25D0/

      DO I=1,3
         GammaPoint(I) = Point1(I)
         LPoint(I)     = Point2(I)
         KPoint(I)     = Point3(I)
         XPoint(I)     = Point4(I)
         UPoint(I)     = Point5(I)
      END DO

      PseudConst(1) = -0.2241D0
      PseudConst(2) =  0.0551D0
      PseudConst(3) =  0.0724D0
      PseudConst(0) =  0.D0
      END



      SUBROUTINE ReadKVecs

C *** The K-vectors of the reciprocal lattice which are involved  ***
C *** in this calculation are read from a file 'KVectors'.        ***


      include "pseudo.glb"

      INTEGER I, J

      OPEN (UNIT=8, FILE='KVectors')
      READ (8,*) ((KVec(I, J), J=1, 3), I=1,NumberKAB)
      END



      SUBROUTINE FillDist

C *** The 'DistMatrix' is filled with the squares of the norms of the ***
C *** relevant K-vectors.                                             ***

      include "pseudo.glb"

      INTEGER I, J, K

      DO I=-3,3
        DO J=-3,3
          DO  K=-3,3
            DistMatrix(I,J,K) = 3*(I*I+J*J+K*K)-
     .                          2*I*J-2*J*K-2*I*K
          END DO
        END DO
      END DO
      END



      SUBROUTINE MakeLine (FirstPoint, SecPoint, SegNum)

C *** A line in the IBZ is scanned. On a series of points on this ***
C *** line, the Hamiltonian is diagonalised.                      ***

      DOUBLE PRECISION FirstPoint(3), SecPoint(3), Delta(3), Length, 
     .                 Step(3), StepLength, kBZ(3)
      INTEGER SegNum

C StepLength is the cartesian distance between subsequent points along a line 
C in the BZ. 
      PARAMETER (StepLength = 0.03D0)

      INTEGER I, L, Num

      DO L=1,3
        Delta(L) = SecPoint(L)-FirstPoint(L)
      END DO
      Length = Delta(1)*Delta(1)+Delta(2)*Delta(2)+Delta(3)*Delta(3)
      Length = SQRT(Length)
      Num = NINT (Length/StepLength)

      DO L=1,3
        Step(L) = Delta(L)/Num
      END DO

      DO L=1,3
        kBZ(L) =  FirstPoint(L)
      END DO
      DO I=1, Num
        CALL Spectrum (kBZ, SegNum, Num)
        DO L=1,3
          kBZ(L) = kBZ(L) + Step(L)
        END DO
      END DO

      END
        


      SUBROUTINE InitSpectrum (kBZ)

C *** For the first point to be calculated, the full Hamiltonian is ***
C *** diagonalised.                                                 ***


      include "pseudo.glb"


      DOUBLE PRECISION kBZ(3), HMatrix(NumberKAB, NumberKAB),
     .                 Diag(NumberKAB)

      INTEGER Level

      CALL FillH (HMatrix, kBZ)
      CALL DiagMat(HMatrix, Diag, NumberKA, NumberKAB)

      DO Level = 1, MaxLevel
        LowLevels(Level) = Diag(Level)      
      END DO
      END




      SUBROUTINE Spectrum (kBZ, SegNum, KNum)

C *** Calculation of the spectrum on subsequent points using Lowdin's ***
C *** perturbation theory. The energy used to build the Lowdin matrix ***
C *** is for every level equal to the energy value found in the       ***
C *** previous calculation.                                           ***
C kBZ is the BZ, I is

      include "pseudo.glb"

      INTEGER Level, I, SegNum, KNum

      DOUBLE PRECISION kBZ(3), HMatrix(NumberKAB, NumberKAB),
     .             UMatrix(NumberKA,  NumberKA),
     .             Diag(NumberKA)

      CALL FillH (HMatrix, kBZ)

      DO Level = 1, MaxLevel
        CALL BuildU (HMatrix, UMatrix, LowLevels(Level))
        CALL DiagMat(UMatrix, Diag, NumberKA, NumberKA)
        LowLevels(Level) = Diag(Level)
      END DO
      WRITE (7, '(3F8.4)') kBZ
      WRITE (7,'(8F8.4)') (LowLevels(I), I=1, 8)
      
      END




      SUBROUTINE SelectLow (Diag, Range, LowLevels, MaxLevel)

C *** This subroutine selects the Range lowest eigenvalues from the 
C spectrum. in Diag, and stores them in the array LowLevels, which
C is MaxLevels long.

      IMPLICIT NONE

      INTEGER I, Level, SecLevel, Range, LCount, MaxLevel

      DOUBLE PRECISION Diag(Range), LowLevels(MaxLevel), Temp

      LowLevels(1) = Diag(1)
      Level = 1
      DO I=2, Range
        Temp = Diag(I)
        IF (Temp .LT. LowLevels(Level)) THEN
          SecLevel = Level
          DO WHILE ((SecLevel.GT.1).AND.(LowLevels(SecLevel-1).GT.Temp))
            SecLevel = SecLevel-1
          END DO
          IF (Level.LT.MaxLevel) Level = Level+1
          DO LCount = Level, SecLevel+1, -1
            LowLevels(LCount) = LowLevels(LCount-1)
          END DO
          LowLevels(SecLevel) = Temp
        END IF
      END DO

      END

               


      SUBROUTINE BuildU (HMatrix, UMatrix, Energy)

C *** In this subroutine, the Lowdin matrix is calculated in lowest ***
C *** order perturbation theory.                                    ***

      include "pseudo.glb"
 

      DOUBLE PRECISION Energy, UMatrix(NumberKA, NumberKA), 
     .               HMatrix(NumberKAB, NumberKAB)
      INTEGER IA, JA, IB

      DO IA = 1,NumberKA
        DO JA = IA,NumberKA
          UMatrix(Ia, Ja) = HMatrix(Ia, Ja)
          DO IB = NumberKA+1, NumberKAB
            UMatrix(Ia, Ja) = UMatrix(Ia, Ja) + 
     .        HMatrix(Ia, Ib)*HMatrix(Ib, Ja)/(Energy-HMatrix(Ib,Ib))
          END DO
          UMatrix(Ja, Ia) = UMatrix (Ia, Ja)
        END DO
      END DO

      END



      INTEGER FUNCTION Code(I)
      IMPLICIT NONE
      INTEGER I
      Code = 0
      IF (I.EQ.3) Code = 1
      IF (I.EQ.8) Code = 2
      IF (I.EQ.11) Code = 3
      END



      SUBROUTINE FillH (HMatrix, kBZ)
C The hamiltonian HMatrix is filled for kBZ in the BZ
C see Eq. (6.51)

      include "pseudo.glb"

      DOUBLE PRECISION PI

      PARAMETER (PI=3.14159265358979D0)

      DOUBLE PRECISION HMatrix(NumberKAB, NumberKAB), kBZ(3),
     .       InProd, KTot(3), NormKTot

      INTEGER I, J, L, Dist, DeltaK(3), Code
      DO I=1, NumberKAB
        DO J=1,I-1
          DO L=1, 3
            DeltaK(L) = KVec(I, L) - KVec(J,L)
          END DO
          HMatrix(I, J) = 0.D0
          Dist = DistMatrix(DeltaK(1), DeltaK(2), DeltaK(3))
          InProd  = PI/4*DBLE(DeltaK(1)+DeltaK(2)+DeltaK(3))
          HMatrix(I,J) = PseudConst(Code(Dist))*COS(InProd)
          HMatrix(J,I) = HMatrix(I,J)
        END DO
        Ktot(1) = kBZ(1) - KVec(I,1) + KVec(I,2) + KVec(I,3)
        Ktot(2) = kBZ(2) + KVec(I,1) - KVec(I,2) + KVec(I,3)
        Ktot(3) = kBZ(3) + KVec(I,1) + KVec(I,2) - KVec(I,3)
        NormKTot = KTot(1)*KTot(1)+KTot(2)*KTot(2)+KTot(3)*KTot(3)
        HMatrix (I,I) = RecLatConst*RecLatConst*NormKTot
      END DO
      END


