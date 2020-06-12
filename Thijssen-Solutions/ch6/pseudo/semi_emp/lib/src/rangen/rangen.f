C******************************************************************
C* This file contains some random generator functions based on    *
C* the shift register principle.                                  *
C*                                                                *
C* This program is written in f77.                                *
C*                                                                *
C* The routine "InitRand(InitJ)" must be called before anything   *
C* else. This initialises the random generator with seed InitJ.   *
C* The function "RealRand()" yields a random double precision,    *
C* uniformly distributed between 0 and 1.                         *
C* The function "InitRand(Num)" yields a random integer with      *
C* uniform distribution between 0 and Num-1                       *
C******************************************************************


      SUBROUTINE Ran(JRand)
      include "rand.glob"
      INTEGER JRand, MaxInt, Factor
      PARAMETER (MaxInt=2147483647, Factor=16807)

      JRand=JRand*factor
      if (JRand.LT.0) then
        JRand = JRand+MaxInt+1
      end if
      end




      SUBROUTINE InitRand(InitJ)
      include "rand.glob"
      INTEGER InitJ, I, JRand

      JRand = InitJ
      DO I=1, 256
        CALL Ran(JRand)
        RanArr(I) = JRand
      END DO
      PI = 4.D0*ATAN(1.D0)
      PCount = 0
      END



      INTEGER FUNCTION RandGen()
      include "rand.glob"
      INTEGER Pos1, Pos2, SCount

      Pos1 = (PCount)
      Pos2 = MOD (PCount+147,251)
      RanArr(PCount+1) = IEOR(RanArr(Pos1+1),RanArr(Pos2+1))
      SCount = PCount
      PCount = MOD(PCount+1, 251)
      RandGen = RanArr(scount+1)
      END




      INTEGER FUNCTION IntRand(Num)
      include "rand.glob"
      INTEGER Num, Randgen

      IntRand = RealFac*Num*Randgen()
      END


      DOUBLE PRECISION FUNCTION RealRand()
      include "rand.glob"

      INTEGER RandGen

      RealRand = RandGen()*RealFac
      END 


      SUBROUTINE ExpRand(R1, R2)
      include "rand.glob"
      
      DOUBLE PRECISION R1, R2, Norm, RealRand, Phi

      Norm = RealRand()
      Phi = RealRand()*2.D0*PI
      Norm = SQRT(-2.D0*LOG(Norm))
      R1 = Norm*cos(Phi)
      R2 = Norm*sin(Phi)
      END






