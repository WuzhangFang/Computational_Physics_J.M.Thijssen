      PROGRAM Scatter

C Program to evaluate the total cross section for H-Kr scattering
C as a function of energy.
C Program described in "Computational Physics", J. M. Thijssen
C Section 2.***
C Program written by J. M. Thijssen, Feb. 1998

      CALL InputParams
      Call CalcScatter

      END      



      SUBROUTINE InputParams

C Asks the user for input parameters

      include "globscatter"
 
      PI = 4.D0*DATAN(1.D0)
      Rm = 3.57D0
      Epsilon = 5.9D0
      RydConst = Rm*Rm*0.48D0

      WRITE (6,*) 'Give max. energy (e.g. 3.5)'
      READ  (5,*) MaxEner
      WRITE (6,*) 'Give min. energy (e.g. 0.1)'
      READ  (5,*) MinEner
      WRITE (6,*) 'Give number of equidistant different energy values'
      WRITE (6,*) '(e.g. 100)'
      READ  (5,*) EnerNum      
      WRITE (6,*) 'Give maximal L-value (6 is recommended)'
      READ  (5,*) LMax      
      WRITE (6,*) 'Give starting radius for integration of radial'
      WRITE (6,*) 'Schrodinger equation (e.g. 0.75)'
      READ  (5,*) Start      
      WRITE (6,*) 'Give Integration step (e.g. 0.02)'
      READ  (5,*) HStep
      WRITE (6,*) 'Give maximal radius for first integration R1',
     .           ' (e.g. 5.0)'
      READ  (5,*) MaxDist    
      MaxI = NINT((MaxDist-Start)/HStep)
      MaxDist = Start+MaxI*HStep
      print *, 'MaxDist is now ', MaxDist

      END





      SUBROUTINE CalcScatter

C Calculates cross sections for a sequence of energies and writes 
C results to a file "sigmadat"

      include "globscatter"

      DOUBLE PRECISION Delta, DeltaE, Ener, SigmaTot

      INTEGER L, I

      OPEN (4, file='sigmadat')

      DeltaE = (MaxEner-MinEner)/EnerNum
      Ener = MinEner

      DO I =1, EnerNum+1
       K = SQRT(RydConst*Ener)
       SigmaTot = 0.D0
       DO L=0, LMax
         CALL CalcDelta (Delta, L, Ener)
         SigmaTot = SigmaTot+4.D0*PI/K/K*(2*L+1)*DSIN(Delta)**2.D0
       ENDDO
       WRITE(6,1000) I, Ener, SigmaTot, TAN(Delta)
       WRITE(4, '(2F12.4)') Ener, SigmaTot
       Ener = Ener + DeltaE
      ENDDO

1000  FORMAT (I6, ' Ener = ', F12.4,'   SigmaTot=',F12.4, 
     .       '  TanDelta = ', F12.4)
      CLOSE(4)
      WRITE (6,*) 'Output written to file ''sigmadat'''
      END



      SUBROUTINE CalcDelta (Delta, L, Ener)

C Calculate the phase shift 
C Delta: Phase shift
C L: angular momentum
C Ener: Energy

      IMPLICIT NONE

      DOUBLE PRECISION GQuotient, Delta, PhaseSft, Ener, SecR

      INTEGER L

      CALL FindGQuotient(L, GQuotient, Ener, SecR)
      
      Delta = PhaseSft (L, GQuotient, SecR)

      END      



      SUBROUTINE FindGQuotient(L, GQuotient, Ener, SecR)

C Find GQuotient defined in Eq. ***, from which the phase shift can be found
C L: Angular momentum
C Ener: Energy
C MaxR: Maximum integration radius for first integration
C SecR: MaxR+1/4 of a wavelength

      include "globscatter"

      DOUBLE PRECISION GQuotient, SecR, Ener, ExpCoeff, F,
     .                Phi1, Phi2, PhiEnd1, PhiEnd2,
     .                QuartWav, Solution(MaxSol), a, b, c,
     .                PhiDeriv

      INTEGER L, SecMaxI

      EXTERNAL F

      ExpCoeff=SQRT(Epsilon*RydConst/25.D0)
      
      Phi1 = EXP(-ExpCoeff*Start**(-5))
      PhiDeriv = 5*ExpCoeff*Start**(-6)*Phi1
      a = HStep*HStep*F(Start+HStep, L, Ener)/12.D0
      b = HStep*HStep*F(Start-HStep, L, Ener)/12.D0
      c = HStep*HStep*F(Start, L, Ener)
      c = (2.D0+5.D0/6.D0*c)*Phi1
      Phi2 = 2*HStep*PhiDeriv*(1-b)+c*(1-2*b)
      Phi2 = Phi2/((1-2*a)*(1-b)+(1-a)*(1-2*b))

      QuartWav = 0.5D0*PI/K
      SecR = MaxDist + QuartWav
      SecMaxI = NINT((SecR-Start)/HStep)
      
      SecR = Start+SecMaxI*HStep
      CALL FillFArr(1, SecMaxI, L, Ener)

      CALL Numerov(HStep, 1, SecMaxI+1, MaxSol, FArr, 
     .            .FALSE., Phi1, Phi2, Solution)
      PhiEnd1 = Solution(MaxI+1)
      PhiEnd2 = Solution(SecMaxI+1)

      GQuotient = PhiEnd2*MaxDist/(PhiEnd1*SecR)

      END



      SUBROUTINE FillFArr(First, Last, L, Ener)
C Fill the array FArr = F(R, L, E) (see below)
C for use in the Numerov routine

      include "globscatter"

      INTEGER I, L, First, Last
      DOUBLE PRECISION R, Ener, F

      EXTERNAL F

      DO I=First, Last
        R = Start+(I-1)*HStep
        FArr(I) = F(R, L, Ener)
      END DO
      END



      DOUBLE PRECISION FUNCTION F (R, L, Ener)

C Function F defined in Eq. ***
C R: radial distance
C L: Angular momentum
C Ener: Energy

      include "globscatter"

      DOUBLE PRECISION R, R2, Ener, R6, R12, R4
       
      INTEGER L

      R2 = R*R
      R4 = R2*R2
      R6 = R4*R2
      R12 = R6*R6

      F = L*(L+1)/R2 + RydConst*Epsilon*(1.D0/R12 - 2.D0/R6)
     .     - RydConst*Ener

      END





      DOUBLE PRECISION FUNCTION PhaseSft (L, GQuotient, SecR)

C Calculates the phase shift Delta_l, using Eq. ***
C GQuotient: defined in Eq. ***
C L: Angular momentum
C MaxR: First radius, see Eq. ***
C SecR: Second radius, see Eq. ***


      include "globscatter"

      DOUBLE PRECISION GQuotient, SecR, Help, SphBesJ, SphBesN

      INTEGER L

      Help = (GQuotient*SphBesJ(L,K*MaxDist)-SphBesJ(L,K*SecR))
      Help = Help/(GQuotient*SphBesN(L,K*MaxDist)-SphBesN(L, K*SecR))
      PhaseSft = ATAN(Help)

      END 
