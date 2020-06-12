       PROGRAM DFT
C Program to calculate the l=0 ground state of the helium
C atom in full DFT, using the Slater X-alpha exchange potential
C and the parametrisation of the Ceperley-Alder correlation potential
C (Parametrisation given by Perdew and Zunger). 
C Program described in "Computational Physics", J. M. Thijssen
C Problem 5.3
C Program written by J. M. Thijssen, July 1998

       INCLUDE 'globDFT'

       DOUBLE PRECISION Energy, Precision, OldEnergy
       PARAMETER (Precision = 1.D-7)

       CALL Initialise

       OldEnergy = -1.D3
       Energy = -0.3D0

C Self-consistency loop:
       DO WHILE(ABS(OldEnergy-Energy).GT.Precision)
         OldEnergy = Energy
         CALL FindBound(Energy, Precision)
         CALL SolveRad(Energy)
         CALL HartCorrect
         CALL ExcCorrect
         CALL CorrCorrect
         print *, 'Eigenvalue', Energy
         print *, 'total energy',2*Energy-HartCorr+ExcCorr+CorrCorr
         CALL CalcHartPot
         CALL CalcExcPot
         CALL CalcCorrPot
       END DO

       END


       SUBROUTINE FindBound(Energy, Precision)
C Find the energy of the bound state.
C The bound state is characterised as having a zero at r=0

       IMPLICIT NONE

       DOUBLE PRECISION Energy, Low, Step, High, Precision,
     .                  NewPrec, PhiMax
       EXTERNAL PhiMax

       Low = -3.D0
       Step = 0.1D0
       CALL FindStep(Low, Step, PhiMax)
C Search for change of sign
       High = Low+Step
       NewPrec = Precision*0.1D0
C Find root by interpolation
       CALL Interpolate (Low, High, Energy, NewPrec, PhiMax)

       END


       SUBROUTINE Initialise 
C Ask for values of input parameters

       INCLUDE 'globDFT'

       WRITE (6,*) 'Give integration step h'
       READ (5, *) h
       WRITE (6,*) 'Give maximum integration radius'
       READ (5, *) MaxRad

C Charge of helium nuclues
       Z = 2.D0

C Number of integration grid points
       MaxI = INT(MaxRad/h)
C Make sure that MaxRad is an integer times h
       MaxRad = MaxI*h

       PI = 4.D0*DATAN(1.D0)
       Third = 1.D0/3.D0

       END
      

       DOUBLE PRECISION FUNCTION PhiMax(Energy)
C PhiMax returns the value of the radial wave function (times r)
C at the nucleus, r=0. This is used in subroutine FindBound
C to find the bound states. 

       INCLUDE 'globDFT'

       DOUBLE PRECISION PhiNext, PhiStart, Energy, R

       CALL FillFArr(Energy)
       PhiStart = MaxRad*EXP(-Z*MaxRad)
       R = MaxRad-h
       PhiNext = R*EXP(-Z*R)
       CALL Numerov (-h, MaxI, 2, MaxSol, FArr, .FALSE., PhiStart, 
     .               PhiNext, RadArr)
       RadArr(1) = 2*RadArr(2)-RadArr(3)+h*h*FArr(2)*RadArr(2)
       PhiMax = RadArr(1)

       END



       SUBROUTINE SolveRad(Energy)
C Solve the radial Schrodinger equation
C The solution is used to determine the charge density

       INCLUDE 'globDFT'

       INTEGER I

       DOUBLE PRECISION PhiNext, PhiStart, NormFac, Energy, R

       CALL FillFArr(Energy)
       PhiStart = MaxRad*EXP(-Z*MaxRad)
       R = MaxRad-h
       PhiNext = R*EXP(-Z*R)
       CALL Numerov (-h, MaxI, 2, MaxSol, FArr, .FALSE., PhiStart, 
     .               PhiNext, RadArr)
       RadArr(1) = 2*RadArr(2)-RadArr(3)+h*h*FArr(2)*RadArr(2)

       DO I=1, MaxI
         ChDens(I) = RadArr(I)**2
       END DO
       CALL CalcInt(ChDens, NormFac, MaxI, h)

       NormFac = 1.D0/NormFac
       ChDens(1) = 0.D0
       DO I=2, MaxI
         R = (I-1)*h
         ChDens(I) = -NormFac*ChDens(I)/R
C Note that ChDens is in fact -n(r)*r*4*pi !!!
       END DO
       END


       SUBROUTINE CalcHartPot
C Calculate the hartree potential. This is done by radially integrating
C Poisson's equation. The correct boundary equations are met by starting
C at Vh(0) = 0, and by adding the homegeneous solution alpha x r to the result
C in order to fix the potential at MaxRad to Z/r. 
C The charge ZScr is that part of the charge of the screening
C electron which lies in the sphere with radius 
C MaxRad (this charge must be close to 1 if MaxRad is large).
C Note that HartPot is r x potential

       INCLUDE 'globDFT'

       INTEGER I
       DOUBLE PRECISION Alpha, ZScr
       
       CALL NumInhom (h, 0, 0.D0, MaxRad, .FALSE., 0.D0, 
     .                h, ChDens, HartPot, MaxSol)
       ZScr = 1.D0-(MaxRad*MaxRad+2*MaxRad+1.D0)*EXP(-MaxRad)
       Alpha = (ZScr-HartPot(MaxI))/(MaxI-1)
       DO I=1, MaxI
         HartPot(I) = HartPot(I)+Alpha*(I-1)
       END DO
       END



       SUBROUTINE CalcExcPot
C Calculate the exchange potential. This potential is given as 
C -[ 3 u^2(r) /(2 \pi^2 r^2) ]^{1/3}


       INCLUDE 'globDFT'
       DOUBLE PRECISION Fac, R

       INTEGER I
       Fac = (1.5D0/(PI*PI))**Third
       Vx(1) = 0.D0
       DO I=2, MaxI
         R = (I-1)*h
         HartPot(I) = 2*HartPot(I)
         Vx(I) = -Fac*(-ChDens(I)/R)**Third
       END DO
       END



       SUBROUTINE CeperAlder(CorrEn, CorrPot, Dens)
 
       INCLUDE 'globDFT'

       DOUBLE PRECISION Dens, A, B, Beta1, Beta2, Gamma, Rs,
     .                  C, D, CorrPot, CorrEn
       DATA   Gamma,    Beta1,  Beta2,   A,      B, 
     .         C,         D
     .     /-0.1423D0,1.0529D0,0.3334D0,0.0311D0,-0.048D0, 
     .       0.0020D0, -0.0116D0/

       Rs = (0.75D0/(PI*Dens))**Third
       IF (Rs.GE. 1.D0) THEN
         CorrEn = Gamma/(1+Beta1*SQRT(Rs)+Beta2*Rs)
         CorrPot = CorrEn*(1+7.D0/6.D0*Beta1*SQRT(Rs)+4/3.D0*Beta2*
     .                     Rs)/
     .             (1+Beta1*SQRT(Rs)+Beta2*Rs)
       ELSE
         CorrEn = A*LOG(Rs)+B+C*Rs*LOG(Rs)+D*Rs
         CorrPot = A*LOG(Rs)+(B-A/3)+2/3.D0*C*Rs*LOG(Rs)
     .            + (2*D-C)*Rs/3.D0
       ENDIF
       END
       
       


       SUBROUTINE CalcCorrPot
       INCLUDE 'globDFT'

       DOUBLE PRECISION Dens, R, CorrEn, CorrPot

       INTEGER I

       Dens = -2*ChDens(2)/h
       CALL CeperAlder(CorrEn, CorrPot, Dens)
       Vc(1) = CorrPot
       Ec(1) = CorrEn

       DO I=2, MaxI
         R = (I-1)*h
C Factor of 2 because of two electrons
         Dens = -2*ChDens(I)/(R*4*PI)
         CALL CeperAlder(CorrEn, CorrPot, Dens)
         Vc(I) = CorrPot
         Ec(I) = CorrEn
       END DO
       END



       SUBROUTINE FillFArr(Energy)
C Fills the array FArr, appearing in the radial equation
C Phi'' = F Phi with the appropriate values. 
       
       INCLUDE 'globDFT'

       DOUBLE PRECISION R, Energy

       INTEGER I

       DO I=2, MaxI
         r = (I-1)*h
         FArr(I) = -2.d0*(Energy+(Z-HartPot(I))/R-Vx(I)-Vc(I))
       END DO
       END

       SUBROUTINE HartCorrect
C Calculate \int d^r r^2 V_h(r) n (r) , which is needed in the calculation
C of the energy. 
       
       INCLUDE 'globDFT'
  
       INTEGER I
       DOUBLE PRECISION TempArr(MaxSol)

       DO I=1, MaxI
         TempArr(I) = -HartPot(I)*ChDens(I)
       END DO
       CALL CalcInt(TempArr, HartCorr, MaxI, h)
       END


       SUBROUTINE ExcCorrect
C Calculate \int d^r r^2 V_h(r) n (r) , which is needed in the calculation
C of the energy. 
       
       INCLUDE 'globDFT'
  
       INTEGER I
       DOUBLE PRECISION TempArr(MaxSol), r, Fac, FourThirds

       Fac = (1.5D0/(PI*PI))**Third
       FourThirds = 4*Third
       DO I=2, MaxI
         R = (I-1)*h
         TempArr(I) = 0.5D0*Fac*(-ChDens(I)/R)**FourThirds*R*R
       END DO
       CALL CalcInt(TempArr, ExcCorr, MaxI, h)
       ExcCorr = ExcCorr
       END


       SUBROUTINE CorrCorrect
C The name is a bit awkward, but it should not be confused with
C CalcExc by a `strict' compiler (6 chars for an identifier)
C Calculate \int d^r r^2 V_h(r) n (r) , which is needed in the calculation
C of the energy. 
       
       INCLUDE 'globDFT'
  
       INTEGER I
       DOUBLE PRECISION TempArr(MaxSol), r

       DO I=2, MaxI
         R = (I-1)*h
         TempArr(I) = -(Ec(I)-Vc(I))*2*ChDens(I)*R
       END DO
       CALL CalcInt(TempArr, CorrCorr, MaxI, h)
       CorrCorr = CorrCorr
       END


       SUBROUTINE CalcInt (Array, Result, MaxI, h)
C Calculate an integral over the radial coordinate grid. A fourth-order
C integration method is used. 

       IMPLICIT NONE
 
       INTEGER I, MaxI

       DOUBLE PRECISION Array(MaxI), Result, Help, h

       Help = (Array(1)*17.d0+Array(2)*59.d0+ 
     .         Array(3)*43.d0+Array(4)*49.d0)/48.D0

       DO I=5, MaxI-4
         Help = Help + Array(I)
       ENDDO

       Help = (Array(MaxI)*17.d0+Array(MaxI-1)*59.d0+ 
     .         Array(MaxI-2)*43.d0+Array(MaxI-3)*49.d0)/48.D0
     .        + Help
       Result = Help*h
       END
