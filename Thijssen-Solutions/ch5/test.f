       PROGRAM test
C Program to calculate the l=0 ground state of the helium
C atom, using the Hartree-potential of only one electron.
C Program described in "Computational Physics", J. M. Thijssen
C Section 5.3.2
C Program written by J. M. Thijssen, July 1998

       INCLUDE 'globDFT'

       DOUBLE PRECISION Energy, Precision, OldEnergy, R
       PARAMETER (Precision = 1.D-7)
       INTEGER I

       CALL Initialise

       OldEnergy = -1.D3
       Energy = -0.5D0
       CALL SolveRad(Energy)
       CALL CalcHartPot

       DO I=1, MaxI
         R=(I-1)*h
         print *, I, HartPot(I), 1-(R+1)*EXP(-2*R)
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
     .               PhiNext, RadArr, Energy)
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
     .               PhiNext, RadArr, Energy)
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



       SUBROUTINE FillFArr(Energy)
C Fills the array FArr, appearing in the radial equation
C Phi'' = F Phi with the appropriate values. 
       
       INCLUDE 'globDFT'

       DOUBLE PRECISION R, Energy

       INTEGER I

       DO I=2, MaxI
         r = (I-1)*h
         FArr(I) = -2.d0*(Energy+(2-HartPot(I))/R)
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
