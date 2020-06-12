       PROGRAM hydro
C Program to calculate the l=0 ground state of the hydrogen
C atom, by radial integration of the Schrodinger equation.
C Program described in "Computational Physics", J. M. Thijssen
C Section 5.3.1
C Program written by J. M. Thijssen, July 1998

       INCLUDE 'globDFT'

       DOUBLE PRECISION Energy, Precision
       PARAMETER (Precision = 1.D-7)

       CALL Initialise

       Energy = -0.3D0

       CALL FindBound(Energy, Precision)
       CALL SolveRad(Energy)
       print *, Energy

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
       Z = 1.D0

C Number of integration grid points
       MaxI = INT(MaxRad/h)
C Make sure that MaxRad is an integer times h
       MaxRad = MaxI*h

       PI = 4.D0*DATAN(1.D0)

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
       PhiMax = 2*RadArr(2)-RadArr(3)+h*h*FArr(2)*RadArr(2)

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


       SUBROUTINE FillFArr(Energy)
C Fills the array FArr, appearing in the radial equation
C Phi'' = F Phi with the appropriate values. 
       
       INCLUDE 'globDFT'

       DOUBLE PRECISION R, Energy

       INTEGER I

       DO I=2, MaxI
         r = (I-1)*h
         FArr(I) = -2.d0*(Energy+Z/R)
       END DO

       END


       SUBROUTINE CalcInt (Array, Result, MaxI, h)
C Calculate an integral over the radial coordinate grid. A fourth-order
C integration method is used. 

       IMPLICIT NONE
 
       INTEGER I, MaxI

       DOUBLE PRECISION Array(MaxI), Result, Help, h

       Help = (Array(1)*17.d0+Array(2)*59.d0+ 
     .         Array(3)*43.d0+Array(4)*49.d0)/48.D0

       DO I=4, MaxI-4
         Help = Help + Array(I)
       ENDDO

       Help = (Array(MaxI)*17.d0+Array(MaxI-1)*59.d0+ 
     .         Array(MaxI-2)*43.d0+Array(MaxI-3)*49.d0)/48.D0
     .        + Help
       Result = Help*h
       END
