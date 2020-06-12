       PROGRAM Exchange
C Program to calculate the l=0 ground state of the helium
C atom in the Hartree-approximation with exchange correction by direct
C integration of the radial Schrodinger equation. 
C Program described in "Computational Physics", J. M. Thijssen
C Section 5.3.3
C Program written by J. M. Thijssen, July 1998

       INCLUDE 'globDFT'

       DOUBLE PRECISION Energy, Precision, OldEnergy
       PARAMETER (Precision = 1.D-9)

       CALL Initialise

       OldEnergy = -1.D3
       Energy = -0.3D0
       DO WHILE(ABS(OldEnergy-Energy).GT.Precision)
         OldEnergy = Energy
         CALL FindBound(Energy, Precision)
         CALL SolveRad(Energy)
         CALL ClcHartCorr
         CALL ClcExcCorr
c         print *, HartCorr*27.2D0, ExcCorr*27.2D0
         print *, 'total energy', 2*Energy-HartCorr-0.5D0*ExcCorr
         CALL CalcHartree
         CALL CalcExc
       END DO

       END


       SUBROUTINE FindBound(Energy, Precision)
C Find the energy of the bound state.
C The bound state is characterised as having a zero at r=0

       IMPLICIT NONE
c       INCLUDE 'globDFT'

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

       PI = 4.D0*DATAN(1.D0)

       END
      

       DOUBLE PRECISION FUNCTION PhiMax(Energy)
C PhiMax returns the value of the radial wave function (times r)
C at the nucleus, r=0. This is used in subroutine FindBound
C to find the bound states. 

       INCLUDE 'globDFT'

       INTEGER I

       DOUBLE PRECISION PhiNext, PhiStart, NormFac, Energy, R

       CALL FillFArr(Energy)
       PhiStart = MaxRad*EXP(-Z*MaxRad)
       R = MaxRad-h
       PhiNext = R*EXP(-Z*R)
       CALL Numerov (-h, 0, MaxRad, MaxI, 3, MaxSol, FArr, PhiStart, 
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
       CALL Numerov (-h, 0, MaxRad, MaxI, 3, MaxSol, FArr, PhiStart, 
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
C Note that ChDens is in fact -n(r)*r !!!
       END DO
       END


       SUBROUTINE CalcHartree
C Calculate the hartree potential. This is done by radially integrating
C Poisson's equation. The correct boundary equations are met by starting
C at Vh(0) = 0, and by adding the homegeneous solution alpha x r to the result
C in order to fix the potential at MaxRad to 1/r. Note that HartPot is r x potential

       INCLUDE 'globDFT'

       INTEGER I
       DOUBLE PRECISION Alpha
       
       CALL NumInhom (h, 0, 0.D0, MaxRad, .FALSE., 0.D0, 
     .                h, ChDens, HartPot, MaxSol)
       Alpha = (1.D0-HartPot(MaxI))/(MaxI-1)
       OPEN (8, file='hartpot')
       DO I=1, MaxI
         HartPot(I) = HartPot(I)+Alpha*(I-1)
         WRITE (8, '(2F12.7)') (I-1)*h, HartPot(I)
       END DO
       close (8)
       END



       SUBROUTINE CalcExc
C Calculate the exchange potential. This potential is given as 
C -[ 3 u^2(r) /(2 \pi^2 r^2) ]^{1/3}


       INCLUDE 'globDFT'
       DOUBLE PRECISION Fac, R

       INTEGER I

       Fac = (1.5D0/(PI*PI))**(0.33333333333D0)
       Vx(1) = 0.D0
       DO I=2, MaxI
         R = (I-1)*h
         HartPot(I) = 2*HartPot(I)
c         WRITE (7, '(2F12.7)') r, ChDens(I)*r
         Vx(I) = -Fac*(-ChDens(I)/R)**0.333333333333D0
c         WRITE (9, '(2F12.7)') r, Vx(I)
       END DO
c       stop
       END


       DOUBLE PRECISION FUNCTION FillFArr(Energy)
C Fills the array FArr, appearing in the radial equation
C Phi'' = F Phi with the appropriate values. 
       
       INCLUDE 'globDFT'

       DOUBLE PRECISION R, Energy

       INTEGER I

       DO I=2, MaxI
         r = (I-1)*h
         FArr(I) = -2.d0*(Energy+(2-HartPot(I))/R-Vx(I))
       END DO

       END

       SUBROUTINE ClcHartCorr
C The name is a bit awkward, but it should not be confused with
C CalcHertree by a `strict' compiler (6 chars for an identifier)
C Calculate \int d^r r^2 V_h(r) n (r) , which is needed in the calculation
C of the energy. 
       
       INCLUDE 'globDFT'
  
       INTEGER I
       DOUBLE PRECISION TempArr(MaxSol), r

       DO I=1, MaxI
         r = (I-1)*h
         TempArr(I) = -HartPot(I)*ChDens(I)
       END DO
       CALL CalcInt(TempArr, HartCorr, MaxI, h)
       END


       SUBROUTINE ClcExcCorr
C The name is a bit awkward, but it should not be confused with
C CalcExc by a `strict' compiler (6 chars for an identifier)
C Calculate \int d^r r^2 V_h(r) n (r) , which is needed in the calculation
C of the energy. 
       
       INCLUDE 'globDFT'
  
       INTEGER I
       DOUBLE PRECISION TempArr(MaxSol), r

       DO I=1, MaxI
         R = (I-1)*h
         TempArr(I) = -Vx(I)*ChDens(I)*R
       END DO
       CALL CalcInt(TempArr, ExcCorr, MaxI, h)
       ExcCorr = 0.5D0*ExcCorr
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
