       SUBROUTINE Numerov (Delta, StartI, EndI, MaxSol, FArr,
     .                     Sing, PhiStart, PhiNext, Solution)

C Integrates the Schrodinger equation. The initial values
C of the solution are PhiStart and PhiNext resectively. The integration
C step is Delta, and the integration steps run from StartI to EndI. StartI may be larger
C than EndI; in that case, integration is performed backward.
C The output values is the solution, stored in the array "Solution".
C Sing determines whether the potential contains a singularity at
C r=0. 
C If there is a singularity at r=0, the value of the Numerov 
C function w at r=0 is taken equal to PhiStart, and not 
C equal to PhiStart/(1-h^2 FArr/12).

C This array is declared with linear size MaxSol.
C Delta is the integration step.
C The equation solved is
C Psi''(R_I) = FArr(I) Psi(R_I)
C FArr must therefore be filled with the appropriate values before
C calling the present routine. In the case of the radial 
C Schrodinger equation, FArr would contain the values
C FArr(I) = 2*(V(R)-E)+L*(L+1)/R**2 for R=R_I.

       IMPLICIT NONE

       INTEGER I, StartI, EndI, MaxSol, IStep

       DOUBLE PRECISION Phi, PhiStart, PhiNext, W, WNext, WPrev,
     .        Delta, Fac, FArr(MaxSol), DeltaSq, Solution(MaxSol)

       LOGICAL Sing

       IF (Delta.LT.0) THEN
         IStep = -1
       ELSE
         IStep = 1
       ENDIF

       DeltaSq = Delta*Delta
       Fac = DeltaSq/12.D0

       IF (Sing) THEN
         WPrev = PhiStart
       ELSE
         WPrev = (1-Fac*FArr(StartI))*PhiStart
         Solution(StartI) = PhiStart
       ENDIF

       Phi = PhiNext
       Solution(StartI+IStep) = PhiNext
       W   = (1-Fac*FArr(StartI+IStep))*Phi

       DO I = StartI+IStep, EndI-IStep, IStep
          WNext = W*2.D0 - WPrev + DeltaSq*Phi*FArr(I)
          WPrev = W
          W     = WNext
          Phi   = W/(1-Fac*FArr(I+IStep))
          Solution(I+IStep) = Phi
       ENDDO
       END 


