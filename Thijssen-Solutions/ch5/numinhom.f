
       SUBROUTINE NumInhom (Delta, L, StartR, MaxR, ZeroSing,
     .               PhiStart, PhiNext, RHS, Solution, MaxSol)

C Integrates an inhomogeneous the Schrodinger equation from StartR to MaxR. The initial values
C of the solution are PhiStart and PhiNext resectively. 
C MaxR is also output as the value nearest to MaxR which is an integer
C times the integration step from StartR.
C The output values is the solution, stored in the array "Solution". 
C This array is declared with linear size MaxSol.
C Delta is the integration step. L is the angular momentum quantum number,
C which is passed to the function F(R,L,E) which returns 
C F(R,L) = V(R)-E+h^2 L(L+1)/(2MR^2)

       IMPLICIT NONE

       INTEGER L, I, MaxI, MaxSol

       DOUBLE PRECISION Phi, PhiStart, PhiNext, 
     .        PhiRMax, PhiRNext, W, WNext, WPrev,
     .        Delta, Fac, R, StartR, MaxR, DeltaSq,
     .        RHS(MaxSol), Solution(MaxSol)

       LOGICAL ZeroSing

       MaxI = NINT((MaxR-StartR)/Delta)
       IF (MaxSol.LT.MaxI) THEN
         WRITE (6,*) MaxI, ' integration points' 
         WRITE (6,*) 'Size of array Solution not large enough'
         WRITE (6,*) 'Program stopped in routine Numerov'
         STOP
       END IF

       DeltaSq = Delta*Delta
       Fac = DeltaSq/12.D0

       R  = StartR
       IF (ZeroSing) THEN
         WPrev = PhiStart
         Solution(1) = 0.D0
       ELSE
         WPrev = PhiStart-Fac*RHS(1)
         Solution(1) = PhiStart
       ENDIF

       R  = StartR+Delta
       Phi = PhiNext
       Solution(2) = PhiNext
       W   = Phi-Fac*RHS(2)
       DO I = 2, MaxI+1
c          print *, rhs(i)
          WNext = W*2.D0 - WPrev + DeltaSq*RHS(I)
          WPrev = W
          W     = WNext
          R = R + Delta
          Phi   = W+Fac*RHS(I+1)
          Solution(I+1) = Phi
       ENDDO
       MaxR = R-Delta
       PhiRNext = Phi 
       PhiRMax = WPrev+Fac*RHS(MaxI)

       END 


