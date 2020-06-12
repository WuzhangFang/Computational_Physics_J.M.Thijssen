       DOUBLE PRECISION FUNCTION SphBesJ (L, X)

C Returns the spherical bessel function j_l(x) as a function of l and x
C Upward recursion is used; therefore reliable for a restricted number of
C l-values


       IMPLICIT NONE

       INTEGER L, HelpL

       DOUBLE PRECISION X, JL, JlMin1, JlMin2, HelpSin

       IF (L.EQ.0) THEN
         SphBesJ = SIN(X)/X
       ELSE IF (L.EQ.1) THEN
         SphBesJ = SIN(X)/X/X-COS(X)/X
       ELSE
         HelpSin =  SIN(X)
         JlMin1 = HelpSin/X
         Jl     = HelpSin/X/X-COS(X)/X
         DO HelpL=2, L
           JlMin2 = JlMin1
           JlMin1 = Jl
           Jl = (2*HelpL-1)/X*JlMin1 - JlMin2
         ENDDO
         SphBesJ = Jl
       END IF

       END

      



       DOUBLE PRECISION FUNCTION SphBesN (L, X)
C Returns the spherical bessel function n_l(x) as a function of l and x
C Upward recursion is used.

       IMPLICIT NONE

       INTEGER L, HelpL

       DOUBLE PRECISION X, NL, NlMin1, NlMin2, HelpCos

       IF (L.EQ.0) THEN
         SphBesN = -COS(X)/X
       ELSE IF (L.EQ.1) THEN
         SphBesN = -COS(X)/X/X-SIN(X)/X
       ELSE
         HelpCos =  COS(X)
         NlMin1 = -HelpCos/X
         Nl     = -HelpCos/X/X-SIN(X)/X
         DO HelpL=2, L
           NlMin2 = NlMin1
           NlMin1 = Nl
           Nl = (2*HelpL-1)/X*NlMin1 - NlMin2
         ENDDO
         SphBesN = Nl
       END IF

       END



