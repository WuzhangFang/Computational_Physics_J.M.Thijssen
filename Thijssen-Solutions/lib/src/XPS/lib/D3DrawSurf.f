      SUBROUTINE D3DrawSurf(Psi, Size, MaxSize)
      IMPLICIT NONE
      INTEGER Size, MaxSize
      DOUBLE PRECISION Psi(MaxSize, MaxSize)
      DOUBLE PRECISION x1, y1, x2, y2, z1, z2, z3, z4, Fac
      INTEGER i, j
  
      Fac = 1.D0
      DO I=1, Size-1
        X1 = Fac*I
        X2 = X1 + Fac
        DO J=1, Size-1
          y1 = Fac*j 
          y2 = y1 + Fac
          z1 = Psi(I, J)
          z2 = Psi(I+1, J)
          z3 = Psi(I, J+1)
          z4 = Psi(I+1, J+1)
          CALL SetFastColor(3)
          CALL D3FillRect(x1, y1, z1, x2, y1, z2, 
     .                    x2, y2, z4, x1, y2, z3)
          CALL SetFastColor(6)
          CALL D3DrawLine (x2, y2, z4, x1, y2, z3)
          CALL D3DrawLine (x2, y2, z4, x2, y1, z2)
          CALL D3DrawLine (x1, y1, z1, x1, y2, z3)
          CALL D3DrawLine (x1, y1, z1, x2, y1, z2) 
        END DO
      END DO
      END 
