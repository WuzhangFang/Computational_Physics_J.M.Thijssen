       SUBROUTINE FindStep (MinX, Step, F)
       IMPLICIT NONE
       DOUBLE PRECISION MinX, Step, FPrev, F

       FPrev = F(MinX)
       DO WHILE (F(MinX+Step)*FPrev.GT.0) 
         MinX = MinX + Step         
       END DO
       END      





       SUBROUTINE Interpolate (MinX, MaxX, IntX, IntWidth, F)

       IMPLICIT NONE

       DOUBLE PRECISION MinX, MaxX, IntX, IntWidth, LeftRes, F,
     .               RightRes, MidX, Wronski, LeftX, RightX, MidRes

       LeftX = MinX
       RightX = MaxX
       LeftRes = F (MinX)
       RightRes = F (MaxX)
       MidRes = 1.d0
       DO WHILE (ABS(MidRes).GT.IntWidth)
         MidX = (RightX*LeftRes-LeftX*RightRes)/
     .          (LeftRes-RightRes)
c         print *, 'te', midx, rightx, leftx, rightres, leftres
         MidRes = F(MidX)
         IF (MidRes*LeftRes.GT.0.d0) THEN
           LeftRes = MidRes
           LeftX = MidX
         ELSE
           RightRes = MidRes
           RightX = MidX
         ENDIF
       ENDDO
       IntX = MidX
       END

