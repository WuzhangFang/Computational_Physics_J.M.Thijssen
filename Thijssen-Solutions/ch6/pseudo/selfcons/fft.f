      SUBROUTINE FFT_Plan(N)
      IMPLICIT NONE
      INCLUDE "fftw_f77.i"

      INTEGER NArr(3), N

      INTEGER PlanTo, PlanFro

      COMMON /FFTPlan/PlanTo, PlanFro

      NArr(1) = N
      NArr(2) = N
      NArr(3) = N
      CALL FFTWND_F77_Create_Plan(PlanTo, 3, NArr, FFTW_Forward, 
     .                             FFTW_Estimate)
      CALL FFTWND_F77_Create_Plan(PlanFro, 3, NArr, FFTW_Backward, 
     .                             FFTW_Estimate)
      END


      SUBROUTINE Forward_FFT(N, InField, OutField)
      IMPLICIT NONE
      INCLUDE "fftw_f77.i"

      DOUBLE COMPLEX InField(0:N-1,0:N-1,0:N-1), 
     .               OutField(0:N-1,0:N-1,0:N-1)

      INTEGER PlanTo, PlanFro, N

      COMMON /FFTPlan/PlanTo, PlanFro

      CALL FFTWND_F77_One(PlanTo, InField, OutField)
      END


      SUBROUTINE Backward_FFT(N, InField, OutField)
      IMPLICIT NONE
      INCLUDE "fftw_f77.i"

      DOUBLE COMPLEX InField(0:N-1,0:N-1,0:N-1), 
     .               OutField(0:N-1,0:N-1,0:N-1)
      INTEGER PlanTo, PlanFro, N

      COMMON /FFTPlan/PlanTo, PlanFro

      CALL FFTWND_F77_One(PlanFro, InField, OutField)
      OutField = OutField/N**3
      END

      SUBROUTINE DestroyPlans_FFT ()
      IMPLICIT NONE
      INCLUDE "fftw_f77.i"

      INTEGER PlanTo, PlanFro

      COMMON /FFTPlan/PlanTo, PlanFro

      CALL FFTWND_F77_Destroy_Plan(PlanTo)
      CALL FFTWND_F77_Destroy_Plan(PlanFro)
      END
