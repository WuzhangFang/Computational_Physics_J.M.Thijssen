MODULE calcul

USE globals
USE utilities

CONTAINS
  SUBROUTINE TotalEnergy
! Calculates the total energy for any density
  

  CALL Density_From_Coeffs()
! Kinetic energy
  Fac = (2*PI/(BoxL*MaxBas))**2
  DO I=0, MaxBas
    DO J=0, MaxBas
      DO K=0, MaxBas
        Fac2 = 0.5D0*(I*I+J*J+K*K)*Fac
          EKin = EKin+SUM(Coeffs_K(:,I,J,K))*Fac2
      END DO
    END DO
  END DO

! Local PP energy                      
