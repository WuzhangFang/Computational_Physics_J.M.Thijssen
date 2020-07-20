MODULE utilities

USE globals

CONTAINS

  SUBROUTINE Gram_Schmidt()
  IMPLICIT NONE
  INTEGER ElecCnt, ElecCnt2
  DOUBLE COMPLEX :: IP

  DO ElecCnt = 1, N_orbitals
    DO ElecCnt2 = 1, ElecCnt-1
      IP = InnerProd(Coeffs_K(ElecCnt2,:,:,:),Coeffs_K(ElecCnt,:,:,:))
      Coeffs_K(ElecCnt,:,:,:) = Coeffs_K(ElecCnt,:,:,:)-&
          IP*Coeffs_K(ElecCnt2,:,:,:)
    END DO
    IP = InnerProd(Coeffs_K(ElecCnt,:,:,:),Coeffs_K(ElecCnt,:,:,:))
    IP = 1/SQRT(IP)
    Coeffs_K(ElecCnt,:,:,:) = Coeffs_K(ElecCnt,:,:,:)*IP
  END DO
  DO ElecCnt = 1, N_orbitals
    DO ElecCnt2 = 1, ElecCnt
      IP = InnerProd(Coeffs_K(ElecCnt,:,:,:),Coeffs_K(ElecCnt2,:,:,:))
    END DO
  END DO
  DO ElecCnt = 1, N_orbitals
    CALL Forward_FFT(GridSize, Coeffs_K(ElecCnt,:,:,:), &
                     Coeffs_R(ElecCnt,:,:,:))
  END DO
  DO ElecCnt = 1, N_orbitals
    IP = InnerProd(Coeffs_R(ElecCnt,:,:,:),Coeffs_R(ElecCnt,:,:,:))
    DO ElecCnt2 = 1, ElecCnt
      IP = InnerProd(Coeffs_R(ElecCnt,:,:,:),Coeffs_R(ElecCnt2,:,:,:))
!      print *, ElecCnt, ElecCnt2, IP
    END DO
  END DO

  END SUBROUTINE Gram_Schmidt



  DOUBLE COMPLEX FUNCTION InnerProd(Arr1, Arr2)

  IMPLICIT NONE

  DOUBLE COMPLEX :: Arr1(0:GridSize-1,0:GridSize-1,0:GridSize-1), &
                    Arr2(0:GridSize-1,0:GridSize-1,0:GridSize-1)
  InnerProd = SUM(CONJG(Arr1)*Arr2)
  END FUNCTION InnerProd  



  SUBROUTINE Density_From_Coeffs
! Calculates the density_K from Coeffs_K
  IMPLICIT NONE
  
  INTEGER :: ElecCnt, I, N

  DO ElecCnt=1, N_orbitals
    CALL Forward_FFT(GridSize, Coeffs_K(ElecCnt,:,:,:), &
                     Coeffs_R(ElecCnt,:,:,:))
  END DO
  Coeffs_R = Coeffs_R/SQRT(Omega)
  Density_R = CMPLX(0.D0)
  DO N = 1, N_orbitals
    Density_R = Density_R + FillFac(N)*Coeffs_R(N,:,:,:)*CONJG(Coeffs_R(N,:,:,:))
  END DO
  CALL Backward_FFT(GridSize, Density_R, Density_K)
  END SUBROUTINE Density_From_Coeffs



  DOUBLE PRECISION FUNCTION R2_Short (I,J,K)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: I, J, K
  INTEGER :: II, JJ, KK

  II = I-INT((2.D0*I)/GridSize)*GridSize
  JJ = J-INT((2.D0*J)/GridSize)*GridSize
  KK = K-INT((2.D0*K)/GridSize)*GridSize

  R2_Short = II*II+JJ*JJ+KK*KK
  END FUNCTION R2_Short




  SUBROUTINE CalcFacs(I,J,K,N,R2,StructFac)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: I, J, K, N
  INTEGER, INTENT(OUT) :: R2
  DOUBLE COMPLEX, INTENT(OUT) :: StructFac
  INTEGER :: II, JJ, KK
  DOUBLE COMPLEX :: Im=(0.d0, 1.d0)

  II = I-INT((2.D0*I)/GridSize)*GridSize
  JJ = J-INT((2.D0*J)/GridSize)*GridSize
  KK = K-INT((2.D0*K)/GridSize)*GridSize

  R2 = II*II+JJ*JJ+KK*KK
  StructFac = EXP(Im*2*PI*(Ions(N,2)*II+Ions(N,3)*JJ+ &
                              Ions(N,4)*KK)/BoxL)
  END SUBROUTINE CalcFacs

END MODULE Utilities
