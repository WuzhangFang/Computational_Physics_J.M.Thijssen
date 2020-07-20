MODULE globals

IMPLICIT NONE

DOUBLE PRECISION :: PI 

!#################################################
! Geometry and basis set size. Values are set during initialisation
!#################################################

DOUBLE PRECISION :: BoxL, Omega	! Linear Box Size
INTEGER :: MaxBas, GridSize, &	! Maximum number of waves along 
	   DiagSize, MatSize	! LINEAR direction, Matrix sizes

!#################################################
! Physical fields
!#################################################

DOUBLE COMPLEX, ALLOCATABLE :: Density_K(:,:,:)
DOUBLE COMPLEX, ALLOCATABLE :: Density_R(:,:,:)

!#################################################
! Physical system: nr of ions, electrons and orbitals
!#################################################

INTEGER :: N_ion, N_electron, N_orbitals

!#################################################
! Wavefunction coefficients
!#################################################

DOUBLE COMPLEX, ALLOCATABLE :: Coeffs_K(:,:,:,:), Coeffs_R(:,:,:,:)

!#################################################
! Kohn-Sham Hamiltonian, stored values of pseudopot.
! core charges and short range part of local pseudopot
!#################################################

DOUBLE COMPLEX, ALLOCATABLE :: H_KS(:,:), Fixed_KS(:,:), CoreCharge(:,:,:), &
                               ShortLocal(:,:,:), NonLocal(:,:,:,:,:), PseudoGrid(:,:,:)

!#################################################
! Data of ions, stored Coulomb potential and kinetic
! term.
!#################################################

DOUBLE PRECISION, ALLOCATABLE :: Ions(:,:), CoulombGrid(:,:,:), &
                                 KinGrid (:,:,:), FillFac(:)
!#################################################
! Connction between linear indices of H_KS and 
! grid positions 
!#################################################


INTEGER, ALLOCATABLE :: GridIndex(:,:,:), GridPos(:,:)

!#################################################
! Cut-off in reciprocal space
!#################################################

DOUBLE PRECISION :: GMax


!#################################################
! Pseudopotential parameters
!#################################################

TYPE Type_PP_Params
  INTEGER :: AtomNum
  INTEGER :: Zion
  INTEGER :: N_nonzero_C_i
  DOUBLE PRECISION :: Xi
  DOUBLE PRECISION :: C(4)
  DOUBLE PRECISION :: MaxL
  DOUBLE PRECISION :: r_s
  DOUBLE PRECISION :: h_1s, h_2s, r_p, h_1p
END TYPE Type_PP_Params

INTEGER :: No_OF_DIFF_IONS !Has to be assigned in main/InCP

TYPE (Type_PP_Params), ALLOCATABLE :: PP_params(:) 


END MODULE Globals
