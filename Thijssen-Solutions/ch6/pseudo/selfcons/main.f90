PROGRAM CarPar

USE globals
USE utilities
USE pseudo
USE KS
USE excorr

CALL InitParams
CALL SCF


CONTAINS

  SUBROUTINE InitParams
! All parameters are read from file InCP and initialised
    IMPLICIT NONE

    DOUBLE PRECISION :: Energy_CutOff
    INTEGER :: I, I2, I3, I5, N

    PI = 4.D0*ATAN(1.D0)
 
    OPEN (8, File="InCP")
    READ (8, *) BoxL
    Omega = BoxL**3
    READ (8, *) Energy_CutOff
    GMax = SQRT(2.D0*Energy_Cutoff)
    GMax = GMax*BoxL*0.5D0/PI
    print *, 'gmax', 4*gmax
    I2 = Ceiling(log(4*Gmax)/Log(2.D0))
    I2 = 2**I2
    I3 = Ceiling(log(4*Gmax)/Log(3.D0))
    I3 = 3**I3
    I5 = Ceiling(log(4*Gmax)/Log(5.D0))
    I5 = 5**I5
    GridSize = MIN(I2,I3,I5)
!    Gridsize = 4
    MaxBas = INT(GMax)
    ALLOCATE (Density_K(0:GridSize-1,0:GridSize-1,0:GridSize-1))
    ALLOCATE (Density_R(0:GridSize-1,0:GridSize-1,0:GridSize-1))
    CALL FFT_Plan(GridSize)
    print *, 'Per direction, the linear index of PWs runs from' 
    print *, 0, '  to ', GridSize-1
    READ (8, *) No_OF_DIFF_IONS
    print *, 'there are ',  No_OF_DIFF_IONS, ' different ions'
    ALLOCATE (PP_Params(No_OF_DIFF_IONS))
    PP_Params(:)%AtomNum = 0
    READ (8, *) N_ion
    ALLOCATE (Ions(N_ion,4))
    DO I=1, N_ion
      READ (8,*) Ions(I, 1), &  ! Atomic Number
                 Ions(I,2), Ions(I,3), Ions(I,4) ! Positions
      CALL Init_PP_Params(NINT(Ions(I,1)))
    END DO
    print *, 'No of ions', N_ion
    print *, 'Ion data', Ions  
    CALL FillGrids()

    READ (8, *) N_electron
    READ (8, *) N_Orbitals
    ALLOCATE (FillFac(N_orbitals))
    DO N = 1, N_orbitals
      READ (8,*) FillFac(N)
    END DO
    ALLOCATE(Coeffs_R(N_orbitals, 0:GridSize-1, 0:GridSize-1, 0:GridSize-1))
    ALLOCATE(Coeffs_K(N_orbitals, 0:GridSize-1, 0:GridSize-1, 0:GridSize-1))
    CALL InitSol()
  END SUBROUTINE InitParams



  SUBROUTINE SCF
! The actual SCF procedure is called
  IMPLICIT NONE

  ALLOCATE (H_KS(MatSize, MatSize))
  ALLOCATE (Fixed_KS(MatSize, MatSize))
  CALL KohnShamLoop()
  END SUBROUTINE SCF


  SUBROUTINE FillGrids()
! Convenient storage of pseudopotential, G^2, G^{-2} etcetera
    IMPLICIT NONE

    INTEGER :: I, J, K, N, R2, Index, I1, J1, K1, PP_Index, AtomNum
    DOUBLE COMPLEX :: StructFac, Im

    MatSize = 0
    DO I=0, GridSize-1
      DO J=0, GridSize-1
        DO K=0, GridSize-1
          R2 = R2_Short(I,J,K)
          IF (R2<Gmax**2) MatSize = MatSize + 1
        END DO
      END DO
    END DO
    ALLOCATE (GridIndex(0:GridSize-1, 0:GridSize-1, 0:GridSize-1))
    ALLOCATE (GridPos(MatSize,3))
    ALLOCATE (PseudoGrid(0:GridSize-1, 0:GridSize-1, 0:GridSize-1))
    ALLOCATE (CoreCharge(0:GridSize-1, 0:GridSize-1, 0:GridSize-1))
    ALLOCATE (ShortLocal(0:GridSize-1, 0:GridSize-1, 0:GridSize-1))
    ALLOCATE (CoulombGrid(0:GridSize-1, 0:GridSize-1, 0:GridSize-1))
    ALLOCATE (KinGrid(0:GridSize-1, 0:GridSize-1, 0:GridSize-1))
    ALLOCATE (NonLocal(0:GridSize-1, 0:GridSize-1, 0:GridSize-1,N_ion,5))
    GridIndex = 0
    Index = 0
    PseudoGrid = 0.D0
    ShortLocal = 0.D0
    CoreCharge = 0.D0
    NonLocal = 0.D0
    Im = CMPLX(0.D0, 1.D0)
    DO I=0, GridSize-1
      DO J=0, GridSize-1
        DO K=0, GridSize-1
          R2 = R2_Short(I,J,K)
          KinGrid(I,J,K) = R2
          IF (R2/=0) THEN
            CoulombGrid(I,J,K) = 1.D0/R2
          ELSE
            CoulombGrid(I,J,K) = 0.D0
          END IF
          DO N = 1, N_ion
            AtomNum = Ions(N,1)
            PP_Index = GetIndexPP(AtomNum)
            I1 = I-INT((2.D0*I)/GridSize)*GridSize
            J1 = J-INT((2.D0*J)/GridSize)*GridSize
            K1 = K-INT((2.D0*K)/GridSize)*GridSize
            CALL CalcFacs(I, J, K, N, R2, StructFac)
            PseudoGrid(I,J,K) = PseudoGrid(I,J,K) + Local_PP(R2,NINT(Ions(N,1)))*StructFac
!            print *, I, J, K, PseudoGrid(I,J,K)
            ShortLocal(I,J,K) = ShortLocal(I,J,K) + Short_PP(R2,NINT(Ions(N,1)))*StructFac
            CoreCharge(I,J,K) = CoreCharge(I,J,K) + CoreDens(R2,NINT(Ions(N,1)))*StructFac
            NonLocal(I,J,K,N,1) = NonLoc(I1,J1,K1,AtomNum,0,0,1)*StructFac
            IF (PP_Params(PP_Index)%MaxL>0) THEN
              NonLocal(I,J,K,N,2) = NonLoc(I1,J1,K1,AtomNum,0,0,2)*StructFac
              NonLocal(I,J,K,N,3) = NonLoc(I1,J1,K1,AtomNum,1,1,1)*StructFac
              NonLocal(I,J,K,N,4) = NonLoc(I1,J1,K1,AtomNum,1,0,1)*StructFac
              NonLocal(I,J,K,N,5) = NonLoc(I1,J1,K1,AtomNum,1,-1,1)*StructFac
            END IF
          END DO
          IF (R2<Gmax**2) THEN
            Index = Index + 1
            GridIndex(I, J, K) = Index
            GridPos(Index,:) = (/I,J,K/)
          END IF
        END DO
      END DO
    END DO
  END SUBROUTINE
          
        
  
  SUBROUTINE InitSol
! Starting solution, Gaussian distribution
  IMPLICIT NONE
  INTEGER I, J, K, ElecCnt, II, JJ, KK
  DOUBLE PRECISION :: Norm, Alpha=0.5D0, R2, X, Y
  
  Coeffs_K = (0.D0, 0.D0)
  DO I=-MaxBas, MaxBas
    DO J=-MaxBas, MaxBas
      DO K=-MaxBas, MaxBas
        R2 = (I*I+J*J+K*K)/(MaxBas+1)**2
        DO ElecCnt = 1, N_orbitals
          Norm = EXP(-Alpha*R2)
          CALL Random_Number(X)
          CALL Random_Number(Y)
          IF (I<0) THEN
            II = I + GridSize
          ELSE
            II = I
          END IF
          IF (J<0) THEN
            JJ = J + GridSize
          ELSE
            JJ = J
          END IF
          IF (K<0) THEN
            KK = K + GridSize
          ELSE
            KK = K
          END IF

          Coeffs_K(ElecCnt, II, JJ, KK) = CMPLX(Norm*X, Norm*Y)
        END DO
      END DO
    END DO  
  END DO
  CALL Gram_Schmidt()
  
  END SUBROUTINE InitSol
          
          
          

END PROGRAM CarPar
    
