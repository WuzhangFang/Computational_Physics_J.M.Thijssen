MODULE KS
! In this module, the Kohn-Sham matrix is set up and diagonalised.
! The density is then fed into a new KS Hamiltonian etcetera
USE pseudo
USE excorr
USE Globals
USE utilities
! USE totener

CONTAINS

  SUBROUTINE FillFixed_KS()
! Fixed part of Kohn-Sham matrix
  IMPLICIT NONE

  INTEGER :: IIndex, JIndex, I1, J1, K1, I2, J2, K2, AtomNum, &
             I, J, K, N, PP_Index
  DOUBLE PRECISION :: PreFac, h_1s, h_2s, h_1p
  DOUBLE COMPLEX, ALLOCATABLE :: TempVec(:,:,:,:,:)
  DOUBLE COMPLEX :: Im, StructFac

  ALLOCATE (TempVec(0:GridSize-1, 0:GridSize-1, 0:GridSize-1,N_ion,5))

  Im = CMPLX(0.D0, 1.D0)

  Fixed_KS=CMPLX(0.D0)

! KINETIC TERM

  PreFac = CMPLX(2*PI*PI/BoxL**2)
  DO IIndex = 2, MatSize
    Fixed_KS(IIndex, IIndex) = PreFac*KinGrid(GridPos(IIndex,1),&
                                              GridPos(IIndex,2),&
                                              GridPos(IIndex,3))
  END DO

  DO N = 1, N_ion
    AtomNum = NINT(Ions(N,1))
    IF (AtomNum > 4) THEN
      PP_Index = GetIndexPP(AtomNum)
      h_1s = PP_Params(PP_Index)%h_1s
      h_2s = PP_Params(PP_Index)%h_2s
      h_1p = PP_Params(PP_Index)%h_1p
      DO IIndex = 1, MatSize
        DO JIndex = 1, MatSize
          I1 = GridPos(IIndex,1)
          J1 = GridPos(IIndex,2)
          K1 = GridPos(IIndex,3)
          I2 = GridPos(JIndex,1)
          J2 = GridPos(JIndex,2)
          K2 = GridPos(JIndex,3)
          Fixed_KS(IIndex, JIndex) = Fixed_KS(IIndex, JIndex) + &
            h_1s*(NonLocal(I1,J1,K1,N,1))*CONJG(NonLocal(I2,J2,K2,N,1))+&
            h_2s*(NonLocal(I1,J1,K1,N,2))*CONJG(NonLocal(I2,J2,K2,N,2))+&
            h_1p*(NonLocal(I1,J1,K1,N,3))*CONJG(NonLocal(I2,J2,K2,N,3))+&
            h_1p*(NonLocal(I1,J1,K1,N,4))*CONJG(NonLocal(I2,J2,K2,N,4))+&
            h_1p*(NonLocal(I1,J1,K1,N,5))*CONJG(NonLocal(I2,J2,K2,N,5))
        END DO
      END DO
    END IF
  END DO   
  END SUBROUTINE FillFixed_KS


  SUBROUTINE Fill_KS()
! Kohn-Sham matrix is filled

  IMPLICIT NONE
  INTEGER :: I1, J1, K1, I2, J2, K2, IIndex, JIndex, N, IT, JT, KT, &
             R2 ! , R2_Short
  DOUBLE COMPLEX, ALLOCATABLE :: TempVec(:,:,:), TempVec2_R(:,:,:), &
                                 TempVec2_K(:,:,:)
  DOUBLE COMPLEX :: PreFac, II

! ALLOCATE STORAGE
  IF (.NOT.ALLOCATED(TempVec)) THEN 
    ALLOCATE (TempVec(0:GridSize-1, 0:GridSize-1, 0:GridSize-1))
  END IF
  IF (.NOT.ALLOCATED(TempVec2_R)) THEN 
    ALLOCATE (TempVec2_R(0:GridSize-1, 0:GridSize-1, 0:GridSize-1))
  END IF

  H_KS=CMPLX(0.D0)

! KINETIC TERM
  H_KS=Fixed_KS

! EXCHANGE CORRELATION 
  TempVec = CMPLX(0.D0) 
  TempVec2_R = CMPLX(0.D0) 
  print *, gridsize
  DO I1 = 0, GridSize-1
    DO J1 = 0, GridSize-1
      DO K1 = 0, GridSize-1
        TempVec2_R(I1, J1, K1) = CMPLX(Vxc(I1, J1, K1))
!        print *, i1, j1, k1, Vxc(I1, J1, K1)
      END DO
    END DO
  END DO
  CALL BackWard_FFT(GridSize, TempVec2_R, TempVec)

! LOCAL PSEUDOPOTENTIAL
  TempVec = TempVec + PseudoGrid

! HARTREE
  PreFac = BoxL**2/PI
  TempVec = TempVec + PreFac*Density_K*CoulombGrid

! ADD POTENTIAL TO H_KS:

  IIndex = 0
  DO IIndex = 1, MatSize
    DO JIndex = 1, MatSize
!    DO JIndex = IIndex, MatSize
      I1 = GridPos(IIndex,1)
      J1 = GridPos(IIndex,2)
      K1 = GridPos(IIndex,3)
      I2 = GridPos(JIndex,1)
      J2 = GridPos(JIndex,2)
      K2 = GridPos(JIndex,3)
      IT = MOD(I1-I2+GridSize, GridSize)
      JT = MOD(J1-J2+GridSize, GridSize)
      KT = MOD(K1-K2+GridSize, GridSize)
      H_KS(IIndex,JIndex) = H_KS(IIndex, JIndex) + &
                            TempVec(IT, JT, KT)
    END DO
  END DO
  END SUBROUTINE Fill_KS



  SUBROUTINE Total_Energy

! Evaluate the energy for a given solution Coeffs_K

  IMPLICIT NONE
  INTEGER :: I1, J1, K1, I2, J2, K2, IIndex, JIndex, R2, N, IT, JT, KT, &
             IElec, M, PP_Index, AtomNum, IOrb
  DOUBLE COMPLEX, ALLOCATABLE :: TempVec(:,:,:), TempVec2_R(:,:,:), &
                                 TempVec2_K(:,:,:)
  DOUBLE COMPLEX :: II, E_KS, n_tot, Norm, Etemp, PreFac, shortloc
  DOUBLE PRECISION :: RPos, h_1s, h_2s, h_1p

  E_KS=CMPLX(0.D0)
  II = (0.D0, -1.D0)
  Norm = CMPLX(0.D0)

! KINETIC ENERGY

  print *, 'nr el. ', N_electron
  PreFac = 2*PI*PI/BoxL**2
  DO N = 1, N_orbitals
    E_KS = E_KS + FillFac(N)*PreFac*&
        SUM(Coeffs_K(N,:,:,:)*CONJG(Coeffs_K(N,:,:,:))*KinGrid)
  END DO
  print '(A23 F15.8)', 'kin', DBLE(e_ks)
  
! SHORT RANGE PART OF LOCAL PP
  
  ETemp = Omega*SUM(ShortLocal*CONJG(Density_K))
  print '(A23 F15.8)', 'pp_sr', DBLE(etemp)
  E_KS=E_KS+Etemp

  Etemp = Omega*SUM(PseudoGrid*CONJG(Density_K))

  print '(A23 F15.8)', 'loc pp', DBLE(etemp)

! EXCHANGE CORRELATION ENERGY
  Etemp = CMPLX(0.D0)
  DO I1 = 0, GridSize-1
    DO J1 = 0, GridSize-1
      DO K1 = 0, GridSize-1
! Evaluate this expression in real space!!!!!!!!!!!!
        Etemp = Etemp+CONJG(Density_R(I1, J1, K1))*&
               CMPLX(epsilon_xc(I1,J1,K1))*Omega/GridSize**3
      END DO
    END DO
  END DO

  print '(A23 F15.8)', 'xc', DBLE(etemp)
  E_KS=E_KS+Etemp

! HARTREE
  PreFac = BoxL**2*Omega/(2*PI)
  ETemp = PreFac*SUM(CoulombGrid*Density_K*CONJG(Density_K))
  print '(A23 F15.8)', 'Hartree energy 1:', DBLE(Etemp)

! Nonlocal PsP
  Etemp = CMPLX(0.D0)
  DO IOrb = 1, N_orbitals
    DO N = 1, N_ion
      AtomNum = NINT(Ions(N,1))
      IF (AtomNum > 4) THEN
        PP_Index = GetIndexPP(AtomNum)
        h_1s = PP_Params(PP_Index)%h_1s
        h_2s = PP_Params(PP_Index)%h_2s
        h_1p = PP_Params(PP_Index)%h_1p
        PreFac = SUM(NonLocal(:,:,:,N,1)*CONJG(Coeffs_K(IOrb,:,:,:)))
        Etemp  = Etemp + FillFac(IOrb)*h_1s*PreFac*CONJG(PreFac)
        IF (PP_Params(PP_Index)%MaxL>0) THEN
          PreFac = SUM(NonLocal(:,:,:,N,2)*CONJG(Coeffs_K(IOrb,:,:,:)))
          Etemp  = Etemp + FillFac(IOrb)*h_2s*PreFac*CONJG(PreFac)
          PreFac = SUM(NonLocal(:,:,:,N,3)*CONJG(Coeffs_K(IOrb,:,:,:)))
          Etemp  = Etemp + FillFac(IOrb)*h_1p*PreFac*CONJG(PreFac)
          PreFac = SUM(NonLocal(:,:,:,N,4)*CONJG(Coeffs_K(IOrb,:,:,:)))
          Etemp  = Etemp + FillFac(IOrb)*h_1p*PreFac*CONJG(PreFac)
          PreFac = SUM(NonLocal(:,:,:,N,5)*CONJG(Coeffs_K(IOrb,:,:,:)))
          Etemp  = Etemp + FillFac(IOrb)*h_1p*PreFac*CONJG(PreFac)
        END IF
      END IF
    END DO  
  END DO
  print '(A23 F15.8)', 'Nonlocal Psp:', DBLE(Etemp)
  E_KS=E_KS+Etemp

! CORE
  PreFac = BoxL**2*Omega/(2*PI)
  ETemp = PreFac*SUM(CoulombGrid*CoreCharge*CONJG(CoreCharge))

  print '(A23 F15.8)', 'Local core energy:', DBLE(Etemp)

  ETemp = PreFac*SUM(CoulombGrid*(CoreCharge+Density_K)*CONJG(CoreCharge+Density_K))

  print '(A23 F15.8)', 'self-energy 1:', DBLE(Etemp)

  E_KS=E_KS+Etemp
  Etemp=CMPLX(0.D0)
  Etemp = Etemp+E_ovrl() 
  print '(A23 F15.8)', 'ovrl:', DBLE(Etemp)
  E_KS = E_KS + Etemp
  DO N=1, N_ion
    E_KS=E_KS-E_Self(NINT(Ions(N,1)))
  END DO
  print '(A23 F15.8)', 'Total energy:', DBLE(E_KS)
  END SUBROUTINE Total_Energy
  



  SUBROUTINE  KohnShamLoop()
! SCF iterations
  IMPLICIT NONE
  INTEGER :: N, LDWork, INFO, Iter, MaxIter=12, Ielec, I, J, K
  DOUBLE COMPLEX, ALLOCATABLE :: Work1(:), OldDens_K(:,:,:)
  DOUBLE PRECISION, ALLOCATABLE :: Work2(:), Backup(:,:), Vec(:)
  DOUBLE PRECISION, ALLOCATABLE :: EigVals(:)
  DOUBLE PRECISION :: Norm, R2

  DiagSize = MatSize
  IF (.NOT.ALLOCATED(Work1)) THEN
    ALLOCATE (Work1(2*MatSize-1))
  END IF
  IF (.NOT.ALLOCATED(Work2)) THEN
    ALLOCATE (Work2(3*MatSize-2))
  END IF
  IF (.NOT.ALLOCATED(EigVals)) THEN
    ALLOCATE (EigVals(MatSize))
  END IF
  IF (.NOT.ALLOCATED(OldDens_K)) THEN
    ALLOCATE (OldDens_K(0:GridSize-1,0:GridSize-1,0:GridSize-1))
  END IF
  Density_K = 0.d0
  Density_K(0,0,0) = 1.d0
  Norm =  Omega*Density_K(0,0,0)/N_electron
  Density_K = Density_K/Norm
  CALL Forward_FFT(GridSize, Density_K, Density_R)
!  CALL Density_From_Coeffs()
  LDWork = 2*MatSize-1
  CALL FillFixed_KS()
  DO Iter = 1, MaxIter
    CALL Fill_KS()
    print *, DiagSize, MatSize
!    do i= 1, MatSize
!      do j=1, MatSize
!        print *, i, j, H_KS(i,J)
!      end do
!    end do
    print *, 'grs', gridsize
    CALL ZHEEV('V', 'U', DiagSize, H_KS, MatSize, EigVals, Work1, LDWork, Work2, INFO)
    print '(A15)', 'eigenvalues'
    print '(F15.8)',  eigvals(1)
    CALL FillNewCoeffs(Coeffs_K)
    OldDens_K = Density_K
    CALL Density_From_Coeffs()
!    do i=0, gridsize-1
!    do j=0, gridsize-1
!    do k=0, gridsize-1
!      print '(I4 I4 I4 F10.5)', i, j, k, dble(density_R(i, j, k))
!    end do
!    end do
!    end do  
    CALL Total_Energy
! Mixing:
!    Density_K = 0.5D0*OldDens_K+0.5D0*Density_K
    Norm =  Omega*Density_K(0,0,0)/N_electron
    Density_K = Density_K/Norm
    CALL Forward_FFT(GridSize, Density_K, Density_R)
!    print *, "checknorm", sum(density_r)*OMEGA/GRIDSIZE**3
  END DO


    
  END SUBROUTINE  KohnShamLoop



  SUBROUTINE FillNewCoeffs(NewCoeffs_K)
! Store the lowest eigenvectors of H_KS in NewCoeffs_K
  IMPLICIT NONE

  DOUBLE COMPLEX, INTENT(INOUT) :: NewCoeffs_K(1:N_orbitals,& 
                     0:GridSize-1, 0:GridSize-1, 0:GridSize-1)

  DOUBLE COMPLEX :: TotDens=(0.d0,0.d0)

  INTEGER :: IElec, I, J, K, LinCnt, R2

!  print *, 'nu!', N_orbitals, N_electron
  DO IElec = 1, N_orbitals
    NewCoeffs_K(IElec,:,:,:) = CMPLX(0.D0)
    DO LinCnt = 1, MatSize
      I = GridPos(LinCnt,1)
      J = GridPos(LinCnt,2)
      K = GridPos(LinCnt,3)
      if (ielec .eq. 1)  print *, LinCnt, i,j,k
      NewCoeffs_K(IElec,I,J,K) = H_KS(LinCnt,IElec)
    END DO
  END DO
  END SUBROUTINE FillNewCoeffs
  
END MODULE KS

