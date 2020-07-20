      SUBROUTINE GenEig(HMatrix, SMatrix, N, MaxN, Diag)

C Driver subroutine for solving the generalised eigenvalue problem using the
C LAPACK subroutine DSYGV. The generalised eigenvalue problem is of the form
C H C = lambda S C, where H is HMatrix, S is SMatrix. These should both be 
C symmetric. 
C Parameters: HMatrix: matrix with leading dimension MaxN
C             SMatrix: matrix with leading dimension MaxN
C             N:       only the NxN blocks of HMatrix and SMatrix are used
C             MaxN:    leading dimension of HMatrix and SMatrix
C             Diag:    array containing the eigenvalues. Its size should be 
C                      at least N.

      IMPLICIT NONE
      INTEGER MaxWork, MaxN, N, ErrInfo
      PARAMETER (MaxWork = 10000)
      DOUBLE PRECISION SMatrix(MaxN, MaxN), HMatrix(MaxN, MaxN),
     .                 Work(MaxWork), Diag(MaxN)

      IF (3*N-1.GT.MaxWork) THEN
        WRITE (6,*) 'Size of array Work too small in subroutine', 
     .              ' Geneig'
        STOP
      ENDIF
      CALL DSYGV(1, 'v', 'u', N, HMatrix, MaxN, SMatrix, MaxN, Diag, 
     .           Work, MaxWork, ErrInfo)
      IF (ErrInfo.LT.0) THEN
        WRITE (6,*) 'Wrong value of argument ', -ErrInfo, 
     .              ' in subroutine call to DSYGV'
        STOP
      ELSE IF ((ErrInfo.GT.0).AND.(ErrInfo.LE.N)) THEN
        WRITE (6,*) 'Convergence failure routine DSYGV'
        STOP
      ELSE IF (ErrInfo.GT.N) THEN
        WRITE (6,*) 'Leading minor of order', ErrInfo-N,
     .              ' is not positive definite'
        STOP
      ENDIF
      END


      SUBROUTINE DiagMat (Matrix, Diag, N, NSize)
C Driver subroutine for matrix diagonalisation
C of a symmetric matrix, argument 'Matrix'. This is 
C a square matrix of size NSizexNSize, but only the 
C upper left block of size NxN is diagonalised.
C Furthermore, only the lower left triangle of this
C NxN block is needed.
C On exit, the eigenvalues are stored in the array Diag,
C and the eigenvectors are the columns of Matrix

      IMPLICIT NONE

      INTEGER N, NSize, MaxWork, ErrInfo

      PARAMETER (MaxWork = 10000)

      DOUBLE PRECISION Matrix(NSize, NSize), Diag(N),
     .                 Work (MaxWork)

      IF (3*N-1.GT.MaxWork) THEN
        WRITE (6,*) 'Size of array Work too small in subroutine', 
     .              ' Geneig'
        STOP
      ENDIF
      CALL DSYEV('v', 'l', N, Matrix, NSize, Diag,
     .           Work, MaxWork, ErrInfo)
      IF (ErrInfo.LT.0) THEN
        WRITE (6, *) 'Argument ', -ErrInfo, ' of subroutine',
     .               ' DSYEV had illegal value'
      ELSE IF (ErrInfo.GT.0) THEN
        WRITE (6, *) 'Convergence failure LAPACK subroutine ',
     .               'DSYEV'
      END IF
      END
      

      SUBROUTINE CalcVMatrix (Smatrix, VMatrix, N, NSize)
C This routine calculates the matrix Us^{-1/2}, where U is the 
C matrix which diagonalises S, and s is the diagonal form of S.

      IMPLICIT NONE

      INTEGER N, NSize, I

      DOUBLE PRECISION SMatrix(NSize, NSize), Diag(10000),
     .                 VMatrix(NSize, NSize)

      DO I=1, N
        CALL DCOPY (N, SMatrix(1, I), 1, VMatrix(1, I), 1)
      END DO

      CALL DiagMat (VMatrix, Diag, N, NSize)
      DO I=1, N
        Diag(I) = 1/SQRT(Diag(I))
      END DO

C BLAS-1 Routine DCAL is used to calculate Us^{-1/2}
      DO I=1, N
        CALL DSCAL (N, Diag(I), VMatrix(1, I), 1)
      END DO

      END


      SUBROUTINE LowdinDiag(VMatrix, HMatrix, Eigen, 
     .                      Diag, N, NSize)
C This routines uses the matrix V=Us^{-1/2}, where U^\dagger S U = 1,
C to construct H'= V^\dagger H V, and then diagonalises H'.
C It is assumed that only the lower left triangle of HMatrix 
C is filled, i.e. H(I,J) = 0 for I<J

      IMPLICIT NONE

      INTEGER N, NSize, I, MaxWork
      PARAMETER (MaxWork = 10000)

      DOUBLE PRECISION VMatrix(NSize, NSize), HMatrix(NSize, NSize),
     .                 Diag(NSize), Eigen(NSize, NSize),
     .                 Work (MaxWork)

      IF (N*N.GT.SQRT(DBLE(MaxWork))) THEN
        WRITE (6,*) 'Not enough workspace in routine "LowdinDiag"' 
      END IF
      DO I=1, N
        CALL DCOPY (N-I+1, HMatrix(I, I), 1, Eigen(I, I), 1)
        CALL DCOPY (I, HMatrix(I, 1), NSize, Eigen(1, I), 1)
      END DO

      CALL DGEMM ('t', 'n', N, N, N, 1.d0, VMatrix, NSize, Eigen, 
     .            NSize, 0.D0, Work, NSize)
      CALL DGEMM ('n', 'n', N, N, N, 1.d0, Work, NSize, VMatrix, 
     .            NSize, 0.D0, Eigen, NSize)
      CALL DiagMat(Eigen, Diag, N, NSize)
      CALL DGEMM ('n', 'n', N, N, N, 1.D0, VMatrix, NSize, Eigen,
     .            NSize, 0.D0, Work, NSize)
      DO I=1, N
        CALL DCOPY (N, Work(N*(I-1)+1), 1, Eigen(1, I), 1)
      END DO
      
      END

      