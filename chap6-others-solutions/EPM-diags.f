c    This matrix diagonalization code was obtained from 
c    www.netlib.org/eispack
c
      subroutine ch(nm,n,ar,ai,w,matz,zr,zi,fv1,fv2,fm1,ierr)
c
      integer i,j,n,nm,ierr,matz
      real ar(nm,n),ai(nm,n),w(n),zr(nm,n),zi(nm,n),
     x       fv1(n),fv2(n),fm1(2,n)
c
c     this subroutine calls the recommended sequence of
c     subroutines from the eigensystem subroutine package (eispack)
c     to find the eigenvalues and eigenvectors (if desired)
c     of a complex hermitian matrix.
c
c     on input-
c
c        nm  must be set to the row dimension of the two-dimensional
c        array parameters as declared in the calling program
c        dimension statement,
c
c        n  is the order of the matrix  a=(ar,ai),
c
c        ar  and  ai  contain the real and imaginary parts,
c        respectively, of the complex hermitian matrix,
c
c        matz  is an integer variable set equal to zero if
c        only eigenvalues are desired,  otherwise it is set to
c        any non-zero integer for both eigenvalues and eigenvectors.
c
c     on output-
c
c        ar and ai contain information about the unitary transformations
c        used in the reduction to tridiagonal form in their full
c        lower triangles.  their
c        strict upper triangles and the diagonal of ar are unaltered,
c
c        w  contains the eigenvalues in ascending order,
c
c        zr  and  zi  contain the real and imaginary parts,
c        respectively, of the eigenvectors if matz is not zero,
c
c        ierr  is an integer output variable set equal to an
c        error completion code described in section 2b of the
c        documentation.  the normal completion code is zero,
c
c        fv1, fv2, and  fm1  are temporary storage arrays.
c
c     questions and comments should be directed to b. s. garbow,
c     applied mathematics division, argonne national laboratory
c
c     ------------------------------------------------------------------
c
c
c
      if (n .le. nm) go to 10
      ierr = 10 * n
      go to 50
c
   10 call  htridi(nm,n,ar,ai,w,fv1,fv2,fm1)
      if (matz .ne. 0) go to 20
c     .......... find eigenvalues only ..........
      call  tqlrat(n,w,fv2,ierr)
      go to 50
c     .......... find both eigenvalues and eigenvectors ..........
   20 do 40 i = 1, n
c
         do 30 j = 1, n
            zr(j,i) = 0.0
   30    continue
c
         zr(i,i) = 1.0
   40 continue
c
      call  tql2(nm,n,w,fv1,zr,ierr)
      if (ierr .ne. 0) go to 50
      call  htribk(nm,n,ar,ai,fm1,n,zr,zi)
   50 return
      end
c
c     ----------------------------------------------------------------
c
      subroutine htridi(nm,n,ar,ai,d,e,e2,tau)
c
      integer i,j,k,l,n,ii,nm,jp1
      real ar(nm,n),ai(nm,n),d(n),e(n),e2(n),tau(2,n)
      real f,g,h,fi,gi,hh,si,scale
c     real sqrt,cabs,abs
c     complex cmplx
c
c     this subroutine is a translation of a complex analogue of
c     the algol procedure tred1, num. math. 11, 181-195(1968)
c     by martin, reinsch, and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
c
c     this subroutine reduces a complex hermitian matrix
c     to a real symmetric tridiagonal matrix using
c     unitary similarity transformations.
c
c     on input-
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement,
c
c        n is the order of the matrix,
c
c        ar and ai contain the real and imaginary parts,
c          respectively, of the complex hermitian input matrix.
c          only the lower triangle of the matrix need be supplied.
c
c     on output-
c
c        ar and ai contain information about the unitary trans-
c          formations used in the reduction in their full lower
c          triangles.  their strict upper triangles and the
c          diagonal of ar are unaltered,
c
c        d contains the diagonal elements of the the tridiagonal matrix,
c
c        e contains the subdiagonal elements of the tridiagonal
c          matrix in its last n-1 positions.  e(1) is set to zero,
c
c        e2 contains the squares of the corresponding elements of e.
c          e2 may coincide with e if the squares are not needed,
c
c        tau contains further information about the transformations.
c
c     arithmetic is real except for the use of the subroutines
c     cabs and cmplx in computing complex absolute values.
c
c     questions and comments should be directed to b. s. garbow,
c     applied mathematics division, argonne national laboratory
c
c     ------------------------------------------------------------------
c
c
c
      tau(1,n) = 1.0
      tau(2,n) = 0.0
c
      do 100 i = 1, n
  100 d(i) = ar(i,i)
c     .......... for i=n step -1 until 1 do -- ..........
      do 300 ii = 1, n
         i = n + 1 - ii
         l = i - 1
         h = 0.0
         scale = 0.0
         if (l .lt. 1) go to 130
c     .......... scale row (algol tol then not needed) ..........
         do 120 k = 1, l
  120    scale = scale + abs(ar(i,k)) + abs(ai(i,k))
c
         if (scale .ne. 0.0) go to 140
         tau(1,l) = 1.0
         tau(2,l) = 0.0
  130    e(i) = 0.0
         e2(i) = 0.0
         go to 290
c
  140    do 150 k = 1, l
            ar(i,k) = ar(i,k) / scale
            ai(i,k) = ai(i,k) / scale
            h = h + ar(i,k) * ar(i,k) + ai(i,k) * ai(i,k)
  150    continue
c
         e2(i) = scale * scale * h
         g = sqrt(h)
         e(i) = scale * g
         f = cabs(cmplx(ar(i,l),ai(i,l)))
c     .......... form next diagonal element of matrix t ..........
         if (f .eq. 0.0) go to 160
         tau(1,l) = (ai(i,l) * tau(2,i) - ar(i,l) * tau(1,i)) / f
         si = (ar(i,l) * tau(2,i) + ai(i,l) * tau(1,i)) / f
         h = h + f * g
         g = 1.0 + g / f
         ar(i,l) = g * ar(i,l)
         ai(i,l) = g * ai(i,l)
         if (l .eq. 1) go to 270
         go to 170
  160    tau(1,l) = -tau(1,i)
         si = tau(2,i)
         ar(i,l) = g
  170    f = 0.0
c
         do 240 j = 1, l
            g = 0.0
            gi = 0.0
c     .......... form element of a*u ..........
            do 180 k = 1, j
               g = g + ar(j,k) * ar(i,k) + ai(j,k) * ai(i,k)
               gi = gi - ar(j,k) * ai(i,k) + ai(j,k) * ar(i,k)
  180       continue
c
            jp1 = j + 1
            if (l .lt. jp1) go to 220
c
            do 200 k = jp1, l
               g = g + ar(k,j) * ar(i,k) - ai(k,j) * ai(i,k)
               gi = gi - ar(k,j) * ai(i,k) - ai(k,j) * ar(i,k)
  200       continue
c     .......... form element of p ..........
  220       e(j) = g / h
            tau(2,j) = gi / h
            f = f + e(j) * ar(i,j) - tau(2,j) * ai(i,j)
  240    continue
c
         hh = f / (h + h)
c     .......... form reduced a ..........
         do 260 j = 1, l
            f = ar(i,j)
            g = e(j) - hh * f
            e(j) = g
            fi = -ai(i,j)
            gi = tau(2,j) - hh * fi
            tau(2,j) = -gi
c
            do 260 k = 1, j
               ar(j,k) = ar(j,k) - f * e(k) - g * ar(i,k)
     x                           + fi * tau(2,k) + gi * ai(i,k)
               ai(j,k) = ai(j,k) - f * tau(2,k) - g * ai(i,k)
     x                           - fi * e(k) - gi * ar(i,k)
  260    continue
c
  270    do 280 k = 1, l
            ar(i,k) = scale * ar(i,k)
            ai(i,k) = scale * ai(i,k)
  280    continue
c
         tau(2,l) = -si
  290    hh = d(i)
         d(i) = ar(i,i)
         ar(i,i) = hh
         ai(i,i) = scale * sqrt(h)
  300 continue
c
      return
      end
c
c     ------------------------------------------------------------------
c
      subroutine tqlrat(n,d,e2,ierr)
c
      integer i,j,l,m,n,ii,l1,mml,ierr
      real d(n),e2(n)
      real b,c,f,g,h,p,r,s,machep
c     real sqrt,abs,sign
c
c     this subroutine is a translation of the algol procedure tqlrat,
c     algorithm 464, comm. acm 16, 689(1973) by reinsch.
c
c     this subroutine finds the eigenvalues of a symmetric
c     tridiagonal matrix by the rational ql method.
c
c     on input-
c
c        n is the order of the matrix,
c
c        d contains the diagonal elements of the input matrix,
c
c        e2 contains the squares of the subdiagonal elements of the
c          input matrix in its last n-1 positions.  e2(1) is arbitrary.
c
c      on output-
c
c        d contains the eigenvalues in ascending order.  if an
c          error exit is made, the eigenvalues are correct and
c          ordered for indices 1,2,...ierr-1, but may not be
c          the smallest eigenvalues,
c
c        e2 has been destroyed,
c
c        ierr is set to
c          zero       for normal return,
c          j          if the j-th eigenvalue has not been
c                     determined after 30 iterations.
c
c     questions and comments should be directed to b. s. garbow,
c     applied mathematics division, argonne national laboratory
c
c     ------------------------------------------------------------------
c
c     .......... machep is a machine dependent parameter specifying
c                the relative precision of floating point arithmetic.
c                ..........
      machep = 2.**(-47)
c
      ierr = 0
      if (n .eq. 1) go to 1001
c
      do 100 i = 2, n
  100 e2(i-1) = e2(i)
c
      f = 0.0
      b = 0.0
      e2(n) = 0.0
c
      do 290 l = 1, n
         j = 0
         h = machep * (abs(d(l)) + sqrt(e2(l)))
         if (b .gt. h) go to 105
         b = h
         c = b * b
c     .......... look for small squared sub-diagonal element ..........
  105    do 110 m = l, n
            if (e2(m) .le. c) go to 120
c     .......... e2(n) is always zero, so there is no exit
c                through the bottom of the loop ..........
  110    continue
c
  120    if (m .eq. l) go to 210
  130    if (j .eq. 30) go to 1000
         j = j + 1
c     .......... form shift ..........
         l1 = l + 1
         s = sqrt(e2(l))
         g = d(l)
         p = (d(l1) - g) / (2.0 * s)
         r = sqrt(p*p+1.0)
         d(l) = s / (p + sign(r,p))
         h = g - d(l)
c
         do 140 i = l1, n
  140    d(i) = d(i) - h
c
         f = f + h
c     .......... rational ql transformation ..........
         g = d(m)
         if (g .eq. 0.0) g = b
         h = g
         s = 0.0
         mml = m - l
c     .......... for i=m-1 step -1 until l do -- ..........
         do 200 ii = 1, mml
            i = m - ii
            p = g * h
            r = p + e2(i)
            e2(i+1) = s * r
            s = e2(i) / r
            d(i+1) = h + s * (h + d(i))
            g = d(i) - e2(i) / g
            if (g .eq. 0.0) g = b
            h = g * p / r
  200    continue
c
c
         e2(l) = s * g
         d(l) = h
c     .......... guard against underflow in convergence test ..........
         if (h .eq. 0.0) go to 210
         if (abs(e2(l)) .le. abs(c/h)) go to 210
         e2(l) = h * e2(l)
         if (e2(l) .ne. 0.0) go to 130
  210    p = d(l) + f
c     ......... order eigenvalues ..........
         if (l .eq. 1) go to 250
c     .......... for i=l step -1 until 2 do -- ..........
         do 230 ii = 2, l
            i = l + 2 - ii
            if (p .ge. d(i-1)) go to 270
            d(i) = d(i-1)
  230    continue
c
  250    i = 1
  270    d(i) = p
  290 continue
c
      go to 1001
c     .......... set error -- no convergence to an
c                eigenvalue after 30 iterations **********
 1000 ierr = l
 1001 return
      end
c
c     ------------------------------------------------------------------
c
      subroutine tql2(nm,n,d,e,z,ierr)
c
      integer i,j,k,l,m,n,ii,l1,nm,mml,ierr
      real d(n),e(n),z(nm,n)
      real b,c,f,g,h,p,r,s,machep
c     real sqrt,abs,sign
c
c     this subroutine is a translation of the algol procedure tql2,
c     num. math. 11, 293-306(1968) by bowdler, martin, reinsch, and
c     wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 227-240(1971).
c
c     this subroutine finds the eigenvalues and eigenvectors
c     of a symmetric tridiagonal matrix by the ql method.
c     the eigenvectors of a full symmetric matrix can also
c     be found if  tred2  has been used to reduce this
c     full matrix to tridiagonal form.
c
c     on input-
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement,
c
c        n is the order of the matrix,
c
c        d contains the diagonal elements of the input matrix,
c
c        e contains the subdiagonal elements of the input matrix
c          in its last n-1 positions.  e(1) is arbitrary,
c
c        z contains the transformation matrix produced in the
c          reduction by  tred2, if performed.  if the eigenvectors
c          of the tridiagonal matrix are desired, z must contain
c          the identity matrix.
c
c      on output-
c
c        d contains the eigenvalues in ascending order.  if an
c          error exit is made, the eigenvalues are correct but
c          unordered for indices 1,2,...,ierr-1,
c
c        e has been destroyed,
c
c        z contains orthonormal eigenvectors of the symmetric
c          tridiagonal (or full) matrix.  if an error exit is made,
c          z contains the eigenvectors associated with the stored
c          eigenvalues,
c
c        ierr is set to
c          zero       for normal return,
c          j          if the j-th eigenvalue has not been
c                     determined after 30 iterations.
c
c     questions and comments should be directed to b. s. garbow,
c     applied mathematics division, argonne national laboratory
c
c     ------------------------------------------------------------------
c
c     .......... machep is a machine dependent parameter specifying
c                the relative precision of floating point arithmetic.
c                ..........
      machep = 2.**(-47)
c
      ierr = 0
      if (n .eq. 1) go to 1001
c
      do 100 i = 2, n
  100 e(i-1) = e(i)
c
      f = 0.0
      b = 0.0
      e(n) = 0.0
c
      do 240 l = 1, n
         j = 0
         h = machep * (abs(d(l)) + abs(e(l)))
         if (b .lt. h) b = h
c     .......... look for small sub-diagonal element ..........
         do 110 m = l, n
            if (abs(e(m)) .le. b) go to 120
c     .......... e(n) is always zero, so there is no exit
c                through the bottom of the loop ..........
  110    continue
c
  120    if (m .eq. l) go to 220
  130    if (j .eq. 30) go to 1000
         j = j + 1
c     ........ form shift ..........
         l1 = l + 1
         g = d(l)
         p = (d(l1) - g) / (2.0 * e(l))
         r = sqrt(p*p+1.0)
         d(l) = e(l) / (p + sign(r,p))
         h = g - d(l)
c
         do 140 i = l1, n
  140    d(i) = d(i) - h
c
         f = f + h
c     .......... ql transformation ..........
         p = d(m)
         c = 1.0
         s = 0.0
         mml = m - l
c     .......... for i=m-1 step -1 until l do -- ..........
         do 200 ii = 1, mml
            i = m - ii
            g = c * e(i)
            h = c * p
            if (abs(p) .lt. abs(e(i))) go to 150
            c = e(i) / p
            r = sqrt(c*c+1.0)
            e(i+1) = s * p * r
            s = c / r
            c = 1.0 / r
            go to 160
  150       c = p / e(i)
            r = sqrt(c*c+1.0)
            e(i+1) = s * e(i) * r
            s = 1.0 / r
            c = c * s
  160       p = c * d(i) - s * g
            d(i+1) = h + s * (c * g + s * d(i))
c     .......... form vector ..........
            do 180 k = 1, n
               h = z(k,i+1)
               z(k,i+1) = s * z(k,i) + c * h
               z(k,i) = c * z(k,i) - s * h
  180       continue
c
  200    continue
c
         e(l) = s * p
         d(l) = c * p
         if (abs(e(l)) .gt. b) go to 130
  220    d(l) = d(l) + f
  240 continue
c     .......... order eigenvalues and eigenvectors ..........
      do 300 ii = 2, n
         i = ii - 1
         k = i
         p = d(i)
c
         do 260 j = ii, n
            if (d(j) .ge. p) go to 260
            k = j
            p = d(j)
  260    continue
c
         if (k .eq. i) go to 300
         d(k) = d(i)
         d(i) = p
c
         do 280 j = 1, n
            p = z(j,i)
            z(j,i) = z(j,k)
            z(j,k) = p
  280    continue
c
  300 continue
c
      go to 1001
c     .......... set error -- no convergence to an
c                eigenvalue after 30 iterations ..........
 1000 ierr = l
 1001 return
      end
c
c     ------------------------------------------------------------------
c
      subroutine htribk(nm,n,ar,ai,tau,m,zr,zi)
c
      integer i,j,k,l,m,n,nm
      real ar(nm,n),ai(nm,n),tau(2,n),zr(nm,m),zi(nm,m)
      real h,s,si
c
c     this subroutine is a translation of a complex analogue of
c     the algol procedure trbak1, num. math. 11, 181-195(1968)
c     by martin, reinsch, and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
c
c     this subroutine forms the eigenvectors of a complex hermitian
c     matrix by back transforming those of the corresponding
c     real symmetric tridiagonal matrix determined by  htridi.
c
c     on input-
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement,
c
c        n is the order of the matrix,
c
c        ar and ai contain information about the unitary trans-
c          formations used in the reduction by  htridi  in their
c          full lower triangles except for the diagonal of ar,
c
c        tau contains further information about the transformations,
c
c        m is the number of eigenvectors to be back transformed,
c
c        zr contains the eigenvectors to be back transformed
c          in its first m columns.
c
c     on output-
c
c        zr and zi contain the real and imaginary parts,
c          respectively, of the transformed eigenvectors
c          in their first m columns.
c
c     note that the last component of each returned vector
c     is real and that vector euclidean norms are preserved.
c
c     questions and comments should be directed to b. s. garbow,
c     applied mathematics division, argonne national laboratory
c
c     ------------------------------------------------------------------
c
      if (m .eq. 0) go to 200
c     .......... transform the eigenvectors of the real symmetric
c                tridiagonal matrix to those of the hermitian
c                tridiagonal matrix. ..........
      do 50 k = 1, n
c
         do 50 j = 1, m
            zi(k,j) = -zr(k,j) * tau(2,k)
            zr(k,j) = zr(k,j) * tau(1,k)
   50 continue
c
      if (n .eq. 1) go to 200
c     .......... recover and apply the householder matrices ..........
      do 140 i = 2, n
         l = i - 1
         h = ai(i,i)
         if (h .eq. 0.0) go to 140
c
         do 130 j = 1, m
            s = 0.0
            si = 0.0
c
            do 110 k = 1, l
               s = s + ar(i,k) * zr(k,j) - ai(i,k) * zi(k,j)
               si = si + ar(i,k) * zi(k,j) + ai(i,k) * zr(k,j)
  110       continue
c     .......... double divisions avoid possible underflow ..........
            s = (s / h) / h
            si = (si / h) / h
c
            do 120 k = 1, l
               zr(k,j) = zr(k,j) - s * ar(i,k) - si * ai(i,k)
               zi(k,j) = zi(k,j) - si * ar(i,k) + s * ai(i,k)
  120       continue
c
  130    continue
c
  140 continue
c
  200 return
      end

    

