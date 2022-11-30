#include "fintrf.h"
! ------------------------------------------------------------------
! sglfitF.f90: block coordinate descent for sg-LASSO MSE regression.
! ------------------------------------------------------------------
      SUBROUTINE sglfit(gamma, ngroups, gindex, nobs,
     & nvars, x, y, pf, dfmax, pmax, nlam, flmin, ulam, eps, peps,
     & isd, intr, maxit, nalam, b0, beta, ibeta, nbeta, alam,
     & npass, jerr)
      
      IMPLICIT NONE
      ! -------- INPUT VARIABLES -------- !
      Integer nobs, nvars, dfmax, pmax, nlam, isd, intr, ngroups, nalam
      Integer maxit, npass, jerr
      !Integer nbeta(nlam), ibeta(pmax)
	  !Integer gindex(ngroups)
      Integer*4  nbeta(nlam), ibeta(pmax)
	  Integer*4  gindex(ngroups)

      Real*8 flmin, eps, peps, gamma
      Real*8 x(nobs, nvars), y(nobs), pf(nvars)
      Real*8 b0(nlam), beta(pmax, nlam)
      Real*8 alam(nlam), ulam(nlam) 


      INTEGER j, l, nk, ierr
      INTEGER, Dimension(:), Allocatable :: ju
      Real*8, Dimension(:), Allocatable :: xmean
	  Real*8, Dimension(:), Allocatable :: xnorm
	  Real*8, Dimension(:), Allocatable :: maj
      
	  nalam=0
      b0=0.d0
      beta=0.d0
      ibeta=0
      nbeta=0
      alam=0.d0
      npass=0
      jerr=0
	  
      ALLOCATE(ju(1:nvars), Stat=ierr)
      jerr = jerr + ierr
      ALLOCATE(xmean(1:nvars), Stat=ierr)
      jerr = jerr + ierr
      ALLOCATE(xnorm(1:nvars), Stat=ierr)
      jerr = jerr + ierr
      ALLOCATE(maj(1:nvars), Stat=ierr)
      jerr = jerr + ierr
      If (jerr /= 0) Return
      Call chkvars(nobs, nvars, x, ju)
      If (maxval(pf) <= 0.0D0) Then
          jerr = 10000
          Return
      End If
      pf = max(0.0D0, pf)
      

      
      Call standard(nobs, nvars, x, ju, isd, intr, xmean, xnorm, maj)
      If (gamma == 1.0D0) Then
          Call lassofitpathF(maj, nobs, nvars, x, y, ju, pf, dfmax,
     & pmax, nlam, flmin, ulam, eps, maxit, nalam, b0, beta, ibeta,
     & nbeta, alam, npass, jerr, intr)
      Else
          Call sglfitpathF(maj, gamma, ngroups, gindex, nobs, nvars,
     & x, y, ju, pf, dfmax, pmax, nlam, flmin, ulam, eps, peps, maxit, 
     & nalam, b0, beta, ibeta, nbeta, alam, npass, jerr, intr)
      End If
      If (jerr > 0) Return
      
      Do l = 1, nalam
          nk = nbeta(l)
          If (isd == 1) Then
              Do j = 1, nk
                  beta(j,l) = beta(j,l)/xnorm(ibeta(j))
              End Do
          End If
          If (intr == 1) Then
              b0(l)=b0(l)-dot_product(beta(1:nk,l),xmean(ibeta(1:nk)))
          End If
      End Do
	  
      DEALLOCATE(ju,xmean,xnorm,maj)
      Return
      End SUBROUTINE sglfit
	  
	  
	  
      SUBROUTINE sglfitpathF(maj, gamma, ngroups, gindex,
     & nobs, nvars, x, y, ju, pf, dfmax, pmax, nlam, flmin,
     & ulam, eps, peps, maxit, nalam, b0, beta, m, nbeta,
     & alam, npass, jerr, intr)
      IMPLICIT NONE
      INTEGER mnl, nobs, nvars, dfmax, pmax, nlam, maxit
      INTEGER nalam, npass, jerr, intr, ngroups
      !INTEGER ju(nvars), m(pmax), nbeta(nlam), gindex(ngroups)
      INTEGER ju(nvars) 
      INTEGER*4 m(pmax), nbeta(nlam), gindex(ngroups)
      Real*8 eps, gamma, bnorm, peps
      Real*8 x(nobs, nvars), y(nobs), maj(nvars)
      Real*8 pf(nvars)
      Real*8 beta(pmax, nlam), b0(nlam)
      Real*8 ulam(nlam), alam(nlam)
      
      INTEGER, Parameter :: mnlam = 6
      Real*8 d, dif, oldb, u, al, flmin
      Real*8, Dimension(:), Allocatable :: b, oldbeta, r
      Real*8 gw
      INTEGER gstart, gend
      INTEGER k, j, l, g, vrg, ctr, ierr, ni, me, pln
      INTEGER gs, gj, skip
      INTEGER, Dimension(:), Allocatable :: mm
	  
	  
      ALLOCATE(b(0:nvars), Stat=jerr)
      ALLOCATE(oldbeta(0:nvars), Stat=ierr)
      jerr = jerr + ierr
      ALLOCATE(mm(1:nvars), Stat=ierr)
      jerr = jerr + ierr
      ALLOCATE(r(1:nobs), Stat=ierr)
      jerr = jerr + ierr
      If (jerr /= 0) Return
      b = 0.0D0
      oldbeta = 0.0D0
      m = 0
      mm = 0
      npass = 0
      ni = npass
      mnl = min(mnlam,nlam)
      ju = 0
      r = y
      Do l = 1, nlam
        al = ulam(l)
        ctr = 0
        pln = 0
        Do
          If (intr == 1) oldbeta(0) = b(0)
          If (ni > 0) oldbeta(m(1:ni)) = b(m(1:ni))
          Do
            npass = npass + 1
            dIf = 0.0D0
            If (intr == 1) oldbeta(0) = b(0)
            If (ni > 0) oldbeta(m(1:ni)) = b(m(1:ni))
            bnorm = sum(b**2)
            bnorm = sqrt(bnorm)
            pln = pln + 1
            Do k = 1, ngroups
              gend = gindex(k)
              If (k == 1) Then
                gstart = 1
              Else
                gstart = gindex(k-1) + 1
              End If
              gs = gend - gstart + 1
              gw = 0.0D0
              Do gj = gstart, gend
                gw = gw + pf(gj)
              End Do
              gw = sqrt(gw)
              skip = 1
              If (pln == 1) Then
                skip = 0
              End If
              If (ju(k) == 1) Then
                skip = 0
              End If
              If (skip == 0) Then
                    Call prox_sgl(gstart, gend, nvars, nobs, x,
     & r, b(1:nvars), al, gamma, pf, peps, gw)
                Do g = gstart, gend
                  If (abs(b(g))>0.0D0) Then
                    If (pln == 1) Then
                      ju(k) = 1
                    End If
                    d = oldbeta(g) - b(g)
                    dif = max(dif, maj(g) * d**2)
                    If (mm(g) == 0) Then
                      ni = ni + 1
                      If (ni > pmax) Exit
                      mm(g) = ni
                      m(ni) = g
                    End If
                  End If
                End Do
              End If
            End Do
            If (ni > pmax) Exit
            If (intr == 1) Then
              oldb = oldbeta(0)
              u = 0.0D0
              DO ! BEGIN GRADIENT DESCENT
                d = SUM(r)/nobs
                IF (d**2 < eps) EXIT
                b(0) = b(0) + d
                r = r - d
              END DO ! END GRADIENT DESCENT
              !b(0) = u
              d = b(0) - oldb
              If (abs(d) > 0.0D0) dif = max(dif, d**2)
            End If
            If (dif < eps) Exit
          End Do
          If (ni > pmax) Exit
          vrg = 1
          If (intr == 1) Then
              If ((b(0) - oldbeta(0))**2 >= eps) vrg = 0
          End If
          Do j = 1, ni
              If ((b(m(j)) - oldbeta(m(j)))**2 >= eps) Then
                  vrg = 0
                  Exit
              End If
          End Do
          If (vrg == 1) Exit
          ctr = ctr + 1
          If (ctr > maxit) Then
              jerr = - l
              Return
          End If
        End Do
        If (ni > pmax) Then
          jerr = - 10000 - l
          Exit
        End If
        If (ni > 0) beta(1:ni,l) = b(m(1:ni))
        nbeta(l) = ni
        b0(l) = b(0)
        alam(l) = al
        nalam = l
        If (l < mnl) Cycle
        If (flmin >= 1.0D0) Cycle
        me = count(abs(beta(1:ni,l)) > 0.0D0)
        If (me > dfmax) Exit
      End Do
      DEALLOCATE(b,oldbeta,mm,r)
	  
      RETURN
      End SUBROUTINE sglfitpathF



      SUBROUTINE lassofitpathF(maj, nobs, nvars, x, y, ju, pf,
     & dfmax, pmax, nlam, flmin, ulam, eps, maxit, nalam, b0, 
     & beta, m, nbeta, alam, npass, jerr, intr)
      
      IMPLICIT NONE
      ! -------- INPUT VARIABLES -------- !
      INTEGER mnl, nobs, nvars, dfmax, pmax, nlam, maxit
      INTEGER nalam, npass, jerr, intr
      !INTEGER ju(nvars), m(pmax), nbeta(nlam)
      INTEGER ju(nvars)
      INTEGER*4 m(pmax), nbeta(nlam)
      Real*8 eps
      Real*8 x(nobs, nvars), y(nobs), maj(nvars)
      Real*8 pf(nvars)
      Real*8 beta(pmax, nlam), b0(nlam)
      Real*8 ulam(nlam), alam(nlam)
      ! -------- LOCAL DECLARATIONS -------- !
      INTEGER, PARAMETER :: mnlam = 6
      Real*8 tmp, d, dif, oldb, u, v, al, flmin
      Real*8, DIMENSION(:), ALLOCATABLE :: b, oldbeta, r
      
      INTEGER k, j, l, vrg, ctr, ierr, ni, me, pln
      INTEGER, DIMENSION(:), ALLOCATABLE :: mm
      ! -------- ALLOCATE VARIABLES -------- !
      ALLOCATE(b(0:nvars), STAT=jerr)
      ALLOCATE(oldbeta(0:nvars), STAT=ierr)
      jerr = jerr + ierr
      ALLOCATE(mm(1:nvars), STAT=ierr)
      jerr = jerr + ierr
      ALLOCATE(r(1:nobs), STAT=ierr)
      jerr = jerr + ierr
      IF (jerr /= 0) RETURN
      ! ---------------- INITIALIZATION ---------------- !
      r = y
      b = 0.0D0
      oldbeta = 0.0D0
      m = 0
      mm = 0
      npass = 0
      ni = npass
      mnl = MIN(mnlam,nlam)
      ju = 0
      
      DO l = 1, nlam
        al = ulam(l)
        ctr = 0
        ! ------------------ OUTER LOOP -------------------- !
        DO
          IF (intr == 1) oldbeta(0) = b(0)
          IF (ni > 0) oldbeta(m(1:ni)) = b(m(1:ni))
          pln = 0
          ! MIDDLE LOOP - !
          DO
            npass = npass + 1
            dif = 0.0D0
            pln = pln + 1
            DO k = 1, nvars
              IF (pln == 1) THEN
                oldb = b(k)
                DO ! BEGIN PROXIMAL COORDINATE DESCENT
                  !u = 0.0D0
                  u = maj(k) * b(k) + DOT_PRODUCT(r,x(:,k))/nobs
                  v = ABS(u) - al * pf(k)
                  IF (v > 0.0D0) THEN
                    tmp = SIGN(v,u)/maj(k)
                  ELSE
                    tmp = 0.0D0
                  END IF
                  d = tmp - b(k)
                  IF (d**2 < eps) EXIT
                  b(k) = tmp
                  r = r - x(:,k) * d
                END DO ! END PROXIMAL GRADIENT DESCENT
                d = b(k) - oldb
                IF (ABS(d) > 0.0D0) THEN
                  dif = dif + maj(k) * d**2
                  IF (mm(k) == 0) THEN
                    ni = ni + 1
                    IF (ni > pmax) EXIT
                    mm(k) = ni
                    m(ni) = k ! RECORD ACTIVE VARIABLES
                  END IF
                END IF
                IF (ABS(b(k))>0.0D0) ju(k) = 1
              ELSE
                IF (ju(k) == 1) THEN
                  oldb = b(k)
                  DO ! BEGIN PROXIMAL GRADIENT DESCENT
                    !u = 0.0D0
                    u = DOT_PRODUCT(r,x(:,k))
                    u = maj(k) * b(k) + u/nobs
                    v = ABS(u) - al * pf(k)
                    IF (v > 0.0D0) THEN
                      tmp = SIGN(v,u)/maj(k)
                    ELSE
                      tmp = 0.0D0
                    END IF
                    d = tmp - b(k)
                    IF (d**2 < eps) EXIT
                    b(k) = tmp
                    r = r - x(:,k) * d
                  END DO ! END PROXIMAL GRADIENT DESCENT
                  d = b(k) - oldb
                  IF (ABS(d) > 0.0D0) THEN
                    dif = MAX(dif, maj(k) * d**2)
                  END IF
                END IF
              END IF
            END DO
            IF (ni > pmax) EXIT
            IF (intr == 1) THEN
              oldb = b(0)
              DO ! BEGIN GRADIENT DESCENT
                d = SUM(r)/nobs
                IF (d**2 < eps) EXIT
                b(0) = b(0) + d
                r = r - d
              END DO ! END GRADIENT DESCENT
              d = b(0) - oldb
              IF (ABS(d) > 0.0D0) dif = MAX(dif, d**2)
            END IF
            IF (dif < eps) EXIT
          END DO ! ----------> END MIDDLE LOOP
          IF (ni > pmax) EXIT
          ! -------------- FINAL CHECK ---------------- !
          vrg = 1
          IF (intr == 1) THEN
            IF ((b(0) - oldbeta(0))**2 >= eps) vrg = 0
          END IF
          DO j = 1, ni
            IF ((b(m(j)) - oldbeta(m(j)))**2 >= eps) THEN
              vrg = 0
              EXIT
            END IF
          END DO
          IF (vrg == 1) EXIT
          ctr = ctr + 1
          IF (ctr > maxit) THEN
            jerr = - l
            RETURN
          END IF
        END DO ! -------> END OUTER LOOP
        ! ----------- FINAL UPDATE & SAVE RESULTS ------------ !
        IF (ni > pmax) THEN
          jerr = - 10000 - l
          EXIT
        END IF
        IF (ni > 0) beta(1:ni,l) = b(m(1:ni))
        nbeta(l) = ni
        b0(l) = b(0)
        alam(l) = al
        nalam = l
        IF (l < mnl) CYCLE
        IF (flmin >= 1.0D0) CYCLE
        me = COUNT(ABS(beta(1:ni,l)) > 0.0D0)
        IF (me > dfmax) EXIT
      END DO ! -------> END LAMBDA LOOP
      DEALLOCATE(b,oldbeta,r,mm)
      RETURN
      END SUBROUTINE lassofitpathF



      SUBROUTINE prox_sgl(gstart, gend, nvars, nobs, x, r,
     & b, al, gamma, pf, peps, gw)
      IMPLICIT NONE
      ! -------- INPUT VARIABLES -------- !
      INTEGER gstart, gend, nvars, nobs
      Real*8 x(nobs,nvars), r(nobs), b(nvars), al
      Real*8 gamma, pf(nvars), peps, gw
      ! -------- LOCAL DECLARATIONS -------- !
      INTEGER :: g
      Real*8 u, scl, tmp, maxg, normg, d 
      Real*8 bold(nvars), vg
      Real*8, Parameter :: big = 9.9D30
      ! -------- BEGIN PROGRAM -------- !
	  Do
       bold(gstart:gend) = b(gstart:gend)
        !--------- LASSO PART ----------!
        Do g = gstart, gend
          u = b(g) + DOT_PRODUCT(x(:,g),r)/nobs
          !S(.) map
          tmp=max(u-al*gamma*pf(g),0.0D0)-max(-u-al*gamma*pf(g),0.0D0)
          b(g) = tmp
        End Do
         !--------- g-LASSO PART ----------!
         ! L2 norm of b_g
         normg = NORM2(b(gstart:gend))
         ! Initialize storage vars
         maxg = 0.0D0
         vg = gw * al * (1.0D0-gamma)/normg
         If (normg .EQ. 0.0D0) Then
           vg = big
         End If
         Do g = gstart, gend
           scl = 1.0D0 - pf(g) * vg
           scl = MAX(scl, 0.0D0)
           ! l_2,1 norm map
           tmp = scl*b(g)
           d = tmp - bold(g)
           r = r - x(:,g)*d
           maxg = MAX(maxg, ABS(d))
           b(g) = tmp
         End Do
      !--------- CHECK CONVERGENCE ----------!
        If (maxg < peps) Exit
	  !Exit
      End Do
      End SUBROUTINE prox_sgl

      SUBROUTINE standard(nobs, nvars, x, ju, isd, intr,
     & xmean, xnorm, maj)     
      IMPLICIT NONE
      ! -------- INPUT VARIABLES -------- !
      INTEGER nobs, nvars, isd, intr, ju(nvars)
      Real*8 x(nobs, nvars), xmean(nvars)
      Real*8xnorm(nvars), maj(nvars)
      ! -------- LOCAL DECLARATIONS -------- !
      INTEGER j
      Real*8 xmsq, xvar
      ! -------- BEGIN PROGRAM -------- !
      IF (intr == 0) THEN
      DO j = 1, nvars
      IF (ju(j) == 1) THEN
      xmean(j) = 0.0D0
      maj(j) = DOT_PRODUCT(x(:,j),x(:,j))/nobs
      IF (isd == 1) THEN
      xmsq = (SUM(x(:,j))/nobs)**2
      xvar = maj(j) - xmsq
      xnorm(j) = SQRT(xvar)
      x(:,j) = x(:,j)/xnorm(j)
      maj(j) = 1.0D0 + xmsq/xvar
      END IF
      END IF
      END DO
      ELSE                   
      DO j = 1, nvars                                  
      IF (ju(j) == 1) THEN                         
      xmean(j) = SUM(x(:,j))/nobs  ! MEAN                        
      x(:,j) = x(:,j) - xmean(j)    
      maj(j) = DOT_PRODUCT(x(:,j),x(:,j))/nobs                                             
      IF (isd == 1) THEN
      xnorm(j) = SQRT(maj(j))  ! STANDARD DEVIATION              
      x(:,j) = x(:,j)/xnorm(j)
      maj(j) = 1.0D0
      END IF                                                        
      END IF                                     
      END DO  
      END IF                           
      END SUBROUTINE standard
      
      ! ---------------------------------------------------------- !
      SUBROUTINE chkvars(nobs, nvars, x, ju)
      IMPLICIT NONE
      ! -------- INPUT VARIABLES -------- !
      INTEGER :: nobs,nvars,ju(nvars)
      Real*8 :: x(nobs,nvars)
      ! -------- LOCAL DECLARATIONS -------- !
      INTEGER :: i,j
      Real*8 :: t
      ! -------- BEGIN PROGRAM -------- ! 
      DO j = 1, nvars
      ju(j) = 0
      t = x(1,j)
      DO i = 2, nobs
      IF (x(i,j) /= t) THEN
      ju(j) = 1
      EXIT
      END IF
      END DO
      END DO
      END SUBROUTINE chkvars



      SUBROUTINE mexFunction(nlhs, plhs, nrhs, prhs)
        IMPLICIT NONE
      mwPointer  mxGetM, mxGetN, mxGetPr
      mwPointer  mxIsNumeric, mxCreateDoubleMatrix
      mwPointer  mxClassIDFromClassName, mxCreateNumericArray
      mwPointer  plhs(*), prhs(*)
      mwPointer  temp_pr
      mwPointer  gamma_pr, ngroups_pr, gindex_pr, nobs_pr
      mwPointer  nvars_pr, x_pr, y_pr, pf_pr, dfmax_pr
      mwPointer  pmax_pr, nlam_pr, flmin_pr, ulam_pr
      mwPointer  eps_pr, peps_pr, isd_pr, intr_pr, maxit_pr
      mwPointer  nalam_pr, b0_pr, beta_pr, ibeta_pr
      mwPointer  nbeta_pr, alam_pr, npass_pr, jerr_pr
      INTEGER nlhs, nrhs
      INTEGER m, n, size, tmp1, tmp2, nlams, q
      INTEGER dims(2), clsid

      Real*8 gamma,flmin, eps, peps
      Integer dfmax, ngroups, nobs, nvars
      Integer  isd
      Integer intr, maxit, nalam
      Integer npass, jerr
      Integer*4 gindex(100000)
      Integer pmax
      Integer ingroups, nlam


      
      Real*8, Dimension(:,:), Allocatable :: x
	  Real*8, Dimension(:), Allocatable :: y
	  Real*8, Dimension(:), Allocatable :: pf
	  Real*8, Dimension(:), Allocatable :: ulam
	  Real*8, Dimension(:), Allocatable :: b0
	  Real*8, Dimension(:,:), Allocatable :: beta
      Integer*4, Dimension(:), Allocatable :: ibeta
	  Integer*4, Dimension(:), Allocatable :: nbeta
	  Real*8, Dimension(:), Allocatable :: alam


      gamma_pr = mxGetPr(prhs(1))
	  call mxCopyPtrToReal8(gamma_pr, gamma, 1)
      
	  ngroups_pr = mxGetPr(prhs(2))
      call mxCopyPtrToInteger8(ngroups_pr, ngroups, 1)
	  
	  gindex_pr = mxGetPr(prhs(3))
	  call mxCopyPtrToInteger4(gindex_pr, gindex, ngroups)
	  
	  nobs_pr = mxGetPr(prhs(4))
      call mxCopyPtrToInteger8(nobs_pr, nobs, 1)
	  
	  nvars_pr = mxGetPr(prhs(5))
      call mxCopyPtrToInteger8(nvars_pr, nvars, 1)
	  
	  m = mxGetM(prhs(6))
          n = mxGetN(prhs(6))
	  
          ALLOCATE(x(m, n))
	  x_pr = mxGetPr(prhs(6))
	  call mxCopyPtrToReal8(x_pr, x, m*n)
	  
	  ALLOCATE(y(m))
	  y_pr = mxGetPr(prhs(7))
	  call mxCopyPtrToReal8(y_pr, y, m)
	  
	  ALLOCATE(pf(n))
	  pf_pr = mxGetPr(prhs(8))
	  call mxCopyPtrToReal8(pf_pr, pf, n)
	  
	  dfmax_pr = mxGetPr(prhs(9))
          call mxCopyPtrToInteger8(dfmax_pr, dfmax, 1)
	  
	  pmax_pr = mxGetPr(prhs(10))
          call mxCopyPtrToInteger8(pmax_pr, pmax, 1)
	  
	  nlam_pr = mxGetPr(prhs(11))
           call mxCopyPtrToInteger8(nlam_pr, nlam, 1)

           
	  flmin_pr = mxGetPr(prhs(12))
	  call mxCopyPtrToReal8(flmin_pr, flmin, 1)
	  
	  tmp1 = mxGetM(prhs(13))
          tmp2 = mxGetN(prhs(13))
          If (tmp1 > tmp2) Then
             q = tmp1
          End If
          If (tmp2 > tmp1) Then
             q = tmp2
          End If
          If (tmp1 == tmp2) Then
             q = 1
          End If
	  
	  ALLOCATE(ulam(q))
	  ulam_pr = mxGetPr(prhs(13))
	  call mxCopyPtrToReal8(ulam_pr, ulam, q)
	  
	  eps_pr = mxGetPr(prhs(14))
	  call mxCopyPtrToReal8(eps_pr, eps, 1)
	  
	  peps_pr = mxGetPr(prhs(15))
	  call mxCopyPtrToReal8(peps_pr, peps, 1)
	  
	  isd_pr = mxGetPr(prhs(16))
          call mxCopyPtrToInteger8(isd_pr, isd, 1)
	  
	  intr_pr = mxGetPr(prhs(17))
          call mxCopyPtrToInteger8(intr_pr, intr, 1)
	  
	  maxit_pr = mxGetPr(prhs(18))
	  call mxCopyPtrToInteger8(maxit_pr, maxit, 1)

   

	  clsid = mxClassIDFromClassName('int32')
          plhs(1) = mxCreateNumericArray(1, 1, clsid, 0)
  
          nalam_pr = mxGetPr(plhs(1))
	  
	  plhs(2) = mxCreateDoubleMatrix(nlam, 1, 0)
	  tmp1 = mxGetM(plhs(2))
	  ALLOCATE(b0(tmp1))
	  b0_pr = mxGetPr(plhs(2))
	  
	  plhs(3) = mxCreateDoubleMatrix(pmax, nlam, 0)
	  m = mxGetM(plhs(3))
          n = mxGetN(plhs(3))
	  ALLOCATE(beta(m, n))
	  beta_pr = mxGetPr(plhs(3))
	  
	  plhs(4) = mxCreateDoubleMatrix(pmax, 1, 0)
	  tmp1 = mxGetM(plhs(4))
	  dims(1) = tmp1
	  dims(2) = 1
	  clsid = mxClassIDFromClassName('int32')
          plhs(4) = mxCreateNumericArray(2, dims, clsid, 0)
	  
	  ALLOCATE(ibeta(tmp1))
	  ibeta_pr = mxGetPr(plhs(4))
	  
	  plhs(5) = mxCreateDoubleMatrix(nlam, 1, 0)
	  tmp1 = mxGetM(plhs(5))
	  dims(1) = tmp1
	  dims(2) = 1
	  clsid = mxClassIDFromClassName('int32')
          plhs(5) = mxCreateNumericArray(2, dims, clsid, 0)
	  

	  tmp1 = mxGetM(plhs(5))
	  ALLOCATE(nbeta(tmp1))
	  nbeta_pr = mxGetPr(plhs(5))
	  
	  plhs(6) = mxCreateDoubleMatrix(nlam, 1, 0)
	  tmp1 = mxGetM(plhs(6))
	  ALLOCATE(alam(tmp1))
	  alam_pr = mxGetPr(plhs(6))
	  
	  clsid = mxClassIDFromClassName('int32')
          plhs(7) = mxCreateNumericArray(1, 1, clsid, 0)
	  
	  npass_pr = mxGetPr(plhs(7))
	  
	  clsid = mxClassIDFromClassName('int32')
          plhs(8) = mxCreateNumericArray(1, 1, clsid, 0)
	  
	  jerr_pr = mxGetPr(plhs(8))
      
        
      Call sglfit(gamma, ngroups, gindex, nobs, nvars, 
     & x, y, pf, dfmax, pmax, nlam, flmin, ulam, eps, peps, isd,
     & intr, maxit, nalam, b0, beta, ibeta, nbeta, alam, 
     & npass, jerr)
	  
      call mxCopyInteger8ToPtr(nalam, nalam_pr, 1)     
      
	  tmp1 = mxGetM(plhs(2))
	  call mxCopyReal8ToPtr(b0, b0_pr, tmp1)
	  
	  m = mxGetM(plhs(3))
          n = mxGetN(plhs(3))
	  call mxCopyReal8ToPtr(beta, beta_pr, m*n)     


   
	  tmp1 = mxGetM(plhs(4))
          call mxCopyInteger4ToPtr(ibeta, ibeta_pr, tmp1)

          tmp1 = mxGetM(plhs(5))




      call mxCopyInteger4ToPtr(nbeta, nbeta_pr, tmp1)
  

   
	  tmp1 = mxGetM(plhs(6))
	  call mxCopyReal8ToPtr(alam, alam_pr, tmp1)
	  
	  call mxCopyInteger8ToPtr(npass, npass_pr, 1)     
	  call mxCopyInteger8ToPtr(jerr, jerr_pr, 1)        
	  
      
      Return
      End