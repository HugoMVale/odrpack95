subroutine dodmn &
   (head, fstitr, prtpen, &
    fcn, n, m, np, nq, job, beta, y, ldy, x, ldx, &
    we, we1, ldwe, ld2we, wd, ldwd, ld2wd, &
    ifixb, ifixx, ldifx, &
    betac, betan, betas, s, delta, deltan, deltas, &
    lower, upper, &
    t, f, fn, fs, fjacb, msgb, fjacd, msgd, &
    ssf, ss, tt, ldtt, stpb, stpd, ldstpd, &
    xplusd, wrk, lwrk, work, lwork, iwork, liwork, info, &
    bound)
!! Iteratively compute least squares solution.
! Routines Called FCN, DACCES, DCOPY, DDOT, DEVJAC, DFLAGS, DNRM2, DODLM,
!                 DODPCR, DODVCV, DUNPAC, DWGHT
! Date Written   860529   (YYMMDD)
! Revision Date  920619   (YYMMDD)

   use odrpack_kinds, only: wp, zero, one
   use odrpack, only: tempret

   logical, intent(in) :: head
      !! The variable designating whether the heading is to be printed (`head = .true.`)
      !! or not (`head = .false.`).
   logical, intent(inout) :: fstitr
      !! The variable designating whether this is the first iteration (`fstitr = .true.`)
      !! or not (`fstitr = .false.`).
   logical, intent(inout) :: prtpen
      !! The value designating whether the penalty parameter is to be printed in the
      !! iteration report (`prtpen = .true.`) or not (`prtpen = .false.`).
   procedure() :: fcn
      !! The user supplied subroutine for evaluating the model.
   integer, intent(in) :: n
      !! The number of observations.
   integer, intent(in) :: m
      !! The number of columns of data in the explanatory variable.
   integer, intent(in) :: np
      !! The number of function parameters.
   integer, intent(in) :: nq
      !! The number of responses per observation.
   integer, intent(in) :: job
      !! The variable controlling problem initialization and computational method.
   real(kind=wp), intent(inout) :: beta(np)
      !! The function parameters.
   real(kind=wp), intent(in) :: y(ldy, nq)
      !! The dependent variable. Unused when the model is implicit.
   integer, intent(in) :: ldy
      !! The leading dimension of array `y`.
   real(kind=wp), intent(in) :: x(ldx, m)
      !! The explanatory variable.
   integer, intent(in) :: ldx
      !! The leading dimension of array `x`.
   real(kind=wp), intent(in) :: we(ldwe, ld2we, nq)
      !! The `epsilon` weights.
   real(kind=wp), intent(in) :: we1(ldwe, ld2we, nq)
      !! The square root of the `epsilon` weights.
   integer, intent(in) :: ldwe
      !! The leading dimension of arrays `we` and `we1`.
   integer, intent(in) :: ld2we
      !! The second dimension of arrays `we` and `we1`.
   real(kind=wp), intent(in) :: wd(ldwd, ld2wd, m)
      !! The `delta` weights.
   integer, intent(in) :: ldwd
      !! The leading dimension of array `wd`.
   integer, intent(in) :: ld2wd
      !! The second dimension of array `wd`.
   integer, intent(in) :: ifixb(np)
      !! The values designating whether the elements of `beta` are fixed at their input values or not.
   integer, intent(in) :: ifixx(ldifx, m)
      !! The values designating whether the elements of `x` are fixed at their input values or not.
   integer, intent(in) :: ldifx
      !! The leading dimension of array `ifixx`.
   real(kind=wp), intent(inout) :: betac(np)
      !! The current estimated values of the unfixed `beta`s.
   real(kind=wp), intent(out) :: betan(np)
      !! The new estimated values of the unfixed `beta`s.
   real(kind=wp), intent(in) :: betas(np)
      !! The saved estimated values of the unfixed `beta`s.
   real(kind=wp), intent(out) :: s(np)
      !! The step for `beta`.
   real(kind=wp), intent(in) :: delta(n, m)
      !! The estimated errors in the explanatory variables.
   real(kind=wp), intent(out) :: deltan(n, m)
      !! The new estimated errors in the explanatory variables.
   real(kind=wp), intent(inout) :: deltas(n, m)
      !! The saved estimated errors in the explanatory variables.
   real(kind=wp), intent(in) :: lower(np)
      !! The lower bound for unfixed `beta`s.
   real(kind=wp), intent(in) :: upper(np)
      !! The upper bound for unfixed `beta`s.
   real(kind=wp), intent(out) :: t(n, m)
      !! The step for `delta`.
   real(kind=wp), intent(inout) :: f(n, nq)
      !! The (weighted) estimated values of `epsilon`.
   real(kind=wp), intent(out) :: fn(n, nq)
      !! The new predicted values from the function.
   real(kind=wp), intent(out) :: fs(n, nq)
      !! The saved predicted values from the function.
   real(kind=wp), intent(out) :: fjacb(n, np, nq)
      !! The Jacobian with respect to `beta`.
   integer, intent(in) :: msgb(nq*np + 1)
      !! The error checking results for the Jacobian with respect to `beta`.
   real(kind=wp), intent(out) :: fjacd(n, m, nq)
      !! The Jacobian with respect to `delta`.
   integer, intent(in) :: msgd(nq*m + 1)
      !! The error checking results for the Jacobian with respect to `delta`.
   real(kind=wp), intent(in) :: ssf(np)
      !! The scaling values used for `beta`.
   real(kind=wp), intent(in) :: ss(np)
      !! The scaling values used for the unfixed `beta`s.
   real(kind=wp), intent(in) :: tt(ldtt, m)
      !! The scaling values used for `delta`.
   integer, intent(in) :: ldtt
      !! The leading dimension of array `tt`.
   real(kind=wp), intent(in) :: stpb(np)
      !! The relative step used for computing finite difference derivatives with respect to each `beta`.
   real(kind=wp), intent(in) :: stpd(ldstpd, m)
      !! The relative step used for computing finite difference derivatives with respect to `delta`.
   integer, intent(in) :: ldstpd
      !! The leading dimension of array `stpd`.
   real(kind=wp), intent(out) :: xplusd(n, m)
      !! The values of `x + delta`.
   real(kind=wp), intent(out) :: wrk(lwrk)
      !! A work array, _equivalenced_ to `wrk1` and `wrk2`.
   integer, intent(in) :: lwrk
      !! The length of vector `wrk`.
   real(kind=wp), intent(inout) :: work(lwork)
      !! The real (kind=wp) workspace.
   integer, intent(in) :: lwork
      !! The length of vector `work`.
   integer, intent(inout) :: iwork(liwork)
      !! The integer workspace.
   integer, intent(in) :: liwork
      !! The length of vector `iwork`.
   integer, intent(inout) :: info
      !! The variable designating why the computations were stopped.
   integer, intent(out) :: bound(np)
      !! The values of the bounds for `beta`.

   ! Local scalars
   real(kind=wp), parameter :: p0001 = 0.00010_wp, &
                               p1 = 0.1_wp, &
                               p25 = 0.25_wp, &
                               p5 = 0.5_wp, &
                               p75 = 0.75_wp
   real(kind=wp) :: actred, actrs, alpha, dirder, eta, olmavg, partol, pnorm, prered, &
                    prers, ratio, rcond, rnorm, rnormn, rnorms, rss, rvar, sstol, tau, &
                    taufac, temp, temp1, temp2, tsnorm

   integer, parameter :: ludflt = 6
   integer :: i, idf, iflag, int2, ipr, ipr1, ipr2, ipr2f, ipr3, irank, istop, istopc, &
              iwrk, j, jpvt, l, looped, lunr, lunrpt, maxit, neta, nfev, niter, njev, &
              nlms, nnzw, npp, npr, npu, omega, qraux, sd, u, vcv, wrk1, wrk2, wrk3, &
              wrk4, wrk5, wrk6
   logical :: access, anajac, cdjac, chkjac, cnvpar, cnvss, didvcv, dovcv, implct, initd, &
              intdbl, isodr, lstep, redoj, restrt

   ! Local arrays
   real(kind=wp) :: loweru(np), upperu(np), wss(3)

   ! External functions
   real(kind=wp), external :: ddot, dnrm2

   ! External subroutines
   external :: dacces, dcopy, devjac, dflags, dodlm, dodpcr, dodvcv, dunpac

   ! Interface blocks
   interface
      subroutine dwght(n, m, wt, ldwt, ld2wt, t, wtt)
         use odrpack_kinds, only: wp
         integer, intent(in) :: n, m, ldwt, ld2wt
         real(kind=wp), intent(in) :: t(:, :), wt(:, :, :)
         real(kind=wp), intent(out) :: wtt(:, :)
      end subroutine
   end interface

   ! Variable Definitions (alphabetically)
   !  ACCESS:  The variable designating whether information is to be accessed from the work
   !           arrays (ACCESS=TRUE) or stored in them (ACCESS=FALSE).
   !  ACTRED:  The actual relative reduction in the sum-of-squares.
   !  ACTRS:   The saved actual relative reduction in the sum-of-squares.
   !  ALPHA:   The Levenberg-Marquardt parameter.
   !  ANAJAC:  The variable designating whether the Jacobians are computed by finite
   !           differences (ANAJAC=FALSE) or not (ANAJAC=TRUE).
   !  BETA:    The function parameters.
   !  BETAC:   The current estimated values of the unfixed BETA'S.
   !  BETAN:   The new estimated values of the unfixed BETA'S.
   !  BETAS:   The saved estimated values of the unfixed BETA'S.
   !  CDJAC:   The variable designating whether the Jacobians are computed by central
   !           differences (cdjac=true) or by forward differences (CDJAC=FALSE).
   !  CHKJAC:  The variable designating whether the user supplied Jacobians are to be
   !           checked (CHKJAC=TRUE) or not (CHKJAC=FALSE).
   !  CNVPAR:  The variable designating whether parameter convergence was attained
   !           (CNVPAR=TRUE) or not (CNVPAR=FALSE).
   !  CNVSS:   The variable designating whether sum-of-squares convergence was attained
   !           (CNVSS=TRUE) or not (CNVSS=FALSE).
   !  DELTA:   The estimated errors in the explanatory variables.
   !  DELTAN:  The new estimated errors in the explanatory variables.
   !  DELTAS:  The saved estimated errors in the explanatory variables.
   !  DIDVCV:  The variable designating whether the covariance matrix was computed
   !           (DIDVCV=TRUE) or not (DIDVCV=FALSE).
   !  DIRDER:  The directional derivative.
   !  DOVCV:   The variable designating whether the covariance matrix should to be
   !           computed (DOVCV=TRUE) or not (DOVCV=FALSE).
   !  ETA:     The relative noise in the function results.
   !  F:       The (weighted) estimated values of EPSILON.
   !  FCN:     The user supplied subroutine for evaluating the model.
   !  FJACB:   The Jacobian with respect to BETA.
   !  FJACD:   The Jacobian with respect to DELTA.
   !  FN:      The new predicted values from the function.
   !  FS:      The saved predicted values from the function.
   !  FSTITR:  The variable designating whether this is the first iteration (FSTITR=TRUE)
   !           or not (FSTITR=FALSE).
   !  HEAD:    The variable designating whether the heading is to be printed (HEAD=TRUE) or
   !           not (HEAD=FALSE).
   !  I:       An indexing variable.
   !  IDF:     The degrees of freedom of the fit, equal to the number of observations with
   !           nonzero weighted derivatives minus the number of parameters being estimated.
   !  IFIXB:   The values designating whether the elements of BETA are fixed at their input
   !           values or not.
   !  IFIXX:   The values designating whether the elements of X are fixed at their input
   !           values or not.
   !  IFLAG:   The variable designating which report is to be printed.
   !  IMPLCT:  The variable designating whether the solution is by implicit ODR (IMPLCT=TRUE)
   !           or explicit ODR (IMPLCT=FALSE).
   !  INFO:    The variable designating why the computations were stopped.
   !  INITD:   The variable designating whether delta is initialized to zero (INITD=TRUE) or
   !           to the values in the first N by M elements of array work (INITD=FALSE).
   !  INT2:    The number of internal doubling steps taken.
   !  INTDBL:  The variable designating whether internal doubling is to be used (INTDBL=TRUE)
   !           or NOT (INTDBL=FALSE).
   !  IPR:     The values designating the length of the printed report.
   !  IPR1:    The value of the 4th digit (from the right) of iprint, which controls the
   !           initial summary report.
   !  IPR2:    The value of the 3rd digit (from the right) of iprint, which controls the
   !           final iteration report.
   !  IPR2F:   The value of the 2nd digit (from the right) of iprint, which controls the
   !           frequency of the iteration reports.
   !  IPR3:    The value of the 1st digit (from the right) of iprint, which controls the final
   !           summary report.
   !  IRANK:   The rank deficiency of the Jacobian wrt BETA.
   !  ISODR:   The variable designating whether the solution is by ODR (ISODR=TRUE) or
   !           OLS (ISODR=FALSE).
   !  ISTOP:   The variable designating whether there are problems computing the function
   !           at the current BETA and DELTA.
   !  ISTOPC:  The variable designating whether the computations were stoped due to some
   !           numerical error within routine  DODSTP.
   !  IWORK:   The integer work space.
   !  IWRK:    An index variable.
   !  J:       An index variable.
   !  JOB:     The variable controling problem initialization and computational method.
   !  JPVT:    The starting location in IWORK of array JPVT.
   !  L:       An index variable.
   !  LDIFX:   The leading dimension of array IFIXX.
   !  LDTT:    The leading dimension of array TT.
   !  LDWD:    The leading dimension of array WD.
   !  LDWE:    The leading dimension of array WE and WE1.
   !  LDX:     The leading dimension of array X.
   !  LDY:     The leading dimension of array Y.
   !  LD2WD:   The second dimension of array WD.
   !  LD2WE:   The second dimension of array WE and WE1.
   !  LIWORK:  The length of vector IWORK.
   !  LOOPED:  A counter used to determine how many times the subloop has been executed,
   !           where if the count becomes large enough the computations will be stopped.
   !  LOWERU:  The lower bound for unfixed BETAs.
   !  LSTEP:   The variable designating whether a successful step has been found (LSTEP=TRUE)
   !           or not (LSTEP=FALSE).
   !  LUDFLT:  The default logical unit number, used for computation reports to the screen.
   !  LUNR:    The logical unit number used for computation reports.
   !  LUNRPT:  The logical unit number used for computation reports.
   !  LWORK:   The length of vector WORK.
   !  LWRK:    The length of vector WRK.
   !  M:       The number of columns of data in the explanatory variable.
   !  MAXIT:   The maximum number of iterations allowed.
   !  MSGB:    The error checking results for the Jacobian wrt BETA.
   !  MSGD:    The error checking results for the Jacobian wrt DELTA.
   !  N:       The number of observations.
   !  NETA:    The number of accurate digits in the function results.
   !  NFEV:    The number of function evaluations.
   !  NITER:   The number of iterations taken.
   !  NJEV:    The number of Jacobian evaluations.
   !  NLMS:    The number of Levenberg-Marquardt steps taken.
   !  NNZW:    The number of nonzero weighted observations.
   !  NP:      The number of function parameters.
   !  NPP:     The number of function parameters being estimated.
   !  NPR:     The number of times the report is to be written.
   !  NPU:     The number of unfixed parameters.
   !  NQ:      The number of responses per observation.
   !  OLMAVG:  The average number of Levenberg-Marquardt steps per iteration.
   !  OMEGA:   The starting location in WORK of array OMEGA.
   !  PARTOL:  The parameter convergence stopping tolerance.
   !  PNORM:   The norm of the scaled estimated parameters.
   !  PRERED:  The predicted relative reduction in the sum-of-squares.
   !  PRERS:   The old predicted relative reduction in the sum-of-squares.
   !  PRTPEN:  The value designating whether the penalty parameter is to be printed in the
   !           iteration report (PRTPEN=TRUE) or not (PRTPEN=FALSE).
   !  QRAUX:   The starting location in array WORK of array QRAUX.
   !  RATIO:   The ratio of the actual relative reduction to the predicted relative reduction
   !           in the sum-of-squares.
   !  RCOND:   The approximate reciprocal condition of FJACB.
   !  REDOJ:   The variable designating whether the Jacobian matrix is to be recomputed for
   !           the computation of the covariance matrix (REDOJ=TRUE) or not (REDOJ=FALSE).
   !  RESTRT:  The variable designating whether the call is a restart (RESTRT=TRUE) or
   !           not (RESTRT=FALSE).
   !  RNORM:   The norm of the weighted errors.
   !  RNORMN:  The new norm of the weighted errors.
   !  RNORMS:  The saved norm of the weighted errors.
   !  RSS:     The residual sum of squares.
   !  RVAR:    The residual variance.
   !  S:       The step for BETA.
   !  SD:      The starting location in array work of array SD.
   !  SS:      The scaling values used for the unfixed BETAS.
   !  SSF:     The scaling values used for BETA.
   !  SSTOL:   The sum-of-squares convergence stopping tolerance.
   !  STPB:    The relative step used for computing finite difference derivatives with
   !           respect to each BETA.
   !  STPD:    The relative step used for computing finite difference derivatives with respect
   !           to DELTA.
   !  T:       The step for DELTA.
   !  TAU:     The trust region diameter.
   !  TAUFAC:  The factor used to compute the initial trust region diameter.
   !  TEMP:    A temporary storage location.
   !  TEMP1:   A temporary storage location.
   !  TEMP2:   A temporary storage location.
   !  TSNORM:  The norm of the scaled step.
   !  TT:      The scaling values used for DELTA.
   !  U:       The starting location in array WORK of array U.
   !  UPPERU:  The upper bound for unfixed BETAs.
   !  VCV:     The starting location in array WORK of array VCV.
   !  WE:      The EPSILON weights.
   !  WE1:     The square root of the EPSILON weights.
   !  WD:      The DELTA weights.
   !  WORK:    The REAL (KIND=wp) work space.
   !  WSS:     The sum-of-squares of the weighted EPSILONS and DELTAS, the sum-of-squares
   !           of the weighted DELTAS, and the sum-of-squares of the weighted EPSILONS.
   !  WRK:     A work array, equivalenced to WRK1 and WRK2
   !  WRK1:    The starting location in array WORK of array WRK1.
   !  WRK2:    The starting location in array WORK of array WRK2.
   !  WRK3:    The starting location in array WORK of array WRK3.
   !  WRK4:    The starting location in array WORK of array WRK4.
   !  WRK5:    The starting location in array WORK of array WRK5.
   !  WRK6:    The starting location in array WORK of array WRK6.
   !  X:       The explanatory variable.
   !  XPLUSD:  The values of X + DELTA.
   !  Y:       The dependent variable. Unused when the model is implicit.

   ! Initialize necessary variables
   call dpack(np, npu, loweru, lower, ifixb)
   call dpack(np, npu, upperu, upper, ifixb)
   call dflags(job, restrt, initd, dovcv, redoj, &
               anajac, cdjac, chkjac, isodr, implct)
   access = .true.
   call dacces(n, m, np, nq, ldwe, ld2we, &
               work, lwork, iwork, liwork, &
               access, isodr, &
               jpvt, omega, u, qraux, sd, vcv, &
               wrk1, wrk2, wrk3, wrk4, wrk5, wrk6, &
               nnzw, npp, &
               job, partol, sstol, maxit, taufac, eta, neta, &
               lunrpt, ipr1, ipr2, ipr2f, ipr3, &
               wss, rvar, idf, &
               tau, alpha, niter, nfev, njev, int2, olmavg, &
               rcond, irank, actrs, pnorm, prers, rnorms, istop)
   rnorm = sqrt(wss(1))

   didvcv = .false.
   intdbl = .false.
   lstep = .true.

   ! Print initial summary if desired
   if (ipr1 .ne. 0 .and. lunrpt .ne. 0) then
      iflag = 1
      if (ipr1 .ge. 3 .and. lunrpt .ne. ludflt) then
         npr = 2
      else
         npr = 1
      end if
      if (ipr1 .ge. 6) then
         ipr = 2
      else
         ipr = 2 - mod(ipr1, 2)
      end if
      lunr = lunrpt
      do i = 1, npr
         call dodpcr(ipr, lunr, &
                     head, prtpen, fstitr, didvcv, iflag, &
                     n, m, np, nq, npp, nnzw, &
                     msgb, msgd, beta, y, ldy, x, ldx, delta, &
                     we, ldwe, ld2we, wd, ldwd, ld2wd, &
                     ifixb, ifixx, ldifx, &
                     lower, upper, &
                     ssf, tt, ldtt, stpb, stpd, ldstpd, &
                     job, neta, taufac, sstol, partol, maxit, &
                     wss, rvar, idf, work(sd), &
                     niter, nfev, njev, actred, prered, &
                     tau, pnorm, alpha, f, rcond, irank, info, istop)
         if (ipr1 .ge. 5) then
            ipr = 2
         else
            ipr = 1
         end if
         lunr = ludflt
      end do
   end if

   ! Stop if initial estimates are exact solution
   if (rnorm .eq. zero) then
      info = 1
      olmavg = zero
      istop = 0
      goto 150
   end if

   ! Stop if number of iterations already equals maximum permitted
   if (restrt .and. &
       (niter .ge. maxit)) then
      istop = 0
      goto 150
   elseif (niter .ge. maxit) then
      info = 4
      istop = 0
      goto 150
   end if

   ! MAIN LOOP
100 continue

   niter = niter + 1
   rnorms = rnorm
   looped = 0

   ! Evaluate jacobian using best estimate of function (FS)
   if ((niter .eq. 1) .and. &
       (anajac .and. chkjac)) then
      istop = 0
   else
      call devjac(fcn, &
                  anajac, cdjac, &
                  n, m, np, nq, &
                  betac, beta, stpb, &
                  ifixb, ifixx, ldifx, &
                  x, ldx, delta, xplusd, stpd, ldstpd, &
                  ssf, tt, ldtt, neta, fs, &
                  t, work(wrk1), work(wrk2), work(wrk3), work(wrk6), &
                  fjacb, isodr, fjacd, we1, ldwe, ld2we, &
                  njev, nfev, istop, info, &
                  lower, upper)
   end if
   if (istop .ne. 0) then
      info = 51000
      goto 200
   elseif (info .eq. 50300) then
      goto 200
   end if

   ! SUB-LOOP for internal doubling or computing new step when old failed
110 continue

   ! Compute steps S and T
   if (looped .gt. 100) then
      info = 60000
      goto 200
   else
      looped = looped + 1
      call dodlm(n, m, np, nq, npp, &
                 f, fjacb, fjacd, &
                 wd, ldwd, ld2wd, ss, tt, ldtt, delta, &
                 alpha, tau, eta, isodr, &
                 work(wrk6), work(omega), &
                 work(u), work(qraux), iwork(jpvt), &
                 s, t, nlms, rcond, irank, &
                 work(wrk1), work(wrk2), work(wrk3), work(wrk4), &
                 work(wrk5), wrk, lwrk, istopc)
   end if
   if (istopc .ne. 0) then
      info = istopc
      goto 200
   end if
   olmavg = olmavg + nlms

   ! Compute BETAN = BETAC + S
   !         DELTAN = DELTA + T
   betan = betac + s
   if (isodr) deltan = delta + t

   ! Project the step wrt the bounds
   do i = 1, npu
      if (loweru(i) .eq. upperu(i)) then
         betan(i) = upperu(i)
         s(i) = upperu(i) - betac(i)
         bound(i) = 3
      elseif (betan(i) .le. loweru(i)) then
         betan(i) = loweru(i)
         s(i) = loweru(i) - betac(i)
         bound(i) = 2
      elseif (betan(i) .ge. upperu(i)) then
         betan(i) = upperu(i)
         s(i) = upperu(i) - betac(i)
         bound(i) = 1
      else
         bound(i) = 0
      end if
   end do

   ! Compute norm of scaled steps S and T (TSNORM)
   call dwght(npp, 1, reshape(ss, [npp, 1, 1]), npp, 1, &
              reshape(s, [npp, 1]), tempret(1:npp, 1:1))
   wrk(1:npp) = tempret(1:npp, 1)
   if (isodr) then
      call dwght(n, m, reshape(tt, [ldtt, 1, m]), ldtt, 1, t, tempret(1:n, 1:m))
      wrk(npp + 1:npp + 1 + n*m - 1) = reshape(tempret(1:n, 1:m), [n*m])
      tsnorm = dnrm2(npp + n*m, wrk, 1)
   else
      tsnorm = dnrm2(npp, wrk, 1)
   end if

   ! Compute scaled predicted reduction
   iwrk = 0
   do l = 1, nq
      do i = 1, n
         iwrk = iwrk + 1
         wrk(iwrk) = ddot(npp, fjacb(i, 1, l), n, s, 1)
         if (isodr) wrk(iwrk) = wrk(iwrk) + ddot(m, fjacd(i, 1, l), n, t(i, 1), n)
      end do
   end do
   if (isodr) then
      call dwght(n, m, wd, ldwd, ld2wd, t, tempret(1:n, 1:m))
      wrk(n*nq + 1:n*nq + 1 + n*m - 1) = reshape(tempret(1:n, 1:m), [n*m])
      temp1 = ddot(n*nq, wrk, 1, wrk, 1) + ddot(n*m, t, 1, wrk(n*nq + 1), 1)
      temp1 = sqrt(temp1)/rnorm
   else
      temp1 = dnrm2(n*nq, wrk, 1)/rnorm
   end if
   temp2 = sqrt(alpha)*tsnorm/rnorm
   prered = temp1**2 + temp2**2/p5

   dirder = -(temp1**2 + temp2**2)

   ! Evaluate predicted values at new point
   call dunpac(np, betan, beta, ifixb)
   xplusd = x(1:n, :) + deltan
   istop = 0
   call fcn(n, m, np, nq, &
            n, m, np, &
            beta, xplusd, &
            ifixb, ifixx, ldifx, &
            002, fn, work(wrk6), work(wrk1), &
            istop)
   if (istop .eq. 0) then
      nfev = nfev + 1
   end if

   if (istop .lt. 0) then
      ! Set INFO to indicate user has stopped the computations in FCN
      info = 51000
      goto 200
   elseif (istop .gt. 0) then
      ! Set norm to indicate step should be rejected
      rnormn = rnorm/(p1*p75)
   else
      ! Compute norm of new weighted EPSILONS and weighted DELTAS (RNORMN)
      if (implct) then
         call dcopy(n*nq, fn, 1, wrk, 1)
      else
         !call dxmy( n, nq, fn, n, y, ldy, wrk, n)
         wrk(1:n*nq) = reshape(fn - y(1:n, :), [n*nq])
      end if
      call dwght(n, nq, we1, ldwe, ld2we, reshape(wrk, [n, nq]), &
                 tempret(1:n, 1:nq))
      wrk(1:n*nq) = reshape(tempret(1:n, 1:nq), [n*nq])
      if (isodr) then
         call dwght(n, m, wd, ldwd, ld2wd, deltan, tempret(1:n, 1:m))
         wrk(n*nq + 1:n*nq + 1 + n*m - 1) = reshape(tempret(1:n, 1:m), [ &
                                                    n*m])
         rnormn = sqrt(ddot(n*nq, wrk, 1, wrk, 1) + &
                       ddot(n*m, deltan, 1, wrk(n*nq + 1), 1))
      else
         rnormn = dnrm2(n*nq, wrk, 1)
      end if
   end if

   ! Compute scaled actual reduction
   if (p1*rnormn .lt. rnorm) then
      actred = one - (rnormn/rnorm)**2
   else
      actred = -one
   end if

   ! Compute ratio of actual reduction to predicted reduction
   if (prered .eq. zero) then
      ratio = zero
   else
      ratio = actred/prered
   end if

   ! Check on lack of reduction in internal doubling case
   if (intdbl .and. (ratio .lt. p0001 .or. rnormn .gt. rnorms)) then
      istop = 0
      tau = tau*p5
      alpha = alpha/p5
      call dcopy(npp, betas, 1, betan, 1)
      call dcopy(n*m, deltas, 1, deltan, 1)
      call dcopy(n*nq, fs, 1, fn, 1)
      actred = actrs
      prered = prers
      rnormn = rnorms
      ratio = p5
   end if

   ! Update step bound
   intdbl = .false.
   if (ratio .lt. p25) then
      if (actred .ge. zero) then
         temp = p5
      else
         temp = p5*dirder/(dirder + p5*actred)
      end if
      if (p1*rnormn .ge. rnorm .or. temp .lt. p1) then
         temp = p1
      end if
      tau = temp*min(tau, tsnorm/p1)
      alpha = alpha/temp
   elseif (alpha .eq. zero) then
      tau = tsnorm/p5

   elseif (ratio .ge. p75 .and. nlms .le. 11) then
      ! Step qualifies for internal doubling
      !  - Update TAU and ALPHA
      !  - Save information for current point

      intdbl = .true.

      tau = tsnorm/p5
      alpha = alpha*p5

      call dcopy(npp, betan, 1, betas, 1)
      call dcopy(n*m, deltan, 1, deltas, 1)
      call dcopy(n*nq, fn, 1, fs, 1)
      actrs = actred
      prers = prered
      rnorms = rnormn
   end if

   ! If internal doubling, skip convergence checks
   if (intdbl .and. tau .gt. zero) then
      int2 = int2 + 1
      goto 110
   end if

   ! Check acceptance
   if (ratio .ge. p0001) then
      call dcopy(n*nq, fn, 1, fs, 1)
      if (implct) then
         call dcopy(n*nq, fs, 1, f, 1)
      else
         !call dxmy( n, nq, fs, n, y, ldy, f, n)
         f = fs - y(1:n, :)
      end if
      call dwght(n, nq, we1, ldwe, ld2we, f, tempret(1:n, 1:nq))
      f(1:n, 1:nq) = tempret(1:n, 1:nq)
      call dcopy(npp, betan, 1, betac, 1)
      call dcopy(n*m, deltan, 1, delta, 1)
      rnorm = rnormn
      call dwght(npp, 1, reshape(ss, [npp, 1, 1]), npp, 1, &
                 reshape(betac, [npp, 1]), tempret(1:npp, 1:1))
      wrk(1:npp) = tempret(1:npp, 1)
      if (isodr) then
         call dwght(n, m, reshape(tt, [ldtt, 1, m]), ldtt, 1, &
                    delta, tempret(1:n, 1:m))
         wrk(npp + 1:npp + 1 + n*m - 1) = reshape(tempret(1:n, 1:m), [n*m])
         pnorm = dnrm2(npp + n*m, wrk, 1)
      else
         pnorm = dnrm2(npp, wrk, 1)
      end if
      lstep = .true.
   else
      lstep = .false.
   end if

   ! Test convergence
   info = 0
   cnvss = rnorm .eq. zero &
           .or. &
           (abs(actred) .le. sstol .and. &
            prered .le. sstol .and. &
            p5*ratio .le. one)
   cnvpar = (tau .le. partol*pnorm) .and. (.not. implct)
   if (cnvss) info = 1
   if (cnvpar) info = 2
   if (cnvss .and. cnvpar) info = 3

   ! Print iteration report
   if (info .ne. 0 .or. lstep) then
      if (ipr2 .ne. 0 .and. ipr2f .ne. 0 .and. lunrpt .ne. 0) then
         if (ipr2f .eq. 1 .or. mod(niter, ipr2f) .eq. 1) then
            iflag = 2
            call dunpac(np, betac, beta, ifixb)
            wss(1) = rnorm*rnorm
            if (ipr2 .ge. 3 .and. lunrpt .ne. ludflt) then
               npr = 2
            else
               npr = 1
            end if
            if (ipr2 .ge. 6) then
               ipr = 2
            else
               ipr = 2 - mod(ipr2, 2)
            end if
            lunr = lunrpt
            do i = 1, npr
               call dodpcr(ipr, lunr, &
                           head, prtpen, fstitr, didvcv, iflag, &
                           n, m, np, nq, npp, nnzw, &
                           msgb, msgd, beta, y, ldy, x, ldx, delta, &
                           we, ldwe, ld2we, wd, ldwd, ld2wd, &
                           ifixb, ifixx, ldifx, &
                           lower, upper, &
                           ssf, tt, ldtt, stpb, stpd, ldstpd, &
                           job, neta, taufac, sstol, partol, maxit, &
                           wss, rvar, idf, work(sd), &
                           niter, nfev, njev, actred, prered, &
                           tau, pnorm, alpha, f, rcond, irank, info, istop)
               if (ipr2 .ge. 5) then
                  ipr = 2
               else
                  ipr = 1
               end if
               lunr = ludflt
            end do
            fstitr = .false.
            prtpen = .false.
         end if
      end if
   end if

   ! Check if finished
   if (info .eq. 0) then

      if (lstep) then
         ! Begin next interation unless a stopping criteria has been met
         if (niter .ge. maxit) then
            info = 4
         else
            goto 100
         end if
      else
         ! Step failed - recompute unless a stopping criteria has been met
         goto 110
      end if

   end if

150 continue

   if (istop .gt. 0) info = info + 100

   ! Store unweighted EPSILONS and X+DELTA to return to user
   if (implct) then
      call dcopy(n*nq, fs, 1, f, 1)
   else
      !call dxmy( n, nq, fs, n, y, ldy, f, n)
      f = fs - y(1:n, :)
   end if
   call dunpac(np, betac, beta, ifixb)
   xplusd = x(1:n, :) + delta

   ! Compute covariance matrix of estimated parameters in upper NP by NP portion
   ! of WORK(VCV) if requested
   if (dovcv .and. istop .eq. 0) then

      ! Re-evaluate Jacobian at final solution, if requested
      ! Otherwise, Jacobian from beginning of last iteration will be used
      ! to compute covariance matrix
      if (redoj) then
         call devjac(fcn, &
                     anajac, cdjac, &
                     n, m, np, nq, &
                     betac, beta, stpb, &
                     ifixb, ifixx, ldifx, &
                     x, ldx, delta, xplusd, stpd, ldstpd, &
                     ssf, tt, ldtt, neta, fs, &
                     t, work(wrk1), work(wrk2), work(wrk3), work(wrk6), &
                     fjacb, isodr, fjacd, we1, ldwe, ld2we, &
                     njev, nfev, istop, info, &
                     lower, upper)

         if (istop .ne. 0) then
            info = 51000
            goto 200
         elseif (info .eq. 50300) then
            goto 200
         end if
      end if

      if (implct) then
         call dwght(n, m, wd, ldwd, ld2wd, delta, tempret(1:n, 1:m))
         wrk(n*nq + 1:n*nq + 1 + n*m - 1) = reshape(tempret(1:n, 1:m), [n*m])
         rss = ddot(n*m, delta, 1, wrk(n*nq + 1), 1)
      else
         rss = rnorm*rnorm
      end if
      if (redoj .or. niter .ge. 1) then
         call dodvcv(n, m, np, nq, npp, &
                     f, fjacb, fjacd, &
                     wd, ldwd, ld2wd, ssf, ss, tt, ldtt, delta, &
                     eta, isodr, &
                     work(vcv), work(sd), &
                     work(wrk6), work(omega), &
                     work(u), work(qraux), iwork(jpvt), &
                     s, t, irank, rcond, rss, idf, rvar, ifixb, &
                     work(wrk1), work(wrk2), work(wrk3), work(wrk4), &
                     work(wrk5), wrk, lwrk, istopc)
         if (istopc .ne. 0) then
            info = istopc
            goto 200
         end if
         didvcv = .true.
      end if

   end if

   ! Set JPVT to indicate dropped, fixed and estimated parameters
200 do i = 0, np - 1
      work(wrk3 + i) = iwork(jpvt + i)
      iwork(jpvt + i) = -2
   end do
   if (redoj .or. niter .ge. 1) then
      do i = 0, npp - 1
         j = int(work(wrk3 + i)) - 1
         if (i .le. npp - irank - 1) then
            iwork(jpvt + j) = 1
         else
            iwork(jpvt + j) = -1
         end if
      end do
      if (npp .lt. np) then
         j = npp - 1
         do i = np - 1, 0, -1
            if (ifixb(i + 1) .eq. 0) then
               iwork(jpvt + i) = 0
            else
               iwork(jpvt + i) = iwork(jpvt + j)
               j = j - 1
            end if
         end do
      end if
   end if

   ! Store various scalars in work arrays for return to user
   if (niter .ge. 1) then
      olmavg = olmavg/niter
   else
      olmavg = zero
   end if

   ! Compute weighted sums of squares for return to user
   call dwght(n, nq, we1, ldwe, ld2we, f, tempret(1:n, 1:nq))
   wrk(1:n*nq) = reshape(tempret(1:n, 1:nq), [n*nq])
   wss(3) = ddot(n*nq, wrk, 1, wrk, 1)
   if (isodr) then
      call dwght(n, m, wd, ldwd, ld2wd, delta, tempret(1:n, 1:m))
      wrk(n*nq + 1:n*nq + 1 + n*m - 1) = reshape(tempret(1:n, 1:m), [n*m])
      wss(2) = ddot(n*m, delta, 1, wrk(n*nq + 1), 1)
   else
      wss(2) = zero
   end if
   wss(1) = wss(2) + wss(3)

   access = .false.
   call dacces(n, m, np, nq, ldwe, ld2we, &
               work, lwork, iwork, liwork, &
               access, isodr, &
               jpvt, omega, u, qraux, sd, vcv, &
               wrk1, wrk2, wrk3, wrk4, wrk5, wrk6, &
               nnzw, npp, &
               job, partol, sstol, maxit, taufac, eta, neta, &
               lunrpt, ipr1, ipr2, ipr2f, ipr3, &
               wss, rvar, idf, &
               tau, alpha, niter, nfev, njev, int2, olmavg, &
               rcond, irank, actrs, pnorm, prers, rnorms, istop)

   ! Encode existance of questionable results into info
   if (info .le. 9 .or. info .ge. 60000) then
      if (msgb(1) .eq. 1 .or. msgd(1) .eq. 1) then
         info = info + 1000
      end if
      if (istop .ne. 0) then
         info = info + 100
      end if
      if (irank .ge. 1) then
         if (npp .gt. irank) then
            info = info + 10
         else
            info = info + 20
         end if
      end if
   end if

   ! Print final summary
   if (ipr3 .ne. 0 .and. lunrpt .ne. 0) then
      iflag = 3

      if (ipr3 .ge. 3 .and. lunrpt .ne. ludflt) then
         npr = 2
      else
         npr = 1
      end if
      if (ipr3 .ge. 6) then
         ipr = 2
      else
         ipr = 2 - mod(ipr3, 2)
      end if
      lunr = lunrpt
      do i = 1, npr
         call dodpcr(ipr, lunr, &
                     head, prtpen, fstitr, didvcv, iflag, &
                     n, m, np, nq, npp, nnzw, &
                     msgb, msgd, beta, y, ldy, x, ldx, delta, &
                     we, ldwe, ld2we, wd, ldwd, ld2wd, &
                     iwork(jpvt), ifixx, ldifx, &
                     lower, upper, &
                     ssf, tt, ldtt, stpb, stpd, ldstpd, &
                     job, neta, taufac, sstol, partol, maxit, &
                     wss, rvar, idf, work(sd), &
                     niter, nfev, njev, actred, prered, &
                     tau, pnorm, alpha, f, rcond, irank, info, istop)
         if (ipr3 .ge. 5) then
            ipr = 2
         else
            ipr = 1
         end if
         lunr = ludflt
      end do
   end if

end subroutine dodmn
