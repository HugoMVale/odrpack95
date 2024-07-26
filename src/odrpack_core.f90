module odrpack_core
!! Core mathematical routines, except drivers, and BLAS/LINPACK.

   use odrpack_kinds, only: wp
   implicit none

   abstract interface
      subroutine fcn_t(n, m, np, nq, ldn, ldm, ldnp, beta, xplusd, ifixb, ifixx, ldifx, &
                       ideval, f, fjacb, fjacd, istop)
      !! User-supplied subroutine for evaluating the model.
         import :: wp
         implicit none
         integer, intent(in) :: n
            !! Number of observations.
         integer, intent(in) :: m
            !! Number of columns of data in the independent variable.
         integer, intent(in) :: np
            !! Number of function parameters.
         integer, intent(in) :: nq
            !! Number of responses per observation.
         integer, intent(in) :: ldn
            !! Leading dimension declarator equal or exceeding `n`.
         integer, intent(in) :: ldm
            !! Leading dimension declarator equal or exceeding `m`.
         integer, intent(in) :: ldnp
            !! Leading dimension declarator equal or exceeding `np`.
         real(wp), intent(in) :: beta(np)
            !! Current values of parameters.
         real(wp), intent(in) :: xplusd(ldn, m)
            !! Current value of explanatory variable, i.e., `x + delta`.
         integer, intent(in) :: ifixb(np)
            !! Indicators for "fixing" parameters (`beta`).
         integer, intent(in) :: ifixx(ldifx, m)
            !! Indicators for "fixing" explanatory variable (`x`).
         integer, intent(in) :: ldifx
            !! Leading dimension of array `ifixx`.
         integer, intent(in) :: ideval
            !! Indicator for selecting computation to be performed.
         real(wp), intent(out) :: f(ldn, nq)
            !! Predicted function values.
         real(wp), intent(out) :: fjacb(ldn, ldnp, nq)
            !! Jacobian with respect to `beta`.
         real(wp), intent(out) :: fjacd(ldn, ldm, nq)
            !! Jacobian with respect to errors `delta`.
         integer, intent(out) :: istop
            !! Stopping condition, with meaning as follows. 0 means current `beta` and
            !! `x+delta` were acceptable and values were computed successfully. 1 means current
            !! `beta` and `x+delta` are not acceptable;  'odrpack' should select values closer
            !! to most recently used values if possible. -1 means current `beta` and `x+delta`
            !! are not acceptable; 'odrpack' should stop.
      end subroutine fcn_t
   end interface

contains

   subroutine dodlm &
      (n, m, np, nq, npp, &
       f, fjacb, fjacd, &
       wd, ldwd, ld2wd, ss, tt, ldtt, delta, &
       alpha2, tau, epsfcn, isodr, &
       tfjacb, omega, u, qraux, jpvt, &
       s, t, nlms, rcond, irank, &
       wrk1, wrk2, wrk3, wrk4, wrk5, wrk, lwrk, tempret, istopc)
   !! Compute Levenberg-Marquardt parameter and steps `s` and `t` using analog of the
   !! trust-region Levenberg-Marquardt algorithm.
      ! Routines Called  DDOT, DNRM2, DODSTP, DSCALE, DWGHT
      ! Date Written   860529   (YYMMDD)
      ! Revision Date  920619   (YYMMDD)

      use odrpack_kinds, only: zero

      integer, intent(in) :: n
         !! The number of observations.
      integer, intent(in) :: m
         !! The number of columns of data in the explanatory variable.
      integer, intent(in) :: np
         !! The number of function parameters.
      integer, intent(in) :: nq
         !! The number of responses per observation.
      integer, intent(in) :: npp
         !! The number of function parameters being estimated.
      real(wp), intent(in) :: f(n, nq)
         !! The (weighted) estimated values of `epsilon`.
      real(wp), intent(in) :: fjacb(n, np, nq)
         !! The Jacobian with respect to `beta`.
      real(wp), intent(in) :: fjacd(n, m, nq)
         !! The Jacobian with respect to `delta`.
      real(wp), intent(in) :: wd(ldwd, ld2wd, m)
         !! The `delta` weights.
      integer, intent(in) :: ldwd
         !! The leading dimension of array `wd`.
      integer, intent(in) :: ld2wd
         !! The second dimension of array `wd`.
      real(wp), intent(in) :: ss(np)
         !! The scaling values used for the unfixed `beta`s.
      real(wp), intent(in) :: tt(ldtt, m)
         !! The scale used for the `delta`s.
      integer, intent(in) :: ldtt
         !! The leading dimension of array `tt`.
      real(wp), intent(in) :: delta(n, m)
         !! The estimated errors in the explanatory variables.
      real(wp), intent(inout) :: alpha2
         !! The current Levenberg-Marquardt parameter.
      real(wp), intent(inout) :: tau
         !! The trust region diameter.
      real(wp), intent(in) :: epsfcn
         !! The function's precision.
      logical, intent(in) :: isodr
         !! The variable designating whether the solution is by ODR (`isodr = .true.`)
         !! or by OLS (`isodr = .false.`).
      real(wp), intent(out) :: tfjacb(n, nq, np)
         !! The array `omega*fjacb`.
      real(wp), intent(out) :: omega(nq, nq)
         !! The array `(I-fjacd*inv(p)*trans(fjacd))**(-1/2)`.
      real(wp), intent(out) :: u(np)
         !! The approximate null vector for tfjacb.
      real(wp), intent(out) :: qraux(np)
         !! The array required to recover the orthogonal part of the Q-R decomposition.
      integer, intent(out) :: jpvt(np)
         !! The pivot vector.
      real(wp), intent(out) :: s(np)
         !! The step for `beta`.
      real(wp), intent(out) :: t(n, m)
         !! The step for `delta`.
      integer, intent(out) :: nlms
         !! The number of Levenberg-Marquardt steps taken.
      real(wp), intent(out) :: rcond
         !! The approximate reciprocal condition of `tfjacb`.
      integer, intent(out) :: irank
         !! The rank deficiency of the Jacobian wrt `beta`.
      real(wp), intent(out) :: wrk1(n, nq, m)
         !! A work array of `(n, nq, m)` elements.
      real(wp), intent(out) :: wrk2(n, nq)
         !! A work array of `(n, nq)` elements.
      real(wp), intent(out) :: wrk3(np)
         !! A work array of `(np)` elements.
      real(wp), intent(out) :: wrk4(m, m)
         !! A work array of `(m, m)` elements.
      real(wp), intent(out) :: wrk5(m)
         !! A work array of `(m)` elements.
      real(wp), intent(out) :: wrk(lwrk)
         !! A work array of `(lwrk)` elements, _equivalenced_ to `wrk1` and `wrk2`.
      integer, intent(in) :: lwrk
         !! The length of vector `wrk`.
      real(wp), intent(inout) :: tempret(:, :)
         !! Temporary work array for holding return values before copying to a lower rank array.
      integer, intent(out) :: istopc
         !! The variable designating whether the computations were stopped due to some other
         !! numerical error detected within subroutine `dodstp`.

      ! Local scalars
      real(wp), parameter :: p001 = 0.001_wp, p1 = 0.1_wp
      real(wp) :: alpha1, alphan, bot, phi1, phi2, sa, top
      integer :: i, iwrk, j, k
      logical :: forvcv

      ! External BLAS/LINPACK procedures
      real(wp), external :: ddot, dnrm2

      ! Variable Definitions (alphabetically)
      !  ALPHAN:  The new Levenberg-Marquardt parameter.
      !  ALPHA1:  The previous Levenberg-Marquardt parameter.
      !  ALPHA2:  The current Levenberg-Marquardt parameter.
      !  BOT:     The lower limit for setting ALPHA.
      !  DELTA:   The estimated errors in the explanatory variables.
      !  EPSFCN:  The function's precision.
      !  F:       The (weighted) estimated values of EPSILON.
      !  FJACB:   The Jacobian with respect to BETA.
      !  FJACD:   The Jacobian with respect to DELTA.
      !  FORVCV:  The variable designating whether this subroutine was called to set up for the
      !           covariance matrix computations (FORVCV=TRUE) or not (FORVCV=FALSE).
      !  I:       An indexing variable.
      !  IRANK:   The rank deficiency of the Jacobian wrt BETA.
      !  ISODR:   The variable designating whether the solution is by ODR (ISODR=TRUE) or
      !           by OLS (ISODR=FALSE).
      !  ISTOPC:  The variable designating whether the computations were stoped due to some
      !           other numerical error detected within subroutine DODSTP.
      !  IWRK:    An indexing variable.
      !  J:       An indexing variable.
      !  K:       An indexing variable.
      !  L:       An indexing variable.
      !  JPVT:    The pivot vector.
      !  LDTT:    The leading dimension of array TT.
      !  LDWD:    The leading dimension of array WD.
      !  LD2WD:   The second dimension of array WD.
      !  LWRK:    The length of vector WRK.
      !  M:       The number of columns of data in the explanatory variable.
      !  N:       The number of observations.
      !  NLMS:    The number of Levenberg-Marquardt steps taken.
      !  NP:      The number of function parameters.
      !  NPP:     The number of function parameters being estimated.
      !  NQ:      The number of responses per observation.
      !  OMEGA:   The array (I-FJACD*INV(P)*trans(FJACD))**(-1/2)  where
      !            P = trans(FJACD)*FJACD + D**2 + ALPHA*TT**2
      !  P001:    The value 0.001E0_wp
      !  P1:      The value 0.1E0_wp
      !  PHI1:    The previous difference between the norm of the scaled step and the trust
      !           region diameter.
      !  PHI2:    The current difference between the norm of the scaled step and the trust region
      !           diameter.
      !  QRAUX:   The array required to recover the orthogonal part of the Q-R decomposition.
      !  RCOND:   The approximate reciprocal condition of TFJACB.
      !  S:       The step for BETA.
      !  SA:      The scalar PHI2*(ALPHA1-ALPHA2)/(PHI1-PHI2).
      !  SS:      The scaling values used for the unfixed BETAS.
      !  T:       The step for DELTA.
      !  TAU:     The trust region diameter.
      !  TFJACB:  The array OMEGA*FJACB.
      !  TOP:     The upper limit for setting ALPHA.
      !  TT:      The scale used for the DELTA'S.
      !  U:       The approximate null vector for TFJACB.
      !  WD:      The DELTA weights.
      !  WRK:     A work array of (LWRK) elements, equivalenced to WRK1 and WRK2.
      !  WRK1:    A work array of (N by NQ by M) elements.
      !  WRK2:    A work array of (N by NQ) elements.
      !  WRK3:    A work array of (NP) elements.
      !  WRK4:    A work array of (M by M) elements.
      !  WRK5:    A work array of (M) elements.

      forvcv = .false.
      istopc = 0

      ! Compute full Gauss-Newton step (ALPHA=0)
      alpha1 = zero
      call dodstp(n, m, np, nq, npp, &
                  f, fjacb, fjacd, &
                  wd, ldwd, ld2wd, ss, tt, ldtt, delta, &
                  alpha1, epsfcn, isodr, &
                  tfjacb, omega, u, qraux, jpvt, &
                  s, t, phi1, irank, rcond, forvcv, &
                  wrk1, wrk2, wrk3, wrk4, wrk5, wrk, lwrk, tempret, istopc)
      if (istopc /= 0) then
         return
      end if

      ! Initialize TAU if necessary
      if (tau < zero) then
         tau = abs(tau)*phi1
      end if

      ! Check if full Gauss-Newton step is optimal
      if ((phi1 - tau) <= p1*tau) then
         nlms = 1
         alpha2 = zero
         return
      end if

      ! Full Gauss-Newton step is outside trust region - find locally constrained optimal step
      phi1 = phi1 - tau

      ! Initialize upper and lower bounds for ALPHA
      bot = zero

      do k = 1, npp
         tfjacb(1:n, 1:nq, k) = fjacb(1:n, k, 1:nq)
         wrk(k) = ddot(n*nq, tfjacb(1, 1, k), 1, f(1, 1), 1)
      end do
      call dscale(npp, 1, ss, npp, wrk, npp, wrk, npp) ! work is input (as t) and output (as sclt)

      if (isodr) then
         call dwght(n, m, wd, ldwd, ld2wd, delta, tempret(1:n, 1:m))
         wrk(npp + 1:npp + 1 + n*m - 1) = reshape(tempret(1:n, 1:m), (/n*m/))
         iwrk = npp
         do j = 1, m
            do i = 1, n
               iwrk = iwrk + 1
               wrk(iwrk) = wrk(iwrk) + ddot(nq, fjacd(i, j, 1), n*m, f(i, 1), n)
            end do
         end do
         call dscale(n, m, tt, ldtt, wrk(npp + 1), n, wrk(npp + 1), n)
         top = dnrm2(npp + n*m, wrk, 1)/tau
      else
         top = dnrm2(npp, wrk, 1)/tau
      end if

      if (alpha2 > top .or. alpha2 == zero) then
         alpha2 = p001*top
      end if

      ! Main loop

      do i = 1, 10

         ! Compute locally constrained steps S and T and PHI(ALPHA) for current value of ALPHA
         call dodstp(n, m, np, nq, npp, &
                     f, fjacb, fjacd, &
                     wd, ldwd, ld2wd, ss, tt, ldtt, delta, &
                     alpha2, epsfcn, isodr, &
                     tfjacb, omega, u, qraux, jpvt, &
                     s, t, phi2, irank, rcond, forvcv, &
                     wrk1, wrk2, wrk3, wrk4, wrk5, wrk, lwrk, tempret, istopc)
         if (istopc /= 0) then
            return
         end if
         phi2 = phi2 - tau

         ! Check whether current step is optimal
         if (abs(phi2) <= p1*tau .or. (alpha2 == bot .and. phi2 < zero)) then
            nlms = i + 1
            return
         end if

         ! Current step is not optimaL

         ! Update bounds for ALPHA and compute new ALPHA
         if (phi1 - phi2 == zero) then
            nlms = 12
            return
         end if
         sa = phi2*(alpha1 - alpha2)/(phi1 - phi2)
         if (phi2 < zero) then
            top = min(top, alpha2)
         else
            bot = max(bot, alpha2)
         end if
         if (phi1*phi2 > zero) then
            bot = max(bot, alpha2 - sa)
         else
            top = min(top, alpha2 - sa)
         end if

         alphan = alpha2 - sa*(phi1 + tau)/tau
         if (alphan >= top .or. alphan <= bot) then
            alphan = max(p001*top, sqrt(top*bot))
         end if

         ! Get ready for next iteration
         alpha1 = alpha2
         alpha2 = alphan
         phi1 = phi2

      end do

      ! Set NLMS to indicate an optimal step could not be found in 10 trys
      nlms = 12

   end subroutine dodlm

   pure subroutine dacces &
      (n, m, np, nq, ldwe, ld2we, &
       work, lwork, iwork, liwork, &
       access, isodr, &
       jpvt, omega, u, qraux, sd, vcv, wrk1, wrk2, wrk3, wrk4, wrk5, wrk6, &
       nnzw, npp, &
       job, partol, sstol, maxit, taufac, eta, neta, &
       lunrpt, ipr1, ipr2, ipr2f, ipr3, &
       wss, rvar, idf, &
       tau, alpha, niter, nfev, njev, int2, olmavg, &
       rcond, irank, actrs, pnorm, prers, rnorms, istop)
   !! Access or store values in the work arrays.
      ! Routines Called  DIWINF, DWINF
      ! Date Written   860529   (YYMMDD)
      ! Revision Date  920619   (YYMMDD)

      integer, intent(in) :: n
         !! The number of observations.
      integer, intent(in) :: m
         !! The number of columns of data in the explanatory variable.
      integer, intent(in) :: np
         !! The number of function parameters.
      integer, intent(in) :: nq
         !! The number of responses per observation.
      integer, intent(in) :: ldwe
         !! The leading dimension of array `we`.
      integer, intent(in) :: ld2we
         !! The second dimension of array `we`.
      real(wp), intent(inout) :: work(lwork)
         !! The real work space.
      integer, intent(in) :: lwork
         !! The length of vector `work`.
      integer, intent(inout) :: iwork(liwork)
         !! The integer work space.
      integer, intent(in) :: liwork
         !! The length of vector `iwork`.
      logical, intent(in) :: access
         !! The variable designating whether information is to be accessed from the work
         !! arrays (`access = .true.`) or stored in them (`access = .false.`).
      logical, intent(in) :: isodr
         !! The variable designating whether the solution is to be found by ODR (`isodr = .true.`)
         !! or by OLS (`isodr = .false.`).
      integer, intent(out) :: jpvt
         !! The pivot vector.
      integer, intent(out) :: omega
         !! The starting location in array `work` of array `omega`.
      integer, intent(out) :: u
         !! The starting location in array `work` of array `u`.
      integer, intent(out) :: qraux
         !! The starting location in array `work` of array `qraux`.
      integer, intent(out) :: sd
         !! The starting location in array `work` of array `sd`.
      integer, intent(out) :: vcv
         !! The starting location in array `work` of array `vcv`.
      integer, intent(out) :: wrk1
         !! The starting location in array `work` of array `wrk1`.
      integer, intent(out) :: wrk2
         !! The starting location in array `work` of array `wrk2`.
      integer, intent(out) :: wrk3
         !! The starting location in array `work` of array `wrk3`.
      integer, intent(out) :: wrk4
         !! The starting location in array `work` of array `wrk4`.
      integer, intent(out) :: wrk5
         !! The starting location in array `work` of array `wrk5`.
      integer, intent(out) :: wrk6
         !! The starting location in array `work` of array `wrk6`.
      integer, intent(out) :: nnzw
         !! The number of nonzero weighted observations.
      integer, intent(out) :: npp
         !! The number of function parameters actually estimated.
      integer, intent(out) :: job
         !! The variable controlling problem initialization and computational method.
      real(wp), intent(inout) :: partol
         !! The parameter convergence stopping tolerance.
      real(wp), intent(inout) :: sstol
         !! The sum-of-squares convergence stopping tolerance.
      integer, intent(out) :: maxit
         !! The maximum number of iterations allowed.
      real(wp), intent(out) :: taufac
         !! The factor used to compute the initial trust region diameter.
      real(wp), intent(out) :: eta
         !! The relative noise in the function results.
      integer, intent(out) :: neta
         !! The number of accurate digits in the function results.
      integer, intent(out) :: lunrpt
          !! The logical unit number used for computation reports.
      integer, intent(out) :: ipr1
         !! The value of the fourth digit (from the right) of `iprint`, which controls the
         !! initial summary report.
      integer, intent(out) :: ipr2
         !! The value of the third digit (from the right) of `iprint`, which controls the
         !! iteration reports.
      integer, intent(out) :: ipr2f
         !! The value of the second digit (from the right) of `iprint`, which controls the
         !! frequency of the iteration reports.
      integer, intent(out) :: ipr3
         !! The value of the first digit (from the right) of `iprint`, which controls the final
         !! summary report.
      real(wp), intent(inout) :: wss(3)
         !! The sum of the squares of the weighted `epsilons` and `deltas`, the sum of the squares
         !! of the weighted `deltas`, and the sum of the squares of the weighted `epsilons`.
      real(wp), intent(inout) :: rvar
         !! The residual variance, i.e. the standard deviation squared.
      integer, intent(inout) :: idf
         !! The degrees of freedom of the fit, equal to the number of observations with nonzero
         !! weighted derivatives minus the number of parameters being estimated.
      real(wp), intent(inout) :: tau
         !! The trust region diameter.
      real(wp), intent(inout) :: alpha
         !! The Levenberg-Marquardt parameter.
      integer, intent(inout) :: niter
         !! The number of iterations taken.
      integer, intent(inout) :: nfev
         !! The number of function evaluations.
      integer, intent(inout) :: njev
         !! The number of Jacobian evaluations.
      integer, intent(inout) :: int2
         !! The number of internal doubling steps.
      real(wp), intent(inout) :: olmavg
         !! The average number of Levenberg-Marquardt steps per iteration.
      real(wp), intent(inout) :: rcond
         !! The approximate reciprocal condition of `fjacb`.
      integer, intent(inout) :: irank
         !! The rank deficiency of the Jacobian wrt `beta`.
      real(wp), intent(inout) :: actrs
         !! The saved actual relative reduction in the sum-of-squares.
      real(wp), intent(inout) :: pnorm
         !! The norm of the scaled estimated parameters.
      real(wp), intent(inout) :: prers
         !! The saved predicted relative reduction in the sum-of-squares.
      real(wp), intent(inout) :: rnorms
         !! The norm of the saved weighted `epsilons` and `deltas`.
      integer, intent(inout) :: istop
         !! The variable designating whether there are problems computing the function at the
         !! current `beta` and `delta`.

      ! Local scalars
      integer :: actrsi, alphai, betaci, betani, betasi, beta0i, boundi, deltai, deltni, deltsi, &
                 diffi, epsi, epsmai, etai, fjacbi, fjacdi, fni, fsi, idfi, int2i, iprini, &
                 iprint, iranki, istopi, jobi, jpvti, ldtti, liwkmn, loweri, luneri, lunrpi, &
                 lwkmn, maxiti, msgb, msgd, netai, nfevi, niteri, njevi, nnzwi, nppi, nrowi, &
                 ntoli, olmavi, omegai, partli, pnormi, prersi, qrauxi, rcondi, rnorsi, rvari, &
                 sdi, si, ssfi, ssi, sstoli, taufci, taui, ti, tti, ui, upperi, vcvi, we1i, &
                 wrk1i, wrk2i, wrk3i, wrk4i, wrk5i, wrk6i, wrk7i, wssi, wssdei, wssepi, xplusi

      ! Variable Definitions (alphabetically)
      !  ACCESS:  The variable designating whether information is to be accessed from the work
      !           arrays (ACCESS=TRUE) or stored in them (ACCESS=FALSE).
      !  ACTRS:   The saved actual relative reduction in the sum-of-squares.
      !  ACTRSI:  The location in array WORK of variable ACTRS.
      !  ALPHA:   The Levenberg-Marquardt parameter.
      !  ALPHAI:  The location in array WORK of variable ALPHA.
      !  BETACI:  The starting location in array WORK of array BETAC.
      !  BETANI:  The starting location in array WORK of array BETAN.
      !  BETASI:  The starting location in array WORK of array BETAS.
      !  BETA0I:  The starting location in array WORK of array BETA0.
      !  DELTAI:  The starting location in array WORK of array DELTA.
      !  DELTNI:  The starting location in array WORK of array DELTAN.
      !  DELTSI:  The starting location in array WORK of array DELTAS.
      !  DIFFI:   The starting location in array WORK of array DIFF.
      !  EPSI:    The starting location in array WORK of array EPS.
      !  EPSMAI:  The location in array WORK of variable EPSMAC.
      !  ETA:     The relative noise in the function results.
      !  ETAI:    The location in array WORK of variable ETA.
      !  FJACBI:  The starting location in array WORK of array FJACB.
      !  FJACDI:  The starting location in array WORK of array FJACD.
      !  FNI:     The starting location in array WORK of array FN.
      !  FSI:     The starting location in array WORK of array FS.
      !  IDF:     The degrees of freedom of the fit, equal to the number of observations with
      !           nonzero weighted derivatives minus the number of parameters being estimated.
      !  IDFI:    The starting location in array IWORK of variable IDF.
      !  INT2:    The number of internal doubling steps.
      !  INT2I:   The location in array IWORK of variable INT2.
      !  IPR1:    The value of the fourth digit (from the right) of IPRINT, which controls the
      !           initial summary report.
      !  IPR2:    The value of the third digit (from the right) of IPRINT, which controls the
      !           iteration reports.
      !  IPR2F:   The value of the second digit (from the right) of IPRINT, which controls the
      !           frequency of the iteration reports.
      !  IPR3:    The value of the first digit (from the right) of IPRINT, which controls the
      !           final summary report.
      !  IPRINI:  The location in array IWORK of variable IPRINT.
      !  IPRINT:  The print control variable.
      !  IRANK:   The rank deficiency of the Jacobian wrt BETA.
      !  IRANKI:  The location in array IWORK of variable IRANK.
      !  ISODR:   The variable designating whether the solution is to be found by ODR (ISODR=TRUE)
      !           or by OLS (ISODR=FALSE).
      !  ISTOP:   The variable designating whether there are problems computing the function at
      !           the current BETA and DELTA.
      !  ISTOPI:  The location in array IWORK of variable ISTOP.
      !  IWORK:   The integer work space.
      !  JOB:     The variable controling problem initialization and computational method.
      !  JOBI:    The location in array IWORK of variable JOB.
      !  JPVT:    The pivot vector.
      !  JPVTI:   The starting location in array IWORK of variable JPVT.
      !  LDTTI:   The starting location in array IWORK of variable LDTT.
      !  LDWE:    The leading dimension of array WE.
      !  LD2WE:   The second dimension of array WE.
      !  LIWORK:  The length of vector IWORK.
      !  LUNERI:  The location in array IWORK of variable LUNERR.
      !  LUNERR:  The logical unit number used for error messages.
      !  LUNRPI:  The location in array IWORK of variable LUNRPT.
      !  LUNRPT:  The logical unit number used for computation reports.
      !  LWKMN:   The minimum acceptable length of array WORK.
      !  LWORK:   The length of vector WORK.
      !  M:       The number of columns of data in the explanatory variable.
      !  MAXIT:   The maximum number of iterations allowed.
      !  MAXITI:  The location in array IWORK of variable MAXIT.
      !  MSGB:    The starting location in array IWORK of array MSGB.
      !  MSGD:    The starting location in array IWORK of array MSGD.
      !  N:       The number of observations.
      !  NETA:    The number of accurate digits in the function results.
      !  NETAI:   The location in array IWORK of variable NETA.
      !  NFEV:    The number of function evaluations.
      !  NFEVI:   The location in array IWORK of variable NFEV.
      !  NITER:   The number of iterations taken.
      !  NITERI:  The location in array IWORK of variable NITER.
      !  NJEV:    The number of Jacobian evaluations.
      !  NJEVI:   The location in array IWORK of variable NJEV.
      !  NNZW:    The number of nonzero weighted observations.
      !  NNZWI:   The location in array IWORK of variable NNZW.
      !  NP:      The number of function parameters.
      !  NPP:     The number of function parameters actually estimated.
      !  NPPI:    The location in array IWORK of variable NPP.
      !  NQ:      The number of responses per observation.
      !  NROWI:   The location in array IWORK of variable NROW.
      !  NTOLI:   The location in array IWORK of variable NTOL.
      !  OLMAVG:  The average number of Levenberg-Marquardt steps per iteration.
      !  OLMAVI:  The location in array WORK of variable OLMAVG.
      !  OMEGA:   The starting location in array WORK of array OMEGA.
      !  OMEGAI:  The starting location in array WORK of array OMEGA.
      !  PARTLI:  The location in array work of variable PARTOL.
      !  PARTOL:  The parameter convergence stopping tolerance.
      !  PNORM:   The norm of the scaled estimated parameters.
      !  PNORMI:  The location in array WORK of variable PNORM.
      !  PRERS:   The saved predicted relative reduction in the sum-of-squares.
      !  PRERSI:  The location in array WORK of variable PRERS.
      !  QRAUX:   The starting location in array WORK of array QRAUX.
      !  QRAUXI:  The starting location in array WORK of array QRAUX.
      !  RCOND:   The approximate reciprocal condition of FJACB.
      !  RCONDI:  The location in array WORK of variable RCOND.
      !  RESTRT:  The variable designating whether the call is a restart (RESTRT=TRUE) or
      !           not (RESTRT=FALSE).
      !  RNORMS:  The norm of the saved weighted EPSILONS and DELTAS.
      !  RNORSI:  The location in array WORK of variable RNORMS.
      !  RVAR:    The residual variance, i.e. standard deviation squared.
      !  RVARI:   The location in array WORK of variable RVAR.
      !  SCLB:    The scaling values used for BETA.
      !  SCLD:    The scaling values used for DELTA.
      !  SD:      The starting location in array WORK of array SD.
      !  SDI:     The starting location in array WORK of array SD.
      !  SI:      The starting location in array WORK of array S.
      !  SSFI:    The starting location in array WORK of array SSF.
      !  SSI:     The starting location in array WORK of array SS.
      !  SSTOL:   The sum-of-squares convergence stopping tolerance.
      !  SSTOLI:  The location in array WORK of variable SSTOL.
      !  TAU:     The trust region diameter.
      !  TAUFAC:  The factor used to compute the initial trust region diameter.
      !  TAUFCI:  The location in array WORK of variable TAUFAC.
      !  TAUI:    the location in array WORK of variable TAU.
      !  TI:      The starting location in array WORK of array T.
      !  TTI:     The starting location in array WORK of array TT.
      !  U:       The starting location in array WORK of array U.
      !  UI:      The starting location in array WORK of array U.
      !  VCV:     The starting location in array WORK of array VCV.
      !  VCVI:    The starting location in array WORK of array VCV.
      !  WE1I:    The starting location in array WORK of array WE1.
      !  WORK:    The REAL (wp) work space.
      !  WRK1:    The starting location in array WORK of array WRK1.
      !  WRK1I:   The starting location in array WORK of array WRK1.
      !  WRK2:    The starting location in array WORK of array WRK2.
      !  WRK2I:   The starting location in array WORK of array WRK2.
      !  WRK3:    The starting location in array WORK of array wrk3.
      !  WRK3I:   The starting location in array WORK of array wrk3.
      !  WRK4:    The starting location in array WORK of array wrk4.
      !  WRK4I:   The starting location in array WORK of array wrk4.
      !  WRK5:    The starting location in array WORK of array wrk5.
      !  WRK5I:   The starting location in array WORK of array wrk5.
      !  WRK6:    The starting location in array WORK of array wrk6.
      !  WRK6I:   The starting location in array WORK of array wrk6.
      !  WRK7I:   The starting location in array WORK of array wrk7.
      !  WSS:     The sum of the squares of the weighted EPSILONS and DELTAS.
      !  WSSI:    The starting location in array WORK of variable WSS(1).
      !  WSSDEI:  The starting location in array WORK of variable WSS(2).
      !  WSSEPI:  The starting location in array WORK of variable WSS(3).
      !  XPLUSI:  The starting location in array WORK of array XPLUSD.

      ! Find starting locations within integer workspace
      call diwinf(m, np, nq, &
                  msgb, msgd, jpvti, istopi, &
                  nnzwi, nppi, idfi, &
                  jobi, iprini, luneri, lunrpi, &
                  nrowi, ntoli, netai, &
                  maxiti, niteri, nfevi, njevi, int2i, iranki, ldtti, &
                  boundi, &
                  liwkmn)

      ! Find starting locations within REAL work space
      call dwinf(n, m, np, nq, ldwe, ld2we, isodr, &
                 deltai, epsi, xplusi, fni, sdi, vcvi, &
                 rvari, wssi, wssdei, wssepi, rcondi, etai, &
                 olmavi, taui, alphai, actrsi, pnormi, rnorsi, prersi, &
                 partli, sstoli, taufci, epsmai, &
                 beta0i, betaci, betasi, betani, si, ssi, ssfi, qrauxi, ui, &
                 fsi, fjacbi, we1i, diffi, &
                 deltsi, deltni, ti, tti, omegai, fjacdi, &
                 wrk1i, wrk2i, wrk3i, wrk4i, wrk5i, wrk6i, wrk7i, &
                 loweri, upperi, &
                 lwkmn)

      if (access) then
         ! Set starting locations for work vectors
         jpvt = jpvti
         omega = omegai
         qraux = qrauxi
         sd = sdi
         vcv = vcvi
         u = ui
         wrk1 = wrk1i
         wrk2 = wrk2i
         wrk3 = wrk3i
         wrk4 = wrk4i
         wrk5 = wrk5i
         wrk6 = wrk6i
         ! Access values from the work vectors
         actrs = work(actrsi)
         alpha = work(alphai)
         eta = work(etai)
         olmavg = work(olmavi)
         partol = work(partli)
         pnorm = work(pnormi)
         prers = work(prersi)
         rcond = work(rcondi)
         wss(1) = work(wssi)
         wss(2) = work(wssdei)
         wss(3) = work(wssepi)
         rvar = work(rvari)
         rnorms = work(rnorsi)
         sstol = work(sstoli)
         tau = work(taui)
         taufac = work(taufci)

         neta = iwork(netai)
         irank = iwork(iranki)
         job = iwork(jobi)
         lunrpt = iwork(lunrpi)
         maxit = iwork(maxiti)
         nfev = iwork(nfevi)
         niter = iwork(niteri)
         njev = iwork(njevi)
         nnzw = iwork(nnzwi)
         npp = iwork(nppi)
         idf = iwork(idfi)
         int2 = iwork(int2i)

         ! Set up print control variables
         iprint = iwork(iprini)
         ipr1 = mod(iprint, 10000)/1000
         ipr2 = mod(iprint, 1000)/100
         ipr2f = mod(iprint, 100)/10
         ipr3 = mod(iprint, 10)

      else
         ! Store values into the work vectors
         work(actrsi) = actrs
         work(alphai) = alpha
         work(olmavi) = olmavg
         work(partli) = partol
         work(pnormi) = pnorm
         work(prersi) = prers
         work(rcondi) = rcond
         work(wssi) = wss(1)
         work(wssdei) = wss(2)
         work(wssepi) = wss(3)
         work(rvari) = rvar
         work(rnorsi) = rnorms
         work(sstoli) = sstol
         work(taui) = tau

         iwork(iranki) = irank
         iwork(istopi) = istop
         iwork(nfevi) = nfev
         iwork(niteri) = niter
         iwork(njevi) = njev
         iwork(idfi) = idf
         iwork(int2i) = int2
      end if

   end subroutine dacces

   real(wp) pure function derstep(itype, k, betak, ssf, stpb, neta) result(derstepr)
   !! Compute step size for center and forward difference calculations.
      ! Routines Called  DHSTEP
      ! Date Written   20040616   (YYYYMMDD)
      ! Revision Date  20040616   (YYYYMMDD)

      use odrpack_kinds, only: zero, one

      integer, intent(in) :: itype
         !! The finite difference method being used, where: `itype = 0` indicates forward
         !! finite differences, and `itype = 1` indicates central finite differences.
      integer, intent(in) :: k
         !! Index into `beta` where `betak` resides.
      real(wp), intent(in) :: betak
         !! The `k`-th function parameter.
      real(wp), intent(in) :: ssf(k)
         !! The scale used for the `beta`s.
      real(wp), intent(in) :: stpb(k)
         !! The relative step used for computing finite difference derivatives with respect
         !! to `beta`.
      integer, intent(in) :: neta
         !! Number of good digits in the function results.

      ! Local scalars
      real(wp) :: typj

      ! Variable definitions (alphabetically)
      !  BETAK:   The K-th function parameter.
      !  ITYPE:   0 - calc foward difference step, 1 - calc center difference step.
      !  K:       Index into beta where BETAK resides.
      !  NETA:    Number of good digits in the function results.
      !  SSF:     The scale used for the BETA'S.
      !  STPB:    The relative step used for computing finite difference derivatives with
      !           respect to BETA.
      !  TYPJ:    The typical size of the J-th unkonwn BETA.

      if (betak == zero) then
         if (ssf(1) < zero) then
            typj = one/abs(ssf(1))
         else
            typj = one/ssf(k)
         end if
      else
         typj = abs(betak)
      end if
      derstepr = sign(one, betak)*typj*dhstep(itype, neta, 1, k, stpb, 1)

   end function derstep

   pure subroutine desubi(n, m, wd, ldwd, ld2wd, alpha, tt, ldtt, i, e)
   !! Compute `e = wd + alpha*tt**2`.
      ! Routines Called (NONE)
      ! Date Written   860529   (YYMMDD)
      ! Revision Date  920304   (YYMMDD)

      use odrpack_kinds, only: zero

      integer, intent(in) :: n
         !! The number of observations.
      integer, intent(in) :: m
         !! The number of columns of data in the independent variable.
      real(wp), intent(in) :: wd(ldwd, ld2wd, m)
         !! The squared `delta` weights.
      integer, intent(in) :: ldwd
         !! The leading dimension of array `wd`.
      integer, intent(in) :: ld2wd
         !! The second dimension of array `wd`.
      real(wp), intent(in) :: alpha
         !! The Levenberg-Marquardt parameter.
      real(wp), intent(in) :: tt(ldtt, m)
         !! The scaling values used for `delta`.
      integer, intent(in) :: ldtt
         !! The leading dimension of array `tt`.
      integer, intent(in) :: i
         !! An indexing variable.
      real(wp), intent(out) :: e(m, m)
         !! The value of the array `e = wd + alpha*tt**2`.

      ! Local scalars
      integer :: j, j1, j2

      ! Variable Definitions (alphabetically)
      !  ALPHA:  The Levenberg-Marquardt parameter.
      !  E:      The value of the array E = WD + ALPHA*TT**2
      !  I:      An indexing variable.
      !  J:      An indexing variable.
      !  J1:     An indexing variable.
      !  J2:     An indexing variable.
      !  LDWD:   The leading dimension of array WD.
      !  LD2WD:  The second dimension of array WD.
      !  M:      The number of columns of data in the independent variable.
      !  N:      The number of observations.
      !  NP:     The number of responses per observation.
      !  TT:     The scaling values used for DELTA.
      !  WD:     The squared DELTA weights, D**2.

      ! N.B. the locations of WD and TT accessed depend on the value of the first element of
      ! each array and the leading dimensions of the multiply subscripted arrays.

      if (n == 0 .or. m == 0) return

      if (wd(1, 1, 1) >= zero) then
         if (ldwd >= n) then
            ! The elements of WD have been individually specified
            if (ld2wd == 1) then
               ! The arrays stored in WD are diagonal
               e = zero
               do j = 1, m
                  e(j, j) = wd(i, 1, j)
               end do
            else
               ! The arrays stored in WD are full positive semidefinite matrices
               do j1 = 1, m
                  do j2 = 1, m
                     e(j1, j2) = wd(i, j1, j2)
                  end do
               end do
            end if
            !
            if (tt(1, 1) > zero) then
               if (ldtt >= n) then
                  do j = 1, m
                     e(j, j) = e(j, j) + alpha*tt(i, j)**2
                  end do
               else
                  do j = 1, m
                     e(j, j) = e(j, j) + alpha*tt(1, j)**2
                  end do
               end if
            else
               do j = 1, m
                  e(j, j) = e(j, j) + alpha*tt(1, 1)**2
               end do
            end if
         else
            ! WD is an M by M matrix
            if (ld2wd == 1) then
               ! The array stored in WD is diagonal
               e = zero
               do j = 1, m
                  e(j, j) = wd(1, 1, j)
               end do
            else
               ! The array stored in WD is a full positive semidefinite matrix
               do j1 = 1, m
                  do j2 = 1, m
                     e(j1, j2) = wd(1, j1, j2)
                  end do
               end do
            end if

            if (tt(1, 1) > zero) then
               if (ldtt >= n) then
                  do j = 1, m
                     e(j, j) = e(j, j) + alpha*tt(i, j)**2
                  end do
               else
                  do j = 1, m
                     e(j, j) = e(j, j) + alpha*tt(1, j)**2
                  end do
               end if
            else
               do j = 1, m
                  e(j, j) = e(j, j) + alpha*tt(1, 1)**2
               end do
            end if
         end if
      else
         ! WD is a diagonal matrix with elements ABS(WD(1,1,1))
         e = zero
         if (tt(1, 1) > zero) then
            if (ldtt >= n) then
               do j = 1, m
                  e(j, j) = abs(wd(1, 1, 1)) + alpha*tt(i, j)**2
               end do
            else
               do j = 1, m
                  e(j, j) = abs(wd(1, 1, 1)) + alpha*tt(1, j)**2
               end do
            end if
         else
            do j = 1, m
               e(j, j) = abs(wd(1, 1, 1)) + alpha*tt(1, 1)**2
            end do
         end if
      end if

   end subroutine desubi

   subroutine detaf &
      (fcn, &
       n, m, np, nq, &
       xplusd, beta, epsmac, nrow, &
       partmp, pv0, &
       ifixb, ifixx, ldifx, &
       istop, nfev, eta, neta, &
       wrk1, wrk2, wrk6, wrk7, &
       info, &
       lower, upper)
   !! Compute noise and number of good digits in function results.
      ! Adapted from STARPAC subroutine ETAFUN.
      ! Routines Called  FCN
      ! Date Written   860529   (YYMMDD)
      ! Revision Date  920619   (YYMMDD)

      use odrpack_kinds, only: zero, one, two, hundred

      procedure(fcn_t) :: fcn
         !! The user-supplied subroutine for evaluating the model.
      integer, intent(in) :: n
         !! The number of observations.
      integer, intent(in) :: m
         !! The number of columns of data in the explanatory variable.
      integer, intent(in) :: np
         !! The number of function parameters.
      integer, intent(in) :: nq
         !! The number of responses per observation.
      real(wp), intent(in) :: xplusd(n, m)
         !! The values of `x + delta`.
      real(wp), intent(in) :: beta(np)
         !! The function parameters.
      real(wp), intent(in) :: epsmac
         !! The value of machine precision.
      integer, intent(in) :: nrow
         !! The row number at which the derivative is to be checked.
      real(wp), intent(out) :: partmp(np)
         !! The model parameters.
      real(wp), intent(in) :: pv0(n, nq)
         !! The original predicted values.
      integer, intent(in) :: ifixb(np)
         !! The values designating whether the elements of `beta` are fixed at their input values or not.
      integer, intent(in) :: ifixx(ldifx, m)
         !! The values designating whether the elements of `x` are fixed at their input values or not.
      integer, intent(in) :: ldifx
         !! The leading dimension of array `ifixx`.
      integer, intent(out) :: istop
         !! The variable designating whether there are problems computing the function at the
         !! current `beta` and `delta`.
      integer, intent(inout) :: nfev
         !! The number of function evaluations.
      real(wp), intent(out) :: eta
         !! The noise in the model results.
      integer, intent(out) :: neta
         !! The number of accurate digits in the model results.
      real(wp), intent(out) :: wrk1(n, m, nq)
         !! A work array of `(n, m, nq)` elements.
      real(wp), intent(out) :: wrk2(n, nq)
         !! A work array of `(n, nq)` elements.
      real(wp), intent(out) :: wrk6(n, np, nq)
         !! A work array of `(n, np, nq)` elements.
      real(wp), intent(out) :: wrk7(-2:2, nq)
         !! A work array of `(5, nq)` elements.
      integer, intent(out) :: info
         !! The variable indicating the status of the computation.
      real(wp), intent(in) :: lower(np)
         !! The lower bound of `beta`.
      real(wp), intent(in) :: upper(np)
         !! The upper bound of `beta`.

      ! Local scalars
      real(wp), parameter :: p1 = 0.1_wp, p2 = 0.2_wp, p5 = 0.5_wp
      real(wp) :: a, b, fac, shift, stp
      integer :: j, k, l, sbk

      ! Local arrays
      real(wp) :: parpts(-2:2, np)

      ! Variable Definitions (ALPHABETICALLY)
      !  A:       Parameters of the local fit.
      !  B:       Parameters of the local fit.
      !  BETA:    The function parameters.
      !  EPSMAC:  The value of machine precision.
      !  ETA:     The noise in the model results.
      !  FAC:     A factor used in the computations.
      !  FCN:     The user supplied subroutine for evaluating the model.
      !  IFIXB:   The values designating whether the elements of BETA are fixed at their input
      !           values or not.
      !  IFIXX:   The values designating whether the elements of X are fixed at their input values
      !           or not.
      !  ISTOP:   The variable designating whether there are problems computing the function at the
      !           current BETA and DELTA.
      !  J:       An index variable.
      !  K:       An index variable.
      !  L:       AN INDEX VARIABLE.
      !  LDIFX:   The leading dimension of array IFIXX.
      !  LOWER:   The lower bound of BETA.
      !  M:       The number of columns of data in the explanatory variable.
      !  N:       The number of observations.
      !  NETA:    The number of accurate digits in the model results.
      !  NFEV:    The number of function evaluations.
      !  NP:      The number of function parameters.
      !  NQ:      The number of responses per observation.
      !  NROW:    The row number at which the derivative is to be checked.
      !  P1:      The value 0.1E0_wp.
      !  P2:      The value 0.2E0_wp.
      !  P5:      The value 0.5E0_wp.
      !  PARPTS:  The points that PARTMP will take on during FCN evaluations.
      !  PARTMP:  The model parameters.
      !  PV0:     The original predicted values.
      !  SHIFT:   When PARPTS cross the parameter bounds they are shifted by SHIFT.
      !  SBK:     The sign of BETA(K).
      !  STP:     A small value used to perturb the parameters.
      !  UPPER:   The upper bound of BETA.
      !  WRK1:    A work array of (N BY M BY NQ) elements.
      !  WRK2:    A work array of (N BY NQ) elements.
      !  WRK6:    A work array of (N BY NP BY NQ) elements.
      !  WRK7:    A work array of (5 BY NQ) elements.
      !  XPLUSD:  The values of X + DELTA.

      stp = hundred*epsmac
      eta = epsmac

      ! Create points to use in calculating FCN for ETA and NETA
      do j = -2, 2
         if (j == 0) then
            parpts(0, :) = beta(:)
         else
            do k = 1, np
               if (ifixb(1) < 0) then
                  parpts(j, k) = beta(k) + j*stp*beta(k)
               elseif (ifixb(k) /= 0) then
                  parpts(j, k) = beta(k) + j*stp*beta(k)
               else
                  parpts(j, k) = beta(k)
               end if
            end do
         end if
      end do

      ! Adjust the points used in calculating FCN to uphold the boundary constraints
      do k = 1, np
         sbk = int(sign(one, parpts(2, k) - parpts(-2, k)))
         if (parpts(sbk*2, k) > upper(k)) then
            shift = upper(k) - parpts(sbk*2, k)
            parpts(sbk*2, k) = upper(k)
            do j = -sbk*2, sbk*1, sbk
               parpts(j, k) = parpts(j, k) + shift
            end do
            if (parpts(-sbk*2, k) < lower(k)) then
               info = 90010
               return
            end if
         end if
         if (parpts(-sbk*2, k) < lower(k)) then
            shift = lower(k) - parpts(-sbk*2, k)
            parpts(-sbk*2, k) = lower(k)
            do j = -sbk*1, sbk*2, sbk
               parpts(j, k) = parpts(j, k) + shift
            end do
            if (parpts(sbk*2, k) > upper(k)) then
               info = 90010
               return
            end if
         end if
      end do

      ! Evaluate FCN for all points in PARPTS
      do j = -2, 2
         if (all(parpts(j, :) == beta(:))) then
            do l = 1, nq
               wrk7(j, l) = pv0(nrow, l)
            end do
         else
            partmp(:) = parpts(j, :)
            istop = 0
            call fcn(n, m, np, nq, &
                     n, m, np, &
                     partmp(:), xplusd, &
                     ifixb, ifixx, ldifx, &
                     003, wrk2, wrk6, wrk1, istop)
            if (istop /= 0) then
               return
            else
               nfev = nfev + 1
            end if
            do l = 1, nq
               wrk7(j, l) = wrk2(nrow, l)
            end do
         end if
      end do

      ! Calculate ETA and NETA
      do l = 1, nq
         a = zero
         b = zero
         do j = -2, 2
            a = a + wrk7(j, l)
            b = b + j*wrk7(j, l)
         end do
         a = p2*a
         b = p1*b
         if ((wrk7(0, l) /= zero) .and. (abs(wrk7(1, l) + wrk7(-1, l)) > hundred*epsmac)) then
            fac = one/abs(wrk7(0, l))
         else
            fac = one
         end if
         do j = -2, 2
            wrk7(j, l) = abs((wrk7(j, l) - (a + j*b))*fac)
            eta = max(wrk7(j, l), eta)
         end do
      end do
      neta = int(max(two, p5 - log10(eta)))

   end subroutine detaf

   subroutine devjac &
      (fcn, &
       anajac, cdjac, &
       n, m, np, nq, &
       betac, beta, stpb, &
       ifixb, ifixx, ldifx, &
       x, ldx, delta, xplusd, stpd, ldstpd, &
       ssf, tt, ldtt, neta, fn, &
       stp, wrk1, wrk2, wrk3, wrk6, tempret, &
       fjacb, isodr, fjacd, we1, ldwe, ld2we, &
       njev, nfev, istop, info, &
       lower, upper)
   !! Compute the weighted Jacobians wrt `beta` and `delta`.
      ! Routines Called  FCN, DDOT, DIFIX, DJACCD, DJACFD, DWGHT, DUNPAC
      ! Date Written   860529   (YYMMDD)
      ! Revision Date  920304   (YYMMDD)

      use odrpack_kinds, only: zero

      procedure(fcn_t) :: fcn
         !! The user-supplied subroutine for evaluating the model.
      logical, intent(in) :: anajac
         !! The variable designating whether the Jacobians are computed by finite differences
         !! (`anajac = .false.`) or not (`anajac = .true.`).
      logical, intent(in) :: cdjac
         !! The variable designating whether the Jacobians are computed by central differences
         !! (`cdjac = .true.`) or by forward differences (`cdjac = .false.`).
      integer, intent(in) :: n
         !! The number of observations.
      integer, intent(in) :: m
         !! The number of columns of data in the independent variable.
      integer, intent(in) :: np
         !! The number of function parameters.
      integer, intent(in) :: nq
         !! The number of responses per observation.
      real(wp), intent(in) :: betac(np)
         !! The current estimated values of the unfixed `beta`s.
      real(wp), intent(out) :: beta(np)
         !! The function parameters.
      real(wp), intent(in) :: stpb(np)
         !! The relative step used for computing finite difference derivatives with respect to `beta`.
      integer, intent(in) :: ifixb(np)
         !! The values designating whether the elements of `beta` are fixed at their input values or not.
      integer, intent(in) :: ifixx(ldifx, m)
         !! The values designating whether the elements of `delta` are fixed at their input values or not.
      integer, intent(in) :: ldifx
         !! The leading dimension of array `ifixx`.
      real(wp), intent(in) :: x(ldx, m)
         !! The independent variable.
      integer, intent(in) :: ldx
         !! The leading dimension of array `x`.
      real(wp), intent(in) :: delta(n, m)
         !! The estimated values of `delta`.
      real(wp), intent(out) :: xplusd(n, m)
         !! The values of `x + delta`.
      real(wp), intent(in) :: stpd(ldstpd, m)
         !! The relative step used for computing finite difference derivatives with respect to `delta`.
      integer, intent(in) :: ldstpd
         !! The leading dimension of array `stpd`.
      real(wp), intent(in) :: ssf(np)
         !! The scale used for the `beta`s.
      real(wp), intent(in) :: tt(ldtt, m)
         !! The scaling values used for `delta`.
      integer, intent(in) :: ldtt
         !! The leading dimension of array `tt`.
      integer, intent(in) :: neta
         !! The number of accurate digits in the function results.
      real(wp), intent(in) :: fn(n, nq)
         !! The predicted values of the function at the current point.
      real(wp), intent(out) :: stp(n)
         !! The step used for computing finite difference derivatives with respect to `delta`.
      real(wp), intent(out) :: wrk1(n, m, nq)
         !! A work array of `(n, m, nq)` elements.
      real(wp), intent(out) :: wrk2(n, nq)
         !! A work array of `(n, nq)` elements.
      real(wp), intent(out) :: wrk3(np)
         !! A work array of `(np)` elements.
      real(wp), intent(out) :: wrk6(n, np, nq)
         !! A work array of `(n, np, nq)` elements.
      real(wp), intent(inout) :: tempret(:, :)
         !! Temporary work array for holding return values before copying to a lower rank array.
      real(wp), intent(out) :: fjacb(n, np, nq)
         !! The Jacobian with respect to `beta`.
      logical, intent(in) :: isodr
         !! The variable designating whether the solution is by ODR (`isodr = .true.`) or
         !! by OLS (`isodr = .false.`).
      real(wp), intent(out) :: fjacd(n, m, nq)
         !! The Jacobian with respect to `delta`.
      real(wp), intent(in) :: we1(ldwe, ld2we, nq)
         !! The square roots of the `epsilon` weights in array `we`.
      integer, intent(in) :: ldwe
         !! The leading dimension of arrays `we` and `we1`.
      integer, intent(in) :: ld2we
         !! The second dimension of arrays `we` and `we1`.
      integer, intent(inout) :: njev
         !! The number of Jacobian evaluations.
      integer, intent(inout) :: nfev
         !! The number of function evaluations.
      integer, intent(out) :: istop
         !! The variable designating that the user wishes the computations stopped.
      integer, intent(out) :: info
         !! The variable designating why the computations were stopped.
      real(wp), intent(in) :: lower(np)
         !! The lower bound of `beta`.
      real(wp), intent(in) :: upper(np)
         !! The upper bound of `beta`.

      ! Local scalars
      integer :: ideval, j, k, k1, l
      logical :: ferror

      ! External BLAS/LINPACK procedures
      real(wp), external :: ddot

      ! Variable Definitions (alphabetically)
      !  ANAJAC:  The variable designating whether the Jacobians are computed by finite differences
      !          (ANAJAC=FALSE) or not (ANAJAC=TRUE).
      !  BETA:    The function parameters.
      !  BETAC:   The current estimated values of the unfixed BETA's.
      !  CDJAC:   The variable designating whether the Jacobians are computed by central differences
      !           (CDJAC=TRUE) or by forward differences (CDJAC=FALSE).
      !  DELTA:   The estimated values of DELTA.
      !  FERROR:  The variable designating whether ODRPACK95 detected nonzero values in array DELTA
      !           in the OLS case, and thus whether the user may have overwritten important information
      !           by computing FJACD in the OLS case.
      !  FCN:     The user-supplied subroutine for evaluating the model.
      !  FJACB:   The Jacobian with respect to BETA.
      !  FJACD:   The Jacobian with respect to DELTA.
      !  FN:      The predicted values of the function at the current point.
      !  IDEVAL:  The variable designating what computations are to be performed by user-supplied
      !           subroutine FCN.
      !  IFIXB:   The values designating whether the elements of BETA are fixed at their input values or not.
      !  IFIXX:   The values designating whether the elements of DELTA are fixed at their input values or not.
      !  INFO:    The variable designating why the computations were stopped.
      !  ISTOP:   The variable designating that the user wishes the computations stopped.
      !  ISODR:   The variable designating whether the solution is by ODR (ISODR=TRUE) or OLS (ISODR=FALSE).
      !  J:       An indexing variable.
      !  K:       An indexing variable.
      !  K1:      An indexing variable.
      !  L:       An indexing variable.
      !  LDIFX:   The leading dimension of array IFIXX.
      !  LDSTPD:  The leading dimension of array STPD.
      !  LDTT:    The leading dimension of array TT.
      !  LDWE:    The leading dimension of arrays WE and WE1.
      !  LDX:     The leading dimension of array X.
      !  LD2WE:   The second dimension of arrays WE and WE1.
      !  M:       The number of columns of data in the independent variable.
      !  N:       The number of observations.
      !  NETA:    The number of accurate digits in the function results.
      !  NFEV:    The number of function evaluations.
      !  NJEV:    The number of Jacobian evaluations.
      !  NP:      The number of function parameters.
      !  NQ:      The number of responses per observation.
      !  SSF:     The scale used for the BETA's.
      !  STP:     The step used for computing finite difference derivatives with respect to DELTA.
      !  STPB:    The relative step used for computing finite difference derivatives with respect to BETA.
      !  STPD:    The relative step used for computing finite difference derivatives with respect to DELTA.
      !  TT:      The scaling values used for DELTA.
      !  WE1:     The square roots of the EPSILON weights in array WE.
      !  WRK1:    A work array of (N by M by NQ) elements.
      !  WRK2:    A work array of (N by NQ) elements.
      !  WRK3:    A work array of (NP) elements.
      !  WRK6:    A work array of (N BY NP BY NQ) elements.
      !  X:       The independent variable.
      !  XPLUSD:  The values of X + DELTA.

      ! Insert current unfixed BETA estimates into BETA
      call dunpac(np, betac, beta, ifixb)

      ! Compute XPLUSD = X + DELTA
      xplusd = x(1:n, :) + delta

      ! Compute the Jacobian wrt the estimated BETAS (FJACB) and the Jacobian wrt DELTA (FJACD)
      istop = 0
      if (isodr) then
         ideval = 110
      else
         ideval = 010
      end if
      if (anajac) then
         call fcn(n, m, np, nq, &
                  n, m, np, &
                  beta, xplusd, &
                  ifixb, ifixx, ldifx, &
                  ideval, wrk2, fjacb, fjacd, &
                  istop)
         if (istop /= 0) then
            return
         else
            njev = njev + 1
         end if
         ! Make sure fixed elements of FJACD are zero
         if (isodr) then
            do l = 1, nq
               call difix(n, m, ifixx, ldifx, fjacd(1, 1, l), n, fjacd(1, 1, l), n)
            end do
         end if
      elseif (cdjac) then
         call djaccd(fcn, &
                     n, m, np, nq, &
                     beta, x, ldx, delta, xplusd, ifixb, ifixx, ldifx, &
                     stpb, stpd, ldstpd, &
                     ssf, tt, ldtt, neta, fn, stp, wrk1, wrk2, wrk3, wrk6, &
                     fjacb, isodr, fjacd, nfev, istop, info, &
                     lower, upper)
      else
         call djacfd(fcn, &
                     n, m, np, nq, &
                     beta, x, ldx, delta, xplusd, ifixb, ifixx, ldifx, &
                     stpb, stpd, ldstpd, &
                     ssf, tt, ldtt, neta, fn, stp, wrk1, wrk2, wrk3, wrk6, &
                     fjacb, isodr, fjacd, nfev, istop, info, &
                     lower, upper)
      end if
      if (istop < 0 .or. info >= 10000) then
         return
      elseif (.not. isodr) then
         ! Try to detect whether the user has computed JFACD within FCN in the OLS case
         ferror = ddot(n*m, delta, 1, delta, 1) /= zero
         if (ferror) then
            info = 50300
            return
         end if
      end if

      ! Weight the Jacobian wrt the estimated BETAS
      if (ifixb(1) < 0) then
         do k = 1, np
            call dwght(n, nq, we1, ldwe, ld2we, fjacb(1:n, k, 1:nq), tempret(1:n, 1:nq))
            fjacb(1:n, k, 1:nq) = tempret(1:n, 1:nq)
         end do
      else
         k1 = 0
         do k = 1, np
            if (ifixb(k) >= 1) then
               k1 = k1 + 1
               call dwght(n, nq, we1, ldwe, ld2we, fjacb(1:n, k, 1:nq), tempret(1:n, 1:nq))
               fjacb(1:n, k1, 1:nq) = tempret(1:n, 1:nq)
            end if
         end do
      end if

      ! Weight the Jacobian's wrt DELTA as appropriate
      if (isodr) then
         do j = 1, m
            call dwght(n, nq, we1, ldwe, ld2we, fjacd(1:n, j, 1:nq), tempret(1:n, 1:nq))
            fjacd(1:n, j, 1:nq) = tempret(1:n, 1:nq)
         end do
      end if

   end subroutine devjac

   subroutine dfctr(oksemi, a, lda, n, info)
   !! Factor the positive (semi)definite matrix `a` using a modified Cholesky factorization.
      ! Adapted from LINPACK subroutine DPOFA.
      ! Routines Called  DDOT
      ! Date Written   910706   (YYMMDD)
      ! Revision Date  920619   (YYMMDD)
      ! References  Dongarra J.J., Bunch J.R., Moler C.B., Stewart G.W., *LINPACK Users Guide*, SIAM, 1979.

      use odrpack_kinds, only: zero, ten

      logical, intent(in) :: oksemi
         !! The indicating whether the factored array can be positive semidefinite
         !! (`oksemi = .true.`) or whether it must be found to be positive definite
         !! (`oksemi = .false.`).
      real(wp), intent(inout) :: a(lda, n)
         !! The array to be factored. Upon return, `a` contains the upper triangular matrix
         !! `r` so that `a = trans(r)*r` where the strict lower triangle is set to zero.
         !! If `info /= 0`, the factorization is not complete.
      integer, intent(in) :: lda
         !! The leading dimension of array `a`.
      integer, intent(in) :: n
         !! The number of rows and columns of data in array `a`.
      integer, intent(out) :: info
         !! An indicator variable, where if:
         !!   `info = 0` then factorization was completed.
         !!   `info = k` signals an error condition. The leading minor of order `k` is not
         !!   positive (semi)definite.

      ! Local scalars
      real(wp) :: xi, s, t
      integer j, k

      ! External BLAS/LINPACK procedures
      real(wp), external :: ddot

      ! Variable Definitions (alphabetically)
      !  A:       The array to be factored.  Upon return, A contains the upper triangular matrix
      !           R so that A = trans(R)*R where the strict lower triangle is set to zero.
      !           If INFO /= 0, the factorization is not complete.
      !  I:       An indexing variable.
      !  INFO:    An idicator variable, where if
      !           INFO = 0  then factorization was completed
      !           INFO = K  signals an error condition.  The leading minor of order  K  is not
      !           positive (semi)definite.
      !  J:       An indexing variable.
      !  LDA:     The leading dimension of array A.
      !  N:       The number of rows and columns of data in array A.
      !  OKSEMI:  The indicating whether the factored array can be positive semidefinite
      !           (OKSEMI=TRUE) or whether it must be found to be positive definite (OKSEMI=FALSE).
      !  XI:      A value used to test for non positive semidefiniteness.

      ! Set relative tolerance for detecting non positive semidefiniteness.
      xi = -ten*epsilon(zero)

      ! Compute factorization, storing in upper triangular portion of A
      do j = 1, n
         info = j
         s = zero
         do k = 1, j - 1
            if (a(k, k) == zero) then
               t = zero
            else
               t = a(k, j) - ddot(k - 1, a(1, k), 1, a(1, j), 1)
               t = t/a(k, k)
            end if
            a(k, j) = t
            s = s + t*t
         end do
         s = a(j, j) - s

         ! ...Exit
         if (a(j, j) < zero .or. s < xi*abs(a(j, j))) then
            return
         elseif (.not. oksemi .and. s <= zero) then
            return
         elseif (s <= zero) then
            a(j, j) = zero
         else
            a(j, j) = sqrt(s)
         end if
      end do
      info = 0

      ! Zero out lower portion of A
      do j = 2, n
         do k = 1, j - 1
            a(j, k) = zero
         end do
      end do

   end subroutine dfctr

   subroutine dfctrw &
      (n, m, nq, npp, &
       isodr, &
       we, ldwe, ld2we, wd, ldwd, ld2wd, &
       wrk0, wrk4, &
       we1, nnzw, info)
   !! Check input parameters, indicating errors found using nonzero values of argument `info` as
   !! described in the ODRPACK95 reference guide.
      ! Routines Called  DFCTR
      ! Date Written   860529   (YYMMDD)
      ! Revision Date  920619   (YYMMDD)

      use odrpack_kinds, only: zero

      integer, intent(in) :: n
         !! The number of observations.
      integer, intent(in) :: m
         !! The number of columns of data in the explanatory variable.
      integer, intent(in) :: nq
         !! The number of responses per observation.
      integer, intent(in) :: npp
         !! The number of function parameters being estimated.
      logical, intent(in) :: isodr
         !! The variable designating whether the solution is by ODR (`isodr = .true`) or
         !! by OLS (`isodr = .false`).
      real(wp), intent(in) :: we(ldwe, ld2we, nq)
         !! The (squared) `epsilon` weights.
      integer, intent(in) :: ldwe
         !! The leading dimension of array `we`.
      integer, intent(in) :: ld2we
         !! The second dimension of array `we`.
      real(wp), intent(in) :: wd(ldwd, ld2wd, m)
         !! The (squared) `delta` weights.
      integer, intent(in) :: ldwd
      ! The leading dimension of array `wd`.
      integer, intent(in) :: ld2wd
         !! The second dimension of array `wd`.
      real(wp), intent(out) :: wrk0(nq, nq)
         !! A work array of `(nq, nq)` elements.
      real(wp), intent(out) :: wrk4(m, m)
         !! A work array of `(m, m)` elements.
      real(wp), intent(out) :: we1(ldwe, ld2we, nq)
         !! The factored `epsilon` weights, such that `trans(we1)*we1 = we`.
      integer, intent(out) :: nnzw
         !! The number of nonzero weighted observations.
      integer, intent(out) :: info
         !! The variable designating why the computations were stopped.

      ! Local scalars
      integer :: i, inf, j, j1, j2, l, l1, l2
      logical :: notzro, exited

      ! Variable Definitions (alphabetically)
      !  I:       An indexing variable.
      !  INFO:    The variable designating why the computations were stopped.
      !  ISODR:   The variable designating whether the solution is by ODR (ISODR=TRUE) or by
      !           OLS (ISODR=FALSE).
      !  J:       An indexing variable.
      !  J1:      An indexing variable.
      !  J2:      An indexing variable.
      !  L:       An indexing variable.
      !  L1:      An indexing variable.
      !  L2:      An indexing variable.
      !  LAST:    The last row of the array to be accessed.
      !  LDWD:    The leading dimension of array WD.
      !  LDWE:    The leading dimension of array WE.
      !  LD2WD:   The second dimension of array WD.
      !  LD2WE:   The second dimension of array WE.
      !  M:       The number of columns of data in the explanatory variable.
      !  N:       The number of observations.
      !  NNZW:    The number of nonzero weighted observations.
      !  NOTZRO:  The variable designating whether a given component of the weight array WE
      !           contains a nonzero element (NOTZRO=FALSE) or not (NOTZRO=TRUE).
      !  NPP:     The number of function parameters being estimated.
      !  NQ:      The number of responses per observations.
      !  WE:      The (squared) EPSILON weights.
      !  WE1:     The factored EPSILON weights, S.T. trans(WE1)*WE1 = WE.
      !  WD:      The (squared) DELTA weights.
      !  WRK0:    A work array of (NQ BY NQ) elements.
      !  WRK4:    A work array of (M BY M) elements.

      ! Check EPSILON weights, and store factorization in WE1

      exited = .false.

      if (we(1, 1, 1) < zero) then
         ! WE contains a scalar
         we1(1, 1, 1) = -sqrt(abs(we(1, 1, 1)))
         nnzw = n
      else
         nnzw = 0
         if (ldwe == 1) then
            if (ld2we == 1) then
               ! WE contains a diagonal matrix
               do l = 1, nq
                  if (we(1, 1, l) > zero) then
                     nnzw = n
                     we1(1, 1, l) = sqrt(we(1, 1, l))
                  elseif (we(1, 1, l) < zero) then
                     info = 30010
                     exited = .true.
                     exit
                  end if
               end do
            else
               ! WE contains a full NQ by NQ semidefinite matrix
               do l1 = 1, nq
                  do l2 = l1, nq
                     wrk0(l1, l2) = we(1, l1, l2)
                  end do
               end do
               call dfctr(.true., wrk0, nq, nq, inf)
               if (inf /= 0) then
                  info = 30010
                  exited = .true.
               else
                  do l1 = 1, nq
                     do l2 = 1, nq
                        we1(1, l1, l2) = wrk0(l1, l2)
                     end do
                     if (we1(1, l1, l1) /= zero) then
                        nnzw = n
                     end if
                  end do
               end if
            end if
         else
            if (ld2we == 1) then
               ! WE contains an array of  diagonal matrix
               do i = 1, n
                  notzro = .false.
                  do l = 1, nq
                     if (we(i, 1, l) > zero) then
                        notzro = .true.
                        we1(i, 1, l) = sqrt(we(i, 1, l))
                     elseif (we(i, 1, l) < zero) then
                        info = 30010
                        exited = .true.
                        exit
                     end if
                  end do
                  if (exited) exit
                  if (notzro) then
                     nnzw = nnzw + 1
                  end if
               end do
            else
               ! WE contains an array of full NQ by NQ semidefinite matrices
               do i = 1, n
                  do l1 = 1, nq
                     do l2 = l1, nq
                        wrk0(l1, l2) = we(i, l1, l2)
                     end do
                  end do
                  call dfctr(.true., wrk0, nq, nq, inf)
                  if (inf /= 0) then
                     info = 30010
                     exited = .true.
                     exit
                  else
                     notzro = .false.
                     do l1 = 1, nq
                        do l2 = 1, nq
                           we1(i, l1, l2) = wrk0(l1, l2)
                        end do
                        if (we1(i, l1, l1) /= zero) then
                           notzro = .true.
                        end if
                     end do
                  end if
                  if (notzro) then
                     nnzw = nnzw + 1
                  end if
               end do
            end if
         end if
      end if

      ! Check for a sufficient number of nonzero EPSILON weights
      if (.not. exited) then
         if (nnzw < npp) info = 30020
      end if

      ! Check DELTA weights
      if (.not. isodr .or. wd(1, 1, 1) < zero) then
         ! Problem is not ODR, or WD contains a scalar
         return
      else
         if (ldwd == 1) then
            if (ld2wd == 1) then
               ! WD contains a diagonal matrix
               do j = 1, m
                  if (wd(1, 1, j) <= zero) then
                     info = max(30001, info + 1)
                     return
                  end if
               end do
            else
               ! WD contains a full M by M positive definite matrix
               do j1 = 1, m
                  do j2 = j1, m
                     wrk4(j1, j2) = wd(1, j1, j2)
                  end do
               end do
               call dfctr(.false., wrk4, m, m, inf)
               if (inf /= 0) then
                  info = max(30001, info + 1)
                  return
               end if
            end if
         else
            if (ld2wd == 1) then
               ! WD contains an array of diagonal matrices
               do i = 1, n
                  do j = 1, m
                     if (wd(i, 1, j) <= zero) then
                        info = max(30001, info + 1)
                        return
                     end if
                  end do
               end do
            else
               ! WD contains an array of full M by M positive definite matrices
               do i = 1, n
                  do j1 = 1, m
                     do j2 = j1, m
                        wrk4(j1, j2) = wd(i, j1, j2)
                     end do
                  end do
                  call dfctr(.false., wrk4, m, m, inf)
                  if (inf /= 0) then
                     info = max(30001, info + 1)
                     return
                  end if
               end do
            end if
         end if
      end if

   end subroutine dfctrw

   pure subroutine dflags(job, restrt, initd, dovcv, redoj, anajac, cdjac, chkjac, isodr, implct)
   !! Set flags indicating conditions specified by `job`.
      ! Routines Called  (NONE)
      ! Date Written   860529   (YYMMDD)
      ! Revision Date  920304   (YYMMDD)

      integer, intent(in) :: job
         !! The variable controlling problem initialization and computational method.
      logical, intent(out) :: restrt
         !! The variable designating whether the call is a restart (`restrt = .true.`)
         !! or not (`restrt = .false.`).
      logical, intent(out) :: initd
         !! The variable designating whether `delta` is to be initialized to zero (`initd = .true.`)
         !! or to the first `n` by `m` elements of array `work` (`initd = .false.`).
      logical, intent(out) :: dovcv
         !! The variable designating whether the covariance matrix is to be computed
         !! (`dovcv = .true.`) or not (`dovcv = .false.`).
      logical, intent(out) :: redoj
         !! The variable designating whether the Jacobian matrix is to be recomputed for the
         !! computation of the covariance matrix (`redoj = .true.`) or not (`redoj = .false.`).
      logical, intent(out) :: anajac
         !! The variable designating whether the Jacobians are computed by finite differences
         !! (`anajac = .false.`) or not (`anajac = .true.`).
      logical, intent(out) :: cdjac
         !! The variable designating whether the Jacobians are computed by central differences
         !! (`cdjac = .true.`) or by forward differences (`cdjac = .false.`).
      logical, intent(out) :: chkjac
         !! The variable designating whether the user-supplied Jacobians are to be checked
         !! (`chkjac = .true.`) or not (`chkjac = .false.`).
      logical, intent(out) :: isodr
         !! The variable designating whether the solution is by ODR (`isodr = .true.`)
         !! or by OLS (`isodr = .false.`).
      logical, intent(out) :: implct
         !! The variable designating whether the solution is by implicit ODR (`implct = .true.`)
         !! or explicit ODR (`implct = .false.`).

      ! Local scalars
      integer :: j

      ! Variable Definitions (alphabetically)
      !  ANAJAC:  The variable designating whether the Jacobians are computed by finite differences
      !           (ANAJAC=FALSE) or not (ANAJAC=TRUE).
      !  CDJAC:   The variable designating whether the Jacobians are computed by central differences
      !           (CDJAC=TRUE) or by forward differences (CDJAC=FALSE).
      !  CHKJAC:  The variable designating whether the user-supplied Jacobians are to be checked
      !           (CHKJAC=TRUE) or not (CHKJAC=FALSE).
      !  DOVCV:   The variable designating whether the covariance matrix is to be computed
      !           (DOVCV=TRUE) or not (DOVCV=FALSE).
      !  IMPLCT:  The variable designating whether the solution is by implicit ODR (IMPLCT=TRUE)
      !           or explicit ODR (IMPLCT=FALSE).
      !  INITD:   The variable designating whether DELTA is to be initialized to zero (INITD=TRUE)
      !           or to the first N by M elements of array WORK (INITD=FALSE).
      !  ISODR:   The variable designating whether the solution is by ODR (ISODR=TRUE) or
      !           by OLS (ISODR=FALSE).
      !  J:       The value of a specific digit of JOB.
      !  JOB:     The variable controling problem initialization and computational method.
      !  REDOJ:   The variable designating whether the Jacobian matrix is to be recomputed for the
      !           computation of the covariance matrix (REDOJ=TRUE) or not (REDOJ=FALSE).
      !  RESTRT:  The variable designating whether the call is a restart (RESTRT=TRUE) or
      !           not (RESTRT=FALSE).

      if (job >= 0) then

         restrt = job >= 10000
         initd = mod(job, 10000)/1000 == 0
         j = mod(job, 1000)/100

         if (j == 0) then
            dovcv = .true.
            redoj = .true.
         elseif (j == 1) then
            dovcv = .true.
            redoj = .false.
         else
            dovcv = .false.
            redoj = .false.
         end if

         j = mod(job, 100)/10
         if (j == 0) then
            anajac = .false.
            cdjac = .false.
            chkjac = .false.
         elseif (j == 1) then
            anajac = .false.
            cdjac = .true.
            chkjac = .false.
         elseif (j == 2) then
            anajac = .true.
            cdjac = .false.
            chkjac = .true.
         else
            anajac = .true.
            cdjac = .false.
            chkjac = .false.
         end if

         j = mod(job, 10)
         if (j == 0) then
            isodr = .true.
            implct = .false.
         elseif (j == 1) then
            isodr = .true.
            implct = .true.
         else
            isodr = .false.
            implct = .false.
         end if

      else

         restrt = .false.
         initd = .true.
         dovcv = .true.
         redoj = .true.
         anajac = .false.
         cdjac = .false.
         chkjac = .false.
         isodr = .true.
         implct = .false.

      end if

   end subroutine dflags

   real(wp) pure function dhstep(itype, neta, i, j, stp, ldstp) result(dhstepr)
   !! Set relative step size for finite difference derivatives.
      ! Routines Called  (NONE)
      ! Date Written   860529   (YYMMDD)
      ! Revision Date  920304   (YYMMDD)

      use odrpack_kinds, only: zero, two, three, ten

      integer, intent(in) :: itype
         !! The finite difference method being used, where: `itype = 0` indicates forward
         !! finite differences, and `itype = 1` indicates central finite differences.
      integer, intent(in) :: neta
         !! The number of good digits in the function results.
      integer, intent(in) :: i
         !! An identifier for selecting user-supplied step sizes.
      integer, intent(in) :: j
         !! An identifier for selecting user-supplied step sizes.
      real(wp), intent(in) :: stp(ldstp, j)
         !! The step size for the finite difference derivative.
      integer, intent(in) :: ldstp
         !! The leading dimension of array `stp`.

      ! Variable Definitions (alphabetically)
      !  I:       An identifier for selecting user supplied step sizes.
      !  ITYPE:   The finite difference method being used, where
      !           ITYPE = 0 indicates forward finite differences, and
      !           ITYPE = 1 indicates central finite differences.
      !  J:       An identifier for selecting user supplied step sizes.
      !  LDSTP:   The leading dimension of array STP.
      !  NETA:    The number of good digits in the function results.
      !  STP:     The step size for the finite difference derivative.

      if (stp(1, 1) <= zero) then
         if (itype == 0) then
            ! Use default forward finite difference step size
            dhstepr = ten**(-abs(neta)/two - two)
         else
            ! Use default central finite difference step size
            dhstepr = ten**(-abs(neta)/three)
         end if

      elseif (ldstp == 1) then
         dhstepr = stp(1, j)
      else
         dhstepr = stp(i, j)
      end if

   end function dhstep

   pure subroutine difix(n, m, ifix, ldifix, t, ldt, tfix, ldtfix)
   !! Set elements of `t` to zero according to `ifix`.
      ! Routines Called  (NONE)
      ! Date Written   910612   (YYMMDD)
      ! Revision Date  920304   (YYMMDD)

      use odrpack_kinds, only: zero

      integer, intent(in) :: n
         !! The number of rows of data in the array.
      integer, intent(in) :: m
         !! The number of columns of data in the array.
      integer, intent(in) :: ifix(ldifix, m)
         !! The array designating whether an element of `t` is to be set to zero.
      integer, intent(in) :: ldifix
         !! The leading dimension of array `ifix`.
      real(wp), intent(in) :: t(ldt, m)
         !! The array being set to zero according to the elements of `ifix`.
      integer, intent(in) :: ldt
         !! The leading dimension of array `t`.
      real(wp), intent(out) :: tfix(ldtfix, m)
         !! The resulting array.
      integer, intent(in) :: ldtfix
         !! The leading dimension of array `tfix`.

      ! Local scalars
      integer :: i, j

      ! Variable Definitions (alphabetically)
      !  I:       An indexing variable.
      !  IFIX:    The array designating whether an element of T is to be set to zero.
      !  J:       an indexing variable.
      !  LDT:     The leading dimension of array T.
      !  LDIFIX:  The leading dimension of array IFIX.
      !  LDTFIX:  The leading dimension of array TFIX.
      !  M:       The number of columns of data in the array.
      !  N:       The number of rows of data in the array.
      !  T:       The array being set to zero according to the elements of IFIX.
      !  TFIX:    The resulting array.

      if (n == 0 .or. m == 0) return

      if (ifix(1, 1) >= zero) then
         if (ldifix >= n) then
            do j = 1, m
               do i = 1, n
                  if (ifix(i, j) == 0) then
                     tfix(i, j) = zero
                  else
                     tfix(i, j) = t(i, j)
                  end if
               end do
            end do
         else
            do j = 1, m
               if (ifix(1, j) == 0) then
                  do i = 1, n
                     tfix(i, j) = zero
                  end do
               else
                  do i = 1, n
                     tfix(i, j) = t(i, j)
                  end do
               end if
            end do
         end if
      end if

   end subroutine difix

   pure subroutine diwinf &
      (m, np, nq, &
       msgbi, msgdi, ifix2i, istopi, &
       nnzwi, nppi, idfi, &
       jobi, iprini, luneri, lunrpi, &
       nrowi, ntoli, netai, &
       maxiti, niteri, nfevi, njevi, int2i, iranki, ldtti, &
       boundi, &
       liwkmn)
   !! Set storage locations within integer work space.
      ! Routines Called  (NONE)
      ! Date Written   860529   (YYMMDD)
      ! Revision Date  920304   (YYMMDD)

      integer, intent(in) :: m
         !! The number of columns of data in the independent variable.
      integer, intent(in) :: np
         !! The number of function parameters.
      integer, intent(in) :: nq
         !! The number of responses per observation.
      integer, intent(out) :: msgbi
         !! The starting location in array `iwork` of array `msgb`.
      integer, intent(out) :: msgdi
         !! The starting location in array `iwork` of array `msgd`.
      integer, intent(out) :: ifix2i
         !! The starting location in array `iwork` of array `ifix2`.
      integer, intent(out) :: istopi
         !! The location in array `iwork` of variable `istop`.
      integer, intent(out) :: nnzwi
         !! The location in array `iwork` of variable `nnzw`.
      integer, intent(out) :: nppi
         !! The location in array `iwork` of variable `npp`.
      integer, intent(out) :: idfi
         !! The location in array `iwork` of variable `idf`.
      integer, intent(out) :: jobi
         !! The location in array `iwork` of variable `job`.
      integer, intent(out) :: iprini
         !! The location in array `iwork` of variable `iprint`.
      integer, intent(out) :: luneri
         !! The location in array `iwork` of variable `lunerr`.
      integer, intent(out) :: lunrpi
         !! The location in array `iwork` of variable `lunrpt`.
      integer, intent(out) :: nrowi
         !! The location in array `iwork` of variable `nrow`.
      integer, intent(out) :: ntoli
         !! The location in array `iwork` of variable `ntol`.
      integer, intent(out) :: netai
         !! The location in array `iwork` of variable `neta`.
      integer, intent(out) :: maxiti
         !! The location in array `iwork` of variable `maxit`.
      integer, intent(out) :: niteri
         !! The location in array `iwork` of variable `niter`.
      integer, intent(out) :: nfevi
         !! The location in array `iwork` of variable `nfev`.
      integer, intent(out) :: njevi
         !! The location in array `iwork` of variable `njev`.
      integer, intent(out) :: int2i
         !! The location in array `iwork` of variable `int2`.
      integer, intent(out) :: iranki
         !! The location in array `iwork` of variable `irank`.
      integer, intent(out) :: ldtti
         !! The location in array `iwork` of variable `ldtt`.
      integer, intent(out) :: boundi
         !! The location in array `iwork` of variable `bound`.
      integer, intent(out) :: liwkmn
         !! The minimum acceptable length of array `iwork`.

      ! Variable Definitions (alphabetically)
      !  IDFI:    The location in array IWORK of variable IDF.
      !  IFIX2I:  The starting location in array IWORK of array IFIX2.
      !  INT2I:   The location in array IWORK of variable INT2.
      !  IPRINI:  The location in array IWORK of variable IPRINT.
      !  IRANKI:  The location in array IWORK of variable IRANK.
      !  ISTOPI:  The location in array IWORK of variable ISTOP.
      !  JOBI:    The location in array IWORK of variable JOB.
      !  LDTTI:   The location in array IWORK of variable LDTT.
      !  LIWKMN:  The minimum acceptable length of array IWORK.
      !  LUNERI:  The location in array IWORK of variable LUNERR.
      !  LUNRPI:  The location in array IWORK of variable LUNRPT.
      !  M:       The number of columns of data in the independent variable.
      !  MAXITI:  The location in array iwork of variable MAXIT.
      !  MSGBI:   The starting location in array IWORK of array MSGB.
      !  MSGDI:   The starting location in array IWORK of array MSGD.
      !  NETAI:   The location in array IWORK of variable NETA.
      !  NFEVI:   The location in array IWORK of variable NFEV.
      !  NITERI:  The location in array IWORK of variabel NITER.
      !  NJEVI:   The location in array IWORK of variable NJEV.
      !  NNZWI:   The location in array IWORK of variable NNZW.
      !  NP:      The number of function parameters.
      !  NPPI:    The location in array IWORK of variable NPP.
      !  NQ:      The number of responses per observation.
      !  NROWI:   The location in array IWORK of variable NROW.
      !  NTOLI:   The location in array IWORK of variable NTOL.

      if (np >= 1 .and. m >= 1) then
         msgbi = 1
         msgdi = msgbi + nq*np + 1
         ifix2i = msgdi + nq*m + 1
         istopi = ifix2i + np
         nnzwi = istopi + 1
         nppi = nnzwi + 1
         idfi = nppi + 1
         jobi = idfi + 1
         iprini = jobi + 1
         luneri = iprini + 1
         lunrpi = luneri + 1
         nrowi = lunrpi + 1
         ntoli = nrowi + 1
         netai = ntoli + 1
         maxiti = netai + 1
         niteri = maxiti + 1
         nfevi = niteri + 1
         njevi = nfevi + 1
         int2i = njevi + 1
         iranki = int2i + 1
         ldtti = iranki + 1
         boundi = ldtti + 1
         liwkmn = boundi + np - 1
      else
         msgbi = 1
         msgdi = 1
         ifix2i = 1
         istopi = 1
         nnzwi = 1
         nppi = 1
         idfi = 1
         jobi = 1
         iprini = 1
         luneri = 1
         lunrpi = 1
         nrowi = 1
         ntoli = 1
         netai = 1
         maxiti = 1
         niteri = 1
         nfevi = 1
         njevi = 1
         int2i = 1
         iranki = 1
         ldtti = 1
         boundi = 1
         liwkmn = 1
      end if

   end subroutine diwinf

   subroutine diniwk &
      (n, m, np, work, lwork, iwork, liwork, &
       x, ldx, ifixx, ldifx, scld, ldscld, &
       beta, sclb, &
       sstol, partol, maxit, taufac, &
       job, iprint, lunerr, lunrpt, &
       lower, upper, &
       epsmai, sstoli, partli, maxiti, taufci, &
       jobi, iprini, luneri, lunrpi, &
       ssfi, tti, ldtti, deltai, &
       loweri, upperi, boundi)
   !! Initialize work vectors as necessary.
      ! Routines Called  DFLAGS, DSCLB, DSCLD, DCOPY
      ! Date Written   860529   (YYMMDD)
      ! Revision Date  920304   (YYMMDD)

      use odrpack_kinds, only: zero, one, two, three

      integer, intent(in) :: n
         !! The number of observations.
      integer, intent(in) :: m
         !! The number of columns of data in the independent variable.
      integer, intent(in) :: np
         !! The number of function parameters.
      real(wp), intent(out) :: work(lwork)
         !! The real work space.
      integer, intent(in) :: lwork
         !! The length of vector `work`.
      integer, intent(out) :: iwork(liwork)
         !! The integer work space.
      integer, intent(in) :: liwork
         !! The length of vector `iwork`.
      real(wp), intent(in) :: x(ldx, m)
         !! The independent variable.
      integer, intent(in) :: ldx
         !! The leading dimension of array `x`.
      integer, intent(in) :: ifixx(ldifx, m)
         !! The values designating whether the elements of `x` are fixed at their input values or not.
      integer, intent(in) :: ldifx
         !! The leading dimension of array `ifixx`.
      real(wp), intent(in) :: scld(ldscld, m)
         !! The scaling values for `delta`.
      integer, intent(in) :: ldscld
         !! The leading dimension of array `scld`.
      real(wp), intent(in) :: beta(np)
         !! The function parameters.
      real(wp), intent(in) :: sclb(np)
         !! The scaling values for `beta`.
      real(wp), intent(in) :: sstol
         !! The sum-of-squares convergence stopping criteria.
      real(wp), intent(in) :: partol
         !! The parameter convergence stopping criteria.
      integer, intent(in) :: maxit
         !! The maximum number of iterations allowed.
      real(wp), intent(in) :: taufac
         !! The factor used to compute the initial trust region diameter.
      integer, intent(in) :: job
         !! The variable controlling problem initialization and computational method.
      integer, intent(in) :: iprint
         !! The print control variable.
      integer, intent(in) :: lunerr
         !! The logical unit number used for error messages.
      integer, intent(in) :: lunrpt
         !! The logical unit number used for computation reports.
      real(wp), intent(in) :: lower(np)
         !! The lower bounds for the function parameters.
      real(wp), intent(in) :: upper(np)
         !! The upper bounds for the function parameters.
      integer, intent(in) :: epsmai
         !! The location in array `work` of variable `epsmac`.
      integer, intent(in) :: sstoli
         !! The location in array `work` of variable `sstol`.
      integer, intent(in) :: partli
         !! The location in array `work` of variable `partol`.
      integer, intent(in) :: maxiti
         !! The location in array `iwork` of variable `maxit`.
      integer, intent(in) :: taufci
         !! The location in array `work` of variable `taufac`.
      integer, intent(in) :: jobi
         !! The location in array `iwork` of variable `job`.
      integer, intent(in) :: iprini
         !! The location in array `iwork` of variable `iprint`.
      integer, intent(in) :: luneri
         !! The location in array `iwork` of variable `lunerr`.
      integer, intent(in) :: lunrpi
         !! The location in array `iwork` of variable `lunrpt`.
      integer, intent(in) :: ssfi
         !! The starting location in array `work` of array `ssf`.
      integer, intent(in) :: tti
         !! The starting location in array `work` of the array `tt`.
      integer, intent(in) :: ldtti
         !! The leading dimension of array `tt`.
      integer, intent(in) :: deltai
         !! The starting location in array `work` of array `delta`.
      integer, intent(in) :: loweri
         !! The starting location in array `iwork` of array `lower`.
      integer, intent(in) :: upperi
         !! The starting location in array `iwork` of array `upper`.
      integer, intent(in) :: boundi
         !! The location in array `iwork` of variable `bound`.

      ! Local scalars
      integer :: i, j, istart
      logical :: anajac, cdjac, chkjac, dovcv, implct, initd, isodr, redoj, restrt

      ! External BLAS/LINPACK procedures
      external :: dcopy

      ! Variable Definitions (alphabetically)
      !  ANAJAC:  The variable designating whether the Jacobians are computed by finite differences
      !           (ANAJAC=FALSE) or not (ANAJAC=TRUE).
      !  BETA:    The function parameters.
      !  CDJAC:   The variable designating whether the Jacobians are computed by central differences
      !           (CDJAC=TRUE) or by forward differences (CDJAC=FALSE).
      !  CHKJAC:  The variable designating whether the user-supplied Jacobians are to be checked
      !           (CHKJAC=TRUE) or not (CHKJAC=FALSE).
      !  DELTAI:  The starting location in array WORK of array DELTA.
      !  DOVCV:   The variable designating whether the covariance matrix is to be computed (DOVCV=TRUE)
      !           or not (DOVCV=FALSE).
      !  EPSMAI:  The location in array WORK of variable EPSMAC.
      !  I:       An indexing variable.
      !  IFIXX:   The values designating whether the elements of X are fixed at their input values or not.
      !  IMPLCT:  The variable designating whether the solution is by implicit ODR (IMPLCT=TRUE) or
      !           explicit ODR (IMPLCT=FALSE).
      !  INITD:   The variable designating whether DELTA is to be initialized to zero (INITD=TRUE) or
      !           to the values in the first N by M elements of array WORK (INITD=FALSE).
      !  IPRINI:  The location in array IWORK of variable IPRINT.
      !  IPRINT:  The print control variable.
      !  ISODR:   The variable designating whether the solution is by ODR (ISODR=TRUE) or by OLS
      !           ISODR=FALSE).
      !  IWORK:   The integer work space.
      !  J:       An indexing variable.
      !  JOB:     The variable controling problem initialization and computational method.
      !  JOBI:    The location in array IWORK of variable JOB.
      !  LDIFX:   The leading dimension of array IFIXX.
      !  LDSCLD:  The leading dimension of array SCLD.
      !  LDTTI:   The leading dimension of array TT.
      !  LDX:     The leading dimension of array X.
      !  LIWORK:  The length of vector IWORK.
      !  LUNERI:  The location in array IWORK of variable LUNERR.
      !  LUNERR:  The logical unit number used for error messages.
      !  LUNRPI:  The location in array iwork of variable LUNRPT.
      !  LUNRPT:  The logical unit number used for computation reports.
      !  LWORK:   The length of vector WORK.
      !  M:       The number of columns of data in the independent variable.
      !  MAXIT:   The maximum number of iterations allowed.
      !  MAXITI:  The location in array IWORK of variable MAXIT.
      !  N:       The number of observations.
      !  NP:      The number of function parameters.
      !  PARTLI:  The location in array work of variable partol.
      !  PARTOL:  The parameter convergence stopping criteria.
      !  REDOJ:   The variable designating whether the Jacobian matrix is to be recomputed for the
      !           computation of the covariance matrix (REDOJ=TRUE) or not (REDOJ=FALSE).
      !  RESTRT:  The variable designating whether the call is a restart (RESTRT=TRUE) or not
      !           (RESTRT=FALSE).
      !  SCLB:    The scaling values for BETA.
      !  SCLD:    The scaling values for DELTA.
      !  SSFI:    The starting location in array WORK of array SSF.
      !  SSTOL:   The sum-of-squares convergence stopping criteria.
      !  SSTOLI:  The location in array WORK of variable SSTOL.
      !  TAUFAC:  The factor used to compute the initial trust region diameter.
      !  TAUFCI:  The location in array WORK of variable TAUFAC.
      !  TTI:     The starting location in array WORK of the ARRAY TT.
      !  WORK:    The REAL (wp) work space.
      !  X:       The independent variable.

      call dflags(job, restrt, initd, dovcv, redoj, anajac, cdjac, chkjac, isodr, implct)

      ! Store value of machine precision in work vector
      work(epsmai) = epsilon(zero)

      ! Set tolerance for stopping criteria based on the change in the parameters (see also
      ! subprogram DODCNT)
      if (partol < zero) then
         work(partli) = work(epsmai)**(two/three)
      else
         work(partli) = min(partol, one)
      end if

      ! Set tolerance for stopping criteria based on the change in the sum of squares of the
      ! weighted observational errors
      if (sstol < zero) then
         work(sstoli) = sqrt(work(epsmai))
      else
         work(sstoli) = min(sstol, one)
      end if

      ! Set factor for computing trust region diameter at first iteration
      if (taufac <= zero) then
         work(taufci) = one
      else
         work(taufci) = min(taufac, one)
      end if

      ! Set maximum number of iterations
      if (maxit < 0) then
         iwork(maxiti) = 50
      else
         iwork(maxiti) = maxit
      end if

      ! Store problem initialization and computational method control variable
      if (job <= 0) then
         iwork(jobi) = 0
      else
         iwork(jobi) = job
      end if

      ! Set print control
      if (iprint < 0) then
         iwork(iprini) = 2001
      else
         iwork(iprini) = iprint
      end if

      ! Set logical unit number for error messages
      iwork(luneri) = lunerr

      ! Set logical unit number for computation reports
      iwork(lunrpi) = lunrpt

      ! Compute scaling for BETA's and DELTA's
      if (sclb(1) <= zero) then
         call dsclb(np, beta, work(ssfi))
      else
         call dcopy(np, sclb, 1, work(ssfi), 1)
      end if
      if (isodr) then
         if (scld(1, 1) <= zero) then
            iwork(ldtti) = n
            call dscld(n, m, x, ldx, work(tti), iwork(ldtti))
         else
            if (ldscld == 1) then
               iwork(ldtti) = 1
               call dcopy(m, scld(1, 1), 1, work(tti), 1)
            else
               iwork(ldtti) = n
               do j = 1, m
                  call dcopy(n, scld(1, j), 1, work(tti + (j - 1)*iwork(ldtti)), 1)
               end do
            end if
         end if
      end if

      ! Initialize DELTA's as necessary
      if (isodr) then
         if (initd) then
            !call dzero( n, m, work( deltai), n)
            work(deltai:deltai + (n*m - 1)) = zero
         else
            if (ifixx(1, 1) >= 0) then
               if (ldifx == 1) then
                  do j = 1, m
                     if (ifixx(1, j) == 0) then
                        istart = deltai + (j - 1)*n
                        work(istart:istart + (n - 1)) = zero
                     end if
                  end do
               else
                  do j = 1, m
                     do i = 1, n
                        if (ifixx(i, j) == 0) then
                           work(deltai - 1 + i + (j - 1)*n) = zero
                        end if
                     end do
                  end do
               end if
            end if
         end if
      else
         !call dzero( n, m, work( deltai), n)
         work(deltai:deltai + (n*m - 1)) = zero
      end if

      ! Copy bounds into WORK
      work(loweri:loweri + np - 1) = lower(1:np)
      work(upperi:upperi + np - 1) = upper(1:np)

      ! Initialize parameters on bounds in IWORK
      iwork(boundi:boundi + np - 1) = 0

   end subroutine diniwk

   subroutine djaccd &
      (fcn, &
       n, m, np, nq, &
       beta, x, ldx, delta, xplusd, ifixb, ifixx, ldifx, &
       stpb, stpd, ldstpd, &
       ssf, tt, ldtt, neta, fn, stp, wrk1, wrk2, wrk3, wrk6, &
       fjacb, isodr, fjacd, nfev, istop, info, &
       lower, upper)
   !! Compute central difference approximations to the Jacobian wrt the estimated `beta`s and
   !! wrt the `delta`s.
      ! Routines Called  FCN, DHSTEP
      ! Date Written   860529   (YYMMDD)
      ! Revision Date  920619   (YYMMDD)

      use odrpack_kinds, only: zero, one

      procedure(fcn_t) :: fcn
         !! The user supplied subroutine for evaluating the model.
      integer, intent(in) :: n
         !! The number of observations.
      integer, intent(in) :: m
         !! The number of columns of data in the explanatory variable.
      integer, intent(in) :: np
         !! The number of function parameters.
      integer, intent(in) :: nq
         !! The number of responses per observation.
      real(wp), intent(inout) :: beta(np)
         !! The function parameters.
      real(wp), intent(in) :: x(ldx, m)
         !! The explanatory variable.
      integer, intent(in) :: ldx
         !! The leading dimension of array `x`.
      real(wp), intent(in) :: delta(n, m)
         !! The estimated errors in the explanatory variables.
      real(wp), intent(inout) :: xplusd(n, m)
         !! The values of `x + delta`.
      integer, intent(in) :: ifixb(np)
         !! The values designating whether the elements of `beta` are fixed at their input values or not.
      integer, intent(in) :: ifixx(ldifx, m)
         !! The values designating whether the elements of `x` are fixed at their input values or not.
      integer, intent(in) :: ldifx
         !! The leading dimension of array `ifixx`.
      real(wp), intent(in) :: stpb(np)
         !! The relative step used for computing finite difference derivatives with respect to each `beta`.
      real(wp), intent(in) :: stpd(ldstpd, m)
         !! The relative step used for computing finite difference derivatives with respect to each `delta`.
      integer, intent(in) :: ldstpd
         !! The leading dimension of array `stpd`.
      real(wp), intent(in) :: ssf(np)
         !! The scaling values used for `beta`.
      real(wp), intent(in) :: tt(ldtt, m)
         !! The scaling values used for `delta`.
      integer, intent(in) :: ldtt
         !! The leading dimension of array `tt`.
      integer, intent(in) :: neta
         !! The number of good digits in the function results.
      real(wp), intent(in) :: fn(n, nq)
         !! The new predicted values from the function. Used when parameter is on a boundary.
      real(wp), intent(out) :: stp(n)
         !! The step used for computing finite difference derivatives with respect to each `delta`.
      real(wp), intent(out) :: wrk1(n, m, nq)
         !! A work array of `(n, m, nq)` elements.
      real(wp), intent(out) :: wrk2(n, nq)
         !! A work array of `(n, nq)` elements.
      real(wp), intent(out) :: wrk3(np)
         !! A work array of `(np)` elements.
      real(wp), intent(out) :: wrk6(n, np, nq)
         !! A work array of `(n, np, nq)` elements.
      real(wp), intent(out) :: fjacb(n, np, nq)
         !! The Jacobian with respect to `beta`.
      logical, intent(in) :: isodr
         !! The variable designating whether the solution is by ODR (`isodr = .true.`) or
         !! by OLS (`isodr = .false.`).
      real(wp), intent(out) :: fjacd(n, m, nq)
         !! The Jacobian with respect to `delta`.
      integer, intent(inout) :: nfev
         !! The number of function evaluations.
      integer, intent(out) :: istop
         !! The variable designating whether there are problems computing the function at the
         !! current `beta` and `delta`.
      integer, intent(out) :: info
         !! The variable designating why the computations were stopped.
      real(wp), intent(in) :: lower(np)
         !! The lower bound on `beta`.
      real(wp), intent(in) :: upper(np)
         !! The upper bound on `beta`.

      ! Local scalars
      real(wp) :: betak, typj
      integer :: i, j, k, l
      logical :: doit, setzro

      ! Variable Definitions (alphabetically)
      !  BETA:    The function parameters.
      !  BETAK:   The K-th function parameter.
      !  DELTA:   The estimated errors in the explanatory variables.
      !  DOIT:    The variable designating whether the derivative wrt a given BETA or DELTA needs
      !           to be computed (DOIT=TRUE) or not (DOIT=FALSE).
      !  FCN:     The user supplied subroutine for evaluating the model.
      !  FJACB:   The Jacobian with respect to BETA.
      !  FJACD:   The Jacobian with respect to DELTA.
      !  FN:      The new predicted values from the function. Used when parameter is on a boundary.
      !  I:       An indexing variable.
      !  IFIXB:   The values designating whether the elements of BETA are fixed at their input
      !           values or not.
      !  IFIXX:   The values designating whether the elements of X are fixed at their input values
      !           or not.
      !  INFO:    The variable designating why the computations were stopped.
      !  ISODR:   The variable designating whether the solution is by ODR (ISODR=TRUE) or by
      !           OLS (ISODR=FALSE).
      !  ISTOP:   The variable designating whether there are problems computing the function at the
      !           current BETA and DELTA.
      !  J:       An indexing variable.
      !  K:       An indexing variable.
      !  L:       An indexing variable.
      !  LDIFX:   The leading dimension of array IFIXX.
      !  LDSTPD:  The leading dimension of array STPD.
      !  LDTT:    The leading dimension of array TT.
      !  LDX:     The leading dimension of array X.
      !  LOWER:   The lower bound on BETA.
      !  M:       The number of columns of data in the explanatory variable.
      !  N:       The number of observations.
      !  NETA:    The number of good digits in the function results.
      !  NFEV:    The number of function evaluations.
      !  NP:      The number of function parameters.
      !  SETZRO:  The variable designating whether the derivative wrt some DELTA needs to be set to
      !           zero (SETZRO=TRUE) or not (SETZRO=FALSE).
      !  SSF:     The scaling values used for BETA.
      !  STP:     The step used for computing finite difference derivatives with respect to each DELTA.
      !  STPB:    the relative step used for computing finite difference derivatives with respect
      !           to each BETA.
      !  STPD:    The relative step used for computing finite difference derivatives with respect
      !           to each DELTA.
      !  TT:      The scaling values used for DELTA.
      !  TYPJ:    The typical size of the J-th unknown BETA or DELTA.
      !  UPPER:   The upper bound on BETA.
      !  X:       The explanatory variable.
      !  XPLUSD:  The values of X + DELTA.
      !  WRK1:    A work array of (N BY M BY NQ) elements.
      !  WRK2:    A work array of (N BY NQ) elements.
      !  WRK3:    A work array of (NP) elements.
      !  WRK6:    A WORK ARRAY OF (N BY NP BY NQ) elements.

      ! Compute the Jacobian wrt the estimated BETAS
      do k = 1, np
         if (ifixb(1) >= 0) then
            if (ifixb(k) == 0) then
               doit = .false.
            else
               doit = .true.
            end if
         else
            doit = .true.
         end if
         if (.not. doit) then
            fjacb(1:n, k, 1:nq) = zero
         else
            betak = beta(k)
            wrk3(k) = betak + derstep(1, k, betak, ssf, stpb, neta)
            wrk3(k) = wrk3(k) - betak

            beta(k) = betak + wrk3(k)
            if (beta(k) > upper(k)) then
               beta(k) = upper(k)
            elseif (beta(k) < lower(k)) then
               beta(k) = lower(k)
            end if
            if (beta(k) - 2*wrk3(k) < lower(k)) then
               beta(k) = lower(k) + 2*wrk3(k)
            elseif (beta(k) - 2*wrk3(k) > upper(k)) then
               beta(k) = upper(k) + 2*wrk3(k)
            end if
            if (beta(k) > upper(k) .or. beta(k) < lower(k)) then
               info = 60001
               return
            end if
            istop = 0
            if (beta(k) == betak) then
               wrk2(1:n, 1:nq) = fn(1:n, 1:nq)
            else
               call fcn(n, m, np, nq, &
                        n, m, np, &
                        beta, xplusd, &
                        ifixb, ifixx, ldifx, &
                        001, wrk2, wrk6, wrk1, &
                        istop)
               if (istop /= 0) then
                  return
               else
                  nfev = nfev + 1
               end if
            end if
            fjacb(1:n, k, 1:nq) = wrk2(1:n, 1:nq)

            beta(k) = beta(k) - 2*wrk3(k)
            if (beta(k) > upper(k)) then
               info = 60001
               return
            end if
            if (beta(k) < lower(k)) then
               info = 60001
               return
            end if
            istop = 0
            if (beta(k) == betak) then
               wrk2(1:n, 1:nq) = fn(1:n, 1:nq)
            else
               call fcn(n, m, np, nq, &
                        n, m, np, &
                        beta, xplusd, &
                        ifixb, ifixx, ldifx, &
                        001, wrk2, wrk6, wrk1, &
                        istop)
               if (istop /= 0) then
                  return
               else
                  nfev = nfev + 1
               end if
            end if

            do l = 1, nq
               do i = 1, n
                  fjacb(i, k, l) = (fjacb(i, k, l) - wrk2(i, l))/(2*wrk3(k))
               end do
            end do
            beta(k) = betak
         end if
      end do

      ! Compute the Jacobian wrt the X's
      if (isodr) then
         do j = 1, m
            if (ifixx(1, 1) < 0) then
               doit = .true.
               setzro = .false.
            elseif (ldifx == 1) then
               if (ifixx(1, j) == 0) then
                  doit = .false.
               else
                  doit = .true.
               end if
               setzro = .false.
            else
               doit = .false.
               setzro = .false.
               do i = 1, n
                  if (ifixx(i, j) /= 0) then
                     doit = .true.
                  else
                     setzro = .true.
                  end if
               end do
            end if
            if (.not. doit) then
               do l = 1, nq
                  fjacd(1:n, j, l) = zero
               end do
            else
               do i = 1, n
                  if (xplusd(i, j) == zero) then
                     if (tt(1, 1) < zero) then
                        typj = one/abs(tt(1, 1))
                     elseif (ldtt == 1) then
                        typj = one/tt(1, j)
                     else
                        typj = one/tt(i, j)
                     end if
                  else
                     typj = abs(xplusd(i, j))
                  end if
                  stp(i) = xplusd(i, j) + sign(one, xplusd(i, j))*typj*dhstep(1, neta, i, j, stpd, ldstpd)
                  stp(i) = stp(i) - xplusd(i, j)
                  xplusd(i, j) = xplusd(i, j) + stp(i)
               end do
               istop = 0
               call fcn(n, m, np, nq, &
                        n, m, np, &
                        beta, xplusd, &
                        ifixb, ifixx, ldifx, &
                        001, wrk2, wrk6, wrk1, &
                        istop)
               if (istop /= 0) then
                  return
               else
                  nfev = nfev + 1
                  do l = 1, nq
                     do i = 1, n
                        fjacd(i, j, l) = wrk2(i, l)
                     end do
                  end do
               end if

               do i = 1, n
                  xplusd(i, j) = x(i, j) + delta(i, j) - stp(i)
               end do
               istop = 0
               call fcn(n, m, np, nq, &
                        n, m, np, &
                        beta, xplusd, &
                        ifixb, ifixx, ldifx, &
                        001, wrk2, wrk6, wrk1, &
                        istop)
               if (istop /= 0) then
                  return
               else
                  nfev = nfev + 1
               end if

               if (setzro) then
                  do i = 1, n
                     if (ifixx(i, j) == 0) then
                        fjacd(i, j, 1:nq) = zero
                     else
                        fjacd(i, j, 1:nq) = (fjacd(i, j, 1:nq) - wrk2(i, 1:nq))/(2*stp(i))
                     end if
                  end do
               else
                  do l = 1, nq
                     fjacd(1:n, j, l) = (fjacd(1:n, j, l) - wrk2(1:n, l))/(2*stp(1:n))
                  end do
               end if
               xplusd(1:n, j) = x(1:n, j) + delta(1:n, j)
            end if
         end do
      end if

   end subroutine djaccd

   subroutine djacfd &
      (fcn, &
       n, m, np, nq, &
       beta, x, ldx, delta, xplusd, ifixb, ifixx, ldifx, &
       stpb, stpd, ldstpd, &
       ssf, tt, ldtt, neta, fn, stp, wrk1, wrk2, wrk3, wrk6, &
       fjacb, isodr, fjacd, nfev, istop, info, &
       lower, upper)
   !! Compute forward difference approximations to the Jacobian wrt the estimated `beta`s and
   !! wrt the `delta`s.
      ! Routines Called  FCN, DHSTEP, DERSTEP
      ! Date Written   860529   (YYMMDD)
      ! Revision Date  920619   (YYMMDD)

      use odrpack_kinds, only: zero, one

      procedure(fcn_t) :: fcn
         !! The user supplied subroutine for evaluating the model.
      integer, intent(in) :: n
         !! The number of observations.
      integer, intent(in) :: m
         !! The number of columns of data in the explanatory variable.
      integer, intent(in) :: np
         !! The number of function parameters.
      integer, intent(in) :: nq
         !! The number of responses per observation.
      real(wp), intent(inout) :: beta(np)
         !! The function parameters.
      real(wp), intent(in) :: x(ldx, m)
         !! The explanatory variable.
      integer, intent(in) :: ldx
         !! The leading dimension of array `x`.
      real(wp), intent(in) :: delta(n, m)
         !! The estimated errors in the explanatory variables.
      real(wp), intent(inout) :: xplusd(n, m)
         !! The values of `x + delta`.
      integer, intent(in) :: ifixb(np)
         !! The values designating whether the elements of `beta` are fixed at their input values or not.
      integer, intent(in) :: ifixx(ldifx, m)
         !! The values designating whether the elements of `x` are fixed at their input values or not.
      integer, intent(in) :: ldifx
         !! The leading dimension of array `ifixx`.
      real(wp), intent(in) :: stpb(np)
         !! The relative step used for computing finite difference derivatives with respect to each `beta`.
      real(wp), intent(in) :: stpd(ldstpd, m)
         !! The relative step used for computing finite difference derivatives with respect to each `delta`.
      integer, intent(in) :: ldstpd
         !! The leading dimension of array `stpd`.
      real(wp), intent(in) :: ssf(np)
         !! The scaling values used for `beta`.
      real(wp), intent(in) :: tt(ldtt, m)
         !! The scaling values used for `delta`.
      integer, intent(in) :: ldtt
         !! The leading dimension of array `tt`.
      integer, intent(in) :: neta
         !! The number of good digits in the function results.
      real(wp), intent(in) :: fn(n, nq)
         !! The new predicted values from the function. Used when parameter is on a boundary.
      real(wp), intent(out) :: stp(n)
         !! The step used for computing finite difference derivatives with respect to each `delta`.
      real(wp), intent(out) :: wrk1(n, m, nq)
         !! A work array of `(n, m, nq)` elements.
      real(wp), intent(out) :: wrk2(n, nq)
         !! A work array of `(n, nq)` elements.
      real(wp), intent(out) :: wrk3(np)
         !! A work array of `(np)` elements.
      real(wp), intent(out) :: wrk6(n, np, nq)
         !! A work array of `(n, np, nq)` elements.
      real(wp), intent(out) :: fjacb(n, np, nq)
         !! The Jacobian with respect to `beta`.
      logical, intent(in) :: isodr
         !! The variable designating whether the solution is by ODR (`isodr = .true.`) or
         !! by OLS (`isodr = .false.`).
      real(wp), intent(out) :: fjacd(n, m, nq)
         !! The Jacobian with respect to `delta`.
      integer, intent(inout) :: nfev
         !! The number of function evaluations.
      integer, intent(out) :: istop
         !! The variable designating whether there are problems computing the function at the
         !! current `beta` and `delta`.
      integer, intent(out) :: info
         !! The variable designating why the computations were stopped.
      real(wp), intent(in) :: lower(np)
         !! The lower bound on `beta`.
      real(wp), intent(in) :: upper(np)
         !! The upper bound on `beta`.

      ! Local scalars
      real(wp) :: betak, step, typj
      integer :: i, j, k, l
      logical :: doit, setzro

      ! Variable Definitions (alphabetically)
      !  BETA:    The function parameters.
      !  BETAK:   The K-th function parameter.
      !  DELTA:   The estimated errors in the explanatory variables.
      !  DOIT:    The variable designating whether the derivative wrt a given BETA or DELTA needs
      !           to be computed (DOIT=TRUE) or not (DOIT=FALSE).
      !  FCN:     The user supplied subroutine for evaluating the model.
      !  FJACB:   The Jacobian with respect to BETA.
      !  FJACD:   The Jacobian with respect to DELTA.
      !  FN:      The new predicted values from the function.
      !  I:       An indexing variable.
      !  IFIXB:   The values designating whether the elements of BETA are fixed at their input
      !           values or not.
      !  IFIXX:   The values designating whether the elements of X are fixed at their input values
      !           or not.
      !  ISODR:   The variable designating whether the solution is by ODR (ISODR=TRUE) or by
      !           OLS (ISODR=FALSE).
      !  ISTOP:   The variable designating whether there are problems computing the function at
      !           the current BETA and DELTA.
      !  J:       An indexing variable.
      !  K:       An indexing variable.
      !  L:       An indexing variable.
      !  LDIFX:   The leading dimension of array IFIXX.
      !  LDSTPD:  The leading dimension of array STPD.
      !  LDTT:    The leading dimension of array TT.
      !  LDX:     The leading dimension of array X.
      !  M:       The number of columns of data in the explanatory variable.
      !  N:       The number of observations.
      !  NETA:    The number of good digits in the function results.
      !  NFEV:    The number of function evaluations.
      !  NP:      The number of function parameters.
      !  SETZRO:  The variable designating whether the derivative wrt some DELTA needs to be set
      !           to zero (SETZRO=TRUE) or not (SETZRO=FALSE).
      !  SSF:     The scale used for the BETA'S.
      !  STP:     The step used for computing finite difference derivatives with respect to DELTA.
      !  STPB:    The relative step used for computing finite difference derivatives with respect
      !           to BETA.
      !  STPD:    The relative step used for computing finite difference derivatives with respect
      !           to DELTA.
      !  TT:      The scaling values used for DELTA.
      !  TYPJ:    The typical size of the J-th unknown BETA or DELTA.
      !  X:       The explanatory variable.
      !  XPLUSD:  The values of X + DELTA.
      !  WRK1:    A work array of (N by M by NQ) elements.
      !  WRK2:    A work array of (N BY NQ) elements.
      !  WRK3:    A work array of (NP) elements.
      !  WRK6:    A work array of (N BY NP BY NQ) elements.

      ! Compute the Jacobian wrt the estimated BETAS
      do k = 1, np
         if (ifixb(1) >= 0) then
            if (ifixb(k) == 0) then
               doit = .false.
            else
               doit = .true.
            end if
         else
            doit = .true.
         end if
         if (.not. doit) then
            do l = 1, nq
               fjacb(1:n, k, l) = zero
            end do
         else
            betak = beta(k)
            step = derstep(0, k, betak, ssf, stpb, neta)
            wrk3(k) = betak + step
            wrk3(k) = wrk3(k) - betak
            beta(k) = betak + wrk3(k)
            if (beta(k) > upper(k)) then
               step = -step
               wrk3(k) = betak + step
               wrk3(k) = wrk3(k) - betak
               beta(k) = betak + wrk3(k)
            end if
            if (beta(k) < lower(k)) then
               step = -step
               wrk3(k) = betak + step
               wrk3(k) = wrk3(k) - betak
               beta(k) = betak + wrk3(k)
               if (beta(k) > upper(k)) then
                  info = 60001
                  return
               end if
            end if
            istop = 0
            call fcn(n, m, np, nq, &
                     n, m, np, &
                     beta, xplusd, &
                     ifixb, ifixx, ldifx, &
                     001, wrk2, wrk6, wrk1, &
                     istop)
            if (istop /= 0) then
               return
            else
               nfev = nfev + 1
            end if
            do l = 1, nq
               do i = 1, n
                  fjacb(i, k, l) = (wrk2(i, l) - fn(i, l))/wrk3(k)
               end do
            end do
            beta(k) = betak
         end if
      end do

      ! Compute the Jacobian wrt the X'S
      if (isodr) then
         do j = 1, m
            if (ifixx(1, 1) < 0) then
               doit = .true.
               setzro = .false.
            elseif (ldifx == 1) then
               if (ifixx(1, j) == 0) then
                  doit = .false.
               else
                  doit = .true.
               end if
               setzro = .false.
            else
               doit = .false.
               setzro = .false.
               do i = 1, n
                  if (ifixx(i, j) /= 0) then
                     doit = .true.
                  else
                     setzro = .true.
                  end if
               end do
            end if
            if (.not. doit) then
               do l = 1, nq
                  fjacd(1:n, j, l) = zero
               end do
            else
               do i = 1, n
                  if (xplusd(i, j) == zero) then
                     if (tt(1, 1) < zero) then
                        typj = one/abs(tt(1, 1))
                     elseif (ldtt == 1) then
                        typj = one/tt(1, j)
                     else
                        typj = one/tt(i, j)
                     end if
                  else
                     typj = abs(xplusd(i, j))
                  end if

                  stp(i) = xplusd(i, j) &
                           + sign(one, xplusd(i, j)) &
                           *typj*dhstep(0, neta, i, j, stpd, ldstpd)
                  stp(i) = stp(i) - xplusd(i, j)
                  xplusd(i, j) = xplusd(i, j) + stp(i)
               end do

               istop = 0
               call fcn(n, m, np, nq, &
                        n, m, np, &
                        beta, xplusd, &
                        ifixb, ifixx, ldifx, &
                        001, wrk2, wrk6, wrk1, &
                        istop)
               if (istop /= 0) then
                  return
               else
                  nfev = nfev + 1
                  do l = 1, nq
                     do i = 1, n
                        fjacd(i, j, l) = wrk2(i, l)
                     end do
                  end do

               end if

               if (setzro) then
                  do i = 1, n
                     if (ifixx(i, j) == 0) then
                        do l = 1, nq
                           fjacd(i, j, l) = zero
                        end do
                     else
                        do l = 1, nq
                           fjacd(i, j, l) = (fjacd(i, j, l) - fn(i, l))/ &
                                            stp(i)
                        end do
                     end if
                  end do
               else
                  do l = 1, nq
                     do i = 1, n
                        fjacd(i, j, l) = (fjacd(i, j, l) - fn(i, l))/stp(i)
                     end do
                  end do
               end if
               do i = 1, n
                  xplusd(i, j) = x(i, j) + delta(i, j)
               end do
            end if
         end do
      end if

   end subroutine djacfd

   subroutine djck &
      (fcn, &
       n, m, np, nq, &
       beta, betaj, xplusd, &
       ifixb, ifixx, ldifx, stpb, stpd, ldstpd, &
       ssf, tt, ldtt, &
       eta, neta, ntol, nrow, isodr, epsmac, &
       pv0i, fjacb, fjacd, &
       msgb, msgd, diff, istop, nfev, njev, &
       wrk1, wrk2, wrk6, &
       interval)
   !! Driver routine for the derivative checking process.
      ! Adapted from STARPAC subroutine DCKCNT.
      ! Routines Called  FCN, DHSTEP, DJCKM
      ! Date Written   860529   (YYMMDD)
      ! Revision Date  920619   (YYMMDD)

      use odrpack_kinds, only: zero, one, p5 => half

      procedure(fcn_t) :: fcn
         !! The user supplied subroutine for evaluating the model.
      integer, intent(in) :: n
         !! The number of observations.
      integer, intent(in) :: m
         !! The number of columns of data in the explanatory variable.
      integer, intent(in) :: np
         !! The number of function parameters.
      integer, intent(in) :: nq
         !! The number of responses per observation.
      real(wp), intent(inout) :: beta(np)
         !! The function parameters.
      real(wp), intent(inout) :: betaj(np)
         !! The function parameters offset such that steps don't cross bounds.
      real(wp), intent(inout) :: xplusd(n, m)
         !! The values of `x + delta`.
      integer, intent(in) :: ifixb(np)
         !! The values designating whether the elements of `beta` are fixed at their input values or not.
      integer, intent(in) :: ifixx(ldifx, m)
         !! The values designating whether the elements of `x` are fixed at their input values or not.
      integer, intent(in) :: ldifx
         !! The leading dimension of array `ifixx`.
      real(wp), intent(in) :: stpb(np)
         !! The step size for finite difference derivatives wrt `beta`.
      real(wp), intent(in) :: stpd(ldstpd, m)
         !! The step size for finite difference derivatives wrt `delta`.
      integer, intent(in) :: ldstpd
         !! The leading dimension of array `stpd`.
      real(wp), intent(in) :: ssf(np)
         !! The scaling values used for `beta`.
      real(wp), intent(in) :: tt(ldtt, m)
         !! The scaling values used for `delta`.
      integer, intent(in) :: ldtt
         !! The leading dimension of array `tt`.
      real(wp), intent(in) :: eta
         !! The relative noise in the function results.
      integer, intent(in) :: neta
         !! The number of reliable digits in the model results.
      integer, intent(out) :: ntol
         !! The number of digits of agreement required between the numerical derivatives and the
         !! user supplied derivatives.
      integer, intent(in) :: nrow
         !! The row number of the explanatory variable array at which the derivative is checked.
      logical, intent(in) :: isodr
         !! The variable designating whether the solution is by ODR (`isodr = .true.`) or
         !! by OLS (`isodr = .false.`).
      real(wp), intent(in) :: epsmac
         !! The value of machine precision.
      real(wp), intent(in) :: pv0i(n, nq)
         !! The predicted values using the user supplied parameter estimates.
      real(wp), intent(out) :: fjacb(n, np, nq)
         !! The Jacobian with respect to `beta`.
      real(wp), intent(out) :: fjacd(n, m, nq)
         !! The Jacobian with respect to `delta`.
      integer, intent(out) :: msgb(1 + nq*np)
         !! The error checking results for the Jacobian wrt `beta`.
      integer, intent(out) :: msgd(1 + nq*m)
         !! The error checking results for the Jacobian wrt `delta`.
      real(wp), intent(out) :: diff(nq, np + m)
         !! The relative differences between the user supplied and finite difference derivatives
         !! for each derivative checked.
      integer, intent(out) :: istop
         !! The variable designating whether there are problems computing the function at the
         !! current `beta` and `delta`
      integer, intent(inout) :: nfev
         !! The number of function evaluations.
      integer, intent(inout) :: njev
         !! The number of Jacobian evaluations.
      real(wp), intent(out) :: wrk1(n, m, nq)
         !! A work array of `(n, m, nq)` elements.
      real(wp), intent(out) :: wrk2(n, nq)
         !! A work array of `(n, nq)` elements.
      real(wp), intent(out) :: wrk6(n, np, nq)
         !! A work array of `(n, np, nq)` elements.
      integer, intent(in) :: interval(np)
         !! Specifies which checks can be performed when checking derivatives based on the
         !! interval of the bound constraints.

      ! Local scalars
      real(wp) :: diffj, h0, hc0, pv, tol, typj
      integer :: ideval, j, lq, msgb1, msgd1
      logical :: isfixd, iswrtb

      ! Local arrays
      real(wp) :: pv0(n, nq)

      ! Variable Definitions (alphabetically)
      !  BETA:     The function parameters.
      !  BETAJ:    The function parameters offset such that steps don't cross bounds.
      !  DIFF:     The relative differences between the user supplied and finite difference
      !            derivatives for each derivative checked.
      !  DIFFJ:    The relative differences between the user supplied and finite difference
      !            derivatives for the derivative being checked.
      !  EPSMAC:   The value of machine precision.
      !  ETA:      The relative noise in the function results.
      !  FCN:      The user supplied subroutine for evaluating the model.
      !  FJACB:    The Jacobian with respect to BETA.
      !  FJACD:    The Jacobian with respect to DELTA.
      !  H0:       The initial relative step size for forward differences.
      !  HC0:      The initial relative step size for central differences.
      !  IDEVAL:   The variable designating what computations are to be performed by user supplied
      !            subroutine FCN.
      !  IFIXB:    The values designating whether the elements of BETA are fixed at their input
      !            values or not.
      !  IFIXX:    The values designating whether the elements of X are fixed at their input values or not.
      !  INTERVAL: Specifies which checks can be performed when checking derivatives based on the
      !            interval of the bound constraints.
      !  ISFIXD:   The variable designating whether the parameter is fixed (ISFIXD=TRUE) or not (ISFIXD=FALSE).
      !  ISTOP:    The variable designating whether there are problems computing the function at the
      !            current BETA and DELTA.
      !  ISODR:    The variable designating whether the solution is by ODR (ISODR=.TRUE.) or by
      !            OLS (ISODR=.FALSE.).
      !  ISWRTB:   The variable designating whether the derivatives wrt BETA (ISWRTB=TRUE) or DELTA
      !            (ISWRTB=FALSE) are being checked.
      !  J:        An index variable.
      !  LDIFX:    The leading dimension of array IFIXX.
      !  LDSTPD:   The leading dimension of array STPD.
      !  LDTT:     The leading dimension of array TT.
      !  LQ:       The response currently being examined.
      !  M:        The number of columns of data in the explanatory variable.
      !  MSGB:     The error checking results for the Jacobian wrt BETA.
      !  MSGB1:    The error checking results for the Jacobian wrt BETA.
      !  MSGD:     The error checking results for the Jacobian wrt DELTA.
      !  MSGD1:    The error checking results for the Jacobian wrt DELTA.
      !  N:        The number of observations.
      !  NETA:     The number of reliable digits in the model results, either set by the user or
      !            computed by DETAF.
      !  NFEV:     The number of function evaluations.
      !  NJEV:     The number of Jacobian evaluations.
      !  NP:       The number of function parameters.
      !  NQ:       The number of responses per observation.
      !  NROW:     The row number of the explanatory variable array at which the derivative is checked.
      !  NTOL:     The number of digits of agreement required between the numerical derivatives and
      !            the user supplied derivatives.
      !  PV:       The scalar in which the predicted value from the model for row NROW is stored.
      !  PV0:      The predicted values using the current parameter estimates (possibly offset from
      !            the user supplied estimates to create distance between parameters and the bounds
      !            on the parameters).
      !  PV0I:     The predicted values using the user supplied parameter estimates.
      !  SSF:      The scaling values used for BETA.
      !  STPB:     The step size for finite difference derivatives wrt BETA.
      !  STPD:     The step size for finite difference derivatives wrt DELTA.
      !  TOL:      The agreement tolerance.
      !  TT:       The scaling values used for DELTA.
      !  TYPJ:     The typical size of the J-th unknown BETA or DELTA.
      !  WRK1:     A work array of (N BY M BY NQ) elements.
      !  WRK2:     A work array of (N BY NQ) elements.
      !  WRK6:     A work array of (N BY NP BY NQ) elements.
      !  XPLUSD:   The values of X + DELTA.

      ! Set tolerance for checking derivatives
      tol = eta**(0.25E0_wp)
      ntol = int(max(one, p5 - log10(tol)))

      ! Compute, if necessary, PV0
      pv0 = pv0i
      if (any(beta(:) /= betaj(:))) then
         istop = 0
         ideval = 001
         call fcn(n, m, np, nq, &
                  n, m, np, &
                  betaj, xplusd, &
                  ifixb, ifixx, ldifx, &
                  ideval, pv0, fjacb, fjacd, &
                  istop)
         if (istop /= 0) then
            return
         else
            njev = njev + 1
         end if
      end if

      ! Compute user-supplied derivative values
      istop = 0
      if (isodr) then
         ideval = 110
      else
         ideval = 010
      end if
      call fcn(n, m, np, nq, &
               n, m, np, &
               betaj, xplusd, &
               ifixb, ifixx, ldifx, &
               ideval, wrk2, fjacb, fjacd, &
               istop)
      if (istop /= 0) then
         return
      else
         njev = njev + 1
      end if

      ! Check derivatives wrt BETA for each response of observation NROW
      msgb1 = 0
      msgd1 = 0

      do lq = 1, nq

         ! Set predicted value of model at current parameter estimates
         pv = pv0(nrow, lq)

         iswrtb = .true.
         do j = 1, np

            if (ifixb(1) < 0) then
               isfixd = .false.
            elseif (ifixb(j) == 0) then
               isfixd = .true.
            else
               isfixd = .false.
            end if

            if (isfixd) then
               msgb(1 + lq + (j - 1)*nq) = -1
            else
               if (beta(j) == zero) then
                  if (ssf(1) < zero) then
                     typj = one/abs(ssf(1))
                  else
                     typj = one/ssf(j)
                  end if
               else
                  typj = abs(beta(j))
               end if

               h0 = dhstep(0, neta, 1, j, stpb, 1)
               hc0 = h0

               ! Check derivative wrt the J-th parameter at the NROW-th row
               if (interval(j) >= 1) then
                  call djckm(fcn, &
                             n, m, np, nq, &
                             betaj, xplusd, &
                             ifixb, ifixx, ldifx, &
                             eta, tol, nrow, epsmac, j, lq, typj, h0, hc0, &
                             iswrtb, pv, fjacb(nrow, j, lq), &
                             diffj, msgb1, msgb(2), istop, nfev, &
                             wrk1, wrk2, wrk6, interval)
                  if (istop /= 0) then
                     msgb(1) = -1
                     return
                  else
                     diff(lq, j) = diffj
                  end if
               else
                  msgb(1 + j) = 9
               end if
            end if

         end do

         ! Check derivatives wrt X for each response of observation NROW
         if (isodr) then
            iswrtb = .false.
            do j = 1, m

               if (ifixx(1, 1) < 0) then
                  isfixd = .false.
               elseif (ldifx == 1) then
                  if (ifixx(1, j) == 0) then
                     isfixd = .true.
                  else
                     isfixd = .false.
                  end if
               else
                  isfixd = .false.
               end if

               if (isfixd) then
                  msgd(1 + lq + (j - 1)*nq) = -1
               else
                  if (xplusd(nrow, j) == zero) then
                     if (tt(1, 1) < zero) then
                        typj = one/abs(tt(1, 1))
                     elseif (ldtt == 1) then
                        typj = one/tt(1, j)
                     else
                        typj = one/tt(nrow, j)
                     end if
                  else
                     typj = abs(xplusd(nrow, j))
                  end if

                  h0 = dhstep(0, neta, nrow, j, stpd, ldstpd)
                  hc0 = dhstep(1, neta, nrow, j, stpd, ldstpd)

                  ! Check derivative wrt the J-th column of DELTA at row NROW
                  call djckm(fcn, &
                             n, m, np, nq, &
                             betaj, xplusd, &
                             ifixb, ifixx, ldifx, &
                             eta, tol, nrow, epsmac, j, lq, typj, h0, hc0, &
                             iswrtb, pv, fjacd(nrow, j, lq), &
                             diffj, msgd1, msgd(2), istop, nfev, &
                             wrk1, wrk2, wrk6, interval)
                  if (istop /= 0) then
                     msgd(1) = -1
                     return
                  else
                     diff(lq, np + j) = diffj
                  end if
               end if

            end do
         end if
      end do

      msgb(1) = msgb1
      msgd(1) = msgd1

   end subroutine djck

   subroutine djckc &
      (fcn, &
       n, m, np, nq, &
       beta, xplusd, ifixb, ifixx, ldifx, &
       eta, tol, nrow, epsmac, j, lq, hc, iswrtb, &
       fd, typj, pvpstp, stp0, &
       pv, d, &
       diffj, msg, istop, nfev, &
       wrk1, wrk2, wrk6)
   !! Check whether high curvature could be the cause of the disagreement between the numerical
   !! and analytic derviatives.
      ! Adapted from STARPAC subroutine DCKCRV.
      ! Routines Called  DJCKF, DPVB, DPVD
      ! Date Written   860529   (YYMMDD)
      ! Revision Date  920619   (YYMMDD)

      use odrpack_kinds, only: one, two, ten

      procedure(fcn_t) :: fcn
         !! The user supplied subroutine for evaluating the model.
      integer, intent(in) :: n
         !! The number of observations.
      integer, intent(in) :: m
         !! The number of columns of data in the explanatory variable.
      integer, intent(in) :: np
         !! The number of function parameters.
      integer, intent(in) :: nq
         !! The number of responses per observation.
      real(wp), intent(inout) :: beta(np)
         !! The function parameters.
      real(wp), intent(inout) :: xplusd(n, m)
         !! The values of `x` + `delta`.
      integer, intent(in) :: ifixb(np)
         !! The values designating whether the elements of `beta` are fixed at their input values or not.
      integer, intent(in) :: ifixx(ldifx, m)
         !! The values designating whether the elements of `x` are fixed at their input values or not.
      integer, intent(in) :: ldifx
         !! The leading dimension of array `ifixx`.
      real(wp), intent(in) :: eta
         !! The relative noise in the model.
      real(wp), intent(in) :: tol
         !! The agreement tolerance.
      integer, intent(in) :: nrow
         !! The row number of the explanatory variable array at which the derivative is to be checked.
      real(wp), intent(in) :: epsmac
         !! The value of machine precision.
      integer, intent(in) :: j
         !! The index of the partial derivative being examined.
      integer, intent(in) :: lq
         !! The response currently being examined.
      real(wp), intent(in) :: hc
         !! The relative step size for central finite differences.
      logical, intent(in) :: iswrtb
         !! The variable designating whether the derivatives wrt `beta` (`iswrtb = .true.`) or
         !! `delta` (`iswrtb = .false.`) are being checked.
      real(wp), intent(out) :: fd
         !! The forward difference derivative wrt the `j`-th parameter.
      real(wp), intent(in) :: typj
         !! The typical size of the `j`-th unknown `beta` or `delta`.
      real(wp), intent(out) :: pvpstp
         !! The predicted value for row `nrow` of the model based on the current parameter estimates
         !! for all but the `j`-th parameter value, which is `beta(j) + stp0`.
      real(wp), intent(in) :: stp0
         !! The initial step size for the finite difference derivative.
      real(wp), intent(in) :: pv
         !! The predicted value of the model for row `nrow`.
      real(wp), intent(in) :: d
         !! The derivative with respect to the `j`-th unknown parameter.
      real(wp), intent(out) :: diffj
         !! The relative differences between the user supplied and finite difference derivatives
         !! for the derivative being checked.
      integer, intent(out) :: msg(nq, j)
         !! The error checking results.
      integer, intent(out) :: istop
         !! The variable designating whether there are problems computing the function at the
         !! current `beta` and `delta`.
      integer, intent(inout) :: nfev
         !! The number of function evaluations.
      real(wp), intent(out) :: wrk1(n, m, nq)
         !! A work array of `(n, m, nq)` elements.
      real(wp), intent(out) :: wrk2(n, nq)
         !! A work array of `(n, nq)` elements.
      real(wp), intent(out) :: wrk6(n, np, nq)
         !! A work array of `(n, np, nq)` elements.

      ! Local scalars
      real(wp), parameter :: p01 = 0.01_wp
      real(wp) :: curve, pvmcrv, pvpcrv, stp, stpcrv

      ! Variable Definitions (alphabetically)
      !  BETA:    The function parameters.
      !  CURVE:   A measure of the curvature in the model.
      !  D:       The derivative with respect to the Jth unknown parameter.
      !  DIFFJ:   The relative differences between the user supplied and finite difference
      !           derivatives for the derivative being checked.
      !  EPSMAC:  The value of machine precision.
      !  ETA:     The relative noise in the model.
      !  FCN:     The user supplied subroutine for evaluating the model.
      !  FD:      The forward difference derivative wrt the Jth parameter.
      !  HC:      The relative step size for central finite differences.
      !  IFIXB:   The values designating whether the elements of BETA are fixed at their input
      !           values or not.
      !  IFIXX:   The values designating whether the elements of X are fixed at their input values
      !           or not.
      !  ISTOP:   The variable designating whether there are problems computing the function at the
      !           current BETA and DELTA.
      !  ISWRTB:  The variable designating whether the derivatives wrt BETA (ISWRTB=TRUE) or
      !           DELTA(ISWRTB=FALSE) are being checked.
      !  J:       The index of the partial derivative being examined.
      !  LDIFX:   The leading dimension of array IFIXX.
      !  LQ:      The response currently being examined.
      !  M:       The number of columns of data in the explanatory variable.
      !  MSG:     The error checking results.
      !  N:       The number of observations.
      !  NFEV:    The number of function evaluations.
      !  NP:      The number of function parameters.
      !  NQ:      The number of responses per observation.
      !  NROW:    The row number of the explanatory variable array at which the derivative is to be
      !           checked.
      !  PV:      The predicted value of the model for row NROW.
      !  PVMCRV:  The predicted value for row  NROW of the model based on the current parameter
      !           estimates for all but the Jth parameter value, which is BETA(J)-STPCRV.
      !  PVPCRV:  The predicted value for row NROW of the model based on the current parameter
      !           estimates for all but the Jth parameter value, which is BETA(J)+STPCRV.
      !  PVPSTP:  The predicted value for row NROW of the model based on the current parameter
      !           estimates for all but the Jth parameter value, which is BETA(J) + STP0.
      !  P01:     The value 0.01E0_wp.
      !  STP0:    The initial step size for the finite difference derivative.
      !  STP:     A step size for the finite difference derivative.
      !  STPCRV:  The step size selected to check for curvature in the model.
      !  TOL:     The agreement tolerance.
      !  TYPJ:    The typical size of the J-th unknown BETA or DELTA.
      !  WRK1:    A work array of (N BY M BY NQ) elements.
      !  WRK2:    A work array of (N BY NQ) elements.
      !  WRK6:    A work array of (N BY NP BY NQ) elements.
      !  XPLUSD:  The values of X + DELTA.

      if (iswrtb) then

         ! Perform central difference computations for derivatives wrt BETA
         stpcrv = (hc*typj*sign(one, beta(j)) + beta(j)) - beta(j)
         call dpvb(fcn, &
                   n, m, np, nq, &
                   beta, xplusd, ifixb, ifixx, ldifx, &
                   nrow, j, lq, stpcrv, &
                   istop, nfev, pvpcrv, &
                   wrk1, wrk2, wrk6)
         if (istop /= 0) then
            return
         end if
         call dpvb(fcn, &
                   n, m, np, nq, &
                   beta, xplusd, ifixb, ifixx, ldifx, &
                   nrow, j, lq, -stpcrv, &
                   istop, nfev, pvmcrv, &
                   wrk1, wrk2, wrk6)
         if (istop /= 0) then
            return
         end if
      else

         ! Perform central difference computations for derivatives wrt DELTA
         stpcrv = (hc*typj*sign(one, xplusd(nrow, j)) + xplusd(nrow, j)) - xplusd(nrow, j)
         call dpvd(fcn, &
                   n, m, np, nq, &
                   beta, xplusd, ifixb, ifixx, ldifx, &
                   nrow, j, lq, stpcrv, &
                   istop, nfev, pvpcrv, &
                   wrk1, wrk2, wrk6)
         if (istop /= 0) then
            return
         end if
         call dpvd(fcn, &
                   n, m, np, nq, &
                   beta, xplusd, ifixb, ifixx, ldifx, &
                   nrow, j, lq, -stpcrv, &
                   istop, nfev, pvmcrv, &
                   wrk1, wrk2, wrk6)
         if (istop /= 0) then
            return
         end if
      end if

      ! Estimate curvature by second derivative of model
      curve = abs((pvpcrv - pv) + (pvmcrv - pv))/(stpcrv*stpcrv)
      curve = curve + eta*(abs(pvpcrv) + abs(pvmcrv) + two*abs(pv))/(stpcrv**2)

      ! Check if finite precision arithmetic could be the culprit.
      call djckf(fcn, &
                 n, m, np, nq, &
                 beta, xplusd, ifixb, ifixx, ldifx, &
                 eta, tol, nrow, j, lq, iswrtb, &
                 fd, typj, pvpstp, stp0, curve, pv, d, &
                 diffj, msg, istop, nfev, &
                 wrk1, wrk2, wrk6)
      if (istop /= 0) then
         return
      end if
      if (msg(lq, j) == 0) then
         return
      end if

      ! Check if high curvature could be the problem.
      stp = two*max(tol*abs(d)/curve, epsmac)
      if (stp < abs(ten*stp0)) then
         stp = min(stp, p01*abs(stp0))
      end if

      if (iswrtb) then
         ! Perform computations for derivatives wrt BETA
         stp = (stp*sign(one, beta(j)) + beta(j)) - beta(j)
         call dpvb(fcn, &
                   n, m, np, nq, &
                   beta, xplusd, ifixb, ifixx, ldifx, &
                   nrow, j, lq, stp, &
                   istop, nfev, pvpstp, &
                   wrk1, wrk2, wrk6)
         if (istop /= 0) then
            return
         end if
      else

         ! Perform computations for derivatives wrt DELTA
         stp = (stp*sign(one, xplusd(nrow, j)) + xplusd(nrow, j)) - xplusd(nrow, j)
         call dpvd(fcn, &
                   n, m, np, nq, &
                   beta, xplusd, ifixb, ifixx, ldifx, &
                   nrow, j, lq, stp, &
                   istop, nfev, pvpstp, &
                   wrk1, wrk2, wrk6)
         if (istop /= 0) then
            return
         end if
      end if

      ! Compute the new numerical derivative
      fd = (pvpstp - pv)/stp
      diffj = min(diffj, abs(fd - d)/abs(d))

      ! Check whether the new numerical derivative is ok
      if (abs(fd - d) <= tol*abs(d)) then
         msg(lq, j) = 0

         ! Check if finite precision may be the culprit (fudge factor = 2)
      elseif (abs(stp*(fd - d)) < two*eta*(abs(pv) + abs(pvpstp)) + &
              curve*(epsmac*typj)**2) then
         msg(lq, j) = 5
      end if

   end subroutine djckc

   subroutine djckf &
      (fcn, &
       n, m, np, nq, &
       beta, xplusd, ifixb, ifixx, ldifx, &
       eta, tol, nrow, j, lq, iswrtb, &
       fd, typj, pvpstp, stp0, curve, pv, d, &
       diffj, msg, istop, nfev, &
       wrk1, wrk2, wrk6)
   !! Check whether finite precision arithmetic could be the cause of the disagreement between
   !! the derivatives.
      ! Adapted from STARPAC subroutine DCKFPA.
      ! Routines Called  DPVB, DPVD
      ! Date Written   860529   (YYMMDD)
      ! Revision Date  920619   (YYMMDD)

      use odrpack_kinds, only: one, two, hundred

      procedure(fcn_t) :: fcn
         !! The user supplied subroutine for evaluating the model.
      integer, intent(in) :: n
         !! The number of observations.
      integer, intent(in) :: m
         !! The number of columns of data in the explanatory variable.
      integer, intent(in) :: np
         !! The number of function parameters.
      integer, intent(in) :: nq
         !! The number of responses per observation.
      real(wp), intent(inout) :: beta(np)
         !! The function parameters.
      real(wp), intent(inout) :: xplusd(n, m)
         !! The values of `x + delta`.
      integer, intent(in) :: ifixb(np)
         !! The values designating whether the elements of `beta` are fixed at their input values or not.
      integer, intent(in) :: ifixx(ldifx, m)
         !! The values designating whether the elements of `x` are fixed at their input values or not.
      integer, intent(in) :: ldifx
         !! The leading dimension of array `ifixx`.
      real(wp), intent(in) :: eta
         !! The relative noise in the model.
      real(wp), intent(in) :: tol
         !! The agreement tolerance.
      integer, intent(in) :: nrow
         !! The row number of the explanatory variable array at which the derivative is to be checked.
      integer, intent(in) :: j
         !! The index of the partial derivative being examined.
      integer, intent(in) :: lq
         !! The response currently being examined.
      logical, intent(in) :: iswrtb
         !! The variable designating whether the derivatives wrt `beta` (`iswrtb = .true.`)
         !! or `delta` (`iswrtb = .false.`) are being checked.
      real(wp), intent(out) :: fd
         !! The forward difference derivative wrt the `j`-th parameter.
      real(wp), intent(in) :: typj
         !! The typical size of the `j`-th unknown `beta` or `delta`.
      real(wp), intent(out) :: pvpstp
         !! The predicted value for row `nrow` of the model based on the current parameter
         !! estimates for all but the `j`-th parameter value, which is `beta(j) + stp0`.
      real(wp), intent(in) :: stp0
         !! The step size for the finite difference derivative.
      real(wp), intent(inout) :: curve
         !! A measure of the curvature in the model.
      real(wp), intent(in) :: pv
         !! The predicted value for row `nrow`.
      real(wp), intent(in) :: d
         !! The derivative with respect to the `j`-th unknown parameter.
      real(wp), intent(out) :: diffj
         !! The relative differences between the user supplied and finite difference derivatives
         !! for the derivative being checked.
      integer, intent(out) :: msg(nq, j)
         !! The error checking results.
      integer, intent(out) :: istop
         !! The variable designating whether there are problems computing the function at the
         !! current `beta` and `delta`.
      integer, intent(inout) :: nfev
         !! The number of function evaluations.
      real(wp), intent(out) :: wrk1(n, m, nq)
         !! A work array of `(n, m, nq)` elements.
      real(wp), intent(out) :: wrk2(n, nq)
         !! A work array of `(n, nq)` elements.
      real(wp), intent(out) :: wrk6(n, np, nq)
         !! A work array of `(n, np, nq)` elements.

      ! Local scalars
      real(wp), parameter :: p1 = 0.1_wp
      real(wp) :: stp
      logical :: large

      ! Variable Definitions (alphabetically)
      !  BETA:    The function parameters.
      !  CURVE:   A measure of the curvature in the model.
      !  D:       The derivative with respect to the Jth unknown parameter.
      !  DIFFJ:   The relative differences between the user supplied and finite difference
      !           derivatives for the derivative being checked.
      !  ETA:     The relative noise in the model.
      !  FCN:     The user supplied subroutine for evaluating the model.
      !  FD:      The forward difference derivative wrt the Jth parameter.
      !  IFIXB:   The values designating whether the elements of BETA are fixed at their input
      !           values or not.
      !  IFIXX:   The values designating whether the elements of X are fixed at their input values
      !           or not.
      !  ISTOP:   The variable designating whether there are problems computing the function at
      !           the current BETA and DELTA.
      !  ISWRTB:  The variable designating whether the derivatives wrt BETA (ISWRTB=TRUE) or
      !           DELTA(ISWRTB=FALSE) are being checked.
      !  J:       The index of the partial derivative being examined.
      !  LARGE:   The value designating whether the recommended increase in the step size would
      !           be greater than TYPJ.
      !  LDIFX:   The leading dimension of array IFIXX.
      !  LQ:      The response currently being examined.
      !  M:       The number of columns of data in the explanatory variable.
      !  MSG:     The error checking results.
      !  N:       The number of observations.
      !  NFEV:    The number of function evaluations.
      !  NP:      The number of function parameters.
      !  NQ:      The number of responses per observation.
      !  NROW:    The row number of the explanatory variable array at which the derivative is to
      !           be checked.
      !  PV:      The predicted value for row NROW.
      !  PVPSTP:  The predicted value for row NROW of the model based on the current parameter
      !           estimates for all but the Jth parameter value, which is BETA(J) + STP0.
      !  P1:      The value 0.1E0_wp.
      !  STP0:    The step size for the finite difference derivative.
      !  TOL:     The agreement tolerance.
      !  TYPJ:    The typical size of the J-th unknown BETA or DELTA.
      !  WRK1:    A work array of (N BY M BY NQ) elements.
      !  WRK2:    A work array of (N BY NQ) elements.
      !  WRK6:    A work array of (N BY NP BY NQ) elements.
      !  XPLUSD:  The values of X + DELTA.

      ! Finite precision arithmetic could be the problem.
      ! Try a larger step size based on estimate of condition error.
      stp = eta*(abs(pv) + abs(pvpstp))/(tol*abs(d))
      if (stp > abs(p1*stp0)) then
         stp = max(stp, hundred*abs(stp0))
      end if
      if (stp > typj) then
         stp = typj
         large = .true.
      else
         large = .false.
      end if

      if (iswrtb) then
         ! Perform computations for derivatives wrt BETA
         stp = (stp*sign(one, beta(j)) + beta(j)) - beta(j)
         call dpvb(fcn, &
                   n, m, np, nq, &
                   beta, xplusd, ifixb, ifixx, ldifx, &
                   nrow, j, lq, stp, &
                   istop, nfev, pvpstp, &
                   wrk1, wrk2, wrk6)
      else
         ! Perform computations for derivatives wrt DELTA
         stp = (stp*sign(one, xplusd(nrow, j)) + xplusd(nrow, j)) - xplusd(nrow, j)
         call dpvd(fcn, &
                   n, m, np, nq, &
                   beta, xplusd, ifixb, ifixx, ldifx, &
                   nrow, j, lq, stp, &
                   istop, nfev, pvpstp, &
                   wrk1, wrk2, wrk6)
      end if
      if (istop /= 0) then
         return
      end if

      fd = (pvpstp - pv)/stp
      diffj = min(diffj, abs(fd - d)/abs(d))

      ! Check for agreement
      if ((abs(fd - d)) <= tol*abs(d)) then
         ! Forward difference quotient and analytic derivatives agree.
         msg(lq, j) = 0

      elseif ((abs(fd - d) <= abs(two*curve*stp)) .or. large) then
         ! Curvature may be the culprit (fudge factor = 2)
         if (large) then
            msg(lq, j) = 4
         else
            msg(lq, j) = 5
         end if
      end if

   end subroutine djckf

   subroutine djckm &
      (fcn, &
       n, m, np, nq, &
       beta, xplusd, ifixb, ifixx, ldifx, &
       eta, tol, nrow, epsmac, j, lq, typj, h0, hc0, &
       iswrtb, pv, d, &
       diffj, msg1, msg, istop, nfev, &
       wrk1, wrk2, wrk6, interval)
   !! Check user supplied analytic derivatives against numerical derivatives.
      ! Adapted from STARPAC subroutine DCKMN.
      ! Routines Called  DJCKC, DJCKZ, DPVB, DPVD
      ! Date Written   860529   (YYMMDD)
      ! Revision Date  920619   (YYMMDD)

      use odrpack_kinds, only: zero, one, two, three, ten, hundred

      procedure(fcn_t) :: fcn
         !! The user supplied subroutine for evaluating the model.
      integer, intent(in) :: n
         !! The number of observations.
      integer, intent(in) :: m
         !! The number of columns of data in the explanatory variable.
      integer, intent(in) :: np
         !! The number of function parameters.
      integer, intent(in) :: nq
         !! The number of responses per observation.
      real(wp), intent(inout) :: beta(np)
         !! The function parameters.
      real(wp), intent(inout) :: xplusd(n, m)
         !! The values of `x + delta`.
      integer, intent(in) :: ifixb(np)
         !! The values designating whether the elements of `beta` are fixed at their input values or not.
      integer, intent(in) :: ifixx(ldifx, m)
         !! The values designating whether the elements of `x` are fixed at their input values or not.
      integer, intent(in) :: ldifx
         !! The leading dimension of array `ifixx`.
      real(wp), intent(in) :: eta
         !! The relative noise in the model.
      real(wp), intent(in) :: tol
         !! The agreement tolerance.
      integer, intent(in) :: nrow
         !! The row number of the explanatory variable array at which the derivative is to be checked.
      real(wp), intent(in) :: epsmac
         !! The value of machine precision.
      integer, intent(in) :: j
         !! The index of the partial derivative being examined.
      integer, intent(in) :: lq
         !! The response currently being examined.
      real(wp), intent(in) :: typj
         !! The typical size of the `j`-th unknown `beta` or `delta`.
      real(wp), intent(in) :: h0
         !! The initial step size for the finite difference derivative.
      real(wp), intent(in) :: hc0
         !! The relative step size for central finite differences.
      logical, intent(in) :: iswrtb
         !! The variable designating whether the derivatives wrt `beta` (`iswrtb = .true.`)
         !! or `delta` (`iswrtb = .false.`) are being checked.
      real(wp), intent(in) :: pv
         !! The predicted value for row `nrow`.
      real(wp), intent(in) :: d
         !! The derivative with respect to the `j`-th unknown parameter.
      real(wp), intent(out) :: diffj
         !! The relative differences between the user supplied and finite difference derivatives
         !! for the derivative being checked.
      integer, intent(out) :: msg1
         !! The first set of error checking results.
      integer, intent(out) :: msg(nq, j)
         !! The error checking results.
      integer, intent(out) :: istop
         !! The variable designating whether there are problems computing the function at the
         !! current `beta` and `delta`.
      integer, intent(inout) :: nfev
         !! The number of function evaluations.
      real(wp), intent(out) :: wrk1(n, m, nq)
         !! A work array of `(n, m, nq)` elements.
      real(wp), intent(out) :: wrk2(n, nq)
         !! A work array of `(n, nq)` elements.
      real(wp), intent(out) :: wrk6(n, np, nq)
         !! A work array of `(n, np, nq)` elements.
      integer, intent(in) :: interval(np)
         !! Specifies which checks can be performed when checking derivatives based on the
         !! interval of the bound constraints.

      ! Local scalars
      real(wp), parameter :: p01 = 0.01_wp, p1 = 0.1_wp
      real(wp), parameter :: big = 1.0E19_wp, tol2 = 5.0E-2_wp
      real(wp) :: fd, h, hc, h1, hc1, pvpstp, stp0
      integer :: i

      ! Variable Definitions (alphabetically)
      !  BETA:     The function parameters.
      !  BIG:      A big value, used to initialize DIFFJ.
      !  D:        The derivative with respect to the Jth unknown parameter.
      !  DIFFJ:    The relative differences between the user supplied and finite difference
      !            derivatives for the derivative being checked.
      !  EPSMAC:   The value of machine precision.
      !  ETA:      The relative noise in the function results.
      !  FCN:      The user supplied subroutine for evaluating the model.
      !  FD:       The forward difference derivative wrt the Jth parameter.
      !  H:        The relative step size for forward differences.
      !  H0:       The initial relative step size for forward differences.
      !  H1:       The default relative step size for forward differences.
      !  HC:       The relative step size for central differences.
      !  HC0:      The initial relative step size for central differences.
      !  HC1:      The default relative step size for central differences.
      !  IFIXB:    The values designating whether the elements of BETA are fixed at their input
      !            values or not.
      !  IFIXX:    The values designating whether the elements of X are fixed at their input values or not.
      !  INTERVAL: Specifies which checks can be performed when checking derivatives based on the
      !            interval of the bound constraints.
      !  ISTOP:    The variable designating whether there are problems computing the function at
      !            the current BETA and DELTA.
      !  ISWRTB:   The variable designating whether the derivatives wrt BETA (ISWRTB=TRUE) or
      !            DELTAS (ISWRTB=FALSE) are being checked.
      !  J:        The index of the partial derivative being examined.
      !  LDIFX:    The leading dimension of array IFIXX.
      !  LQ:       The response currently being examined.
      !  M:        The number of columns of data in the explanatory variable.
      !  MSG:      The error checking results.
      !  MSG1:     The error checking results summary.
      !  N:        The number of observations.
      !  NFEV:     The number of function evaluations.
      !  NP:       The number of function parameters.
      !  NQ:       The number of responses per observation.
      !  NROW:     The row number of the explanatory variable array at which the derivative is to
      !            be checked.
      !  PV:       The predicted value from the model for row NROW.
      !  PVPSTP:   The predicted value for row NROW of the model using the current parameter
      !            estimates for all but the Jth parameter value, which is BETA(J) + STP0.
      !  P01:     The value 0.01E0_wp.
      !  P1:      The value 0.1E0_wp.
      !  STP0:    The initial step size for the finite difference derivative.
      !  TOL:     The agreement tolerance.
      !  TOL2:    A minimum agreement tolerance.
      !  TYPJ:    The typical size of the J-th unknown BETA or DELTA.
      !  WRK1:    A work array of (N BY M BY NQ) elements.
      !  WRK2:    A work array of (N BY NQ) elements.
      !  WRK6:    A work array of (N BY NP BY NQ) elements.
      !  XPLUSD:  The values of X + DELTA.

      ! Calculate the Jth partial derivative using forward difference quotients and decide if it
      ! agrees with user supplied values

      h1 = sqrt(eta)
      hc1 = eta**(one/three)

      msg(lq, j) = 7
      diffj = big

      do i = 1, 3

         if (i == 1) then
            ! Try initial relative step size
            h = h0
            hc = hc0
         elseif (i == 2) then
            ! Try larger relative step size
            h = max(ten*h1, min(hundred*h0, one))
            hc = max(ten*hc1, min(hundred*hc0, one))
         elseif (i == 3) then
            ! Try smaller relative step size
            h = min(p1*h1, max(p01*h, two*epsmac))
            hc = min(p1*hc1, max(p01*hc, two*epsmac))
         end if

         if (iswrtb) then
            ! Perform computations for derivatives wrt BETA
            stp0 = (h*typj*sign(one, beta(j)) + beta(j)) - beta(j)
            call dpvb(fcn, &
                      n, m, np, nq, &
                      beta, xplusd, ifixb, ifixx, ldifx, &
                      nrow, j, lq, stp0, &
                      istop, nfev, pvpstp, &
                      wrk1, wrk2, wrk6)
         else
            ! Perform computations for derivatives wrt DELTA
            stp0 = (h*typj*sign(one, xplusd(nrow, j)) + xplusd(nrow, j)) - xplusd(nrow, j)
            call dpvd(fcn, &
                      n, m, np, nq, &
                      beta, xplusd, ifixb, ifixx, ldifx, &
                      nrow, j, lq, stp0, &
                      istop, nfev, pvpstp, &
                      wrk1, wrk2, wrk6)
         end if
         if (istop /= 0) then
            return
         end if

         fd = (pvpstp - pv)/stp0

         ! Check for agreement

         ! Numerical and analytic derivatives agree
         if (abs(fd - d) <= tol*abs(d)) then

            ! Set relative difference for derivative checking report
            if ((d == zero) .or. (fd == zero)) then
               diffj = abs(fd - d)
            else
               diffj = abs(fd - d)/abs(d)
            end if

            ! Set message flag
            if (d == zero) then
               ! JTH analytic and numerical derivatives are both zero.
               msg(lq, j) = 1

            else
               ! JTH analytic and numerical derivatives are both nonzero.
               msg(lq, j) = 0
            end if

         else
            ! Numerical and analytic derivatives disagree.  Check why
            if ((d == zero) .or. (fd == zero)) then
               if (interval(j) >= 10 .or. .not. iswrtb) then
                  call djckz(fcn, &
                             n, m, np, nq, &
                             beta, xplusd, ifixb, ifixx, ldifx, &
                             nrow, epsmac, j, lq, iswrtb, &
                             tol, d, fd, typj, pvpstp, stp0, pv, &
                             diffj, msg, istop, nfev, &
                             wrk1, wrk2, wrk6)
               else
                  msg(lq, j) = 8
               end if
            else
               if (interval(j) >= 100 .or. .not. iswrtb) then
                  call djckc(fcn, &
                             n, m, np, nq, &
                             beta, xplusd, ifixb, ifixx, ldifx, &
                             eta, tol, nrow, epsmac, j, lq, hc, iswrtb, &
                             fd, typj, pvpstp, stp0, pv, d, &
                             diffj, msg, istop, nfev, &
                             wrk1, wrk2, wrk6)
               else
                  msg(lq, j) = 8
               end if
            end if

            if (msg(lq, j) <= 2) then
               exit
            end if

         end if
      end do

      ! Set summary flag to indicate questionable results
      if ((msg(lq, j) >= 7) .and. (diffj <= tol2)) then
         msg(lq, j) = 6
      end if
      if ((msg(lq, j) >= 1) .and. (msg(lq, j) <= 6)) then
         msg1 = max(msg1, 1)
      elseif (msg(lq, j) >= 7) then
         msg1 = 2
      end if

   end subroutine djckm

   subroutine djckz &
      (fcn, &
       n, m, np, nq, &
       beta, xplusd, ifixb, ifixx, ldifx, &
       nrow, epsmac, j, lq, iswrtb, &
       tol, d, fd, typj, pvpstp, stp0, pv, &
       diffj, msg, istop, nfev, &
       wrk1, wrk2, wrk6)
   !! Recheck the derivatives in the case where the finite difference derivative disagrees with
   !! the analytic derivative and the analytic derivative is zero.
      ! Adapted from STARPAC subroutine DCKZRO.
      ! Routines Called  DPVB, DPVD
      ! Date Written   860529   (YYMMDD)
      ! Revision Date  920619   (YYMMDD)

      use odrpack_kinds, only: zero, one, two, three

      procedure(fcn_t) :: fcn
         !! The user supplied subroutine for evaluating the model.
      integer, intent(in) :: n
         !! The number of observations.
      integer, intent(in) :: m
         !! The number of columns of data in the explanatory variable.
      integer, intent(in) :: np
         !! The number of function parameters.
      integer, intent(in) :: nq
         !! The number of responses per observation.
      real(wp), intent(inout) :: beta(np)
         !! The function parameters.
      real(wp), intent(inout) :: xplusd(n, m)
         !! The values of `x + delta`.
      integer, intent(in) :: ifixb(np)
         !! The values designating whether the elements of `beta` are fixed at their input values or not.
      integer, intent(in) :: ifixx(ldifx, m)
         !! The values designating whether the elements of `x` are fixed at their input values or not.
      integer, intent(in) :: ldifx
         !! The leading dimension of array `ifixx`.
      integer, intent(in) :: nrow
         !! The row number of the explanatory variable array at which the derivative is to be checked.
      real(wp), intent(in) :: epsmac
         !! The value of machine precision.
      integer, intent(in) :: j
         !! The index of the partial derivative being examined.
      integer, intent(in) :: lq
         !! The response currently being examined.
      logical, intent(in) :: iswrtb
         !! The variable designating whether the derivatives wrt `beta` (`iswrtb = .true.`)
         !! or `delta` (`iswrtb = .false.`) are being checked.
      real(wp), intent(in) :: tol
         !! The agreement tolerance.
      real(wp), intent(in) :: d
         !! The derivative with respect to the `j`-th unknown parameter.
      real(wp), intent(in) :: fd
         !! The forward difference derivative wrt the `j`-th parameter.
      real(wp), intent(in) :: typj
         !! The typical size of the `j`-th unknown `beta` or `delta`.
      real(wp), intent(in) :: pvpstp
         !! The predicted value for row `nrow` of the model using the current parameter estimates
         !! for all but the `j`-th parameter value, which is `beta(j) + stp0`.
      real(wp), intent(in) :: stp0
         !! The initial step size for the finite difference derivative.
      real(wp), intent(in) :: pv
         !! The predicted value from the model for row `nrow`.
      real(wp), intent(out) :: diffj
         !! The relative differences between the user supplied and finite difference derivatives
         !! for the derivative being checked.
      integer, intent(out) :: msg(nq, j)
         !! The error checking results.
      integer, intent(out) :: istop
         !! The variable designating whether there are problems computing the function at the
         !! current `beta` and `delta`.
      integer, intent(inout) :: nfev
         !! The number of function evaluations.
      real(wp), intent(out) :: wrk1(n, m, nq)
         !! A work array of `(n, m, nq)` elements.
      real(wp), intent(out) :: wrk2(n, nq)
         !! A work array of `(n, nq)` elements.
      real(wp), intent(out) :: wrk6(n, np, nq)
         !! A work array of `(n, np, nq)` elements.

      ! Local scalars
      real(wp) :: cd, pvmstp

      ! Variable Definitions (alphabetically)
      !  BETA:    The function parameters.
      !  CD:      The central difference derivative wrt the Jth parameter.
      !  D:       The derivative with respect to the Jth unknown parameter.
      !  DIFFJ:   The relative differences between the user supplied and finite difference
      !           derivatives for the derivative being checked.
      !  EPSMAC:  The value of machine precision.
      !  FCN:     The user supplied subroutine for evaluating the model.
      !  FD:      The forward difference derivative wrt the Jth parameter.
      !  IFIXB:   The values designating whether the elements of BETA are fixed at their input
      !           values or not.
      !  IFIXX:   The values designating whether the elements of X are fixed at their input values
      !           or not.
      !  ISTOP:   The variable designating whether there are problems computing the function at the
      !           current BETA and DELTA.
      !  ISWRTB:  The variable designating whether the derivatives wrt BETA (ISWRTB=TRUE) or
      !           X (ISWRTB=FALSE) are being checked.
      !  J:       The index of the partial derivative being examined.
      !  LDIFX:   The leading dimension of array IFIXX.
      !  LQ:      The response currently being examined.
      !  M:       The number of columns of data in the explanatory variable.
      !  MSG:     The error checking results.
      !  N:       The number of observations.
      !  NFEV:    The number of function evaluations.
      !  NP:      The number of function parameters.
      !  NQ:      The number of responses per observation.
      !  NROW:    The row number of the explanatory variable array at which the derivative is to be
      !           checked.
      !  PV:      The predicted value from the model for row NROW.
      !  PVMSTP:  The predicted value for row NROW of the model using the current parameter
      !           estimates for all but the Jth parameter value, which is BETA(J) - STP0.
      !  PVPSTP:  The predicted value for row NROW of the model using the current parameter
      !           estimates for all but the JTH parameter value, which is BETA(J) + STP0.
      !  STP0:    The initial step size for the finite difference derivative.
      !  TOL:     The agreement tolerance.
      !  TYPJ:    The typical size of the J-th unknown BETA or DELTA.
      !  WRK1:    A work array of (N BY M BY NQ) elements.
      !  WRK2:    A work array of (N BY NQ) elements.
      !  WRK6:    A work array of (N BY NP BY NQ) elements.
      !  XPLUSD:  The values of X + DELTA.

      ! Recalculate numerical derivative using central difference and step size of 2*STP0
      if (iswrtb) then
         ! Perform computations for derivatives wrt BETA
         call dpvb(fcn, &
                   n, m, np, nq, &
                   beta, xplusd, ifixb, ifixx, ldifx, &
                   nrow, j, lq, -stp0, &
                   istop, nfev, pvmstp, &
                   wrk1, wrk2, wrk6)
      else
         ! Perform computations for derivatives wrt DELTA
         call dpvd(fcn, &
                   n, m, np, nq, &
                   beta, xplusd, ifixb, ifixx, ldifx, &
                   nrow, j, lq, -stp0, &
                   istop, nfev, pvmstp, &
                   wrk1, wrk2, wrk6)
      end if

      if (istop /= 0) then
         return
      end if

      cd = (pvpstp - pvmstp)/(two*stp0)
      diffj = min(abs(cd - d), abs(fd - d))

      ! Check for agreement
      if (diffj <= tol*abs(d)) then
         ! Finite difference and analytic derivatives now agree
         if (d == zero) then
            msg(lq, j) = 1
         else
            msg(lq, j) = 0
         end if
      elseif (diffj*typj <= abs(pv*epsmac**(one/three))) then
         ! Derivatives are both close to zero
         msg(lq, j) = 2
      else
         ! Derivatives are not both close to zero
         msg(lq, j) = 3
      end if

   end subroutine djckz

   pure subroutine dodchk &
      (n, m, np, nq, &
       isodr, anajac, implct, &
       beta, ifixb, &
       ldx, ldifx, ldscld, ldstpd, ldwe, ld2we, ldwd, ld2wd, &
       ldy, &
       lwork, lwkmn, liwork, liwkmn, &
       sclb, scld, stpb, stpd, &
       info, &
       lower, upper)
   !! Check input parameters, indicating errors found using nonzero values of argument `info`.
      ! Routines Called  (None)
      ! Date Written   860529   (YYMMDD)
      ! Revision Date  920619   (YYMMDD)

      use odrpack_kinds, only: zero

      integer, intent(in) :: n
         !! The number of observations.
      integer, intent(in) :: m
         !! The number of columns of data in the explanatory variable.
      integer, intent(in) :: np
         !! The number of function parameters.
      integer, intent(in) :: nq
         !! The number of responses per observation.
      logical, intent(in) :: isodr
         !! The variable designating whether the solution is by ODR (`isodr = .true.`) or
         !! by OLS (`isodr = .false.`).
      logical, intent(in) :: anajac
         !! The variable designating whether the Jacobians are computed by finite differences
         !! (`anajac = .false.`) or not (`anajac = .true.`).
      logical, intent(in) :: implct
         !! The variable designating whether the solution is by implicit ODR (`implct = .true.`)
         !! or explicit ODR (`implct = .false.`).
      real(wp), intent(in) :: beta(np)
         !! The function parameters.
      integer, intent(in) :: ifixb(np)
         !! The values designating whether the elements of `beta` are fixed at their input values or not.
      integer, intent(in) :: ldx
         !! The leading dimension of array `x`.
      integer, intent(in) :: ldifx
         !! The leading dimension of array `ifixx`.
      integer, intent(in) :: ldscld
         !! The leading dimension of array `scld`.
      integer, intent(in) :: ldstpd
         !! The leading dimension of array `stpd`.
      integer, intent(in) :: ldwe
      !! The leading dimension of array `we`.
      integer, intent(in) :: ld2we
         !! The second dimension of array `we`.
      integer, intent(in) :: ldwd
         !! The leading dimension of array `wd`.
      integer, intent(in) :: ld2wd
         !! The second dimension of array `wd`.
      integer, intent(in) :: ldy
         !! The leading dimension of array `y`.
      integer, intent(in) :: lwork
         !! The length of vector `work`.
      integer, intent(in) :: lwkmn
         !! The minimum acceptable length of array `work`.
      integer, intent(in) :: liwork
         !! The length of vector `iwork`.
      integer, intent(in) :: liwkmn
         !! The minimum acceptable length of array `iwork`.
      real(wp), intent(in) :: sclb(np)
         !! The scaling values for `beta`.
      real(wp), intent(in) :: scld(ldscld, m)
         !! The scaling value for `delta`.
      real(wp), intent(in) :: stpb(np)
         !! The step for the finite difference derivative wrt `beta`.
      real(wp), intent(in) :: stpd(ldstpd, m)
         !! The step for the finite difference derivative wrt `delta`.
      integer, intent(out) :: info
         !! The variable designating why the computations were stopped.
      real(wp), intent(in) :: lower(np)
         !! The lower bound on `beta`.
      real(wp), intent(in) :: upper(np)
         !! The upper bound on `beta`.

      ! Local scalars
      integer :: last, npp

      ! Variable Definitions (alphabetically)
      !  ANAJAC:  The variable designating whether the Jacobians are computed by finite
      !           differences (ANAJAC=FALSE) or not (ANAJAC=TRUE).
      !  I:       An indexing variable.
      !  IFIXB:   The values designating whether the elements of BETA are fixed at their input
      !           values or not.
      !  IMPLCT:  The variable designating whether the solution is by implicit ODR (IMPLCT=TRUE)
      !           or explicit ODR (IMPLCT=FALSE).
      !  INFO:    The variable designating why the computations were stopped.
      !  ISODR:   The variable designating whether the solution is by ODR (ISODR=TRUE) or
      !           by OLS (ISODR=FALSE).
      !  J:       An indexing variable.
      !  K:       An indexing variable.
      !  LAST:    The last row of the array to be accessed.
      !  LDIFX:   The leading dimension of array IFIXX.
      !  LDSCLD:  The leading dimension of array SCLD.
      !  LDSTPD:  The leading dimension of array STPD.
      !  LDWD:    The leading dimension of array WD.
      !  LDWE:    The leading dimension of array WE.
      !  LDX:     The leading dimension of array X.
      !  LDY:     The leading dimension of array X.
      !  LD2WD:   The second dimension of array WD.
      !  LD2WE:   The second dimension of array WE.
      !  LIWKMN:  The minimum acceptable length of array IWORK.
      !  LIWORK:  The length of vector IWORK.
      !  LWKMN:   The minimum acceptable length of array WORK.
      !  LWORK:   The length of vector WORK.
      !  M:       The number of columns of data in the explanatory variable.
      !  N:       The number of observations.
      !  NP:      The number of function parameters.
      !  NPP:     The number of function parameters being estimated.
      !  NQ:      The number of responses per observations.
      !  SCLB:    The scaling values for BETA.
      !  SCLD:    The scaling value for DELTA.
      !  STPB:    The step for the finite difference derivitive wrt BETA.
      !  STPD:    The step for the finite difference derivitive wrt DELTA.

      ! Find actual number of parameters being estimated
      if ((np <= 0) .or. (ifixb(1) < 0)) then
         npp = np
      else
         npp = count(ifixb(1:np) /= 0)
      end if

      ! Check problem specification parameters
      if ((n <= 0) .or. (m <= 0) .or. (npp <= 0 .or. npp > n) .or. (nq <= 0)) then
         info = 10000
         if (n <= 0) then
            info = info + 1000
         end if
         if (m <= 0) then
            info = info + 100
         end if
         if (npp <= 0 .or. npp > n) then
            info = info + 10
         end if
         if (nq <= 0) then
            info = info + 1
         end if
         return
      end if

      ! Check dimension specification parameters
      if ((.not. implct .and. (ldy < n)) .or. &
          (ldx < n) .or. &
          ((ldwe /= 1) .and. (ldwe < n)) .or. &
          ((ld2we /= 1) .and. (ld2we < nq)) .or. &
          (isodr .and. ((ldwd /= 1) .and. (ldwd < n))) .or. &
          (isodr .and. ((ld2wd /= 1) .and. (ld2wd < m))) .or. &
          (isodr .and. ((ldifx /= 1) .and. (ldifx < n))) .or. &
          (isodr .and. ((ldstpd /= 1) .and. (ldstpd < n))) .or. &
          (isodr .and. ((ldscld /= 1) .and. (ldscld < n))) .or. &
          (lwork < lwkmn) .or. &
          (liwork < liwkmn)) then

         info = 20000

         if (.not. implct .and. ldy < n) then
            info = info + 1000
         end if

         if (ldx < n) then
            info = info + 2000
         end if

         if ((ldwe /= 1 .and. ldwe < n) .or. (ld2we /= 1 .and. ld2we < nq)) then
            info = info + 100
         end if

         if (isodr .and. &
             ((ldwd /= 1 .and. ldwd < n) .or. (ld2wd /= 1 .and. ld2wd < m))) then
            info = info + 200
         end if

         if (isodr .and. (ldifx /= 1 .and. ldifx < n)) then
            info = info + 10
         end if

         if (isodr .and. (ldstpd /= 1 .and. ldstpd < n)) then
            info = info + 20
         end if

         if (isodr .and. &
             (ldscld /= 1 .and. ldscld < n)) then
            info = info + 40
         end if

         if (lwork < lwkmn) then
            info = info + 1
         end if

         if (liwork < liwkmn) then
            info = info + 2
         end if

      end if

      ! Check DELTA scaling
      if (isodr .and. scld(1, 1) > zero) then
         if (ldscld >= n) then
            last = n
         else
            last = 1
         end if
         if (any(scld(1:last, 1:m) <= zero)) then
            info = 30200
         end if
      end if

      ! Check BETA scaling
      if (sclb(1) > zero) then
         if (any(sclb(1:np) <= zero)) then
            if (info == 0) then
               info = 30100
            else
               info = info + 100
            end if
         end if
      end if

      ! Check DELTA finite difference step sizes
      if (anajac .and. isodr .and. stpd(1, 1) > zero) then
         if (ldstpd >= n) then
            last = n
         else
            last = 1
         end if
         if (any(stpd(1:last, 1:m) <= zero)) then
            if (info == 0) then
               info = 32000
            else
               info = info + 2000
            end if
         end if
      end if

      ! Check BETA finite difference step sizes
      if (anajac .and. stpb(1) > zero) then
         if (any(stpb(1:np) <= zero)) then
            if (info == 0) then
               info = 31000
            else
               info = info + 1000
            end if
         end if
      end if

      !  Check bounds
      if (any(upper(1:np) < lower(1:np))) then
         if (info == 0) then
            info = 91000
         end if
      end if

      if (any( &
          ((upper(1:np) < beta(1:np)) .or. (lower(1:np) > beta(1:np))) &
          .and. .not. (upper(1:np) < lower(1:np)))) then
         if (info >= 90000) then
            info = info + 100
         else
            info = 90100
         end if
      end if

   end subroutine dodchk

   subroutine dodstp &
      (n, m, np, nq, npp, &
       f, fjacb, fjacd, &
       wd, ldwd, ld2wd, ss, tt, ldtt, delta, &
       alpha, epsfcn, isodr, &
       tfjacb, omega, u, qraux, kpvt, &
       s, t, phi, irank, rcond, forvcv, &
       wrk1, wrk2, wrk3, wrk4, wrk5, wrk, lwrk, tempret, istopc)
   !! Compute locally constrained steps `s` and `t`, and `phi(alpha)`. 
      ! @note: This is one of the most time-consuming subroutines in ODRPACK (~25% of total).
      ! Routines Called  IDAMAX, DCHEX, DESUBI, DFCTR, DNRM2, DQRDC, DQRSL, DROT,
      !                  DROTG, DSOLVE, DTRCO, DTRSL, DVEVTR, DWGHT
      ! Date Written   860529   (YYMMDD)
      ! Revision Date  920619   (YYMMDD)

      use odrpack_kinds, only: zero, one

      integer, intent(in) :: n
         !! The number of observations.
      integer, intent(in) :: m
         !! The number of columns of data in the explanatory variable.
      integer, intent(in) :: np
         !! The number of function parameters.
      integer, intent(in) :: nq
         !! The number of responses per observation.
      integer, intent(in) :: npp
         !! The number of function parameters being estimated.
      real(wp), intent(in) :: f(n, nq)
         !! The (weighted) estimated values of `epsilon`.
      real(wp), intent(in) :: fjacb(n, np, nq)
         !! The Jacobian with respect to `beta`.
      real(wp), intent(in) :: fjacd(n, m, nq)
         !! The Jacobian with respect to `delta`.
      real(wp), intent(in) :: wd(ldwd, ld2wd, m)
         !! The (squared) `delta` weights.
      integer, intent(in) :: ldwd
         !! The leading dimension of array `wd`.
      integer, intent(in) :: ld2wd
         !! The second dimension of array `wd`.
      real(wp), intent(in) :: ss(np)
         !! The scaling values for the unfixed `beta`s.
      real(wp), intent(in) :: tt(ldtt, m)
         !! The scaling values for `delta`.
      integer, intent(in) :: ldtt
         !! The leading dimension of array `tt`.
      real(wp), intent(in) :: delta(n, m)
         !! The estimated errors in the explanatory variables.
      real(wp), intent(in) :: alpha
         !! The Levenberg-Marquardt parameter.
      real(wp), intent(in) :: epsfcn
         !! The function's precision.
      logical, intent(in) :: isodr
         !! The variable designating whether the solution is by ODR (`isodr = .true.`) or
         !! by OLS (`isodr = .false.`).
      real(wp), intent(out) :: tfjacb(n, nq, np)
         !! The array `omega*fjacb`.
      real(wp), intent(out) :: omega(nq, nq)
         !! The array defined such that:
         !! `omega*trans(omega) = inv(I + fjacd*inv(e)*trans(fjacd))
         !! = (I - fjacd*inv(p)*trans(fjacd))`
         !! where `e = d**2 + alpha*tt**2` and
         !! `p = trans(fjacd)*fjacd + d**2 + alpha*tt**2`.
      real(wp), intent(out) :: u(np)
         !! The approximate null vector for `tfjacb`.
      real(wp), intent(out) :: qraux(np)
         !! The array required to recover the orthogonal part of the Q-R decomposition.
      integer, intent(out) :: kpvt(np)
         !! The pivot vector.
      real(wp), intent(out) :: s(np)
         !! The step for `beta`.
      real(wp), intent(out) :: t(n, m)
         !! The step for `delta`.
      real(wp), intent(out) :: phi
         !! The difference between the norm of the scaled step and the trust region diameter.
      integer, intent(out) :: irank
         !! The rank deficiency of the Jacobian wrt `beta`.
      real(wp), intent(out) :: rcond
         !! The approximate reciprocal condition number of `tfjacb`.
      logical, intent(in) :: forvcv
         !! The variable designating whether this subroutine was called to set up for the
         !! covariance matrix computations (`forvcv = .true.`) or not (`forvcv = .false.`).
      real(wp), intent(out) :: wrk1(n, nq, m)
         !! A work array of `(n, nq, m)` elements.
      real(wp), intent(out) :: wrk2(n, nq)
         !! A work array of `(n, nq)` elements.
      real(wp), intent(out) :: wrk3(np)
         !! A work array of `(np)` elements.
      real(wp), intent(out) :: wrk4(m, m)
         !! A work array of `(m, m)` elements.
      real(wp), intent(out) :: wrk5(m)
         !! A work array of `(m)` elements.
      real(wp), intent(out) :: wrk(lwrk)
         !! A work array of `(lwrk)` elements, _equivalenced_ to `wrk1` and `wrk2`.
      integer, intent(in) :: lwrk
         !! The length of vector `wrk`.
      real(wp), intent(inout) :: tempret(:, :)
         !! Temporary work array for holding return values before copying to a lower rank array.
      integer, intent(inout) :: istopc
         !! The variable designating whether the computations were stopped due to a numerical
         !! error within subroutine `dodstp`.

      ! Local scalars
      real(wp) :: co, si, temp
      integer :: i, imax, inf, ipvt, j, k, k1, k2, kp, l
      logical :: elim

      ! Local arrays
      real(wp) :: dum(2)

      ! External BLAS/LINPACK procedures
      real(wp), external :: dnrm2
      integer, external :: idamax
      external :: dchex, dqrdc, dqrsl, drot, drotg, dtrco, dtrsl

      ! Variable definitions (alphabetically)
      !  ALPHA:   The Levenberg-Marquardt parameter.
      !  CO:      The cosine from the plane rotation.
      !  DELTA:   The estimated errors in the explanatory variables.
      !  DUM:     A dummy array.
      !  ELIM:    The variable designating whether columns of the Jacobian
      !           wrt BETA have been eliminated (ELIM=TRUE) or not (ELIM=FALSE).
      !  EPSFCN:  The function's precision.
      !  F:       The (weighted) estimated values of EPSILON.
      !  FJACB:   The Jacobian with respect to BETA.
      !  FJACD:   The Jacobian with respect to DELTA.
      !  FORVCV:  The variable designating whether this subroutine was called to set up for the
      !           covariance matrix computations (FORVCV=TRUE) or not (FORVCV=FALSE).
      !  I:       An indexing variable.
      !  IMAX:    The index of the element of U having the largest absolute value.
      !  INF:     The return code from LINPACK routines.
      !  IPVT:    The variable designating whether pivoting is to be done.
      !  IRANK:   The rank deficiency of the Jacobian wrt BETA.
      !  ISODR:   The variable designating whether the solution is by ODR (ISODR=TRUE) or by
      !           OLS (ISODR=FALSE).
      !  ISTOPC:  The variable designating whether the computations were stoped due to a numerical
      !           error within subroutine DODSTP.
      !  J:       An indexing variable.
      !  K:       An indexing variable.
      !  K1:      An indexing variable.
      !  K2:      An indexing variable.
      !  KP:      The rank of the Jacobian wrt BETA.
      !  KPVT:    The pivot vector.
      !  L:       An indexing variable.
      !  LDTT:    The leading dimension of array TT.
      !  LDWD:    The leading dimension of array WD.
      !  LD2WD:   The second dimension of array WD.
      !  LWRK:    The length of vector WRK.
      !  M:       The number of columns of data in the explanatory variable.
      !  N:       The number of observations.
      !  NP:      The number of function parameters.
      !  NPP:     The number of function parameters being estimated.
      !  OMEGA:   The array defined S.T.
      !           OMEGA*trans(OMEGA) = inv(I+FJACD*inv(E)*trans(FJACD))
      !           = (I-FJACD*inv(P)*trans(FJACD))
      !           where E = D**2 + ALPHA*TT**2
      !           P = trans(FJACD)*FJACD + D**2 + ALPHA*TT**2
      !  PHI:     The difference between the norm of the scaled step and the trust region diameter.
      !  QRAUX:   The array required to recover the orthogonal part of the Q-R decomposition.
      !  RCOND:   The approximate reciprocal condition number of TFJACB.
      !  S:       The step for BETA.
      !  SI:      The sine from the plane rotation.
      !  SS:      The scaling values for the unfixed BETAS.
      !  T:       The step for DELTA.
      !  TEMP:    A temporary storage LOCATION.
      !  TFJACB:  The array OMEGA*FJACB.
      !  TT:      The scaling values for DELTA.
      !  U:       The approximate null vector for TFJACB.
      !  WD:      The (squared) DELTA weights.
      !  WRK:     A work array of (LWRK) elements, equivalenced to WRK1 and WRK2.
      !  WRK1:    A work array of (N by NQ by M) elements.
      !  WRK2:    A work array of (N by NQ) elements.
      !  WRK3:    A work array of (NP) elements.
      !  WRK4:    A work array of (M by M) elements.
      !  WRK5:    A work array of (M) elements.

      ! Compute loop parameters which depend on weight structure

      !  Set up KPVT if ALPHA = 0
      if (alpha == zero) then
         kp = npp
         do k = 1, np
            kpvt(k) = k
         end do
      else
         if (npp >= 1) then
            kp = npp - irank
         else
            kp = npp
         end if
      end if

      if (isodr) then
         ! T = WD * DELTA = D*G2
         call dwght(n, m, wd, ldwd, ld2wd, delta, t)

         do i = 1, n
            !  Compute WRK4, such that TRANS(WRK4)*WRK4 = E = (D**2 + ALPHA*TT**2)
            call desubi(n, m, wd, ldwd, ld2wd, alpha, tt, ldtt, i, wrk4)
            call dfctr(.false., wrk4, m, m, inf)
            if (inf /= 0) then
               istopc = 60000
               return
            end if
            ! Compute OMEGA, such that
            ! trans(OMEGA)*OMEGA = I+FJACD*inv(E)*trans(FJACD)
            ! inv(trans(OMEGA)*OMEGA) = I-FJACD*inv(P)*trans(FJACD)
            call dvevtr(m, nq, i, fjacd, n, m, wrk4, m, wrk1, n, nq, omega, nq, wrk5)
            do l = 1, nq
               omega(l, l) = one + omega(l, l)
            end do
            call dfctr(.false., omega, nq, nq, inf)
            if (inf /= 0) then
               istopc = 60000
               return
            end if
            ! Compute WRK1 = trans(FJACD)*(I-FJACD*inv(P)*trans(JFACD))
            !              = trans(FJACD)*inv(trans(OMEGA)*OMEGA)
            do j = 1, m
               do l = 1, nq
                  wrk1(i, l, j) = fjacd(i, j, l)
               end do
               call dsolve(nq, omega, nq, wrk1(i, 1:nq, j), 4)
               call dsolve(nq, omega, nq, wrk1(i, 1:nq, j), 2)
            end do

            ! Compute WRK5 = inv(E)*D*G2
            do j = 1, m
               wrk5(j) = t(i, j)
            end do
            call dsolve(m, wrk4, m, wrk5, 4)
            call dsolve(m, wrk4, m, wrk5, 2)

            ! Compute TFJACB = inv(trans(OMEGA))*FJACB
            do k = 1, kp
               do l = 1, nq
                  tfjacb(i, l, k) = fjacb(i, kpvt(k), l)
               end do
               call dsolve(nq, omega, nq, tfjacb(i, 1:nq, k), 4)
               do l = 1, nq
                  if (ss(1) > zero) then
                     tfjacb(i, l, k) = tfjacb(i, l, k)/ss(kpvt(k))
                  else
                     tfjacb(i, l, k) = tfjacb(i, l, k)/abs(ss(1))
                  end if
               end do
            end do
            ! Compute WRK2 = (V*inv(E)*D**2*G2 - G1)
            do l = 1, nq
               wrk2(i, l) = zero
               do j = 1, m
                  wrk2(i, l) = wrk2(i, l) + fjacd(i, j, l)*wrk5(j)
               end do
               wrk2(i, l) = wrk2(i, l) - f(i, l)
            end do

            ! Compute WRK2 = inv(trans(OMEGA))*(V*inv(E)*D**2*G2 - G1)
            call dsolve(nq, omega, nq, wrk2(i, 1:nq), 4)
         end do

      else
         do i = 1, n
            do l = 1, nq
               do k = 1, kp
                  tfjacb(i, l, k) = fjacb(i, kpvt(k), l)
                  if (ss(1) > zero) then
                     tfjacb(i, l, k) = tfjacb(i, l, k)/ss(kpvt(k))
                  else
                     tfjacb(i, l, k) = tfjacb(i, l, k)/abs(ss(1))
                  end if
               end do
               wrk2(i, l) = -f(i, l)
            end do
         end do
      end if

      ! Compute S
      ! Do QR factorization (with column pivoting of TFJACB if ALPHA = 0)
      if (alpha == zero) then
         ipvt = 1
         do k = 1, np
            kpvt(k) = 0
         end do
      else
         ipvt = 0
      end if

      call dqrdc(tfjacb, n*nq, n*nq, kp, qraux, kpvt, wrk3, ipvt)
      call dqrsl(tfjacb, n*nq, n*nq, kp, qraux, wrk2, dum, wrk2, dum, dum, dum, 1000, inf)
      if (inf /= 0) then
         istopc = 60000
         return
      end if

      ! Eliminate alpha part using givens rotations
      if (alpha /= zero) then
         s(1:npp) = zero
         do k1 = 1, kp
            wrk3(1:kp) = zero
            wrk3(k1) = sqrt(alpha)
            do k2 = k1, kp
               call drotg(tfjacb(k2, 1, k2), wrk3(k2), co, si)
               if (kp - k2 >= 1) then
                  call drot(kp - k2, tfjacb(k2, 1, k2 + 1), n*nq, wrk3(k2 + 1), 1, co, si)
               end if
               temp = co*wrk2(k2, 1) + si*s(kpvt(k1))
               s(kpvt(k1)) = -si*wrk2(k2, 1) + co*s(kpvt(k1))
               wrk2(k2, 1) = temp
            end do
         end do
      end if

      ! Compute solution - eliminate variables if necessary
      if (npp >= 1) then
         if (alpha == zero) then
            kp = npp
            elim = .true.
            do while (elim .and. kp >= 1)
               ! Estimate RCOND - U will contain approx null vector
               call dtrco(tfjacb, n*nq, kp, rcond, u, 1)
               if (rcond <= epsfcn) then
                  elim = .true.
                  imax = idamax(kp, u, 1)
                  ! IMAX is the column to remove - use DCHEX and fix KPVT
                  if (imax /= kp) then
                     call dchex(tfjacb, n*nq, kp, imax, kp, wrk2, n*nq, 1, qraux, wrk3, 2)
                     k = kpvt(imax)
                     do i = imax, kp - 1
                        kpvt(i) = kpvt(i + 1)
                     end do
                     kpvt(kp) = k
                  end if
                  kp = kp - 1
               else
                  elim = .false.
               end if
            end do
            irank = npp - kp
         end if
      end if

      if (forvcv) return

      ! Backsolve and unscramble
      if (npp >= 1) then
         do i = kp + 1, npp
            wrk2(i, 1) = zero
         end do
         if (kp >= 1) then
            call dtrsl(tfjacb, n*nq, kp, wrk2, 01, inf)
            if (inf /= 0) then
               istopc = 60000
               return
            end if
         end if
         do i = 1, npp
            if (ss(1) > zero) then
               s(kpvt(i)) = wrk2(i, 1)/ss(kpvt(i))
            else
               s(kpvt(i)) = wrk2(i, 1)/abs(ss(1))
            end if
         end do
      end if

      if (isodr) then
         !  NOTE: T and WRK1 have been initialized above,
         !          where T    = WD * DELTA = D*G2
         !          WRK1 = trans(FJACD)*(I-FJACD*inv(P)*trans(JFACD))
         do i = 1, n
            ! Compute WRK4, such that trans(WRK4)*WRK4 = E = (D**2 + ALPHA*TT**2)
            call desubi(n, m, wd, ldwd, ld2wd, alpha, tt, ldtt, i, wrk4)
            call dfctr(.false., wrk4, m, m, inf)
            if (inf /= 0) then
               istopc = 60000
               return
            end if

            ! Compute WRK5 = inv(E)*D*G2
            do j = 1, m
               wrk5(j) = t(i, j)
            end do
            call dsolve(m, wrk4, m, wrk5, 4)
            call dsolve(m, wrk4, m, wrk5, 2)

            do l = 1, nq
               wrk2(i, l) = f(i, l)
               do k = 1, npp
                  wrk2(i, l) = wrk2(i, l) + fjacb(i, k, l)*s(k)
               end do
               do j = 1, m
                  wrk2(i, l) = wrk2(i, l) - fjacd(i, j, l)*wrk5(j)
               end do
            end do

            do j = 1, m
               wrk5(j) = zero
               do l = 1, nq
                  wrk5(j) = wrk5(j) + wrk1(i, l, j)*wrk2(i, l)
               end do
               t(i, j) = -(wrk5(j) + t(i, j))
            end do
            call dsolve(m, wrk4, m, t(i, 1:m), 4)
            call dsolve(m, wrk4, m, t(i, 1:m), 2)
         end do

      end if

      ! Compute PHI(ALPHA) from scaled S and T
      call dwght(npp, 1, reshape(ss, [npp, 1, 1]), npp, 1, reshape(s, [npp, 1]), tempret(1:npp, 1:1))
      wrk(1:npp) = tempret(1:npp, 1)
      if (isodr) then
         call dwght(n, m, reshape(tt, [ldtt, 1, m]), ldtt, 1, t, tempret(1:n, 1:m))
         wrk(npp + 1:npp + 1 + n*m - 1) = reshape(tempret(1:n, 1:m), [n*m])
         phi = dnrm2(npp + n*m, wrk, 1)
      else
         phi = dnrm2(npp, wrk, 1)
      end if

   end subroutine dodstp

   subroutine dodvcv &
      (n, m, np, nq, npp, &
       f, fjacb, fjacd, &
       wd, ldwd, ld2wd, ssf, ss, tt, ldtt, delta, &
       epsfcn, isodr, &
       vcv, sd, &
       wrk6, omega, u, qraux, jpvt, &
       s, t, irank, rcond, rss, idf, rvar, ifixb, &
       wrk1, wrk2, wrk3, wrk4, wrk5, wrk, lwrk, tempret, istopc)
   !! Compute covariance matrix of estimated parameters.
      ! Routines Called  DPODI, DODSTP
      ! Date Written   901207   (YYMMDD)
      ! Revision Date  920619   (YYMMDD)

      use odrpack_kinds, only: zero

      integer, intent(in) :: n
         !! The number of observations.
      integer, intent(in) :: m
         !! The number of columns of data in the explanatory variable.
      integer, intent(in) :: np
         !! The number of function parameters.
      integer, intent(in) :: nq
         !! The number of responses per observation.
      integer, intent(in) :: npp
         !! The number of function parameters being estimated.
      real(wp), intent(in) :: f(n, nq)
         !! The (weighted) estimated values of `epsilon`.
      real(wp), intent(in) :: fjacb(n, np, nq)
         !! The Jacobian with respect to `beta`.
      real(wp), intent(in) :: fjacd(n, m, nq)
         !! The Jacobian with respect to `delta`.
      real(wp), intent(in) :: wd(ldwd, ld2wd, m)
         !! The `delta` weights.
      integer, intent(in) :: ldwd
         !! The leading dimension of array `wd`.
      integer, intent(in) :: ld2wd
         !! The second dimension of array `wd`.
      real(wp), intent(in) :: ssf(np)
         !! The scaling values used for `beta`.
      real(wp), intent(in) :: ss(np)
         !! The scaling values for the unfixed `beta`s.
      real(wp), intent(in) :: tt(ldtt, m)
         !! The scaling values for `delta`.
      integer, intent(in) :: ldtt
         !! The leading dimension of array `tt`.
      real(wp), intent(in) :: delta(n, m)
         !! The estimated errors in the explanatory variables.
      real(wp), intent(in) :: epsfcn
         !! The function's precision.
      logical, intent(in) :: isodr
         !! The variable designating whether the solution is by ODR (`isodr = .true.`) or
         !! by OLS (`isodr = .false.`).
      real(wp), intent(out) :: vcv(np, np)
         !! The covariance matrix of the estimated `beta`s.
      real(wp), intent(out) :: sd(np)
         !! The standard deviations of the estimated `beta`s.
      real(wp), intent(out) :: wrk6(n*nq, np)
         !! A work array of `(n*nq, np)` elements.
      real(wp), intent(out) :: omega(nq, nq)
         !! The array defined such that `omega*trans(omega) = inv(I + fjacd*inv(e)*trans(fjacd))
         !! = (I - fjacd*inv(p)*trans(fjacd))`.
      real(wp), intent(out) :: u(np)
         !! The approximate null vector for `fjacb`.
      real(wp), intent(out) :: qraux(np)
         !! The array required to recover the orthogonal part of the Q-R decomposition.
      integer, intent(out) :: jpvt(np)
         !! The pivot vector.
      real(wp), intent(out) :: s(np)
         !! The step for `beta`.
      real(wp), intent(out) :: t(n, m)
         !! The step for `delta`.
      integer, intent(out) :: irank
         !! The rank deficiency of the Jacobian wrt `beta`.
      real(wp), intent(out) :: rcond
         !! The approximate reciprocal condition of `fjacb`.
      real(wp), intent(inout) :: rss
         !! The residual sum of squares.
      integer, intent(out) :: idf
         !! The degrees of freedom of the fit, equal to the number of observations with nonzero
         !! weighted derivatives minus the number of parameters being estimated.
      real(wp), intent(out) :: rvar
         !! The residual variance.
      integer, intent(in) :: ifixb(np)
         !! The values designating whether the elements of `beta` are fixed at their input values or not.
      real(wp), intent(out) :: wrk1(n, nq, m)
         !! A work array of `(n, nq, m)` elements.
      real(wp), intent(out) :: wrk2(n, nq)
         !! A work array of `(n, nq)` elements.
      real(wp), intent(out) :: wrk3(np)
         !! A work array of `(np)` elements.
      real(wp), intent(out) :: wrk4(m, m)
         !! A work array of `(m, m)` elements.
      real(wp), intent(out) :: wrk5(m)
         !! A work array of `(m)` elements.
      real(wp), intent(out) :: wrk(lwrk)
         !! A work array of `(lwrk)` elements, _equivalenced_ to `wrk1` and `wrk2`.
      integer, intent(in) :: lwrk
         !! The length of vector `lwrk`.
      real(wp), intent(inout) :: tempret(:, :)
         !! Temporary work array for holding return values before copying to a lower rank array.
      integer, intent(out) :: istopc
         !! The variable designating whether the computations were stoped due to a numerical
         !! error within subroutine `dodstp`.

      ! Local scalars
      real(wp) :: temp
      integer :: i, iunfix, j, junfix, kp
      logical :: forvcv

      ! External BLAS/LINPACK procedures
      external :: dpodi

      ! Variable definitions (alphabetically)
      !  DELTA:   The estimated errors in the explanatory variables.
      !  EPSFCN:  The function's precision.
      !  F:       The (weighted) estimated values of EPSILON.
      !  FJACB:   The Jacobian with respect to BETA.
      !  FJACD:   The Jacobian with respect to DELTA.
      !  FORVCV:  The variable designating whether subroutine DODSTP is called to set up for the
      !           covariance matrix computations (FORVCV=TRUE) or not (FORVCV=FALSE).
      !  I:       An indexing variable.
      !  IDF:     The degrees of freedom of the fit, equal to the number of observations with nonzero
      !           weighted derivatives minus the number of parameters being estimated.
      !  IFIXB:   The values designating whether the elements of BETA are fixed at their input values
      !           or not.
      !  IMAX:    The index of the element of U having the largest absolute value.
      !  IRANK:   The rank deficiency of the Jacobian wrt BETA.
      !  ISODR:   The variable designating whether the solution is by ODR (ISODR=TRUE) or by
      !           OLS (ISODR=FALSE).
      !  ISTOPC:  The variable designating whether the computations were stoped due to a numerical
      !           error within subroutine DODSTP.
      !  IUNFIX:  The index of the next unfixed parameter.
      !  J:       An indexing variable.
      !  JPVT:    The pivot vector.
      !  JUNFIX:  The index of the next unfixed parameter.
      !  KP:      The rank of the Jacobian wrt BETA.
      !  L:       An indexing variable.
      !  LDTT:    The leading dimension of array TT.
      !  LDWD:    The leading dimension of array WD.
      !  LD2WD:   The second dimension of array WD.
      !  LWRK:    The length of vector WRK.
      !  M:       The number of columns of data in the explanatory variable.
      !  N:       The number of observations.
      !  NP:      The number of function parameters.
      !  NPP:     The number of function parameters being estimated.
      !  NQ:      The number of responses per observation.
      !  OMEGA:   The array defined S.T.
      !           OMEGA*trans(OMEGA) = inv(I+FJACD*inv(E)*trans(FJACD))
      !           = (I-FJACD*inv(P)*trans(FJACD))
      !           where E = D**2 + ALPHA*TT**2 and P = trans(FJACD)*FJACD + D**2 + ALPHA*TT**2
      !  QRAUX:   The array required to recover the orthogonal part of the Q-R decomposition.
      !  RCOND:   The approximate reciprocal condition of FJACB.
      !  RSS:     The residual sum of squares.
      !  RVAR:    The residual variance.
      !  S:       The step for BETA.
      !  SD:      The standard deviations of the estimated BETAS.
      !  SS:      The scaling values for the unfixed BETAS.
      !  SSF:     The scaling values used for BETA.
      !  T:       The step for DELTA.
      !  TEMP:    A temporary storage location
      !  TT:      The scaling values for DELTA.
      !  U:       The approximate null vector for FJACB.
      !  VCV:     The covariance matrix of the estimated BETAS.
      !  WD:      The DELTA weights.
      !  WRK:     A work array of (LWRK) elements, equivalenced to WRK1 and WRK2.
      !  WRK1:    A work array of (N by NQ by M) elements.
      !  WRK2:    A work array of (N by NQ) elements.
      !  WRK3:    A work array of (NP) elements.
      !  WRK4:    A work array of (M by M) elements.
      !  WRK5:    A work array of (M) elements.
      !  WRK6:    A work array of (N*NQ by P) elements.

      forvcv = .true.
      istopc = 0

      call dodstp(n, m, np, nq, npp, &
                  f, fjacb, fjacd, &
                  wd, ldwd, ld2wd, ss, tt, ldtt, delta, &
                  zero, epsfcn, isodr, &
                  wrk6, omega, u, qraux, jpvt, &
                  s, t, temp, irank, rcond, forvcv, &
                  wrk1, wrk2, wrk3, wrk4, wrk5, wrk, lwrk, tempret, istopc)
      if (istopc /= 0) then
         return
      end if
      kp = npp - irank
      call dpodi(wrk6, n*nq, kp, wrk3, 1)

      idf = 0
      do i = 1, n
         if (any(fjacb(i, :, :) /= zero)) then
            idf = idf + 1
            cycle
         end if
         if (isodr) then
            if (any(fjacd(i, :, :) /= zero)) then
               idf = idf + 1
               cycle
            end if
         end if
      end do

      if (idf > kp) then
         idf = idf - kp
         rvar = rss/idf
      else
         idf = 0
         rvar = rss
      end if

      ! Store variances in SD, restoring original order
      do i = 1, np
         sd(i) = zero
      end do
      do i = 1, kp
         sd(jpvt(i)) = wrk6(i, i)
      end do
      if (np > npp) then
         junfix = npp
         do j = np, 1, -1
            if (ifixb(j) == 0) then
               sd(j) = zero
            else
               sd(j) = sd(junfix)
               junfix = junfix - 1
            end if
         end do
      end if

      ! Store covariance matrix in VCV, restoring original order
      do i = 1, np
         do j = 1, i
            vcv(i, j) = zero
         end do
      end do
      do i = 1, kp
         do j = i + 1, kp
            if (jpvt(i) > jpvt(j)) then
               vcv(jpvt(i), jpvt(j)) = wrk6(i, j)
            else
               vcv(jpvt(j), jpvt(i)) = wrk6(i, j)
            end if
         end do
      end do
      if (np > npp) then
         iunfix = npp
         do i = np, 1, -1
            if (ifixb(i) == 0) then
               do j = i, 1, -1
                  vcv(i, j) = zero
               end do
            else
               junfix = npp
               do j = np, 1, -1
                  if (ifixb(j) == 0) then
                     vcv(i, j) = zero
                  else
                     vcv(i, j) = vcv(iunfix, junfix)
                     junfix = junfix - 1
                  end if
               end do
               iunfix = iunfix - 1
            end if
         end do
      end if

      do i = 1, np
         vcv(i, i) = sd(i)
         sd(i) = sqrt(rvar*sd(i))
         do j = 1, i
            vcv(j, i) = vcv(i, j)
         end do
      end do

      ! Unscale standard errors and covariance matrix
      do i = 1, np
         if (ssf(1) > zero) then
            sd(i) = sd(i)/ssf(i)
         else
            sd(i) = sd(i)/abs(ssf(1))
         end if
         do j = 1, np
            if (ssf(1) > zero) then
               vcv(i, j) = vcv(i, j)/(ssf(i)*ssf(j))
            else
               vcv(i, j) = vcv(i, j)/(ssf(1)*ssf(1))
            end if
         end do
      end do

   end subroutine dodvcv

   subroutine dpack(n2, n1, v1, v2, ifix)
   !! Select the unfixed elements of `v2` and return them in `v1`.
      ! Routines Called  DCOPY
      ! Date Written   860529   (YYMMDD)
      ! Revision Date  920304   (YYMMDD)

      integer, intent(in) :: n2
         !! The number of items in `v2`.
      integer, intent(out) :: n1
         !! The number of items in `v1`.
      real(wp), intent(out) :: v1(n2)
         !! The vector of the unfixed items from `v2`.
      real(wp), intent(in) :: v2(n2)
         !! The vector of the fixed and unfixed items from which the unfixed elements are to be extracted.
      integer, intent(in) :: ifix(n2)
         !! The values designating whether the elements of `v2` are fixed at their input values or not.

      ! Local scalars
      integer :: i

      ! External subroutines
      external :: dcopy

      ! Variable definitions (alphabetically)
      !  I:       An indexing variable.
      !  IFIX:    The values designating whether the elements of V2 are fixed at their input
      !           values or not.
      !  N1:      The number of items in V1.
      !  N2:      The number of items in V2.
      !  V1:      The vector of the unfixed items from V2.
      !  V2:      The vector of the fixed and unfixed items from which the unfixed elements are
      !           to be extracted.

      n1 = 0
      if (ifix(1) >= 0) then
         do i = 1, n2
            if (ifix(i) /= 0) then
               n1 = n1 + 1
               v1(n1) = v2(i)
            end if
         end do
      else
         n1 = n2
         call dcopy(n2, v2, 1, v1, 1)
      end if

   end subroutine dpack

   real(wp) pure function dppnml(p) result(dppnmlr)
   !! Compute the percent point function value for the normal (Gaussian) distribution with
   !!  mean 0 and standard deviation 1, and with probability density function:
   !!
   !!       `f(x) = (1/sqrt(2*pi))*exp(-x^2/2)`
   !!
      ! Adapted from DATAPAC subroutine TPPF, with modifications to facilitate conversion to
      ! real(wp) automatically.
      ! Routines Called  (NONE)
      ! Date Written   901207   (YYMMDD)
      ! Revision Date  920304   (YYMMDD)
      !***Author  Filliben, James J.,
      !       Statistical Engineering Division
      !       National Bureau of Standards
      !       Washington, D. C. 20234
      !       (Original Version--June      1972.
      !       (Updated         --September 1975,
      !       November  1975, AND
      !       October   1976.
      !***Description
      !       --The coding as presented below is essentially
      !       identical to that presented by Odeh and Evans
      !       as Algortihm 70 of Applied Statistics.
      !       --As pointed out by Odeh and Evans in Applied
      !       Statistics, their algorithm representes a
      !       substantial improvement over the previously employed
      !       Hastings approximation for the normal percent point
      !       function, with accuracy improving from 4.5*(10**-4)
      !       to 1.5*(10**-8).
      !***References  Odeh and Evans, the Percentage Points of the Normal
      !       Distribution, Algortihm 70, Applied Statistics, 1974,
      !       Pages 96-97.
      !       Evans, Algorithms for Minimal Degree Polynomial and
      !       Rational Approximation, M. Sc. Thesis, 1972,
      !       University of Victoria, B. C., Canada.
      !       Hastings, Approximations for Digital Computers, 1955,
      !       Pages 113, 191, 192.
      !       National Bureau of Standards Applied Mathematics
      !       Series 55, 1964, Page 933, Formula 26.2.23.
      !       Filliben, Simple and Robust Linear Estimation of the
      !       Location Parameter of a Symmetric Distribution
      !       (Unpublished Ph.D. Dissertation, Princeton
      !       University), 1969, Pages 21-44, 229-231.
      !       Filliben, "The Percent Point Function",
      !       (Unpublished Manuscript), 1970, Pages 28-31.
      !       Johnson and Kotz, Continuous Univariate Distributions,
      !       Volume 1, 1970, Pages 40-111.
      !       Kelley Statistical Tables, 1948.
      !       Owen, Handbook of Statistical Tables, 1962, Pages 3-16.
      !       Pearson and Hartley, Biometrika Tables for
      !       Statisticians, Volume 1, 1954, Pages 104-113.

      use odrpack_kinds, only: zero, half, one, two

      real(wp), intent(in) :: p
         !! The probability at which the percent point is to be evaluated. `p` must lie between
         !! 0.0 and 1.0, exclusive.

      ! Local scalars
      real(wp), parameter :: p0 = -0.322232431088E0_wp, &
                             p1 = -1.0E0_wp, &
                             p2 = -0.342242088547E0_wp, &
                             p3 = -0.204231210245E-1_wp, &
                             p4 = -0.453642210148E-4_wp, &
                             q0 = 0.993484626060E-1_wp, &
                             q1 = 0.588581570495E0_wp, &
                             q2 = 0.531103462366E0_wp, &
                             q3 = 0.103537752850E0_wp, &
                             q4 = 0.38560700634E-2_wp
      real(wp) :: aden, anum, r, t

      ! Variable Definitions (alphabetically)
      !  ADEN:    A value used in the approximation.
      !  ANUM:    A value used in the approximation.
      !  P:       The probability at which the percent point is to be evaluated. P must be between
      !           0.0E0_wp and 1.0E0_wp, exclusive.
      !  P0:      A parameter used in the approximation.
      !  P1:      A parameter used in the approximation.
      !  P2:      A parameter used in the approximation.
      !  P3:      A parameter used in the approximation.
      !  P4:      A parameter used in the approximation.
      !  Q0:      A parameter used in the approximation.
      !  Q1:      A parameter used in the approximation.
      !  Q2:      A parameter used in the approximation.
      !  Q3:      A parameter used in the approximation.
      !  Q4:      A parameter used in the approximation.
      !  R:       The probability at which the percent point is evaluated.
      !  T:       A value used in the approximation.

      if (p == half) then
         dppnmlr = zero
      else
         r = p
         if (p > half) r = one - r
         t = sqrt(-two*log(r))
         anum = ((((t*p4 + p3)*t + p2)*t + p1)*t + p0)
         aden = ((((t*q4 + q3)*t + q2)*t + q1)*t + q0)
         dppnmlr = t + (anum/aden)
         if (p < half) dppnmlr = -dppnmlr
      end if

   end function dppnml

   real(wp) pure function dppt(p, idf) result(dpptr)
   !! Compute the percent point function value for the student's T distribution with `idf`
   !! degrees of freedom.
      ! Adapted from DATAPAC subroutine TPPF, with modifications to facilitate conversion to
      ! real(wp) automatically.
      ! Routines Called  DPPNML
      ! Date Written   901207   (YYMMDD)
      ! Revision Date  920304   (YYMMDD)
      !***Author  Filliben, James J.,
      !       Statistical Engineering Division
      !       National Bureau of Standards
      !       Washington, D. C. 20234
      !       (Original Version--October   1975.)
      !       (Updated         --November  1975.)
      !***Description
      !       --For IDF = 1 AND IDF = 2, the percent point function
      !       for the T distribution exists in simple closed form
      !       and so the computed percent points are exact.
      !       --For IDF between 3 and 6, inclusively, the approximation
      !       is augmented by 3 iterations of Newton's method to
      !       improve the accuracy, especially for P near 0 or 1.
      !***References  National Bureau of Standards Applied Mathmatics
      !       Series 55, 1964, Page 949, Formula 26.7.5.
      !       Johnson and Kotz, Continuous Univariate Distributions,
      !       Volume 2, 1970, Page 102, Formula 11.
      !       Federighi, "Extended Tables of the Percentage Points
      !       of Student"S T Distribution, Journal of the American
      !       Statistical Association, 1969, Pages 683-688.
      !       Hastings and Peacock, Statistical Distributions, A
      !       Handbook for Students and Practitioners, 1975,
      !       Pages 120-123.

      use odrpack_kinds, only: pi, zero, half, one, two, three, eight, fiftn

      real(wp), intent(in) :: p
         !! The probability at which the percent point is to be evaluated. `p` must lie between
         !! 0.0 and 1.0, exclusive.
      integer, intent(in) :: idf
         !! The (positive integer) degrees of freedom.

      ! Local scalars
      real(wp), parameter :: b21 = 4.0E0_wp, &
                             b31 = 96.0E0_wp, &
                             b32 = 5.0E0_wp, &
                             b33 = 16.0E0_wp, &
                             b34 = 3.0E0_wp, &
                             b41 = 384.0E0_wp, &
                             b42 = 3.0E0_wp, &
                             b43 = 19.0E0_wp, &
                             b44 = 17.0E0_wp, &
                             b45 = -15.0E0_wp, &
                             b51 = 9216.0E0_wp, &
                             b52 = 79.0E0_wp, &
                             b53 = 776.0E0_wp, &
                             b54 = 1482.0E0_wp, &
                             b55 = -1920.0E0_wp, &
                             b56 = -945.0E0_wp
      real(wp) :: arg, c, con, d1, d3, d5, d7, d9, df, ppfn, s, term1, term2, term3, &
                  term4, term5, z
      integer :: ipass, maxit

      ! Variable definitions (alphabetically)
      !  ARG:    A value used in the approximation.
      !  B21:    A parameter used in the approximation.
      !  B31:    A parameter used in the approximation.
      !  B32:    A parameter used in the approximation.
      !  B33:    A parameter used in the approximation.
      !  B34:    A parameter used in the approximation.
      !  B41:    A parameter used in the approximation.
      !  B42:    A parameter used in the approximation.
      !  B43:    A parameter used in the approximation.
      !  B44:    A parameter used in the approximation.
      !  B45:    A parameter used in the approximation.
      !  B51:    A parameter used in the approximation.
      !  B52:    A parameter used in the approximation.
      !  B53:    A parameter used in the approximation.
      !  B54:    A parameter used in the approximation.
      !  B55:    A parameter used in the approximation.
      !  B56:    A parameter used in the approximation.
      !  C:      A value used in the approximation.
      !  CON:    A value used in the approximation.
      !  DF:     The degrees of freedom.
      !  D1:     A value used in the approximation.
      !  D3:     A value used in the approximation.
      !  D5:     A value used in the approximation.
      !  D7:     A value used in the approximation.
      !  D9:     A value used in the approximation.
      !  IDF:    The (positive integer) degrees of freedom.
      !  IPASS:  A value used in the approximation.
      !  MAXIT:  The maximum number of iterations allowed for the approx.
      !  P:      The probability at which the percent point is to be
      !          evaluated.  P must lie between 0.0DO and 1.0E0_wp, exclusive.
      !  PPFN:   The normal percent point value.
      !  S:      A value used in the approximation.
      !  TERM1:  A value used in the approximation.
      !  TERM2:  A value used in the approximation.
      !  TERM3:  A value used in the approximation.
      !  TERM4:  A value used in the approximation.
      !  TERM5:  A value used in the approximation.
      !  Z:      A value used in the approximation.

      df = idf
      maxit = 5

      if (idf <= 0) then
         !Treat the IDF < 1 case
         dpptr = zero

      elseif (idf == 1) then
         !Treat the IDF = 1 (Cauchy) case
         arg = pi*p
         dpptr = -cos(arg)/sin(arg)

      elseif (idf == 2) then
         !  Treat the IDF = 2 case
         term1 = sqrt(two)/two
         term2 = two*p - one
         term3 = sqrt(p*(one - p))
         dpptr = term1*term2/term3

      elseif (idf >= 3) then
         ! Treat the IDF greater than or equal to 3 case
         ppfn = dppnml(p)
         d1 = ppfn
         d3 = ppfn**3
         d5 = ppfn**5
         d7 = ppfn**7
         d9 = ppfn**9
         term1 = d1
         term2 = (one/b21)*(d3 + d1)/df
         term3 = (one/b31)*(b32*d5 + b33*d3 + b34*d1)/(df**2)
         term4 = (one/b41)*(b42*d7 + b43*d5 + b44*d3 + b45*d1)/(df**3)
         term5 = (one/b51)*(b52*d9 + b53*d7 + b54*d5 + b55*d3 + b56*d1)/(df**4)
         dpptr = term1 + term2 + term3 + term4 + term5

         if (idf == 3) then
            ! Augment the results for the IDF = 3 case
            con = pi*(p - half)
            arg = dpptr/sqrt(df)
            z = atan(arg)
            do ipass = 1, maxit
               s = sin(z)
               c = cos(z)
               z = z - (z + s*c - con)/(two*c**2)
            end do
            dpptr = sqrt(df)*s/c

         elseif (idf == 4) then
            ! Augment the results for the IDF = 4 case
            con = two*(p - half)
            arg = dpptr/sqrt(df)
            z = atan(arg)
            do ipass = 1, maxit
               s = sin(z)
               c = cos(z)
               z = z - ((one + half*c**2)*s - con)/((one + half)*c**3)
            end do
            dpptr = sqrt(df)*s/c

         elseif (idf == 5) then
            ! Augment the results for the IDF = 5 case
            con = pi*(p - half)
            arg = dpptr/sqrt(df)
            z = atan(arg)
            do ipass = 1, maxit
               s = sin(z)
               c = cos(z)
               z = z - (z + (c + (two/three)*c**3)*s - con)/ &
                   ((eight/three)*c**4)
            end do
            dpptr = sqrt(df)*s/c

         elseif (idf == 6) then
            !  Augment the results for the IDF = 6 case
            con = two*(p - half)
            arg = dpptr/sqrt(df)
            z = atan(arg)
            do ipass = 1, maxit
               s = sin(z)
               c = cos(z)
               z = z - ((one + half*c**2 + (three/eight)*c**4)*s - con)/ &
                   ((fiftn/eight)*c**5)
            end do
            dpptr = sqrt(df)*s/c
         end if
      end if

   end function dppt

   subroutine dpvb &
      (fcn, &
       n, m, np, nq, &
       beta, xplusd, ifixb, ifixx, ldifx, &
       nrow, j, lq, stp, &
       istop, nfev, pvb, &
       wrk1, wrk2, wrk6)
   !! Compute the `nrow`-th function value using `beta(j) + stp`.
      ! Routines Called  FCN
      ! Date Written   860529   (YYMMDD)
      ! Revision Date  920304   (YYMMDD)

      procedure(fcn_t) :: fcn
         !! The user-supplied subroutine for evaluating the model.
      integer, intent(in) :: n
         !! The number of observations.
      integer, intent(in) :: m
         !! The number of columns of data in the independent variable.
      integer, intent(in) :: np
         !! The number of function parameters.
      integer, intent(in) :: nq
         !! The number of responses per observation.
      real(wp), intent(inout) :: beta(np)
         !! The function parameters.
      real(wp), intent(in) :: xplusd(n, m)
         !! The values of `x + delta`.
      integer, intent(in) :: ifixb(np)
         !! The values designating whether the elements of `beta` are fixed at their input values or not.
      integer, intent(in) :: ifixx(ldifx, m)
         !! The values designating whether the elements of `x` are fixed at their input values or not.
      integer, intent(in) :: ldifx
         !! The leading dimension of array `ifixx`.
      integer, intent(in) :: nrow
         !! The row number of the independent variable array at which the derivative is to be checked.
      integer, intent(in) :: j
         !! The index of the partial derivative being examined.
      integer, intent(in) :: lq
         !! The response currently being examined.
      real(wp), intent(in) :: stp
         !! The step size for the finite difference derivative.
      integer, intent(out) :: istop
         !! The variable designating whether there are problems computing the function at the
         !! current `beta` and `delta`.
      integer, intent(inout) :: nfev
         !! The number of function evaluations.
      real(wp), intent(out) :: pvb
         !! The function value for the selected observation & response.
      real(wp), intent(out) :: wrk1(n, m, nq)
         !! Work array.
      real(wp), intent(out) :: wrk2(n, nq)
         !! Work array.
      real(wp), intent(out) :: wrk6(n, np, nq)
         !! Work array.

      ! Local scalars
      real(wp) :: betaj

      ! Variable Definitions (alphabetically)
      !  BETA:    The function parameters.
      !  BETAJ:   The current estimate of the jth parameter.
      !  FCN:     The user-supplied subroutine for evaluating the model.
      !  IFIXB:   The values designating whether the elements of BETA are fixed at their input
      !           values or not.
      !  IFIXX:   The values designating whether the elements of X are fixed at their input
      !           values or not.
      !  ISTOP:   The variable designating whether there are problems computing the function at
      !           the current BETA and DELTA.
      !  J:       The index of the partial derivative being examined.
      !  LDIFX:   The leading dimension of array IFIXX.
      !  LQ:      The response currently being examined.
      !  M:       The number of columns of data in the independent variable.
      !  N:       The number of observations.
      !  NFEV:    The number of function evaluations.
      !  NP:      The number of function parameters.
      !  NQ:      The number of responses per observation.
      !  NROW:    The row number of the independent variable array at which the derivative is
      !           to be checked.
      !  PVB:     The function value for the selected observation & response.
      !  STP:     The step size for the finite difference derivative.
      !  XPLUSD:  The values of X + DELTA.

      betaj = beta(j)
      beta(j) = beta(j) + stp
      istop = 0
      call fcn(n, m, np, nq, &
               n, m, np, &
               beta, xplusd, &
               ifixb, ifixx, ldifx, &
               003, wrk2, wrk6, wrk1, &
               istop)
      if (istop == 0) then
         nfev = nfev + 1
      else
         return
      end if
      beta(j) = betaj

      pvb = wrk2(nrow, lq)

   end subroutine dpvb

   subroutine dpvd &
      (fcn, &
       n, m, np, nq, &
       beta, xplusd, ifixb, ifixx, ldifx, &
       nrow, j, lq, stp, &
       istop, nfev, pvd, &
       wrk1, wrk2, wrk6)
   !! Compute `nrow`-th function value using `x(nrow, j) + delta(nrow, j) + stp`.
      ! Routines Called FCN
      ! Date Written   860529   (YYMMDD)
      ! Revision Date  920304   (YYMMDD)

      procedure(fcn_t) :: fcn
         !! The user-supplied subroutine for evaluating the model.
      integer, intent(in) :: n
         !! The number of observations.
      integer, intent(in) :: m
         !! The number of columns of data in the independent variable.
      integer, intent(in) :: np
         !! The number of function parameters.
      integer, intent(in) :: nq
         !! The number of responses per observation.
      real(wp), intent(in) :: beta(np)
         !! The function parameters.
      real(wp), intent(inout) :: xplusd(n, m)
         !! The values of `x + delta`.
      integer, intent(in) :: ifixb(np)
         !! The values designating whether the elements of `beta` are fixed at their input values or not.
      integer, intent(in) :: ifixx(ldifx, m)
         !! The values designating whether the elements of `x` are fixed at their input values or not.
      integer, intent(in) :: ldifx
         !! The leading dimension of array `ifixx`.
      integer, intent(in) :: nrow
         !! The row number of the independent variable array at which the derivative is to be checked.
      integer, intent(in) :: j
         !! The index of the partial derivative being examined.
      integer, intent(in) :: lq
         !! The response currently being examined.
      real(wp), intent(in) :: stp
         !! The step size for the finite difference derivative.
      integer, intent(out) :: istop
         !! The variable designating whether there are problems computing the function at the
         !! current `beta` and `delta`.
      integer, intent(inout) :: nfev
         !! The number of function evaluations.
      real(wp), intent(out) :: pvd
         !! The function value for the selected observation & response.
      real(wp), intent(out) :: wrk1(n, m, nq)
         !! Work array.
      real(wp), intent(out) :: wrk2(n, nq)
         !! Work array.
      real(wp), intent(out) :: wrk6(n, np, nq)
         !! Work array.

      ! Local scalars
      real(wp) :: xpdj

      ! Variable Definitions (alphabetically)
      !  BETA:    The function parameters.
      !  FCN:     The user-supplied subroutine for evaluating the model.
      !  IFIXB:   The values designating whether the elements of BETA are
      !           fixed at their input values or not.
      !  IFIXX:   The values designating whether the elements of X are fixed at their input values or not.
      !  ISTOP:   The variable designating whether there are problems
      !           computing the function at the current BETA and DELTA.
      !  J:       The index of the partial derivative being examined.
      !  LDIFX:   The leading dimension of array IFIXX.
      !  LQ:      The response currently being examined.
      !  M:       The number of columns of data in the independent variable.
      !  N:       The number of observations.
      !  NFEV:    The number of function evaluations.
      !  NP:      The number of function parameters.
      !  NQ:      The number of responses per observation.
      !  NROW:    The row number of the independent variable array at which the derivative is to be checked.
      !  PVD:     The function value for the selected observation & response.
      !  STP:     The step size for the finite difference derivative.
      !  XPDJ:    The (NROW,J)th element of XPLUSD.
      !  XPLUSD:  The values of X + DELTA.

      xpdj = xplusd(nrow, j)
      xplusd(nrow, j) = xplusd(nrow, j) + stp
      istop = 0
      call fcn(n, m, np, nq, &
               n, m, np, &
               beta, xplusd, &
               ifixb, ifixx, ldifx, &
               003, wrk2, wrk6, wrk1, &
               istop)
      if (istop == 0) then
         nfev = nfev + 1
      else
         return
      end if
      xplusd(nrow, j) = xpdj

      pvd = wrk2(nrow, lq)

   end subroutine dpvd

   pure subroutine dscale(n, m, scl, ldscl, t, ldt, sclt, ldsclt)
   !! Scale `t` by the inverse of `scl`, i.e., compute `t/scl`.
      ! Routines Called (NONE)
      ! Date Written   860529   (YYMMDD)
      ! Revision Date  920304   (YYMMDD)

      use odrpack_kinds, only: zero, one

      integer, intent(in) :: n
         !! The number of rows of data in `t`.
      integer, intent(in) :: m
         !! The number of columns of data in `t`.
      real(wp), intent(in) :: scl(ldscl, m)
         !! The scale values.
      integer, intent(in) :: ldscl
         !! The leading dimension of array `scl`.
      real(wp), intent(in) :: t(ldt, m)
         !! The array to be inversely scaled by `scl`.
      integer, intent(in) :: ldt
         !! The leading dimension of array `t`.
      real(wp), intent(out) :: sclt(ldsclt, m)
         !! The inversely scaled matrix.
      integer, intent(in) :: ldsclt
         !! The leading dimension of array `sclt`.

      ! Local scalars
      real(wp) :: temp
      integer :: i, j

      ! Variable Definitions (alphabetically)
      !  I:       An indexing variable.
      !  J:       An indexing variable.
      !  LDSCL:   The leading dimension of array SCL.
      !  LDSCLT:  The leading dimension of array SCLT.
      !  LDT:     The leading dimension of array T.
      !  M:       The number of columns of data in T.
      !  N:       The number of rows of data in T.
      !  SCL:     The scale values.
      !  SCLT:    The inversely scaled matrix.
      !  T:       The array to be inversely scaled by SCL.
      !  TEMP:    A temporary scalar.

      if (n == 0 .or. m == 0) return

      if (scl(1, 1) >= zero) then
         if (ldscl >= n) then
            do j = 1, m
               do i = 1, n
                  sclt(i, j) = t(i, j)/scl(i, j)
               end do
            end do
         else
            do j = 1, m
               temp = one/scl(1, j)
               do i = 1, n
                  sclt(i, j) = t(i, j)*temp
               end do
            end do
         end if
      else
         temp = one/abs(scl(1, 1))
         do j = 1, m
            do i = 1, n
               sclt(i, j) = t(i, j)*temp
            end do
         end do
      end if

   end subroutine dscale

   pure subroutine dsclb(np, beta, ssf)
   !! Select scaling values for `beta` according to the algorithm given in the ODRPACK95
   !! reference guide.
      ! Routines Called (NONE)
      ! Date Written   860529   (YYMMDD)
      ! Revision Date  920304   (YYMMDD)

      use odrpack_kinds, only: zero, one, ten

      integer, intent(in) :: np
         !! The number of function parameters.
      real(wp), intent(in) :: beta(np)
         !! The function parameters.
      real(wp), intent(out) :: ssf(np)
         !! The scaling values for `beta`.

      ! Local scalars
      real(wp) :: bmax, bmin
      integer :: k
      logical ::bigdif

      ! Variable Definitions (alphabetically)
      !  BETA:    The function parameters.
      !  BIGDIF:  The variable designating whether there is a significant difference in the
      !           magnitudes of the nonzero elements of BETA (BIGDIF=.TRUE.) or not (BIGDIF=.FALSE.).
      !  BMAX:    The largest nonzero magnitude.
      !  BMIN:    The smallest nonzero magnitude.
      !  K:       An indexing variable.
      !  NP:      The number of function parameters.
      !  SSF:     The scaling values for BETA.

      bmax = abs(beta(1))
      do k = 2, np
         bmax = max(bmax, abs(beta(k)))
      end do

      if (bmax == zero) then
         !  All input values of BETA are zero
         ssf(1:np) = one
      else
         !  Some of the input values are nonzero
         bmin = bmax
         do k = 1, np
            if (beta(k) /= zero) then
               bmin = min(bmin, abs(beta(k)))
            end if
         end do
         bigdif = log10(bmax) - log10(bmin) >= one
         do k = 1, np
            if (beta(k) == zero) then
               ssf(k) = ten/bmin
            else
               if (bigdif) then
                  ssf(k) = one/abs(beta(k))
               else
                  ssf(k) = one/bmax
               end if
            end if
         end do

      end if

   end subroutine dsclb

   pure subroutine dscld(n, m, x, ldx, tt, ldtt)
   !! Select scaling values for `delta` according to the algorithm given in the ODRPACK95
   !! reference guide.
      ! Routines Called (NONE)
      ! Date Written   860529   (YYMMDD)
      ! Revision Date  920304   (YYMMDD)

      use odrpack_kinds, only: zero, one, ten

      integer, intent(in) :: n
         !! The number of observations.
      integer, intent(in) :: m
         !! The number of columns of data in the independent variable.
      real(wp), intent(in) :: x(ldx, m)
         !! The independent variable.
      integer, intent(in) :: ldx
         !! The leading dimension of array `x`.
      real(wp), intent(out) :: tt(ldtt, m)
         !! The scaling values for `delta`.
      integer, intent(in) :: ldtt
         !! The leading dimension of array `tt`.

      ! Local scalars
      real(wp) :: xmax, xmin
      integer :: i, j
      logical :: bigdif

      ! Variable Definitions (alphabetically)
      !  BIGDIF:  The variable designating whether there is a significant
      !           difference in the magnitudes of the nonzero elements of
      !           X (BIGDIF=.TRUE.) or not (BIGDIF=.FALSE.).
      !  I:       An indexing variable.
      !  J:       An indexing variable.
      !  LDTT:    The leading dimension of array TT.
      !  LDX:     The leading dimension of array X.
      !  M:       The number of columns of data in the independent variable.
      !  N:       The number of observations.
      !  TT:      THE SCALING VALUES FOR DELTA.
      !  X:       The independent variable.
      !  XMAX:    The largest nonzero magnitude.
      !  XMIN:    THE SMALLEST NONZERO MAGNITUDE.

      do j = 1, m
         xmax = abs(x(1, j))
         do i = 2, n
            xmax = max(xmax, abs(x(i, j)))
         end do
         if (xmax == zero) then
            !  All input values of X(I,J), I=1,...,N, are zero
            do i = 1, n
               tt(i, j) = one
            end do
         else
            !  Some of the input values are nonzero
            xmin = xmax
            do i = 1, n
               if (x(i, j) /= zero) then
                  xmin = min(xmin, abs(x(i, j)))
               end if
            end do
            bigdif = log10(xmax) - log10(xmin) >= one
            do i = 1, n
               if (x(i, j) /= zero) then
                  if (bigdif) then
                     tt(i, j) = one/abs(x(i, j))
                  else
                     tt(i, j) = one/xmax
                  end if
               else
                  tt(i, j) = ten/xmin
               end if
            end do
         end if
      end do

   end subroutine dscld

   pure subroutine dsetn(n, m, x, ldx, nrow)
   !! Select the row at which the derivative will be checked.
      ! Routines Called  (None)
      ! Date Written   860529   (YYMMDD)
      ! Revision Date  920304   (YYMMDD)

      use odrpack_kinds, only: zero

      integer, intent(in) :: n
         !! The number of observations.
      integer, intent(in) :: m
         !! The number of columns of data in the independent variable.
      real(wp), intent(in) :: x(ldx, m)
         !! The independent variable.
      integer, intent(in) :: ldx
         !! The leading dimension of array `x`.
      integer, intent(inout) :: nrow
         !! The selected row number of the independent variable.

      ! Local scalars
      integer :: i

      ! Variable Definitions (alphabetically)
      !  I:       An index variable.
      !  J:       An index variable.
      !  LDX:     The leading dimension of array X.
      !  M:       The number of columns of data in the independent variable.
      !  N:       The number of observations.
      !  NROW:    The selected row number of the independent variable.
      !  X:       The independent variable.

      if ((nrow >= 1) .and. (nrow <= n)) return

      ! Select first row of independent variables which contains no zeros
      ! if there is one, otherwise first row is used.
      nrow = 1
      do i = 1, n
         if (all(x(i, :) /= zero)) then
            nrow = i
            return
         end if
      end do

   end subroutine dsetn

   subroutine dsolve(n, t, ldt, b, job)
   !! Solve systems of the form:
   !!
   !!  `t * x = b  or  trans(t) * x = b`
   !!
   !! where `t` is an upper or lower triangular matrix of order `n`, and the solution `x`
   !! overwrites the RHS `b`.
   !! Adapted from LINPACK subroutine DTRSL.
   !!
   !! References:
   !! * Dongarra J.J., Bunch J.R., Moler C.B., Stewart G.W., *LINPACK Users Guide*, SIAM, 1979.
      ! @note: This is one of the most time-consuming subroutines in ODRPACK (~25% of total).
      ! Routines Called  DAXPY, DDOT
      ! Date Written   920220   (YYMMDD)
      ! Revision Date  920619   (YYMMDD)

      use odrpack_kinds, only: zero

      integer, intent(in) :: n
         !! The number of rows and columns of data in array `t`.
      real(wp), intent(in) :: t(ldt, n)
         !! The upper or lower tridiagonal system.
      integer, intent(in) :: ldt
         !! The leading dimension of array `t`.
      real(wp), intent(inout) :: b(n)
         !! On input: the right hand side; On exit: the solution.
      integer, intent(in) :: job
         !! What kind of system is to be solved:
         !!   1: Solve `t * x = b`, where `t` is lower triangular,
         !!   2: Solve `t * x = b`, where `t` is upper triangular,
         !!   3: Solve `trans(t) * x = b`, where `t` is lower triangular,
         !!   4: Solve `trans(t) * x = b`, where `t` is upper triangular.

      ! Local scalars
      real(wp) :: temp
      integer :: j1, j, jn

      ! External BLAS/LINPACK procedures
      real(wp), external :: ddot
      external :: daxpy

      ! Variable Definitions (alphabetically)
      !  B:       On input:  the right hand side;  On exit:  the solution
      !  J1:      The first nonzero entry in T.
      !  J:       An indexing variable.
      !  JN:      The last nonzero entry in T.
      !  JOB:     What kind of system is to be solved, where if JOB is
      !           1   Solve T*X=B, T lower triangular,
      !           2   Solve T*X=B, T upper triangular,
      !           3   Solve trans(T)*X=B, T lower triangular,
      !           4   Solve trans(T)*X=B, T upper triangular.
      !  LDT:     The leading dimension of array T.
      !  N:       The number of rows and columns of data in array T.
      !  T:       The upper or lower tridiagonal system.

      !  Find first nonzero diagonal entry in T
      j1 = 0
      do j = 1, n
         if (j1 == 0 .and. t(j, j) /= zero) then
            j1 = j
         elseif (t(j, j) == zero) then
            b(j) = zero
         end if
      end do
      if (j1 == 0) return

      ! Find last nonzero diagonal entry in T
      jn = 0
      do j = n, j1, -1
         if (jn == 0 .and. t(j, j) /= zero) then
            jn = j
         elseif (t(j, j) == zero) then
            b(j) = zero
         end if
      end do

      if (job == 1) then
         ! Solve T*X=B for T lower triangular
         b(j1) = b(j1)/t(j1, j1)
         do j = j1 + 1, jn
            temp = -b(j - 1)
            call daxpy(jn - j + 1, temp, t(j, j - 1), 1, b(j), 1)
            if (t(j, j) /= zero) then
               b(j) = b(j)/t(j, j)
            else
               b(j) = zero
            end if
         end do

      elseif (job == 2) then
         ! Solve T*X=B for T upper triangular.
         b(jn) = b(jn)/t(jn, jn)
         do j = jn - 1, j1, -1
            temp = -b(j + 1)
            call daxpy(j, temp, t(1, j + 1), 1, b(1), 1)
            if (t(j, j) /= zero) then
               b(j) = b(j)/t(j, j)
            else
               b(j) = zero
            end if
         end do

      elseif (job == 3) then
         ! Solve trans(T)*X=B for T lower triangular.
         b(jn) = b(jn)/t(jn, jn)
         do j = jn - 1, j1, -1
            b(j) = b(j) - ddot(jn - j + 1, t(j + 1, j), 1, b(j + 1), 1)
            if (t(j, j) /= zero) then
               b(j) = b(j)/t(j, j)
            else
               b(j) = zero
            end if
         end do

      elseif (job == 4) then
         ! Solve trans(T)*X=B for T upper triangular
         b(j1) = b(j1)/t(j1, j1)
         do j = j1 + 1, jn
            b(j) = b(j) - ddot(j - 1, t(1, j), 1, b(1), 1)
            if (t(j, j) /= zero) then
               b(j) = b(j)/t(j, j)
            else
               b(j) = zero
            end if
         end do
      end if

   end subroutine dsolve

   subroutine dunpac(n2, v1, v2, ifix)
   !! Copy the elements of `v1` into the locations of `v2` which are unfixed.
      ! Routines Called  DCOPY
      ! Date Written   860529   (YYMMDD)
      ! Revision Date  920304   (YYMMDD)

      integer, intent(in) :: n2
         !! The number of items in `v2`.
      real(wp), intent(in) :: v1(n2)
         !! The vector of the unfixed items.
      real(wp), intent(out) :: v2(n2)
         !! The vector of the fixed and unfixed items into which the elements of `v1` are to
         !! be inserted.
      integer, intent(in) :: ifix(n2)
         !! The values designating whether the elements of `v2` are fixed at their input values or not.

      ! Local scalars
      integer :: i, n1

      ! External BLAS/LINPACK procedures
      external :: dcopy

      ! Variable Definitions (alphabetically)
      !  I:       An indexing variable.
      !  IFIX:    The values designating whether the elements of V2 are fixed at their input
      !           values or not.
      !  N1:      The number of items in V1.
      !  N2:      The number of items in V2.
      !  V1:      The vector of the unfixed items.
      !  V2:      The vector of the fixed and unfixed items into which the elements of V1 are
      !           to be inserted.

      n1 = 0
      if (ifix(1) >= 0) then
         do i = 1, n2
            if (ifix(i) /= 0) then
               n1 = n1 + 1
               v2(i) = v1(n1)
            end if
         end do
      else
         n1 = n2
         call dcopy(n2, v1, 1, v2, 1)
      end if

   end subroutine dunpac

   subroutine dvevtr &
      (m, nq, indx, &
       v, ldv, ld2v, e, lde, ve, ldve, ld2ve, vev, ldvev, &
       wrk5)
   !! Compute `v*e*trans(v)` for the (`indx`)th `m` by `nq` array in `v`.
      ! Routines Called  DSOLVE
      ! Date Written   910613   (YYMMDD)
      ! Revision Date  920304   (YYMMDD)

      use odrpack_kinds, only: zero

      integer, intent(in) :: m
         !! The number of columns of data in the independent variable.
      integer, intent(in) :: nq
         !! The number of responses per observation.
      integer, intent(in) :: indx
         !! The row in `v` in which the `m` by `nq` array is stored.
      integer, intent(in) :: ldv
         !! The leading dimension of array `v`.
      integer, intent(in) :: ld2v
         !! The second dimension of array `v`.
      integer, intent(in) :: lde
         !! The leading dimension of array `e`.
      integer, intent(in) :: ldve
         !! The leading dimension of array `ve`.
      integer, intent(in) :: ldvev
         !! The leading dimension of array `vev`.
      integer, intent(in) :: ld2ve
         !! The second dimension of array `ve`.
      real(wp), intent(in) :: v(ldv, ld2v, nq)
         !! An array of `nq` by `m` matrices.
      real(wp), intent(in) :: e(lde, m)
         !! The `m` by `m` matrix of the factors, so `ete = (d**2 + alpha*t**2)`.
      real(wp), intent(out) :: ve(ldve, ld2ve, m)
         !! The `nq` by `m` array `ve = v * inv(e)`.
      real(wp), intent(out) :: vev(ldvev, nq)
         !! The `nq` by `nq` array `vev = v * inv(ete) * trans(v)`.
      real(wp), intent(out) :: wrk5(m)
         !! An `m` work vector.

      ! Local scalars
      integer :: j, l1, l2

      ! Variable Definitions (alphabetically)
      !  INDX:    The row in V in which the M by NQ array is stored.
      !  J:       An indexing variable.
      !  LDE:     The leading dimension of array E.
      !  LDV:     The leading dimension of array V.
      !  LDVE:    The leading dimension of array VE.
      !  LDVEV:   The leading dimension of array VEV.
      !  LD2V:    The second dimension of array V.
      !  L1:      An indexing variable.
      !  L2:      An indexing variable.
      !  M:       The number of columns of data in the independent variable.
      !  NQ:      The number of responses per observation.
      !  E:       The M by M matrix of the factors so ETE = (D**2 + ALPHA*T**2).
      !  V:       An array of NQ by M matrices.
      !  VE:      The NQ by M array VE = V * inv(E)
      !  VEV:     The NQ by NQ array VEV = V * inv(ETE) * trans(V).
      !  WRK5:    An M work vector.

      if (nq == 0 .or. m == 0) return

      do l1 = 1, nq
         do j = 1, m
            wrk5(j) = v(indx, j, l1)
         end do
         call dsolve(m, e, lde, wrk5, 4)
         do j = 1, m
            ve(indx, l1, j) = wrk5(j)
         end do
      end do

      do l1 = 1, nq
         do l2 = 1, l1
            vev(l1, l2) = zero
            do j = 1, m
               vev(l1, l2) = vev(l1, l2) + ve(indx, l1, j)*ve(indx, l2, j)
            end do
            vev(l2, l1) = vev(l1, l2)
         end do
      end do

   end subroutine dvevtr

   pure subroutine dwght(n, m, wt, ldwt, ld2wt, t, wtt)
   !! Scale matrix `t` using `wt`, i.e., compute `wtt = wt*t`.
      ! Routines Called  (NONE)
      ! Date Written   860529   (YYMMDD)
      ! Revision Date  920304   (YYMMDD)

      use odrpack_kinds, only: zero

      integer, intent(in) :: n
         !! The number of rows of data in `t`.
      integer, intent(in) :: m
         !! The number of columns of data in `t`.
      integer, intent(in) :: ldwt
         !! The leading dimension of array `wt`.
      integer, intent(in) :: ld2wt
         !! The second dimension of array `wt`.
      real(wp), intent(in) :: wt(:, :, :)
         !! The weights.
      real(wp), intent(in) :: t(:, :)
         !! The array being scaled by `wt`.
      real(wp), intent(out) :: wtt(:, :)
         !! The results of weighting array `t` by `wt`. Array `wtt` can be the same as `t` only if
         !! the arrays in `wt` are upper triangular with zeros below the diagonal.

      ! Local scalars
      real(wp) :: temp
      integer :: i, j, k

      ! Variable Definitions (alphabetically)
      !  I:       An indexing variable.
      !  J:       An indexing variable.
      !  K:       An indexing variable.
      !  LDWT:    The leading dimension of array WT.
      !  LD2WT:   The second dimension of array WT.
      !  M:       The number of columns of data in T.
      !  N:       The number of rows of data in T.
      !  T:       The array being scaled by WT.
      !  TEMP:    A temporary scalar.
      !  WT:      The weights.
      !  WTT:     The results of weighting array T by WT. Array WTT can be the same as T only if
      !           the arrays in WT are upper triangular with zeros below the diagonal.

      if (n == 0 .or. m == 0) return

      if (wt(1, 1, 1) >= zero) then
         if (ldwt >= n) then
            if (ld2wt >= m) then
               ! WT is an N-array of M by M matrices
               do i = 1, n
                  do j = 1, m
                     temp = zero
                     do k = 1, m
                        temp = temp + wt(i, j, k)*t(i, k)
                     end do
                     wtt(i, j) = temp
                  end do
               end do
            else
               ! WT is an N-array of diagonal matrices
               do i = 1, n
                  do j = 1, m
                     wtt(i, j) = wt(i, 1, j)*t(i, j)
                  end do
               end do
            end if
         else
            if (ld2wt >= m) then
               ! WT is an M by M matrix
               do i = 1, n
                  do j = 1, m
                     temp = zero
                     do k = 1, m
                        temp = temp + wt(1, j, k)*t(i, k)
                     end do
                     wtt(i, j) = temp
                  end do
               end do
            else
               ! WT is a diagonal matrice
               do i = 1, n
                  do j = 1, m
                     wtt(i, j) = wt(1, 1, j)*t(i, j)
                  end do
               end do
            end if
         end if
      else
         ! WT is a scalar
         do j = 1, m
            do i = 1, n
               wtt(i, j) = abs(wt(1, 1, 1))*t(i, j)
            end do
         end do
      end if

   end subroutine dwght

   pure subroutine dwinf &
      (n, m, np, nq, ldwe, ld2we, isodr, &
       deltai, epsi, xplusi, fni, sdi, vcvi, &
       rvari, wssi, wssdei, wssepi, rcondi, etai, &
       olmavi, taui, alphai, actrsi, pnormi, rnorsi, prersi, &
       partli, sstoli, taufci, epsmai, &
       beta0i, betaci, betasi, betani, si, ssi, ssfi, qrauxi, ui, &
       fsi, fjacbi, we1i, diffi, &
       deltsi, deltni, ti, tti, omegai, fjacdi, &
       wrk1i, wrk2i, wrk3i, wrk4i, wrk5i, wrk6i, wrk7i, &
       loweri, upperi, &
       lwkmn)
   !! Set storage locations within REAL (wp) work space.
      ! Routines Called  (NONE)
      ! Date Written   860529   (YYMMDD)
      ! Revision Date  920619   (YYMMDD)
      integer, intent(in) :: n
         !! The number of observations.
      integer, intent(in) :: m
         !! The number of columns of data in the explanatory variable.
      integer, intent(in) :: np
         !! The number of function parameters.
      integer, intent(in) :: nq
         !! The number of responses per observation.
      integer, intent(in) :: ldwe
         !! The leading dimension of array `we`.
      integer, intent(in) :: ld2we
         !! The second dimension of array `we`.
      logical, intent(in) :: isodr
         !! The variable designating whether the solution is by ODR (`isodr`=.true.) or by OLS (`isodr`=.false.).
      integer, intent(out) :: deltai
         !! The starting location in array `work` of array `delta`.
      integer, intent(out) :: epsi
         !! The starting location in array `work` of array `eps`.
      integer, intent(out) :: xplusi
         !! The starting location in array `work` of array `xplusd`.
      integer, intent(out) :: fni
         !! The starting location in array `work` of array `fn`.
      integer, intent(out) :: sdi
         !! The starting location in array `work` of array `sd`.
      integer, intent(out) :: vcvi
         !! The starting location in array `work` of array `vcv`.
      integer, intent(out) :: rvari
         !! The location in array `work` of variable `rvar`.
      integer, intent(out) :: wssi
         !! The location in array `work` of variable `wss`.
      integer, intent(out) :: wssdei
         !! The location in array `work` of variable `wssdel`.
      integer, intent(out) :: wssepi
         !! The location in array `work` of variable `wsseps`.
      integer, intent(out) :: rcondi
         !! The location in array `work` of variable `rcondi`.
      integer, intent(out) :: etai
         !! The location in array `work` of variable `eta`.
      integer, intent(out) :: olmavi
         !! The location in array `work` of variable `olmavg`.
      integer, intent(out) :: taui
         !! The location in array `work` of variable `tau`.
      integer, intent(out) :: alphai
         !! The location in array `work` of variable `alpha`.
      integer, intent(out) :: actrsi
         !! The location in array `work` of variable `actrs`.
      integer, intent(out) :: pnormi
         !! The location in array `work` of variable `pnorm`.
      integer, intent(out) :: rnorsi
         !! The location in array `work` of variable `rnorms`.
      integer, intent(out) :: prersi
         !! The location in array `work` of variable `prers`.
      integer, intent(out) :: partli
         !! The location in array `work` of variable `partol`.
      integer, intent(out) :: sstoli
         !! The location in array `work` of variable `sstol`.
      integer, intent(out) :: taufci
         !! The location in array `work` of variable `taufac`.
      integer, intent(out) :: epsmai
         !! The location in array `work` of variable `epsmac`.
      integer, intent(out) :: beta0i
         !! The starting location in array `work` of array `beta0`.
      integer, intent(out) :: betaci
         !! The starting location in array `work` of array `betac`.
      integer, intent(out) :: betasi
         !! The starting location in array `work` of array `betas`.
      integer, intent(out) :: betani
         !! The starting location in array `work` of array `betan`.
      integer, intent(out) :: si
         !! The starting location in array `work` of array `s`.
      integer, intent(out) :: ssi
         !! The starting location in array `work` of array `ss`.
      integer, intent(out) :: ssfi
         !! The starting location in array `work` of array `ssf`.
      integer, intent(out) :: qrauxi
         !! The starting location in array `work` of array `qraux`.
      integer, intent(out) :: ui
         !! The starting location in array `work` of array `u`.
      integer, intent(out) :: fsi
         !! The starting location in array `work` of array `fs`.
      integer, intent(out) :: fjacbi
         !! The starting location in array `work` of array `fjacb`.
      integer, intent(out) :: we1i
         !! The starting location in array `work` of array `we1`.
      integer, intent(out) :: diffi
         !! The starting location in array `work` of array `diff`.
      integer, intent(out) :: deltsi
         !! The starting location in array `work` of array `deltas`.
      integer, intent(out) :: deltni
         !! The starting location in array `work` of array `deltan`.
      integer, intent(out) :: ti
         !! The starting location in array `work` of array `t`.
      integer, intent(out) :: tti
         !! The starting location in array `work` of array `tt`.
      integer, intent(out) :: omegai
         !! The starting location in array `work` of array `omega`.
      integer, intent(out) :: fjacdi
         !! The starting location in array `work` of array `fjacd`.
      integer, intent(out) :: wrk1i
         !! The starting location in array `work` of array `wrk1`.
      integer, intent(out) :: wrk2i
         !! The starting location in array `work` of array `wrk2`.
      integer, intent(out) :: wrk3i
         !! The starting location in array `work` of array `wrk3`.
      integer, intent(out) :: wrk4i
         !! The starting location in array `work` of array `wrk4`.
      integer, intent(out) :: wrk5i
         !! The starting location in array `work` of array `wrk5`.
      integer, intent(out) :: wrk6i
         !! The starting location in array `work` of array `wrk6`.
      integer, intent(out) :: wrk7i
         !! The starting location in array `work` of array `wrk7`.
      integer, intent(out) :: loweri
         !! The starting location in array `work` of array `lower`.
      integer, intent(out) :: upperi
         !! The starting location in array `work` of array `upper`.
      integer, intent(out) :: lwkmn
         !! The minimum acceptable length of vector `work`.

      ! Local scalars
      integer :: next

      ! Variable Definitions (alphabetically)
      !  ACTRSI:  The location in array WORK of variable ACTRS.
      !  ALPHAI:  The location in array WORK of variable ALPHA.
      !  BETACI:  The starting location in array WORK of array BETAC.
      !  BETANI:  The starting location in array WORK of array BETAN.
      !  BETASI:  The starting location in array WORK of array BETAS.
      !  BETA0I:  The starting location in array WORK of array BETA0.
      !  DELTAI:  The starting location in array WORK of array DELTA.
      !  DELTNI:  The starting location in array WORK of array DELTAN.
      !  DELTSI:  The starting location in array WORK of array DELTAS.
      !  DIFFI:   The starting location in array WORK of array DIFF.
      !  EPSI:    The starting location in array WORK of array EPS.
      !  EPSMAI:  The location in array WORK of variable EPSMAC.
      !  ETAI:    The location in array WORK of variable ETA.
      !  FJACBI:  The starting location in array WORK of array FJACB.
      !  FJACDI:  The starting location in array WORK of array FJACD.
      !  FNI:     The starting location in array WORK of array FN.
      !  FSI:     The starting location in array WORK of array FS.
      !  ISODR:   The variable designating whether the solution is by ODR (ISODR=TRUE) or
      !           by OLS (ISODR=FALSE).
      !  LDWE:    The leading dimension of array WE.
      !  LD2WE:   The second dimension of array WE.
      !  LWKMN:   The minimum acceptable length of vector work.
      !  M:       The number of columns of data in the explanatory variable.
      !  N:       The number of observations.
      !  NEXT:    The next available location with WORK.
      !  NP:      The number of function parameters.
      !  NQ:      The number of responses per observation.
      !  OLMAVI:  The location in array WORK of variable OLMAVG.
      !  OMEGAI:  The starting location in array WORK of array OMEGA.
      !  PARTLI:  The location in array WORK of variable PARTOL.
      !  PNORMI:  The location in array WORK of variable PNORM.
      !  PRERSI:  The location in array WORK of variable PRERS.
      !  QRAUXI:  The starting location in array WORK of array QRAUX.
      !  RCONDI:  The location in array WORK of variable RCONDI.
      !  RNORSI:  The location in array WORK of variable RNORMS.
      !  RVARI:   The location in array WORK of variable RVAR.
      !  SDI:     The starting location in array WORK of array SD.
      !  SI:      The starting location in array WORK of array S.
      !  SSFI:    The starting location in array WORK of array SSF.
      !  SSI:     The starting location in array WORK of array SS.
      !  SSTOLI:  The location in array WORK of variable SSTOL.
      !  TAUFCI:  The location in array WORK of variable TAUFAC.
      !  TAUI:    The location in array WORK of variable TAU.
      !  TI:      The starting location in array WORK of array T.
      !  TTI:     The starting location in array WORK of array TT.
      !  UI:      The starting location in array WORK of array U.
      !  VCVI:    The starting location in array WORK of array VCV.
      !  WE1I:    The starting location in array WORK of array WE1.
      !  WRK1I:   The starting location in array WORK of array WRK1.
      !  WRK2I:   The starting location in array WORK of array WRK2.
      !  WRK3I:   The starting location in array WORK of array WRK3.
      !  WRK4I:   The starting location in array WORK of array WRK4.
      !  WRK5I:   The starting location in array WORK of array WRK5.
      !  WRK6I:   The starting location in array WORK of array WRK6.
      !  WRK7I:   The starting location in array WORK of array WRK7.
      !  WSSI:    The location in array WORK of variable WSS.
      !  WSSDEI:  The location in array WORK of variable WSSDEL.
      !  WSSEPI:  The location in array work of variable WSSEPS.
      !  XPLUSI:  The starting location in array WORK of array XPLUSD.

      if (n >= 1 .and. m >= 1 .and. np >= 1 .and. nq >= 1 .and. &
          ldwe >= 1 .and. ld2we >= 1) then

         deltai = 1
         epsi = deltai + n*m
         xplusi = epsi + n*nq
         fni = xplusi + n*m
         sdi = fni + n*nq
         vcvi = sdi + np
         rvari = vcvi + np*np

         wssi = rvari + 1
         wssdei = wssi + 1
         wssepi = wssdei + 1
         rcondi = wssepi + 1
         etai = rcondi + 1
         olmavi = etai + 1

         taui = olmavi + 1
         alphai = taui + 1
         actrsi = alphai + 1
         pnormi = actrsi + 1
         rnorsi = pnormi + 1
         prersi = rnorsi + 1
         partli = prersi + 1
         sstoli = partli + 1
         taufci = sstoli + 1
         epsmai = taufci + 1
         beta0i = epsmai + 1

         betaci = beta0i + np
         betasi = betaci + np
         betani = betasi + np
         si = betani + np
         ssi = si + np
         ssfi = ssi + np
         qrauxi = ssfi + np
         ui = qrauxi + np
         fsi = ui + np

         fjacbi = fsi + n*nq

         we1i = fjacbi + n*np*nq

         diffi = we1i + ldwe*ld2we*nq

         next = diffi + nq*(np + m)

         if (isodr) then
            deltsi = next
            deltni = deltsi + n*m
            ti = deltni + n*m
            tti = ti + n*m
            omegai = tti + n*m
            fjacdi = omegai + nq*nq
            wrk1i = fjacdi + n*m*nq
            next = wrk1i + n*m*nq
         else
            deltsi = deltai
            deltni = deltai
            ti = deltai
            tti = deltai
            omegai = deltai
            fjacdi = deltai
            wrk1i = deltai
         end if

         wrk2i = next
         wrk3i = wrk2i + n*nq
         wrk4i = wrk3i + np
         wrk5i = wrk4i + m*m
         wrk6i = wrk5i + m
         wrk7i = wrk6i + n*nq*np
         loweri = wrk7i + 5*nq
         upperi = loweri + np
         next = upperi + np

         lwkmn = next
      else
         deltai = 1
         epsi = 1
         xplusi = 1
         fni = 1
         sdi = 1
         vcvi = 1
         rvari = 1
         wssi = 1
         wssdei = 1
         wssepi = 1
         rcondi = 1
         etai = 1
         olmavi = 1
         taui = 1
         alphai = 1
         actrsi = 1
         pnormi = 1
         rnorsi = 1
         prersi = 1
         partli = 1
         sstoli = 1
         taufci = 1
         epsmai = 1
         beta0i = 1
         betaci = 1
         betasi = 1
         betani = 1
         si = 1
         ssi = 1
         ssfi = 1
         qrauxi = 1
         fsi = 1
         ui = 1
         fjacbi = 1
         we1i = 1
         diffi = 1
         deltsi = 1
         deltni = 1
         ti = 1
         tti = 1
         fjacdi = 1
         omegai = 1
         wrk1i = 1
         wrk2i = 1
         wrk3i = 1
         wrk4i = 1
         wrk5i = 1
         wrk6i = 1
         wrk7i = 1
         loweri = 1
         upperi = 1
         lwkmn = 1
      end if

   end subroutine dwinf

   pure subroutine mbfb(np, beta, lower, upper, ssf, stpb, neta, eta, interval)
   !! Ensure range of bounds is large enough for derivative checking.
   !! Move beta away from bounds so that derivatives can be calculated.
      ! ROUTINES CALLED  DHSTEP
      ! DATE WRITTEN   20040624   (YYYYMMDD)
      ! REVISION DATE  20040624   (YYYYMMDD)

      use odrpack_kinds, only: zero, one, three, ten, hundred

      integer, intent(in) :: np
         !! The number of parameters `np`.
      real(wp), intent(inout) :: beta(np)
         !! Function parameters.
      real(wp), intent(in) :: lower(np)
         !! !! Lower bound on `beta`.
      real(wp), intent(in) :: upper(np)
         !! Upper bound on `beta`.
      real(wp), intent(in) :: ssf(np)
         !! The scale used for the `beta`s.
      real(wp), intent(in) :: stpb(np)
         !! The relative step used for computing finite difference derivatives with respect to `beta`.
      integer, intent(in) :: neta
         !! Number of good digits in the function results.
      real(wp), intent(in) :: eta
         !! The relative noise in the function results.
      integer, intent(out) :: interval(np)
         !! Specifies which difference methods and step sizes are supported by the current
         !! interval `upper-lower`.

      ! Local scalars
      integer :: k
      real(wp) :: h, h0, h1, hc, hc0, hc1, stpr, stpl, typj

      ! VARIABLE DEFINITIONS (ALPHABETICALLY)
      !  BETA:     BETA for the jacobian checker.  BETA will be moved far enough from the bounds so
      !            that the derivative checker may proceed.
      !  H:        Relative step size for forward differences.
      !  H0:       Initial relative step size for forward differences.
      !  H1:       Default relative step size for forward differences.
      !  HC:       Relative step size for center differences.
      !  HC0:      Initial relative step size for center differences.
      !  HC1:      Default relative step size for center differences.
      !  INTERVAL: Specifies which difference methods and step sizes are supported by the current
      !            interval UPPER-LOWER.
      !  K:        Index variable for BETA.
      !  NETA:     Number of good digits in the function results.
      !  SSF:      The scale used for the BETA'S.
      !  STPB:     The relative step used for computing finite difference derivatives with respect
      !            to BETA.
      !  STPL:     Maximum step to the left of BETA (-) the derivative checker will use.
      !  STPR:     Maximum step to the right of BETA (+) the derivative checker will use.
      !  TYPJ:     The typical size of the J-th unkonwn BETA.

      interval(:) = 111
      do k = 1, np
         h0 = dhstep(0, neta, 1, k, stpb, 1)
         hc0 = h0
         h1 = sqrt(eta)
         hc1 = eta**(one/three)
         h = max(ten*h1, min(hundred*h0, one))
         hc = max(ten*hc1, min(hundred*hc0, one))
         if (beta(k) == zero) then
            if (ssf(1) < zero) then
               typj = one/abs(ssf(1))
            else
               typj = one/ssf(k)
            end if
         else
            typj = abs(beta(k))
         end if
         stpr = (h*typj*sign(one, beta(k)) + beta(k)) - beta(k)
         stpl = (hc*typj*sign(one, beta(k)) + beta(k)) - beta(k)
         ! Check outer interval
         if (lower(k) + 2*abs(stpl) > upper(k)) then
            if (interval(k) >= 100) then
               interval(k) = interval(k) - 100
            end if
         elseif (beta(k) + stpl > upper(k) .or. beta(k) - stpl > upper(k)) then
            beta(k) = upper(k) - abs(stpl)
         elseif (beta(k) + stpl < lower(k) .or. beta(k) - stpl < lower(k)) then
            beta(k) = lower(k) + abs(stpl)
         end if
         ! Check middle interval
         if (lower(k) + 2*abs(stpr) > upper(k)) then
            if (mod(interval(k), 100) >= 10) then
               interval(k) = interval(k) - 10
            end if
         elseif (beta(k) + stpr > upper(k) .or. beta(k) - stpr > upper(k)) then
            beta(k) = upper(k) - abs(stpr)
         elseif (beta(k) + stpr < lower(k) .or. beta(k) - stpr < lower(k)) then
            beta(k) = lower(k) + abs(stpr)
         end if
         ! Check inner interval
         if (lower(k) + abs(stpr) > upper(k)) then
            interval(k) = 0
         elseif (beta(k) + stpr > upper(k)) then
            beta(k) = upper(k) - stpr
         elseif (beta(k) + stpr < lower(k)) then
            beta(k) = lower(k) - stpr
         end if
      end do

   end subroutine mbfb

end module odrpack_core
