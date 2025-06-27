module odrpack_core
!! Core mathematical routines, except drivers, and BLAS/LINPACK.

   use odrpack_kinds, only: wp
   implicit none

   abstract interface
      subroutine fcn_t(beta, xplusd, ifixb, ifixx, ideval, f, fjacb, fjacd, istop)
      !! User-supplied subroutine for evaluating the model.
         import :: wp
         implicit none
         real(wp), intent(in) :: beta(:)
            !! Current values of parameters.
         real(wp), intent(in) :: xplusd(:, :)
            !! Current value of explanatory variable, i.e., `x + delta`. Shape is `(n, m)`.
         integer, intent(in) :: ifixb(:)
            !! Indicators for "fixing" parameters (`beta`).
         integer, intent(in) :: ifixx(:, :)
            !! Indicators for "fixing" explanatory variable (`x`). Shape is `(ldifx, m)`.
         integer, intent(in) :: ideval
            !! Indicator for selecting computation to be performed.
         real(wp), intent(out) :: f(:, :)
            !! Predicted function values. Shape is `(n, q)`.
         real(wp), intent(out) :: fjacb(:, :, :)
            !! Jacobian with respect to `beta`. Shape is `(n, np, q)`.
         real(wp), intent(out) :: fjacd(:, :, :)
            !! Jacobian with respect to errors `delta`. Shape is `(n, m, q)`.
         integer, intent(out) :: istop
            !! Stopping condition, with meaning as follows. 0 means current `beta` and
            !! `x + delta` were acceptable and values were computed successfully. 1 means current
            !! `beta` and `x + delta` are not acceptable;  'odrpack' should select values closer
            !! to most recently used values if possible. -1 means current `beta` and `x + delta`
            !! are not acceptable; 'odrpack' should stop.
      end subroutine fcn_t
   end interface

contains

   subroutine trust_region &
      (n, m, np, q, npp, &
       f, fjacb, fjacd, &
       wd, ldwd, ld2wd, ss, tt, ldtt, delta, &
       alpha2, tau, epsfcn, isodr, &
       tfjacb, omega, u, qraux, jpvt, &
       s, t, nlms, rcond, irank, &
       wrk1, wrk2, wrk3, wrk4, wrk5, wrk, lwrk, istopc)
   !! Compute Levenberg-Marquardt parameter and steps `s` and `t` using analog of the
   !! trust-region Levenberg-Marquardt algorithm.

      use odrpack_kinds, only: zero
      use blas_interfaces, only: ddot, dnrm2

      integer, intent(in) :: n
         !! Number of observations.
      integer, intent(in) :: m
         !! Number of columns of data in the explanatory variable.
      integer, intent(in) :: np
         !! Number of function parameters.
      integer, intent(in) :: q
         !! Number of responses per observation.
      integer, intent(in) :: npp
         !! Number of function parameters being estimated.
      real(wp), intent(in) :: f(n, q)
         !! Weighted estimated values of `epsilon`.
      real(wp), intent(in) :: fjacb(n, np, q)
         !! Jacobian with respect to `beta`.
      real(wp), intent(in) :: fjacd(n, m, q)
         !! Jacobian with respect to `delta`.
      real(wp), intent(in) :: wd(ldwd, ld2wd, m)
         !! `delta` weights.
      integer, intent(in) :: ldwd
         !! Leading dimension of array `wd`.
      integer, intent(in) :: ld2wd
         !! Second dimension of array `wd`.
      real(wp), intent(in) :: ss(np)
         !! Scaling values used for the unfixed `beta`s.
      real(wp), intent(in) :: tt(ldtt, m)
         !! Scale used for the `delta`s.
      integer, intent(in) :: ldtt
         !! Leading dimension of array `tt`.
      real(wp), intent(in) :: delta(n, m)
         !! Estimated errors in the explanatory variables.
      real(wp), intent(inout) :: alpha2
         !! Current Levenberg-Marquardt parameter.
      real(wp), intent(inout) :: tau
         !! Trust region diameter.
      real(wp), intent(in) :: epsfcn
         !! Function's precision.
      logical, intent(in) :: isodr
         !! Variable designating whether the solution is by ODR (`isodr = .true.`)
         !! or by OLS (`isodr = .false.`).
      real(wp), intent(out) :: tfjacb(n, q, np)
         !! Array `omega*fjacb`.
      real(wp), intent(out) :: omega(q, q)
         !! Array `(I-fjacd*inv(p)*trans(fjacd))**(-1/2)`.
      real(wp), intent(out) :: u(np)
         !! Approximate null vector for `tfjacb`.
      real(wp), intent(out) :: qraux(np)
         !! Array required to recover the orthogonal part of the Q-R decomposition.
      integer, intent(out) :: jpvt(np)
         !! Pivot vector.
      real(wp), intent(out) :: s(np)
         !! Step for `beta`.
      real(wp), intent(out) :: t(n, m)
         !! Step for `delta`.
      integer, intent(out) :: nlms
         !! Number of Levenberg-Marquardt steps taken.
      real(wp), intent(out) :: rcond
         !! Approximate reciprocal condition of `tfjacb`.
      integer, intent(out) :: irank
         !! Aank deficiency of the Jacobian wrt `beta`.
      real(wp), intent(out) :: wrk1(n, q, m)
         !! Work array of `(n, q, m)` elements.
      real(wp), intent(out) :: wrk2(n, q)
         !! Work array of `(n, q)` elements.
      real(wp), intent(out) :: wrk3(np)
         !! Work array of `(np)` elements.
      real(wp), intent(out) :: wrk4(m, m)
         !! Work array of `(m, m)` elements.
      real(wp), intent(out) :: wrk5(m)
         !! Work array of `(m)` elements.
      real(wp), intent(out) :: wrk(lwrk)
         !! Work array of `(lwrk)` elements, _equivalenced_ to `wrk1` and `wrk2`.
      integer, intent(in) :: lwrk
         !! Length of vector `wrk`.
      integer, intent(out) :: istopc
         !! Variable designating whether the computations were stopped due to some other
         !! numerical error detected within subroutine `dodstp`.

      ! Local scalars
      real(wp), parameter :: p001 = 0.001_wp, p1 = 0.1_wp
      real(wp) :: alpha1, alphan, bot, phi1, phi2, sa, top
      integer :: i, iwrk, j, k
      logical :: forvcv

      ! Variable Definitions (alphabetically)
      !  ALPHAN:  The new Levenberg-Marquardt parameter.
      !  ALPHA1:  The previous Levenberg-Marquardt parameter.
      !  BOT:     The lower limit for setting ALPHA.
      !  FORVCV:  The variable designating whether this subroutine was called to set up for the
      !           covariance matrix computations (FORVCV=TRUE) or not (FORVCV=FALSE).
      !  I:       An indexing variable.
      !  IWRK:    An indexing variable.
      !  J:       An indexing variable.
      !  K:       An indexing variable.
      !  PHI1:    The previous difference between the norm of the scaled step and the trust
      !           region diameter.
      !  PHI2:    The current difference between the norm of the scaled step and the trust region
      !           diameter.
      !  SA:      The scalar PHI2*(ALPHA1-ALPHA2)/(PHI1-PHI2).
      !  TOP:     The upper limit for setting ALPHA.

      forvcv = .false.
      istopc = 0

      ! Compute full Gauss-Newton step (ALPHA=0)
      alpha1 = zero
      call lcstep(n, m, np, q, npp, &
                  f, fjacb, fjacd, &
                  wd, ldwd, ld2wd, ss, tt, ldtt, delta, &
                  alpha1, epsfcn, isodr, &
                  tfjacb, omega, u, qraux, jpvt, &
                  s, t, phi1, irank, rcond, forvcv, &
                  wrk1, wrk2, wrk3, wrk4, wrk5, wrk, lwrk, istopc)
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
         tfjacb(1:n, 1:q, k) = fjacb(1:n, k, 1:q)
         wrk(k) = ddot(n*q, tfjacb(1, 1, k), 1, f(1, 1), 1)
      end do

      call scale_vec(npp, 1, ss, npp, wrk, npp, wrk, npp) ! work is input (as t) and output (as sclt)

      if (isodr) then
         call scale_mat(n, m, wd, ldwd, ld2wd, delta, wrk(npp + 1:npp + 1 + n*m - 1))
         iwrk = npp
         do j = 1, m
            do i = 1, n
               iwrk = iwrk + 1
               wrk(iwrk) = wrk(iwrk) + ddot(q, fjacd(i, j, 1), n*m, f(i, 1), n)
            end do
         end do
         call scale_vec(n, m, tt, ldtt, wrk(npp + 1), n, wrk(npp + 1), n)
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
         call lcstep(n, m, np, q, npp, &
                     f, fjacb, fjacd, &
                     wd, ldwd, ld2wd, ss, tt, ldtt, delta, &
                     alpha2, epsfcn, isodr, &
                     tfjacb, omega, u, qraux, jpvt, &
                     s, t, phi2, irank, rcond, forvcv, &
                     wrk1, wrk2, wrk3, wrk4, wrk5, wrk, lwrk, istopc)
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

   end subroutine trust_region

   pure subroutine access_workspace &
      (n, m, np, q, ldwe, ld2we, &
       rwork, lrwork, iwork, liwork, &
       access, isodr, &
       jpvt, omega, u, qraux, sd, vcv, wrk1, wrk2, wrk3, wrk4, wrk5, wrk6, &
       nnzw, npp, &
       job, partol, sstol, maxit, taufac, eta, neta, &
       lunrpt, ipr1, ipr2, ipr2f, ipr3, &
       wss, rvar, idf, &
       tau, alpha, niter, nfev, njev, int2, olmavg, &
       rcond, irank, actrs, pnorm, prers, rnorms, istop)
   !! Access or store values in the work arrays.

      integer, intent(in) :: n
         !! Number of observations.
      integer, intent(in) :: m
         !! Number of columns of data in the explanatory variable.
      integer, intent(in) :: np
         !! Number of function parameters.
      integer, intent(in) :: q
         !! Number of responses per observation.
      integer, intent(in) :: ldwe
         !! Leading dimension of array `we`.
      integer, intent(in) :: ld2we
         !! Second dimension of array `we`.
      real(wp), intent(inout) :: rwork(lrwork)
         !! Real work space.
      integer, intent(in) :: lrwork
         !! Length of vector `rwork`.
      integer, intent(inout) :: iwork(liwork)
         !! Integer work space.
      integer, intent(in) :: liwork
         !! Length of vector `iwork`.
      logical, intent(in) :: access
         !! Variable designating whether information is to be accessed from the work
         !! arrays (`access = .true.`) or stored in them (`access = .false.`).
      logical, intent(in) :: isodr
         !! Variable designating whether the solution is to be found by ODR (`isodr = .true.`)
         !! or by OLS (`isodr = .false.`).
      integer, intent(out) :: jpvt
         !! Pivot vector.
      integer, intent(out) :: omega
         !! Starting location in array `rwork` of array `omega`.
      integer, intent(out) :: u
         !! Starting location in array `rwork` of array `u`.
      integer, intent(out) :: qraux
         !! Starting location in array `rwork` of array `qraux`.
      integer, intent(out) :: sd
         !! Starting location in array `rwork` of array `sd`.
      integer, intent(out) :: vcv
         !! Starting location in array `rwork` of array `vcv`.
      integer, intent(out) :: wrk1
         !! Starting location in array `rwork` of array `wrk1`.
      integer, intent(out) :: wrk2
         !! Starting location in array `rwork` of array `wrk2`.
      integer, intent(out) :: wrk3
         !! Starting location in array `rwork` of array `wrk3`.
      integer, intent(out) :: wrk4
         !! Starting location in array `rwork` of array `wrk4`.
      integer, intent(out) :: wrk5
         !! Starting location in array `rwork` of array `wrk5`.
      integer, intent(out) :: wrk6
         !! Starting location in array `rwork` of array `wrk6`.
      integer, intent(out) :: nnzw
         !! Number of nonzero weighted observations.
      integer, intent(out) :: npp
         !! Number of function parameters actually estimated.
      integer, intent(out) :: job
         !! Variable controlling problem initialization and computational method.
      real(wp), intent(inout) :: partol
         !! Parameter convergence stopping tolerance.
      real(wp), intent(inout) :: sstol
         !! Sum-of-squares convergence stopping tolerance.
      integer, intent(out) :: maxit
         !! Maximum number of iterations allowed.
      real(wp), intent(out) :: taufac
         !! Factor used to compute the initial trust region diameter.
      real(wp), intent(out) :: eta
         !! Relative noise in the function results.
      integer, intent(out) :: neta
         !! Number of accurate digits in the function results.
      integer, intent(out) :: lunrpt
          !! Logical unit number used for computation reports.
      integer, intent(out) :: ipr1
         !! Value of the fourth digit (from the right) of `iprint`, which controls the
         !! initial summary report.
      integer, intent(out) :: ipr2
         !! Value of the third digit (from the right) of `iprint`, which controls the
         !! iteration reports.
      integer, intent(out) :: ipr2f
         !! Value of the second digit (from the right) of `iprint`, which controls the
         !! frequency of the iteration reports.
      integer, intent(out) :: ipr3
         !! Value of the first digit (from the right) of `iprint`, which controls the final
         !! summary report.
      real(wp), intent(inout) :: wss(3)
         !! Sum of the squares of the weighted `epsilons` and `deltas`, the sum of the squares
         !! of the weighted `deltas`, and the sum of the squares of the weighted `epsilons`.
      real(wp), intent(inout) :: rvar
         !! Residual variance, i.e. the standard deviation squared.
      integer, intent(inout) :: idf
         !! Degrees of freedom of the fit, equal to the number of observations with nonzero
         !! weighted derivatives minus the number of parameters being estimated.
      real(wp), intent(inout) :: tau
         !! Trust region diameter.
      real(wp), intent(inout) :: alpha
         !! Levenberg-Marquardt parameter.
      integer, intent(inout) :: niter
         !! Number of iterations taken.
      integer, intent(inout) :: nfev
         !! Number of function evaluations.
      integer, intent(inout) :: njev
         !! Number of Jacobian evaluations.
      integer, intent(inout) :: int2
         !! Number of internal doubling steps.
      real(wp), intent(inout) :: olmavg
         !! Average number of Levenberg-Marquardt steps per iteration.
      real(wp), intent(inout) :: rcond
         !! Approximate reciprocal condition of `fjacb`.
      integer, intent(inout) :: irank
         !! Rank deficiency of the Jacobian wrt `beta`.
      real(wp), intent(inout) :: actrs
         !! Saved actual relative reduction in the sum-of-squares.
      real(wp), intent(inout) :: pnorm
         !! Norm of the scaled estimated parameters.
      real(wp), intent(inout) :: prers
         !! Saved predicted relative reduction in the sum-of-squares.
      real(wp), intent(inout) :: rnorms
         !! Norm of the saved weighted `epsilons` and `deltas`.
      integer, intent(inout) :: istop
         !! Variable designating whether there are problems computing the function at the
         !! current `beta` and `delta`.

      ! Local scalars
      integer :: actrsi, alphai, betaci, betani, betasi, beta0i, boundi, deltai, deltni, &
                 deltsi, diffi, epsi, epsmai, etai, fjacbi, fjacdi, fni, fsi, idfi, int2i, &
                 iprini, iprint, iranki, istopi, jobi, jpvti, ldtti, liwkmn, loweri, luneri, &
                 lunrpi, lrwkmn, maxiti, msgb, msgd, netai, nfevi, niteri, njevi, nnzwi, nppi, &
                 nrowi, ntoli, olmavi, omegai, partli, pnormi, prersi, qrauxi, rcondi, rnorsi, &
                 rvari, sdi, si, ssfi, ssi, sstoli, taufci, taui, ti, tti, ui, upperi, vcvi, &
                 we1i, wrk1i, wrk2i, wrk3i, wrk4i, wrk5i, wrk6i, wrk7i, wssi, wssdei, wssepi, &
                 xplusi

      ! Variable Definitions (alphabetically)
      !  ACTRSI:  The location in array RWORK of variable ACTRS.
      !  ALPHAI:  The location in array RWORK of variable ALPHA.
      !  BETACI:  The starting location in array RWORK of array BETAC.
      !  BETANI:  The starting location in array RWORK of array BETAN.
      !  BETASI:  The starting location in array RWORK of array BETAS.
      !  BETA0I:  The starting location in array RWORK of array BETA0.
      !  DELTAI:  The starting location in array RWORK of array DELTA.
      !  DELTNI:  The starting location in array RWORK of array DELTAN.
      !  DELTSI:  The starting location in array RWORK of array DELTAS.
      !  DIFFI:   The starting location in array RWORK of array DIFF.
      !  EPSI:    The starting location in array RWORK of array EPS.
      !  EPSMAI:  The location in array RWORK of variable EPSMAC.
      !  ETAI:    The location in array RWORK of variable ETA.
      !  FJACBI:  The starting location in array RWORK of array FJACB.
      !  FJACDI:  The starting location in array RWORK of array FJACD.
      !  FNI:     The starting location in array RWORK of array FN.
      !  FSI:     The starting location in array RWORK of array FS.
      !  IDFI:    The starting location in array IWORK of variable IDF.
      !  INT2I:   The location in array IWORK of variable INT2.
      !  IPRINI:  The location in array IWORK of variable IPRINT.
      !  IPRINT:  The print control variable.
      !  IRANKI:  The location in array IWORK of variable IRANK.
      !  ISTOPI:  The location in array IWORK of variable ISTOP.
      !  JOBI:    The location in array IWORK of variable JOB.
      !  JPVTI:   The starting location in array IWORK of variable JPVT.
      !  LDTTI:   The starting location in array IWORK of variable LDTT.
      !  LUNERI:  The location in array IWORK of variable LUNERR.
      !  LUNRPI:  The location in array IWORK of variable LUNRPT.
      !  LRWKMN:  The minimum acceptable length of array RWORK.
      !  MAXITI:  The location in array IWORK of variable MAXIT.
      !  MSGB:    The starting location in array IWORK of array MSGB.
      !  MSGD:    The starting location in array IWORK of array MSGD.
      !  NETAI:   The location in array IWORK of variable NETA.
      !  NFEVI:   The location in array IWORK of variable NFEV.
      !  NITERI:  The location in array IWORK of variable NITER.
      !  NJEVI:   The location in array IWORK of variable NJEV.
      !  NNZWI:   The location in array IWORK of variable NNZW.
      !  NPPI:    The location in array IWORK of variable NPP.
      !  NROWI:   The location in array IWORK of variable NROW.
      !  NTOLI:   The location in array IWORK of variable NTOL.
      !  OLMAVI:  The location in array RWORK of variable OLMAVG.
      !  OMEGAI:  The starting location in array RWORK of array OMEGA.
      !  PARTLI:  The location in array work of variable PARTOL.
      !  PNORMI:  The location in array RWORK of variable PNORM.
      !  PRERSI:  The location in array RWORK of variable PRERS.
      !  QRAUXI:  The starting location in array RWORK of array QRAUX.
      !  RCONDI:  The location in array RWORK of variable RCOND.
      !  RNORSI:  The location in array RWORK of variable RNORMS.
      !  RVARI:   The location in array RWORK of variable RVAR.
      !  SDI:     The starting location in array RWORK of array SD.
      !  SSFI:    The starting location in array RWORK of array SSF.
      !  SSI:     The starting location in array RWORK of array SS.
      !  SSTOLI:  The location in array RWORK of variable SSTOL.
      !  TAUFCI:  The location in array RWORK of variable TAUFAC.
      !  TAUI:    the location in array RWORK of variable TAU.
      !  TI:      The starting location in array RWORK of array T.
      !  TTI:     The starting location in array RWORK of array TT.
      !  UI:      The starting location in array RWORK of array U.
      !  VCVI:    The starting location in array RWORK of array VCV.
      !  WE1I:    The starting location in array RWORK of array WE1.
      !  WRK1I:   The starting location in array RWORK of array WRK1.
      !  WRK2I:   The starting location in array RWORK of array WRK2.
      !  WRK3I:   The starting location in array RWORK of array wrk3.
      !  WRK4I:   The starting location in array RWORK of array wrk4.
      !  WRK5I:   The starting location in array RWORK of array wrk5.
      !  WRK6I:   The starting location in array RWORK of array wrk6.
      !  WRK7I:   The starting location in array RWORK of array wrk7.
      !  WSSI:    The starting location in array RWORK of variable WSS(1).
      !  WSSDEI:  The starting location in array RWORK of variable WSS(2).
      !  WSSEPI:  The starting location in array RWORK of variable WSS(3).
      !  XPLUSI:  The starting location in array RWORK of array XPLUSD.

      ! Find starting locations within integer workspace
      call loc_iwork(m, q, np, &
                     msgb, msgd, jpvti, istopi, &
                     nnzwi, nppi, idfi, &
                     jobi, iprini, luneri, lunrpi, &
                     nrowi, ntoli, netai, &
                     maxiti, niteri, nfevi, njevi, int2i, iranki, ldtti, &
                     boundi, &
                     liwkmn)

      ! Find starting locations within REAL work space
      call loc_rwork(n, m, q, np, ldwe, ld2we, isodr, &
                     deltai, epsi, xplusi, fni, sdi, vcvi, &
                     rvari, wssi, wssdei, wssepi, rcondi, etai, &
                     olmavi, taui, alphai, actrsi, pnormi, rnorsi, prersi, &
                     partli, sstoli, taufci, epsmai, &
                     beta0i, betaci, betasi, betani, si, ssi, ssfi, qrauxi, ui, &
                     fsi, fjacbi, we1i, diffi, &
                     deltsi, deltni, ti, tti, omegai, fjacdi, &
                     wrk1i, wrk2i, wrk3i, wrk4i, wrk5i, wrk6i, wrk7i, &
                     loweri, upperi, &
                     lrwkmn)

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
         actrs = rwork(actrsi)
         alpha = rwork(alphai)
         eta = rwork(etai)
         olmavg = rwork(olmavi)
         partol = rwork(partli)
         pnorm = rwork(pnormi)
         prers = rwork(prersi)
         rcond = rwork(rcondi)
         wss(1) = rwork(wssi)
         wss(2) = rwork(wssdei)
         wss(3) = rwork(wssepi)
         rvar = rwork(rvari)
         rnorms = rwork(rnorsi)
         sstol = rwork(sstoli)
         tau = rwork(taui)
         taufac = rwork(taufci)

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
         rwork(actrsi) = actrs
         rwork(alphai) = alpha
         rwork(olmavi) = olmavg
         rwork(partli) = partol
         rwork(pnormi) = pnorm
         rwork(prersi) = prers
         rwork(rcondi) = rcond
         rwork(wssi) = wss(1)
         rwork(wssdei) = wss(2)
         rwork(wssepi) = wss(3)
         rwork(rvari) = rvar
         rwork(rnorsi) = rnorms
         rwork(sstoli) = sstol
         rwork(taui) = tau

         iwork(iranki) = irank
         iwork(istopi) = istop
         iwork(nfevi) = nfev
         iwork(niteri) = niter
         iwork(njevi) = njev
         iwork(idfi) = idf
         iwork(int2i) = int2
      end if

   end subroutine access_workspace

   real(wp) pure function derstep(itype, k, betak, ssf, stpb, neta) result(res)
   !! Compute step size for center and forward difference calculations.

      use odrpack_kinds, only: zero, one

      integer, intent(in) :: itype
         !! Finite difference method being used, where: `itype = 0` indicates forward
         !! finite differences, and `itype = 1` indicates central finite differences.
      integer, intent(in) :: k
         !! Index into `beta` where `betak` resides.
      real(wp), intent(in) :: betak
         !! `k`-th function parameter.
      real(wp), intent(in) :: ssf(k)
         !! Scale used for the `beta`s.
      real(wp), intent(in) :: stpb(k)
         !! Relative step used for computing finite difference derivatives with respect
         !! to `beta`.
      integer, intent(in) :: neta
         !! Number of good digits in the function results.

      ! Local scalars
      real(wp) :: typj

      ! Variable definitions (alphabetically)
      !  TYPJ:    The typical size of the J-th unkonwn BETA.

      if (betak == zero) then
         if (ssf(1) < zero) then
            typj = 1/abs(ssf(1))
         else
            typj = 1/ssf(k)
         end if
      else
         typj = abs(betak)
      end if
      res = sign(one, betak)*typj*hstep(itype, neta, 1, k, stpb, 1)

   end function derstep

   pure subroutine esubi(n, m, wd, ldwd, ld2wd, alpha, tt, ldtt, i, e)
   !! Compute `e = wd + alpha*tt**2`.

      use odrpack_kinds, only: zero

      integer, intent(in) :: n
         !! Number of observations.
      integer, intent(in) :: m
         !! Number of columns of data in the independent variable.
      real(wp), intent(in) :: wd(ldwd, ld2wd, m)
         !! Squared `delta` weights.
      integer, intent(in) :: ldwd
         !! Leading dimension of array `wd`.
      integer, intent(in) :: ld2wd
         !! Second dimension of array `wd`.
      real(wp), intent(in) :: alpha
         !! Levenberg-Marquardt parameter.
      real(wp), intent(in) :: tt(ldtt, m)
         !! Scaling values used for `delta`.
      integer, intent(in) :: ldtt
         !! Leading dimension of array `tt`.
      integer, intent(in) :: i
         !! Indexing variable.
      real(wp), intent(out) :: e(m, m)
         !! Value of the array `e = wd + alpha*tt**2`.

      ! Local scalars
      integer :: j

      ! Variable Definitions (alphabetically)
      !  J:      An indexing variable.

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
               e = wd(i, :, :)
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
               e = wd(1, :, :)
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

   end subroutine esubi

   subroutine eta_fcn &
      (fcn, &
       n, m, np, q, &
       xplusd, beta, epsmac, nrow, &
       partmp, pv0, &
       ifixb, ifixx, ldifx, &
       istop, nfev, eta, neta, &
       fjacd, f, fjacb, wrk7, &
       info, &
       lower, upper)
   !! Compute noise and number of good digits in function results.
      ! Adapted from STARPAC subroutine ETAFUN.

      use odrpack_kinds, only: zero, one, two

      procedure(fcn_t) :: fcn
         !! User-supplied subroutine for evaluating the model.
      integer, intent(in) :: n
         !! Number of observations.
      integer, intent(in) :: m
         !! Number of columns of data in the explanatory variable.
      integer, intent(in) :: np
         !! Number of function parameters.
      integer, intent(in) :: q
         !! Number of responses per observation.
      real(wp), intent(in) :: xplusd(n, m)
         !! Values of `x + delta`.
      real(wp), intent(in) :: beta(np)
         !! Function parameters.
      real(wp), intent(in) :: epsmac
         !! Value of machine precision.
      integer, intent(in) :: nrow
         !! Row number at which the derivative is to be checked.
      real(wp), intent(out) :: partmp(np)
         !! Model parameters.
      real(wp), intent(in) :: pv0(n, q)
         !! Original predicted values.
      integer, intent(in) :: ifixb(np)
         !! Values designating whether the elements of `beta` are fixed at their input values or not.
      integer, intent(in) :: ifixx(ldifx, m)
         !! Values designating whether the elements of `x` are fixed at their input values or not.
      integer, intent(in) :: ldifx
         !! Leading dimension of array `ifixx`.
      integer, intent(out) :: istop
         !! Variable designating whether there are problems computing the function at the
         !! current `beta` and `delta`.
      integer, intent(inout) :: nfev
         !! Number of function evaluations.
      real(wp), intent(out) :: eta
         !! Noise in the model results.
      integer, intent(out) :: neta
         !! Number of accurate digits in the model results.
      real(wp), intent(out) :: fjacd(n, m, q)
         !! Jacobian wrt `delta`.
      real(wp), intent(out) :: f(n, q)
         !! Function values.
      real(wp), intent(out) :: fjacb(n, np, q)
         !! Jacobian wrt `beta`.
      real(wp), intent(out) :: wrk7(-2:2, q)
         !! Work array of `(5, q)` elements.
      integer, intent(out) :: info
         !! Variable indicating the status of the computation.
      real(wp), intent(in) :: lower(np)
         !! Lower bound of `beta`.
      real(wp), intent(in) :: upper(np)
         !! Upper bound of `beta`.

      ! Local scalars
      real(wp), parameter :: p1 = 0.1_wp, p2 = 0.2_wp, p5 = 0.5_wp
      real(wp) :: a, b, fac, shift, stp
      integer :: j, k, sbk

      ! Local arrays
      real(wp) :: parpts(-2:2, np)

      ! Variable Definitions (ALPHABETICALLY)
      !  A:       Parameters of the local fit.
      !  B:       Parameters of the local fit.
      !  FAC:     A factor used in the computations.
      !  J:       An index variable.
      !  K:       An index variable.
      !  SHIFT:   When PARPTS cross the parameter bounds they are shifted by SHIFT.
      !  SBK:     The sign of BETA(K).
      !  STP:     A small value used to perturb the parameters.

      stp = 100*epsmac
      eta = epsmac

      ! Create points to use in calculating FCN for ETA and NETA
      do j = -2, 2
         if (j == 0) then
            parpts(0, :) = beta
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
         if (all(parpts(j, :) == beta)) then
            wrk7(j, :) = pv0(nrow, :)
         else
            partmp = parpts(j, :)
            istop = 0
            call fcn(partmp, xplusd, ifixb, ifixx, 003, f, fjacb, fjacd, istop)
            if (istop /= 0) then
               return
            else
               nfev = nfev + 1
            end if
            wrk7(j, :) = f(nrow, :)
         end if
      end do

      ! Calculate ETA and NETA
      do k = 1, q
         a = zero
         b = zero
         do j = -2, 2
            a = a + wrk7(j, k)
            b = b + j*wrk7(j, k)
         end do
         a = p2*a
         b = p1*b
         if ((wrk7(0, k) /= zero) .and. (abs(wrk7(1, k) + wrk7(-1, k)) > 100*epsmac)) then
            fac = 1/abs(wrk7(0, k))
         else
            fac = one
         end if
         do j = -2, 2
            wrk7(j, k) = abs((wrk7(j, k) - (a + j*b))*fac)
            eta = max(wrk7(j, k), eta)
         end do
      end do
      neta = int(max(two, p5 - log10(eta)))

   end subroutine eta_fcn

   subroutine eval_jac &
      (fcn, &
       anajac, cdjac, &
       n, m, np, q, &
       betac, beta, stpb, &
       ifixb, ifixx, ldifx, &
       x, delta, xplusd, stpd, ldstpd, &
       ssf, tt, ldtt, neta, fn, &
       stp, wrk1, wrk2, wrk3, wrk6, tempret, &
       fjacb, isodr, fjacd, we1, ldwe, ld2we, &
       njev, nfev, istop, info, &
       lower, upper)
   !! Compute the weighted Jacobians wrt `beta` and `delta`.

      use odrpack_kinds, only: zero
      use blas_interfaces, only: ddot

      procedure(fcn_t) :: fcn
         !! User-supplied subroutine for evaluating the model.
      logical, intent(in) :: anajac
         !! Variable designating whether the Jacobians are computed by finite differences
         !! (`anajac = .false.`) or not (`anajac = .true.`).
      logical, intent(in) :: cdjac
         !! Variable designating whether the Jacobians are computed by central differences
         !! (`cdjac = .true.`) or by forward differences (`cdjac = .false.`).
      integer, intent(in) :: n
         !! Number of observations.
      integer, intent(in) :: m
         !! Number of columns of data in the independent variable.
      integer, intent(in) :: np
         !! Number of function parameters.
      integer, intent(in) :: q
         !! Number of responses per observation.
      real(wp), intent(in) :: betac(np)
         !! Current estimated values of the unfixed `beta`s.
      real(wp), intent(out) :: beta(np)
         !! Function parameters.
      real(wp), intent(in) :: stpb(np)
         !! Relative step used for computing finite difference derivatives with respect to `beta`.
      integer, intent(in) :: ifixb(np)
         !! Values designating whether the elements of `beta` are fixed at their input values or not.
      integer, intent(in) :: ifixx(ldifx, m)
         !! Values designating whether the elements of `delta` are fixed at their input values or not.
      integer, intent(in) :: ldifx
         !! Leading dimension of array `ifixx`.
      real(wp), intent(in) :: x(n, m)
         !! Independent variable.
      real(wp), intent(in) :: delta(n, m)
         !! Estimated values of `delta`.
      real(wp), intent(out) :: xplusd(n, m)
         !! Values of `x + delta`.
      real(wp), intent(in) :: stpd(ldstpd, m)
         !! Relative step used for computing finite difference derivatives with respect to `delta`.
      integer, intent(in) :: ldstpd
         !! Leading dimension of array `stpd`.
      real(wp), intent(in) :: ssf(np)
         !! Scale used for the `beta`s.
      real(wp), intent(in) :: tt(ldtt, m)
         !! Scaling values used for `delta`.
      integer, intent(in) :: ldtt
         !! Leading dimension of array `tt`.
      integer, intent(in) :: neta
         !! Number of accurate digits in the function results.
      real(wp), intent(in) :: fn(n, q)
         !! Predicted values of the function at the current point.
      real(wp), intent(out) :: stp(n)
         !! Step used for computing finite difference derivatives with respect to `delta`.
      real(wp), intent(out) :: wrk1(n, m, q)
         !! Work array of `(n, m, q)` elements.
      real(wp), intent(out) :: wrk2(n, q)
         !! Work array of `(n, q)` elements.
      real(wp), intent(out) :: wrk3(np)
         !! Work array of `(np)` elements.
      real(wp), intent(out) :: wrk6(n, np, q)
         !! Work array of `(n, np, q)` elements.
      real(wp), intent(inout) :: tempret(:, :)
         !! Temporary work array for holding return values before copying to a lower rank array.
      real(wp), intent(out) :: fjacb(n, np, q)
         !! Jacobian with respect to `beta`.
      logical, intent(in) :: isodr
         !! Variable designating whether the solution is by ODR (`isodr = .true.`) or
         !! by OLS (`isodr = .false.`).
      real(wp), intent(out) :: fjacd(n, m, q)
         !! Jacobian with respect to `delta`.
      real(wp), intent(in) :: we1(ldwe, ld2we, q)
         !! Square roots of the `epsilon` weights in array `we`.
      integer, intent(in) :: ldwe
         !! Leading dimension of arrays `we` and `we1`.
      integer, intent(in) :: ld2we
         !! Second dimension of arrays `we` and `we1`.
      integer, intent(inout) :: njev
         !! Number of Jacobian evaluations.
      integer, intent(inout) :: nfev
         !! Number of function evaluations.
      integer, intent(out) :: istop
         !! Variable designating that the user wishes the computations stopped.
      integer, intent(out) :: info
         !! Variable designating why the computations were stopped.
      real(wp), intent(in) :: lower(np)
         !! Lower bound of `beta`.
      real(wp), intent(in) :: upper(np)
         !! Upper bound of `beta`.

      ! Local scalars
      integer :: ideval, j, j1
      logical :: ferror

      ! Variable Definitions (alphabetically)
      !  FERROR:  The variable designating whether ODRPACK95 detected nonzero values in array DELTA
      !           in the OLS case, and thus whether the user may have overwritten important information
      !           by computing FJACD in the OLS case.
      !  IDEVAL:  The variable designating what computations are to be performed by user-supplied
      !           subroutine FCN.
      !  J:       An indexing variable.
      !  J1:      An indexing variable.

      ! Insert current unfixed BETA estimates into BETA
      call unpack_vec(np, betac, beta, ifixb)

      ! Compute XPLUSD = X + DELTA
      xplusd = x + delta

      ! Compute the Jacobian wrt the estimated BETAS (FJACB) and the Jacobian wrt DELTA (FJACD)
      istop = 0
      if (isodr) then
         ideval = 110
      else
         ideval = 010
      end if
      if (anajac) then
         call fcn(beta, xplusd, ifixb, ifixx, ideval, wrk2, fjacb, fjacd, istop)
         if (istop /= 0) then
            return
         else
            njev = njev + 1
         end if
         ! Make sure fixed elements of FJACD are zero
         if (isodr) then
            do j = 1, q
               call set_ifix(n, m, ifixx, ldifx, fjacd(1, 1, j), n, fjacd(1, 1, j), n)
            end do
         end if
      elseif (cdjac) then
         call jac_cdiff(fcn, &
                        n, m, np, q, &
                        beta, x, delta, xplusd, ifixb, ifixx, ldifx, &
                        stpb, stpd, ldstpd, &
                        ssf, tt, ldtt, neta, fn, stp, wrk1, wrk2, wrk3, wrk6, &
                        fjacb, isodr, fjacd, nfev, istop, info, &
                        lower, upper)
      else
         call jac_fwdiff(fcn, &
                         n, m, np, q, &
                         beta, x, delta, xplusd, ifixb, ifixx, ldifx, &
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
         do j = 1, np
            call scale_mat(n, q, we1, ldwe, ld2we, &
                           reshape(fjacb(:, j, :), [n, q]), &
                           tempret(1:n, 1:q))
            fjacb(:, j, :) = tempret(1:n, 1:q)
         end do
      else
         j1 = 0
         do j = 1, np
            if (ifixb(j) >= 1) then
               j1 = j1 + 1
               call scale_mat(n, q, we1, ldwe, ld2we, &
                              reshape(fjacb(:, j, :), [n, q]), &
                              tempret(1:n, 1:q))
               fjacb(:, j1, :) = tempret(1:n, 1:q)
            end if
         end do
      end if

      ! Weight the Jacobian's wrt DELTA as appropriate
      if (isodr) then
         do j = 1, m
            call scale_mat(n, q, we1, ldwe, ld2we, &
                           reshape(fjacd(:, j, :), [n, q]), &
                           tempret(1:n, 1:q))
            fjacd(:, j, :) = tempret(1:n, 1:q)
         end do
      end if

   end subroutine eval_jac

   pure subroutine fctr(oksemi, a, lda, n, info)
   !! Factor the positive (semi)definite matrix `a` using a modified Cholesky factorization.
      ! Adapted from LINPACK subroutine DPOFA.

      use odrpack_kinds, only: zero

      logical, intent(in) :: oksemi
         !! Flag indicating whether the factored array can be positive semidefinite
         !! (`oksemi = .true.`) or whether it must be found to be positive definite
         !! (`oksemi = .false.`).
      real(wp), intent(inout) :: a(lda, n)
         !! Array to be factored. Upon return, `a` contains the upper triangular matrix
         !! `r` so that `a = trans(r)*r` where the strict lower triangle is set to zero.
         !! If `info /= 0`, the factorization is not complete.
      integer, intent(in) :: lda
         !! Leading dimension of array `a`.
      integer, intent(in) :: n
         !! Number of rows and columns of data in array `a`.
      integer, intent(out) :: info
         !! Output flag.
         !! `0`: Factorization was completed.
         !! `k`: Error condition. The leading minor of order `k` is not positive (semi)definite.

      ! Local scalars
      real(wp) :: xi, s, t
      integer j, k

      ! Variable Definitions (alphabetically)
      !  J:       An indexing variable.
      !  K:       An indexing variable.
      !  XI:      A value used to test for non positive semidefiniteness.

      ! Set relative tolerance for detecting non positive semidefiniteness.
      xi = -10*epsilon(zero)

      ! Compute factorization, storing in upper triangular portion of A
      do j = 1, n
         info = j
         s = zero
         do k = 1, j - 1
            if (a(k, k) == zero) then
               t = zero
            else
               t = a(k, j) - dot_product(a(1:k - 1, k), a(1:k - 1, j))
               t = t/a(k, k)
            end if
            a(k, j) = t
            s = s + t**2
         end do
         s = a(j, j) - s

         ! Exit
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
      do j = 1, n - 1
         a(j + 1:n, j) = zero
      end do

   end subroutine fctr

   pure subroutine fctrw &
      (n, m, q, npp, &
       isodr, &
       we, ldwe, ld2we, wd, ldwd, ld2wd, &
       wrk0, wrk4, &
       we1, nnzw, info)
   !! Check input parameters, indicating errors found using nonzero values of argument `info`
   !! as described in the ODRPACK95 reference guide.

      use odrpack_kinds, only: zero

      integer, intent(in) :: n
         !! Number of observations.
      integer, intent(in) :: m
         !! Number of columns of data in the explanatory variable.
      integer, intent(in) :: q
         !! Number of responses per observation.
      integer, intent(in) :: npp
         !! Number of function parameters being estimated.
      logical, intent(in) :: isodr
         !! Variable designating whether the solution is by ODR (`isodr = .true`) or
         !! by OLS (`isodr = .false`).
      real(wp), intent(in) :: we(ldwe, ld2we, q)
         !! Squared `epsilon` weights.
      integer, intent(in) :: ldwe
         !! Leading dimension of array `we`.
      integer, intent(in) :: ld2we
         !! Second dimension of array `we`.
      real(wp), intent(in) :: wd(ldwd, ld2wd, m)
         !! Squared `delta` weights.
      integer, intent(in) :: ldwd
         !! Leading dimension of array `wd`.
      integer, intent(in) :: ld2wd
         !! Second dimension of array `wd`.
      real(wp), intent(out) :: wrk0(q, q)
         !! Work array of `(q, q)` elements.
      real(wp), intent(out) :: wrk4(m, m)
         !! Work array of `(m, m)` elements.
      real(wp), intent(out) :: we1(ldwe, ld2we, q)
         !! Factored `epsilon` weights, such that `trans(we1)*we1 = we`.
      integer, intent(out) :: nnzw
         !! Number of nonzero weighted observations.
      integer, intent(out) :: info
         !! Variable designating why the computations were stopped.

      ! Local scalars
      integer :: i, finfo, j
      logical :: notzero, exited

      ! Variable Definitions (alphabetically)
      !  I:        An indexing variable.
      !  J:        An indexing variable.
      !  L:        An indexing variable.
      !  NOTZERO:  The variable designating whether a given component of the weight array WE
      !            contains a nonzero element (NOTZRO=FALSE) or not (NOTZRO=TRUE).

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
               do j = 1, q
                  if (we(1, 1, j) > zero) then
                     nnzw = n
                     we1(1, 1, j) = sqrt(we(1, 1, j))
                  elseif (we(1, 1, j) < zero) then
                     info = 30010
                     exited = .true.
                     exit
                  end if
               end do
            else
               ! WE contains a full Q by Q semidefinite matrix
               do j = 1, q
                  wrk0(1:j, j) = we(1, 1:j, j)
               end do
               call fctr(.true., wrk0, q, q, finfo)
               if (finfo /= 0) then
                  info = 30010
                  exited = .true.
               else
                  do j = 1, q
                     we1(1, :, j) = wrk0(:, j)
                     if (we1(1, j, j) /= zero) then
                        nnzw = n
                     end if
                  end do
               end if
            end if
         else
            if (ld2we == 1) then
               ! WE contains an array of  diagonal matrix
               do i = 1, n
                  notzero = .false.
                  do j = 1, q
                     if (we(i, 1, j) > zero) then
                        notzero = .true.
                        we1(i, 1, j) = sqrt(we(i, 1, j))
                     elseif (we(i, 1, j) < zero) then
                        info = 30010
                        exited = .true.
                        exit
                     end if
                  end do
                  if (exited) exit
                  if (notzero) then
                     nnzw = nnzw + 1
                  end if
               end do
            else
               ! WE contains an array of full Q by Q semidefinite matrices
               do i = 1, n
                  do j = 1, q
                     wrk0(1:j, j) = we(i, 1:j, j)
                  end do
                  call fctr(.true., wrk0, q, q, finfo)
                  if (finfo /= 0) then
                     info = 30010
                     exited = .true.
                     exit
                  else
                     notzero = .false.
                     do j = 1, q
                        we1(i, :, j) = wrk0(:, j)
                        if (we1(i, j, j) /= zero) then
                           notzero = .true.
                        end if
                     end do
                  end if
                  if (notzero) then
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
               if (any(wd(1, 1, :) <= zero)) then
                  info = max(30001, info + 1)
                  return
               end if
            else
               ! WD contains a full M by M positive definite matrix
               do j = 1, m
                  wrk4(1:j, j) = wd(1, 1:j, j)
               end do
               call fctr(.false., wrk4, m, m, finfo)
               if (finfo /= 0) then
                  info = max(30001, info + 1)
                  return
               end if
            end if
         else
            if (ld2wd == 1) then
               ! WD contains an array of diagonal matrices
               if (any(wd(:, 1, :) <= zero)) then
                  info = max(30001, info + 1)
                  return
               end if
            else
               ! WD contains an array of full M by M positive definite matrices
               do i = 1, n
                  do j = 1, m
                     wrk4(1:j, j) = wd(i, 1:j, j)
                  end do
                  call fctr(.false., wrk4, m, m, finfo)
                  if (finfo /= 0) then
                     info = max(30001, info + 1)
                     return
                  end if
               end do
            end if
         end if
      end if

   end subroutine fctrw

   pure subroutine set_flags( &
      job, restrt, initd, dovcv, redoj, anajac, cdjac, chkjac, isodr, implct)
   !! Set flags indicating conditions specified by `job`.

      integer, intent(in) :: job
         !! Variable controlling problem initialization and computational method.
      logical, intent(out) :: restrt
         !! Variable designating whether the call is a restart (`restrt = .true.`)
         !! or not (`restrt = .false.`).
      logical, intent(out) :: initd
         !! Variable designating whether `delta` is to be initialized to zero (`initd = .true.`)
         !! or to the first `n` by `m` elements of array `rwork` (`initd = .false.`).
      logical, intent(out) :: dovcv
         !! Variable designating whether the covariance matrix is to be computed
         !! (`dovcv = .true.`) or not (`dovcv = .false.`).
      logical, intent(out) :: redoj
         !! Variable designating whether the Jacobian matrix is to be recomputed for the
         !! computation of the covariance matrix (`redoj = .true.`) or not (`redoj = .false.`).
      logical, intent(out) :: anajac
         !! Variable designating whether the Jacobians are computed by finite differences
         !! (`anajac = .false.`) or not (`anajac = .true.`).
      logical, intent(out) :: cdjac
         !! Variable designating whether the Jacobians are computed by central differences
         !! (`cdjac = .true.`) or by forward differences (`cdjac = .false.`).
      logical, intent(out) :: chkjac
         !! Variable designating whether the user-supplied Jacobians are to be checked
         !! (`chkjac = .true.`) or not (`chkjac = .false.`).
      logical, intent(out) :: isodr
         !! Variable designating whether the solution is by ODR (`isodr = .true.`)
         !! or by OLS (`isodr = .false.`).
      logical, intent(out) :: implct
         !! Variable designating whether the solution is by implicit ODR (`implct = .true.`)
         !! or explicit ODR (`implct = .false.`).

      ! Local scalars
      integer :: j

      ! Variable Definitions (alphabetically)
      !  J:       The value of a specific digit of JOB.

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

   end subroutine set_flags

   real(wp) pure function hstep(itype, neta, i, j, stp, ldstp) result(res)
   !! Set relative step size for finite difference derivatives.

      use odrpack_kinds, only: zero, two, three, ten

      integer, intent(in) :: itype
         !! Finite difference method being used, where: `itype = 0` indicates forward
         !! finite differences, and `itype = 1` indicates central finite differences.
      integer, intent(in) :: neta
         !! Number of good digits in the function results.
      integer, intent(in) :: i
         !! Identifier for selecting user-supplied step sizes.
      integer, intent(in) :: j
         !! Identifier for selecting user-supplied step sizes.
      real(wp), intent(in) :: stp(ldstp, j)
         !! Step size for the finite difference derivative.
      integer, intent(in) :: ldstp
         !! Leading dimension of array `stp`.

      if (stp(1, 1) <= zero) then
         if (itype == 0) then
            ! Use default forward finite difference step size
            res = ten**(-abs(neta)/two - two)
         else
            ! Use default central finite difference step size
            res = ten**(-abs(neta)/three)
         end if
      elseif (ldstp == 1) then
         res = stp(1, j)
      else
         res = stp(i, j)
      end if

   end function hstep

   pure subroutine set_ifix(n, m, ifix, ldifix, t, ldt, tfix, ldtfix)
   !! Set elements of `t` to zero according to `ifix`.

      use odrpack_kinds, only: zero

      integer, intent(in) :: n
         !! Number of rows of data in the array.
      integer, intent(in) :: m
         !! Number of columns of data in the array.
      integer, intent(in) :: ifix(ldifix, m)
         !! Array designating whether an element of `t` is to be set to zero.
      integer, intent(in) :: ldifix
         !! Leading dimension of array `ifix`.
      real(wp), intent(in) :: t(ldt, m)
         !! Array being set to zero according to the elements of `ifix`.
      integer, intent(in) :: ldt
         !! Leading dimension of array `t`.
      real(wp), intent(out) :: tfix(ldtfix, m)
         !! Resulting array.
      integer, intent(in) :: ldtfix
         !! Leading dimension of array `tfix`.

      ! Local scalars
      integer :: j

      ! Variable Definitions (alphabetically)
      !  J:       an indexing variable.

      if (n == 0 .or. m == 0) return

      if (ifix(1, 1) >= zero) then
         if (ldifix == n) then
            where (ifix(:, :) == 0)
               tfix(1:n, :) = zero
            else where
               tfix(1:n, :) = t(1:n, :)
            end where
         else
            do j = 1, m
               if (ifix(1, j) == 0) then
                  tfix(1:n, j) = zero
               else
                  tfix(1:n, j) = t(1:n, j)
               end if
            end do
         end if
      end if

   end subroutine set_ifix

   pure subroutine loc_iwork &
      (m, q, np, &
       msgbi, msgdi, ifix2i, istopi, &
       nnzwi, nppi, idfi, &
       jobi, iprini, luneri, lunrpi, &
       nrowi, ntoli, netai, &
       maxiti, niteri, nfevi, njevi, int2i, iranki, ldtti, &
       boundi, &
       liwkmn)
   !! Get storage locations within integer work space.

      integer, intent(in) :: m
         !! Number of columns of data in the independent variable.
      integer, intent(in) :: q
         !! Number of responses per observation.
      integer, intent(in) :: np
         !! Number of function parameters.
      integer, intent(out) :: msgbi
         !! Starting location in array `iwork` of array `msgb`.
      integer, intent(out) :: msgdi
         !! Starting location in array `iwork` of array `msgd`.
      integer, intent(out) :: ifix2i
         !! Starting location in array `iwork` of array `ifix2`.
      integer, intent(out) :: istopi
         !! Location in array `iwork` of variable `istop`.
      integer, intent(out) :: nnzwi
         !! Location in array `iwork` of variable `nnzw`.
      integer, intent(out) :: nppi
         !! Location in array `iwork` of variable `npp`.
      integer, intent(out) :: idfi
         !! Location in array `iwork` of variable `idf`.
      integer, intent(out) :: jobi
         !! Location in array `iwork` of variable `job`.
      integer, intent(out) :: iprini
         !! Location in array `iwork` of variable `iprint`.
      integer, intent(out) :: luneri
         !! Location in array `iwork` of variable `lunerr`.
      integer, intent(out) :: lunrpi
         !! Location in array `iwork` of variable `lunrpt`.
      integer, intent(out) :: nrowi
         !! Location in array `iwork` of variable `nrow`.
      integer, intent(out) :: ntoli
         !! Location in array `iwork` of variable `ntol`.
      integer, intent(out) :: netai
         !! Location in array `iwork` of variable `neta`.
      integer, intent(out) :: maxiti
         !! Location in array `iwork` of variable `maxit`.
      integer, intent(out) :: niteri
         !! Location in array `iwork` of variable `niter`.
      integer, intent(out) :: nfevi
         !! Location in array `iwork` of variable `nfev`.
      integer, intent(out) :: njevi
         !! Location in array `iwork` of variable `njev`.
      integer, intent(out) :: int2i
         !! Location in array `iwork` of variable `int2`.
      integer, intent(out) :: iranki
         !! Location in array `iwork` of variable `irank`.
      integer, intent(out) :: ldtti
         !! Location in array `iwork` of variable `ldtt`.
      integer, intent(out) :: boundi
         !! Location in array `iwork` of variable `bound`.
      integer, intent(out) :: liwkmn
         !! Minimum acceptable length of array `iwork`.

      if (np >= 1 .and. m >= 1) then
         msgbi = 1
         msgdi = msgbi + q*np + 1
         ifix2i = msgdi + q*m + 1
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

   end subroutine loc_iwork

   pure subroutine loc_rwork &
      (n, m, q, np, ldwe, ld2we, isodr, &
       deltai, epsi, xplusi, fni, sdi, vcvi, &
       rvari, wssi, wssdei, wssepi, rcondi, etai, &
       olmavi, taui, alphai, actrsi, pnormi, rnorsi, prersi, &
       partli, sstoli, taufci, epsmai, &
       beta0i, betaci, betasi, betani, si, ssi, ssfi, qrauxi, ui, &
       fsi, fjacbi, we1i, diffi, &
       deltsi, deltni, ti, tti, omegai, fjacdi, &
       wrk1i, wrk2i, wrk3i, wrk4i, wrk5i, wrk6i, wrk7i, &
       loweri, upperi, &
       lrwkmn)
   !! Get storage locations within real work space.

      integer, intent(in) :: n
         !! Number of observations.
      integer, intent(in) :: m
         !! Number of columns of data in the explanatory variable.
      integer, intent(in) :: q
         !! Number of responses per observation.
      integer, intent(in) :: np
         !! Number of function parameters.
      integer, intent(in) :: ldwe
         !! Leading dimension of array `we`.
      integer, intent(in) :: ld2we
         !! Second dimension of array `we`.
      logical, intent(in) :: isodr
         !! Variable designating whether the solution is by ODR (`isodr`=.true.) or by OLS (`isodr`=.false.).
      integer, intent(out) :: deltai
         !! Starting location in array `rwork` of array `delta`.
      integer, intent(out) :: epsi
         !! Starting location in array `rwork` of array `eps`.
      integer, intent(out) :: xplusi
         !! Starting location in array `rwork` of array `xplusd`.
      integer, intent(out) :: fni
         !! Starting location in array `rwork` of array `fn`.
      integer, intent(out) :: sdi
         !! Starting location in array `rwork` of array `sd`.
      integer, intent(out) :: vcvi
         !! Starting location in array `rwork` of array `vcv`.
      integer, intent(out) :: rvari
         !! Location in array `rwork` of variable `rvar`.
      integer, intent(out) :: wssi
         !! Location in array `rwork` of variable `wss`.
      integer, intent(out) :: wssdei
         !! Location in array `rwork` of variable `wssdel`.
      integer, intent(out) :: wssepi
         !! Location in array `rwork` of variable `wsseps`.
      integer, intent(out) :: rcondi
         !! Location in array `rwork` of variable `rcondi`.
      integer, intent(out) :: etai
         !! Location in array `rwork` of variable `eta`.
      integer, intent(out) :: olmavi
         !! Location in array `rwork` of variable `olmavg`.
      integer, intent(out) :: taui
         !! Location in array `rwork` of variable `tau`.
      integer, intent(out) :: alphai
         !! Location in array `rwork` of variable `alpha`.
      integer, intent(out) :: actrsi
         !! Location in array `rwork` of variable `actrs`.
      integer, intent(out) :: pnormi
         !! Location in array `rwork` of variable `pnorm`.
      integer, intent(out) :: rnorsi
         !! Location in array `rwork` of variable `rnorms`.
      integer, intent(out) :: prersi
         !! Location in array `rwork` of variable `prers`.
      integer, intent(out) :: partli
         !! Location in array `rwork` of variable `partol`.
      integer, intent(out) :: sstoli
         !! Location in array `rwork` of variable `sstol`.
      integer, intent(out) :: taufci
         !! Location in array `rwork` of variable `taufac`.
      integer, intent(out) :: epsmai
         !! Location in array `rwork` of variable `epsmac`.
      integer, intent(out) :: beta0i
         !! Starting location in array `rwork` of array `beta0`.
      integer, intent(out) :: betaci
         !! Starting location in array `rwork` of array `betac`.
      integer, intent(out) :: betasi
         !! Starting location in array `rwork` of array `betas`.
      integer, intent(out) :: betani
         !! Starting location in array `rwork` of array `betan`.
      integer, intent(out) :: si
         !! Starting location in array `rwork` of array `s`.
      integer, intent(out) :: ssi
         !! Starting location in array `rwork` of array `ss`.
      integer, intent(out) :: ssfi
         !! Starting location in array `rwork` of array `ssf`.
      integer, intent(out) :: qrauxi
         !! Starting location in array `rwork` of array `qraux`.
      integer, intent(out) :: ui
         !! Starting location in array `rwork` of array `u`.
      integer, intent(out) :: fsi
         !! Starting location in array `rwork` of array `fs`.
      integer, intent(out) :: fjacbi
         !! Starting location in array `rwork` of array `fjacb`.
      integer, intent(out) :: we1i
         !! Starting location in array `rwork` of array `we1`.
      integer, intent(out) :: diffi
         !! Starting location in array `rwork` of array `diff`.
      integer, intent(out) :: deltsi
         !! Starting location in array `rwork` of array `deltas`.
      integer, intent(out) :: deltni
         !! Starting location in array `rwork` of array `deltan`.
      integer, intent(out) :: ti
         !! Starting location in array `rwork` of array `t`.
      integer, intent(out) :: tti
         !! Starting location in array `rwork` of array `tt`.
      integer, intent(out) :: omegai
         !! Starting location in array `rwork` of array `omega`.
      integer, intent(out) :: fjacdi
         !! Starting location in array `rwork` of array `fjacd`.
      integer, intent(out) :: wrk1i
         !! Starting location in array `rwork` of array `wrk1`.
      integer, intent(out) :: wrk2i
         !! Starting location in array `rwork` of array `wrk2`.
      integer, intent(out) :: wrk3i
         !! Starting location in array `rwork` of array `wrk3`.
      integer, intent(out) :: wrk4i
         !! Starting location in array `rwork` of array `wrk4`.
      integer, intent(out) :: wrk5i
         !! Starting location in array `rwork` of array `wrk5`.
      integer, intent(out) :: wrk6i
         !! Starting location in array `rwork` of array `wrk6`.
      integer, intent(out) :: wrk7i
         !! Starting location in array `rwork` of array `wrk7`.
      integer, intent(out) :: loweri
         !! Starting location in array `rwork` of array `lower`.
      integer, intent(out) :: upperi
         !! Starting location in array `rwork` of array `upper`.
      integer, intent(out) :: lrwkmn
         !! Minimum acceptable length of vector `rwork`.

      ! Local scalars
      integer :: next

      ! Variable Definitions (alphabetically)
      !  NEXT:    The next available location with RWORK.

      if (n >= 1 .and. m >= 1 .and. np >= 1 .and. q >= 1 .and. ldwe >= 1 &
          .and. ld2we >= 1) then

         deltai = 1
         epsi = deltai + n*m
         xplusi = epsi + n*q
         fni = xplusi + n*m
         sdi = fni + n*q
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

         fjacbi = fsi + n*q

         we1i = fjacbi + n*np*q

         diffi = we1i + ldwe*ld2we*q

         next = diffi + q*(np + m)

         if (isodr) then
            deltsi = next
            deltni = deltsi + n*m
            ti = deltni + n*m
            tti = ti + n*m
            omegai = tti + n*m
            fjacdi = omegai + q**2
            wrk1i = fjacdi + n*m*q
            next = wrk1i + n*m*q
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
         wrk3i = wrk2i + n*q
         wrk4i = wrk3i + np
         wrk5i = wrk4i + m*m
         wrk6i = wrk5i + m
         wrk7i = wrk6i + n*q*np
         loweri = wrk7i + 5*q
         upperi = loweri + np
         next = upperi + np

         lrwkmn = next
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
         lrwkmn = 1
      end if

   end subroutine loc_rwork

   pure subroutine init_work &
      (n, m, np, rwork, lrwork, iwork, liwork, &
       x, ifixx, ldifx, scld, ldscld, &
       beta, sclb, &
       sstol, partol, maxit, taufac, &
       job, iprint, lunerr, lunrpt, &
       lower, upper, &
       epsmai, sstoli, partli, maxiti, taufci, &
       jobi, iprini, luneri, lunrpi, &
       ssfi, tti, ldtti, deltai, &
       loweri, upperi, boundi)
   !! Initialize work vectors as necessary.

      use odrpack_kinds, only: zero, one, two, three

      integer, intent(in) :: n
         !! Number of observations.
      integer, intent(in) :: m
         !! Number of columns of data in the independent variable.
      integer, intent(in) :: np
         !! Number of function parameters.
      real(wp), intent(out) :: rwork(lrwork)
         !! Real work space.
      integer, intent(in) :: lrwork
         !! Length of vector `rwork`.
      integer, intent(out) :: iwork(liwork)
         !! Integer work space.
      integer, intent(in) :: liwork
         !! Length of vector `iwork`.
      real(wp), intent(in) :: x(n, m)
         !! Independent variable.
      integer, intent(in) :: ifixx(ldifx, m)
         !! Values designating whether the elements of `x` are fixed at their input values or not.
      integer, intent(in) :: ldifx
         !! Leading dimension of array `ifixx`.
      real(wp), intent(in) :: scld(ldscld, m)
         !! Scaling values for `delta`.
      integer, intent(in) :: ldscld
         !! Leading dimension of array `scld`.
      real(wp), intent(in) :: beta(np)
         !! Function parameters.
      real(wp), intent(in) :: sclb(np)
         !! Scaling values for `beta`.
      real(wp), intent(in) :: sstol
         !! Sum-of-squares convergence stopping criteria.
      real(wp), intent(in) :: partol
         !! Parameter convergence stopping criteria.
      integer, intent(in) :: maxit
         !! Maximum number of iterations allowed.
      real(wp), intent(in) :: taufac
         !! Factor used to compute the initial trust region diameter.
      integer, intent(in) :: job
         !! Variable controlling problem initialization and computational method.
      integer, intent(in) :: iprint
         !! Print control variable.
      integer, intent(in) :: lunerr
         !! Logical unit number used for error messages.
      integer, intent(in) :: lunrpt
         !! Logical unit number used for computation reports.
      real(wp), intent(in) :: lower(np)
         !! Lower bounds for the function parameters.
      real(wp), intent(in) :: upper(np)
         !! Upper bounds for the function parameters.
      integer, intent(in) :: epsmai
         !! Location in array `rwork` of variable `epsmac`.
      integer, intent(in) :: sstoli
         !! Location in array `rwork` of variable `sstol`.
      integer, intent(in) :: partli
         !! Location in array `rwork` of variable `partol`.
      integer, intent(in) :: maxiti
         !! Location in array `iwork` of variable `maxit`.
      integer, intent(in) :: taufci
         !! Location in array `rwork` of variable `taufac`.
      integer, intent(in) :: jobi
         !! Location in array `iwork` of variable `job`.
      integer, intent(in) :: iprini
         !! Location in array `iwork` of variable `iprint`.
      integer, intent(in) :: luneri
         !! Location in array `iwork` of variable `lunerr`.
      integer, intent(in) :: lunrpi
         !! Location in array `iwork` of variable `lunrpt`.
      integer, intent(in) :: ssfi
         !! Starting location in array `rwork` of array `ssf`.
      integer, intent(in) :: tti
         !! Starting location in array `rwork` of the array `tt`.
      integer, intent(in) :: ldtti
         !! Leading dimension of array `tt`.
      integer, intent(in) :: deltai
         !! Starting location in array `rwork` of array `delta`.
      integer, intent(in) :: loweri
         !! Starting location in array `iwork` of array `lower`.
      integer, intent(in) :: upperi
         !! Starting location in array `iwork` of array `upper`.
      integer, intent(in) :: boundi
         !! Location in array `iwork` of variable `bound`.

      ! Local scalars
      integer :: i, j, istart
      logical :: anajac, cdjac, chkjac, dovcv, implct, initd, isodr, redoj, restrt

      ! Variable Definitions (alphabetically)
      !  ANAJAC:  The variable designating whether the Jacobians are computed by finite differences
      !           (ANAJAC=FALSE) or not (ANAJAC=TRUE).
      !  CDJAC:   The variable designating whether the Jacobians are computed by central differences
      !           (CDJAC=TRUE) or by forward differences (CDJAC=FALSE).
      !  CHKJAC:  The variable designating whether the user-supplied Jacobians are to be checked
      !           (CHKJAC=TRUE) or not (CHKJAC=FALSE).
      !  DOVCV:   The variable designating whether the covariance matrix is to be computed (DOVCV=TRUE)
      !           or not (DOVCV=FALSE).
      !  I:       An indexing variable.
      !  IMPLCT:  The variable designating whether the solution is by implicit ODR (IMPLCT=TRUE) or
      !           explicit ODR (IMPLCT=FALSE).
      !  INITD:   The variable designating whether DELTA is to be initialized to zero (INITD=TRUE) or
      !           to the values in the first N by M elements of array RWORK (INITD=FALSE).
      !  ISODR:   The variable designating whether the solution is by ODR (ISODR=TRUE) or by OLS
      !           ISODR=FALSE).
      !  J:       An indexing variable.
      !  REDOJ:   The variable designating whether the Jacobian matrix is to be recomputed for the
      !           computation of the covariance matrix (REDOJ=TRUE) or not (REDOJ=FALSE).
      !  RESTRT:  The variable designating whether the call is a restart (RESTRT=TRUE) or not
      !           (RESTRT=FALSE).

      call set_flags(job, restrt, initd, dovcv, redoj, anajac, cdjac, chkjac, isodr, implct)

      ! Store value of machine precision in work vector
      rwork(epsmai) = epsilon(zero)

      ! Set tolerance for stopping criteria based on the change in the parameters (see also
      ! subprogram DODCNT)
      if (partol < zero) then
         rwork(partli) = rwork(epsmai)**(two/three)
      else
         rwork(partli) = min(partol, one)
      end if

      ! Set tolerance for stopping criteria based on the change in the sum of squares of the
      ! weighted observational errors
      if (sstol < zero) then
         rwork(sstoli) = sqrt(rwork(epsmai))
      else
         rwork(sstoli) = min(sstol, one)
      end if

      ! Set factor for computing trust region diameter at first iteration
      if (taufac <= zero) then
         rwork(taufci) = one
      else
         rwork(taufci) = min(taufac, one)
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
         call scale_beta(np, beta, rwork(ssfi))
      else
         rwork(ssfi:ssfi + (np - 1)) = sclb
      end if

      if (isodr) then
         if (scld(1, 1) <= zero) then
            iwork(ldtti) = n
            call scale_delta(n, m, x, rwork(tti), iwork(ldtti))
         else
            if (ldscld == 1) then
               iwork(ldtti) = 1
               rwork(tti:tti + (m - 1)) = scld(1:m, 1)
            else
               iwork(ldtti) = n
               do j = 1, m
                  istart = tti + (j - 1)*iwork(ldtti)
                  rwork(istart:istart + (n - 1)) = scld(1:n, j)
               end do
            end if
         end if
      end if

      ! Initialize DELTA's as necessary
      if (isodr) then
         if (initd) then
            rwork(deltai:deltai + (n*m - 1)) = zero
         else
            if (ifixx(1, 1) >= 0) then
               if (ldifx == 1) then
                  do j = 1, m
                     if (ifixx(1, j) == 0) then
                        istart = deltai + (j - 1)*n
                        rwork(istart:istart + (n - 1)) = zero
                     end if
                  end do
               else
                  do j = 1, m
                     do i = 1, n
                        if (ifixx(i, j) == 0) then
                           rwork(deltai - 1 + i + (j - 1)*n) = zero
                        end if
                     end do
                  end do
               end if
            end if
         end if
      else
         rwork(deltai:deltai + (n*m - 1)) = zero
      end if

      ! Copy bounds into RWORK
      rwork(loweri:loweri + np - 1) = lower
      rwork(upperi:upperi + np - 1) = upper

      ! Initialize parameters on bounds in IWORK
      iwork(boundi:boundi + np - 1) = 0

   end subroutine init_work

   subroutine jac_cdiff &
      (fcn, &
       n, m, np, q, &
       beta, x, delta, xplusd, ifixb, ifixx, ldifx, &
       stpb, stpd, ldstpd, &
       ssf, tt, ldtt, neta, fn, stp, wrk1, wrk2, wrk3, wrk6, &
       fjacb, isodr, fjacd, nfev, istop, info, &
       lower, upper)
   !! Compute central difference approximations to the Jacobian wrt the estimated `beta`s and
   !! wrt the `delta`s.

      use odrpack_kinds, only: zero, one

      procedure(fcn_t) :: fcn
         !! User supplied subroutine for evaluating the model.
      integer, intent(in) :: n
         !! Number of observations.
      integer, intent(in) :: m
         !! Number of columns of data in the explanatory variable.
      integer, intent(in) :: np
         !! Number of function parameters.
      integer, intent(in) :: q
         !! Number of responses per observation.
      real(wp), intent(inout) :: beta(np)
         !! Function parameters.
      real(wp), intent(in) :: x(n, m)
         !! Explanatory variable.
      real(wp), intent(in) :: delta(n, m)
         !! Estimated errors in the explanatory variables.
      real(wp), intent(inout) :: xplusd(n, m)
         !! Values of `x + delta`.
      integer, intent(in) :: ifixb(np)
         !! Values designating whether the elements of `beta` are fixed at their input values or not.
      integer, intent(in) :: ifixx(ldifx, m)
         !! Values designating whether the elements of `x` are fixed at their input values or not.
      integer, intent(in) :: ldifx
         !! Leading dimension of array `ifixx`.
      real(wp), intent(in) :: stpb(np)
         !! Relative step used for computing finite difference derivatives with respect to each `beta`.
      real(wp), intent(in) :: stpd(ldstpd, m)
         !! Relative step used for computing finite difference derivatives with respect to each `delta`.
      integer, intent(in) :: ldstpd
         !! Leading dimension of array `stpd`.
      real(wp), intent(in) :: ssf(np)
         !! Scaling values used for `beta`.
      real(wp), intent(in) :: tt(ldtt, m)
         !! Scaling values used for `delta`.
      integer, intent(in) :: ldtt
         !! Leading dimension of array `tt`.
      integer, intent(in) :: neta
         !! Number of good digits in the function results.
      real(wp), intent(in) :: fn(n, q)
         !! New predicted values from the function. Used when parameter is on a boundary.
      real(wp), intent(out) :: stp(n)
         !! Step used for computing finite difference derivatives with respect to each `delta`.
      real(wp), intent(out) :: wrk1(n, m, q)
         !! A work array of `(n, m, q)` elements.
      real(wp), intent(out) :: wrk2(n, q)
         !! A work array of `(n, q)` elements.
      real(wp), intent(out) :: wrk3(np)
         !! A work array of `(np)` elements.
      real(wp), intent(out) :: wrk6(n, np, q)
         !! A work array of `(n, np, q)` elements.
      real(wp), intent(out) :: fjacb(n, np, q)
         !! Jacobian with respect to `beta`.
      logical, intent(in) :: isodr
         !! Variable designating whether the solution is by ODR (`isodr = .true.`) or
         !! by OLS (`isodr = .false.`).
      real(wp), intent(out) :: fjacd(n, m, q)
         !! Jacobian with respect to `delta`.
      integer, intent(inout) :: nfev
         !! Number of function evaluations.
      integer, intent(out) :: istop
         !! Variable designating whether there are problems computing the function at the
         !! current `beta` and `delta`.
      integer, intent(out) :: info
         !! Variable designating why the computations were stopped.
      real(wp), intent(in) :: lower(np)
         !! Lower bound on `beta`.
      real(wp), intent(in) :: upper(np)
         !! Upper bound on `beta`.

      ! Local scalars
      real(wp) :: betak, typj
      integer :: i, j, k, l
      logical :: doit, setzero

      ! Variable Definitions (alphabetically)
      !  BETAK:    The K-th function parameter.
      !  DOIT:     The variable designating whether the derivative wrt a given BETA or DELTA
      !            needs to be computed (DOIT=TRUE) or not (DOIT=FALSE).
      !  I:        An indexing variable.
      !  J:        An indexing variable.
      !  K:        An indexing variable.
      !  L:        An indexing variable.
      !  SETZERO:  The variable designating whether the derivative wrt some DELTA needs to be
      !            set to zero (SETZRO=TRUE) or not (SETZRO=FALSE).
      !  TYPJ:     The typical size of the J-th unknown BETA or DELTA.

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
            fjacb(:, k, :) = zero
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
               wrk2 = fn
            else
               call fcn(beta, xplusd, ifixb, ifixx, 001, wrk2, wrk6, wrk1, istop)
               if (istop /= 0) then
                  return
               else
                  nfev = nfev + 1
               end if
            end if
            fjacb(:, k, :) = wrk2

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
               wrk2 = fn
            else
               call fcn(beta, xplusd, ifixb, ifixx, 001, wrk2, wrk6, wrk1, istop)
               if (istop /= 0) then
                  return
               else
                  nfev = nfev + 1
               end if
            end if

            fjacb(:, k, :) = (fjacb(:, k, :) - wrk2)/(2*wrk3(k))
            beta(k) = betak

         end if
      end do

      ! Compute the Jacobian wrt the X's
      if (isodr) then
         do j = 1, m
            if (ifixx(1, 1) < 0) then
               doit = .true.
               setzero = .false.
            elseif (ldifx == 1) then
               if (ifixx(1, j) == 0) then
                  doit = .false.
               else
                  doit = .true.
               end if
               setzero = .false.
            else
               doit = any(ifixx(:, j) /= 0)
               setzero = any(ifixx(:, j) == 0)
            end if

            if (.not. doit) then
               fjacd(:, j, :) = zero
            else
               do i = 1, n

                  if (xplusd(i, j) == zero) then
                     if (tt(1, 1) < zero) then
                        typj = 1/abs(tt(1, 1))
                     elseif (ldtt == 1) then
                        typj = 1/tt(1, j)
                     else
                        typj = 1/tt(i, j)
                     end if
                  else
                     typj = abs(xplusd(i, j))
                  end if

                  stp(i) = xplusd(i, j) &
                           + sign(one, xplusd(i, j))*typj*hstep(1, neta, i, j, stpd, ldstpd)
                  stp(i) = stp(i) - xplusd(i, j)
                  xplusd(i, j) = xplusd(i, j) + stp(i)

               end do

               istop = 0
               call fcn(beta, xplusd, ifixb, ifixx, 001, wrk2, wrk6, wrk1, istop)
               if (istop /= 0) then
                  return
               else
                  nfev = nfev + 1
                  fjacd(:, j, :) = wrk2
               end if

               xplusd(:, j) = x(:, j) + delta(:, j) - stp

               istop = 0
               call fcn(beta, xplusd, ifixb, ifixx, 001, wrk2, wrk6, wrk1, istop)
               if (istop /= 0) then
                  return
               else
                  nfev = nfev + 1
               end if

               if (setzero) then
                  do i = 1, n
                     if (ifixx(i, j) == 0) then
                        fjacd(i, j, :) = zero
                     else
                        fjacd(i, j, :) = (fjacd(i, j, :) - wrk2(i, :))/(2*stp(i))
                     end if
                  end do
               else
                  do l = 1, q
                     fjacd(:, j, l) = (fjacd(:, j, l) - wrk2(:, l))/(2*stp(:))
                  end do
               end if

               xplusd(:, j) = x(:, j) + delta(:, j)

            end if

         end do
      end if

   end subroutine jac_cdiff

   subroutine jac_fwdiff &
      (fcn, &
       n, m, np, q, &
       beta, x, delta, xplusd, ifixb, ifixx, ldifx, &
       stpb, stpd, ldstpd, &
       ssf, tt, ldtt, neta, fn, stp, wrk1, wrk2, wrk3, wrk6, &
       fjacb, isodr, fjacd, nfev, istop, info, &
       lower, upper)
   !! Compute forward difference approximations to the Jacobian wrt the estimated `beta`s and
   !! wrt the `delta`s.

      use odrpack_kinds, only: zero, one

      procedure(fcn_t) :: fcn
         !! User supplied subroutine for evaluating the model.
      integer, intent(in) :: n
         !! Number of observations.
      integer, intent(in) :: m
         !! Number of columns of data in the explanatory variable.
      integer, intent(in) :: np
         !! Number of function parameters.
      integer, intent(in) :: q
         !! Number of responses per observation.
      real(wp), intent(inout) :: beta(np)
         !! Function parameters.
      real(wp), intent(in) :: x(n, m)
         !! Explanatory variable.
      real(wp), intent(in) :: delta(n, m)
         !! Estimated errors in the explanatory variables.
      real(wp), intent(inout) :: xplusd(n, m)
         !! Values of `x + delta`.
      integer, intent(in) :: ifixb(np)
         !! Values designating whether the elements of `beta` are fixed at their input values or not.
      integer, intent(in) :: ifixx(ldifx, m)
         !! Values designating whether the elements of `x` are fixed at their input values or not.
      integer, intent(in) :: ldifx
         !! Leading dimension of array `ifixx`.
      real(wp), intent(in) :: stpb(np)
         !! Relative step used for computing finite difference derivatives with respect to each `beta`.
      real(wp), intent(in) :: stpd(ldstpd, m)
         !! Relative step used for computing finite difference derivatives with respect to each `delta`.
      integer, intent(in) :: ldstpd
         !! Leading dimension of array `stpd`.
      real(wp), intent(in) :: ssf(np)
         !! Scaling values used for `beta`.
      real(wp), intent(in) :: tt(ldtt, m)
         !! Scaling values used for `delta`.
      integer, intent(in) :: ldtt
         !! Leading dimension of array `tt`.
      integer, intent(in) :: neta
         !! Number of good digits in the function results.
      real(wp), intent(in) :: fn(n, q)
         !! New predicted values from the function. Used when parameter is on a boundary.
      real(wp), intent(out) :: stp(n)
         !! Step used for computing finite difference derivatives with respect to each `delta`.
      real(wp), intent(out) :: wrk1(n, m, q)
         !! A work array of `(n, m, q)` elements.
      real(wp), intent(out) :: wrk2(n, q)
         !! A work array of `(n, q)` elements.
      real(wp), intent(out) :: wrk3(np)
         !! A work array of `(np)` elements.
      real(wp), intent(out) :: wrk6(n, np, q)
         !! A work array of `(n, np, q)` elements.
      real(wp), intent(out) :: fjacb(n, np, q)
         !! Jacobian with respect to `beta`.
      logical, intent(in) :: isodr
         !! Variable designating whether the solution is by ODR (`isodr = .true.`) or
         !! by OLS (`isodr = .false.`).
      real(wp), intent(out) :: fjacd(n, m, q)
         !! Jacobian with respect to `delta`.
      integer, intent(inout) :: nfev
         !! Number of function evaluations.
      integer, intent(out) :: istop
         !! Variable designating whether there are problems computing the function at the
         !! current `beta` and `delta`.
      integer, intent(out) :: info
         !! Variable designating why the computations were stopped.
      real(wp), intent(in) :: lower(np)
         !! Lower bound on `beta`.
      real(wp), intent(in) :: upper(np)
         !! Upper bound on `beta`.

      ! Local scalars
      real(wp) :: betak, step, typj
      integer :: i, j, k, l
      logical :: doit, setzero

      ! Variable Definitions (alphabetically)
      !  BETAK:   The K-th function parameter.
      !  DOIT:    The variable designating whether the derivative wrt a given BETA or DELTA needs
      !           to be computed (DOIT=TRUE) or not (DOIT=FALSE).
      !  I:       An indexing variable.
      !  J:       An indexing variable.
      !  K:       An indexing variable.
      !  L:       An indexing variable.
      !  SETZRO:  The variable designating whether the derivative wrt some DELTA needs to be set
      !           to zero (SETZRO=TRUE) or not (SETZRO=FALSE).
      !  TYPJ:    The typical size of the J-th unknown BETA or DELTA.

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
            fjacb(:, k, :) = zero
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
            call fcn(beta, xplusd, ifixb, ifixx, 001, wrk2, wrk6, wrk1, istop)
            if (istop /= 0) then
               return
            else
               nfev = nfev + 1
            end if
            fjacb(:, k, :) = (wrk2 - fn)/wrk3(k)
            beta(k) = betak
         end if
      end do

      ! Compute the Jacobian wrt the X'S
      if (isodr) then
         do j = 1, m

            if (ifixx(1, 1) < 0) then
               doit = .true.
               setzero = .false.
            elseif (ldifx == 1) then
               if (ifixx(1, j) == 0) then
                  doit = .false.
               else
                  doit = .true.
               end if
               setzero = .false.
            else
               doit = any(ifixx(:, j) /= 0)
               setzero = any(ifixx(:, j) == 0)
            end if

            if (.not. doit) then
               fjacd(:, j, :) = zero
            else
               do i = 1, n

                  if (xplusd(i, j) == zero) then
                     if (tt(1, 1) < zero) then
                        typj = 1/abs(tt(1, 1))
                     elseif (ldtt == 1) then
                        typj = 1/tt(1, j)
                     else
                        typj = 1/tt(i, j)
                     end if
                  else
                     typj = abs(xplusd(i, j))
                  end if

                  stp(i) = xplusd(i, j) &
                           + sign(one, xplusd(i, j))*typj*hstep(0, neta, i, j, stpd, ldstpd)
                  stp(i) = stp(i) - xplusd(i, j)
                  xplusd(i, j) = xplusd(i, j) + stp(i)

               end do

               istop = 0
               call fcn(beta, xplusd, ifixb, ifixx, 001, wrk2, wrk6, wrk1, istop)
               if (istop /= 0) then
                  return
               else
                  nfev = nfev + 1
                  fjacd(:, j, :) = wrk2
               end if

               if (setzero) then
                  do i = 1, n
                     if (ifixx(i, j) == 0) then
                        fjacd(i, j, :) = zero
                     else
                        fjacd(i, j, :) = (fjacd(i, j, :) - fn(i, :))/stp(i)
                     end if
                  end do
               else
                  do l = 1, q
                     fjacd(:, j, l) = (fjacd(:, j, l) - fn(:, l))/stp(:)
                  end do
               end if

               xplusd(:, j) = x(:, j) + delta(:, j)

            end if

         end do
      end if

   end subroutine jac_fwdiff

   subroutine check_jac &
      (fcn, &
       n, m, np, q, &
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

      use odrpack_kinds, only: zero, one, p5 => half

      procedure(fcn_t) :: fcn
         !! User supplied subroutine for evaluating the model.
      integer, intent(in) :: n
         !! Number of observations.
      integer, intent(in) :: m
         !! Number of columns of data in the explanatory variable.
      integer, intent(in) :: np
         !! Number of function parameters.
      integer, intent(in) :: q
         !! Number of responses per observation.
      real(wp), intent(inout) :: beta(np)
         !! Function parameters.
      real(wp), intent(inout) :: betaj(np)
         !! Function parameters offset such that steps don't cross bounds.
      real(wp), intent(inout) :: xplusd(n, m)
         !! Values of `x + delta`.
      integer, intent(in) :: ifixb(np)
         !! Values designating whether the elements of `beta` are fixed at their input values or not.
      integer, intent(in) :: ifixx(ldifx, m)
         !! Values designating whether the elements of `x` are fixed at their input values or not.
      integer, intent(in) :: ldifx
         !! Leading dimension of array `ifixx`.
      real(wp), intent(in) :: stpb(np)
         !! Step size for finite difference derivatives wrt `beta`.
      real(wp), intent(in) :: stpd(ldstpd, m)
         !! Step size for finite difference derivatives wrt `delta`.
      integer, intent(in) :: ldstpd
         !! Leading dimension of array `stpd`.
      real(wp), intent(in) :: ssf(np)
         !! Scaling values used for `beta`.
      real(wp), intent(in) :: tt(ldtt, m)
         !! Scaling values used for `delta`.
      integer, intent(in) :: ldtt
         !! Leading dimension of array `tt`.
      real(wp), intent(in) :: eta
         !! Relative noise in the function results.
      integer, intent(in) :: neta
         !! Number of reliable digits in the model results.
      integer, intent(out) :: ntol
         !! Number of digits of agreement required between the numerical derivatives and the
         !! user supplied derivatives.
      integer, intent(in) :: nrow
         !! Row number of the explanatory variable array at which the derivative is checked.
      logical, intent(in) :: isodr
         !! Variable designating whether the solution is by ODR (`isodr = .true.`) or
         !! by OLS (`isodr = .false.`).
      real(wp), intent(in) :: epsmac
         !! Value of machine precision.
      real(wp), intent(in) :: pv0i(n, q)
         !! Predicted values using the user supplied parameter estimates.
      real(wp), intent(out) :: fjacb(n, np, q)
         !! Jacobian with respect to `beta`.
      real(wp), intent(out) :: fjacd(n, m, q)
         !! Jacobian with respect to `delta`.
      integer, intent(out) :: msgb(1 + q*np)
         !! Error checking results for the Jacobian wrt `beta`.
      integer, intent(out) :: msgd(1 + q*m)
         !! Error checking results for the Jacobian wrt `delta`.
      real(wp), intent(out) :: diff(q, np + m)
         !! Relative differences between the user supplied and finite difference derivatives
         !! for each derivative checked.
      integer, intent(out) :: istop
         !! Variable designating whether there are problems computing the function at the
         !! current `beta` and `delta`
      integer, intent(inout) :: nfev
         !! Number of function evaluations.
      integer, intent(inout) :: njev
         !! Number of Jacobian evaluations.
      real(wp), intent(out) :: wrk1(n, m, q)
         !! A work array of `(n, m, q)` elements.
      real(wp), intent(out) :: wrk2(n, q)
         !! A work array of `(n, q)` elements.
      real(wp), intent(out) :: wrk6(n, np, q)
         !! A work array of `(n, np, q)` elements.
      integer, intent(in) :: interval(np)
         !! Specifies which checks can be performed when checking derivatives based on the
         !! interval of the bound constraints.

      ! Local scalars
      real(wp) :: diffj, h0, hc0, pv, tol, typj
      integer :: ideval, j, lq, msgb1, msgd1
      logical :: isfixd, iswrtb

      ! Local arrays
      real(wp) :: pv0(n, q)

      ! Variable Definitions (alphabetically)
      !  DIFFJ:    The relative differences between the user supplied and finite difference
      !            derivatives for the derivative being checked.
      !  H0:       The initial relative step size for forward differences.
      !  HC0:      The initial relative step size for central differences.
      !  IDEVAL:   The variable designating what computations are to be performed by user supplied
      !            subroutine FCN.
      !  ISFIXD:   The variable designating whether the parameter is fixed (ISFIXD=TRUE) or not (ISFIXD=FALSE).
      !  ISWRTB:   The variable designating whether the derivatives wrt BETA (ISWRTB=TRUE) or DELTA
      !            (ISWRTB=FALSE) are being checked.
      !  J:        An index variable.
      !  LQ:       The response currently being examined.
      !  MSGB1:    The error checking results for the Jacobian wrt BETA.
      !  MSGD1:    The error checking results for the Jacobian wrt DELTA.
      !  PV:       The scalar in which the predicted value from the model for row NROW is stored.
      !  PV0:      The predicted values using the current parameter estimates (possibly offset from
      !            the user supplied estimates to create distance between parameters and the bounds
      !            on the parameters).
      !  TOL:      The agreement tolerance.
      !  TYPJ:     The typical size of the J-th unknown BETA or DELTA.

      ! Set tolerance for checking derivatives
      tol = eta**(0.25E0_wp)
      ntol = int(max(one, p5 - log10(tol)))

      ! Compute, if necessary, PV0
      pv0 = pv0i
      if (any(beta /= betaj)) then
         istop = 0
         ideval = 001
         call fcn(betaj, xplusd, ifixb, ifixx, ideval, pv0, fjacb, fjacd, istop)
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
      call fcn(betaj, xplusd, ifixb, ifixx, ideval, wrk2, fjacb, fjacd, istop)
      if (istop /= 0) then
         return
      else
         njev = njev + 1
      end if

      ! Check derivatives wrt BETA for each response of observation NROW
      msgb1 = 0
      msgd1 = 0

      do lq = 1, q

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
               msgb(1 + lq + (j - 1)*q) = -1
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

               h0 = hstep(0, neta, 1, j, stpb, 1)
               hc0 = h0

               ! Check derivative wrt the J-th parameter at the NROW-th row
               if (interval(j) >= 1) then
                  call check_jac_value(fcn, &
                                       n, m, np, q, &
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
                  msgd(1 + lq + (j - 1)*q) = -1
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

                  h0 = hstep(0, neta, nrow, j, stpd, ldstpd)
                  hc0 = hstep(1, neta, nrow, j, stpd, ldstpd)

                  ! Check derivative wrt the J-th column of DELTA at row NROW
                  call check_jac_value(fcn, &
                                       n, m, np, q, &
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

   end subroutine check_jac

   subroutine check_jac_curv &
      (fcn, &
       n, m, np, q, &
       beta, xplusd, ifixb, ifixx, ldifx, &
       eta, tol, nrow, epsmac, j, lq, hc, iswrtb, &
       fd, typj, pvpstp, stp0, &
       pv, d, &
       diffj, msg, istop, nfev, &
       wrk1, wrk2, wrk6)
   !! Check whether high curvature could be the cause of the disagreement between the numerical
   !! and analytic derviatives.
      ! Adapted from STARPAC subroutine DCKCRV.

      use odrpack_kinds, only: one, two, ten

      procedure(fcn_t) :: fcn
         !! User supplied subroutine for evaluating the model.
      integer, intent(in) :: n
         !! Number of observations.
      integer, intent(in) :: m
         !! Number of columns of data in the explanatory variable.
      integer, intent(in) :: np
         !! Number of function parameters.
      integer, intent(in) :: q
         !! Number of responses per observation.
      real(wp), intent(inout) :: beta(np)
         !! Function parameters.
      real(wp), intent(inout) :: xplusd(n, m)
         !! Values of `x` + `delta`.
      integer, intent(in) :: ifixb(np)
         !! Values designating whether the elements of `beta` are fixed at their input values or not.
      integer, intent(in) :: ifixx(ldifx, m)
         !! Values designating whether the elements of `x` are fixed at their input values or not.
      integer, intent(in) :: ldifx
         !! Leading dimension of array `ifixx`.
      real(wp), intent(in) :: eta
         !! Relative noise in the model.
      real(wp), intent(in) :: tol
         !! Agreement tolerance.
      integer, intent(in) :: nrow
         !! Row number of the explanatory variable array at which the derivative is to be checked.
      real(wp), intent(in) :: epsmac
         !! Value of machine precision.
      integer, intent(in) :: j
         !! Index of the partial derivative being examined.
      integer, intent(in) :: lq
         !! Response currently being examined.
      real(wp), intent(in) :: hc
         !! Relative step size for central finite differences.
      logical, intent(in) :: iswrtb
         !! Variable designating whether the derivatives wrt `beta` (`iswrtb = .true.`) or
         !! `delta` (`iswrtb = .false.`) are being checked.
      real(wp), intent(out) :: fd
         !! Forward difference derivative wrt the `j`-th parameter.
      real(wp), intent(in) :: typj
         !! Typical size of the `j`-th unknown `beta` or `delta`.
      real(wp), intent(out) :: pvpstp
         !! Predicted value for row `nrow` of the model based on the current parameter estimates
         !! for all but the `j`-th parameter value, which is `beta(j) + stp0`.
      real(wp), intent(in) :: stp0
         !! Initial step size for the finite difference derivative.
      real(wp), intent(in) :: pv
         !! Predicted value of the model for row `nrow`.
      real(wp), intent(in) :: d
         !! Derivative with respect to the `j`-th unknown parameter.
      real(wp), intent(out) :: diffj
         !! Relative differences between the user supplied and finite difference derivatives
         !! for the derivative being checked.
      integer, intent(out) :: msg(q, j)
         !! Error checking results.
      integer, intent(out) :: istop
         !! Variable designating whether there are problems computing the function at the
         !! current `beta` and `delta`.
      integer, intent(inout) :: nfev
         !! Number of function evaluations.
      real(wp), intent(out) :: wrk1(n, m, q)
         !! A work array of `(n, m, q)` elements.
      real(wp), intent(out) :: wrk2(n, q)
         !! A work array of `(n, q)` elements.
      real(wp), intent(out) :: wrk6(n, np, q)
         !! A work array of `(n, np, q)` elements.

      ! Local scalars
      real(wp), parameter :: p01 = 0.01_wp
      real(wp) :: curve, pvmcrv, pvpcrv, stp, stpcrv

      ! Variable Definitions (alphabetically)
      !  CURVE:   A measure of the curvature in the model.
      !  PVMCRV:  The predicted value for row  NROW of the model based on the current parameter
      !           estimates for all but the Jth parameter value, which is BETA(J)-STPCRV.
      !  PVPCRV:  The predicted value for row NROW of the model based on the current parameter
      !           estimates for all but the Jth parameter value, which is BETA(J)+STPCRV.
      !  STP:     A step size for the finite difference derivative.
      !  STPCRV:  The step size selected to check for curvature in the model.

      if (iswrtb) then

         ! Perform central difference computations for derivatives wrt BETA
         stpcrv = (hc*typj*sign(one, beta(j)) + beta(j)) - beta(j)
         call fpvb(fcn, &
                   n, m, np, q, &
                   beta, xplusd, ifixb, ifixx, ldifx, &
                   nrow, j, lq, stpcrv, &
                   istop, nfev, pvpcrv, &
                   wrk1, wrk2, wrk6)
         if (istop /= 0) then
            return
         end if
         call fpvb(fcn, &
                   n, m, np, q, &
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
         call fpvd(fcn, &
                   n, m, np, q, &
                   beta, xplusd, ifixb, ifixx, ldifx, &
                   nrow, j, lq, stpcrv, &
                   istop, nfev, pvpcrv, &
                   wrk1, wrk2, wrk6)
         if (istop /= 0) then
            return
         end if
         call fpvd(fcn, &
                   n, m, np, q, &
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
      call check_jac_fp(fcn, &
                        n, m, np, q, &
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
         call fpvb(fcn, &
                   n, m, np, q, &
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
         call fpvd(fcn, &
                   n, m, np, q, &
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

   end subroutine check_jac_curv

   subroutine check_jac_fp &
      (fcn, &
       n, m, np, q, &
       beta, xplusd, ifixb, ifixx, ldifx, &
       eta, tol, nrow, j, lq, iswrtb, &
       fd, typj, pvpstp, stp0, curve, pv, d, &
       diffj, msg, istop, nfev, &
       wrk1, wrk2, wrk6)
   !! Check whether finite precision arithmetic could be the cause of the disagreement between
   !! the derivatives.
      ! Adapted from STARPAC subroutine DCKFPA.

      use odrpack_kinds, only: one, two, hundred

      procedure(fcn_t) :: fcn
         !! User supplied subroutine for evaluating the model.
      integer, intent(in) :: n
         !! Number of observations.
      integer, intent(in) :: m
         !! Number of columns of data in the explanatory variable.
      integer, intent(in) :: np
         !! Number of function parameters.
      integer, intent(in) :: q
         !! Number of responses per observation.
      real(wp), intent(inout) :: beta(np)
         !! Function parameters.
      real(wp), intent(inout) :: xplusd(n, m)
         !! Values of `x + delta`.
      integer, intent(in) :: ifixb(np)
         !! Values designating whether the elements of `beta` are fixed at their input values or not.
      integer, intent(in) :: ifixx(ldifx, m)
         !! Values designating whether the elements of `x` are fixed at their input values or not.
      integer, intent(in) :: ldifx
         !! Leading dimension of array `ifixx`.
      real(wp), intent(in) :: eta
         !! Relative noise in the model.
      real(wp), intent(in) :: tol
         !! Agreement tolerance.
      integer, intent(in) :: nrow
         !! Row number of the explanatory variable array at which the derivative is to be checked.
      integer, intent(in) :: j
         !! Index of the partial derivative being examined.
      integer, intent(in) :: lq
         !! Response currently being examined.
      logical, intent(in) :: iswrtb
         !! Variable designating whether the derivatives wrt `beta` (`iswrtb = .true.`)
         !! or `delta` (`iswrtb = .false.`) are being checked.
      real(wp), intent(out) :: fd
         !! Forward difference derivative wrt the `j`-th parameter.
      real(wp), intent(in) :: typj
         !! Typical size of the `j`-th unknown `beta` or `delta`.
      real(wp), intent(out) :: pvpstp
         !! Predicted value for row `nrow` of the model based on the current parameter
         !! estimates for all but the `j`-th parameter value, which is `beta(j) + stp0`.
      real(wp), intent(in) :: stp0
         !! Step size for the finite difference derivative.
      real(wp), intent(inout) :: curve
         !! A measure of the curvature in the model.
      real(wp), intent(in) :: pv
         !! Predicted value for row `nrow`.
      real(wp), intent(in) :: d
         !! Derivative with respect to the `j`-th unknown parameter.
      real(wp), intent(out) :: diffj
         !! Relative differences between the user supplied and finite difference derivatives
         !! for the derivative being checked.
      integer, intent(out) :: msg(q, j)
         !! Error checking results.
      integer, intent(out) :: istop
         !! Variable designating whether there are problems computing the function at the
         !! current `beta` and `delta`.
      integer, intent(inout) :: nfev
         !! Number of function evaluations.
      real(wp), intent(out) :: wrk1(n, m, q)
         !! A work array of `(n, m, q)` elements.
      real(wp), intent(out) :: wrk2(n, q)
         !! A work array of `(n, q)` elements.
      real(wp), intent(out) :: wrk6(n, np, q)
         !! A work array of `(n, np, q)` elements.

      ! Local scalars
      real(wp), parameter :: p1 = 0.1_wp
      real(wp) :: stp
      logical :: large

      ! Variable Definitions (alphabetically)
      !  LARGE:   The value designating whether the recommended increase in the step size would
      !           be greater than TYPJ.
      !  STP:     A step size for the finite difference derivative.

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
         call fpvb(fcn, &
                   n, m, np, q, &
                   beta, xplusd, ifixb, ifixx, ldifx, &
                   nrow, j, lq, stp, &
                   istop, nfev, pvpstp, &
                   wrk1, wrk2, wrk6)
      else
         ! Perform computations for derivatives wrt DELTA
         stp = (stp*sign(one, xplusd(nrow, j)) + xplusd(nrow, j)) - xplusd(nrow, j)
         call fpvd(fcn, &
                   n, m, np, q, &
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

   end subroutine check_jac_fp

   subroutine check_jac_value &
      (fcn, &
       n, m, np, q, &
       beta, xplusd, ifixb, ifixx, ldifx, &
       eta, tol, nrow, epsmac, j, lq, typj, h0, hc0, &
       iswrtb, pv, d, &
       diffj, msg1, msg, istop, nfev, &
       wrk1, wrk2, wrk6, interval)
   !! Check user supplied analytic derivatives against numerical derivatives.
      ! Adapted from STARPAC subroutine DCKMN.

      use odrpack_kinds, only: zero, one, two, three, ten, hundred

      procedure(fcn_t) :: fcn
         !! User supplied subroutine for evaluating the model.
      integer, intent(in) :: n
         !! Number of observations.
      integer, intent(in) :: m
         !! Number of columns of data in the explanatory variable.
      integer, intent(in) :: np
         !! Number of function parameters.
      integer, intent(in) :: q
         !! Number of responses per observation.
      real(wp), intent(inout) :: beta(np)
         !! Function parameters.
      real(wp), intent(inout) :: xplusd(n, m)
         !! Values of `x + delta`.
      integer, intent(in) :: ifixb(np)
         !! Values designating whether the elements of `beta` are fixed at their input values or not.
      integer, intent(in) :: ifixx(ldifx, m)
         !! Values designating whether the elements of `x` are fixed at their input values or not.
      integer, intent(in) :: ldifx
         !! Leading dimension of array `ifixx`.
      real(wp), intent(in) :: eta
         !! Relative noise in the model.
      real(wp), intent(in) :: tol
         !! Agreement tolerance.
      integer, intent(in) :: nrow
         !! Row number of the explanatory variable array at which the derivative is to be checked.
      real(wp), intent(in) :: epsmac
         !! Value of machine precision.
      integer, intent(in) :: j
         !! Index of the partial derivative being examined.
      integer, intent(in) :: lq
         !! Response currently being examined.
      real(wp), intent(in) :: typj
         !! Typical size of the `j`-th unknown `beta` or `delta`.
      real(wp), intent(in) :: h0
         !! Initial step size for the finite difference derivative.
      real(wp), intent(in) :: hc0
         !! Relative step size for central finite differences.
      logical, intent(in) :: iswrtb
         !! Variable designating whether the derivatives wrt `beta` (`iswrtb = .true.`)
         !! or `delta` (`iswrtb = .false.`) are being checked.
      real(wp), intent(in) :: pv
         !! Predicted value for row `nrow`.
      real(wp), intent(in) :: d
         !! Derivative with respect to the `j`-th unknown parameter.
      real(wp), intent(out) :: diffj
         !! Relative differences between the user supplied and finite difference derivatives
         !! for the derivative being checked.
      integer, intent(out) :: msg1
         !! First set of error checking results.
      integer, intent(out) :: msg(q, j)
         !! Error checking results.
      integer, intent(out) :: istop
         !! Variable designating whether there are problems computing the function at the
         !! current `beta` and `delta`.
      integer, intent(inout) :: nfev
         !! Number of function evaluations.
      real(wp), intent(out) :: wrk1(n, m, q)
         !! A work array of `(n, m, q)` elements.
      real(wp), intent(out) :: wrk2(n, q)
         !! A work array of `(n, q)` elements.
      real(wp), intent(out) :: wrk6(n, np, q)
         !! A work array of `(n, np, q)` elements.
      integer, intent(in) :: interval(np)
         !! Specifies which checks can be performed when checking derivatives based on the
         !! interval of the bound constraints.

      ! Local scalars
      real(wp), parameter :: p01 = 0.01_wp, p1 = 0.1_wp
      real(wp), parameter :: big = 1.0E19_wp, tol2 = 5.0E-2_wp
      real(wp) :: fd, h, hc, h1, hc1, pvpstp, stp0
      integer :: i

      ! Variable Definitions (alphabetically)
      !  BIG:      A big value, used to initialize DIFFJ.
      !  FD:       The forward difference derivative wrt the Jth parameter.
      !  H:        The relative step size for forward differences.
      !  H1:       The default relative step size for forward differences.
      !  HC:       The relative step size for central differences.
      !  HC1:      The default relative step size for central differences.
      !  PVPSTP:   The predicted value for row NROW of the model using the current parameter
      !            estimates for all but the Jth parameter value, which is BETA(J) + STP0.
      !  STP0:     The initial step size for the finite difference derivative.
      !  TOL2:     A minimum agreement tolerance.

      ! Calculate the Jth partial derivative using forward difference quotients and decide if
      ! it agrees with user supplied values

      h1 = sqrt(eta)
      hc1 = eta**(1/three)

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
            call fpvb(fcn, &
                      n, m, np, q, &
                      beta, xplusd, ifixb, ifixx, ldifx, &
                      nrow, j, lq, stp0, &
                      istop, nfev, pvpstp, &
                      wrk1, wrk2, wrk6)
         else
            ! Perform computations for derivatives wrt DELTA
            stp0 = (h*typj*sign(one, xplusd(nrow, j)) + xplusd(nrow, j)) - xplusd(nrow, j)
            call fpvd(fcn, &
                      n, m, np, q, &
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
                  call check_jac_zero(fcn, &
                                      n, m, np, q, &
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
                  call check_jac_curv(fcn, &
                                      n, m, np, q, &
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

   end subroutine check_jac_value

   subroutine check_jac_zero &
      (fcn, &
       n, m, np, q, &
       beta, xplusd, ifixb, ifixx, ldifx, &
       nrow, epsmac, j, lq, iswrtb, &
       tol, d, fd, typj, pvpstp, stp0, pv, &
       diffj, msg, istop, nfev, &
       wrk1, wrk2, wrk6)
   !! Recheck the derivatives in the case where the finite difference derivative disagrees with
   !! the analytic derivative and the analytic derivative is zero.
      ! Adapted from STARPAC subroutine DCKZRO.

      use odrpack_kinds, only: zero, three

      procedure(fcn_t) :: fcn
         !! User supplied subroutine for evaluating the model.
      integer, intent(in) :: n
         !! Number of observations.
      integer, intent(in) :: m
         !! Number of columns of data in the explanatory variable.
      integer, intent(in) :: np
         !! Number of function parameters.
      integer, intent(in) :: q
         !! Number of responses per observation.
      real(wp), intent(inout) :: beta(np)
         !! Function parameters.
      real(wp), intent(inout) :: xplusd(n, m)
         !! Values of `x + delta`.
      integer, intent(in) :: ifixb(np)
         !! Values designating whether the elements of `beta` are fixed at their input values or not.
      integer, intent(in) :: ifixx(ldifx, m)
         !! Values designating whether the elements of `x` are fixed at their input values or not.
      integer, intent(in) :: ldifx
         !! Leading dimension of array `ifixx`.
      integer, intent(in) :: nrow
         !! Row number of the explanatory variable array at which the derivative is to be checked.
      real(wp), intent(in) :: epsmac
         !! Value of machine precision.
      integer, intent(in) :: j
         !! Index of the partial derivative being examined.
      integer, intent(in) :: lq
         !! Response currently being examined.
      logical, intent(in) :: iswrtb
         !! Variable designating whether the derivatives wrt `beta` (`iswrtb = .true.`)
         !! or `delta` (`iswrtb = .false.`) are being checked.
      real(wp), intent(in) :: tol
         !! Agreement tolerance.
      real(wp), intent(in) :: d
         !! Derivative with respect to the `j`-th unknown parameter.
      real(wp), intent(in) :: fd
         !! Forward difference derivative wrt the `j`-th parameter.
      real(wp), intent(in) :: typj
         !! Typical size of the `j`-th unknown `beta` or `delta`.
      real(wp), intent(in) :: pvpstp
         !! Predicted value for row `nrow` of the model using the current parameter estimates
         !! for all but the `j`-th parameter value, which is `beta(j) + stp0`.
      real(wp), intent(in) :: stp0
         !! Initial step size for the finite difference derivative.
      real(wp), intent(in) :: pv
         !! Predicted value from the model for row `nrow`.
      real(wp), intent(out) :: diffj
         !! Relative differences between the user supplied and finite difference derivatives
         !! for the derivative being checked.
      integer, intent(out) :: msg(q, j)
         !! Error checking results.
      integer, intent(out) :: istop
         !! Variable designating whether there are problems computing the function at the
         !! current `beta` and `delta`.
      integer, intent(inout) :: nfev
         !! Number of function evaluations.
      real(wp), intent(out) :: wrk1(n, m, q)
         !! A work array of `(n, m, q)` elements.
      real(wp), intent(out) :: wrk2(n, q)
         !! A work array of `(n, q)` elements.
      real(wp), intent(out) :: wrk6(n, np, q)
         !! A work array of `(n, np, q)` elements.

      ! Local scalars
      real(wp) :: cd, pvmstp

      ! Variable Definitions (alphabetically)
      !  CD:      The central difference derivative wrt the Jth parameter.
      !  PVMSTP:  The predicted value for row NROW of the model using the current parameter
      !           estimates for all but the Jth parameter value, which is BETA(J) - STP0.

      ! Recalculate numerical derivative using central difference and step size of 2*STP0
      if (iswrtb) then
         ! Perform computations for derivatives wrt BETA
         call fpvb(fcn, &
                   n, m, np, q, &
                   beta, xplusd, ifixb, ifixx, ldifx, &
                   nrow, j, lq, -stp0, &
                   istop, nfev, pvmstp, &
                   wrk1, wrk2, wrk6)
      else
         ! Perform computations for derivatives wrt DELTA
         call fpvd(fcn, &
                   n, m, np, q, &
                   beta, xplusd, ifixb, ifixx, ldifx, &
                   nrow, j, lq, -stp0, &
                   istop, nfev, pvmstp, &
                   wrk1, wrk2, wrk6)
      end if

      if (istop /= 0) then
         return
      end if

      cd = (pvpstp - pvmstp)/(2*stp0)
      diffj = min(abs(cd - d), abs(fd - d))

      ! Check for agreement
      if (diffj <= tol*abs(d)) then
         ! Finite difference and analytic derivatives now agree
         if (d == zero) then
            msg(lq, j) = 1
         else
            msg(lq, j) = 0
         end if
      elseif (diffj*typj <= abs(pv*epsmac**(1/three))) then
         ! Derivatives are both close to zero
         msg(lq, j) = 2
      else
         ! Derivatives are not both close to zero
         msg(lq, j) = 3
      end if

   end subroutine check_jac_zero

   pure subroutine check_inputs &
      (n, m, np, q, &
       isodr, anajac, &
       beta, ifixb, &
       ldifx, ldscld, ldstpd, ldwe, ld2we, ldwd, ld2wd, &
       lrwork, lrwkmn, liwork, liwkmn, &
       sclb, scld, stpb, stpd, &
       info, &
       lower, upper)
   !! Check input parameters, indicating errors found using nonzero values of argument `info`.

      use odrpack_kinds, only: zero

      integer, intent(in) :: n
         !! Number of observations.
      integer, intent(in) :: m
         !! Number of columns of data in the explanatory variable.
      integer, intent(in) :: np
         !! Number of function parameters.
      integer, intent(in) :: q
         !! Number of responses per observation.
      logical, intent(in) :: isodr
         !! Variable designating whether the solution is by ODR (`isodr = .true.`) or
         !! by OLS (`isodr = .false.`).
      logical, intent(in) :: anajac
         !! Variable designating whether the Jacobians are computed by finite differences
         !! (`anajac = .false.`) or not (`anajac = .true.`).
      real(wp), intent(in) :: beta(np)
         !! Function parameters.
      integer, intent(in) :: ifixb(np)
         !! Values designating whether the elements of `beta` are fixed at their input values or not.
      integer, intent(in) :: ldifx
         !! Leading dimension of array `ifixx`.
      integer, intent(in) :: ldscld
         !! Leading dimension of array `scld`.
      integer, intent(in) :: ldstpd
         !! Leading dimension of array `stpd`.
      integer, intent(in) :: ldwe
      !! Leading dimension of array `we`.
      integer, intent(in) :: ld2we
         !! Second dimension of array `we`.
      integer, intent(in) :: ldwd
         !! Leading dimension of array `wd`.
      integer, intent(in) :: ld2wd
         !! Second dimension of array `wd`.
      integer, intent(in) :: lrwork
         !! Length of vector `rwork`.
      integer, intent(in) :: lrwkmn
         !! Minimum acceptable length of array `rwork`.
      integer, intent(in) :: liwork
         !! Length of vector `iwork`.
      integer, intent(in) :: liwkmn
         !! Minimum acceptable length of array `iwork`.
      real(wp), intent(in) :: sclb(np)
         !! Scaling values for `beta`.
      real(wp), intent(in) :: scld(ldscld, m)
         !! Scaling value for `delta`.
      real(wp), intent(in) :: stpb(np)
         !! Step for the finite difference derivative wrt `beta`.
      real(wp), intent(in) :: stpd(ldstpd, m)
         !! Step for the finite difference derivative wrt `delta`.
      integer, intent(out) :: info
         !! Variable designating why the computations were stopped.
      real(wp), intent(in) :: lower(np)
         !! Lower bound on `beta`.
      real(wp), intent(in) :: upper(np)
         !! Upper bound on `beta`.

      ! Local scalars
      integer :: last, npp

      ! Variable Definitions (alphabetically)
      !  LAST:    The last row of the array to be accessed.
      !  NPP:     The number of function parameters being estimated.

      ! Find actual number of parameters being estimated
      if ((np <= 0) .or. (ifixb(1) < 0)) then
         npp = np
      else
         npp = count(ifixb /= 0)
      end if

      ! Check problem specification parameters
      if ((n <= 0) .or. (m <= 0) .or. (npp <= 0 .or. npp > n) .or. (q <= 0)) then
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
         if (q <= 0) then
            info = info + 1
         end if
         return
      end if

      ! Check dimension specification parameters
      if (((ldwe /= 1) .and. (ldwe < n)) .or. &
          ((ld2we /= 1) .and. (ld2we < q)) .or. &
          (isodr .and. ((ldwd /= 1) .and. (ldwd < n))) .or. &
          (isodr .and. ((ld2wd /= 1) .and. (ld2wd < m))) .or. &
          (isodr .and. ((ldifx /= 1) .and. (ldifx < n))) .or. &
          (isodr .and. ((ldstpd /= 1) .and. (ldstpd < n))) .or. &
          (isodr .and. ((ldscld /= 1) .and. (ldscld < n))) .or. &
          (lrwork < lrwkmn) .or. &
          (liwork < liwkmn)) then

         info = 20000

         if ((ldwe /= 1 .and. ldwe < n) .or. (ld2we /= 1 .and. ld2we < q)) then
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

         if (lrwork < lrwkmn) then
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
         if (any(scld(1:last, :) <= zero)) then
            info = 30200
         end if
      end if

      ! Check BETA scaling
      if (sclb(1) > zero) then
         if (any(sclb <= zero)) then
            if (info == 0) then
               info = 30100
            else
               info = info + 100
            end if
         end if
      end if

      ! Check DELTA finite difference step sizes
      if ((.not. anajac) .and. isodr .and. stpd(1, 1) > zero) then
         if (ldstpd >= n) then
            last = n
         else
            last = 1
         end if
         if (any(stpd(1:last, :) <= zero)) then
            if (info == 0) then
               info = 32000
            else
               info = info + 2000
            end if
         end if
      end if

      ! Check BETA finite difference step sizes
      if ((.not. anajac) .and. stpb(1) > zero) then
         if (any(stpb <= zero)) then
            if (info == 0) then
               info = 31000
            else
               info = info + 1000
            end if
         end if
      end if

      !  Check bounds
      if (any(upper < lower)) then
         if (info == 0) then
            info = 91000
         end if
      end if

      if (any( &
          ((upper < beta) .or. (lower > beta)) &
          .and. .not. (upper < lower))) then
         if (info >= 90000) then
            info = info + 100
         else
            info = 90100
         end if
      end if

   end subroutine check_inputs

   subroutine lcstep &
      (n, m, np, q, npp, &
       f, fjacb, fjacd, &
       wd, ldwd, ld2wd, ss, tt, ldtt, delta, &
       alpha, epsfcn, isodr, &
       tfjacb, omega, u, qraux, kpvt, &
       s, t, phi, irank, rcond, forvcv, &
       wrk1, wrk2, wrk3, wrk4, wrk5, wrk, lwrk, istopc)
   !! Compute locally constrained steps `s` and `t`, and `phi(alpha)`.
      ! @note: This is one of the most time-consuming subroutines in ODRPACK (~25% of total).

      use odrpack_kinds, only: zero, one
      use linpack, only: dchex, dqrdc, dqrsl, dtrco, dtrsl
      use blas_interfaces, only: dnrm2, drot, drotg

      integer, intent(in) :: n
         !! Number of observations.
      integer, intent(in) :: m
         !! Number of columns of data in the explanatory variable.
      integer, intent(in) :: np
         !! Number of function parameters.
      integer, intent(in) :: q
         !! Number of responses per observation.
      integer, intent(in) :: npp
         !! Number of function parameters being estimated.
      real(wp), intent(in) :: f(n, q)
         !! Weighted estimated values of `epsilon`.
      real(wp), intent(in) :: fjacb(n, np, q)
         !! Jacobian with respect to `beta`.
      real(wp), intent(in) :: fjacd(n, m, q)
         !! Jacobian with respect to `delta`.
      real(wp), intent(in) :: wd(ldwd, ld2wd, m)
         !! Squared `delta` weights.
      integer, intent(in) :: ldwd
         !! Leading dimension of array `wd`.
      integer, intent(in) :: ld2wd
         !! Second dimension of array `wd`.
      real(wp), intent(in) :: ss(np)
         !! Scaling values for the unfixed `beta`s.
      real(wp), intent(in) :: tt(ldtt, m)
         !! Scaling values for `delta`.
      integer, intent(in) :: ldtt
         !! Leading dimension of array `tt`.
      real(wp), intent(in) :: delta(n, m)
         !! Estimated errors in the explanatory variables.
      real(wp), intent(in) :: alpha
         !! Levenberg-Marquardt parameter.
      real(wp), intent(in) :: epsfcn
         !! Function's precision.
      logical, intent(in) :: isodr
         !! Variable designating whether the solution is by ODR (`isodr = .true.`) or
         !! by OLS (`isodr = .false.`).
      real(wp), intent(out) :: tfjacb(n, q, np)
         !! Array `omega*fjacb`.
      real(wp), intent(out) :: omega(q, q)
         !! Array defined such that:
         !! `omega*trans(omega) = inv(I + fjacd*inv(e)*trans(fjacd))
         !! = (I - fjacd*inv(p)*trans(fjacd))`
         !! where `e = d**2 + alpha*tt**2` and
         !! `p = trans(fjacd)*fjacd + d**2 + alpha*tt**2`.
      real(wp), intent(out) :: u(np)
         !! Approximate null vector for `tfjacb`.
      real(wp), intent(out) :: qraux(np)
         !! Array required to recover the orthogonal part of the Q-R decomposition.
      integer, intent(out) :: kpvt(np)
         !! Pivot vector.
      real(wp), intent(out) :: s(np)
         !! Step for `beta`.
      real(wp), intent(out) :: t(n, m)
         !! Step for `delta`.
      real(wp), intent(out) :: phi
         !! Difference between the norm of the scaled step and the trust region diameter.
      integer, intent(out) :: irank
         !! Rank deficiency of the Jacobian wrt `beta`.
      real(wp), intent(out) :: rcond
         !! Approximate reciprocal condition number of `tfjacb`.
      logical, intent(in) :: forvcv
         !! Variable designating whether this subroutine was called to set up for the
         !! covariance matrix computations (`forvcv = .true.`) or not (`forvcv = .false.`).
      real(wp), intent(out) :: wrk1(n, q, m)
         !! A work array of `(n, q, m)` elements.
      real(wp), intent(out) :: wrk2(n, q)
         !! A work array of `(n, q)` elements.
      real(wp), intent(out) :: wrk3(np)
         !! A work array of `(np)` elements.
      real(wp), intent(out) :: wrk4(m, m)
         !! A work array of `(m, m)` elements.
      real(wp), intent(out) :: wrk5(m)
         !! A work array of `(m)` elements.
      real(wp), intent(out) :: wrk(lwrk)
         !! A work array of `(lwrk)` elements, _equivalenced_ to `wrk1` and `wrk2`.
      integer, intent(in) :: lwrk
         !! Length of vector `wrk`.
      integer, intent(inout) :: istopc
         !! Variable designating whether the computations were stopped due to a numerical
         !! error within subroutine `dodstp`.

      ! Local scalars
      real(wp) :: co, si, temp
      integer :: i, imax, inf, ipvt, j, k, k1, k2, kp, l
      logical :: elim

      ! Local arrays
      real(wp) :: dum(2)

      ! Variable definitions (alphabetically)
      !  CO:      The cosine from the plane rotation.
      !  DUM:     A dummy array.
      !  ELIM:    The variable designating whether columns of the Jacobian
      !  I:       An indexing variable.
      !  IMAX:    The index of the element of U having the largest absolute value.
      !  INF:     The return code from LINPACK routines.
      !  IPVT:    The variable designating whether pivoting is to be done.
      !  J:       An indexing variable.
      !  K:       An indexing variable.
      !  K1:      An indexing variable.
      !  K2:      An indexing variable.
      !  KP:      The rank of the Jacobian wrt BETA.
      !  L:       An indexing variable.
      !  SI:      The sine from the plane rotation.
      !  TEMP:    A temporary storage LOCATION.

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
         call scale_mat(n, m, wd, ldwd, ld2wd, delta, t)

         do i = 1, n

            !  Compute WRK4, such that TRANS(WRK4)*WRK4 = E = (D**2 + ALPHA*TT**2)
            call esubi(n, m, wd, ldwd, ld2wd, alpha, tt, ldtt, i, wrk4)
            call fctr(.false., wrk4, m, m, inf)
            if (inf /= 0) then
               istopc = 60000
               return
            end if
            ! Compute OMEGA, such that
            ! trans(OMEGA)*OMEGA = I+FJACD*inv(E)*trans(FJACD)
            ! inv(trans(OMEGA)*OMEGA) = I-FJACD*inv(P)*trans(FJACD)
            call vevtr(m, q, i, fjacd, n, m, wrk4, m, wrk1, n, q, omega, q, wrk5)
            do l = 1, q
               omega(l, l) = one + omega(l, l)
            end do
            call fctr(.false., omega, q, q, inf)
            if (inf /= 0) then
               istopc = 60000
               return
            end if
            ! Compute WRK1 = trans(FJACD)*(I-FJACD*inv(P)*trans(JFACD))
            !              = trans(FJACD)*inv(trans(OMEGA)*OMEGA)
            do j = 1, m
               wrk1(i, :, j) = fjacd(i, j, :)
               call solve_trl(q, omega, q, wrk1(i, 1:q, j), 4)
               call solve_trl(q, omega, q, wrk1(i, 1:q, j), 2)
            end do

            ! Compute WRK5 = inv(E)*D*G2
            wrk5 = t(i, :)
            call solve_trl(m, wrk4, m, wrk5, 4)
            call solve_trl(m, wrk4, m, wrk5, 2)

            ! Compute TFJACB = inv(trans(OMEGA))*FJACB
            do k = 1, kp
               tfjacb(i, :, k) = fjacb(i, kpvt(k), :)
               call solve_trl(q, omega, q, tfjacb(i, 1:q, k), 4)
               if (ss(1) > zero) then
                  tfjacb(i, :, k) = tfjacb(i, :, k)/ss(kpvt(k))
               else
                  tfjacb(i, :, k) = tfjacb(i, :, k)/abs(ss(1))
               end if
            end do

            ! Compute WRK2 = (V*inv(E)*D**2*G2 - G1)
            do l = 1, q
               wrk2(i, l) = dot_product(fjacd(i, :, l), wrk5) - f(i, l)
            end do

            ! Compute WRK2 = inv(trans(OMEGA))*(V*inv(E)*D**2*G2 - G1)
            call solve_trl(q, omega, q, wrk2(i, 1:q), 4)

         end do

      else
         do l = 1, q
            do i = 1, n
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
         kpvt = 0
      else
         ipvt = 0
      end if

      call dqrdc(tfjacb, n*q, n*q, kp, qraux, kpvt, wrk3, ipvt)
      call dqrsl(tfjacb, n*q, n*q, kp, qraux, wrk2, dum, wrk2, dum, dum, dum, 1000, inf)
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
                  call drot(kp - k2, tfjacb(k2, 1, k2 + 1), n*q, wrk3(k2 + 1), 1, co, si)
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
               call dtrco(tfjacb, n*q, kp, rcond, u, 1)
               if (rcond <= epsfcn) then
                  elim = .true.
                  imax = maxloc(u(1:kp), dim=1)
                  ! IMAX is the column to remove - use DCHEX and fix KPVT
                  if (imax /= kp) then
                     call dchex(tfjacb, n*q, kp, imax, kp, wrk2, n*q, 1, qraux, wrk3, 2)
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
            call dtrsl(tfjacb, n*q, kp, wrk2, 01, inf)
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
            call esubi(n, m, wd, ldwd, ld2wd, alpha, tt, ldtt, i, wrk4)
            call fctr(.false., wrk4, m, m, inf)
            if (inf /= 0) then
               istopc = 60000
               return
            end if

            ! Compute WRK5 = inv(E)*D*G2
            wrk5 = t(i, :)
            call solve_trl(m, wrk4, m, wrk5, 4)
            call solve_trl(m, wrk4, m, wrk5, 2)

            do l = 1, q
               wrk2(i, l) = f(i, l) &
                            + dot_product(fjacb(i, 1:npp, l), s(1:npp)) &
                            - dot_product(fjacd(i, :, l), wrk5)
            end do

            do j = 1, m
               wrk5(j) = dot_product(wrk1(i, :, j), wrk2(i, :))
               t(i, j) = -(wrk5(j) + t(i, j))
            end do

            call solve_trl(m, wrk4, m, t(i, 1:m), 4)
            call solve_trl(m, wrk4, m, t(i, 1:m), 2)

         end do

      end if

      ! Compute PHI(ALPHA) from scaled S and T
      call scale_mat(npp, 1, ss, npp, 1, s, wrk(1:npp))
      if (isodr) then
         call scale_mat(n, m, tt, ldtt, 1, t, wrk(npp + 1:npp + 1 + n*m - 1))
         phi = dnrm2(npp + n*m, wrk, 1)
      else
         phi = dnrm2(npp, wrk, 1)
      end if

   end subroutine lcstep

   subroutine vcv_beta &
      (n, m, np, q, npp, &
       f, fjacb, fjacd, &
       wd, ldwd, ld2wd, ssf, ss, tt, ldtt, delta, &
       epsfcn, isodr, &
       vcv, sd, &
       wrk6, omega, u, qraux, jpvt, &
       s, t, irank, rcond, rss, idf, rvar, ifixb, &
       wrk1, wrk2, wrk3, wrk4, wrk5, wrk, lwrk, istopc)
   !! Compute covariance matrix of estimated parameters.

      use odrpack_kinds, only: zero
      use linpack, only: dpodi

      integer, intent(in) :: n
         !! Number of observations.
      integer, intent(in) :: m
         !! Number of columns of data in the explanatory variable.
      integer, intent(in) :: np
         !! Number of function parameters.
      integer, intent(in) :: q
         !! Number of responses per observation.
      integer, intent(in) :: npp
         !! Number of function parameters being estimated.
      real(wp), intent(in) :: f(n, q)
         !! Weighted estimated values of `epsilon`.
      real(wp), intent(in) :: fjacb(n, np, q)
         !! Jacobian with respect to `beta`.
      real(wp), intent(in) :: fjacd(n, m, q)
         !! Jacobian with respect to `delta`.
      real(wp), intent(in) :: wd(ldwd, ld2wd, m)
         !! `delta` weights.
      integer, intent(in) :: ldwd
         !! Leading dimension of array `wd`.
      integer, intent(in) :: ld2wd
         !! Second dimension of array `wd`.
      real(wp), intent(in) :: ssf(np)
         !! Scaling values used for `beta`.
      real(wp), intent(in) :: ss(np)
         !! Scaling values for the unfixed `beta`s.
      real(wp), intent(in) :: tt(ldtt, m)
         !! Scaling values for `delta`.
      integer, intent(in) :: ldtt
         !! Leading dimension of array `tt`.
      real(wp), intent(in) :: delta(n, m)
         !! Estimated errors in the explanatory variables.
      real(wp), intent(in) :: epsfcn
         !! Function's precision.
      logical, intent(in) :: isodr
         !! Variable designating whether the solution is by ODR (`isodr = .true.`) or
         !! by OLS (`isodr = .false.`).
      real(wp), intent(out) :: vcv(np, np)
         !! Covariance matrix of the estimated `beta`s.
      real(wp), intent(out) :: sd(np)
         !! Standard deviations of the estimated `beta`s.
      real(wp), intent(out) :: wrk6(n*q, np)
         !! A work array of `(n*q, np)` elements.
      real(wp), intent(out) :: omega(q, q)
         !! Array defined such that `omega*trans(omega) = inv(I + fjacd*inv(e)*trans(fjacd))
         !! = (I - fjacd*inv(p)*trans(fjacd))`.
      real(wp), intent(out) :: u(np)
         !! Approximate null vector for `fjacb`.
      real(wp), intent(out) :: qraux(np)
         !! Array required to recover the orthogonal part of the Q-R decomposition.
      integer, intent(out) :: jpvt(np)
         !! Pivot vector.
      real(wp), intent(out) :: s(np)
         !! Step for `beta`.
      real(wp), intent(out) :: t(n, m)
         !! Step for `delta`.
      integer, intent(out) :: irank
         !! Rank deficiency of the Jacobian wrt `beta`.
      real(wp), intent(out) :: rcond
         !! Approximate reciprocal condition of `fjacb`.
      real(wp), intent(inout) :: rss
         !! Residual sum of squares.
      integer, intent(out) :: idf
         !! Degrees of freedom of the fit, equal to the number of observations with nonzero
         !! weighted derivatives minus the number of parameters being estimated.
      real(wp), intent(out) :: rvar
         !! Residual variance.
      integer, intent(in) :: ifixb(np)
         !! Values designating whether the elements of `beta` are fixed at their input values or not.
      real(wp), intent(out) :: wrk1(n, q, m)
         !! A work array of `(n, q, m)` elements.
      real(wp), intent(out) :: wrk2(n, q)
         !! A work array of `(n, q)` elements.
      real(wp), intent(out) :: wrk3(np)
         !! A work array of `(np)` elements.
      real(wp), intent(out) :: wrk4(m, m)
         !! A work array of `(m, m)` elements.
      real(wp), intent(out) :: wrk5(m)
         !! A work array of `(m)` elements.
      real(wp), intent(out) :: wrk(lwrk)
         !! A work array of `(lwrk)` elements, _equivalenced_ to `wrk1` and `wrk2`.
      integer, intent(in) :: lwrk
         !! Length of vector `lwrk`.
      integer, intent(out) :: istopc
         !! Variable designating whether the computations were stoped due to a numerical
         !! error within subroutine `dodstp`.

      ! Local scalars
      real(wp) :: temp
      integer :: i, iunfix, j, junfix, kp
      logical :: forvcv

      ! Variable definitions (alphabetically)
      !  FORVCV:  The variable designating whether subroutine DODSTP is called to set up for the
      !           covariance matrix computations (FORVCV=TRUE) or not (FORVCV=FALSE).
      !  I:       An indexing variable.
      !  IUNFIX:  The index of the next unfixed parameter.
      !  J:       An indexing variable.
      !  JUNFIX:  The index of the next unfixed parameter.
      !  KP:      The rank of the Jacobian wrt BETA.
      !  TEMP:    A temporary storage location

      forvcv = .true.
      istopc = 0

      call lcstep(n, m, np, q, npp, &
                  f, fjacb, fjacd, &
                  wd, ldwd, ld2wd, ss, tt, ldtt, delta, &
                  zero, epsfcn, isodr, &
                  wrk6, omega, u, qraux, jpvt, &
                  s, t, temp, irank, rcond, forvcv, &
                  wrk1, wrk2, wrk3, wrk4, wrk5, wrk, lwrk, istopc)
      if (istopc /= 0) then
         return
      end if

      kp = npp - irank
      call dpodi(wrk6, n*q, kp, wrk3, 1)

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
      sd = zero
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
      vcv = zero
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

   end subroutine vcv_beta

   pure subroutine pack_vec(n2, n1, v1, v2, ifix)
   !! Select the unfixed elements of `v2` and return them in `v1`.

      integer, intent(in) :: n2
         !! Number of items in `v2`.
      integer, intent(out) :: n1
         !! Number of items in `v1`.
      real(wp), intent(out) :: v1(n2)
         !! Vector of the unfixed items from `v2`.
      real(wp), intent(in) :: v2(n2)
         !! Vector of the fixed and unfixed items from which the unfixed elements are to be extracted.
      integer, intent(in) :: ifix(n2)
         !! Values designating whether the elements of `v2` are fixed at their input values or not.

      ! Local scalars
      integer :: i

      ! Variable definitions (alphabetically)
      !  I:       An indexing variable.

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
         v1 = v2
      end if

   end subroutine pack_vec

   pure subroutine unpack_vec(n2, v1, v2, ifix)
   !! Copy the elements of `v1` into the locations of `v2` which are unfixed.

      integer, intent(in) :: n2
         !! Number of items in `v2`.
      real(wp), intent(in) :: v1(n2)
         !! Vector of the unfixed items.
      real(wp), intent(out) :: v2(n2)
         !! Vector of the fixed and unfixed items into which the elements of `v1` are to
         !! be inserted.
      integer, intent(in) :: ifix(n2)
         !! Values designating whether the elements of `v2` are fixed at their input values
         !! or not.

      ! Local scalars
      integer :: i, n1

      ! Variable Definitions (alphabetically)
      !  I:       An indexing variable.
      !  N1:      The number of items in V1.

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
         v2 = v1
      end if

   end subroutine unpack_vec

   real(wp) pure function ppf_normal(p) result(res)
   !! Compute the percent point function value for the normal (Gaussian) distribution with
   !!  mean 0 and standard deviation 1, and with probability density function:
   !!
   !!       `f(x) = (1/sqrt(2*pi))*exp(-x^2/2)`
   !!
      ! Adapted from DATAPAC subroutine TPPF, with modifications to facilitate conversion to
      ! real(wp) automatically.
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

      use odrpack_kinds, only: zero, half, one

      real(wp), intent(in) :: p
         !! Probability at which the percent point is to be evaluated. `p` must lie between
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
         res = zero
      else
         r = p
         if (p > half) r = one - r
         t = sqrt(-2*log(r))
         anum = ((((t*p4 + p3)*t + p2)*t + p1)*t + p0)
         aden = ((((t*q4 + q3)*t + q2)*t + q1)*t + q0)
         res = t + (anum/aden)
         if (p < half) res = -res
      end if

   end function ppf_normal

   real(wp) pure function ppf_tstudent(p, idf) result(res)
   !! Compute the percent point function value for the student's T distribution with `idf`
   !! degrees of freedom.
      ! Adapted from DATAPAC subroutine TPPF, with modifications to facilitate conversion to
      ! real(wp) automatically.
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
         !! Probability at which the percent point is to be evaluated. `p` must lie between
         !! 0.0 and 1.0, exclusive.
      integer, intent(in) :: idf
         !! Degrees of freedom.

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
      !  IPASS:  A value used in the approximation.
      !  MAXIT:  The maximum number of iterations allowed for the approx.
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
         res = zero

      elseif (idf == 1) then
         !Treat the IDF = 1 (Cauchy) case
         arg = pi*p
         res = -1/tan(arg)

      elseif (idf == 2) then
         !  Treat the IDF = 2 case
         term1 = sqrt(two)/2
         term2 = 2*p - one
         term3 = sqrt(p*(one - p))
         res = term1*term2/term3

      elseif (idf >= 3) then
         ! Treat the IDF greater than or equal to 3 case
         ppfn = ppf_normal(p)
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
         res = term1 + term2 + term3 + term4 + term5

         if (idf == 3) then
            ! Augment the results for the IDF = 3 case
            con = pi*(p - half)
            arg = res/sqrt(df)
            z = atan(arg)
            do ipass = 1, maxit
               s = sin(z)
               c = cos(z)
               z = z - (z + s*c - con)/(2*c**2)
            end do
            res = sqrt(df)*s/c

         elseif (idf == 4) then
            ! Augment the results for the IDF = 4 case
            con = 2*(p - half)
            arg = res/sqrt(df)
            z = atan(arg)
            do ipass = 1, maxit
               s = sin(z)
               c = cos(z)
               z = z - ((one + half*c**2)*s - con)/((one + half)*c**3)
            end do
            res = sqrt(df)*s/c

         elseif (idf == 5) then
            ! Augment the results for the IDF = 5 case
            con = pi*(p - half)
            arg = res/sqrt(df)
            z = atan(arg)
            do ipass = 1, maxit
               s = sin(z)
               c = cos(z)
               z = z - (z + (c + (two/three)*c**3)*s - con)/((eight/three)*c**4)
            end do
            res = sqrt(df)*s/c

         elseif (idf == 6) then
            !  Augment the results for the IDF = 6 case
            con = 2*(p - half)
            arg = res/sqrt(df)
            z = atan(arg)
            do ipass = 1, maxit
               s = sin(z)
               c = cos(z)
               z = z - ((one + half*c**2 + (three/eight)*c**4)*s - con)/((fiftn/eight)*c**5)
            end do
            res = sqrt(df)*s/c
         end if
      end if

   end function ppf_tstudent

   subroutine fpvb &
      (fcn, &
       n, m, np, q, &
       beta, xplusd, ifixb, ifixx, ldifx, &
       nrow, j, lq, stp, &
       istop, nfev, pvb, &
       wrk1, wrk2, wrk6)
   !! Compute the `nrow`-th function value using `beta(j) + stp`.

      procedure(fcn_t) :: fcn
         !! User-supplied subroutine for evaluating the model.
      integer, intent(in) :: n
         !! Number of observations.
      integer, intent(in) :: m
         !! Number of columns of data in the independent variable.
      integer, intent(in) :: np
         !! Number of function parameters.
      integer, intent(in) :: q
         !! Number of responses per observation.
      real(wp), intent(inout) :: beta(np)
         !! Function parameters.
      real(wp), intent(in) :: xplusd(n, m)
         !! Values of `x + delta`.
      integer, intent(in) :: ifixb(np)
         !! Values designating whether the elements of `beta` are fixed at their input values or not.
      integer, intent(in) :: ifixx(ldifx, m)
         !! Values designating whether the elements of `x` are fixed at their input values or not.
      integer, intent(in) :: ldifx
         !! Leading dimension of array `ifixx`.
      integer, intent(in) :: nrow
         !! Row number of the independent variable array at which the derivative is to be checked.
      integer, intent(in) :: j
         !! Index of the partial derivative being examined.
      integer, intent(in) :: lq
         !! Response currently being examined.
      real(wp), intent(in) :: stp
         !! Step size for the finite difference derivative.
      integer, intent(out) :: istop
         !! Variable designating whether there are problems computing the function at the
         !! current `beta` and `delta`.
      integer, intent(inout) :: nfev
         !! Number of function evaluations.
      real(wp), intent(out) :: pvb
         !! Function value for the selected observation & response.
      real(wp), intent(out) :: wrk1(n, m, q)
         !! Work array.
      real(wp), intent(out) :: wrk2(n, q)
         !! Work array.
      real(wp), intent(out) :: wrk6(n, np, q)
         !! Work array.

      ! Local scalars
      real(wp) :: betaj

      ! Variable Definitions (alphabetically)
      !  BETAJ:   The current estimate of the jth parameter.

      betaj = beta(j)
      beta(j) = beta(j) + stp

      istop = 0
      call fcn(beta, xplusd, ifixb, ifixx, 003, wrk2, wrk6, wrk1, istop)
      if (istop == 0) then
         nfev = nfev + 1
      else
         return
      end if

      beta(j) = betaj
      pvb = wrk2(nrow, lq)

   end subroutine fpvb

   subroutine fpvd &
      (fcn, &
       n, m, np, q, &
       beta, xplusd, ifixb, ifixx, ldifx, &
       nrow, j, lq, stp, &
       istop, nfev, pvd, &
       fjacd, f, fjacb)
   !! Compute `nrow`-th function value using `x(nrow, j) + delta(nrow, j) + stp`.

      procedure(fcn_t) :: fcn
         !! User-supplied subroutine for evaluating the model.
      integer, intent(in) :: n
         !! Number of observations.
      integer, intent(in) :: m
         !! Number of columns of data in the independent variable.
      integer, intent(in) :: np
         !! Number of function parameters.
      integer, intent(in) :: q
         !! Number of responses per observation.
      real(wp), intent(in) :: beta(np)
         !! Function parameters.
      real(wp), intent(inout) :: xplusd(n, m)
         !! Values of `x + delta`.
      integer, intent(in) :: ifixb(np)
         !! Values designating whether the elements of `beta` are fixed at their input values or not.
      integer, intent(in) :: ifixx(ldifx, m)
         !! Values designating whether the elements of `x` are fixed at their input values or not.
      integer, intent(in) :: ldifx
         !! Leading dimension of array `ifixx`.
      integer, intent(in) :: nrow
         !! Row number of the independent variable array at which the derivative is to be checked.
      integer, intent(in) :: j
         !! Index of the partial derivative being examined.
      integer, intent(in) :: lq
         !! Response currently being examined.
      real(wp), intent(in) :: stp
         !! Step size for the finite difference derivative.
      integer, intent(out) :: istop
         !! Variable designating whether there are problems computing the function at the
         !! current `beta` and `delta`.
      integer, intent(inout) :: nfev
         !! Number of function evaluations.
      real(wp), intent(out) :: pvd
         !! Function value for the selected observation & response.
      real(wp), intent(out) :: fjacd(n, m, q)
         !! Jacobian wrt delta.
      real(wp), intent(out) :: f(n, q)
         !! Predicted function values.
      real(wp), intent(out) :: fjacb(n, np, q)
         !! Jocobian wrt beta.

      ! Local scalars
      real(wp) :: xpdj

      ! Variable Definitions (alphabetically)
      !  XPDJ:    The (NROW,J)th element of XPLUSD.

      xpdj = xplusd(nrow, j)
      xplusd(nrow, j) = xplusd(nrow, j) + stp

      istop = 0
      call fcn(beta, xplusd, ifixb, ifixx, 003, f, fjacb, fjacd, istop)
      if (istop == 0) then
         nfev = nfev + 1
      else
         return
      end if

      xplusd(nrow, j) = xpdj
      pvd = f(nrow, lq)

   end subroutine fpvd

   pure subroutine scale_vec(n, m, scl, ldscl, t, ldt, sclt, ldsclt)
   !! Scale `t` by the inverse of `scl`, i.e., compute `t/scl`.

      use odrpack_kinds, only: zero

      integer, intent(in) :: n
         !! Number of rows of data in `t`.
      integer, intent(in) :: m
         !! Number of columns of data in `t`.
      real(wp), intent(in) :: scl(ldscl, m)
         !! Scale values.
      integer, intent(in) :: ldscl
         !! Leading dimension of array `scl`.
      real(wp), intent(in) :: t(ldt, m)
         !! Array to be inversely scaled by `scl`.
      integer, intent(in) :: ldt
         !! Leading dimension of array `t`.
      real(wp), intent(out) :: sclt(ldsclt, m)
         !! Inversely scaled matrix.
      integer, intent(in) :: ldsclt
         !! Leading dimension of array `sclt`.

      ! Local scalars
      integer :: j

      ! Variable Definitions (alphabetically)
      !  J:       An indexing variable.

      if (n == 0 .or. m == 0) return

      if (scl(1, 1) >= zero) then
         if (ldscl >= n) then
            sclt(1:n, 1:m) = t(1:n, 1:m)/scl(1:n, 1:m)
         else
            do j = 1, m
               sclt(1:n, j) = t(1:n, j)/scl(1, j)
            end do
         end if
      else
         sclt(1:n, 1:m) = t(1:n, 1:m)/abs(scl(1, 1))
      end if

   end subroutine scale_vec

   pure subroutine scale_beta(np, beta, ssf)
   !! Select scaling values for `beta` according to the algorithm given in the ODRPACK95
   !! reference guide.

      use odrpack_kinds, only: zero, one

      integer, intent(in) :: np
         !! Number of function parameters.
      real(wp), intent(in) :: beta(np)
         !! Function parameters.
      real(wp), intent(out) :: ssf(np)
         !! Scaling values for `beta`.

      ! Local scalars
      real(wp) :: bmax, bmin
      integer :: k
      logical ::bigdif

      ! Variable Definitions (alphabetically)
      !  BIGDIF:  The variable designating whether there is a significant difference in the
      !           magnitudes of the nonzero elements of BETA (BIGDIF=.TRUE.) or not (BIGDIF=.FALSE.).
      !  BMAX:    The largest nonzero magnitude.
      !  BMIN:    The smallest nonzero magnitude.
      !  K:       An indexing variable.

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
               ssf(k) = 10/bmin
            else
               if (bigdif) then
                  ssf(k) = 1/abs(beta(k))
               else
                  ssf(k) = 1/bmax
               end if
            end if
         end do

      end if

   end subroutine scale_beta

   pure subroutine scale_delta(n, m, x, tt, ldtt)
   !! Select scaling values for `delta` according to the algorithm given in the ODRPACK95
   !! reference guide.

      use odrpack_kinds, only: zero, one

      integer, intent(in) :: n
         !! Number of observations.
      integer, intent(in) :: m
         !! Number of columns of data in the independent variable.
      real(wp), intent(in) :: x(n, m)
         !! Independent variable.
      real(wp), intent(out) :: tt(ldtt, m)
         !! Scaling values for `delta`.
      integer, intent(in) :: ldtt
         !! Leading dimension of array `tt`.

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
      !  XMAX:    The largest nonzero magnitude.
      !  XMIN:    THE SMALLEST NONZERO MAGNITUDE.

      do j = 1, m
         xmax = abs(x(1, j))
         do i = 2, n
            xmax = max(xmax, abs(x(i, j)))
         end do
         if (xmax == zero) then
            !  All input values of X(I,J), I=1,...,N, are zero
            tt(1:n, j) = one
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
                     tt(i, j) = 1/abs(x(i, j))
                  else
                     tt(i, j) = 1/xmax
                  end if
               else
                  tt(i, j) = 10/xmin
               end if
            end do
         end if
      end do

   end subroutine scale_delta

   pure subroutine select_row(n, m, x, nrow)
   !! Select the row at which the derivative will be checked.

      use odrpack_kinds, only: zero

      integer, intent(in) :: n
         !! Number of observations.
      integer, intent(in) :: m
         !! Number of columns of data in the independent variable.
      real(wp), intent(in) :: x(n, m)
         !! Independent variable.
      integer, intent(inout) :: nrow
         !! Selected row number of the independent variable.

      ! Local scalars
      integer :: i

      ! Variable Definitions (alphabetically)
      !  I:       An index variable.
      !  J:       An index variable.
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

   end subroutine select_row

   pure subroutine solve_trl(n, t, ldt, b, job)
   !! Solve systems of the form:
   !!
   !!  `t * x = b  or  trans(t) * x = b`
   !!
   !! where `t` is an upper or lower triangular matrix of order `n`, and the solution `x`
   !! overwrites the RHS `b`.
      ! Adapted from LINPACK subroutine DTRSL.
      ! @note: This is one of the most time-consuming subroutines in ODRPACK (~25% of total).

      use odrpack_kinds, only: zero

      integer, intent(in) :: n
         !! Number of rows and columns of data in array `t`.
      real(wp), intent(in) :: t(ldt, n)
         !! Upper or lower tridiagonal system.
      integer, intent(in) :: ldt
         !! Leading dimension of array `t`.
      real(wp), intent(inout) :: b(:)
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

      select case (job)
      case (1)
         ! Solve T*X=B for T lower triangular
         b(j1) = b(j1)/t(j1, j1)
         do j = j1 + 1, jn
            temp = -b(j - 1)
            b(j:jn) = b(j:jn) + temp*t(j:jn, j - 1)
            if (t(j, j) /= zero) then
               b(j) = b(j)/t(j, j)
            else
               b(j) = zero
            end if
         end do

      case (2)
         ! Solve T*X=B for T upper triangular.
         b(jn) = b(jn)/t(jn, jn)
         do j = jn - 1, j1, -1
            temp = -b(j + 1)
            b(1:j) = b(1:j) + temp*t(1:j, j + 1)
            if (t(j, j) /= zero) then
               b(j) = b(j)/t(j, j)
            else
               b(j) = zero
            end if
         end do

      case (3)
         ! Solve trans(T)*X=B for T lower triangular.
         b(jn) = b(jn)/t(jn, jn)
         do j = jn - 1, j1, -1
            b(j) = b(j) - dot_product(t(j + 1:jn + 1, j), b(j + 1:jn + 1))
            if (t(j, j) /= zero) then
               b(j) = b(j)/t(j, j)
            else
               b(j) = zero
            end if
         end do

      case (4)
         ! Solve trans(T)*X=B for T upper triangular
         b(j1) = b(j1)/t(j1, j1)
         do j = j1 + 1, jn
            b(j) = b(j) - dot_product(t(1:j - 1, j), b(1:j - 1))
            if (t(j, j) /= zero) then
               b(j) = b(j)/t(j, j)
            else
               b(j) = zero
            end if
         end do
      case default
         error stop "Invalid value of JOB."
      end select

   end subroutine solve_trl

   pure subroutine vevtr &
      (m, q, indx, v, ldv, ld2v, e, lde, ve, ldve, ld2ve, vev, ldvev, wrk5)
   !! Compute `v*e*trans(v)` for the (`indx`)th `m` by `q` array in `v`.

      integer, intent(in) :: m
         !! Number of columns of data in the independent variable.
      integer, intent(in) :: q
         !! Number of responses per observation.
      integer, intent(in) :: indx
         !! Row in `v` in which the `m` by `q` array is stored.
      integer, intent(in) :: ldv
         !! Leading dimension of array `v`.
      integer, intent(in) :: ld2v
         !! Second dimension of array `v`.
      integer, intent(in) :: lde
         !! Leading dimension of array `e`.
      integer, intent(in) :: ldve
         !! Leading dimension of array `ve`.
      integer, intent(in) :: ldvev
         !! Leading dimension of array `vev`.
      integer, intent(in) :: ld2ve
         !! Second dimension of array `ve`.
      real(wp), intent(in) :: v(ldv, ld2v, q)
         !! An array of `q` by `m` matrices.
      real(wp), intent(in) :: e(lde, m)
         !! Matrix of the factors, so `ete = (d**2 + alpha*t**2)`.
      real(wp), intent(out) :: ve(ldve, ld2ve, m)
         !! Array `ve = v * inv(e)`.
      real(wp), intent(out) :: vev(ldvev, q)
         !! Array `vev = v * inv(ete) * trans(v)`.
      real(wp), intent(out) :: wrk5(m)
         !! Work vector.

      ! Local scalars
      integer :: l1, l2

      ! Variable Definitions (alphabetically)
      !  J:       An indexing variable.
      !  L1:      An indexing variable.
      !  L2:      An indexing variable.

      if (q == 0 .or. m == 0) return

      do l1 = 1, q
         wrk5 = v(indx, 1:m, l1)
         call solve_trl(m, e, lde, wrk5, 4)
         ve(indx, l1, :) = wrk5
      end do

      do l1 = 1, q
         do l2 = 1, l1
            vev(l1, l2) = dot_product(ve(indx, l1, :), ve(indx, l2, :))
            vev(l2, l1) = vev(l1, l2)
         end do
      end do

   end subroutine vevtr

   pure subroutine scale_mat(n, m, wt, ldwt, ld2wt, t, wtt)
   !! Scale matrix `t` using `wt`, i.e., compute `wtt = wt*t`.

      use odrpack_kinds, only: zero

      integer, intent(in) :: n
         !! Number of rows of data in `t`.
      integer, intent(in) :: m
         !! Number of columns of data in `t`.
      real(wp), intent(in), target :: wt(..)
         !! Array of shape conformable to `(ldwt,ld2wt,m)` holding the weights.
      integer, intent(in) :: ldwt
         !! Leading dimension of array `wt`.
      integer, intent(in) :: ld2wt
         !! Second dimension of array `wt`.
      real(wp), intent(in), target :: t(..)
         !! Array of shape conformable to `(n,m)` being scaled by `wt`.
      real(wp), intent(out), target :: wtt(..)
         !! Array of shape conformable to `(n,m)` holding the result of weighting array `t` by
         !! array `wt`. Array `wtt` can be the same as `t` only if the arrays in `wt` are upper
         !! triangular with zeros below the diagonal.

      ! Local scalars
      integer :: i, j
      real(wp), pointer :: wt_(:, :, :), t_(:, :), wtt_(:, :)

      ! Variable Definitions (alphabetically)
      !  I:       An indexing variable.
      !  J:       An indexing variable.

      if (n == 0 .or. m == 0) return

      select rank (wt)
      rank (1); wt_(1:ldwt, 1:ld2wt, 1:m) => wt
      rank (2); wt_(1:ldwt, 1:ld2wt, 1:m) => wt
      rank (3); wt_(1:ldwt, 1:ld2wt, 1:m) => wt
      rank default; error stop "Invalid rank of `wt`."
      end select

      select rank (t)
      rank (1); t_(1:n, 1:m) => t
      rank (2); t_(1:n, 1:m) => t
      rank default; error stop "Invalid rank of `t`."
      end select

      select rank (wtt)
      rank (1); wtt_(1:n, 1:m) => wtt
      rank (2); wtt_(1:n, 1:m) => wtt
      rank default; error stop "Invalid rank of `wtt`."
      end select

      if (wt_(1, 1, 1) >= zero) then
         if (ldwt >= n) then
            if (ld2wt >= m) then
               ! WT is an N-array of M by M matrices
               do j = 1, m
                  do i = 1, n
                     wtt_(i, j) = dot_product(wt_(i, j, :), t_(i, :))
                  end do
               end do
            else
               ! WT is an N-array of diagonal matrices
               do j = 1, m
                  wtt_(:, j) = wt_(:, 1, j)*t_(:, j)
               end do
            end if
         else
            if (ld2wt >= m) then
               ! WT is an M by M matrix
               do j = 1, m
                  do i = 1, n
                     wtt_(i, j) = dot_product(wt_(1, j, :), t_(i, :))
                  end do
               end do
            else
               ! WT is a diagonal matrix
               do j = 1, m
                  wtt_(:, j) = wt_(1, 1, j)*t_(:, j)
               end do
            end if
         end if
      else
         ! WT is a scalar
         wtt_ = abs(wt_(1, 1, 1))*t_
      end if

   end subroutine scale_mat

   pure subroutine move_beta( &
      np, beta, lower, upper, ssf, stpb, neta, eta, interval)
   !! Ensure range of bounds is large enough for derivative checking.
   !! Move beta away from bounds so that derivatives can be calculated.

      use odrpack_kinds, only: zero, one, three, ten, hundred

      integer, intent(in) :: np
         !! Number of parameters `np`.
      real(wp), intent(inout) :: beta(np)
         !! Function parameters.
      real(wp), intent(in) :: lower(np)
         !! !! Lower bound on `beta`.
      real(wp), intent(in) :: upper(np)
         !! Upper bound on `beta`.
      real(wp), intent(in) :: ssf(np)
         !! Scale used for the `beta`s.
      real(wp), intent(in) :: stpb(np)
         !! Relative step used for computing finite difference derivatives with respect to `beta`.
      integer, intent(in) :: neta
         !! Number of good digits in the function results.
      real(wp), intent(in) :: eta
         !! Relative noise in the function results.
      integer, intent(out) :: interval(np)
         !! Specifies which difference methods and step sizes are supported by the current
         !! interval `upper-lower`.

      ! Local scalars
      integer :: k
      real(wp) :: h, h0, h1, hc, hc0, hc1, stpl, stpr, typj

      ! VARIABLE DEFINITIONS (ALPHABETICALLY)
      !  H:        Relative step size for forward differences.
      !  H0:       Initial relative step size for forward differences.
      !  H1:       Default relative step size for forward differences.
      !  HC:       Relative step size for center differences.
      !  HC0:      Initial relative step size for center differences.
      !  HC1:      Default relative step size for center differences.
      !  K:        Index variable for BETA.
      !  STPL:     Maximum step to the left of BETA (-) the derivative checker will use.
      !  STPR:     Maximum step to the right of BETA (+) the derivative checker will use.
      !  TYPJ:     The typical size of the J-th unkonwn BETA.

      interval = 111
      do k = 1, np
         h0 = hstep(0, neta, 1, k, stpb, 1)
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

   end subroutine move_beta

end module odrpack_core
