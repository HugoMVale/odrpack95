subroutine dacces &
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
! Routines Called  DIWINF,DWINF
! Date Written   860529   (YYMMDD)
! Revision Date  920619   (YYMMDD)

   use odrpack_kinds, only: wp

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
   real(kind=wp), intent(inout) :: work(lwork)
      !! The `real` (kind=wp) work space.
   integer, intent(in) :: lwork
      !! The length of vector `work`.
   integer, intent(inout) :: iwork(liwork)
      !! The integer work space.
   integer, intent(in) :: liwork
      !! The length of vector `iwork`.
   logical, intent(in) :: access
      !! The variable designating whether information is to be accessed from the work
      !! arrays (`access`=TRUE) or stored in them (`access`=FALSE).
   logical, intent(in) :: isodr
      !! The variable designating whether the solution is to be found by ODR (`isodr`=TRUE) or
      !! by OLS (`isodr`=FALSE).
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
   real(kind=wp), intent(out) :: partol
      !! The parameter convergence stopping tolerance.
   real(kind=wp), intent(out) :: sstol
      !! The sum-of-squares convergence stopping tolerance.
   integer, intent(out) :: maxit
      !! The maximum number of iterations allowed.
   real(kind=wp), intent(out) :: taufac
      !! The factor used to compute the initial trust region diameter.
   real(kind=wp), intent(out) :: eta
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
   real(kind=wp), intent(out) :: wss(3)
      !! The sum of the squares of the weighted `epsilons` and `deltas`, the sum of the squares
      !! of the weighted `deltas`, and the sum of the squares of the weighted `epsilons`.
   real(kind=wp), intent(out) :: rvar
      !! The residual variance, i.e. standard deviation squared.
   integer, intent(out) :: idf
      !! The degrees of freedom of the fit, equal to the number of observations with nonzero
      !! weighted derivatives minus the number of parameters being estimated.
   real(kind=wp), intent(out) :: tau
      !! The trust region diameter.
   real(kind=wp), intent(out) :: alpha
      !! The Levenberg-Marquardt parameter.
   integer, intent(out) :: niter
      !! The number of iterations taken.
   integer, intent(out) :: nfev
      !! The number of function evaluations.
   integer, intent(out) :: njev
      !! The number of Jacobian evaluations.
   integer, intent(out) :: int2
      !! The number of internal doubling steps.
   real(kind=wp), intent(out) :: olmavg
      !! The average number of Levenberg-Marquardt steps per iteration.
   real(kind=wp), intent(out) :: rcond
      !! The approximate reciprocal condition of `fjacb`.
   integer, intent(out) :: irank
      !! The rank deficiency of the Jacobian wrt `beta`.
   real(kind=wp), intent(out) :: actrs
      !! The saved actual relative reduction in the sum-of-squares.
   real(kind=wp), intent(out) :: pnorm
      !! The norm of the scaled estimated parameters.
   real(kind=wp), intent(out) :: prers
      !! The saved predicted relative reduction in the sum-of-squares.
   real(kind=wp), intent(out) :: rnorms
      !! The norm of the saved weighted `epsilons` and `deltas`.
   integer, intent(out) :: istop
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

   ! External subroutines
   external :: diwinf, dwinf

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
   !  IDF:     The degrees of freedom of the fit, equal to the number of observations with nonzero
   !           weighted derivatives minus the number of parameters being estimated.
   !  IDFI:    The starting location in array IWORK of variable IDF.
   !  INT2:    The number of internal doubling steps.
   !  INT2I:   The location in array IWORK of variable INT2.
   !  IPR1:    The value of the fourth digit (from the right) of IPRINT, which controls the
   !           initial summary report.
   !  IPR2:    The value of the third digit (from the right) of IPRINT, which controls the
   !           iteration reports.
   !  IPR2F:   The value of the second digit (from the right) of IPRINT, which controls the
   !           frequency of the iteration reports.
   !  IPR3:    The value of the first digit (from the right) of IPRINT, which controls the final
   !           summary report.
   !  IPRINI:  The location in array IWORK of variable IPRINT.
   !  IPRINT:  The print control variable.
   !  IRANK:   The rank deficiency of the Jacobian wrt BETA.
   !  IRANKI:  The location in array IWORK of variable IRANK.
   !  ISODR:   The variable designating whether the solution is to be found by ODR (ISODR=TRUE) or
   !           by OLS (ISODR=FALSE).
   !  ISTOP:   The variable designating whether there are problems computing the function at the
   !           current BETA and DELTA.
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
   !  RESTRT:  The variable designating whether the call is a restart (RESTRT=TRUE) or not (RESTRT=FALSE).
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
   !  WORK:    The REAL (KIND=wp) work space.
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
   !  WSS:     The sum of the squares of the weighted EPSILONS and DELTAS, the sum of the squares
   !           of the weighted DELTAS, and the sum of the squares of the weighted EPSILONS.
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
!! Set storage locations within REAL (KIND=wp) work space
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
   !  ISODR:   The variable designating whether the solution is by ODR (ISODR=TRUE) or by OLS (ISODR=FALSE).
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

   if (n .ge. 1 .and. m .ge. 1 .and. np .ge. 1 .and. nq .ge. 1 .and. &
       ldwe .ge. 1 .and. ld2we .ge. 1) then

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

pure subroutine dwght(n, m, wt, ldwt, ld2wt, t, wtt)
!! Scale matrix T using WT, i.e., compute WTT = WT*T
! Routines Called  (NONE)
! Date Written   860529   (YYMMDD)
! Revision Date  920304   (YYMMDD)

   use odrpack_kinds, only: wp, zero

   integer, intent(in) :: n
      !! The number of rows of data in `t`.
   integer, intent(in) :: m
      !! The number of columns of data in `t`.
   integer, intent(in) :: ldwt
      !! The leading dimension of array `wt`.
   integer, intent(in) :: ld2wt
      !! The second dimension of array `wt`.
   real(kind=wp), intent(in) :: wt(:, :, :)
      !! The weights.
   real(kind=wp), intent(in) :: t(:, :)
      !! The array being scaled by `wt`.
   real(kind=wp), intent(out) :: wtt(:, :)
      !! The results of weighting array `t` by `wt`. Array `wtt` can be the same as `t` only if
      !! the arrays in `wt` are upper triangular with zeros below the diagonal.

   ! Local scalars
   real(kind=wp) :: temp
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

   if (n .eq. 0 .or. m .eq. 0) return

   if (wt(1, 1, 1) .ge. zero) then
      if (ldwt .ge. n) then
         if (ld2wt .ge. m) then
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
         if (ld2wt .ge. m) then
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

subroutine dvevtr &
   (m, nq, indx, &
    v, ldv, ld2v, e, lde, ve, ldve, ld2ve, vev, ldvev, &
    wrk5)
!! Compute  V*E*trans(V) for the (INDX)TH M by NQ array in V
! Routines Called  DSOLVE
! Date Written   910613   (YYMMDD)
! Revision Date  920304   (YYMMDD)

   use odrpack_kinds, only: wp, zero

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
   real(kind=wp), intent(in) :: v(ldv, ld2v, nq)
      !! An array of `nq` by `m` matrices.
   real(kind=wp), intent(in) :: e(lde, m)
      !! The `m` by `m` matrix of the factors so `ete = (d**2 + alpha*t**2)`.
   real(kind=wp), intent(out) :: ve(ldve, ld2ve, m)
      !! The `nq` by `m` array `ve = v * inv(e)`.
   real(kind=wp), intent(out) :: vev(ldvev, nq)
      !! The `nq` by `nq` array `vev = v * inv(ete) * trans(v)`.
   real(kind=wp), intent(out) :: wrk5(m)
      !! An `m` work vector.

   ! Local scalars
   integer :: j, l1, l2

   ! External subroutines
   external :: dsolve

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

   if (nq .eq. 0 .or. m .eq. 0) return

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

subroutine dunpac(n2, v1, v2, ifix)
!! Copy the elements of `v1` into the locations of `v2` which are unfixed
! Routines Called  DCOPY
! Date Written   860529   (YYMMDD)
! Revision Date  920304   (YYMMDD)

   use odrpack_kinds, only: wp

   integer, intent(in) :: n2
      !! The number of items in `v2`.
   real(kind=wp), intent(in) :: v1(n2)
      !! The vector of the unfixed items.
   real(kind=wp), intent(out) :: v2(n2)
      !! The vector of the fixed and unfixed items into which the elements of `v1` are to be inserted.
   integer, intent(in) :: ifix(n2)
      !! The values designating whether the elements of `v2` are fixed at their input values or not.

   ! Local scalars
   integer :: i, n1

   ! External subroutines
   external :: dcopy

   ! Variable Definitions (alphabetically)
   !  I:       An indexing variable.
   !  IFIX:    The values designating whether the elements of V2 are fixed at their input values or not.
   !  N1:      The number of items in V1.
   !  N2:      The number of items in V2.
   !  V1:      The vector of the unfixed items.
   !  V2:      The vector of the fixed and unfixed items into which the elements of V1 are to be inserted.

   n1 = 0
   if (ifix(1) .ge. 0) then
      do i = 1, n2
         if (ifix(i) .ne. 0) then
            n1 = n1 + 1
            v2(i) = v1(n1)
         end if
      end do
   else
      n1 = n2
      call dcopy(n2, v1, 1, v2, 1)
   end if

end subroutine dunpac

subroutine dsolve(n, t, ldt, b, job)
!! Solve systems of the form:
!!  T * X = B  or  trans(T) * X = B
!! where T is an upper or lower triangular matrix of order N, and the solution X overwrites
!! the RHS B. (adapted from LINPACK subroutine DTRSL)
!! References:
!!  Dongarra J.J., Bunch J.R., Moler C.B., Stewart G.W., *LINPACK Users Guide*, SIAM, 1979.
! Routines Called  DAXPY,DDOT
! Date Written   920220   (YYMMDD)
! Revision Date  920619   (YYMMDD)

   use odrpack_kinds, only: wp, zero

   integer, intent(in) :: n
      !! The number of rows and columns of data in array `t`.
   real(kind=wp), intent(in) :: t(ldt, n)
      !! The upper or lower tridiagonal system.
   integer, intent(in) :: ldt
      !! The leading dimension of array `t`.
   real(kind=wp), intent(inout) :: b(n)
      !! On input: the right hand side; On exit: the solution.
   integer, intent(in) :: job
      !! What kind of system is to be solved:
      !!   1   Solve T*X=B, T lower triangular,
      !!   2   Solve T*X=B, T upper triangular,
      !!   3   Solve trans(T)*X=B, T lower triangular,
      !!   4   Solve trans(T)*X=B, T upper triangular.

   ! Local scalars
   real(kind=wp) :: temp
   integer :: j1, j, jn

   ! External functions
   real(kind=wp), external :: ddot

   ! External subroutines
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
      if (j1 .eq. 0 .and. t(j, j) .ne. zero) then
   !!! FPT - 3087 REAL or COMPLEX quantity tested for exact equality/inequality
         j1 = j
      elseif (t(j, j) .eq. zero) then
   !!! FPT - 3087 REAL or COMPLEX quantity tested for exact equality/inequality
         b(j) = zero
      end if
   end do
   if (j1 .eq. 0) return
   !  Find last nonzero diagonal entry in T
   jn = 0
   do j = n, j1, -1
      if (jn .eq. 0 .and. t(j, j) .ne. zero) then
   !!! FPT - 3087 REAL or COMPLEX quantity tested for exact equality/inequality
         jn = j
      elseif (t(j, j) .eq. zero) then
   !!! FPT - 3087 REAL or COMPLEX quantity tested for exact equality/inequality
         b(j) = zero
      end if
   end do

   if (job .eq. 1) then
      !  Solve T*X=B for T lower triangular
      b(j1) = b(j1)/t(j1, j1)
      do j = j1 + 1, jn
         temp = -b(j - 1)
         call daxpy(jn - j + 1, temp, t(j, j - 1), 1, b(j), 1)
         if (t(j, j) .ne. zero) then
   !!! FPT - 3087 REAL or COMPLEX quantity tested for exact equality/inequality
            b(j) = b(j)/t(j, j)
         else
            b(j) = zero
         end if
      end do

   elseif (job .eq. 2) then
      !  Solve T*X=B for T upper triangular.
      b(jn) = b(jn)/t(jn, jn)
      do j = jn - 1, j1, -1
         temp = -b(j + 1)
         call daxpy(j, temp, t(1, j + 1), 1, b(1), 1)
         if (t(j, j) .ne. zero) then
   !!! FPT - 3087 REAL or COMPLEX quantity tested for exact equality/inequality
            b(j) = b(j)/t(j, j)
         else
            b(j) = zero
         end if
      end do

   elseif (job .eq. 3) then
      !  Solve trans(T)*X=B for T lower triangular.
      b(jn) = b(jn)/t(jn, jn)
      do j = jn - 1, j1, -1
         b(j) = b(j) - ddot(jn - j + 1, t(j + 1, j), 1, b(j + 1), 1)
         if (t(j, j) .ne. zero) then
   !!! FPT - 3087 REAL or COMPLEX quantity tested for exact equality/inequality
            b(j) = b(j)/t(j, j)
         else
            b(j) = zero
         end if
      end do

   elseif (job .eq. 4) then
      !  Solve trans(T)*X=B for T upper triangular.
      b(j1) = b(j1)/t(j1, j1)
      do j = j1 + 1, jn
         b(j) = b(j) - ddot(j - 1, t(1, j), 1, b(1), 1)
         if (t(j, j) .ne. zero) then
   !!! FPT - 3087 REAL or COMPLEX quantity tested for exact equality/inequality
            b(j) = b(j)/t(j, j)
         else
            b(j) = zero
         end if
      end do
   end if

end subroutine dsolve

pure subroutine dsetn(n, m, x, ldx, nrow)
!! Select the row at which the derivative will be checked
! Routines Called  (None)
! Date Written   860529   (YYMMDD)
! Revision Date  920304   (YYMMDD)

   use odrpack_kinds, only: wp, zero

   integer, intent(in) :: n
      !! The number of observations.
   integer, intent(in) :: m
      !! The number of columns of data in the independent variable.
   real(kind=wp), intent(in) :: x(ldx, m)
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

   if ((nrow .ge. 1) .and. (nrow .le. n)) return

   ! Select first row of independent variables which contains no zeros
   ! if there is one, otherwise first row is used.
   nrow = 1
   do i = 1, n
      if (all(x(i, :) .ne. zero)) then
         nrow = i
         return
      end if
   end do
!--------------------------------^---------------------------------------------
!!! FPT - 3087 REAL or COMPLEX quantity tested for exact equality/inequality
!------------------------------------------------------------------------------

end subroutine dsetn

pure subroutine dscld(n, m, x, ldx, tt, ldtt)
!! Select scaling values for DELTA according to the algorithm given in the ODRPACK95 reference guide
! Routines Called  (None)
! Date Written   860529   (YYMMDD)
! Revision Date  920304   (YYMMDD)

   use odrpack_kinds, only: wp, zero, one, ten

   integer, intent(in) :: n
      !! The number of observations.
   integer, intent(in) :: m
      !! The number of columns of data in the independent variable.
   real(kind=wp), intent(in) :: x(ldx, m)
      !! The independent variable.
   integer, intent(in) :: ldx
      !! The leading dimension of array `x`.
   real(kind=wp), intent(out) :: tt(ldtt, m)
      !! The scaling values for `delta`.
   integer, intent(in) :: ldtt
      !! The leading dimension of array `tt`.

   ! Local scalars
   real(kind=wp) :: xmax, xmin
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
      if (xmax .eq. zero) then
!-------------------------^----------------------------------------------------
!!! FPT - 3087 REAL or COMPLEX quantity tested for exact equality/inequality
!------------------------------------------------------------------------------
         !  All input values of X(I,J), I=1,...,N, are zero
         do i = 1, n
            tt(i, j) = one
         end do
      else
         !  Some of the input values are nonzero
         xmin = xmax
         do i = 1, n
            if (x(i, j) .ne. zero) then
!-----------------------------------^------------------------------------------
!!! FPT - 3087 REAL or COMPLEX quantity tested for exact equality/inequality
!------------------------------------------------------------------------------
               xmin = min(xmin, abs(x(i, j)))
            end if
         end do
         bigdif = log10(xmax) - log10(xmin) .ge. one
         do i = 1, n
            if (x(i, j) .ne. zero) then
!-----------------------------------^------------------------------------------
!!! FPT - 3087 REAL or COMPLEX quantity tested for exact equality/inequality
!------------------------------------------------------------------------------
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

pure subroutine dsclb(np, beta, ssf)
!! Select scaling values for BETA according to the algorithm given in the ODRPACK95 reference guide
! Routines Called  (NONE)
! Date Written   860529   (YYMMDD)
! Revision Date  920304   (YYMMDD)

   use odrpack_kinds, only: wp, zero, one, ten

   integer, intent(in) :: np
      !! The number of function parameters.
   real(kind=wp), intent(in) :: beta(np)
      !! The function parameters.
   real(kind=wp), intent(out) :: ssf(np)
      !! The scaling values for `beta`.

   ! Local scalars
   real(kind=wp) :: bmax, bmin
   integer :: k
   logical ::bigdif

   ! Variable Definitions (alphabetically)
   !  BETA:    The function parameters.
   !  BIGDIF:  The variable designating whether there is a significant
   !           difference in the magnitudes of the nonzero elements of
   !  BETA (BIGDIF=.TRUE.) or not (BIGDIF=.FALSE.).
   !  BMAX:    The largest nonzero magnitude.
   !  BMIN:    The smallest nonzero magnitude.
   !  K:       An indexing variable.
   !  NP:      The number of function parameters.
   !  SSF:     The scaling values for BETA.

   bmax = abs(beta(1))
   do k = 2, np
      bmax = max(bmax, abs(beta(k)))
   end do

   if (bmax .eq. zero) then
!------------------------------------------------------------------------------
!!! FPT - 3087 REAL or COMPLEX quantity tested for exact equality/inequality
!------------------------------------------------------------------------------
      !  All input values of BETA are zero
      ssf(1:np) = one
   else
      !  Some of the input values are nonzero
      bmin = bmax
      do k = 1, np
         if (beta(k) .ne. zero) then
!--------------------------------^---------------------------------------------
!!! FPT - 3087 REAL or COMPLEX quantity tested for exact equality/inequality
!------------------------------------------------------------------------------
            bmin = min(bmin, abs(beta(k)))
         end if
      end do
      bigdif = log10(bmax) - log10(bmin) .ge. one
      do k = 1, np
         if (beta(k) .eq. zero) then
!--------------------------------^---------------------------------------------
!!! FPT - 3087 REAL or COMPLEX quantity tested for exact equality/inequality
!------------------------------------------------------------------------------
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

pure subroutine dscale(n, m, scl, ldscl, t, ldt, sclt, ldsclt)
!! Scale T by the inverse of SCL, I.E., compute T/SCL
! Routines Called  (NONE)
! Date Written   860529   (YYMMDD)
! Revision Date  920304   (YYMMDD)

   use odrpack_kinds, only: wp, zero, one

   integer, intent(in) :: n
      !! The number of rows of data in `t`.
   integer, intent(in) :: m
      !! The number of columns of data in `t`.
   real(kind=wp), intent(in) :: scl(ldscl, m)
      !! The scale values.
   integer, intent(in) :: ldscl
      !! The leading dimension of array `scl`.
   real(kind=wp), intent(in) :: t(ldt, m)
      !! The array to be inversely scaled by `scl`.
   integer, intent(in) :: ldt
      !! The leading dimension of array `t`.
   real(kind=wp), intent(out) :: sclt(ldsclt, m)
      !! The inversely scaled matrix.
   integer, intent(in) :: ldsclt
      !! The leading dimension of array `sclt`.

   ! Local scalars
   real(kind=wp) :: temp
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

   if (n .eq. 0 .or. m .eq. 0) return

   if (scl(1, 1) .ge. zero) then
      if (ldscl .ge. n) then
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

subroutine dpvd &
   (fcn, &
    n, m, np, nq, &
    beta, xplusd, ifixb, ifixx, ldifx, &
    nrow, j, lq, stp, &
    istop, nfev, pvd, &
    wrk1, wrk2, wrk6)
!! Compute NROW-th function value using X(NROW,J) + DELTA(NROW,J) + STP
! Routines Called  FCN
! Date Written   860529   (YYMMDD)
! Revision Date  920304   (YYMMDD)

   use odrpack_kinds, only: wp

   external :: fcn
      !! The user-supplied subroutine for evaluating the model.
   integer, intent(in) :: n
      !! The number of observations.
   integer, intent(in) :: m
      !! The number of columns of data in the independent variable.
   integer, intent(in) :: np
      !! The number of function parameters.
   integer, intent(in) :: nq
      !! The number of responses per observation.
   real(kind=wp), intent(in) :: beta(np)
      !! The function parameters.
   real(kind=wp), intent(inout) :: xplusd(n, m)
      !! The values of X + DELTA.
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
   real(kind=wp), intent(in) :: stp
      !! The step size for the finite difference derivative.
   integer, intent(out) :: istop
      !! The variable designating whether there are problems computing the function at the current `beta` and `delta`.
   integer, intent(out) :: nfev
      !! The number of function evaluations.
   real(kind=wp), intent(out) :: pvd
      !! The function value for the selected observation & response.
   real(kind=wp), intent(out) :: wrk1(n, m, nq)
      !! Work array.
   real(kind=wp), intent(out) :: wrk2(n, nq)
      !! Work array.
   real(kind=wp), intent(out) :: wrk6(n, np, nq)
      !! Work array.

   ! Local scalars
   real(kind=wp) :: xpdj

   ! Routine names used as subprogram arguments
   !  FCN:     The user-supplied subroutine for evaluating the model.
   !
   ! Variable Definitions (alphabetically)
   !  BETA:    The function parameters.
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
   if (istop .eq. 0) then
      nfev = nfev + 1
   else
      return
   end if
   xplusd(nrow, j) = xpdj

   pvd = wrk2(nrow, lq)

end subroutine dpvd

subroutine dpvb &
   (fcn, &
    n, m, np, nq, &
    beta, xplusd, ifixb, ifixx, ldifx, &
    nrow, j, lq, stp, &
    istop, nfev, pvb, &
    wrk1, wrk2, wrk6)
!! Compute the NROW-th function value using BETA(J) + STP
! Routines Called  FCN
! Date Written   860529   (YYMMDD)
! Revision Date  920304   (YYMMDD)

   use odrpack_kinds, only: wp

   external :: fcn
      !! The user-supplied subroutine for evaluating the model.
   integer, intent(in) :: n
      !! The number of observations.
   integer, intent(in) :: m
      !! The number of columns of data in the independent variable.
   integer, intent(in) :: np
      !! The number of function parameters.
   integer, intent(in) :: nq
      !! The number of responses per observation.
   real(kind=wp), intent(inout) :: beta(np)
      !! The function parameters.
   real(kind=wp), intent(in) :: xplusd(n, m)
      !! The values of X + DELTA.
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
   real(kind=wp), intent(in) :: stp
      !! The step size for the finite difference derivative.
   integer, intent(out) :: istop
      !! The variable designating whether there are problems computing the function at the current `beta` and `delta`.
   integer, intent(inout) :: nfev
      !! The number of function evaluations.
   real(kind=wp), intent(out) :: pvb
      !! The function value for the selected observation & response.
   real(kind=wp), intent(out) :: wrk1(n, m, nq)
      !! Work array.
   real(kind=wp), intent(out) :: wrk2(n, nq)
      !! Work array.
   real(kind=wp), intent(out) :: wrk6(n, np, nq)
      !! Work array.

   ! Local scalars
   real(kind=wp) :: betaj

   ! Routine names used as subprogram arguments
   !  FCN:     The user-supplied subroutine for evaluating the model.
   ! Variable Definitions (alphabetically)
   !  BETA:    The function parameters.
   !  BETAJ:   The current estimate of the jth parameter.
   !  IFIXB:   The values designating whether the elements of BETA are fixed at their input values or not.
   !  IFIXX:   The values designating whether the elements of X are fixed at their input values or not.
   !  ISTOP:   The variable designating whether there are problems computing the function at the current BETA and DELTA.
   !  J:       The index of the partial derivative being examined.
   !  LDIFX:   The leading dimension of array IFIXX.
   !  LQ:      The response currently being examined.
   !  M:       The number of columns of data in the independent variable.
   !  N:       The number of observations.
   !  NFEV:    The number of function evaluations.
   !  NP:      The number of function parameters.
   !  NQ:      The number of responses per observation.
   !  NROW:    The row number of the independent variable array at which the derivative is to be checked.
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
   if (istop .eq. 0) then
      nfev = nfev + 1
   else
      return
   end if
   beta(j) = betaj

   pvb = wrk2(nrow, lq)

end subroutine dpvb

real(kind=wp) function dppt(p, idf) result(dpptr)
!! Compute the percent point function value for the student's T distribution with IDF degrees
!! of freedom. (Adapted from DATAPAC subroutine TPPF, with modifications to facilitate
!! conversion to REAL (KIND=wp) automatically)
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

   use odrpack_kinds, only: wp, pi, zero, half, one, two, three, eight, fiftn

   real(kind=wp), intent(in) :: p
      !! The probability at which the percent point is to be evaluated. `p` must lie between
      !! 0.0 and 1.0E0_wp, exclusive.
   integer, intent(in) :: idf
      !! The (positive integer) degrees of freedom.

   ! Local scalars
   real(kind=wp) :: arg, b21, b31, b32, b33, b34, b41, b42, b43, b44, b45, &
                    b51, b52, b53, b54, b55, b56, c, con, d1, d3, d5, d7, d9, df, &
                    ppfn, s, term1, term2, term3, term4, term5, z
   integer :: ipass, maxit

   ! External functions
   real(kind=wp), external :: dppnml

   ! Data statements
   data &
      b21 &
      /4.0E0_wp/
   data &
      b31, b32, b33, b34 &
      /96.0E0_wp, 5.0E0_wp, 16.0E0_wp, 3.0E0_wp/
   data &
      b41, b42, b43, b44, b45 &
      /384.0E0_wp, 3.0E0_wp, 19.0E0_wp, 17.0E0_wp, -15.0E0_wp/
   data &
      b51, b52, b53, b54, b55, b56 &
      /9216.0E0_wp, 79.0E0_wp, 776.0E0_wp, 1482.0E0_wp, -1920.0E0_wp, &
      -945.0E0_wp/

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

   if (idf .le. 0) then
      !Treat the IDF < 1 case
      dpptr = zero

   elseif (idf .eq. 1) then
      !Treat the IDF = 1 (Cauchy) case
      arg = pi*p
      dpptr = -cos(arg)/sin(arg)

   elseif (idf .eq. 2) then
      !  Treat the IDF = 2 case
      term1 = sqrt(two)/two
      term2 = two*p - one
      term3 = sqrt(p*(one - p))
      dpptr = term1*term2/term3

   elseif (idf .ge. 3) then
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

      if (idf .eq. 3) then
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

      elseif (idf .eq. 4) then
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

      elseif (idf .eq. 5) then
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

      elseif (idf .eq. 6) then
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

real(kind=wp) function dppnml(p) result(dppnmlr)
!! Compute the percent point function value for the normal (Gaussian) distribution with mean 0
!! and standard deviation 1, and with probability density function
!!       F(X) = (1/SQRT(2*PI))*EXP(-X*X/2).
!! (Adapted from DATAPAC subroutine TPPF, with modifications to facilitate conversion to
!! REAL (KIND=wp) automatically)
! Routines Called  (None)
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

   use odrpack_kinds, only: wp, zero, half, one, two

   real(kind=wp), intent(in) :: p
      !! The probability at which the percent point is to be evaluated. `p` must lie between
      !! 0.0 and 1.0E0_wp, exclusive.

   ! Local scalars
   real(kind=wp) :: aden, anum, p0, p1, p2, p3, p4, q0, q1, q2, q3, q4, r, t

   ! Data statements
   data &
      p0, p1, p2, p3, p4 &
      /-0.322232431088E0_wp, -1.0E0_wp, -0.342242088547E0_wp, &
      -0.204231210245E-1_wp, -0.453642210148E-4_wp/
   data &
      q0, q1, q2, q3, q4 &
      /0.993484626060E-1_wp, 0.588581570495E0_wp, &
      0.531103462366E0_wp, 0.103537752850E0_wp, 0.38560700634E-2_wp/

   ! Variable Definitions (alphabetically)
   !  ADEN:    A value used in the approximation.
   !  ANUM:    A value used in the approximation.
   !  P:       The probability at which the percent point is to be evaluated.  P must be between
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

   if (p .eq. half) then
!-------------------^----------------------------------------------------------
!!! FPT - 3087 REAL or COMPLEX quantity tested for exact equality/inequality
!------------------------------------------------------------------------------
      dppnmlr = zero
   else
      r = p
      if (p .gt. half) r = one - r
      t = sqrt(-two*log(r))
      anum = ((((t*p4 + p3)*t + p2)*t + p1)*t + p0)
      aden = ((((t*q4 + q3)*t + q2)*t + q1)*t + q0)
      dppnmlr = t + (anum/aden)
      if (p .lt. half) dppnmlr = -dppnmlr
   end if

end function dppnml

subroutine dpack(n2, n1, v1, v2, ifix)
!! Select the unfixed elements of V2 and return them in V1.
! Routines Called  DCOPY
! Date Written   860529   (YYMMDD)
! Revision Date  920304   (YYMMDD)

   use odrpack_kinds, only: wp

   integer, intent(in) :: n2
      !! The number of items in V2.
   integer, intent(out) :: n1
      !! The number of items in V1.
   real(kind=wp), intent(out) :: v1(n2)
      !! The vector of the unfixed items from V2.
   real(kind=wp), intent(in) :: v2(n2)
      !! The vector of the fixed and unfixed items from which the unfixed elements are to be extracted.
   integer, intent(in) :: ifix(n2)
      !! The values designating whether the elements of V2 are fixed at their input values or not.

   ! Local scalars
   integer :: i

   ! External subroutines
   external :: dcopy

   ! Variable definitions (alphabetically)
   !  I:       An indexing variable.
   !  IFIX:    The values designating whether the elements of V2 are fixed at their input values or not.
   !  N1:      The number of items in V1.
   !  N2:      The number of items in V2.
   !  V1:      The vector of the unfixed items from V2.
   !  V2:      The vector of the fixed and unfixed items from which the unfixed elements are to be extracted.

   n1 = 0
   if (ifix(1) .ge. 0) then
      do i = 1, n2
         if (ifix(i) .ne. 0) then
            n1 = n1 + 1
            v1(n1) = v2(i)
         end if
      end do
   else
      n1 = n2
      call dcopy(n2, v2, 1, v1, 1)
   end if

end subroutine dpack

subroutine dodvcv &
   (n, m, np, nq, npp, &
    f, fjacb, fjacd, &
    wd, ldwd, ld2wd, ssf, ss, tt, ldtt, delta, &
    epsfcn, isodr, &
    vcv, sd, &
    wrk6, omega, u, qraux, jpvt, &
    s, t, irank, rcond, rss, idf, rvar, ifixb, &
    wrk1, wrk2, wrk3, wrk4, wrk5, wrk, lwrk, istopc)
!! Compute covariance matrix of estimated parameters.
! Routines Called  DPODI,DODSTP
! Date Written   901207   (YYMMDD)
! Revision Date  920619   (YYMMDD)

   use odrpack_kinds, only: wp, zero

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
   real(kind=wp), intent(in) :: f(n, nq)
      !! The (weighted) estimated values of `epsilon`.
   real(kind=wp), intent(in) :: fjacb(n, np, nq)
      !! The Jacobian with respect to `beta`.
   real(kind=wp), intent(in) :: fjacd(n, m, nq)
      !! The Jacobian with respect to `delta`.
   real(kind=wp), intent(in) :: wd(ldwd, ld2wd, m)
      !! The `delta` weights.
   integer, intent(in) :: ldwd
      !! The leading dimension of array `wd`.
   integer, intent(in) :: ld2wd
      !! The second dimension of array `wd`.
   real(kind=wp), intent(in) :: ssf(np)
      !! The scaling values used for `beta`.
   real(kind=wp), intent(in) :: ss(np)
      !! The scaling values for the unfixed `betas`.
   real(kind=wp), intent(in) :: tt(ldtt, m)
      !! The scaling values for `delta`.
   integer, intent(in) :: ldtt
      !! The leading dimension of array `tt`.
   real(kind=wp), intent(in) :: delta(n, m)
      !! The estimated errors in the explanatory variables.
   real(kind=wp), intent(in) :: epsfcn
      !! The function's precision.
   logical, intent(in) :: isodr
      !! The variable designating whether the solution is by ODR (`isodr=.TRUE.`) or by OLS (`isodr=.FALSE.`).
   real(kind=wp), intent(out) :: vcv(np, np)
      !! The covariance matrix of the estimated `betas`.
   real(kind=wp), intent(out) :: sd(np)
      !! The standard deviations of the estimated `betas`.
   real(kind=wp), intent(out) :: wrk6(n*nq, np)
      !! A work array of `(n*nq by np)` elements.
   real(kind=wp), intent(out) :: omega(nq, nq)
      !! The array defined such that `omega*trans(omega) = inv(I+fjacd*inv(e)*trans(fjacd)) = (I-fjacd*inv(p)*trans(fjacd))`.
   real(kind=wp), intent(out) :: u(np)
      !! The approximate null vector for `fjacb`.
   real(kind=wp), intent(out) :: qraux(np)
      !! The array required to recover the orthogonal part of the Q-R decomposition.
   integer, intent(out) :: jpvt(np)
      !! The pivot vector.
   real(kind=wp), intent(out) :: s(np)
      !! The step for `beta`.
   real(kind=wp), intent(out) :: t(n, m)
      !! The step for `delta`.
   integer, intent(out) :: irank
      !! The rank deficiency of the Jacobian wrt `beta`.
   real(kind=wp), intent(out) :: rcond
      !! The approximate reciprocal condition of `fjacb`.
   real(kind=wp), intent(out) :: rss
      !! The residual sum of squares.
   integer, intent(out) :: idf
      !! The degrees of freedom of the fit, equal to the number of observations with nonzero weighted derivatives minus the number of parameters being estimated.
   real(kind=wp), intent(out) :: rvar
      !! The residual variance.
   integer, intent(in) :: ifixb(np)
      !! The values designating whether the elements of `beta` are fixed at their input values or not.
   real(kind=wp), intent(out) :: wrk1(n, nq, m)
      !! A work array of `(n by nq by m)` elements.
   real(kind=wp), intent(out) :: wrk2(n, nq)
      !! A work array of `(n by nq)` elements.
   real(kind=wp), intent(out) :: wrk3(np)
      !! A work array of `(np)` elements.
   real(kind=wp), intent(out) :: wrk4(m, m)
      !! A work array of `(m by m)` elements.
   real(kind=wp), intent(out) :: wrk5(m)
      !! A work array of `(m)` elements.
   real(kind=wp), intent(out) :: wrk(lwrk)
      !! A work array of `(lwrk)` elements, equivalenced to `wrk1` and `wrk2`.
   integer, intent(in) :: lwrk
      !! The length of vector `lwrk`.
   integer, intent(out) :: istopc
      !! The variable designating whether the computations were stoped due to a numerical error within subroutine `dodstp`.

   ! Local scalars
   real(kind=wp) :: temp
   integer :: i, iunfix, j, junfix, kp
   logical :: forvcv

   ! External subroutines
   external :: dpodi, dodstp

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
               wrk1, wrk2, wrk3, wrk4, wrk5, wrk, lwrk, istopc)
   if (istopc .ne. 0) then
      return
   end if
   kp = npp - irank
   call dpodi(wrk6, n*nq, kp, wrk3, 1)

   idf = 0
   do i = 1, n
      if (any(fjacb(i, :, :) .ne. zero)) then
         idf = idf + 1
         cycle
      end if
      if (isodr) then
         if (any(fjacd(i, :, :) .ne. zero)) then
            idf = idf + 1
            cycle
         end if
      end if
   end do
!---------------------------------------------^--------------------------------
!!! FPT - 3087 REAL or COMPLEX quantity tested for exact equality/inequality
!------------------------------------------------------------------------------

   if (idf .gt. kp) then
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
   if (np .gt. npp) then
      junfix = npp
      do j = np, 1, -1
         if (ifixb(j) .eq. 0) then
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
         if (jpvt(i) .gt. jpvt(j)) then
            vcv(jpvt(i), jpvt(j)) = wrk6(i, j)
         else
            vcv(jpvt(j), jpvt(i)) = wrk6(i, j)
         end if
      end do
   end do
   if (np .gt. npp) then
      iunfix = npp
      do i = np, 1, -1
         if (ifixb(i) .eq. 0) then
            do j = i, 1, -1
               vcv(i, j) = zero
            end do
         else
            junfix = npp
            do j = np, 1, -1
               if (ifixb(j) .eq. 0) then
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
      if (ssf(1) .gt. zero) then
         sd(i) = sd(i)/ssf(i)
      else
         sd(i) = sd(i)/abs(ssf(1))
      end if
      do j = 1, np
         if (ssf(1) .gt. zero) then
            vcv(i, j) = vcv(i, j)/(ssf(i)*ssf(j))
         else
            vcv(i, j) = vcv(i, j)/(ssf(1)*ssf(1))
         end if
      end do
   end do

end subroutine dodvcv

subroutine dodstp &
   (n, m, np, nq, npp, &
    f, fjacb, fjacd, &
    wd, ldwd, ld2wd, ss, tt, ldtt, delta, &
    alpha, epsfcn, isodr, &
    tfjacb, omega, u, qraux, kpvt, &
    s, t, phi, irank, rcond, forvcv, &
    wrk1, wrk2, wrk3, wrk4, wrk5, wrk, lwrk, istopc)
!! Compute locally constrained steps S and T, and PHI(ALPHA).
! Routines Called  IDAMAX,DCHEX,DESUBI,DFCTR,DNRM2,DQRDC,DQRSL,DROT,
!                  DROTG,DSOLVE,DTRCO,DTRSL,DVEVTR,DWGHT
! Date Written   860529   (YYMMDD)
! Revision Date  920619   (YYMMDD)

   use odrpack_kinds, only: wp, zero, one
   use odrpack, only: tempret

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
   real(kind=wp), intent(in) :: f(n, nq)
      !! The (weighted) estimated values of `epsilon`.
   real(kind=wp), intent(in) :: fjacb(n, np, nq)
      !! The Jacobian with respect to `beta`.
   real(kind=wp), intent(in) :: fjacd(n, m, nq)
      !! The Jacobian with respect to `delta`.
   real(kind=wp), intent(in) :: wd(ldwd, ld2wd, m)
      !! The (squared) `delta` weights.
   integer, intent(in) :: ldwd
      !! The leading dimension of array `wd`.
   integer, intent(in) :: ld2wd
      !! The second dimension of array `wd`.
   real(kind=wp), intent(in) :: ss(np)
      !! The scaling values for the unfixed `betas`.
   real(kind=wp), intent(in) :: tt(ldtt, m)
      !! The scaling values for `delta`.
   integer, intent(in) :: ldtt
      !! The leading dimension of array `tt`.
   real(kind=wp), intent(in) :: delta(n, m)
      !! The estimated errors in the explanatory variables.
   real(kind=wp), intent(in) :: alpha
      !! The Levenberg-Marquardt parameter.
   real(kind=wp), intent(in) :: epsfcn
      !! The function's precision.
   logical, intent(in) :: isodr
      !! The variable designating whether the solution is by ODR (`isodr`=TRUE) or by OLS (`isodr`=FALSE).
   real(kind=wp), intent(out) :: tfjacb(n, nq, np)
      !! The array `omega`*`fjacb`.
   real(kind=wp), intent(out) :: omega(nq, nq)
      !! The array defined S.T.
      !! `omega`*trans(`omega`) = inv(I+`fjacd`*inv(E)*trans(`fjacd`))
      !! = (I-`fjacd`*inv(P)*trans(`fjacd`))
      !! where E = D**2 + `alpha`*`tt`**2
      !! P = trans(`fjacd`)*`fjacd` + D**2 + `alpha`*`tt`**2
   real(kind=wp), intent(out) :: u(np)
      !! The approximate null vector for `tfjacb`.
   real(kind=wp), intent(out) :: qraux(np)
      !! The array required to recover the orthogonal part of the Q-R decomposition.
   integer, intent(out) :: kpvt(np)
      !! The pivot vector.
   real(kind=wp), intent(out) :: s(np)
      !! The step for `beta`.
   real(kind=wp), intent(out) :: t(n, m)
      !! The step for `delta`.
   real(kind=wp), intent(out) :: phi
      !! The difference between the norm of the scaled step and the trust region diameter.
   integer, intent(out) :: irank
      !! The rank deficiency of the Jacobian wrt `beta`.
   real(kind=wp), intent(out) :: rcond
      !! The approximate reciprocal condition number of `tfjacb`.
   logical, intent(in) :: forvcv
      !! The variable designating whether this subroutine was called to set up for the covariance matrix computations (`forvcv`=TRUE) or not (`forvcv`=FALSE).
   real(kind=wp), intent(out) :: wrk1(n, nq, m)
      !! A work array of (`n` by `nq` by `m`) elements.
   real(kind=wp), intent(out) :: wrk2(n, nq)
      !! A work array of (`n` by `nq`) elements.
   real(kind=wp), intent(out) :: wrk3(np)
      !! A work array of (`np`) elements.
   real(kind=wp), intent(out) :: wrk4(m, m)
      !! A work array of (`m` by `m`) elements.
   real(kind=wp), intent(out) :: wrk5(m)
      !! A work array of (`m`) elements.
   real(kind=wp), intent(out) :: wrk(lwrk)
      !! A work array of (`lwrk`) elements, equivalenced to `wrk1` and `wrk2`.
   integer, intent(in) :: lwrk
      !! The length of vector `wrk`.
   integer, intent(out) :: istopc
      !! The variable designating whether the computations were stopped due to a numerical error within subroutine `dodstp`.

   ! Local scalars
   real(kind=wp) :: co, si, temp
   integer :: i, imax, inf, ipvt, j, k, k1, k2, kp, l
   logical :: elim

   ! Local arrays
   real(kind=wp) :: dum(2)

   ! External functions
   real(kind=wp), external :: dnrm2
   integer, external :: idamax

   ! External subroutines
   external :: dchex, desubi, dfctr, dqrdc, dqrsl, drot, drotg, dsolve, dtrco, dtrsl, dvevtr

   ! Interface blocks
   interface
      subroutine dwght &
         (n, m, wt, ldwt, ld2wt, t, wtt)
         use odrpack_kinds, only: wp
         integer &
            ldwt, ld2wt, m, n
         real(kind=wp) &
            t(:, :), wt(:, :, :), wtt(:, :)
      end subroutine
   end interface

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
   !                error within subroutine DODSTP.
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
   if (alpha .eq. zero) then
!-----------------------^------------------------------------------------------
!!! FPT - 3087 REAL or COMPLEX quantity tested for exact equality/inequality
!------------------------------------------------------------------------------
      kp = npp
      do k = 1, np
         kpvt(k) = k
      end do
   else
      if (npp .ge. 1) then
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
         if (inf .ne. 0) then
            istopc = 60000
            return
         end if
         !  Compute OMEGA, such that
         !             trans(OMEGA)*OMEGA = I+FJACD*inv(E)*trans(FJACD)
         !             inv(trans(OMEGA)*OMEGA) = I-FJACD*inv(P)*trans(FJACD)
         call dvevtr(m, nq, i, fjacd, n, m, wrk4, m, wrk1, n, nq, omega, nq, wrk5)
         do l = 1, nq
            omega(l, l) = one + omega(l, l)
         end do
         call dfctr(.false., omega, nq, nq, inf)
         if (inf .ne. 0) then
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
               if (ss(1) .gt. zero) then
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
               if (ss(1) .gt. zero) then
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
   if (alpha .eq. zero) then
!-----------------------^------------------------------------------------------
!!! FPT - 3087 REAL or COMPLEX quantity tested for exact equality/inequality
!------------------------------------------------------------------------------
      ipvt = 1
      do k = 1, np
         kpvt(k) = 0
      end do
   else
      ipvt = 0
   end if

   call dqrdc(tfjacb, n*nq, n*nq, kp, qraux, kpvt, wrk3, ipvt)
   call dqrsl(tfjacb, n*nq, n*nq, kp, &
              qraux, wrk2, dum, wrk2, dum, dum, dum, 1000, inf)
   if (inf .ne. 0) then
      istopc = 60000
      return
   end if

   ! Eliminate alpha part using givens rotations
   if (alpha .ne. zero) then
!-----------------------^------------------------------------------------------
!!! FPT - 3087 REAL or COMPLEX quantity tested for exact equality/inequality
!------------------------------------------------------------------------------
      s(1:npp) = zero
      do k1 = 1, kp
         wrk3(1:kp) = zero
         wrk3(k1) = sqrt(alpha)
         do k2 = k1, kp
            call drotg(tfjacb(k2, 1, k2), wrk3(k2), co, si)
            if (kp - k2 .ge. 1) then
               call drot(kp - k2, tfjacb(k2, 1, k2 + 1), n*nq, &
                         wrk3(k2 + 1), 1, co, si)
            end if
            temp = co*wrk2(k2, 1) + si*s(kpvt(k1))
            s(kpvt(k1)) = -si*wrk2(k2, 1) + co*s(kpvt(k1))
            wrk2(k2, 1) = temp
         end do
      end do
   end if

   ! Compute solution - eliminate variables if necessary
   if (npp .ge. 1) then
      if (alpha .eq. zero) then
!--------------------------^---------------------------------------------------
!!! FPT - 3087 REAL or COMPLEX quantity tested for exact equality/inequality
!------------------------------------------------------------------------------
         kp = npp
         elim = .true.
         do while (elim .and. kp .ge. 1)
            ! Estimate RCOND - U will contain approx null vector
            call dtrco(tfjacb, n*nq, kp, rcond, u, 1)
            if (rcond .le. epsfcn) then
               elim = .true.
               imax = idamax(kp, u, 1)
               ! IMAX is the column to remove - use DCHEX and fix KPVT
               if (imax .ne. kp) then
                  call dchex(tfjacb, n*nq, kp, imax, kp, wrk2, n*nq, 1, &
                             qraux, wrk3, 2)
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
   if (npp .ge. 1) then
      do i = kp + 1, npp
         wrk2(i, 1) = zero
      end do
      if (kp .ge. 1) then
         call dtrsl(tfjacb, n*nq, kp, wrk2, 01, inf)
         if (inf .ne. 0) then
            istopc = 60000
            return
         end if
      end if
      do i = 1, npp
         if (ss(1) .gt. zero) then
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
         if (inf .ne. 0) then
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
