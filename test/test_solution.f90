module test_solution_m

   use odrpack_kinds, only: wp
   implicit none

   integer, parameter :: maxn = 50, maxm = 3, maxq = 2, maxnp = 10
   real(wp), allocatable :: lower(:), upper(:)
   integer :: setno

contains

   impure subroutine test_odr(ntest, tstfac, passed, lunrpt, lunerr, lunsum)
   !! Test solutions delivered with [[odr]].

      use odrpack, only: odr, workspace_dimensions

      integer, intent(in) :: ntest
      real(wp), intent(in) :: tstfac
      logical, intent(inout) :: passed
      integer, intent(in) :: lunrpt
      integer, intent(in) :: lunerr
      integer, intent(in) :: lunsum

      ! Parameters
      real(wp), parameter :: base = radix(1.0_wp)
      integer, parameter :: ntests = 23

      ! Local scalars
      integer :: i, info, iprint, itest, job, l, ldstpd, lun, m, maxit, msg, n, ndigit, &
                 np, q, ldwd, ld2wd, ldwe, ld2we, lrwork, liwork, wssi, wssdei, wssepi
      real(wp) :: bnrm, epsmac, ewrt, ewrt2, hundrd, one, p01, p2, partol, sstol, &
                  taufac, three, tsttol, two, wss, wssdel, wsseps, zero
      logical :: failed, fails, isodr, short, wflat
      character(len=80) :: title

      ! Local arrays
      integer :: idpymp(ntests)
      integer, allocatable :: ifixb(:), ifixx(:, :), iwork(:)
      real(wp) :: dpymp(2, ntests)
      real(wp), allocatable :: delta(:, :), beta(:), x(:, :), y(:, :), sclb(:), scld(:, :), &
                               stpb(:), stpd(:, :), we(:, :, :), wd(:, :, :), beta_last(:), &
                               rwork(:)

      ! Data statements
      data &
         zero, p01, p2, one, two, three, hundrd &
         /0.0E0_wp, 0.01E0_wp, 0.2E0_wp, 1.0E0_wp, 2.0E0_wp, 3.0E0_wp, &
         100.0E0_wp/
      data &
         (dpymp(i, 1), i=1, 2) &
         /2.762733195780256808978449342964E+04_wp, &
         7.532639569022918943695104672512E-04_wp/
      data &
         (dpymp(i, 2), i=1, 2) &
         /2.762732630143673024399942947263E+04_wp, &
         7.538467722687131506874279314940E-04_wp/
      data &
         (dpymp(i, 3), i=1, 2) &
         /1.069944100000000027940905194068E+09_wp, &
         1.212808593256056359629660672046E-05_wp/
      data &
         (dpymp(i, 4), i=1, 2) &
         /1.069944100000000026623461142867E+09_wp, &
         5.452084633790606017572015067556E-07_wp/
      data &
         (dpymp(i, 5), i=1, 2) &
         /1.426988156377258617521571734503E+00_wp, &
         1.084728687127432219753903919409E+00_wp/
      data &
         (dpymp(i, 6), i=1, 2) &
         /4.261321829513978871872508874025E+00_wp, &
         1.477967210398420733565424329280E-02_wp/
      data &
         (dpymp(i, 7), i=1, 2) &
         /4.261272307142888464011486769858E+00_wp, &
         1.477966125465374336804138554559E-02_wp/
      data &
         (dpymp(i, 8), i=1, 2) &
         /4.371487317909745009110272283622E+01_wp, &
         1.144419474408286067112233592550E-03_wp/
      data &
         (dpymp(i, 9), i=1, 2) &
         /3.099048849376848610380977303924E+00_wp, &
         8.824708863783850023783338218501E-02_wp/
      data &
         (dpymp(i, 10), i=1, 2) &
         /9.469917836739932584221023234527E+00_wp, &
         4.205389215588104651198536809880E-01_wp/
      data &
         (dpymp(i, 11), i=1, 2) &
         /3.950949253027682207109233363651E+01_wp, &
         6.651838750834910819636881506915E+01_wp/
      data &
         (dpymp(i, 12), i=1, 2) &
         /3.950949253027682207109233363651E+01_wp, &
         6.651838750834910819636881506915E+01_wp/
      data &
         (dpymp(i, 13), i=1, 2) &
         /1.414213562373095000000000000000E+00_wp, &
         5.250825926608277346013642256883E-26_wp/
      data &
         (dpymp(i, 14), i=1, 2) &
         /1.414213562373095000000000000000E+00_wp, &
         8.159081600696301507018019048968E-26_wp/
      data &
         (dpymp(i, 15), i=1, 2) &
         /1.486588477064952451556223422813E+00_wp, &
         1.841690442255357083922717720270E+03_wp/
      data &
         (dpymp(i, 16), i=1, 2) &
         /2.001224625073357401561224833131E+02_wp, &
         0.000000000000000000000000000000E+00_wp/
      data &
         (dpymp(i, 17), i=1, 2) &
         /2.000099997500125000000000000000E+02_wp, &
         0.000000000000000000000000000000E+00_wp/
      data &
         (dpymp(i, 18), i=1, 2) &
         /1.414213562373095000000000000000E+00_wp, &
         5.816277809383742531415846947805E-26_wp/
      data &
         (dpymp(i, 19), i=1, 2) &
         /2.000624902374255782433465356007E+02_wp, &
         4.568236947482152283374593507328E+30_wp/
      data &
         (dpymp(i, 20), i=1, 2) &
         /2.000624902374255782433465356007E+02_wp, &
         1.848525209410256939008831977844E+05_wp/
      data &
         (dpymp(i, 21), i=1, 2) &
         /2.000624902374255782433465356007E+02_wp, &
         1.848525209410256939008831977844E+05_wp/
      data &
         (dpymp(i, 22), i=1, 2) &
         /2.731300056749532689792659000000E+00_wp, &
         3.378975642596100806258619000000E+05_wp/
      data &
         (dpymp(i, 23), i=1, 2) &
         /2.675757304209387399396291584708E+00_wp, &
         5.174484505019630309341494012187E-02_wp/
      data &
         (idpymp(i), i=1, 23) &
         /1, 1, 3, 1, 1, 4, 1, 1, 2, 1, 1023, 40100, 2, 2, 3, 90100, 91000, 2, 90010, &
         90020, 90010, 21, 1/

      ! Routine names used as subprogram arguments
      !  odrxf:  The user-supplied routine for evaluating the model.

      ! Variable definitions (alphabetically)
      !  BASE:    The base of floating point numbers on the current machine
      !  BETA:    The function parameters.
      !  BNRM:    The norm of BETA.
      !  DELTA:   The error in the X data.
      !  DPYMP:   The floating point results from a cray YMP using REAL (wp).
      !  EPSMAC:  The value of machine precision.
      !  EWRT:    A temporary variable for the denominator of the relative error
      !           calculations (error with respect to).
      !  EWRT2:   A temporary variable for the denominator of the relative error
      !           calculations (error with respect to).
      !  FAILED:  The variable designating whether the results of all of the
      !           demonstration runs agreed with those from the cray YMP
      !           using REAL (wp) (FAILED=FALSE) or whether some of
      !           the tests disagreed (FAILED=TRUE).
      !  FAILS:   The variable designating whether the results of an
      !           individual demonstration run agreed with those from the
      !           cray YMP using REAL (wp) (FAILS=FALSE) or disagree (FAILS=TRUE).
      !  I:       An index variable.
      !  IDPYMP:  The integer results from a cray YMP using REAL (wp).
      !  IFIXB:   The values designating whether the elements of BETA are
      !           fixed at their input values or not.
      !  IFIXX:   The values designating whether the elements of DELTA are
      !           fixed at their input values or not.
      !  INFO:    The variable designating why the computations stopped.
      !  IPRINT:  The print control variable.
      !  ISODR:   The variable designating whether the solution is by odr
      !           (ISODR=TRUE) or by ols (ISODR=FALSE).
      !  ITEST:   The number of the current test being run.
      !  IWORK:   The integer work space.
      !  J:       An index variable.
      !  JOB:     The variable controlling problem initialization and computational method.
      !  LDWD:    The leading dimension of array WD.
      !  LDWE:    The leading dimension of array WE.
      !  LD2WD:   The second dimension of array WD.
      !  LD2WE:   The second dimension of array WE.
      !  LIWORK:  The length of vector IWORK.
      !  LUN:     The logical unit number currently being used.
      !  LUNERR:  The logical unit number used for error messages.
      !  LUNRPT:  The logical unit number used for computation reports.
      !  LUNSUM:  The logical unit number used for a summary report.
      !  LRWORK:  The length of vector RWORK.
      !  M:       The number of columns of data in the explanatory variable.
      !  MAXIT:   The maximum number of iterations allowed.
      !  MSG:     The variable designating which message is to be printed as
      !           a result of the comparison with the cray YMP or x86 (Linux) results.
      !  N:       The number of observations.
      !  NDIGIT:  The number of accurate digits in the function results, as supplied by the user.
      !  NP:      The number of function parameters.
      !  NTEST:   The number of tests to be run.
      !  NTESTS:  The number of different tests available.
      !  PASSED:  The variable designating whether the results of all of the
      !           demonstration runs agreed with those from the cray YMP
      !           using REAL (wp) (PASSED=TRUE), or whether some of
      !           the results disagreed (PASSED=FALSE).
      !  PARTOL:  The parameter convergence stopping criteria.
      !  SCLB:    The scaling values for BETA.
      !  SCLD:    The scaling values for DELTA.
      !  SETNO:   The number of the data set being analyzed.
      !  SHORT:   The variable designating whether ODRPACK95 is invoked by the
      !           short-call (SHORT=.TRUE.) or the long-call (SHORT=.FALSE.).
      !  SSTOL:   The sum-of-squares convergence stopping tolerance.
      !  TAUFAC:  The factor used to compute the initial trust region diameter.
      !  TITLE:   The reference for the data set being analyzed.
      !  TSTFAC:  The user-supplied factor for scaling the test tolerances
      !           used to check for agreement between computed results and
      !           results obtained using REAL (wp) version on cray YMP.
      !  TSTTOL:  The test tolerance used in checking computed values for
      !           purposes of determining proper installation.
      !  WD:      The DELTA weights.
      !  WE:      The EPSILON weights.
      !  RWORK:   The REAL work space.
      !  WRK:     The REAL work space for computing test results.
      !  WSS:     The sum of the squared weighted errors.
      !  WSSDEL:  The sum of the squared weighted errors in X.
      !  WSSEPS:  The sum of the squared weighted errors in Y.
      !  X:       The explanatory variable.
      !  Y:       The response variable.

      ! Initialize test tolerance
      if (tstfac > one) then
         tsttol = tstfac
      else
         tsttol = one
      end if

      ! Initialize machine precision
      epsmac = base**(1 - digits(base))

      ! Initialize miscellaneous variables used in the exercise procedure
      failed = .false.

      ! Begin exercising ODRPACK95
      test: do itest = 1, ntest

         ! Set control values to invoke default values
         ndigit = -1
         taufac = -one
         sstol = -one
         partol = -one
         maxit = -1
         wflat = .true.
         short = .false.
         iprint = 2112
         ! iprint = 6616

         if (itest == 1) then

            ! Test simple ODR problem with analytic derivatives.

            lun = lunrpt
            do i = 1, 2
               write (lun, 1001) itest
               write (lun, 1010)
               lun = lunsum
            end do

            setno = 5
            call set_inputs( &
               wflat, &
               title, n, m, np, q, ldwd, ld2wd, ldwe, ld2we, ldstpd, &
               x, y, beta, lower, upper, delta, ifixb, ifixx, wd, we, sclb, scld, stpb, stpd)
            job = 00020
            short = .true.

         elseif (itest == 2) then

            ! Test simple OLS problem with forward difference derivatives.

            lun = lunrpt
            do i = 1, 2
               write (lun, 1001) itest
               write (lun, 1020)
               lun = lunsum
            end do

            setno = 5
            call set_inputs( &
               wflat, &
               title, n, m, np, q, ldwd, ld2wd, ldwe, ld2we, ldstpd, &
               x, y, beta, lower, upper, delta, ifixb, ifixx, wd, we, sclb, scld, stpb, stpd)
            job = 00002
            short = .true.

         elseif (itest == 3) then

            ! Test parameter fixing capabilities for poorly scaled OLS problem
            ! with analytic derivatives. (derivative checking turned off.)

            lun = lunrpt
            do i = 1, 2
               write (lun, 1001) itest
               write (lun, 1030)
               lun = lunsum
            end do

            setno = 3
            call set_inputs( &
               wflat, &
               title, n, m, np, q, ldwd, ld2wd, ldwe, ld2we, ldstpd, &
               x, y, beta, lower, upper, delta, ifixb, ifixx, wd, we, sclb, scld, stpb, stpd)
            ifixb = [1, 1, 1, 0, 1, 0, 0, 0, 0]
            job = 00042

         elseif (itest == 4) then

            ! Test weighting capabilities for odr problem with analytic derivatives.
            ! Also shows solution of poorly scaled odr problem.
            ! (derivative checking turned off.)
            ! N.B., this run continues from where test 3 left off.

            lun = lunrpt
            do i = 1, 2
               write (lun, 1001) itest
               write (lun, 1040)
               lun = lunsum
            end do

            setno = 3
            wflat = .false.
            beta_last = beta
            call set_inputs( &
               wflat, &
               title, n, m, np, q, ldwd, ld2wd, ldwe, ld2we, ldstpd, &
               x, y, beta, lower, upper, delta, ifixb, ifixx, wd, we, sclb, scld, stpb, stpd)
            beta = beta_last
            do i = 1, n
               wd(i, 1, 1) = (p01/abs(x(i, 1)))**2
               we(i, 1, 1) = one
            end do
            we(28, 1, 1) = zero
            ifixb = [1, 1, 1, 0, 1, 1, 1, 0, 0]
            job = 00030
            iprint = 2232

         elseif (itest == 5) then

            ! Test DELTA initialization capabilities and user-supplied scaling
            ! and use of istop to restrict parameter values for ODR problem with analytic
            ! derivatives.

            lun = lunrpt
            do i = 1, 2
               write (lun, 1001) itest
               write (lun, 1050)
               lun = lunsum
            end do

            setno = 1
            call set_inputs( &
               wflat, &
               title, n, m, np, q, ldwd, ld2wd, ldwe, ld2we, ldstpd, &
               x, y, beta, lower, upper, delta, ifixb, ifixx, wd, we, sclb, scld, stpb, stpd)
            job = 01020
            scld(:, 1) = two
            sclb = [p2, one]
            do i = 20, 21
               delta(i, 1) = beta(1)/y(i, 1) + beta(2) - x(i, 1)
            end do

         elseif (itest == 6) then

            ! Test stiff stopping conditions for unscaled ODR problem  with analytic
            ! derivatives.

            lun = lunrpt
            do i = 1, 2
               write (lun, 1001) itest
               write (lun, 1060)
               lun = lunsum
            end do

            setno = 4
            call set_inputs( &
               wflat, &
               title, n, m, np, q, ldwd, ld2wd, ldwe, ld2we, ldstpd, &
               x, y, beta, lower, upper, delta, ifixb, ifixx, wd, we, sclb, scld, stpb, stpd)
            job = 00020
            sstol = hundrd*epsmac
            partol = epsmac
            maxit = 2

         elseif (itest == 7) then

            ! Test restart for unscaled ODR problem with analytic derivatives.

            lun = lunrpt
            do i = 1, 2
               write (lun, 1001) itest
               write (lun, 1070)
               lun = lunsum
            end do

            setno = 4
            job = 20220
            sstol = hundrd*epsmac
            partol = epsmac
            maxit = 50

         elseif (itest == 8) then

            ! Test use of TAUFAC to restrict first step
            ! for ODR problem with central difference derivatives.

            lun = lunrpt
            do i = 1, 2
               write (lun, 1001) itest
               write (lun, 1080)
               lun = lunsum
            end do

            setno = 6
            call set_inputs( &
               wflat, &
               title, n, m, np, q, ldwd, ld2wd, ldwe, ld2we, ldstpd, &
               x, y, beta, lower, upper, delta, ifixb, ifixx, wd, we, sclb, scld, stpb, stpd)
            job = 00210
            taufac = p01

         elseif (itest == 9) then

            ! Test implicit ODR problem
            ! with forward finite difference derivatives
            ! and covariance matrix constructed with recomputed derivatives.

            lun = lunrpt
            do i = 1, 2
               write (lun, 1001) itest
               write (lun, 1090)
               lun = lunsum
            end do

            setno = 7
            call set_inputs( &
               wflat, &
               title, n, m, np, q, ldwd, ld2wd, ldwe, ld2we, ldstpd, &
               x, y, beta, lower, upper, delta, ifixb, ifixx, wd, we, sclb, scld, stpb, stpd)
            job = 00001
            partol = epsmac**(one/three)

         elseif (itest == 10) then

            ! Test multiresponse ODR problem
            ! with central difference derivatives ,
            ! DELTA initialized to nonzero values,
            ! variable fixing, and weighting.

            lun = lunrpt
            do i = 1, 2
               write (lun, 1001) itest
               write (lun, 1100)
               lun = lunsum
            end do
            setno = 8
            wflat = .false.
            call set_inputs( &
               wflat, &
               title, n, m, np, q, ldwd, ld2wd, ldwe, ld2we, ldstpd, &
               x, y, beta, lower, upper, delta, ifixb, ifixx, wd, we, sclb, scld, stpb, stpd)

            do i = 1, n
               ! Initialize DELTA, and specify first decade of frequencies as fixed
               if (x(i, 1) < 100.0E0_wp) then
                  delta(i, 1) = 0.0E0_wp
                  ifixx(i, 1) = 0
               elseif (x(i, 1) <= 150.0E0_wp) then
                  delta(i, 1) = 0.0E0_wp
                  ifixx(i, 1) = 1
               elseif (x(i, 1) <= 1000.0E0_wp) then
                  delta(i, 1) = 25.0E0_wp
                  ifixx(i, 1) = 1
               elseif (x(i, 1) <= 10000.0E0_wp) then
                  delta(i, 1) = 560.0E0_wp
                  ifixx(i, 1) = 1
               elseif (x(i, 1) <= 100000.0E0_wp) then
                  delta(i, 1) = 9500.0E0_wp
                  ifixx(i, 1) = 1
               else
                  delta(i, 1) = 144000.0E0_wp
                  ifixx(i, 1) = 1
               end if

               ! Set weights
               if (x(i, 1) == 100.0E0_wp .or. x(i, 1) == 150.0E0_wp) then
                  we(i, 1, 1) = 0.0E0_wp
                  we(i, 1, 2) = 0.0E0_wp
                  we(i, 2, 1) = 0.0E0_wp
                  we(i, 2, 2) = 0.0E0_wp
               else
                  we(i, 1, 1) = 559.6E0_wp
                  we(i, 1, 2) = -1634.0E0_wp
                  we(i, 2, 1) = -1634.0E0_wp
                  we(i, 2, 2) = 8397.0E0_wp
               end if
               wd(i, 1, 1) = (1.0E-4_wp)/(x(i, 1)**2)
            end do
            job = 00210

         elseif (itest == 11) then

            ! Test detection of incorrect derivatives

            lun = lunrpt
            do i = 1, 2
               write (lun, 1001) itest
               write (lun, 1110)
               lun = lunsum
            end do

            setno = 6
            call set_inputs( &
               wflat, &
               title, n, m, np, q, ldwd, ld2wd, ldwe, ld2we, ldstpd, &
               x, y, beta, lower, upper, delta, ifixb, ifixx, wd, we, sclb, scld, stpb, stpd)
            job = 00022

         elseif (itest == 12) then

            ! Test detection of incorrect derivatives

            lun = lunrpt
            do i = 1, 2
               write (lun, 1001) itest
               write (lun, 1120)
               lun = lunsum
            end do

            setno = 6
            call set_inputs( &
               wflat, &
               title, n, m, np, q, ldwd, ld2wd, ldwe, ld2we, ldstpd, &
               x, y, beta, lower, upper, delta, ifixb, ifixx, wd, we, sclb, scld, stpb, stpd)
            job = 00020

         elseif (itest == 13) then

            ! Test bounded ODR problem where
            ! parameters start on bound, move away, hit bound, move away, find minimum.

            lun = lunrpt
            do i = 1, 2
               write (lun, 1001) itest
               write (lun, 1010)
               lun = lunsum
            end do

            setno = 9
            call set_inputs( &
               wflat, &
               title, n, m, np, q, ldwd, ld2wd, ldwe, ld2we, ldstpd, &
               x, y, beta, lower, upper, delta, ifixb, ifixx, wd, we, sclb, scld, stpb, stpd)
            job = 00000
            maxit = 100
            beta = [200.0_wp, 5.0_wp]
            lower = [0.1_wp, 0.0_wp]
            upper = [200.0_wp, 5.0_wp]

         elseif (itest == 14) then

            ! Test bounded ODR problem where bounds are never hit.

            lun = lunrpt
            do i = 1, 2
               write (lun, 1001) itest
               write (lun, 1010)
               lun = lunsum
            end do

            setno = 9
            call set_inputs( &
               wflat, &
               title, n, m, np, q, ldwd, ld2wd, ldwe, ld2we, ldstpd, &
               x, y, beta, lower, upper, delta, ifixb, ifixx, wd, we, sclb, scld, stpb, stpd)
            job = 00000
            maxit = 100
            lower = [0.0_wp, 0.0_wp]
            upper = [400.0_wp, 6.0_wp]

         elseif (itest == 15) then

            ! Test bounded ODR problem where minimum is on boundary.

            lun = lunrpt
            do i = 1, 2
               write (lun, 1001) itest
               write (lun, 1010)
               lun = lunsum
            end do

            setno = 9
            call set_inputs( &
               wflat, &
               title, n, m, np, q, ldwd, ld2wd, ldwe, ld2we, ldstpd, &
               x, y, beta, lower, upper, delta, ifixb, ifixx, wd, we, sclb, scld, stpb, stpd)
            job = 00000
            maxit = 1000
            beta = [200.0_wp, 3.0_wp]
            lower = [1.1_wp, 0.0_wp]
            upper = [400.0_wp, 6.0_wp]
            tsttol = 500.0_wp

         elseif (itest == 16) then

            ! Test bounded ODR problem where initial BETA is outside bounds.

            lun = lunrpt
            do i = 1, 2
               write (lun, 1001) itest
               write (lun, 1010)
               lun = lunsum
            end do

            setno = 9
            call set_inputs( &
               wflat, &
               title, n, m, np, q, ldwd, ld2wd, ldwe, ld2we, ldstpd, &
               x, y, beta, lower, upper, delta, ifixb, ifixx, wd, we, sclb, scld, stpb, stpd)
            job = 00000
            maxit = 1000
            beta = [200.0_wp, 7.0_wp]
            lower = [1.1_wp, 0.0_wp]
            upper = [200.0_wp, 5.0_wp]

         elseif (itest == 17) then

            ! Test bounded ODR problem where bounds are ill defined.

            lun = lunrpt
            do i = 1, 2
               write (lun, 1001) itest
               write (lun, 1010)
               lun = lunsum
            end do

            setno = 9
            call set_inputs( &
               wflat, &
               title, n, m, np, q, ldwd, ld2wd, ldwe, ld2we, ldstpd, &
               x, y, beta, lower, upper, delta, ifixb, ifixx, wd, we, sclb, scld, stpb, stpd)
            job = 00000
            maxit = 1000
            beta = [200.0_wp, 2.0_wp]
            lower = [10.0_wp, 0.0_wp]
            upper = [2.0_wp, 5.0_wp]

         elseif (itest == 18) then

            ! Test bounded ODR problem using centered differences where
            ! parameters start on bound, move away, hit bound, move away, find minimum.

            lun = lunrpt
            do i = 1, 2
               write (lun, 1001) itest
               write (lun, 1010)
               lun = lunsum
            end do

            setno = 9
            call set_inputs( &
               wflat, &
               title, n, m, np, q, ldwd, ld2wd, ldwe, ld2we, ldstpd, &
               x, y, beta, lower, upper, delta, ifixb, ifixx, wd, we, sclb, scld, stpb, stpd)
            job = 00010
            maxit = 100
            beta = [200.0_wp, 5.0_wp]
            lower = [0.1_wp, 0.0_wp]
            upper = [200.0_wp, 5.0_wp]

         elseif (itest == 19) then

            ! Test bounded ODR problem when bounds are too small.
            ! Parameters start on bound.

            lun = lunrpt
            do i = 1, 2
               write (lun, 1001) itest
               write (lun, 1010)
               lun = lunsum
            end do

            setno = 9
            call set_inputs( &
               wflat, &
               title, n, m, np, q, ldwd, ld2wd, ldwe, ld2we, ldstpd, &
               x, y, beta, lower, upper, delta, ifixb, ifixx, wd, we, sclb, scld, stpb, stpd)
            job = 00010
            maxit = 100
            beta = [200.0_wp, 5.0_wp]
            upper(1) = 200.0_wp
            lower(2) = 5.0_wp
            lower(1) = upper(1) - 400*upper(1)*epsmac + upper(1)*epsmac
            upper(2) = lower(2) + 400*lower(2)*epsmac - lower(2)*epsmac

         elseif (itest == 20) then

            ! Test bounded ODR problem when bounds are just big enough for ndigit
            ! calculation but too small for difference calculation.
            ! Parameters start on bound.

            lun = lunrpt
            do i = 1, 2
               write (lun, 1001) itest
               write (lun, 1010)
               lun = lunsum
            end do

            setno = 9
            call set_inputs( &
               wflat, &
               title, n, m, np, q, ldwd, ld2wd, ldwe, ld2we, ldstpd, &
               x, y, beta, lower, upper, delta, ifixb, ifixx, wd, we, sclb, scld, stpb, stpd)
            job = 00000
            maxit = 100
            beta = [-200.0_wp, -5.0_wp]
            upper(1) = -200.0_wp
            lower(2) = -5.0_wp
            lower(1) = upper(1) + 400*upper(1)*epsmac
            upper(2) = lower(2) - 400*lower(2)*epsmac

         elseif (itest == 21) then

            ! Test bounded ODR problem when bounds are too small for derivative
            ! step sizes using forward differences. Parameters start on bound.

            lun = lunrpt
            do i = 1, 2
               write (lun, 1001) itest
               write (lun, 1010)
               lun = lunsum
            end do

            setno = 9
            call set_inputs( &
               wflat, &
               title, n, m, np, q, ldwd, ld2wd, ldwe, ld2we, ldstpd, &
               x, y, beta, lower, upper, delta, ifixb, ifixx, wd, we, sclb, scld, stpb, stpd)
            job = 00000
            maxit = 100
            beta = [-200.0_wp, -5.0_wp]
            upper(1) = -200.0_wp
            lower(2) = -5.0_wp
            lower(1) = upper(1) + upper(1)*epsmac
            upper(2) = lower(2) - lower(2)*epsmac

         elseif (itest == 22) then

            ! Test bounded ODR problem when first parameter is fixed and second is bounded.
            ! However, set the bounds on the first parameter to exclude the correct value
            ! of the second parameter. This will exercise the packing and unpacking of
            ! parameters and ensure that bounds and fixed parameters can be mixed.

            lun = lunrpt
            do i = 1, 2
               write (lun, 1001) itest
               write (lun, 1010)
               lun = lunsum
            end do

            setno = 10
            call set_inputs( &
               wflat, &
               title, n, m, np, q, ldwd, ld2wd, ldwe, ld2we, ldstpd, &
               x, y, beta, lower, upper, delta, ifixb, ifixx, wd, we, sclb, scld, stpb, stpd)
            job = 00010
            maxit = 100
            beta = [2.5_wp, 1.5_wp]
            lower = [2.5_wp, 1.1_wp]
            upper = [10.0_wp, 5.0_wp]
            ifixb = [0, 1]

         elseif (itest == 23) then

            ! Similar to test 22 but without bounds.

            lun = lunrpt
            do i = 1, 2
               write (lun, 1001) itest
               write (lun, 1010)
               lun = lunsum
            end do

            setno = 10
            call set_inputs( &
               wflat, &
               title, n, m, np, q, ldwd, ld2wd, ldwe, ld2we, ldstpd, &
               x, y, beta, lower, upper, delta, ifixb, ifixx, wd, we, sclb, scld, stpb, stpd)
            job = 00010
            maxit = 100
            beta = [2.5_wp, 1.5_wp]
            ifixb = [0, 1]

         end if

         ! Allocate work arrays
         if (job < 10000) then
            if (allocated(iwork)) deallocate (iwork)
            if (allocated(rwork)) deallocate (rwork)
            isodr = (job < 0 .or. mod(job, 10) <= 1)
            call workspace_dimensions(n, m, q, np, isodr, lrwork, liwork)
            allocate (iwork(liwork), rwork(lrwork))
         end if

         ! Compute solution

         write (lunrpt, 2200) title
         write (lunsum, 2200) title

         if (short) then
            call odr(fcn=fcn, &
                     n=n, m=m, q=q, np=np, &
                     beta=beta, &
                     y=y, x=x, &
                     delta=delta, &
                     job=job, &
                     iprint=iprint, lunerr=lunerr, lunrpt=lunrpt, &
                     rwork=rwork, iwork=iwork, &
                     info=info)
         else
            call odr(fcn=fcn, &
                     n=n, m=m, q=q, np=np, &
                     beta=beta, &
                     y=y, x=x, &
                     delta=delta, &
                     we=we, wd=wd, &
                     ifixb=ifixb, ifixx=ifixx, &
                     job=job, ndigit=ndigit, taufac=taufac, &
                     sstol=sstol, partol=partol, maxit=maxit, &
                     iprint=iprint, lunerr=lunerr, lunrpt=lunrpt, &
                     stpb=stpb, stpd=stpd, &
                     sclb=sclb, scld=scld, &
                     rwork=rwork, iwork=iwork, &
                     lower=lower, upper=upper, &
                     info=info)
         end if

         ! Compare results with those obtained on the CRAY YMP or the Intel Xeon running
         ! Linux using REAL(8) version of ODRPACK95

         call loc_wsse(n, m, np, q, ldwe, ld2we, isodr, wssi, wssdei, wssepi)

         wssdel = rwork(wssdei)
         wsseps = rwork(wssepi)
         wss = rwork(wssi)

         bnrm = norm2(beta)

         if (sstol < zero) then
            sstol = sqrt(epsmac)
         else
            sstol = min(sstol, one)
         end if

         if (partol < zero) then
            partol = epsmac**(two/three)
         else
            partol = min(partol, one)
         end if

         if (info >= 10000) then
            if (idpymp(itest) == info) then
               fails = .false.
               msg = 1
            else
               fails = .true.
               msg = 3
            end if

         elseif (mod(info, 10) == 1) then
            fails = abs(wss - dpymp(2, itest)) > &
                    dpymp(2, itest)*sstol*tsttol
            msg = 2

         elseif (mod(info, 10) == 2) then
            fails = abs(bnrm - dpymp(1, itest)) > &
                    dpymp(1, itest)*partol*tsttol
            msg = 2

         elseif (mod(info, 10) == 3) then
            fails = (abs(wss - dpymp(2, itest)) > &
                     dpymp(2, itest)*sstol*tsttol) &
                    .and. &
                    (abs(bnrm - dpymp(1, itest)) > &
                     dpymp(1, itest)*partol*tsttol)
            msg = 2

         elseif ((mod(info, 10) == 4) .and. &
                 (idpymp(itest) == 4)) then
            fails = .false.
            msg = 1

         elseif (info == idpymp(itest)) then
            fails = .true.
            msg = 4
         else
            fails = .true.
            msg = 3
         end if

         failed = failed .or. fails

         lun = lunrpt
         do l = 1, 2
            write (lun, 3100)
            write (lun, 3210) &
               ' CRAY YMP OR X86 RESULT = ', &
               dpymp(1, itest), dpymp(2, itest), idpymp(itest)
            write (lun, 3210) ' NEW TEST RESULT      = ', &
               bnrm, wss, info
            write (lun, 3220) ' DIFFERENCE           = ', &
               abs(dpymp(1, itest) - bnrm), abs(dpymp(2, itest) - wss)
            ewrt = abs(dpymp(1, itest))
            ewrt2 = abs(dpymp(2, itest))
            if (ewrt == zero) then
               ewrt = one
            end if
            if (ewrt2 == zero) then
               ewrt2 = one
            end if
            write (lun, 3220) ' RELATIVE ERROR       = ', &
               abs(dpymp(1, itest) - bnrm)/ewrt, &
               abs(dpymp(2, itest) - wss)/ewrt2

            if (msg == 1) then
               write (lun, 3310)
            elseif (msg == 2) then
               if (fails) then
                  write (lun, 3320)
               else
                  write (lun, 3330)
               end if
            elseif (msg == 3) then
               write (lun, 3340)
            elseif (msg == 4) then
               write (lun, 3350)
            end if

            lun = lunsum
         end do
      end do test

      if (failed) then
         write (lunrpt, 4100)
         write (lunsum, 4100)
         passed = .false.
      else
         write (lunrpt, 4200)
         write (lunsum, 4200)
         passed = .true.
      end if

      !  Format statements

1001  format(' Example ', I2,/)
1010  format(' Test simple ODR problem with analytic derivatives.')
1020  format(' Test simple OLS problem with finite difference derivatives.')
1030  format(' Test parameter fixing capabilities for poorly scaled OLS problem with'/ &
             ' analytic derivatives.')
1040  format(' Test weighting capabilities for ODR problem with analytic derivatives'/ &
             ' Also shows solution of poorly scaled ODR problem.'/ &
             ' (derivative checking turned off.)')
1050  format(' Test DELTA initialization capabilities and use of ISTOP to restrict'/ &
             ' parameter values for ODR problem with analytic derivatives.')
1060  format(' Test stiff stopping conditions for unscaled ODR problem with'/ &
             ' analytic derivatives.')
1070  format(' Test restart for unscaled ODR problem with analytic derivatives.')
1080  format(' Test use of TAUFAC to restrict first step for ODR problem with'/ &
             ' finite difference derivatives.')
1090  format(' Test implicit model for OLS problem.')
1100  format(' Test multiresponse model for ODR problem with finite difference derivatives.')
1110  format(' Test detection of questionable analytic derivatives for OLS problem.')
1120  format(' Test detection of incorrect analytic derivatives for ODR problem with'/ &
             ' analytic derivatives.')
2200  format(' Data Set Reference: ', A80)
3100  format &
         (/' Comparison of new results with', &
           ' REAL (wp) Cray YMP or Intel X86 (Linux) '/ &
           ' Result:'// &
           '                         Norm of BETA', &
           '        Sum of Squared WTD OBS Errors  INFO')
3210  format &
         (/A25/1P, 2E37.30, I6)
3220  format &
         (/A25, 1P, D12.5, 25X, D12.5, I6)
3310  format &
         (/' *** Stopping conditions', &
           ' show convergence not attained. ***'/ &
           '        no further comparisons made between results.'//)
3320  format &
         (//' *** WARNING ***', &
           ' results do not agree to within stopping tolerance. ***'//)
3330  format &
         (//' *** Results agree to within stopping tolerance. ***'//)
3340  format &
         (//' *** WARNING ***', &
           ' stopping conditions do not agree. ***'//)
3350  format &
         (//' *** WARNING ***', &
           ' unexpected stopping condition.', &
           '  please contact package authors. ***'//)
4100  format &
         (/// &
           ' *** Summary:', &
           ' one or more tests do not agree with expected results. ***')
4200  format &
         (/// &
           ' *** Summary:', &
           ' all tests agree with expected results. ***')

   end subroutine test_odr

   impure subroutine set_inputs( &
      wflat, &
      title, n, m, np, q, ldwd, ld2wd, ldwe, ld2we, ldstpd, &
      x, y, beta, lower, upper, delta, ifixb, ifixx, wd, we, sclb, scld, stpb, stpd)
   !! Set up input data for odrpack exerciser.

      logical, intent(in) :: wflat
      character(len=80), intent(out) :: title
      integer, intent(out) :: n
      integer, intent(out) :: m
      integer, intent(out) :: np
      integer, intent(out) :: q
      integer, intent(out) :: ldwd
      integer, intent(out) :: ld2wd
      integer, intent(out) :: ldwe
      integer, intent(out) :: ld2we
      integer, intent(out) :: ldstpd
      real(wp), intent(out), allocatable :: x(:, :)
      real(wp), intent(out), allocatable :: y(:, :)
      real(wp), intent(out), allocatable :: beta(:)
      real(wp), intent(out), allocatable :: lower(:)
      real(wp), intent(out), allocatable :: upper(:)
      real(wp), intent(out), allocatable :: delta(:, :)
      integer, intent(out), allocatable :: ifixb(:)
      integer, intent(out), allocatable :: ifixx(:, :)
      real(wp), intent(out), allocatable :: wd(:, :, :)
      real(wp), intent(out), allocatable :: we(:, :, :)
      real(wp), intent(out), allocatable :: sclb(:)
      real(wp), intent(out), allocatable :: scld(:, :)
      real(wp), intent(out), allocatable :: stpb(:)
      real(wp), intent(out), allocatable :: stpd(:, :)

      ! Parameters
      integer, parameter :: maxset = 16

      ! Local scalars
      integer :: k

      ! Local arrays
      real(wp) :: bdata(maxnp, maxset), xdata(maxn, maxm, maxset), ydata(maxn, maxq, maxset)
      integer :: mdata(maxset), ndata(maxset), npdata(maxset), qdata(maxset)
      character(len=80) :: tdata(maxset)

      ! Data statements
      data &
         tdata(1) &
         /' BOGGS, BYRD AND SCHNABEL, 1985, EXAMPLE 1'/
      data &
         ndata(1), mdata(1), npdata(1), qdata(1) &
         /40, 1, 2, 1/
      data &
         (bdata(k, 1), k=1, 2) &
         /1.0E+0_wp, 1.0E+0_wp/
      data &
         ydata(1, 1, 1), xdata(1, 1, 1) &
         /-0.119569795672791172E+1_wp, -0.213701920211315155E-1_wp/
      data &
         ydata(2, 1, 1), xdata(2, 1, 1) &
         /-0.128023349509594288E+1_wp, 0.494813247025012969E-1_wp/
      data &
         ydata(3, 1, 1), xdata(3, 1, 1) &
         /-0.125270693343174591E+1_wp, 0.127889194935560226E+0_wp/
      data &
         ydata(4, 1, 1), xdata(4, 1, 1) &
         /-0.996698267935287383E+0_wp, 0.128615394085645676E+0_wp/
      data &
         ydata(5, 1, 1), xdata(5, 1, 1) &
         /-0.104681033065801934E+1_wp, 0.232544285655021667E+0_wp/
      data &
         ydata(6, 1, 1), xdata(6, 1, 1) &
         /-0.146724952092847308E+1_wp, 0.268151108026504516E+0_wp/
      data &
         ydata(7, 1, 1), xdata(7, 1, 1) &
         /-0.123366891873487528E+1_wp, 0.309041029810905456E+0_wp/
      data &
         ydata(8, 1, 1), xdata(8, 1, 1) &
         /-0.165665097907185554E+1_wp, 0.405991539210081099E+0_wp/
      data &
         ydata(9, 1, 1), xdata(9, 1, 1) &
         /-0.168476460930907119E+1_wp, 0.376611424833536147E+0_wp/
      data &
         ydata(10, 1, 1), xdata(10, 1, 1) &
         /-0.198571971169224491E+1_wp, 0.475875890851020811E+0_wp/
      data &
         ydata(11, 1, 1), xdata(11, 1, 1) &
         /-0.195691696638051344E+1_wp, 0.499246935397386550E+0_wp/
      data &
         ydata(12, 1, 1), xdata(12, 1, 1) &
         /-0.211871342665769836E+1_wp, 0.536615037024021147E+0_wp/
      data &
         ydata(13, 1, 1), xdata(13, 1, 1) &
         /-0.268642932558671020E+1_wp, 0.581830765902996060E+0_wp/
      data &
         ydata(14, 1, 1), xdata(14, 1, 1) &
         /-0.281123260058024347E+1_wp, 0.684512710422277446E+0_wp/
      data &
         ydata(15, 1, 1), xdata(15, 1, 1) &
         /-0.328704486581785920E+1_wp, 0.660219819694757458E+0_wp/
      data &
         ydata(16, 1, 1), xdata(16, 1, 1) &
         /-0.423062993461887032E+1_wp, 0.766990323960781092E+0_wp/
      data &
         ydata(17, 1, 1), xdata(17, 1, 1) &
         /-0.512043906552226903E+1_wp, 0.808270426690578456E+0_wp/
      data &
         ydata(18, 1, 1), xdata(18, 1, 1) &
         /-0.731032616379005535E+1_wp, 0.897410020083189004E+0_wp/
      data &
         ydata(19, 1, 1), xdata(19, 1, 1) &
         /-0.109002759485608993E+2_wp, 0.959199774116277687E+0_wp/
      data &
         ydata(20, 1, 1), xdata(20, 1, 1) &
         /-0.251810238510370206E+2_wp, 0.914675474762916558E+0_wp/
      data &
         ydata(21, 1, 1), xdata(21, 1, 1) &
         /0.100123028650879944E+3_wp, 0.997759691476821892E+0_wp/
      data &
         ydata(22, 1, 1), xdata(22, 1, 1) &
         /0.168225085871915048E+2_wp, 0.107136870384216308E+1_wp/
      data &
         ydata(23, 1, 1), xdata(23, 1, 1) &
         /0.894830510866913009E+1_wp, 0.108033321037888526E+1_wp/
      data &
         ydata(24, 1, 1), xdata(24, 1, 1) &
         /0.645853815227747004E+1_wp, 0.116064198672771453E+1_wp/
      data &
         ydata(25, 1, 1), xdata(25, 1, 1) &
         /0.498218564760117328E+1_wp, 0.119080889359116553E+1_wp/
      data &
         ydata(26, 1, 1), xdata(26, 1, 1) &
         /0.382971664718710476E+1_wp, 0.129418875187635420E+1_wp/
      data &
         ydata(27, 1, 1), xdata(27, 1, 1) &
         /0.344116492497344184E+1_wp, 0.135594148099422453E+1_wp/
      data &
         ydata(28, 1, 1), xdata(28, 1, 1) &
         /0.276840496973858949E+1_wp, 0.135302808716893195E+1_wp/
      data &
         ydata(29, 1, 1), xdata(29, 1, 1) &
         /0.259521665196956666E+1_wp, 0.137994666010141371E+1_wp/
      data &
         ydata(30, 1, 1), xdata(30, 1, 1) &
         /0.205996022794557661E+1_wp, 0.147630019545555113E+1_wp/
      data &
         ydata(31, 1, 1), xdata(31, 1, 1) &
         /0.197939614345337836E+1_wp, 0.153450708076357840E+1_wp/
      data &
         ydata(32, 1, 1), xdata(32, 1, 1) &
         /0.156739340562905589E+1_wp, 0.152805351451039313E+1_wp/
      data &
         ydata(33, 1, 1), xdata(33, 1, 1) &
         /0.159032057073028366E+1_wp, 0.157147316247224806E+1_wp/
      data &
         ydata(34, 1, 1), xdata(34, 1, 1) &
         /0.173102268158937949E+1_wp, 0.166649596005678175E+1_wp/
      data &
         ydata(35, 1, 1), xdata(35, 1, 1) &
         /0.155512561664824758E+1_wp, 0.166505665838718412E+1_wp/
      data &
         ydata(36, 1, 1), xdata(36, 1, 1) &
         /0.149635994944133260E+1_wp, 0.175214128553867338E+1_wp/
      data &
         ydata(37, 1, 1), xdata(37, 1, 1) &
         /0.147487601463073568E+1_wp, 0.180567992463707922E+1_wp/
      data &
         ydata(38, 1, 1), xdata(38, 1, 1) &
         /0.117244575233306998E+1_wp, 0.184624404296278952E+1_wp/
      data &
         ydata(39, 1, 1), xdata(39, 1, 1) &
         /0.910931336069172580E+0_wp, 0.195568727388978002E+1_wp/
      data &
         ydata(40, 1, 1), xdata(40, 1, 1) &
         /0.126172980914513272E+1_wp, 0.199326394036412237E+1_wp/

      data &
         tdata(2) &
         /' BOGGS, BYRD AND SCHNABEL, 1985, EXAMPLE 2'/
      data &
         ndata(2), mdata(2), npdata(2), qdata(2) &
         /50, 2, 3, 1/
      data &
         (bdata(k, 2), k=1, 3) &
         /-1.0E+0_wp, 1.0E+0_wp, 1.0E+0_wp/
      data &
         ydata(1, 1, 2), xdata(1, 1, 2), xdata(1, 2, 2) &
         /0.680832777217942900E+0_wp, &
         0.625474598833994800E-1_wp, 0.110179064209783100E+0_wp/
      data &
         ydata(2, 1, 2), xdata(2, 1, 2), xdata(2, 2, 2) &
         /0.122183594595302200E+1_wp, &
         0.202500343620642400E+0_wp, -0.196140862891327600E-1_wp/
      data &
         ydata(3, 1, 2), xdata(3, 1, 2), xdata(3, 2, 2) &
         /0.118958678734608200E+1_wp, &
         0.164943738599876500E+0_wp, 0.166514874750996600E+0_wp/
      data &
         ydata(4, 1, 2), xdata(4, 1, 2), xdata(4, 2, 2) &
         /0.146982623764094600E+1_wp, &
         0.304874137610506100E+0_wp, 0.612908688041490500E-2_wp/
      data &
         ydata(5, 1, 2), xdata(5, 1, 2), xdata(5, 2, 2) &
         /0.167775338189355300E+1_wp, &
         0.532727445580665100E+0_wp, 0.938248787552444600E-1_wp/
      data &
         ydata(6, 1, 2), xdata(6, 1, 2), xdata(6, 2, 2) &
         /0.202485721906026200E+1_wp, &
         0.508823707598910200E+0_wp, 0.499605775020505400E-2_wp/
      data &
         ydata(7, 1, 2), xdata(7, 1, 2), xdata(7, 2, 2) &
         /0.258912851935938800E+1_wp, &
         0.704227041878554000E+0_wp, 0.819354849092326200E-1_wp/
      data &
         ydata(8, 1, 2), xdata(8, 1, 2), xdata(8, 2, 2) &
         /0.366894203254154800E+1_wp, &
         0.592077736111512000E+0_wp, 0.127113960672389100E-1_wp/
      data &
         ydata(9, 1, 2), xdata(9, 1, 2), xdata(9, 2, 2) &
         /0.574609583351347300E+1_wp, &
         0.104940945646421600E+1_wp, 0.258095243658316100E-1_wp/
      data &
         ydata(10, 1, 2), xdata(10, 1, 2), xdata(10, 2, 2) &
         /0.127676424026489300E+2_wp, &
         0.979382517558619200E+0_wp, 0.124280755181027900E+0_wp/
      data &
         ydata(11, 1, 2), xdata(11, 1, 2), xdata(11, 2, 2) &
         /0.123473079693623100E+1_wp, &
         0.637870453165538700E-1_wp, 0.304856401137196400E+0_wp/
      data &
         ydata(12, 1, 2), xdata(12, 1, 2), xdata(12, 2, 2) &
         /0.142256120864082800E+1_wp, &
         0.176123312906025700E+0_wp, 0.262387028078896900E+0_wp/
      data &
         ydata(13, 1, 2), xdata(13, 1, 2), xdata(13, 2, 2) &
         /0.169889534013024700E+1_wp, &
         0.310965082300263000E+0_wp, 0.226430765474758800E+0_wp/
      data &
         ydata(14, 1, 2), xdata(14, 1, 2), xdata(14, 2, 2) &
         /0.173485577901204400E+1_wp, &
         0.311394269116782100E+0_wp, 0.271375840410281800E+0_wp/
      data &
         ydata(15, 1, 2), xdata(15, 1, 2), xdata(15, 2, 2) &
         /0.277761263972834600E+1_wp, &
         0.447076126190612500E+0_wp, 0.255000858902618300E+0_wp/
      data &
         ydata(16, 1, 2), xdata(16, 1, 2), xdata(16, 2, 2) &
         /0.339163324662617300E+1_wp, &
         0.384786230998211100E+0_wp, 0.154958003178364000E+0_wp/
      data &
         ydata(17, 1, 2), xdata(17, 1, 2), xdata(17, 2, 2) &
         /0.589615137312147500E+1_wp, &
         0.649093176450780500E+0_wp, 0.258301685463773200E+0_wp/
      data &
         ydata(18, 1, 2), xdata(18, 1, 2), xdata(18, 2, 2) &
         /0.124415625214576800E+2_wp, &
         0.685612005372525500E+0_wp, 0.107391260603228600E+0_wp/
      data &
         ydata(19, 1, 2), xdata(19, 1, 2), xdata(19, 2, 2) &
         /-0.498491739153861600E+2_wp, &
         0.968747139425088400E+0_wp, 0.151932526135740700E+0_wp/
      data &
         ydata(20, 1, 2), xdata(20, 1, 2), xdata(20, 2, 2) &
         /-0.832795509000618600E+1_wp, &
         0.869789367989532900E+0_wp, 0.625507500586400000E-1_wp/
      data &
         ydata(21, 1, 2), xdata(21, 1, 2), xdata(21, 2, 2) &
         /0.184934617774239900E+1_wp, &
         -0.465309930332736600E-2_wp, 0.546795662595375200E+0_wp/
      data &
         ydata(22, 1, 2), xdata(22, 1, 2), xdata(22, 2, 2) &
         /0.175192979176839200E+1_wp, &
         0.604753397196646000E-2_wp, 0.230905749473922700E+0_wp/
      data &
         ydata(23, 1, 2), xdata(23, 1, 2), xdata(23, 2, 2) &
         /0.253949381238535800E+1_wp, &
         0.239418809621756000E+0_wp, 0.190752069681170700E+0_wp/
      data &
         ydata(24, 1, 2), xdata(24, 1, 2), xdata(24, 2, 2) &
         /0.373500774928501700E+1_wp, &
         0.456662468911699800E+0_wp, 0.328870615170984400E+0_wp/
      data &
         ydata(25, 1, 2), xdata(25, 1, 2), xdata(25, 2, 2) &
         /0.548408128950331000E+1_wp, &
         0.371115320522079500E+0_wp, 0.439978556640660500E+0_wp/
      data &
         ydata(26, 1, 2), xdata(26, 1, 2), xdata(26, 2, 2) &
         /0.125256880521774300E+2_wp, &
         0.586442107042503000E+0_wp, 0.490689043752286700E+0_wp/
      data &
         ydata(27, 1, 2), xdata(27, 1, 2), xdata(27, 2, 2) &
         /-0.493587797164916600E+2_wp, &
         0.579796274973298000E+0_wp, 0.521860998203383100E+0_wp/
      data &
         ydata(28, 1, 2), xdata(28, 1, 2), xdata(28, 2, 2) &
         /-0.801158974965412700E+1_wp, &
         0.805008094903899900E+0_wp, 0.292283538955391600E+0_wp/
      data &
         ydata(29, 1, 2), xdata(29, 1, 2), xdata(29, 2, 2) &
         /-0.437399487061934100E+1_wp, &
         0.637242340835710000E+0_wp, 0.402261740352486000E+0_wp/
      data &
         ydata(30, 1, 2), xdata(30, 1, 2), xdata(30, 2, 2) &
         /-0.297800103425979600E+1_wp, &
         0.982132817936118700E+0_wp, 0.392546836419047000E+0_wp/
      data &
         ydata(31, 1, 2), xdata(31, 1, 2), xdata(31, 2, 2) &
         /0.271811057454661300E+1_wp, &
         -0.223515657121262700E-1_wp, 0.650479019708978800E+0_wp/
      data &
         ydata(32, 1, 2), xdata(32, 1, 2), xdata(32, 2, 2) &
         /0.377035865613392400E+1_wp, &
         0.136081427545033600E+0_wp, 0.753020101897661800E+0_wp/
      data &
         ydata(33, 1, 2), xdata(33, 1, 2), xdata(33, 2, 2) &
         /0.560111053917143100E+1_wp, &
         0.145367053019870600E+0_wp, 0.611153532003093100E+0_wp/
      data &
         ydata(34, 1, 2), xdata(34, 1, 2), xdata(34, 2, 2) &
         /0.128152376174926800E+2_wp, &
         0.308221919576435500E+0_wp, 0.455217283290423900E+0_wp/
      data &
         ydata(35, 1, 2), xdata(35, 1, 2), xdata(35, 2, 2) &
         /-0.498709177732467200E+2_wp, &
         0.432658769133528300E+0_wp, 0.678607663414113000E+0_wp/
      data &
         ydata(36, 1, 2), xdata(36, 1, 2), xdata(36, 2, 2) &
         /-0.815797696908314300E+1_wp, &
         0.477785501079980300E+0_wp, 0.536178207572157000E+0_wp/
      data &
         ydata(37, 1, 2), xdata(37, 1, 2), xdata(37, 2, 2) &
         /-0.440240491195158600E+1_wp, &
         0.727986827616619000E+0_wp, 0.668497920573493900E+0_wp/
      data &
         ydata(38, 1, 2), xdata(38, 1, 2), xdata(38, 2, 2) &
         /-0.276723957061767500E+1_wp, &
         0.745950385588265100E+0_wp, 0.786077589007263700E+0_wp/
      data &
         ydata(39, 1, 2), xdata(39, 1, 2), xdata(39, 2, 2) &
         /-0.223203667288734800E+1_wp, &
         0.732537503527113500E+0_wp, 0.582625164046828400E+0_wp/
      data &
         ydata(40, 1, 2), xdata(40, 1, 2), xdata(40, 2, 2) &
         /-0.169728270310622000E+1_wp, &
         0.967352361433846300E+0_wp, 0.460779396016832800E+0_wp/
      data &
         ydata(41, 1, 2), xdata(41, 1, 2), xdata(41, 2, 2) &
         /0.551015652153227000E+1_wp, &
         0.129761784310891100E-1_wp, 0.700009537931860000E+0_wp/
      data &
         ydata(42, 1, 2), xdata(42, 1, 2), xdata(42, 2, 2) &
         /0.128036180496215800E+2_wp, &
         0.170163243950629700E+0_wp, 0.853131830764348700E+0_wp/
      data &
         ydata(43, 1, 2), xdata(43, 1, 2), xdata(43, 2, 2) &
         /-0.498257683396339000E+2_wp, &
         0.162768461906274000E+0_wp, 0.865315129048175000E+0_wp/
      data &
         ydata(44, 1, 2), xdata(44, 1, 2), xdata(44, 2, 2) &
         /-0.877334550221761900E+1_wp, &
         0.222914807946165800E+0_wp, 0.797511758502094500E+0_wp/
      data &
         ydata(45, 1, 2), xdata(45, 1, 2), xdata(45, 2, 2) &
         /-0.453820192156867600E+1_wp, &
         0.402910095604624900E+0_wp, 0.761492958727023100E+0_wp/
      data &
         ydata(46, 1, 2), xdata(46, 1, 2), xdata(46, 2, 2) &
         /-0.297499315738677900E+1_wp, &
         0.233770812593443200E+0_wp, 0.896000095844223500E+0_wp/
      data &
         ydata(47, 1, 2), xdata(47, 1, 2), xdata(47, 2, 2) &
         /-0.212743255978538900E+1_wp, &
         0.646528693486914700E+0_wp, 0.968574333700755700E+0_wp/
      data &
         ydata(48, 1, 2), xdata(48, 1, 2), xdata(48, 2, 2) &
         /-0.209703205365401000E+1_wp, &
         0.802811658568969400E+0_wp, 0.904866450476711600E+0_wp/
      data &
         ydata(49, 1, 2), xdata(49, 1, 2), xdata(49, 2, 2) &
         /-0.155287292042086200E+1_wp, &
         0.837137859891222900E+0_wp, 0.835684424990021900E+0_wp/
      data &
         ydata(50, 1, 2), xdata(50, 1, 2), xdata(50, 2, 2) &
         /-0.161356673770480700E+1_wp, &
         0.103165980756526600E+1_wp, 0.793902191912346100E+0_wp/

      data &
         tdata(3) &
         /' BOGGS, BYRD AND SCHNABEL, 1985, EXAMPLE 3'/
      data &
         ndata(3), mdata(3), npdata(3), qdata(3) &
         /44, 1, 9, 1/
      data &
         (bdata(k, 3), k=1, 9) &
         /0.281887509408440189E-5_wp, &
         -0.231290549212363845E-2_wp, 0.583035555572801965E+1_wp, &
         0.000000000000000000E+0_wp, 0.406910776203121026E+8_wp, &
         0.138001105225000000E-2_wp, 0.596038513209999999E-1_wp, &
         0.670582099359999998E+1_wp, 0.106994410000000000E+10_wp/
      data &
         ydata(1, 1, 3), xdata(1, 1, 3) &
         /0.988227696721327788E+0_wp, 0.25E-8_wp/
      data &
         ydata(2, 1, 3), xdata(2, 1, 3) &
         /0.988268083998559958E+0_wp, 0.64E-8_wp/
      data &
         ydata(3, 1, 3), xdata(3, 1, 3) &
         /0.988341022958438831E+0_wp, 1.0E-8_wp/
      data &
         ydata(4, 1, 3), xdata(4, 1, 3) &
         /0.988380557606306446E+0_wp, 0.9E-7_wp/
      data &
         ydata(5, 1, 3), xdata(5, 1, 3) &
         /0.988275062411751338E+0_wp, 1.0E-6_wp/
      data &
         ydata(6, 1, 3), xdata(6, 1, 3) &
         /0.988326680176446987E+0_wp, 0.4E-5_wp/
      data &
         ydata(7, 1, 3), xdata(7, 1, 3) &
         /0.988306058860433439E+0_wp, 0.9E-5_wp/
      data &
         ydata(8, 1, 3), xdata(8, 1, 3) &
         /0.988292880079125555E+0_wp, 0.16E-4_wp/
      data &
         ydata(9, 1, 3), xdata(9, 1, 3) &
         /0.988305279259496905E+0_wp, 0.36E-4_wp/
      data &
         ydata(10, 1, 3), xdata(10, 1, 3) &
         /0.988278142019574202E+0_wp, 0.64E-4_wp/
      data &
         ydata(11, 1, 3), xdata(11, 1, 3) &
         /0.988224953369819946E+0_wp, 1.0E-4_wp/
      data &
         ydata(12, 1, 3), xdata(12, 1, 3) &
         /0.988111989169778223E+0_wp, 0.144E-3_wp/
      data &
         ydata(13, 1, 3), xdata(13, 1, 3) &
         /0.988045627103840613E+0_wp, 0.225E-3_wp/
      data &
         ydata(14, 1, 3), xdata(14, 1, 3) &
         /0.987913715667047655E+0_wp, 0.400E-3_wp/
      data &
         ydata(15, 1, 3), xdata(15, 1, 3) &
         /0.987841994238525678E+0_wp, 0.625E-3_wp/
      data &
         ydata(16, 1, 3), xdata(16, 1, 3) &
         /0.987638450432434270E+0_wp, 0.900E-3_wp/
      data &
         ydata(17, 1, 3), xdata(17, 1, 3) &
         /0.987587364331771395E+0_wp, 0.1225E-2_wp/
      data &
         ydata(18, 1, 3), xdata(18, 1, 3) &
         /0.987576264149633684E+0_wp, 0.1600E-2_wp/
      data &
         ydata(19, 1, 3), xdata(19, 1, 3) &
         /0.987539209110983643E+0_wp, 0.2025E-2_wp/
      data &
         ydata(20, 1, 3), xdata(20, 1, 3) &
         /0.987621143807705698E+0_wp, 0.25E-2_wp/
      data &
         ydata(21, 1, 3), xdata(21, 1, 3) &
         /0.988023229785526217E+0_wp, 0.36E-2_wp/
      data &
         ydata(22, 1, 3), xdata(22, 1, 3) &
         /0.988558376710994197E+0_wp, 0.49E-2_wp/
      data &
         ydata(23, 1, 3), xdata(23, 1, 3) &
         /0.989304775352439885E+0_wp, 0.64E-2_wp/
      data &
         ydata(24, 1, 3), xdata(24, 1, 3) &
         /0.990210452265710472E+0_wp, 0.81E-2_wp/
      data &
         ydata(25, 1, 3), xdata(25, 1, 3) &
         /0.991095950592263900E+0_wp, 1.00E-2_wp/
      data &
         ydata(26, 1, 3), xdata(26, 1, 3) &
         /0.991475677297119272E+0_wp, 0.11025E-1_wp/
      data &
         ydata(27, 1, 3), xdata(27, 1, 3) &
         /0.991901306250746771E+0_wp, 0.12100E-1_wp/
      data &
         ydata(28, 1, 3), xdata(28, 1, 3) &
         /0.992619222425303263E+0_wp, 0.14400E-1_wp/
      data &
         ydata(29, 1, 3), xdata(29, 1, 3) &
         /0.993617037631973475E+0_wp, 0.16900E-1_wp/
      data &
         ydata(30, 1, 3), xdata(30, 1, 3) &
         /0.994727321698030676E+0_wp, 0.19600E-1_wp/
      data &
         ydata(31, 1, 3), xdata(31, 1, 3) &
         /0.996523114720326189E+0_wp, 0.25600E-1_wp/
      data &
         ydata(32, 1, 3), xdata(32, 1, 3) &
         /0.998036909563764020E+0_wp, 0.32400E-1_wp/
      data &
         ydata(33, 1, 3), xdata(33, 1, 3) &
         /0.999151968626971372E+0_wp, 0.40000E-1_wp/
      data &
         ydata(34, 1, 3), xdata(34, 1, 3) &
         /0.100017083706131769E+1_wp, 0.50625E-1_wp/
      data &
         ydata(35, 1, 3), xdata(35, 1, 3) &
         /0.100110046382923523E+1_wp, 0.75625E-1_wp/
      data &
         ydata(36, 1, 3), xdata(36, 1, 3) &
         /0.100059103180404652E+1_wp, 0.12250E+0_wp/
      data &
         ydata(37, 1, 3), xdata(37, 1, 3) &
         /0.999211829791257561E+0_wp, 0.16000E+0_wp/
      data &
         ydata(38, 1, 3), xdata(38, 1, 3) &
         /0.994711451526761862E+0_wp, 0.25000E+0_wp/
      data &
         ydata(39, 1, 3), xdata(39, 1, 3) &
         /0.989844132928847109E+0_wp, 0.33640E+0_wp/
      data &
         ydata(40, 1, 3), xdata(40, 1, 3) &
         /0.987234104554490439E+0_wp, 0.38440E+0_wp/
      data &
         ydata(41, 1, 3), xdata(41, 1, 3) &
         /0.980928240178404887E+0_wp, 0.49E+0_wp/
      data &
         ydata(42, 1, 3), xdata(42, 1, 3) &
         /0.970888680366055576E+0_wp, 0.64E+0_wp/
      data &
         ydata(43, 1, 3), xdata(43, 1, 3) &
         /0.960043769857327398E+0_wp, 0.81E+0_wp/
      data &
         ydata(44, 1, 3), xdata(44, 1, 3) &
         /0.947277159259551068E+0_wp, 1.00E+0_wp/

      data &
         tdata(4) &
         /' HIMMELBLAU, 1970, EXAMPLE 6.2-4, PAGE 188'/
      data &
         ndata(4), mdata(4), npdata(4), qdata(4) &
         /13, 2, 3, 1/
      data &
         (bdata(k, 4), k=1, 3) &
         /3.0E+0_wp, 3.0E+0_wp, -0.5E+0_wp/
      data &
         ydata(1, 1, 4), xdata(1, 1, 4), xdata(1, 2, 4) &
         /2.93E+0_wp, 0.0E+0_wp, 0.0E+0_wp/
      data &
         ydata(2, 1, 4), xdata(2, 1, 4), xdata(2, 2, 4) &
         /1.95E+0_wp, 0.0E+0_wp, 1.0E+0_wp/
      data &
         ydata(3, 1, 4), xdata(3, 1, 4), xdata(3, 2, 4) &
         /0.81E+0_wp, 0.0E+0_wp, 2.0E+0_wp/
      data &
         ydata(4, 1, 4), xdata(4, 1, 4), xdata(4, 2, 4) &
         /0.58E+0_wp, 0.0E+0_wp, 3.0E+0_wp/
      data &
         ydata(5, 1, 4), xdata(5, 1, 4), xdata(5, 2, 4) &
         /5.90E+0_wp, 1.0E+0_wp, 0.0E+0_wp/
      data &
         ydata(6, 1, 4), xdata(6, 1, 4), xdata(6, 2, 4) &
         /4.74E+0_wp, 1.0E+0_wp, 1.0E+0_wp/
      data &
         ydata(7, 1, 4), xdata(7, 1, 4), xdata(7, 2, 4) &
         /4.18E+0_wp, 1.0E+0_wp, 2.0E+0_wp/
      data &
         ydata(8, 1, 4), xdata(8, 1, 4), xdata(8, 2, 4) &
         /4.05E+0_wp, 1.0E+0_wp, 2.0E+0_wp/
      data &
         ydata(9, 1, 4), xdata(9, 1, 4), xdata(9, 2, 4) &
         /9.03E+0_wp, 2.0E+0_wp, 0.0E+0_wp/
      data &
         ydata(10, 1, 4), xdata(10, 1, 4), xdata(10, 2, 4) &
         /7.85E+0_wp, 2.0E+0_wp, 1.0E+0_wp/
      data &
         ydata(11, 1, 4), xdata(11, 1, 4), xdata(11, 2, 4) &
         /7.22E+0_wp, 2.0E+0_wp, 2.0E+0_wp/
      data &
         ydata(12, 1, 4), xdata(12, 1, 4), xdata(12, 2, 4) &
         /8.50E+0_wp, 2.5E+0_wp, 2.0E+0_wp/
      data &
         ydata(13, 1, 4), xdata(13, 1, 4), xdata(13, 2, 4) &
         /9.81E+0_wp, 2.9E+0_wp, 1.8E+0_wp/

      data &
         tdata(5) &
         /' DRAPER AND SMITH, 1981, EXERCISE I, PAGE 521-522'/
      data &
         ndata(5), mdata(5), npdata(5), qdata(5) &
         /8, 2, 2, 1/
      data &
         (bdata(k, 5), k=1, 2) &
         /0.01155E+0_wp, 5000.0E+0_wp/
      data &
         ydata(1, 1, 5), xdata(1, 1, 5), xdata(1, 2, 5) &
         /0.912E+0_wp, 109.0E+0_wp, 600.0E+0_wp/
      data &
         ydata(2, 1, 5), xdata(2, 1, 5), xdata(2, 2, 5) &
         /0.382E+0_wp, 65.0E+0_wp, 640.0E+0_wp/
      data &
         ydata(3, 1, 5), xdata(3, 1, 5), xdata(3, 2, 5) &
         /0.397E+0_wp, 1180.0E+0_wp, 600.0E+0_wp/
      data &
         ydata(4, 1, 5), xdata(4, 1, 5), xdata(4, 2, 5) &
         /0.376E+0_wp, 66.0E+0_wp, 640.0E+0_wp/
      data &
         ydata(5, 1, 5), xdata(5, 1, 5), xdata(5, 2, 5) &
         /0.342E+0_wp, 1270.0E+0_wp, 600.0E+0_wp/
      data &
         ydata(6, 1, 5), xdata(6, 1, 5), xdata(6, 2, 5) &
         /0.358E+0_wp, 69.0E+0_wp, 640.0E+0_wp/
      data &
         ydata(7, 1, 5), xdata(7, 1, 5), xdata(7, 2, 5) &
         /0.348E+0_wp, 1230.0E+0_wp, 600.0E+0_wp/
      data &
         ydata(8, 1, 5), xdata(8, 1, 5), xdata(8, 2, 5) &
         /0.376E+0_wp, 68.0E+0_wp, 640.0E+0_wp/

      data &
         tdata(6) &
         /' POWELL AND MACDONALD, 1972, TABLES 7 AND 8, PAGES 153-154'/
      data &
         ndata(6), mdata(6), npdata(6), qdata(6) &
         /14, 1, 3, 1/
      data &
         (bdata(k, 6), k=1, 3) &
         /25.0E+0_wp, 30.0E+0_wp, 6.0E+0_wp/
      data &
         ydata(1, 1, 6), xdata(1, 1, 6) &
         /26.38E+0_wp, 1.0E+0_wp/
      data &
         ydata(2, 1, 6), xdata(2, 1, 6) &
         /25.79E+0_wp, 2.0E+0_wp/
      data &
         ydata(3, 1, 6), xdata(3, 1, 6) &
         /25.29E+0_wp, 3.0E+0_wp/
      data &
         ydata(4, 1, 6), xdata(4, 1, 6) &
         /24.86E+0_wp, 4.0E+0_wp/
      data &
         ydata(5, 1, 6), xdata(5, 1, 6) &
         /24.46E+0_wp, 5.0E+0_wp/
      data &
         ydata(6, 1, 6), xdata(6, 1, 6) &
         /24.10E+0_wp, 6.0E+0_wp/
      data &
         ydata(7, 1, 6), xdata(7, 1, 6) &
         /23.78E+0_wp, 7.0E+0_wp/
      data &
         ydata(8, 1, 6), xdata(8, 1, 6) &
         /23.50E+0_wp, 8.0E+0_wp/
      data &
         ydata(9, 1, 6), xdata(9, 1, 6) &
         /23.24E+0_wp, 9.0E+0_wp/
      data &
         ydata(10, 1, 6), xdata(10, 1, 6) &
         /23.00E+0_wp, 10.0E+0_wp/
      data &
         ydata(11, 1, 6), xdata(11, 1, 6) &
         /22.78E+0_wp, 11.0E+0_wp/
      data &
         ydata(12, 1, 6), xdata(12, 1, 6) &
         /22.58E+0_wp, 12.0E+0_wp/
      data &
         ydata(13, 1, 6), xdata(13, 1, 6) &
         /22.39E+0_wp, 13.0E+0_wp/
      data &
         ydata(14, 1, 6), xdata(14, 1, 6) &
         /22.22E+0_wp, 14.0E+0_wp/

      data &
         tdata(7) &
         /' FULLER, 1987, TABLE 3.2.10, PAGES 244-245'/
      data &
         ndata(7), mdata(7), npdata(7), qdata(7) &
         /20, 2, 5, 1/
      data &
         (bdata(k, 7), k=1, 5) &
         /-1.0E+0_wp, -3.0E+0_wp, 0.09E+0_wp, 0.02E+0_wp, 0.08E+0_wp/
      data &
         ydata(1, 1, 7), xdata(1, 1, 7), xdata(1, 2, 7) &
         /0.0E+0_wp, 0.50E+0_wp, -0.12E+0_wp/
      data &
         ydata(2, 1, 7), xdata(2, 1, 7), xdata(2, 2, 7) &
         /0.0E+0_wp, 1.20E+0_wp, -0.60E+0_wp/
      data &
         ydata(3, 1, 7), xdata(3, 1, 7), xdata(3, 2, 7) &
         /0.0E+0_wp, 1.60E+0_wp, -1.00E+0_wp/
      data &
         ydata(4, 1, 7), xdata(4, 1, 7), xdata(4, 2, 7) &
         /0.0E+0_wp, 1.86E+0_wp, -1.40E+0_wp/
      data &
         ydata(5, 1, 7), xdata(5, 1, 7), xdata(5, 2, 7) &
         /0.0E+0_wp, 2.12E+0_wp, -2.54E+0_wp/
      data &
         ydata(6, 1, 7), xdata(6, 1, 7), xdata(6, 2, 7) &
         /0.0E+0_wp, 2.36E+0_wp, -3.36E+0_wp/
      data &
         ydata(7, 1, 7), xdata(7, 1, 7), xdata(7, 2, 7) &
         /0.0E+0_wp, 2.44E+0_wp, -4.00E+0_wp/
      data &
         ydata(8, 1, 7), xdata(8, 1, 7), xdata(8, 2, 7) &
         /0.0E+0_wp, 2.36E+0_wp, -4.75E+0_wp/
      data &
         ydata(9, 1, 7), xdata(9, 1, 7), xdata(9, 2, 7) &
         /0.0E+0_wp, 2.06E+0_wp, -5.25E+0_wp/
      data &
         ydata(10, 1, 7), xdata(10, 1, 7), xdata(10, 2, 7) &
         /0.0E+0_wp, 1.74E+0_wp, -5.64E+0_wp/
      data &
         ydata(11, 1, 7), xdata(11, 1, 7), xdata(11, 2, 7) &
         /0.0E+0_wp, 1.34E+0_wp, -5.97E+0_wp/
      data &
         ydata(12, 1, 7), xdata(12, 1, 7), xdata(12, 2, 7) &
         /0.0E+0_wp, 0.90E+0_wp, -6.32E+0_wp/
      data &
         ydata(13, 1, 7), xdata(13, 1, 7), xdata(13, 2, 7) &
         /0.0E+0_wp, -0.28E+0_wp, -6.44E+0_wp/
      data &
         ydata(14, 1, 7), xdata(14, 1, 7), xdata(14, 2, 7) &
         /0.0E+0_wp, -0.78E+0_wp, -6.44E+0_wp/
      data &
         ydata(15, 1, 7), xdata(15, 1, 7), xdata(15, 2, 7) &
         /0.0E+0_wp, -1.36E+0_wp, -6.41E+0_wp/
      data &
         ydata(16, 1, 7), xdata(16, 1, 7), xdata(16, 2, 7) &
         /0.0E+0_wp, -1.90E+0_wp, -6.25E+0_wp/
      data &
         ydata(17, 1, 7), xdata(17, 1, 7), xdata(17, 2, 7) &
         /0.0E+0_wp, -2.50E+0_wp, -5.88E+0_wp/
      data &
         ydata(18, 1, 7), xdata(18, 1, 7), xdata(18, 2, 7) &
         /0.0E+0_wp, -2.88E+0_wp, -5.50E+0_wp/
      data &
         ydata(19, 1, 7), xdata(19, 1, 7), xdata(19, 2, 7) &
         /0.0E+0_wp, -3.18E+0_wp, -5.24E+0_wp/
      data &
         ydata(20, 1, 7), xdata(20, 1, 7), xdata(20, 2, 7) &
         /0.0E+0_wp, -3.44E+0_wp, -4.86E+0_wp/

      data &
         tdata(8) &
         /' BATES AND WATTS, 1988, TABLE A1.13, PAGES 280-281'/
      data &
         ndata(8), mdata(8), npdata(8), qdata(8) &
         /23, 1, 5, 2/
      data &
         (bdata(k, 8), k=1, 5) &
         /4.0E+0_wp, 2.0E+0_wp, 7.0E+0_wp, 0.40E+0_wp, 0.50E+0_wp/
      data &
         ydata(1, 1, 8), ydata(1, 2, 8), xdata(1, 1, 8) &
         /4.220E+0_wp, 0.136E+0_wp, 30.0E+0_wp/
      data &
         ydata(2, 1, 8), ydata(2, 2, 8), xdata(2, 1, 8) &
         /4.167E+0_wp, 0.167E+0_wp, 50.0E+0_wp/
      data &
         ydata(3, 1, 8), ydata(3, 2, 8), xdata(3, 1, 8) &
         /4.132E+0_wp, 0.188E+0_wp, 70.0E+0_wp/
      data &
         ydata(4, 1, 8), ydata(4, 2, 8), xdata(4, 1, 8) &
         /4.038E+0_wp, 0.212E+0_wp, 100.0E+0_wp/
      data &
         ydata(5, 1, 8), ydata(5, 2, 8), xdata(5, 1, 8) &
         /4.019E+0_wp, 0.236E+0_wp, 150.0E+0_wp/
      data &
         ydata(6, 1, 8), ydata(6, 2, 8), xdata(6, 1, 8) &
         /3.956E+0_wp, 0.257E+0_wp, 200.0E+0_wp/
      data &
         ydata(7, 1, 8), ydata(7, 2, 8), xdata(7, 1, 8) &
         /3.884E+0_wp, 0.276E+0_wp, 300.0E+0_wp/
      data &
         ydata(8, 1, 8), ydata(8, 2, 8), xdata(8, 1, 8) &
         /3.784E+0_wp, 0.297E+0_wp, 500.0E+0_wp/
      data &
         ydata(9, 1, 8), ydata(9, 2, 8), xdata(9, 1, 8) &
         /3.713E+0_wp, 0.309E+0_wp, 700.0E+0_wp/
      data &
         ydata(10, 1, 8), ydata(10, 2, 8), xdata(10, 1, 8) &
         /3.633E+0_wp, 0.311E+0_wp, 1000.0E+0_wp/
      data &
         ydata(11, 1, 8), ydata(11, 2, 8), xdata(11, 1, 8) &
         /3.540E+0_wp, 0.314E+0_wp, 1500.0E+0_wp/
      data &
         ydata(12, 1, 8), ydata(12, 2, 8), xdata(12, 1, 8) &
         /3.433E+0_wp, 0.311E+0_wp, 2000.0E+0_wp/
      data &
         ydata(13, 1, 8), ydata(13, 2, 8), xdata(13, 1, 8) &
         /3.358E+0_wp, 0.305E+0_wp, 3000.0E+0_wp/
      data &
         ydata(14, 1, 8), ydata(14, 2, 8), xdata(14, 1, 8) &
         /3.258E+0_wp, 0.289E+0_wp, 5000.0E+0_wp/
      data &
         ydata(15, 1, 8), ydata(15, 2, 8), xdata(15, 1, 8) &
         /3.193E+0_wp, 0.277E+0_wp, 7000.0E+0_wp/
      data &
         ydata(16, 1, 8), ydata(16, 2, 8), xdata(16, 1, 8) &
         /3.128E+0_wp, 0.255E+0_wp, 10000.0E+0_wp/
      data &
         ydata(17, 1, 8), ydata(17, 2, 8), xdata(17, 1, 8) &
         /3.059E+0_wp, 0.240E+0_wp, 15000.0E+0_wp/
      data &
         ydata(18, 1, 8), ydata(18, 2, 8), xdata(18, 1, 8) &
         /2.984E+0_wp, 0.218E+0_wp, 20000.0E+0_wp/
      data &
         ydata(19, 1, 8), ydata(19, 2, 8), xdata(19, 1, 8) &
         /2.934E+0_wp, 0.202E+0_wp, 30000.0E+0_wp/
      data &
         ydata(20, 1, 8), ydata(20, 2, 8), xdata(20, 1, 8) &
         /2.876E+0_wp, 0.182E+0_wp, 50000.0E+0_wp/
      data &
         ydata(21, 1, 8), ydata(21, 2, 8), xdata(21, 1, 8) &
         /2.838E+0_wp, 0.168E+0_wp, 70000.0E+0_wp/
      data &
         ydata(22, 1, 8), ydata(22, 2, 8), xdata(22, 1, 8) &
         /2.798E+0_wp, 0.153E+0_wp, 100000.0E+0_wp/
      data &
         ydata(23, 1, 8), ydata(23, 2, 8), xdata(23, 1, 8) &
         /2.759E+0_wp, 0.139E+0_wp, 150000.0E+0_wp/

      data &
         tdata(9) &
         /' ZWOLAK, WATSON, AND TYSON, 2004.'/
      data &
         ndata(9), mdata(9), npdata(9), qdata(9) &
         /4, 1, 2, 1/
      data &
         (bdata(k, 9), k=1, 2) &
         /200.0_wp, 5.0_wp/
      data &
         ydata(1, 1, 9), xdata(1, 1, 9) &
         /2.718281828459045_wp, 1.0_wp/
      data &
         ydata(2, 1, 9), xdata(2, 1, 9) &
         /7.389056098930650_wp, 2.0_wp/
      data &
         ydata(3, 1, 9), xdata(3, 1, 9) &
         /148.4131591025766_wp, 5.0_wp/
      data &
         ydata(4, 1, 9), xdata(4, 1, 9) &
         /403.4287934927353_wp, 6.0_wp/

      data &
         tdata(10) &
         /' ZWOLAK, WATSON, AND TYSON, 2005.'/
      data &
         ndata(10), mdata(10), npdata(10), qdata(10) &
         /4, 1, 2, 1/
      data &
         (bdata(k, 10), k=1, 2) &
         /200.0_wp, 5.0_wp/
      data &
         ydata(1, 1, 10), xdata(1, 1, 10) &
         /2.718281828459045_wp, 1.0_wp/
      data &
         ydata(2, 1, 10), xdata(2, 1, 10) &
         /7.389056098930650_wp, 2.0_wp/
      data &
         ydata(3, 1, 10), xdata(3, 1, 10) &
         /148.4131591025766_wp, 5.0_wp/
      data &
         ydata(4, 1, 10), xdata(4, 1, 10) &
         /403.4287934927353_wp, 6.0_wp/

      ! Variable definitions (alphabetically)
      !  BDATA:   The function parameter for each data set.
      !  BETA:    The function parameters.
      !  I:       An indexing variable.
      !  J:       An indexing variable.
      !  L:       An indexing variable.
      !  LDX:     The leading dimension of array X.
      !  M:       The number of columns of data in the explanatory variable.
      !  MDATA:   The number of columns of data in the explanatory variable in each data set.
      !  N:       The number of observations.
      !  NDATA:   The number of observations per data set.
      !  NP:      The number of function parameters.
      !  NPDATA:  The number of function parameters in each data set.
      !  QDATA:   The number of responses per observation in each data set.
      !  SETNO:   The number of the data set being analyzed.
      !  TDATA:   The reference for the each of the data sets.
      !  TITLE:   The reference for the data set being analyzed.
      !  X:       The explanatory variables.
      !  XDATA:   The explanatory variables for each data set.
      !  Y:       The response variable.
      !  YDATA:   The response variables for each data set.

      title = tdata(setno)

      n = ndata(setno)
      m = mdata(setno)
      np = npdata(setno)
      q = qdata(setno)

      if (wflat) then
         ldstpd = 1
         ldwd = 1
         ld2wd = 1
         ldwe = 1
         ld2we = 1
      else
         ldstpd = n
         ldwd = n
         ld2wd = m
         ldwe = n
         ld2we = q
      end if

      if (allocated(x)) deallocate (x)
      if (allocated(y)) deallocate (y)
      if (allocated(beta)) deallocate (beta)
      if (allocated(lower)) deallocate (lower)
      if (allocated(upper)) deallocate (upper)
      if (allocated(delta)) deallocate (delta)
      if (allocated(ifixb)) deallocate (ifixb)
      if (allocated(ifixx)) deallocate (ifixx)
      if (allocated(wd)) deallocate (wd)
      if (allocated(we)) deallocate (we)
      if (allocated(sclb)) deallocate (sclb)
      if (allocated(scld)) deallocate (scld)
      if (allocated(stpb)) deallocate (stpb)
      if (allocated(stpd)) deallocate (stpd)

      allocate (x(n, m), y(n, q), beta(np), lower(np), upper(np), delta(n, m), &
                ifixb(np), ifixx(n, m), wd(ldwd, ld2wd, m), we(ldwe, ld2we, q), &
                sclb(np), scld(n, m), stpb(np), stpd(ldstpd, m))

      x = xdata(1:n, 1:m, setno)
      y = ydata(1:n, 1:q, setno)
      beta = bdata(1:np, setno)
      lower = -huge(1.0_wp)
      upper = huge(1.0_wp)
      delta = 0.0_wp
      ifixb = -1
      ifixx = -1
      sclb = -1.0_wp
      scld = -1.0_wp
      stpb = -1.0_wp
      stpd = -1.0_wp
      we = -1.0_wp
      wd = -1.0_wp

   end subroutine set_inputs

   pure subroutine fcn(beta, xplusd, ifixb, ifixx, ideval, f, fjacb, fjacd, istop)
   !! Compute model function and jacobian for odrpack exerciser.

      use odrpack_kinds, only: zero, one

      real(wp), intent(in) :: beta(:)
      real(wp), intent(in) :: xplusd(:, :)
      integer, intent(in) :: ifixb(:)
      integer, intent(in) :: ifixx(:, :)
      integer, intent(in) :: ideval
      real(wp), intent(out) :: f(:, :)
      real(wp), intent(out) :: fjacb(:, :, :)
      real(wp), intent(out) :: fjacd(:, :, :)
      integer, intent(out) :: istop

      ! Local scalars
      real(wp) :: ctheta, fac1, fac2, fac3, fac4, freq, omega, phi, pi, r, stheta, theta
      integer :: i, j, k, n

      ! Variable definitions (alphabetically)
      !  BETA:    Current values of parameters
      !  F:       Predicted function values
      !  FAC1:    A factors or terms used in computing the jacobians.
      !  FAC2:    A factors or terms used in computing the jacobians.
      !  FAC3:    A factors or terms used in computing the jacobians.
      !  FAC4:    A factors or terms used in computing the jacobians.
      !  FJACB:   Jacobian with respect to BETA
      !  FJACD:   Jacobian with respect to errors DELTA
      !  IDEVAL:  Indicator for selecting computation to be performed
      !  IFIXB:   Indicators for "fixing" parameters (BETA)
      !  IFIXX:   Indicators for "fixing" explanatory variable (X)
      !  LDIFX:   Leading dimension of array IFIXX
      !  ISTOP:   Stopping condition, where
      !           0 means current BETA and X+DELTA were
      !           acceptable and values were computed successfully
      !           1 means current BETA and X+DELTA are
      !           not acceptable;  ODRPACK95 should select
      !           values closer to most recently used values
      !           -1 means current BETA and X+DELTA are not acceptable; ODRPACK95 should stop
      !  LDN:     Leading dimension declarator equal or exceeding N
      !  LDM:     Leading dimension declarator equal or exceeding M
      !  LDNP:    Leading dimension declarator equal or exceeding NP
      !  M:       The number of columns of data in the explanatory variable.
      !  N:       The number of observations.
      !  NP:      The number of function parameters.
      !  Q:       The number of responses per observation.
      !  SETNO:   The number of the data set being analyzed.
      !  XPLUSD:  Current value of explanatory variable, i.e., X + DELTA

      !  Check for BETA outside bounds.  Return with error if BETA outside bounds.

      if (any(lower > beta)) then
         istop = -1
         return
      end if

      if (any(upper < beta)) then
         istop = -2
         return
      end if

      n = ubound(xplusd, 1)

      if (setno == 1) then

         !  Setno. 1:  Boggs, Byrd and Schnabel, 1985, example 1

         if (beta(1) <= 1.01E0_wp) then
            istop = 0

            if (mod(ideval, 10) /= 0) then
               do i = 1, n
                  f(i, 1) = beta(1)/(xplusd(i, 1) - beta(2))
               end do
            end if

            if (mod(ideval/10, 10) /= 0) then
               do i = 1, n
                  fjacb(i, 1, 1) = 1/(xplusd(i, 1) - beta(2))
                  fjacb(i, 2, 1) = beta(1)*(xplusd(i, 1) - beta(2))**(-2)
               end do
            end if

            if (mod(ideval/100, 10) /= 0) then
               do i = 1, n
                  fjacd(i, 1, 1) = -beta(1)*(xplusd(i, 1) - beta(2))**(-2)
               end do
            end if

         else
            istop = 1
         end if

      elseif (setno == 2) then

         !  Setno. 2:  Boggs, Byrd and Schnabel, 1985, example 2

         istop = 0

         do i = 1, n
            fac1 = (beta(2)*xplusd(i, 1) + beta(3)*xplusd(i, 2) - one)

            if (mod(ideval, 10) /= 0) then
               f(i, 1) = beta(1)/fac1
            end if

            if (mod(ideval/10, 10) /= 0) then
               fjacb(i, 1, 1) = 1/fac1
               fjacb(i, 2, 1) = -beta(1)*(fac1**(-2))*xplusd(i, 1)
               fjacb(i, 3, 1) = -beta(1)*(fac1**(-2))*xplusd(i, 2)
            end if

            if (mod(ideval/100, 10) /= 0) then
               fjacd(i, 1, 1) = -beta(1)*(fac1**(-2))*beta(2)
               fjacd(i, 2, 1) = -beta(1)*(fac1**(-2))*beta(3)
            end if
         end do

      elseif (setno == 3) then

         !  Setno. 3:  Boggs, Byrd and Schnabel, 1985, example 3

         istop = 0

         if (mod(ideval, 10) /= 0) then
            do i = 1, n
               f(i, 1) = zero
               do j = 1, 4
                  f(i, 1) = f(i, 1) + beta(j)/(xplusd(i, 1) + beta(j + 5))
               end do
               f(i, 1) = f(i, 1) + beta(5)
            end do
         end if

         if (mod(ideval/10, 10) /= 0) then
            do i = 1, n
               fjacb(i, 5, 1) = one
               do k = 1, 4
                  fjacb(i, k, 1) = 1/(xplusd(i, 1) + beta(k + 5))
                  fjacb(i, k + 5, 1) = -beta(k)*(xplusd(i, 1) + beta(k + 5))**(-2)
               end do
            end do
         end if

         if (mod(ideval/100, 10) /= 0) then
            do i = 1, n
               fjacd(i, 1, 1) = zero
               do k = 4, 1, -1
                  fjacd(i, 1, 1) = fjacd(i, 1, 1) - beta(k)*(xplusd(i, 1) + beta(k + 5))**(-2)
               end do
            end do
         end if

      elseif (setno == 4) then

         !  Setno. 4:  Himmelblau, 1970, example 6.2-4, page 188

         istop = 0

         if (mod(ideval, 10) /= 0) then
            do i = 1, n
               f(i, 1) = beta(1)*xplusd(i, 1) + &
                         beta(2)*exp(beta(3)*xplusd(i, 2))
            end do
         end if

         if (mod(ideval/10, 10) /= 0) then
            do i = 1, n
               fjacb(i, 1, 1) = xplusd(i, 1)
               fjacb(i, 2, 1) = exp(beta(3)*xplusd(i, 2))
               fjacb(i, 3, 1) = beta(2)*exp(beta(3)*xplusd(i, 2))*xplusd(i, 2)
            end do
         end if

         if (mod(ideval/100, 10) /= 0) then
            do i = 1, n
               fjacd(i, 1, 1) = beta(1)
               fjacd(i, 2, 1) = beta(2)*exp(beta(3)*xplusd(i, 2))*beta(3)
            end do
         end if

      elseif (setno == 5) then

         !  Setno. 5:  Draper and Smith, 1981, exercise i, page 521-522

         istop = 0

         if (mod(ideval, 10) /= 0) then
            do i = 1, n
               f(i, 1) = exp(-beta(1)*xplusd(i, 1)* &
                             exp(-beta(2)*(1/xplusd(i, 2) - 1/620.0E0_wp)))
            end do
         end if

         if (mod(ideval/10, 10) /= 0) then
            do i = 1, n
               fac1 = 1/xplusd(i, 2) - 1/620.0E0_wp
               fac2 = exp(-beta(2)*fac1)
               fac3 = beta(1)*xplusd(i, 1)
               fac4 = exp(-fac3*fac2)

               fjacb(i, 1, 1) = -fac4*xplusd(i, 1)*fac2
               fjacb(i, 2, 1) = fac4*fac3*fac2*fac1

               if (mod(ideval/100, 10) /= 0) then
                  fjacd(i, 1, 1) = -fac4*beta(1)*fac2
                  fjacd(i, 2, 1) = -fac4*fac3*fac2* &
                                   beta(2)/xplusd(i, 2)**2
               end if
            end do
         end if

      elseif (setno == 6) then

         !  Setno. 6:  Powell and Macdonald, 1972, tables 7 and 8, page 153-154
         !          N.B.  this derivative is intentionally coded incorrectly

         istop = 0

         if (mod(ideval, 10) /= 0) then
            do i = 1, n
               f(i, 1) = beta(1)*(one + beta(3)*xplusd(i, 1)/beta(2))**(-1/beta(3))
            end do
         end if

         if (mod(ideval/10, 10) /= 0) then
            do i = 1, n
               fjacb(i, 1, 1) = zero
               fjacb(i, 2, 1) = zero
               fjacb(i, 3, 1) = zero

               if (mod(ideval/100, 10) /= 0) then
                  fjacd(i, 1, 1) = xplusd(i, 1)
               end if
            end do
         end if

      elseif (setno == 7) then

         !  Setno. 7:  Fuller, 1987, table 3.2.10, pages 244-245
         !          N.B.  this derivative is intentionally coded incorrectly

         istop = 0

         if (mod(ideval, 10) /= 0) then
            do i = 1, n
               f(i, 1) = beta(3)*(xplusd(i, 1) - beta(1))**2 + &
                         2*beta(4)*(xplusd(i, 1) - beta(1))* &
                         (xplusd(i, 2) - beta(2)) + &
                         beta(5)*(xplusd(i, 2) - beta(2))**2 - 1.0E0_wp
            end do
         end if

         if (mod(ideval/10, 10) /= 0) then
            do i = 1, n
               fjacb(i, 1, 1) = zero
               fjacb(i, 2, 1) = zero
               fjacb(i, 3, 1) = zero
               fjacb(i, 4, 1) = zero
               fjacb(i, 5, 1) = zero

               if (mod(ideval/100, 10) /= 0) then
                  fjacd(i, 1, 1) = zero
                  fjacd(i, 2, 1) = zero
               end if
            end do
         end if

      elseif (setno == 8) then

         !  Setno. 8:  Bates and Watts, 1988, table A1.13, pages 280-281
         !          N.B.  This derivative is intentionally coded incorrectly

         do i = 1, n
            if (xplusd(i, 1) < 0.0E0_wp) then
               istop = 1
               return
            end if
         end do
         istop = 0

         if (mod(ideval, 10) /= 0) then
            pi = 3.141592653589793238462643383279E0_wp
            theta = pi*beta(4)*0.5E0_wp
            ctheta = cos(theta)
            stheta = sin(theta)
            do i = 1, n
               freq = xplusd(i, 1)
               omega = (2.0E0_wp*pi*freq*exp(-beta(3)))**beta(4)
               phi = atan2((omega*stheta), (1 + omega*ctheta))
               r = (beta(1) - beta(2))* &
                   sqrt((1 + omega*ctheta)**2 + (omega*stheta)**2)** &
                   (-beta(5))
               f(i, 1) = beta(2) + r*cos(beta(5)*phi)
               f(i, 2) = r*sin(beta(5)*phi)
            end do
         end if

         if (mod(ideval/10, 10) /= 0) then
            do i = 1, n
               fjacb(i, 1, 1) = zero
               fjacb(i, 2, 1) = zero
               fjacb(i, 3, 1) = zero
               fjacb(i, 4, 1) = zero
               fjacb(i, 5, 1) = zero

               fjacb(i, 1, 2) = zero
               fjacb(i, 2, 2) = zero
               fjacb(i, 3, 2) = zero
               fjacb(i, 4, 2) = zero
               fjacb(i, 5, 2) = zero

               if (mod(ideval/100, 10) /= 0) then
                  fjacd(i, 1, 1) = zero
                  fjacd(i, 1, 2) = zero
               end if
            end do
         end if

      elseif (setno == 9) then

         !  Setno. 9:  Zwolak, Watson, and Tyson, 2004.

         istop = 0

         if (mod(ideval, 10) /= 0) then
            do i = 1, n
               f(i, 1) = beta(1)*exp(beta(2)*xplusd(i, 1))
            end do
         end if

         if (mod(ideval/10, 10) /= 0) then
            do i = 1, n
               fjacb(i, 1, 1) = exp(beta(2)*xplusd(i, 1))
               fjacb(i, 2, 1) = beta(1)*xplusd(i, 1)*exp(beta(2)*xplusd(i, 1))
            end do
         end if

         if (mod(ideval/100, 10) /= 0) then
            do i = 1, n
               fjacd(i, 1, 1) = beta(1)*beta(2)*exp(beta(2)*xplusd(i, 1))
            end do
         end if

      elseif (setno == 10) then

         !  Setno. 10:  Zwolak, Watson, and Tyson, 2005.

         istop = 0

         if (mod(ideval, 10) /= 0) then
            do i = 1, n
               f(i, 1) = beta(1)/2.0_wp*exp(beta(2)*xplusd(i, 1))
            end do
         end if

         if (mod(ideval/10, 10) /= 0) then
            do i = 1, n
               fjacb(i, 1, 1) = exp(beta(2)*xplusd(i, 1))/2.0_wp
               fjacb(i, 2, 1) = beta(1)/2.0_wp*xplusd(i, 1)*exp(beta(2)*xplusd(i, 1))
            end do
         end if

         if (mod(ideval/100, 10) /= 0) then
            do i = 1, n
               fjacd(i, 1, 1) = beta(1)/2.0_wp*beta(2)*exp(beta(2)*xplusd(i, 1))
            end do
         end if

      end if

   end subroutine fcn

   subroutine loc_wsse(n, m, np, q, ldwe, ld2we, isodr, wssi, wssdeli, wssepsi)
   !! Locations of the weighted sum of squares results in the `rwork` array.

      use odrpack_core, only: loc_rwork

      integer, intent(in) :: n, m, np, q, ldwe, ld2we
      logical, intent(in) :: isodr
      integer, intent(out) :: wssi, wssdeli, wssepsi

      integer :: deltai, epsi, xplusdi, fni, sdi, vcvi, rvari, &
                 rcondi, etai, olmavgi, taui, alphai, actrsi, pnormi, rnormsi, prersi, partoli, &
                 sstoli, taufaci, epsmaci, beta0i, betaci, betasi, betani, si, ssi, ssfi, &
                 qrauxi, ui, fsi, fjacbi, we1i, diffi, deltasi, deltani, ti, tti, omegai, &
                 fjacdi, wrk1i, wrk2i, wrk3i, wrk4i, wrk5i, wrk6i, wrk7i, loweri, upperi, &
                 lwkmin

      call loc_rwork(n, m, q, np, ldwe, ld2we, isodr, &
                     deltai, epsi, xplusdi, fni, sdi, vcvi, &
                     rvari, wssi, wssdeli, wssepsi, rcondi, etai, &
                     olmavgi, taui, alphai, actrsi, pnormi, rnormsi, prersi, &
                     partoli, sstoli, taufaci, epsmaci, &
                     beta0i, betaci, betasi, betani, si, ssi, ssfi, qrauxi, ui, &
                     fsi, fjacbi, we1i, diffi, &
                     deltasi, deltani, ti, tti, omegai, fjacdi, &
                     wrk1i, wrk2i, wrk3i, wrk4i, wrk5i, wrk6i, wrk7i, &
                     loweri, upperi, &
                     lwkmin)

   end subroutine loc_wsse

end module test_solution_m

program test_solution
!! Solution tests for [[odrpack]].

   use odrpack_kinds, only: wp, one
   use test_solution_m, only: test_odr
   implicit none

   ! Local scalars
   real(wp) :: tstfac
   integer :: ntest, lunerr, lunrpt, lunsum
   logical :: passed

   ! Variable declarations (alphabetically)
   !  LUNERR:  The logical unit number used for error messages.
   !  LUNRPT:  The logical unit number used for computation reports.
   !  LUNSUM:  The logical unit number used for a summary report listing
   !           only the test comparisons and not the odrpack generated reports.
   !  NTEST:   The number of tests to be run.
   !  PASSED:  The variable designating whether the results of all of the tests agree with
   !           those from the cray ymp using double precision (PASSED=TRUE), or whether some
   !           of the results disagreed (PASSED=FALSE).
   !  TSTFAC:  The user-supplied factor for scaling the test tolerances used to check for
   !           agreement between computed results and results obtained using REAL (wp) version
   !           on cray YMP. Values of TSTFAC greater than one increase the test tolerances,
   !           making the tests easier to pass and allowing small discrepancies between the
   !           computed and expected results to be automatically discounted.

   ! Set up necessary files
   open (newunit=lunrpt, file='./test/test_solution_report.txt')
   open (newunit=lunerr, file='./test/test_solution_error.txt')
   open (newunit=lunsum, file='./test/test_solution_summary.txt')

   ! Exercise REAL (wp) version of ODRPACK95
   tstfac = one
   ntest = 23
   passed = .false.
   call test_odr(ntest, tstfac, passed, lunrpt, lunerr, lunsum)

   if (passed) then
      stop "Solution tests passed. See report."
   else
      error stop "Solution tests failed. See report."
   end if

   close (lunrpt)
   close (lunerr)
   close (lunsum)

end program test_solution
