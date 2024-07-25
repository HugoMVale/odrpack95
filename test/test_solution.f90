program dtest
!***Begin Prologue  DTEST
!***Refer to ODR
!***Routines Called  DODRX
!***Date Written   861229   (YYMMDD)
!***Revision Date  920619   (YYMMDD)
!***Purpose  EXERCISE FEATURES OF ODRPACK95 SOFTWARE
!***End Prologue  DTEST

!...Used modules
   use odrpack_kinds, only: wp
   implicit none

!...Scalars in common
   integer :: ntest

!...Local scalars
   real(wp) :: tstfac
   integer :: lunerr, lunrpt, lunsum
   logical :: passed

!...External subroutines
   external :: dodrx

!...Common blocks
   common/tstset/ntest

!***Variable declarations (alphabetically)
!
!       LUNERR:  The logical unit number used for error messages.
!       LUNRPT:  The logical unit number used for computation reports.
!       LUNSUM:  The logical unit number used for a summary report listing
!       only the test comparisons and not the odrpack generated
!       reports.
!       NTEST:   The number of tests to be run.
!       PASSED:  The variable designating whether the results of all of the
!       tests agree with those from the cray ymp using double
!       precision (PASSED=TRUE), or whether some of the results
!       disagreed (PASSED=FALSE).
!       TSTFAC:  The user-supplied factor for scaling the test tolerances
!       used to check for agreement between computed results and
!       results obtained using REAL (wp) version on cray
!       YMP.  Values of TSTFAC greater than one increase the
!       test tolerances, making the tests easier to pass and
!       allowing small discrepancies between the computed and
!       expected results to be automatically discounted.
!
!
!***First executable statement  TEST
!
!
!  Set up necessary files
!
!  NOTE:  ODRPACK95 generates computation and error reports on
!       logical unit 6 by default;
!       logical unit 'LUNSUM' used to summarize results of comparisons
!       from exercise routine DODRX.

   lunrpt = 18
   lunerr = 18
   lunsum = 19

   open (unit=lunrpt, file='./test/test_solution_report.txt')
   open (unit=lunerr, file='./test/test_solution_report.txt')
   open (unit=lunsum, file='./test/test_solution_summary.txt')

!  Exercise REAL (wp) version of ODRPACK95
!  (test reports generated on file 'RESULTS' and
!       summarized in file 'SUMMARY')

   ntest = 23
   tstfac = 1.0E0_wp
   call dodrx(tstfac, passed, lunsum)

end program dtest

subroutine dodrx(tstfac, passed, lunsum)
!***Begin Prologue  DODRX
!***Refer to ODR
!***Routines Called  DDOT,DNRM2,ODR,DODRXD,
!       DODRXF,DODRXW,DWGHT
!***Date Written   860529   (YYMMDD)
!***Revision Date  920619   (YYMMDD)
!***Purpose  Exercise features of ODRPACK95 software
!***End Prologue  DODRX

!...Used modules
   use odrpack, only: odr
   use odrpack_kinds, only: wp
   use odrpack_core, only: dwght
   implicit none

!...Parameters
   real(wp), parameter :: base = radix(1.0E0_wp)
   integer, parameter :: maxn = 50, maxm = 3, maxnp = 10, maxnq = 2, ntests = 23, &
                         ldwe = maxn, ld2we = maxnq, ldwd = maxn, ld2wd = maxm, &
                         lwork = 18 + 11*maxnp + maxnp**2 + maxm + maxm**2 + &
                         4*maxn*maxnq + 6*maxn*maxm + 2*maxn*maxnq*maxnp + &
                         2*maxn*maxnq*maxm + maxnq**2 + &
                         5*maxnq + maxnq*(maxnp + maxm) + ldwe*ld2we*maxnq, &
                         liwork = 20 + maxnp + maxnq*(maxnp + maxm)

!...Scalar arguments
   real(wp) :: tstfac
   integer :: lunsum
   logical :: passed

!...Scalars in common
   integer :: ntest, setno

!...Local scalars
   integer :: i, info, iprint, itest, job, l, ldifx, ldscld, ldstpd, ldwd1, ldwe1, &
              ldx, ldy, ld2wd1, ld2we1, liwmin, lun, lunerr, lunrpt, lwmin, &
              m, maxit, msg, n, ndigit, np, nq
   real(wp) :: bnrm, epsmac, ewrt, ewrt2, hundrd, one, p01, p2, partol, sstol, &
               taufac, three, tsttol, two, wss, wssdel, wsseps, zero
   logical :: failed, fails, isodr, short
   character(len=80) :: title

!...Arrays in common
   real(wp) :: lower(maxnp), upper(maxnp)

!...Local arrays
   real(wp) :: beta(maxnp), dpymp(2, ntests), &
               sclb(maxnp), scld(maxn, maxm), &
               stpb(maxnp), stpd(maxn, maxm), &
               we(maxn, maxnq, maxnq), wd(maxn, maxm, maxm), &
               wrk(maxn*maxm + maxn*maxnq), x(maxn, maxm), y(maxn, maxnq), &
               tempretl(maxn, maxm)
   real(wp), pointer :: work(:)
   real(wp), allocatable :: delta(:, :)
   integer :: idpymp(ntests), ifixb(maxnp), ifixx(maxn, maxm)
   integer, pointer :: iwork(:)

!...External functions
   real(wp), external :: ddot, dnrm2

!...External subroutines
   external :: dodrxd, dodrxf, dodrxw

!...Common blocks
   common/setid/setno
   common/tstset/ntest
   common/bounds/lower, upper

!...Data statements
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

!...Routine names used as subprogram arguments
!       DODRXF:  The user-supplied routine for evaluating the model.

!...Variable definitions (alphabetically)
!       BASE:    The base of floating point numbers on the current machine
!       BETA:    The function parameters.
!       BNRM:    The norm of BETA.
!       DELTA:   The error in the X data.
!       DPYMP:   The floating point results from a cray YMP using
!       REAL (wp).
!       EPSMAC:  The value of machine precision.
!       EWRT:    A temporary variable for the denominator of the relative error
!       calculations (error with respect to).
!       EWRT2:   A temporary variable for the denominator of the relative error
!       calculations (error with respect to).
!       FAILED:  The variable designating whether the results of all of the
!       demonstration runs agreed with those from the cray YMP
!       using REAL (wp) (FAILED=FALSE) or whether some of
!       the tests disagreed (FAILED=TRUE).
!       FAILS:   The variable designating whether the results of an
!       individual demonstration run agreed with those from the
!       cray YMP using REAL (wp) (FAILS=FALSE) or
!       disagree (FAILS=TRUE).
!       HUNDRD:  The value 100.0E0_wp.
!       I:       An index variable.
!       IDPYMP:  The integer results from a cray YMP using
!       REAL (wp).
!       IFIXB:   The values designating whether the elements of BETA are
!       fixed at their input values or not.
!       IFIXX:   The values designating whether the elements of DELTA are
!       fixed at their input values or not.
!       INFO:    The variable designating why the computations stopped.
!       IPRINT:  The print control variable.
!       ISODR:   The variable designating whether the solution is by odr
!       (ISODR=TRUE) or by ols (ISODR=FALSE).
!       ITEST:   The number of the current test being run.
!       IWORK:   The integer work space.
!       J:       An index variable.
!       JOB:     The variable controlling problem initialization and
!       computational method.
!       LDIFX:   The leading dimension of array IFIXX.
!       LDSCLD:  The leading dimension of array SCLD.
!       LDWD:    The leading dimension of array WD.
!       LDWD1:   The leading dimension of array WD as passed to ODRPACK95.
!       LDWE:    The leading dimension of array WE.
!       LDWE1:   The leading dimension of array WE as passed to ODRPACK95.
!       LDX:     The leading dimension of array X.
!       LDY:     The leading dimension of array Y.
!       LD2WD:   The second dimension of array WD.
!       LD2WD1:  The second dimension of array WD as passed to ODRPACK95.
!       LD2WE:   The second dimension of array WE.
!       LD2WE1:  The second dimension of array WE as passed to ODRPACK95.
!       LIWKMN:  The minimum acceptable length of array IWORK.
!       LIWMIN:  The minimum length of vector IWORK for a given problem.
!       LIWORK:  The length of vector IWORK.
!       LUN:     The logical unit number currently being used.
!       LUNERR:  The logical unit number used for error messages.
!       LUNRPT:  The logical unit number used for computation reports.
!       LUNSUM:  The logical unit number used for a summary report.
!       LWKMN:   The minimum acceptable length of array WORK.
!       LWMIN:   The minimum length of vector WORK for a given problem.
!       LWORK:   The length of vector WORK.
!       M:       The number of columns of data in the explanatory variable.
!       MAXIT:   The maximum number of iterations allowed.
!       MSG:     The variable designating which message is to be printed as
!       a result of the comparison with the cray YMP or x86 (Linux)
!       results.
!       N:       The number of observations.
!       NDIGIT:  The number of accurate digits in the function results, as
!       supplied by the user.
!       NP:      The number of function parameters.
!       NTEST:   The number of tests to be run.
!       NTESTS:  The number of different tests available.
!       ONE:     The value 1.0E0_wp.
!       PASSED:  The variable designating whether the results of all of the
!       demonstration runs agreed with those from the cray YMP
!       using REAL (wp) (PASSED=TRUE), or whether some of
!       the results disagreed (PASSED=FALSE).
!       P01:     The value 0.01E0_wp.
!       P2:      The value 0.2E0_wp.
!       PARTOL:  The parameter convergence stopping criteria.
!       SCLB:    The scaling values for BETA.
!       SCLD:    The scaling values for DELTA.
!       SETNO:   The number of the data set being analyzed.
!       SHORT:   The variable designating whether ODRPACK95 is invoked by the
!       short-call (SHORT=.TRUE.) or the long-call (SHORT=.FALSE.).
!       SSTOL:   The sum-of-squares convergence stopping tolerance.
!       TAUFAC:  The factor used to compute the initial trust region
!       diameter.
!       THREE:   The value 3.0E0_wp.
!       TITLE:   The reference for the data set being analyzed.
!       TSTFAC:  The user-supplied factor for scaling the test tolerances
!       used to check for agreement between computed results and
!       results obtained using REAL (wp) version on cray
!       YMP.
!       TSTTOL:  The test tolerance used in checking computed values for
!       purposes of determining proper installation.
!       TWO:     The value 2.0E0_wp.
!       WD:      The DELTA weights.
!       WE:      The EPSILON weights.
!       WORK:    The REAL (wp) work space.
!       WRK:     The REAL (wp) work space for computing test results.
!       WSS:     The sum of the squared weighted errors.
!       WSSDEL:  The sum of the squared weighted errors in X.
!       WSSEPS:  The sum of the squared weighted errors in Y.
!       X:       The explanatory variable.
!       Y:       The response variable.
!       ZERO:    The value 0.0E0_wp.

!***First executable statement  DODRX

!  Allocate work arrays and DELTA

   allocate (delta(maxn, maxm), iwork(liwork), work(lwork))

!  Set logical units for error and computation reports

   lunerr = 18
   lunrpt = 18

!  Initialize test tolerance

   if (tstfac .gt. one) then
      tsttol = tstfac
   else
      tsttol = one
   end if

!  Initialize machine precision

   epsmac = base**(1 - digits(base))

!  Initialize leading dimension of X

   ldx = maxn
   ldy = maxn

!  Initialize miscellaneous variables used in the exercise procedure

   failed = .false.
   short = .true.
   isodr = .true.
   n = 0

!  Begin exercising ODRPACK95

   test: do itest = 1, ntest

!  Set control values to invoke default values

      we(1, 1, 1) = -one
      ldwe1 = ldwe
      ld2we1 = ld2we
      wd(1, 1, 1) = -one
      ldwd1 = ldwd
      ld2wd1 = ld2wd

      ifixb(1) = -1
      ifixx(1, 1) = -1
      ldifx = maxn

      ndigit = -1
      taufac = -one

      sstol = -one
      partol = -one
      maxit = -1

      iprint = 2112
      ! iprint = 6616

      stpb(1) = -one
      stpd(1, 1) = -one
      ldstpd = 1

      sclb(1) = -one
      scld(1, 1) = -one
      ldscld = 1

      upper(:) = huge(one)
      lower(:) = -huge(one)

      if (itest .eq. 1) then

!  Test simple odr problem with analytic derivatives.

         lun = lunrpt
         write (lun, 1000)
         do i = 1, 2
            write (lun, 1001) itest
            write (lun, 1010)
            lun = lunsum
         end do
         setno = 5
         call dodrxd(title, n, m, np, nq, ldx, x, ldy, y, beta)
         work = zero
         delta = zero
         job = 00020
         short = .true.
         isodr = .true.

      elseif (itest .eq. 2) then

!  Test simple ols problem with forward difference derivatives.

         lun = lunrpt
         write (lun, 1000)
         do i = 1, 2
            write (lun, 1001) itest
            write (lun, 1020)
            lun = lunsum
         end do
         setno = 5
         call dodrxd(title, n, m, np, nq, ldx, x, ldy, y, beta)
         work = zero
         delta = zero
         job = 00002
         short = .true.
         isodr = .false.

      elseif (itest .eq. 3) then

!  Test parameter fixing capabilities for poorly scaled ols problem
!  with analytic derivatives.
!  (derivative checking turned off.)

         lun = lunrpt
         write (lun, 1000)
         do i = 1, 2
            write (lun, 1001) itest
            write (lun, 1030)
            lun = lunsum
         end do
         setno = 3
         call dodrxd(title, n, m, np, nq, ldx, x, ldy, y, beta)
         work = zero
         delta = zero
         ifixb(1) = 1
         ifixb(2) = 1
         ifixb(3) = 1
         ifixb(4) = 0
         ifixb(5) = 1
         ifixb(6) = 0
         ifixb(7) = 0
         ifixb(8) = 0
         ifixb(9) = 0
         job = 00042
         short = .false.
         isodr = .false.

      elseif (itest .eq. 4) then

!  Test weighting capabilities for odr problem with
!  analytic derivatives.
!  Also shows solution of poorly scaled odr problem.
!  (derivative checking turned off.)
!  N.B., this run continues from where test 3 left off.

         lun = lunrpt
         write (lun, 1000)
         do i = 1, 2
            write (lun, 1001) itest
            write (lun, 1040)
            lun = lunsum
         end do
         setno = 3
         work = zero
         delta = zero
         ldwd1 = ldwd
         ldwe1 = ldwe
         ld2wd1 = ld2wd
         ld2we1 = ld2we
         do i = 1, n
            wd(i, 1, 1) = (p01/abs(x(i, 1)))**2
            we(i, 1, 1) = one
         end do
         we(28, 1, 1) = zero
         ifixb(1) = 1
         ifixb(2) = 1
         ifixb(3) = 1
         ifixb(4) = 0
         ifixb(5) = 1
         ifixb(6) = 1
         ifixb(7) = 1
         ifixb(8) = 0
         ifixb(9) = 0
         job = 00030
         iprint = 2232
         short = .false.
         isodr = .true.

      elseif (itest .eq. 5) then

!  Test DELTA initialization capabilities and user-supplied scaling
!  and use of istop to restrict parameter values
!  for odr problem with analytic derivatives.

         lun = lunrpt
         write (lun, 1000)
         do i = 1, 2
            write (lun, 1001) itest
            write (lun, 1050)
            lun = lunsum
         end do
         setno = 1
         call dodrxd(title, n, m, np, nq, ldx, x, ldy, y, beta)
         work = zero
         delta = zero
         job = 01020
         ldscld = 1
         scld(1, 1) = two
         sclb(1) = p2
         sclb(2) = one
         ldwe1 = 1
         ld2we1 = 1
         we(1, 1, 1) = -one
         ldwd1 = 1
         ld2wd1 = 1
         wd(1, 1, 1) = -one
         do i = 20, 21
            delta(i, 1) = beta(1)/y(i, 1) + beta(2) - x(i, 1)
         end do
         short = .false.
         isodr = .true.

      elseif (itest .eq. 6) then

!  Test stiff stopping conditions for unscaled odr problem  with analytic derivatives.

         lun = lunrpt
         write (lun, 1000)
         do i = 1, 2
            write (lun, 1001) itest
            write (lun, 1060)
            lun = lunsum
         end do
         setno = 4
         call dodrxd(title, n, m, np, nq, ldx, x, ldy, y, beta)
         work = zero
         delta = zero
         job = 00020
         sstol = hundrd*epsmac
         partol = epsmac
         maxit = 2
         short = .false.
         isodr = .true.

      elseif (itest .eq. 7) then

!  Test restart for unscaled odr problem with analytic derivatives.

         lun = lunrpt
         write (lun, 1000)
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
         short = .false.
         isodr = .true.

      elseif (itest .eq. 8) then

!  Test use of TAUFAC to restrict first step
!  for odr problem with central difference derivatives.

         lun = lunrpt
         write (lun, 1000)
         do i = 1, 2
            write (lun, 1001) itest
            write (lun, 1080)
            lun = lunsum
         end do
         setno = 6
         call dodrxd(title, n, m, np, nq, ldx, x, ldy, y, beta)
         work = zero
         delta = zero
         job = 00210
         taufac = p01
         short = .false.
         isodr = .true.

      elseif (itest .eq. 9) then

!  Test implicit odr problem
!  with forward finite difference derivatives
!  and covariance matrix constructed with recomputed derivatives.

         lun = lunrpt
         write (lun, 1000)
         do i = 1, 2
            write (lun, 1001) itest
            write (lun, 1090)
            lun = lunsum
         end do
         setno = 7
         call dodrxd(title, n, m, np, nq, ldx, x, ldy, y, beta)
         work = zero
         delta = zero
         job = 00001
         partol = epsmac**(one/three)
         short = .true.
         isodr = .true.

      elseif (itest .eq. 10) then

!  Test multiresponse odr problem
!  with central difference derivatives ,
!  DELTA initialized to nonzero values,
!  variable fixing, and weighting.

         lun = lunrpt
         write (lun, 1000)
         do i = 1, 2
            write (lun, 1001) itest
            write (lun, 1100)
            lun = lunsum
         end do
         setno = 8
         call dodrxd(title, n, m, np, nq, ldx, x, ldy, y, beta)
         work = zero
         delta = zero

         ldwd1 = ldwd
         ldwe1 = ldwe
         ld2wd1 = ld2wd
         ld2we1 = ld2we
         do i = 1, n
!  Initialize DELTA, and specify first decade of frequencies as fixed
            if (x(i, 1) .lt. 100.0E0_wp) then
               delta(i, 1) = 0.0E0_wp
               ifixx(i, 1) = 0
            elseif (x(i, 1) .le. 150.0E0_wp) then
               delta(i, 1) = 0.0E0_wp
               ifixx(i, 1) = 1
            elseif (x(i, 1) .le. 1000.0E0_wp) then
               delta(i, 1) = 25.0E0_wp
               ifixx(i, 1) = 1
            elseif (x(i, 1) .le. 10000.0E0_wp) then
               delta(i, 1) = 560.0E0_wp
               ifixx(i, 1) = 1
            elseif (x(i, 1) .le. 100000.0E0_wp) then
               delta(i, 1) = 9500.0E0_wp
               ifixx(i, 1) = 1
            else
               delta(i, 1) = 144000.0E0_wp
               ifixx(i, 1) = 1
            end if

!  Set weights
            if (x(i, 1) .eq. 100.0E0_wp .or. x(i, 1) .eq. 150.0E0_wp) then
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
         short = .false.
         isodr = .true.

      elseif (itest .eq. 11) then

!  Test detection of incorrect derivatives

         lun = lunrpt
         write (lun, 1000)
         do i = 1, 2
            write (lun, 1001) itest
            write (lun, 1110)
            lun = lunsum
         end do
         setno = 6
         call dodrxd(title, n, m, np, nq, ldx, x, ldy, y, beta)
         work = zero
         delta = zero
         job = 00022
         short = .false.
         isodr = .false.

      elseif (itest .eq. 12) then

!  Test detection of incorrect derivatives

         lun = lunrpt
         write (lun, 1000)
         do i = 1, 2
            write (lun, 1001) itest
            write (lun, 1120)
            lun = lunsum
         end do
         setno = 6
         call dodrxd(title, n, m, np, nq, ldx, x, ldy, y, beta)
         work = zero
         delta = zero
         job = 00020
         short = .false.
         isodr = .true.

      elseif (itest .eq. 13) then

!  Test bounded odr problem where
!  parameters start on bound, move away, hit bound, move away, find minimum.

         lun = lunrpt
         write (lun, 1000)
         do i = 1, 2
            write (lun, 1001) itest
            write (lun, 1010)
            lun = lunsum
         end do
         setno = 9
         call dodrxd(title, n, m, np, nq, ldx, x, ldy, y, beta)
         work = zero
         delta = zero
         job = 00000
         short = .false.
         isodr = .true.
         maxit = 100
         beta(1:2) = (/200.0_wp, 5.0_wp/)
         lower(1:2) = (/0.1_wp, 0.0_wp/)
         upper(1:2) = (/200.0_wp, 5.0_wp/)

      elseif (itest .eq. 14) then

!  Test bounded odr problem where
!  bounds are never hit.

         lun = lunrpt
         write (lun, 1000)
         do i = 1, 2
            write (lun, 1001) itest
            write (lun, 1010)
            lun = lunsum
         end do
         setno = 9
         call dodrxd(title, n, m, np, nq, ldx, x, ldy, y, beta)
         work = zero
         delta = zero
         job = 00000
         short = .false.
         isodr = .true.
         maxit = 100
         lower(1:2) = (/0.0_wp, 0.0_wp/)
         upper(1:2) = (/400.0_wp, 6.0_wp/)

      elseif (itest .eq. 15) then

!  Test bounded odr problem where
!  minimum is on boundary.

         lun = lunrpt
         write (lun, 1000)
         do i = 1, 2
            write (lun, 1001) itest
            write (lun, 1010)
            lun = lunsum
         end do
         setno = 9
         call dodrxd(title, n, m, np, nq, ldx, x, ldy, y, beta)
         work = zero
         delta = zero
         job = 00000
         short = .false.
         isodr = .true.
         maxit = 1000
         beta(1:2) = (/200.0_wp, 3.0_wp/)
         lower(1:2) = (/1.1_wp, 0.0_wp/)
         upper(1:2) = (/400.0_wp, 6.0_wp/)
         tsttol = 500.0_wp

      elseif (itest .eq. 16) then

!  Test bounded odr problem where
!  initial BETA is outside bounds.

         lun = lunrpt
         write (lun, 1000)
         do i = 1, 2
            write (lun, 1001) itest
            write (lun, 1010)
            lun = lunsum
         end do
         setno = 9
         call dodrxd(title, n, m, np, nq, ldx, x, ldy, y, beta)
         work = zero
         delta = zero
         job = 00000
         short = .false.
         isodr = .true.
         maxit = 1000
         beta(1:2) = (/200.0_wp, 7.0_wp/)
         lower(1:2) = (/1.1_wp, 0.0_wp/)
         upper(1:2) = (/200.0_wp, 5.0_wp/)

      elseif (itest .eq. 17) then

!  Test bounded odr problem where bounds are ill defined.

         lun = lunrpt
         write (lun, 1000)
         do i = 1, 2
            write (lun, 1001) itest
            write (lun, 1010)
            lun = lunsum
         end do
         setno = 9
         call dodrxd(title, n, m, np, nq, ldx, x, ldy, y, beta)
         work = zero
         delta = zero
         job = 00000
         short = .false.
         isodr = .true.
         maxit = 1000
         beta(1:2) = (/200.0_wp, 2.0_wp/)
         lower(1:2) = (/10.0_wp, 0.0_wp/)
         upper(1:2) = (/2.0_wp, 5.0_wp/)

      elseif (itest .eq. 18) then

!  Test bounded odr problem using centered differences where
!  parameters start on bound, move away, hit bound, move away, find minimum.

         lun = lunrpt
         write (lun, 1000)
         do i = 1, 2
            write (lun, 1001) itest
            write (lun, 1010)
            lun = lunsum
         end do
         setno = 9
         call dodrxd(title, n, m, np, nq, ldx, x, ldy, y, beta)
         work = zero
         delta = zero
         job = 00010
         short = .false.
         isodr = .true.
         maxit = 100
         beta(1:2) = (/200.0_wp, 5.0_wp/)
         lower(1:2) = (/0.1_wp, 0.0_wp/)
         upper(1:2) = (/200.0_wp, 5.0_wp/)

      elseif (itest .eq. 19) then

!  Test bounded odr problem when bounds are too small.
!  Parameters start on bound.

         lun = lunrpt
         write (lun, 1000)
         do i = 1, 2
            write (lun, 1001) itest
            write (lun, 1010)
            lun = lunsum
         end do
         setno = 9
         call dodrxd(title, n, m, np, nq, ldx, x, ldy, y, beta)
         work = zero
         delta = zero
         job = 00010
         short = .false.
         isodr = .true.
         maxit = 100
         beta(1:2) = (/200.0_wp, 5.0_wp/)
         upper(1) = 200.0_wp
         lower(2) = 5.0_wp
         lower(1) = upper(1) - 400*upper(1)*epsmac &
                    + upper(1)*epsmac

         upper(2) = lower(2) + 400*lower(2)*epsmac &
                    - lower(2)*epsmac

      elseif (itest .eq. 20) then

!  Test bounded odr problem when bounds are just big enough for ndigit
!  calculation but too small for difference calculation.
!  Parameters start on bound.

         lun = lunrpt
         write (lun, 1000)
         do i = 1, 2
            write (lun, 1001) itest
            write (lun, 1010)
            lun = lunsum
         end do
         setno = 9
         call dodrxd(title, n, m, np, nq, ldx, x, ldy, y, beta)
         work = zero
         delta = zero
         job = 00000
         short = .false.
         isodr = .true.
         maxit = 100
         beta(1:2) = (/-200.0_wp, -5.0_wp/)
         upper(1) = -200.0_wp
         lower(2) = -5.0_wp
         lower(1) = upper(1) + 400*upper(1)*epsmac
         upper(2) = lower(2) - 400*lower(2)*epsmac

      elseif (itest .eq. 21) then

!  Test bounded odr problem when bounds are too small for derivative
!  step sizes using forward differences.  Parameters start on bound.

         lun = lunrpt
         write (lun, 1000)
         do i = 1, 2
            write (lun, 1001) itest
            write (lun, 1010)
            lun = lunsum
         end do
         setno = 9
         call dodrxd(title, n, m, np, nq, ldx, x, ldy, y, beta)
         work = zero
         delta = zero
         job = 00000
         short = .false.
         isodr = .true.
         maxit = 100
         beta(1:2) = (/-200.0_wp, -5.0_wp/)
         upper(1) = -200.0_wp
         lower(2) = -5.0_wp
         lower(1) = upper(1) + upper(1)*epsmac
         upper(2) = lower(2) - lower(2)*epsmac

      elseif (itest .eq. 22) then

!  Test bounded odr problem when first parameter is fixed and second is bounded.
!  However, set the bounds on the first parameter to exclude the correct value
!  of the second parameter.  This will exercise the packing and unpacking of
!  parameters and ensure that bounds and fixed parameters can be mixed.

         lun = lunrpt
         write (lun, 1000)
         do i = 1, 2
            write (lun, 1001) itest
            write (lun, 1010)
            lun = lunsum
         end do
         setno = 10
         call dodrxd(title, n, m, np, nq, ldx, x, ldy, y, beta)
         work = zero
         delta = zero
         job = 00010
         short = .false.
         isodr = .true.
         maxit = 100
         beta(1:2) = (/2.5_wp, 1.5_wp/)
         lower(1:2) = (/2.5_wp, 1.1_wp/)
         upper(1:2) = (/10.0_wp, 5.0_wp/)
         ifixb(1:2) = (/0, 1/)

      elseif (itest .eq. 23) then

!  Similar to test 22 but without bounds.

         lun = lunrpt
         write (lun, 1000)
         do i = 1, 2
            write (lun, 1001) itest
            write (lun, 1010)
            lun = lunsum
         end do
         setno = 10
         call dodrxd(title, n, m, np, nq, ldx, x, ldy, y, beta)
         work = zero
         delta = zero
         job = 00010
         short = .false.
         isodr = .true.
         maxit = 100
         beta(1:2) = (/2.5_wp, 1.5_wp/)
         lower(1:2) = -huge(1.0_wp)
         upper(1:2) = huge(1.0_wp)
         ifixb(1:2) = (/0, 1/)

      end if

      call dodrxw &
         (n, m, np, nq, ldwe1, ld2we1, isodr, liwmin, lwmin)

!  Compute solution

      write (lunrpt, 2200) title
      write (lunsum, 2200) title
      if (short) then
         call odr(fcn=dodrxf, &
                  n=n, m=m, np=np, nq=nq, &
                  beta=beta, &
                  y=y, x=x, &
                  delta=delta, &
                  we=we(1:ldwe1, 1:ld2we1, :), &
                  wd=wd(1:ldwd1, 1:ld2wd1, :), &
                  job=job, &
                  iprint=iprint, lunerr=lunerr, lunrpt=lunrpt, &
                  work=work, iwork=iwork, &
                  info=info)
      else
         call odr(fcn=dodrxf, &
                  n=n, m=m, np=np, nq=nq, &
                  beta=beta, &
                  y=y, x=x, &
                  delta=delta, &
                  we=we(1:ldwe1, 1:ld2we1, :), &
                  wd=wd(1:ldwd1, 1:ld2wd1, :), &
                  ifixb=ifixb, ifixx=ifixx(1:ldifx, :), &
                  job=job, ndigit=ndigit, taufac=taufac, &
                  sstol=sstol, partol=partol, maxit=maxit, &
                  iprint=iprint, lunerr=lunerr, lunrpt=lunrpt, &
                  stpb=stpb, stpd=stpd(1:ldstpd, :), &
                  sclb=sclb, scld=scld(1:ldscld, :), &
                  work=work, iwork=iwork, &
                  lower=lower(1:np), upper=upper(1:np), &
                  info=info)
      end if

!  Compare results with those obtained on the cray ymp or the intel xeon running
!  Linux using REAL (wp) version of ODRPACK95

      bnrm = dnrm2(np, beta, 1)
      call dwght(n, m, wd, ldwd1, ld2wd1, reshape(work(1:n*m), (/n, m/)), &
                 tempretl(1:n, 1:m))
      wrk(1:n*m) = reshape(tempretl(1:n, 1:m), (/n*m/))
      wssdel = ddot(n*m, work(1:n*m), 1, wrk(1), 1)
      call dwght(n, nq, we, ldwe1, ld2we1, &
                 reshape(work(n*m + 1:n*m + 1 + n*nq - 1), (/n, nq/)), &
                 tempretl(1:n, 1:nq))
      wrk(n*m + 1:n*m + 1 + n*nq - 1) = reshape(tempretl(1:n, 1:nq), (/n*nq/))
      wsseps = ddot(n*nq, work(n*m + 1:n*m + 1 + n*nq - 1), 1, &
                    wrk(n*m + 1:n*m + 1 + n*nq - 1), 1)
      wss = wsseps + wssdel

      if (sstol .lt. zero) then
         sstol = sqrt(epsmac)
      else
         sstol = min(sstol, one)
      end if

      if (partol .lt. zero) then
         partol = epsmac**(two/three)
      else
         partol = min(partol, one)
      end if

      if (info .ge. 10000) then
         if (idpymp(itest) .eq. info) then
            fails = .false.
            msg = 1
         else
            fails = .true.
            msg = 3
         end if

      elseif (mod(info, 10) .eq. 1) then
         fails = abs(wss - dpymp(2, itest)) .gt. &
                 dpymp(2, itest)*sstol*tsttol
         msg = 2

      elseif (mod(info, 10) .eq. 2) then
         fails = abs(bnrm - dpymp(1, itest)) .gt. &
                 dpymp(1, itest)*partol*tsttol
         msg = 2

      elseif (mod(info, 10) .eq. 3) then
         fails = (abs(wss - dpymp(2, itest)) .gt. &
                  dpymp(2, itest)*sstol*tsttol) &
                 .and. &
                 (abs(bnrm - dpymp(1, itest)) .gt. &
                  dpymp(1, itest)*partol*tsttol)
         msg = 2

      elseif ((mod(info, 10) .eq. 4) .and. &
              (idpymp(itest) .eq. 4)) then
         fails = .false.
         msg = 1

      elseif (info .eq. idpymp(itest)) then
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
         if (ewrt .eq. zero) then
!---------------------------^--------------------------------------------------
!!! FPT - 3087 REAL or COMPLEX quantity tested for exact equality/inequality
!------------------------------------------------------------------------------
            ewrt = one
         end if
         if (ewrt2 .eq. zero) then
!----------------------------^-------------------------------------------------
!!! FPT - 3087 REAL or COMPLEX quantity tested for exact equality/inequality
!------------------------------------------------------------------------------
            ewrt2 = one
         end if
         write (lun, 3220) ' RELATIVE ERROR       = ', &
            abs(dpymp(1, itest) - bnrm)/ewrt, &
            abs(dpymp(2, itest) - wss)/ewrt2

         if (msg .eq. 1) then
            write (lun, 3310)
         elseif (msg .eq. 2) then
            if (fails) then
               write (lun, 3320)
            else
               write (lun, 3330)
            end if
         elseif (msg .eq. 3) then
            write (lun, 3340)
         elseif (msg .eq. 4) then
            write (lun, 3350)
         end if

         lun = lunsum
      end do
   end do test

   write (lunrpt, 1000)
   if (failed) then
      write (lunrpt, 4100)
      write (lunsum, 4100)
      passed = .false.
   else
      write (lunrpt, 4200)
      write (lunsum, 4200)
      passed = .true.
   end if

   if (passed) then
      stop "Solution tests passed. See report."
   else
      error stop "Solution tests failed. See report."
   end if

!  Format statements

1000 format('1')
1001 format(' Example ', I2/)
1010 format(' Test simple odr problem'/ &
          ' with analytic derivatives', &
          ' using ODR.')
1020 format(' Test simple OLS problem'/ &
          ' with finite difference derivatives', &
          ' using ODR.')
1030 format(' Test parameter fixing capabilities', &
          ' for poorly scaled OLS problem'/ &
          ' with analytic derivatives', &
          ' using ODR.')
1040 format(' Test weighting capabilities', &
          ' for ODR problem'/ &
          ' with analytic derivatives', &
          ' using ODR. '/ &
          ' also shows solution of poorly scaled', &
          ' ODR problem.'/ &
          ' (derivative checking turned off.)')
1050 format(' Test DELTA initialization capabilities'/ &
          ' and use of ISTOP to restrict parameter values', &
          ' for ODR problem'/ &
          ' with analytic derivatives', &
          ' using ODR.')
1060 format(' Test stiff stopping conditions', &
          ' for unscaled ODR problem'/ &
          ' with analytic derivatives', &
          ' using ODR.')
1070 format(' Test restart', &
          ' for unscaled ODR problem'/ &
          ' with analytic derivatives', &
          ' using ODR.')
1080 format(' Test use of TAUFAC to restrict first step', &
          ' for ODR problem'/ &
          ' with finite difference derivatives', &
          ' using ODR.')
1090 format(' Test implicit model', &
          ' for OLS problem'/ &
          ' using ODR.')
1100 format(' Test multiresponse model', &
          ' for ODR problem'/ &
          ' with finite difference derivatives', &
          ' using ODR.')
1110 format(' Test detection of questionable analytic derivatives', &
          ' for OLS problem'/ &
          ' using ODR.')
1120 format(' Test detection of incorrect analytic derivatives', &
          ' for ODR problem'/ &
          ' with analytic derivatives', &
          ' using ODR.')
2200 format(' Data Set Reference: ', A80)
3100 format &
      (/' Comparison of new results with', &
        ' REAL (wp) Cray YMP or Intel X86 (Linux) '/ &
        ' Result:'// &
        '                         Norm of BETA', &
        '        Sum of Squared WTD OBS Errors  INFO')
3210 format &
      (/A25/1P, 2E37.30, I6)
3220 format &
      (/A25, 1P, D12.5, 25X, D12.5, I6)
3310 format &
      (/' *** Stopping conditions', &
        ' show convergence not attained. ***'/ &
        '        no further comparisons made between results.'//)
3320 format &
      (//' *** WARNING ***', &
        ' results do not agree to within stopping tolerance. ***'//)
3330 format &
      (//' *** Results agree to within stopping tolerance. ***'//)
3340 format &
      (//' *** WARNING ***', &
        ' stopping conditions do not agree. ***'//)
3350 format &
      (//' *** WARNING ***', &
        ' unexpected stopping condition.', &
        '  please contact package authors. ***'//)
4100 format &
      (/// &
        ' *** Summary:', &
        ' one or more tests do not agree with expected results. ***')
4200 format &
      (/// &
        ' *** Summary:', &
        ' all tests agree with expected results. ***')

end

subroutine dodrxd &
   (title, n, m, np, nq, ldx, x, ldy, y, beta)
!***Begin Prologue  DODRXD
!***Refer to  ODR
!***Routines Called  (NONE)
!***Date Written   860529   (YYMMDD)
!***Revision Date  920619   (YYMMDD)
!***Purpose  Set up data for ODRPACK95 exerciser
!***End Prologue  DODRXD
!
!...Used modules
   use odrpack_kinds, only: wp
   implicit none

!...Parameters
   integer, parameter :: maxn = 50, maxm = 3, maxnp = 10, maxnq = 3, maxset = 16

!...Scalar arguments
   integer :: ldx, ldy, m, n, np, nq
   character(len=80) :: title

!...Array arguments
   real(wp) :: beta(*), x(ldx, *), y(ldy, *)

!...Scalars in common
   integer :: setno

!...Local scalars
   integer :: i, j, k, l

!...Local arrays
   real(wp) :: bdata(maxnp, maxset), xdata(maxn, maxm, maxset), ydata(maxn, maxnq, maxset)
   integer :: mdata(maxset), ndata(maxset), npdata(maxset), nqdata(maxset)
   character(len=80) :: tdata(maxset)

!...Common blocks
   common/setid/setno

!...Data statements
   data &
      tdata(1) &
      /' BOGGS, BYRD AND SCHNABEL, 1985, EXAMPLE 1'/
   data &
      ndata(1), mdata(1), npdata(1), nqdata(1) &
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
      ndata(2), mdata(2), npdata(2), nqdata(2) &
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
      ndata(3), mdata(3), npdata(3), nqdata(3) &
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
      ndata(4), mdata(4), npdata(4), nqdata(4) &
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
      ndata(5), mdata(5), npdata(5), nqdata(5) &
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
      ndata(6), mdata(6), npdata(6), nqdata(6) &
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
      ndata(7), mdata(7), npdata(7), nqdata(7) &
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
      ndata(8), mdata(8), npdata(8), nqdata(8) &
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
      ndata(9), mdata(9), npdata(9), nqdata(9) &
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
      ndata(10), mdata(10), npdata(10), nqdata(10) &
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

!...Variable definitions (alphabetically)
!       BDATA:   The function parameter for each data set.
!       BETA:    The function parameters.
!       I:       An indexing variable.
!       J:       An indexing variable.
!       L:       An indexing variable.
!       LDX:     The leading dimension of array X.
!       M:       The number of columns of data in the explanatory variable.
!       MDATA:   The number of columns of data in the explanatory variable
!       in each data set.
!       N:       The number of observations.
!       NDATA:   The number of observations per data set.
!       NP:      The number of function parameters.
!       NPDATA:  The number of function parameters in each data set.
!       NQDATA:  The number of responses per observation in each data set.
!       SETNO:   The number of the data set being analyzed.
!       TDATA:   The reference for the each of the data sets.
!       TITLE:   The reference for the data set being analyzed.
!       X:       The explanatory variables.
!       XDATA:   The explanatory variables for each data set.
!       Y:       The response variable.
!       YDATA:   The response variables for each data set.
!
!
!***First executable statement  DODRXD

   title = tdata(setno)

   n = ndata(setno)
   m = mdata(setno)
   np = npdata(setno)
   nq = nqdata(setno)

   do l = 1, nq
      do i = 1, n
         y(i, l) = ydata(i, l, setno)
      end do
   end do

   do j = 1, m
      do i = 1, n
         x(i, j) = xdata(i, j, setno)
      end do
   end do

   do k = 1, np
      beta(k) = bdata(k, setno)
   end do

end

subroutine dodrxf &
   (n, m, np, nq, &
    ldn, ldm, ldnp, &
    beta, xplusd, &
    ifixb, ifixx, ldifx, &
    ideval, f, fjacb, fjacd, &
    istop)
!***Begin Prologue  DODRXF
!***Refer to  ODR
!***Routines Called  (NONE)
!***Date Written   860529   (YYMMDD)
!***Revision Date  920619   (YYMMDD)
!***Purpose  Compute jacobian matricies for ODRPACK95 exerciser
!***End Prologue  DODRXF

!...Used modules
   use odrpack_kinds, only: wp, zero, one
   implicit none

!...Scalar arguments
   integer :: ideval, istop, ldifx, ldm, ldn, ldnp, m, n, np, nq

!...Array arguments
   real(wp) :: beta(np), f(ldn, nq), fjacb(ldn, ldnp, nq), fjacd(ldn, ldm, nq), xplusd(ldn, m)
   integer :: ifixb(np), ifixx(ldifx, m)

!...Scalar parameters
   integer, parameter :: maxnp = 10

!...Scalars in common
   integer :: setno

!...Arrays in common
   real(wp) :: lower(maxnp), upper(maxnp)

!...Local scalars
   real(wp) :: ctheta, fac1, fac2, fac3, fac4, freq, omega, phi, pi, r, stheta, theta
   integer :: i, j, k

!...Common blocks
   common/setid/setno
   common/bounds/lower, upper

!...Variable definitions (alphabetically)
!       BETA:    Current values of parameters
!       F:       Predicted function values
!       FAC1:    A factors or terms used in computing the jacobians.
!       FAC2:    A factors or terms used in computing the jacobians.
!       FAC3:    A factors or terms used in computing the jacobians.
!       FAC4:    A factors or terms used in computing the jacobians.
!       FJACB:   Jacobian with respect to BETA
!       FJACD:   Jacobian with respect to errors DELTA
!       IDEVAL:  Indicator for selecting computation to be performed
!       IFIXB:   Indicators for "fixing" parameters (BETA)
!       IFIXX:   Indicators for "fixing" explanatory variable (X)
!       LDIFX:   Leading dimension of array IFIXX
!       ISTOP:   Stopping condition, where
!       0 means current BETA and X+DELTA were
!       acceptable and values were computed successfully
!       1 means current BETA and X+DELTA are
!       not acceptable;  ODRPACK95 should select
!       values closer to most recently used values
!       -1 means current BETA and X+DELTA are
!       not acceptable;  ODRPACK95 should stop
!       LDN:     Leading dimension declarator equal or exceeding N
!       LDM:     Leading dimension declarator equal or exceeding M
!       LDNP:    Leading dimension declarator equal or exceeding NP
!       M:       The number of columns of data in the explanatory variable.
!       N:       The number of observations.
!       NP:      The number of function parameters.
!       NQ:      The number of responses per observation.
!       ONE:     The value 1.0E0_wp.
!       SETNO:   The number of the data set being analyzed.
!       XPLUSD:  Current value of explanatory variable, i.e., X + DELTA
!       zero:    The value 0.0E0_wp.
!
!
!***First executable statement  DODRXF

!  Check for BETA outside bounds.  Return with error if BETA outside bounds.

   if (any(lower(1:np) .gt. beta(1:np))) then
      istop = -1
      return
   end if

   if (any(upper(1:np) .lt. beta(1:np))) then
      istop = -2
      return
   end if

   if (setno .eq. 1) then

!  Setno. 1:  Boggs, Byrd and Schnabel, 1985, example 1

      if (beta(1) .le. 1.01E0_wp) then
         istop = 0

         if (mod(ideval, 10) .ne. 0) then
            do i = 1, n
               f(i, 1) = beta(1)/(xplusd(i, 1) - beta(2))
            end do
         end if

         if (mod(ideval/10, 10) .ne. 0) then
            do i = 1, n
               fjacb(i, 1, 1) = one/(xplusd(i, 1) - beta(2))
               fjacb(i, 2, 1) = beta(1)*(xplusd(i, 1) - beta(2))**(-2)
            end do
         end if

         if (mod(ideval/100, 10) .ne. 0) then
            do i = 1, n
               fjacd(i, 1, 1) = -beta(1)*(xplusd(i, 1) - beta(2))**(-2)
            end do
         end if

      else
         istop = 1
      end if

   elseif (setno .eq. 2) then

!  Setno. 2:  Boggs, Byrd and Schnabel, 1985, example 2

      istop = 0

      do i = 1, n
         fac1 = (beta(2)*xplusd(i, 1) + beta(3)*xplusd(i, 2) - one)

         if (mod(ideval, 10) .ne. 0) then
            f(i, 1) = beta(1)/fac1
         end if

         if (mod(ideval/10, 10) .ne. 0) then
            fjacb(i, 1, 1) = one/fac1
            fjacb(i, 2, 1) = -beta(1)*(fac1**(-2))*xplusd(i, 1)
            fjacb(i, 3, 1) = -beta(1)*(fac1**(-2))*xplusd(i, 2)
         end if

         if (mod(ideval/100, 10) .ne. 0) then
            fjacd(i, 1, 1) = -beta(1)*(fac1**(-2))*beta(2)
            fjacd(i, 2, 1) = -beta(1)*(fac1**(-2))*beta(3)
         end if
      end do

   elseif (setno .eq. 3) then

!  Setno. 3:  Boggs, Byrd and Schnabel, 1985, example 3

      istop = 0

      if (mod(ideval, 10) .ne. 0) then
         do i = 1, n
            f(i, 1) = zero
            do j = 1, 4
               f(i, 1) = f(i, 1) + beta(j)/(xplusd(i, 1) + beta(j + 5))
            end do
            f(i, 1) = f(i, 1) + beta(5)
         end do
      end if

      if (mod(ideval/10, 10) .ne. 0) then
         do i = 1, n
            fjacb(i, 5, 1) = one
            do k = 1, 4
               fjacb(i, k, 1) = one/(xplusd(i, 1) + beta(k + 5))
               fjacb(i, k + 5, 1) = -beta(k)*(xplusd(i, 1) + beta(k + 5))**(-2)
            end do
         end do
      end if

      if (mod(ideval/100, 10) .ne. 0) then
         do i = 1, n
            fjacd(i, 1, 1) = zero
            do k = 4, 1, -1
               fjacd(i, 1, 1) = fjacd(i, 1, 1) - beta(k)*(xplusd(i, 1) + beta(k + 5))**(-2)
            end do
         end do
      end if

   elseif (setno .eq. 4) then

!  Setno. 4:  Himmelblau, 1970, example 6.2-4, page 188

      istop = 0

      if (mod(ideval, 10) .ne. 0) then
         do i = 1, n
            f(i, 1) = beta(1)*xplusd(i, 1) + &
                      beta(2)*exp(beta(3)*xplusd(i, 2))
         end do
      end if

      if (mod(ideval/10, 10) .ne. 0) then
         do i = 1, n
            fjacb(i, 1, 1) = xplusd(i, 1)
            fjacb(i, 2, 1) = exp(beta(3)*xplusd(i, 2))
            fjacb(i, 3, 1) = beta(2)*exp(beta(3)*xplusd(i, 2))*xplusd(i, 2)
         end do
      end if

      if (mod(ideval/100, 10) .ne. 0) then
         do i = 1, n
            fjacd(i, 1, 1) = beta(1)
            fjacd(i, 2, 1) = beta(2)*exp(beta(3)*xplusd(i, 2))*beta(3)
         end do
      end if

   elseif (setno .eq. 5) then

!  Setno. 5:  Draper and Smith, 1981, exercise i, page 521-522

      istop = 0

      if (mod(ideval, 10) .ne. 0) then
         do i = 1, n
            f(i, 1) = exp(-beta(1)*xplusd(i, 1)* &
                          exp(-beta(2)*(one/xplusd(i, 2) - one/620.0E0_wp)))
         end do
      end if

      if (mod(ideval/10, 10) .ne. 0) then
         do i = 1, n
            fac1 = one/xplusd(i, 2) - one/620.0E0_wp
            fac2 = exp(-beta(2)*fac1)
            fac3 = beta(1)*xplusd(i, 1)
            fac4 = exp(-fac3*fac2)

            fjacb(i, 1, 1) = -fac4*xplusd(i, 1)*fac2
            fjacb(i, 2, 1) = fac4*fac3*fac2*fac1

            if (mod(ideval/100, 10) .ne. 0) then
               fjacd(i, 1, 1) = -fac4*beta(1)*fac2
               fjacd(i, 2, 1) = -fac4*fac3*fac2* &
                                beta(2)/xplusd(i, 2)**2
            end if
         end do
      end if

   elseif (setno .eq. 6) then

!  Setno. 6:  Powell and Macdonald, 1972, tables 7 and 8, page 153-154
!          N.B.  this derivative is intentionally coded incorrectly

      istop = 0

      if (mod(ideval, 10) .ne. 0) then
         do i = 1, n
            f(i, 1) = beta(1)* &
                      (one + beta(3)*xplusd(i, 1)/beta(2))**(-one/beta(3))
         end do
      end if

      if (mod(ideval/10, 10) .ne. 0) then
         do i = 1, n
            fjacb(i, 1, 1) = zero
            fjacb(i, 2, 1) = zero
            fjacb(i, 3, 1) = zero

            if (mod(ideval/100, 10) .ne. 0) then
               fjacd(i, 1, 1) = xplusd(i, 1)
            end if
         end do
      end if

   elseif (setno .eq. 7) then

!  Setno. 7:  Fuller, 1987, table 3.2.10, pages 244-245
!          N.B.  this derivative is intentionally coded incorrectly

      istop = 0

      if (mod(ideval, 10) .ne. 0) then
         do i = 1, n
            f(i, 1) = beta(3)*(xplusd(i, 1) - beta(1))**2 + &
                      2*beta(4)*(xplusd(i, 1) - beta(1))* &
                      (xplusd(i, 2) - beta(2)) + &
                      beta(5)*(xplusd(i, 2) - beta(2))**2 - 1.0E0_wp
         end do
      end if

      if (mod(ideval/10, 10) .ne. 0) then
         do i = 1, n
            fjacb(i, 1, 1) = zero
            fjacb(i, 2, 1) = zero
            fjacb(i, 3, 1) = zero
            fjacb(i, 4, 1) = zero
            fjacb(i, 5, 1) = zero

            if (mod(ideval/100, 10) .ne. 0) then
               fjacd(i, 1, 1) = zero
               fjacd(i, 2, 1) = zero
            end if
         end do
      end if

   elseif (setno .eq. 8) then

!  Setno. 8:  Bates and Watts, 1988, table A1.13, pages 280-281
!          N.B.  This derivative is intentionally coded incorrectly

      do i = 1, n
         if (xplusd(i, 1) .lt. 0.0E0_wp) then
            istop = 1
            return
         end if
      end do
      istop = 0

      if (mod(ideval, 10) .ne. 0) then
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

      if (mod(ideval/10, 10) .ne. 0) then
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

            if (mod(ideval/100, 10) .ne. 0) then
               fjacd(i, 1, 1) = zero
               fjacd(i, 1, 2) = zero
            end if
         end do
      end if

   elseif (setno .eq. 9) then

!  Setno. 9:  Zwolak, Watson, and Tyson, 2004.

      istop = 0

      if (mod(ideval, 10) .ne. 0) then
         do i = 1, n
            f(i, 1) = beta(1)*exp(beta(2)*xplusd(i, 1))
         end do
      end if

      if (mod(ideval/10, 10) .ne. 0) then
         do i = 1, n
            fjacb(i, 1, 1) = exp(beta(2)*xplusd(i, 1))
            fjacb(i, 2, 1) = beta(1)*xplusd(i, 1)*exp(beta(2)*xplusd(i, 1))
         end do
      end if

      if (mod(ideval/100, 10) .ne. 0) then
         do i = 1, n
            fjacd(i, 1, 1) = beta(1)*beta(2)*exp(beta(2)*xplusd(i, 1))
         end do
      end if

   elseif (setno .eq. 10) then

!  Setno. 10:  Zwolak, Watson, and Tyson, 2005.

      istop = 0

      if (mod(ideval, 10) .ne. 0) then
         do i = 1, n
            f(i, 1) = beta(1)/2.0_wp*exp(beta(2)*xplusd(i, 1))
         end do
      end if

      if (mod(ideval/10, 10) .ne. 0) then
         do i = 1, n
            fjacb(i, 1, 1) = exp(beta(2)*xplusd(i, 1))/2.0_wp
            fjacb(i, 2, 1) = beta(1)/2.0_wp*xplusd(i, 1)*exp(beta(2)*xplusd(i, 1))
         end do
      end if

      if (mod(ideval/100, 10) .ne. 0) then
         do i = 1, n
            fjacd(i, 1, 1) = beta(1)/2.0_wp*beta(2)*exp(beta(2)*xplusd(i, 1))
         end do
      end if

   end if

end

subroutine dodrxw(maxn, maxm, maxnp, maxnq, ldwe, ld2we, isodr, liwmin, lwmin)
!***Begin Prologue  DODRXW
!***Refer to  ODR
!***Routines Called  (NONE)
!***Date Written   890205   (YYMMDD)
!***Revision Date  920619   (YYMMDD)
!***Purpose  Compute minimum lengths for work vectors
!***Routines Called  NONE
!***End Prologue  DODRXW
!
!...Used modules
   implicit none

!...Scalar arguments
   integer :: ldwe, ld2we, liwmin, lwmin, maxn, maxm, maxnp, maxnq
   logical :: isodr

!...Variable definitions (alphabetically)
!       ISODR:   The variable designating whether the solution is by odr
!       (ISODR=TRUE) or by ols (ISODR=FALSE).
!       LDWE:    The leading dimension of array WE.
!       LD2WE:   The second dimension of array WE.
!       LIWMIN:  The minimum length of vector IWORK for a given problem.
!       LWMIN:   The minimum length of vector WORK for a given problem.
!       MAXM:    The number of columns in the explanatory variable.
!       MAXN:    The number of observations.
!       MAXNP:   The number of function parameters.
!       MAXNQ:   The number of responses per observation.

!***First executable statement  DODRXW

   liwmin = 20 + maxnp + maxnq*(maxnp + maxm)
   if (isodr) then
      lwmin = 18 + 11*maxnp + maxnp**2 + maxm + maxm**2 + &
              4*maxn*maxnq + 6*maxn*maxm + 2*maxn*maxnq*maxnp + &
              2*maxn*maxnq*maxm + maxnq**2 + &
              5*maxnq + maxnq*(maxnp + maxm) + ldwe*ld2we*maxnq
   else
      lwmin = 18 + 11*maxnp + maxnp**2 + maxm + maxm**2 + &
              4*maxn*maxnq + 2*maxn*maxm + 2*maxn*maxnq*maxnp + &
              5*maxnq + maxnq*(maxnp + maxm) + ldwe*ld2we*maxnq
   end if

end
