*DTEST
      PROGRAM DTEST
C***Begin Prologue  DTEST
C***Refer to ODR
C***Routines Called  DODRX
C***Date Written   861229   (YYMMDD)
C***Revision Date  920619   (YYMMDD)
C***Purpose  EXERCISE FEATURES OF ODRPACK95 SOFTWARE
C***End Prologue  DTEST

C...Used modules
      USE REAL_PRECISION

C...Scalars in common
      INTEGER
     &   NTEST

C...Local scalars
      REAL (KIND=R8)
     &   TSTFAC
      INTEGER
     &   LUNERR,LUNRPT,LUNSUM
      LOGICAL
     &   PASSED

C...External subroutines
      EXTERNAL
     &   DODRX

C...Common blocks
      COMMON /TSTSET/ NTEST

C***Variable declarations (alphabetically)

C   LUNERR:  The logical unit number used for error messages.
C   LUNRPT:  The logical unit number used for computation reports.
C   LUNSUM:  The logical unit number used for a summary report listing
C            only the test comparisons and not the odrpack generated
C            reports.
C   NTEST:   The number of tests to be run.
C   PASSED:  The variable designating whether the results of all of the 
C            tests agree with those from the cray ymp using double 
C            precision (PASSED=TRUE), or whether some of the results
C            disagreed (PASSED=FALSE).
C   TSTFAC:  The user-supplied factor for scaling the test tolerances 
C            used to check for agreement between computed results and 
C            results obtained using REAL (KIND=R8) version on cray
C            YMP.  Values of TSTFAC greater than one increase the 
C            test tolerances, making the tests easier to pass and 
C            allowing small discrepancies between the computed and 
C            expected results to be automatically discounted.


C***First executable statement  TEST


C  Set up necessary files

C  NOTE:  ODRPACK95 generates computation and error reports on
C         logical unit 6 by default;
C         logical unit 'LUNSUM' used to summarize results of comparisons
C         from exercise routine DODRX.

      LUNRPT = 18
      LUNERR = 18
      LUNSUM = 19

      OPEN(UNIT=LUNRPT,FILE='REPORT')
      OPEN(UNIT=LUNERR,FILE='REPORT')
      OPEN(UNIT=LUNSUM,FILE='SUMMARY')

C  Exercise REAL (KIND=R8) version of ODRPACK95
C  (test reports generated on file 'RESULTS' and
C   summarized in file 'SUMMARY')

      NTEST = 23
      TSTFAC = 1.0E0_R8
      CALL DODRX(TSTFAC,PASSED,LUNSUM)

      END
*DODRX
      SUBROUTINE DODRX
     &   (TSTFAC,PASSED,LUNSUM)
C***Begin Prologue  DODRX
C***Refer to ODR
C***Routines Called  DDOT,DNRM2,ODR,DODRXD,
C                    DODRXF,DODRXW,DWGHT,DZERO
C***Date Written   860529   (YYMMDD)
C***Revision Date  920619   (YYMMDD)
C***Purpose  Exercise features of ODRPACK95 software
C***End Prologue  DODRX

C...Used modules
      USE ODRPACK95
      USE REAL_PRECISION

C...Parameters
      INTEGER
     &   LDWD,LDWE,LD2WD,LD2WE,LIWORK,LWORK,MAXN,MAXM,MAXNP,MAXNQ,NTESTS
      REAL (KIND=R8)
     &   BASE
      PARAMETER
     &   (MAXN=50, MAXM=3, MAXNP=10, MAXNQ=2, NTESTS=23,
     &   LDWE=MAXN, LD2WE=MAXNQ, LDWD=MAXN, LD2WD=MAXM,
     &   LWORK = 18 + 11*MAXNP + MAXNP**2 + MAXM + MAXM**2 + 
     &           4*MAXN*MAXNQ + 6*MAXN*MAXM + 2*MAXN*MAXNQ*MAXNP + 
     &           2*MAXN*MAXNQ*MAXM + MAXNQ**2 +
     &           5*MAXNQ + MAXNQ*(MAXNP+MAXM) + LDWE*LD2WE*MAXNQ,
     &   LIWORK = 20+MAXNP+MAXNQ*(MAXNP+MAXM),
     &   BASE = RADIX(1.0E0_R8))

C...Scalar arguments
      REAL (KIND=R8)
     &   TSTFAC
      INTEGER
     &   LUNSUM
      LOGICAL
     &   PASSED

C...Scalars in common
      INTEGER
     &   NTEST,SETNO

C...Local scalars
      INTEGER   
     &   I,INFO,IPRINT,ITEST,JOB,L,LDIFX,LDSCLD,LDSTPD,LDWD1,LDWE1,
     &   LDX,LDY,LD2WD1,LD2WE1,LIWMIN,LUN,LUNERR,LUNRPT,LWMIN,
     &   M,MAXIT,MSG,N,NDIGIT,NP,NQ
      REAL (KIND=R8)
     &   BNRM,EPSMAC,EWRT,EWRT2,HUNDRD,ONE,P01,P2,PARTOL,SSTOL,
     &   TAUFAC,THREE,TSTTOL,TWO,WSS,WSSDEL,WSSEPS,ZERO
      LOGICAL
     &   FAILED,FAILS,ISODR,SHORT
      CHARACTER TITLE*80

C...Arrays in common
      REAL (KIND=R8)
     &   LOWER(MAXNP),UPPER(MAXNP)

C...Local arrays
      REAL (KIND=R8)
     &   BETA(MAXNP),DELTA(:,:),DPYMP(2,NTESTS),
     &   SCLB(MAXNP),SCLD(MAXN,MAXM),
     &   STPB(MAXNP),STPD(MAXN,MAXM),
     &   WE(MAXN,MAXNQ,MAXNQ),WD(MAXN,MAXM,MAXM),WORK(:),
     &   WRK(MAXN*MAXM+MAXN*MAXNQ),X(MAXN,MAXM),Y(MAXN,MAXNQ),
     &   TEMPRETL(MAXN,MAXM)
      INTEGER
     &   IDPYMP(NTESTS),IFIXB(MAXNP),IFIXX(MAXN,MAXM),IWORK(:)

C...Pointers
      POINTER
     &   DELTA,IWORK,WORK

C...External functions
      REAL (KIND=R8)
     &   DDOT,DNRM2
      EXTERNAL
     &   DDOT,DNRM2

C...External subroutines
      EXTERNAL
     &   DODRXD,DODRXF,DODRXW,DZERO

C...Intrinsic functions
      INTRINSIC
     &   ABS,MOD

C...Common blocks
      COMMON /SETID/SETNO
      COMMON /TSTSET/ NTEST
      COMMON /BOUNDS/ LOWER,UPPER

C...Data statements
      DATA
     &   ZERO,P01,P2,ONE,TWO,THREE,HUNDRD
     &   /0.0E0_R8,0.01E0_R8,0.2E0_R8,1.0E0_R8,2.0E0_R8,3.0E0_R8,
     &   100.0E0_R8/

      DATA
     &   (DPYMP(I,1),I=1,2)
     &   /2.762733195780256808978449342964E+04_R8,
     &    7.532639569022918943695104672512E-04_R8/
      DATA
     &   (DPYMP(I,2),I=1,2)
     &   /2.762732630143673024399942947263E+04_R8,
     &    7.538467722687131506874279314940E-04_R8/
      DATA
     &   (DPYMP(I,3),I=1,2)
     &   /1.069944100000000027940905194068E+09_R8,
     &    1.212808593256056359629660672046E-05_R8/
      DATA
     &   (DPYMP(I,4),I=1,2)
     &   /1.069944100000000026623461142867E+09_R8,
     &    5.452084633790606017572015067556E-07_R8/
      DATA
     &   (DPYMP(I,5),I=1,2)
     &   /1.426988156377258617521571734503E+00_R8,
     &    1.084728687127432219753903919409E+00_R8/
      DATA
     &   (DPYMP(I,6),I=1,2)
     &   /4.261321829513978871872508874025E+00_R8,
     &    1.477967210398420733565424329280E-02_R8/
      DATA
     &   (DPYMP(I,7),I=1,2)
     &   /4.261272307142888464011486769858E+00_R8,
     &    1.477966125465374336804138554559E-02_R8/
      DATA
     &   (DPYMP(I,8),I=1,2)
     &   /4.371487317909745009110272283622E+01_R8,
     &    1.144419474408286067112233592550E-03_R8/
      DATA
     &   (DPYMP(I,9),I=1,2)
     &   /3.099048849376848610380977303924E+00_R8,
     &    8.824708863783850023783338218501E-02_R8/
      DATA
     &   (DPYMP(I,10),I=1,2)
     &   /9.469917836739932584221023234527E+00_R8,
     &    4.205389215588104651198536809880E-01_R8/
      DATA
     &   (DPYMP(I,11),I=1,2)
     &   /3.950949253027682207109233363651E+01_R8,
     &    6.651838750834910819636881506915E+01_R8/
      DATA
     &   (DPYMP(I,12),I=1,2)
     &   /3.950949253027682207109233363651E+01_R8,
     &    6.651838750834910819636881506915E+01_R8/
      DATA
     &   (DPYMP(I,13),I=1,2)
     &   /1.414213562373095000000000000000E+00_R8,
     &    5.250825926608277346013642256883E-26_R8/
      DATA
     &   (DPYMP(I,14),I=1,2)
     &   /1.414213562373095000000000000000E+00_R8,
     &    8.159081600696301507018019048968E-26_R8/
      DATA
     &   (DPYMP(I,15),I=1,2)
     &   /1.486588477064952451556223422813E+00_R8,
     &    1.841690442255357083922717720270E+03_R8/
      DATA
     &   (DPYMP(I,16),I=1,2)
     &   /2.001224625073357401561224833131E+02_R8,
     &    0.000000000000000000000000000000E+00_R8/
      DATA
     &   (DPYMP(I,17),I=1,2)
     &   /2.000099997500125000000000000000E+02_R8,
     &    0.000000000000000000000000000000E+00_R8/
      DATA
     &   (DPYMP(I,18),I=1,2)
     &   /1.414213562373095000000000000000E+00_R8,
     &    5.816277809383742531415846947805E-26_R8/
      DATA
     &   (DPYMP(I,19),I=1,2)
     &   /2.000624902374255782433465356007E+02_R8,
     &    4.568236947482152283374593507328E+30_R8/

      DATA
     &   (DPYMP(I,20),I=1,2)
     &   /2.000624902374255782433465356007E+02_R8,
     &    1.848525209410256939008831977844E+05_R8/

      DATA
     &   (DPYMP(I,21),I=1,2)
     &   /2.000624902374255782433465356007E+02_R8,
     &    1.848525209410256939008831977844E+05_R8/

      DATA
     &   (DPYMP(I,22),I=1,2)
     &   /2.731300056749532689792659000000E+00_R8,
     &    3.378975642596100806258619000000E+05_R8/

      DATA
     &   (DPYMP(I,23),I=1,2)
     &   /2.675757304209387399396291584708E+00_R8,
     &    5.174484505019630309341494012187E-02_R8/

      DATA
     &   (IDPYMP(I),I=1,23)
     &   /1,1,3,1,1,4,1,1,2,1,1023,40100,2,2,3,90100,91000,2,90010,
     &    90020,90010,21,1/

C...Interface blocks
      INTERFACE
      SUBROUTINE DWGHT
     &   (N,M,WT,LDWT,LD2WT,T,WTT)
      USE REAL_PRECISION
      INTEGER
     &   LDWT,LD2WT,M,N
      REAL (KIND=R8)
     &   T(:,:),WT(:,:,:),WTT(:,:)
      END SUBROUTINE
      END INTERFACE

C...Routine names used as subprogram arguments
C   DODRXF:  The user-supplied routine for evaluating the model.

C...Variable definitions (alphabetically)
C   BASE:    The base of floating point numbers on the current machine
C   BETA:    The function parameters.
C   BNRM:    The norm of BETA.
C   DELTA:   The error in the X data.
C   DPYMP:   The floating point results from a cray YMP using
C            REAL (KIND=R8).
C   EPSMAC:  The value of machine precision.
C   EWRT:    A temporary variable for the denominator of the relative error 
C            calculations (error with respect to).
C   EWRT2:   A temporary variable for the denominator of the relative error 
C            calculations (error with respect to).
C   FAILED:  The variable designating whether the results of all of the
C            demonstration runs agreed with those from the cray YMP
C            using REAL (KIND=R8) (FAILED=FALSE) or whether some of
C            the tests disagreed (FAILED=TRUE).
C   FAILS:   The variable designating whether the results of an 
C            individual demonstration run agreed with those from the
C            cray YMP using REAL (KIND=R8) (FAILS=FALSE) or 
C            disagree (FAILS=TRUE).
C   HUNDRD:  The value 100.0E0_R8.
C   I:       An index variable.
C   IDPYMP:  The integer results from a cray YMP using
C            REAL (KIND=R8).
C   IFIXB:   The values designating whether the elements of BETA are 
C            fixed at their input values or not.
C   IFIXX:   The values designating whether the elements of DELTA are 
C            fixed at their input values or not.
C   INFO:    The variable designating why the computations stopped.
C   IPRINT:  The print control variable.
C   ISODR:   The variable designating whether the solution is by odr 
C            (ISODR=TRUE) or by ols (ISODR=FALSE).
C   ITEST:   The number of the current test being run.
C   IWORK:   The integer work space.
C   J:       An index variable.
C   JOB:     The variable controlling problem initialization and 
C            computational method.
C   LDIFX:   The leading dimension of array IFIXX.
C   LDSCLD:  The leading dimension of array SCLD.
C   LDWD:    The leading dimension of array WD.
C   LDWD1:   The leading dimension of array WD as passed to ODRPACK95.
C   LDWE:    The leading dimension of array WE.
C   LDWE1:   The leading dimension of array WE as passed to ODRPACK95.
C   LDX:     The leading dimension of array X.
C   LDY:     The leading dimension of array Y.
C   LD2WD:   The second dimension of array WD.
C   LD2WD1:  The second dimension of array WD as passed to ODRPACK95.
C   LD2WE:   The second dimension of array WE.
C   LD2WE1:  The second dimension of array WE as passed to ODRPACK95.
C   LIWKMN:  The minimum acceptable length of array IWORK.
C   LIWMIN:  The minimum length of vector IWORK for a given problem.
C   LIWORK:  The length of vector IWORK.
C   LUN:     The logical unit number currently being used.
C   LUNERR:  The logical unit number used for error messages.
C   LUNRPT:  The logical unit number used for computation reports.
C   LUNSUM:  The logical unit number used for a summary report.
C   LWKMN:   The minimum acceptable length of array WORK.
C   LWMIN:   The minimum length of vector WORK for a given problem.
C   LWORK:   The length of vector WORK.
C   M:       The number of columns of data in the explanatory variable.
C   MAXIT:   The maximum number of iterations allowed.
C   MSG:     The variable designating which message is to be printed as
C            a result of the comparison with the cray YMP or x86 (Linux) 
C            results.
C   N:       The number of observations.
C   NDIGIT:  The number of accurate digits in the function results, as
C            supplied by the user.
C   NP:      The number of function parameters.
C   NTEST:   The number of tests to be run.
C   NTESTS:  The number of different tests available.
C   ONE:     The value 1.0E0_R8.
C   PASSED:  The variable designating whether the results of all of the
C            demonstration runs agreed with those from the cray YMP
C            using REAL (KIND=R8) (PASSED=TRUE), or whether some of
C            the results disagreed (PASSED=FALSE).
C   P01:     The value 0.01E0_R8.
C   P2:      The value 0.2E0_R8.
C   PARTOL:  The parameter convergence stopping criteria.
C   SCLB:    The scaling values for BETA.
C   SCLD:    The scaling values for DELTA.
C   SETNO:   The number of the data set being analyzed.
C   SHORT:   The variable designating whether ODRPACK95 is invoked by the
C            short-call (SHORT=.TRUE.) or the long-call (SHORT=.FALSE.).
C   SSTOL:   The sum-of-squares convergence stopping tolerance.
C   TAUFAC:  The factor used to compute the initial trust region 
C            diameter.
C   THREE:   The value 3.0E0_R8.
C   TITLE:   The reference for the data set being analyzed.
C   TSTFAC:  The user-supplied factor for scaling the test tolerances
C            used to check for agreement between computed results and
C            results obtained using REAL (KIND=R8) version on cray 
C            YMP.
C   TSTTOL:  The test tolerance used in checking computed values for
C            purposes of determining proper installation.
C   TWO:     The value 2.0E0_R8.
C   WD:      The DELTA weights.
C   WE:      The EPSILON weights.
C   WORK:    The REAL (KIND=R8) work space.
C   WRK:     The REAL (KIND=R8) work space for computing test results.
C   WSS:     The sum of the squared weighted errors.
C   WSSDEL:  The sum of the squared weighted errors in X.
C   WSSEPS:  The sum of the squared weighted errors in Y.
C   X:       The explanatory variable.
C   Y:       The response variable.
C   ZERO:    The value 0.0E0_R8.


C***First executable statement  DODRX


C  Allocate work arrays and DELTA

      ALLOCATE(DELTA(MAXN,MAXM),IWORK(LIWORK),WORK(LWORK))

C  Set logical units for error and computation reports

      LUNERR = 18
      LUNRPT = 18

C  Initialize test tolerance

      IF (TSTFAC.GT.ONE) THEN
         TSTTOL = TSTFAC
      ELSE
         TSTTOL = ONE
      END IF

C  Initialize machine precision

      EPSMAC = BASE**(1-DIGITS(BASE))

C  Initialize leading dimension of X

      LDX = MAXN
      LDY = MAXN

C  Initialize miscellaneous variables used in the exercise procedure

      FAILED = .FALSE.
      SHORT = .TRUE.
      ISODR = .TRUE.
      N = 0

C  Begin exercising ODRPACK95

      DO 400 ITEST=1,NTEST

C  Set control values to invoke default values

         WE(1,1,1) = -ONE
         LDWE1 = LDWE
         LD2WE1 = LD2WE
         WD(1,1,1) = -ONE
         LDWD1 = LDWD
         LD2WD1 = LD2WD

         IFIXB(1) = -1
         IFIXX(1,1) = -1
         LDIFX = MAXN

         NDIGIT = -1
         TAUFAC = -ONE

         SSTOL = -ONE
         PARTOL = -ONE
         MAXIT = -1

         IPRINT = 2112
C        IPRINT = 6616

         STPB(1) = -ONE
         STPD(1,1) = -ONE
         LDSTPD = 1

         SCLB(1) = -ONE
         SCLD(1,1) = -ONE
         LDSCLD = 1

         UPPER(:) = HUGE(ONE)
         LOWER(:) = -HUGE(ONE)


         IF (ITEST.EQ.1) THEN

C  Test simple odr problem 
C  with analytic derivatives.

            LUN = LUNRPT
            WRITE (LUN,1000)
            DO 10 I=1,2
               WRITE (LUN,1001) ITEST
               WRITE (LUN,1010)
               LUN = LUNSUM
   10       CONTINUE
            SETNO = 5
            CALL DODRXD(TITLE,N,M,NP,NQ,LDX,X,LDY,Y,BETA)
            CALL DZERO(LWORK,1,WORK,LWORK)
            DELTA(:,:) = ZERO
            JOB = 00020
            SHORT = .TRUE.
            ISODR = .TRUE.

         ELSE IF (ITEST.EQ.2) THEN

C  Test simple ols problem 
C  with forward difference derivatives.

            LUN = LUNRPT
            WRITE (LUN,1000)
            DO 20 I=1,2
               WRITE (LUN,1001) ITEST
               WRITE (LUN,1020)
               LUN = LUNSUM
   20       CONTINUE
            SETNO = 5
            CALL DODRXD(TITLE,N,M,NP,NQ,LDX,X,LDY,Y,BETA)
            CALL DZERO(LWORK,1,WORK,LWORK)
            DELTA(:,:) = ZERO
            JOB = 00002
            SHORT = .TRUE.
            ISODR = .FALSE.

         ELSE IF (ITEST.EQ.3) THEN

C  Test parameter fixing capabilities for poorly scaled ols problem
C  with analytic derivatives.
C  (derivative checking turned off.)

            LUN = LUNRPT
            WRITE (LUN,1000)
            DO 30 I=1,2
               WRITE (LUN,1001) ITEST
               WRITE (LUN,1030)
               LUN = LUNSUM
   30       CONTINUE
            SETNO = 3
            CALL DODRXD(TITLE,N,M,NP,NQ,LDX,X,LDY,Y,BETA)
            CALL DZERO(LWORK,1,WORK,LWORK)
            DELTA(:,:) = ZERO
            IFIXB(1) = 1
            IFIXB(2) = 1
            IFIXB(3) = 1
            IFIXB(4) = 0
            IFIXB(5) = 1
            IFIXB(6) = 0
            IFIXB(7) = 0
            IFIXB(8) = 0
            IFIXB(9) = 0
            JOB = 00042
            SHORT = .FALSE.
            ISODR = .FALSE.

         ELSE IF (ITEST.EQ.4) THEN

C  Test weighting capabilities for odr problem with
C  analytic derivatives.
C  Also shows solution of poorly scaled odr problem.
C  (derivative checking turned off.)
C  N.B., this run continues from where test 3 left off.

            LUN = LUNRPT
            WRITE (LUN,1000)
            DO 40 I=1,2
               WRITE (LUN,1001) ITEST
               WRITE (LUN,1040)
               LUN = LUNSUM
   40       CONTINUE
            SETNO = 3
            CALL DZERO(LWORK,1,WORK,LWORK)
            DELTA(:,:) = ZERO
            LDWD1 = LDWD
            LDWE1 = LDWE
            LD2WD1 = LD2WD
            LD2WE1 = LD2WE
            DO 45 I=1,N
               WD(I,1,1)   = (P01/ABS(X(I,1)))**2
               WE(I,1,1) = ONE
   45       CONTINUE
            WE(28,1,1) = ZERO
            IFIXB(1) = 1
            IFIXB(2) = 1
            IFIXB(3) = 1
            IFIXB(4) = 0
            IFIXB(5) = 1
            IFIXB(6) = 1
            IFIXB(7) = 1
            IFIXB(8) = 0
            IFIXB(9) = 0
            JOB = 00030
            IPRINT = 2232
            SHORT = .FALSE.
            ISODR = .TRUE.

         ELSE IF (ITEST.EQ.5) THEN

C  Test DELTA initialization capabilities and user-supplied scaling
C  and use of istop to restrict parameter values
C  for odr problem with analytic derivatives.

            LUN = LUNRPT
            WRITE (LUN,1000)
            DO 50 I=1,2
               WRITE (LUN,1001) ITEST
               WRITE (LUN,1050)
               LUN = LUNSUM
   50       CONTINUE
            SETNO = 1
            CALL DODRXD(TITLE,N,M,NP,NQ,LDX,X,LDY,Y,BETA)
            CALL DZERO(LWORK,1,WORK,LWORK)
            DELTA(:,:) = ZERO
            JOB = 01020
            LDSCLD = 1
            SCLD(1,1) = TWO
            SCLB(1) = P2
            SCLB(2) = ONE
            LDWE1 = 1
            LD2WE1 = 1
            WE(1,1,1) = -ONE
            LDWD1  = 1
            LD2WD1 = 1
            WD(1,1,1) = -ONE
            DO 55 I=20,21
               DELTA(I,1) = BETA(1)/Y(I,1) + BETA(2) - X(I,1)
   55       CONTINUE
            SHORT = .FALSE.
            ISODR = .TRUE.

         ELSE IF (ITEST.EQ.6) THEN

C  Test stiff stopping conditions for unscaled odr problem
C  with analytic derivatives.

            LUN = LUNRPT
            WRITE (LUN,1000)
            DO 60 I=1,2
               WRITE (LUN,1001) ITEST
               WRITE (LUN,1060)
               LUN = LUNSUM
   60       CONTINUE
            SETNO = 4
            CALL DODRXD(TITLE,N,M,NP,NQ,LDX,X,LDY,Y,BETA)
            CALL DZERO(LWORK,1,WORK,LWORK)
            DELTA(:,:) = ZERO
            JOB = 00020
            SSTOL = HUNDRD*EPSMAC
            PARTOL = EPSMAC
            MAXIT = 2
            SHORT = .FALSE.
            ISODR = .TRUE.

         ELSE IF (ITEST.EQ.7) THEN

C  Test restart for unscaled odr problem
C  with analytic derivatives.

            LUN = LUNRPT
            WRITE (LUN,1000)
            DO 70 I=1,2
               WRITE (LUN,1001) ITEST
               WRITE (LUN,1070)
               LUN = LUNSUM
   70       CONTINUE
            SETNO = 4
            JOB = 20220
            SSTOL = HUNDRD*EPSMAC
            PARTOL = EPSMAC
            MAXIT = 50
            SHORT = .FALSE.
            ISODR = .TRUE.

         ELSE IF (ITEST.EQ.8) THEN

C  Test use of TAUFAC to restrict first step
C  for odr problem with central difference derivatives.

            LUN = LUNRPT
            WRITE (LUN,1000)
            DO 80 I=1,2
               WRITE (LUN,1001) ITEST
               WRITE (LUN,1080)
               LUN = LUNSUM
   80       CONTINUE
            SETNO = 6
            CALL DODRXD(TITLE,N,M,NP,NQ,LDX,X,LDY,Y,BETA)
            CALL DZERO(LWORK,1,WORK,LWORK)
            DELTA(:,:) = ZERO
            JOB = 00210
            TAUFAC = P01
            SHORT = .FALSE.
            ISODR = .TRUE.

         ELSE IF (ITEST.EQ.9) THEN

C  Test implicit odr problem
C  with forward finite difference derivatives
C  and covariance matrix constructed with recomputed derivatives.


            LUN = LUNRPT
            WRITE (LUN,1000)
            DO 90 I=1,2
               WRITE (LUN,1001) ITEST
               WRITE (LUN,1090)
               LUN = LUNSUM
   90       CONTINUE
            SETNO = 7
            CALL DODRXD(TITLE,N,M,NP,NQ,LDX,X,LDY,Y,BETA)
            CALL DZERO(LWORK,1,WORK,LWORK)
            DELTA(:,:) = ZERO
            JOB = 00001
            PARTOL = EPSMAC**(ONE/THREE)
            SHORT = .TRUE.
            ISODR = .TRUE.

         ELSE IF (ITEST.EQ.10) THEN

C  Test multiresponse odr problem 
C  with central difference derivatives ,
C  DELTA initialized to nonzero values, 
C  variable fixing,  and weighting.

            LUN = LUNRPT
            WRITE (LUN,1000)
            DO 100 I=1,2
               WRITE (LUN,1001) ITEST
               WRITE (LUN,1100)
               LUN = LUNSUM
  100       CONTINUE
            SETNO = 8
            CALL DODRXD(TITLE,N,M,NP,NQ,LDX,X,LDY,Y,BETA)
            CALL DZERO(LWORK,1,WORK,LWORK)
            DELTA(:,:) = ZERO

            LDWD1 = LDWD
            LDWE1 = LDWE
            LD2WD1 = LD2WD
            LD2WE1 = LD2WE
            DO 105 I=1,N
C  Initialize DELTA, and specify first decade of frequencies as fixed
               IF (X(I,1).LT.100.0E0_R8) THEN
                  DELTA(I,1) = 0.0E0_R8
                  IFIXX(I,1) = 0
               ELSE IF (X(I,1).LE.150.0E0_R8) THEN
                  DELTA(I,1) = 0.0E0_R8
                  IFIXX(I,1) = 1
               ELSE IF (X(I,1).LE.1000.0E0_R8) THEN
                  DELTA(I,1) = 25.0E0_R8
                  IFIXX(I,1) = 1
               ELSE IF (X(I,1).LE.10000.0E0_R8) THEN
                  DELTA(I,1) = 560.0E0_R8
                  IFIXX(I,1) = 1
               ELSE IF (X(I,1).LE.100000.0E0_R8) THEN
                  DELTA(I,1) = 9500.0E0_R8
                  IFIXX(I,1) = 1
               ELSE
                  DELTA(I,1) = 144000.0E0_R8
                  IFIXX(I,1) = 1
               END IF

C  Set weights
               IF (X(I,1).EQ.100.0E0_R8 .OR. X(I,1).EQ.150.0E0_R8) THEN
                  WE(I,1,1) = 0.0E0_R8
                  WE(I,1,2) = 0.0E0_R8
                  WE(I,2,1) = 0.0E0_R8
                  WE(I,2,2) = 0.0E0_R8
               ELSE
                  WE(I,1,1) =   559.6E0_R8
                  WE(I,1,2) = -1634.0E0_R8
                  WE(I,2,1) = -1634.0E0_R8
                  WE(I,2,2) =  8397.0E0_R8
               END IF
               WD(I,1,1)    =  (1.0E-4_R8)/(X(I,1)**2)
  105       CONTINUE
            JOB = 00210
            SHORT = .FALSE.
            ISODR = .TRUE.
 
         ELSE IF (ITEST.EQ.11) THEN

C  Test detection of incorrect derivatives

            LUN = LUNRPT
            WRITE (LUN,1000)
            DO 110 I=1,2
               WRITE (LUN,1001) ITEST
               WRITE (LUN,1110)
               LUN = LUNSUM
  110       CONTINUE
            SETNO = 6
            CALL DODRXD(TITLE,N,M,NP,NQ,LDX,X,LDY,Y,BETA)
            CALL DZERO(LWORK,1,WORK,LWORK)
            DELTA(:,:) = ZERO
            JOB = 00022
            SHORT = .FALSE.
            ISODR = .FALSE.

         ELSE IF (ITEST.EQ.12) THEN

C  Test detection of incorrect derivatives

            LUN = LUNRPT
            WRITE (LUN,1000)
            DO 120 I=1,2
               WRITE (LUN,1001) ITEST
               WRITE (LUN,1120)
               LUN = LUNSUM
  120       CONTINUE
            SETNO = 6
            CALL DODRXD(TITLE,N,M,NP,NQ,LDX,X,LDY,Y,BETA)
            CALL DZERO(LWORK,1,WORK,LWORK)
            DELTA(:,:) = ZERO
            JOB = 00020
            SHORT = .FALSE.
            ISODR = .TRUE.

         ELSE IF (ITEST.EQ.13) THEN

C  Test bounded odr problem where
C  parameters start on bound, move away, hit bound, move away, find minimum.

            LUN = LUNRPT
            WRITE (LUN,1000)
            DO I=1,2
               WRITE (LUN,1001) ITEST
               WRITE (LUN,1010)
               LUN = LUNSUM
            END DO
            SETNO = 9
            CALL DODRXD(TITLE,N,M,NP,NQ,LDX,X,LDY,Y,BETA)
            CALL DZERO(LWORK,1,WORK,LWORK)
            DELTA(:,:) = ZERO
            JOB = 00000
            SHORT = .FALSE.
            ISODR = .TRUE.
            MAXIT = 100
            BETA(1:2)  = (/ 200.0_R8, 5.0_R8 /)
            LOWER(1:2) = (/ 0.1_R8, 0.0_R8 /)
            UPPER(1:2) = (/ 200.0_R8, 5.0_R8 /)

         ELSE IF (ITEST.EQ.14) THEN

C  Test bounded odr problem where
C  bounds are never hit.

            LUN = LUNRPT
            WRITE (LUN,1000)
            DO I=1,2
               WRITE (LUN,1001) ITEST
               WRITE (LUN,1010)
               LUN = LUNSUM
            END DO
            SETNO = 9
            CALL DODRXD(TITLE,N,M,NP,NQ,LDX,X,LDY,Y,BETA)
            CALL DZERO(LWORK,1,WORK,LWORK)
            DELTA(:,:) = ZERO
            JOB = 00000
            SHORT = .FALSE.
            ISODR = .TRUE.
            MAXIT = 100
            LOWER(1:2) = (/ 0.0_R8, 0.0_R8 /)
            UPPER(1:2) = (/ 400.0_R8, 6.0_R8 /)

         ELSE IF (ITEST.EQ.15) THEN

C  Test bounded odr problem where
C  minimum is on boundary.

            LUN = LUNRPT
            WRITE (LUN,1000)
            DO I=1,2
               WRITE (LUN,1001) ITEST
               WRITE (LUN,1010)
               LUN = LUNSUM
            END DO
            SETNO = 9
            CALL DODRXD(TITLE,N,M,NP,NQ,LDX,X,LDY,Y,BETA)
            CALL DZERO(LWORK,1,WORK,LWORK)
            DELTA(:,:) = ZERO
            JOB = 00000
            SHORT = .FALSE.
            ISODR = .TRUE.
            MAXIT = 1000
            BETA(1:2) = (/ 200.0_R8, 3.0_R8 /)
            LOWER(1:2) = (/ 1.1_R8, 0.0_R8 /)
            UPPER(1:2) = (/ 400.0_R8, 6.0_R8 /)
            TSTTOL = 500.0_R8

         ELSE IF (ITEST.EQ.16) THEN

C  Test bounded odr problem where
C  initial BETA is outside bounds.

            LUN = LUNRPT
            WRITE (LUN,1000)
            DO I=1,2
               WRITE (LUN,1001) ITEST
               WRITE (LUN,1010)
               LUN = LUNSUM
            END DO
            SETNO = 9
            CALL DODRXD(TITLE,N,M,NP,NQ,LDX,X,LDY,Y,BETA)
            CALL DZERO(LWORK,1,WORK,LWORK)
            DELTA(:,:) = ZERO
            JOB = 00000
            SHORT = .FALSE.
            ISODR = .TRUE.
            MAXIT = 1000
            BETA(1:2) = (/ 200.0_R8, 7.0_R8 /)
            LOWER(1:2) = (/ 1.1_R8, 0.0_R8 /)
            UPPER(1:2) = (/ 200.0_R8, 5.0_R8 /)

         ELSE IF (ITEST.EQ.17) THEN

C  Test bounded odr problem where
C  bounds are ill defined.

            LUN = LUNRPT
            WRITE (LUN,1000)
            DO I=1,2
               WRITE (LUN,1001) ITEST
               WRITE (LUN,1010)
               LUN = LUNSUM
            END DO
            SETNO = 9
            CALL DODRXD(TITLE,N,M,NP,NQ,LDX,X,LDY,Y,BETA)
            CALL DZERO(LWORK,1,WORK,LWORK)
            DELTA(:,:) = ZERO
            JOB = 00000
            SHORT = .FALSE.
            ISODR = .TRUE.
            MAXIT = 1000
            BETA(1:2) = (/ 200.0_R8, 2.0_R8 /)
            LOWER(1:2) = (/ 10.0_R8, 0.0_R8 /)
            UPPER(1:2) = (/ 2.0_R8, 5.0_R8 /)

         ELSE IF (ITEST.EQ.18) THEN

C  Test bounded odr problem using centered differences where
C  parameters start on bound, move away, hit bound, move away, find minimum.

            LUN = LUNRPT
            WRITE (LUN,1000)
            DO I=1,2
               WRITE (LUN,1001) ITEST
               WRITE (LUN,1010)
               LUN = LUNSUM
            END DO
            SETNO = 9
            CALL DODRXD(TITLE,N,M,NP,NQ,LDX,X,LDY,Y,BETA)
            CALL DZERO(LWORK,1,WORK,LWORK)
            DELTA(:,:) = ZERO
            JOB = 00010
            SHORT = .FALSE.
            ISODR = .TRUE.
            MAXIT = 100
            BETA(1:2)  = (/ 200.0_R8, 5.0_R8 /)
            LOWER(1:2) = (/ 0.1_R8, 0.0_R8 /)
            UPPER(1:2) = (/ 200.0_R8, 5.0_R8 /)

         ELSE IF (ITEST.EQ.19) THEN

C  Test bounded odr problem when bounds are too small.
C  Parameters start on bound.

            LUN = LUNRPT
            WRITE (LUN,1000)
            DO I=1,2
               WRITE (LUN,1001) ITEST
               WRITE (LUN,1010)
               LUN = LUNSUM
            END DO
            SETNO = 9
            CALL DODRXD(TITLE,N,M,NP,NQ,LDX,X,LDY,Y,BETA)
            CALL DZERO(LWORK,1,WORK,LWORK)
            DELTA(:,:) = ZERO
            JOB = 00010
            SHORT = .FALSE.
            ISODR = .TRUE.
            MAXIT = 100
            BETA(1:2)  = (/ 200.0_R8, 5.0_R8 /)
            UPPER(1) = 200.0_R8
            LOWER(2) = 5.0_R8
            LOWER(1) = UPPER(1) - 400*UPPER(1)*EPSMAC 
     &                 + UPPER(1)*EPSMAC

            UPPER(2) = LOWER(2) + 400*LOWER(2)*EPSMAC 
     &                   - LOWER(2)*EPSMAC

         ELSE IF (ITEST.EQ.20) THEN

C  Test bounded odr problem when bounds are just big enough for ndigit
C  calculation but too small for difference calculation.
C  Parameters start on bound.

            LUN = LUNRPT
            WRITE (LUN,1000)
            DO I=1,2
               WRITE (LUN,1001) ITEST
               WRITE (LUN,1010)
               LUN = LUNSUM
            END DO
            SETNO = 9
            CALL DODRXD(TITLE,N,M,NP,NQ,LDX,X,LDY,Y,BETA)
            CALL DZERO(LWORK,1,WORK,LWORK)
            DELTA(:,:) = ZERO
            JOB = 00000
            SHORT = .FALSE.
            ISODR = .TRUE.
            MAXIT = 100
            BETA(1:2)  = (/ -200.0_R8, -5.0_R8 /)
            UPPER(1) = -200.0_R8
            LOWER(2) = -5.0_R8
            LOWER(1) = UPPER(1) + 400*UPPER(1)*EPSMAC
            UPPER(2) = LOWER(2) - 400*LOWER(2)*EPSMAC

         ELSE IF (ITEST.EQ.21) THEN

C  Test bounded odr problem when bounds are too small for derivative
C  step sizes using forward differences.  Parameters start on bound.

            LUN = LUNRPT
            WRITE (LUN,1000)
            DO I=1,2
               WRITE (LUN,1001) ITEST
               WRITE (LUN,1010)
               LUN = LUNSUM
            END DO
            SETNO = 9
            CALL DODRXD(TITLE,N,M,NP,NQ,LDX,X,LDY,Y,BETA)
            CALL DZERO(LWORK,1,WORK,LWORK)
            DELTA(:,:) = ZERO
            JOB = 00000
            SHORT = .FALSE.
            ISODR = .TRUE.
            MAXIT = 100
            BETA(1:2)  = (/ -200.0_R8, -5.0_R8 /)
            UPPER(1) = -200.0_R8
            LOWER(2) = -5.0_R8
            LOWER(1) = UPPER(1) + UPPER(1)*EPSMAC 
            UPPER(2) = LOWER(2) - LOWER(2)*EPSMAC 

         ELSE IF (ITEST.EQ.22) THEN

C  Test bounded odr problem when first parameter is fixed and second is bounded.
C  However, set the bounds on the first parameter to exclude the correct value
C  of the second parameter.  This will exercise the packing and unpacking of
C  parameters and ensure that bounds and fixed parameters can be mixed.

            LUN = LUNRPT
            WRITE (LUN,1000)
            DO I=1,2
               WRITE (LUN,1001) ITEST
               WRITE (LUN,1010)
               LUN = LUNSUM
            END DO
            SETNO = 10
            CALL DODRXD(TITLE,N,M,NP,NQ,LDX,X,LDY,Y,BETA)
            CALL DZERO(LWORK,1,WORK,LWORK)
            DELTA(:,:) = ZERO
            JOB = 00010
            SHORT = .FALSE.
            ISODR = .TRUE.
            MAXIT = 100
            BETA(1:2)  = (/ 2.5_R8, 1.5_R8 /)
            LOWER(1:2) = (/ 2.5_R8, 1.1_R8 /)
            UPPER(1:2) = (/ 10.0_R8, 5.0_R8 /)
            IFIXB(1:2) = (/ 0, 1 /)

         ELSE IF (ITEST.EQ.23) THEN

C  Similar to test 22 but without bounds.

            LUN = LUNRPT
            WRITE (LUN,1000)
            DO I=1,2
               WRITE (LUN,1001) ITEST
               WRITE (LUN,1010)
               LUN = LUNSUM
            END DO
            SETNO = 10
            CALL DODRXD(TITLE,N,M,NP,NQ,LDX,X,LDY,Y,BETA)
            CALL DZERO(LWORK,1,WORK,LWORK)
            DELTA(:,:) = ZERO
            JOB = 00010
            SHORT = .FALSE.
            ISODR = .TRUE.
            MAXIT = 100
            BETA(1:2)  = (/ 2.5_R8, 1.5_R8 /)
            LOWER(1:2) = -HUGE(1.0_R8)
            UPPER(1:2) = HUGE(1.0_R8)
            IFIXB(1:2) = (/ 0, 1 /)

         END IF

         CALL DODRXW
     &      (N,M,NP,NQ,LDWE1,LD2WE1,ISODR,LIWMIN,LWMIN)

C  Compute solution 

         WRITE (LUNRPT,2200) TITLE
         WRITE (LUNSUM,2200) TITLE
         IF (SHORT) THEN
            CALL ODR(FCN=DODRXF,
     &               N=N,M=M,NP=NP,NQ=NQ,
     &               BETA=BETA,
     &               Y=Y,X=X,
     &               DELTA=DELTA,
     &               WE=WE(1:LDWE1,1:LD2WE1,:),
     &               WD=WD(1:LDWD1,1:LD2WD1,:),
     &               JOB=JOB,
     &               IPRINT=IPRINT,LUNERR=LUNERR,LUNRPT=LUNRPT,
     &               WORK=WORK,IWORK=IWORK,
     &               INFO=INFO)
         ELSE
            CALL ODR(FCN=DODRXF,
     &               N=N,M=M,NP=NP,NQ=NQ,
     &               BETA=BETA,
     &               Y=Y,X=X,
     &               DELTA=DELTA,
     &               WE=WE(1:LDWE1,1:LD2WE1,:),
     &               WD=WD(1:LDWD1,1:LD2WD1,:),
     &               IFIXB=IFIXB,IFIXX=IFIXX(1:LDIFX,:),
     &               JOB=JOB,NDIGIT=NDIGIT,TAUFAC=TAUFAC,
     &               SSTOL=SSTOL,PARTOL=PARTOL,MAXIT=MAXIT,
     &               IPRINT=IPRINT,LUNERR=LUNERR,LUNRPT=LUNRPT,
     &               STPB=STPB,STPD=STPD(1:LDSTPD,:),
     &               SCLB=SCLB,SCLD=SCLD(1:LDSCLD,:),
     &               WORK=WORK,IWORK=IWORK,
     &               LOWER=LOWER(1:NP),UPPER=UPPER(1:NP),
     &               INFO=INFO)
         END IF

C  Compare results with those obtained on the cray ymp or the intel xeon running
C  Linux using REAL (KIND=R8) version of ODRPACK95

         BNRM = DNRM2(NP,BETA,1)
         CALL DWGHT(N,M,WD,LDWD1,LD2WD1,RESHAPE(WORK(1:N*M),(/N,M/)),
     &      TEMPRETL(1:N,1:M))
         WRK(1:N*M) = RESHAPE(TEMPRETL(1:N,1:M),(/N*M/))
         WSSDEL = DDOT(N*M,WORK(1:N*M),1,WRK(1),1)
         CALL DWGHT(N,NQ,WE,LDWE1,LD2WE1,
     &      RESHAPE(WORK(N*M+1:N*M+1+N*NQ-1),(/N,NQ/)),
     &      TEMPRETL(1:N,1:NQ))
         WRK(N*M+1:N*M+1+N*NQ-1) = RESHAPE(TEMPRETL(1:N,1:NQ),(/N*NQ/))
         WSSEPS = DDOT(N*NQ,WORK(N*M+1:N*M+1+N*NQ-1),1,
     &      WRK(N*M+1:N*M+1+N*NQ-1),1)
         WSS = WSSEPS + WSSDEL

         IF (SSTOL.LT.ZERO) THEN
            SSTOL = SQRT(EPSMAC)
         ELSE
            SSTOL = MIN(SSTOL, ONE)
         END IF

         IF (PARTOL.LT.ZERO) THEN
            PARTOL = EPSMAC**(TWO/THREE)
         ELSE
            PARTOL = MIN(PARTOL, ONE)
         END IF

         IF (INFO.GE.10000) THEN
            IF (IDPYMP(ITEST).EQ.INFO) THEN
                  FAILS = .FALSE.
                  MSG = 1
            ELSE
                  FAILS = .TRUE.
                  MSG = 3
            END IF

         ELSE IF (MOD(INFO,10).EQ.1) THEN
                  FAILS = ABS(WSS-DPYMP(2,ITEST)).GT.
     &                    DPYMP(2,ITEST)*SSTOL*TSTTOL
                  MSG = 2

         ELSE IF (MOD(INFO,10).EQ.2) THEN
                  FAILS = ABS(BNRM-DPYMP(1,ITEST)).GT.
     &                    DPYMP(1,ITEST)*PARTOL*TSTTOL
                  MSG = 2

         ELSE IF (MOD(INFO,10).EQ.3) THEN
                  FAILS = (ABS(WSS-DPYMP(2,ITEST)).GT.
     &                     DPYMP(2,ITEST)*SSTOL*TSTTOL)
     &                    .AND.
     &                    (ABS(BNRM-DPYMP(1,ITEST)).GT.
     &                     DPYMP(1,ITEST)*PARTOL*TSTTOL)
                  MSG = 2

         ELSE IF ((MOD(INFO,10).EQ.4) .AND. (IDPYMP(ITEST).EQ.4)) THEN
                  FAILS = .FALSE.
                  MSG = 1

         ELSE IF (INFO.EQ.IDPYMP(ITEST)) THEN
                  FAILS = .TRUE.
                  MSG = 4
         ELSE
                  FAILS = .TRUE.
                  MSG = 3
         END IF

         FAILED = FAILED .OR. FAILS

         LUN = LUNRPT
         DO 300 L=1,2
            WRITE (LUN,3100)
            WRITE (LUN,3210) 
     &         ' CRAY YMP OR X86 RESULT = ',
     &         DPYMP(1,ITEST),DPYMP(2,ITEST),IDPYMP(ITEST)
            WRITE (LUN,3210) ' NEW TEST RESULT      = ',
     &         BNRM,WSS,INFO
            WRITE (LUN,3220) ' DIFFERENCE           = ',
     &         ABS(DPYMP(1,ITEST)-BNRM),ABS(DPYMP(2,ITEST)-WSS)
            EWRT  = ABS(DPYMP(1,ITEST))
            EWRT2 = ABS(DPYMP(2,ITEST))
            IF (EWRT.EQ.ZERO) THEN
               EWRT = ONE
            END IF
            IF (EWRT2.EQ.ZERO) THEN
               EWRT2 = ONE
            END IF
            WRITE (LUN,3220) ' RELATIVE ERROR       = ',
     &         ABS(DPYMP(1,ITEST)-BNRM)/EWRT,
     &         ABS(DPYMP(2,ITEST)-WSS)/EWRT2

            IF (MSG.EQ.1) THEN
               WRITE (LUN,3310)
            ELSE IF (MSG.EQ.2) THEN
               IF (FAILS) THEN
                  WRITE (LUN,3320)
               ELSE
                  WRITE (LUN,3330)
               END IF
            ELSE IF (MSG.EQ.3) THEN
               WRITE (LUN,3340)
            ELSE IF (MSG.EQ.4) THEN
               WRITE (LUN,3350)
            END IF

            LUN = LUNSUM
  300    CONTINUE
  400 CONTINUE

      WRITE (LUNRPT,1000)
      IF (FAILED) THEN
         WRITE (LUNRPT,4100)
         WRITE (LUNSUM,4100)
         PASSED = .FALSE.
      ELSE
         WRITE (LUNRPT,4200)
         WRITE (LUNSUM,4200)
         PASSED = .TRUE.
      END IF

C  Format statements

 1000 FORMAT('1')
 1001 FORMAT(' Example ', I2/)
 1010 FORMAT(' Test simple odr problem'/
     &       ' with analytic derivatives',
     &       ' using ODR.')
 1020 FORMAT(' Test simple OLS problem'/
     &       ' with finite difference derivatives',
     &       ' using ODR.')
 1030 FORMAT(' Test parameter fixing capabilities',
     &       ' for poorly scaled OLS problem'/
     &       ' with analytic derivatives',
     &       ' using ODR.')
 1040 FORMAT(' Test weighting capabilities',
     &       ' for ODR problem'/
     &       ' with analytic derivatives',
     &       ' using ODR. '/
     &       ' also shows solution of poorly scaled',
     &       ' ODR problem.'/
     &       ' (derivative checking turned off.)')
 1050 FORMAT(' Test DELTA initialization capabilities'/
     &       ' and use of ISTOP to restrict parameter values',
     &       ' for ODR problem'/
     &       ' with analytic derivatives',
     &       ' using ODR.')
 1060 FORMAT(' Test stiff stopping conditions',
     &       ' for unscaled ODR problem'/
     &       ' with analytic derivatives',
     &       ' using ODR.')
 1070 FORMAT(' Test restart',
     &       ' for unscaled ODR problem'/
     &       ' with analytic derivatives',
     &       ' using ODR.')
 1080 FORMAT(' Test use of TAUFAC to restrict first step',
     &       ' for ODR problem'/
     &       ' with finite difference derivatives',
     &       ' using ODR.')
 1090 FORMAT(' Test implicit model',
     &       ' for OLS problem'/
     &       ' using ODR.')
 1100 FORMAT(' Test multiresponse model',
     &       ' for ODR problem'/
     &       ' with finite difference derivatives',
     &       ' using ODR.')
 1110 FORMAT(' Test detection of questionable analytic derivatives',
     &       ' for OLS problem'/
     &       ' using ODR.')
 1120 FORMAT(' Test detection of incorrect analytic derivatives',
     &       ' for ODR problem'/
     &       ' with analytic derivatives',
     &       ' using ODR.')
 2200 FORMAT (' Data Set Reference: ', A80)
 3100 FORMAT
     &   (/' Comparison of new results with',
     &     ' REAL (KIND=R8) Cray YMP or Intel X86 (Linux) '/
     &     ' Result:'//
     &     '                         Norm of BETA',
     &     '        Sum of Squared WTD OBS Errors  INFO')
 3210 FORMAT
     &   (/A25/1P,2E37.30,I6)
 3220 FORMAT
     &   (/A25,1P,D12.5,25X,D12.5,I6)
 3310 FORMAT
     &   (/' *** Stopping conditions',
     &         ' show convergence not attained. ***'/
     &     '        no further comparisons made between results.'//)
 3320 FORMAT
     &   (//' *** WARNING ***',
     &      ' results do not agree to within stopping tolerance. ***'//)
 3330 FORMAT
     &   (//' *** Results agree to within stopping tolerance. ***'//)
 3340 FORMAT
     &   (//' *** WARNING ***',
     &     ' stopping conditions do not agree. ***'//)
 3350 FORMAT
     &   (//' *** WARNING ***',
     &      ' unexpected stopping condition.',
     &      '  please contact package authors. ***'//)
 4100 FORMAT
     &   (///
     &   ' *** Summary:',
     &   ' one or more tests do not agree with expected results. ***')
 4200 FORMAT
     &   (///
     &   ' *** Summary:',
     &   ' all tests agree with expected results. ***')

      END
*DODRXD
      SUBROUTINE DODRXD
     &   (TITLE,N,M,NP,NQ,LDX,X,LDY,Y,BETA)
C***Begin Prologue  DODRXD
C***Refer to  ODR
C***Routines Called  (NONE)
C***Date Written   860529   (YYMMDD)
C***Revision Date  920619   (YYMMDD)
C***Purpose  Set up data for ODRPACK95 exerciser
C***End Prologue  DODRXD

C...Used modules
      USE REAL_PRECISION

C...Parameters
      INTEGER
     &    MAXN,MAXM,MAXNP,MAXNQ,MAXSET
      PARAMETER
     &    (MAXN=50,MAXM=3,MAXNP=10,MAXNQ=3,MAXSET=16)

C...Scalar arguments
      INTEGER
     &   LDX,LDY,M,N,NP,NQ
      CHARACTER TITLE*80

C...Array arguments
      REAL (KIND=R8)
     &   BETA(*),X(LDX,*),Y(LDY,*)

C...Scalars in common
      INTEGER
     &   SETNO

C...Local scalars
      INTEGER
     &   I,J,K,L

C...Local arrays
      REAL (KIND=R8)
     &   BDATA(MAXNP,MAXSET),XDATA(MAXN,MAXM,MAXSET),
     &   YDATA(MAXN,MAXNQ,MAXSET)
      INTEGER
     &   MDATA(MAXSET),NDATA(MAXSET),NPDATA(MAXSET),NQDATA(MAXSET)
      CHARACTER TDATA(MAXSET)*80

C...Common blocks
      COMMON /SETID/SETNO

C...Data statements
      DATA
     &   TDATA(1)
     &   /' BOGGS, BYRD AND SCHNABEL, 1985, EXAMPLE 1'/
      DATA
     &   NDATA(1), MDATA(1), NPDATA(1), NQDATA(1)
     &   /40, 1, 2, 1/
      DATA
     &   (BDATA(K,1),K=1,2)
     &   /1.0E+0_R8, 1.0E+0_R8/
      DATA
     &   YDATA( 1,1,1), XDATA( 1,1,1)
     &   /-0.119569795672791172E+1_R8, -0.213701920211315155E-1_R8/
      DATA
     &   YDATA( 2,1,1), XDATA( 2,1,1)
     &   /-0.128023349509594288E+1_R8, 0.494813247025012969E-1_R8/
      DATA
     &   YDATA( 3,1,1), XDATA( 3,1,1)
     &   /-0.125270693343174591E+1_R8, 0.127889194935560226E+0_R8/
      DATA
     &   YDATA( 4,1,1), XDATA( 4,1,1)
     &   /-0.996698267935287383E+0_R8, 0.128615394085645676E+0_R8/
      DATA
     &   YDATA( 5,1,1), XDATA( 5,1,1)
     &   /-0.104681033065801934E+1_R8, 0.232544285655021667E+0_R8/
      DATA
     &   YDATA( 6,1,1), XDATA( 6,1,1)
     &   /-0.146724952092847308E+1_R8, 0.268151108026504516E+0_R8/
      DATA
     &   YDATA( 7,1,1), XDATA( 7,1,1)
     &   /-0.123366891873487528E+1_R8, 0.309041029810905456E+0_R8/
      DATA
     &   YDATA( 8,1,1), XDATA( 8,1,1)
     &   /-0.165665097907185554E+1_R8, 0.405991539210081099E+0_R8/
      DATA
     &   YDATA( 9,1,1), XDATA( 9,1,1)
     &   /-0.168476460930907119E+1_R8, 0.376611424833536147E+0_R8/
      DATA
     &   YDATA(10,1,1), XDATA(10,1,1)
     &   /-0.198571971169224491E+1_R8, 0.475875890851020811E+0_R8/
      DATA
     &   YDATA(11,1,1), XDATA(11,1,1)
     &   /-0.195691696638051344E+1_R8, 0.499246935397386550E+0_R8/
      DATA
     &   YDATA(12,1,1), XDATA(12,1,1)
     &   /-0.211871342665769836E+1_R8, 0.536615037024021147E+0_R8/
      DATA
     &   YDATA(13,1,1), XDATA(13,1,1)
     &   /-0.268642932558671020E+1_R8, 0.581830765902996060E+0_R8/
      DATA
     &   YDATA(14,1,1), XDATA(14,1,1)
     &   /-0.281123260058024347E+1_R8, 0.684512710422277446E+0_R8/
      DATA
     &   YDATA(15,1,1), XDATA(15,1,1)
     &   /-0.328704486581785920E+1_R8, 0.660219819694757458E+0_R8/
      DATA
     &   YDATA(16,1,1), XDATA(16,1,1)
     &   /-0.423062993461887032E+1_R8, 0.766990323960781092E+0_R8/
      DATA
     &   YDATA(17,1,1), XDATA(17,1,1)
     &   /-0.512043906552226903E+1_R8, 0.808270426690578456E+0_R8/
      DATA
     &   YDATA(18,1,1), XDATA(18,1,1)
     &   /-0.731032616379005535E+1_R8, 0.897410020083189004E+0_R8/
      DATA
     &   YDATA(19,1,1), XDATA(19,1,1)
     &   /-0.109002759485608993E+2_R8, 0.959199774116277687E+0_R8/
      DATA
     &   YDATA(20,1,1), XDATA(20,1,1)
     &   /-0.251810238510370206E+2_R8, 0.914675474762916558E+0_R8/
      DATA
     &   YDATA(21,1,1), XDATA(21,1,1)
     &   /0.100123028650879944E+3_R8, 0.997759691476821892E+0_R8/
      DATA
     &   YDATA(22,1,1), XDATA(22,1,1)
     &   /0.168225085871915048E+2_R8, 0.107136870384216308E+1_R8/
      DATA
     &   YDATA(23,1,1), XDATA(23,1,1)
     &   /0.894830510866913009E+1_R8, 0.108033321037888526E+1_R8/
      DATA
     &   YDATA(24,1,1), XDATA(24,1,1)
     &   /0.645853815227747004E+1_R8, 0.116064198672771453E+1_R8/
      DATA
     &   YDATA(25,1,1), XDATA(25,1,1)
     &   /0.498218564760117328E+1_R8, 0.119080889359116553E+1_R8/
      DATA
     &   YDATA(26,1,1), XDATA(26,1,1)
     &   /0.382971664718710476E+1_R8, 0.129418875187635420E+1_R8/
      DATA
     &   YDATA(27,1,1), XDATA(27,1,1)
     &   /0.344116492497344184E+1_R8, 0.135594148099422453E+1_R8/
      DATA
     &   YDATA(28,1,1), XDATA(28,1,1)
     &   /0.276840496973858949E+1_R8, 0.135302808716893195E+1_R8/
      DATA
     &   YDATA(29,1,1), XDATA(29,1,1)
     &   /0.259521665196956666E+1_R8, 0.137994666010141371E+1_R8/
      DATA
     &   YDATA(30,1,1), XDATA(30,1,1)
     &   /0.205996022794557661E+1_R8, 0.147630019545555113E+1_R8/
      DATA
     &   YDATA(31,1,1), XDATA(31,1,1)
     &   /0.197939614345337836E+1_R8, 0.153450708076357840E+1_R8/
      DATA
     &   YDATA(32,1,1), XDATA(32,1,1)
     &   /0.156739340562905589E+1_R8, 0.152805351451039313E+1_R8/
      DATA
     &   YDATA(33,1,1), XDATA(33,1,1)
     &   /0.159032057073028366E+1_R8, 0.157147316247224806E+1_R8/
      DATA
     &   YDATA(34,1,1), XDATA(34,1,1)
     &   /0.173102268158937949E+1_R8, 0.166649596005678175E+1_R8/
      DATA
     &   YDATA(35,1,1), XDATA(35,1,1)
     &   /0.155512561664824758E+1_R8, 0.166505665838718412E+1_R8/
      DATA
     &   YDATA(36,1,1), XDATA(36,1,1)
     &   /0.149635994944133260E+1_R8, 0.175214128553867338E+1_R8/
      DATA
     &   YDATA(37,1,1), XDATA(37,1,1)
     &   /0.147487601463073568E+1_R8, 0.180567992463707922E+1_R8/
      DATA
     &   YDATA(38,1,1), XDATA(38,1,1)
     &   /0.117244575233306998E+1_R8, 0.184624404296278952E+1_R8/
      DATA
     &   YDATA(39,1,1), XDATA(39,1,1)
     &   /0.910931336069172580E+0_R8, 0.195568727388978002E+1_R8/
      DATA
     &   YDATA(40,1,1), XDATA(40,1,1)
     &   /0.126172980914513272E+1_R8, 0.199326394036412237E+1_R8/

      DATA
     &   TDATA(2)
     &   /' BOGGS, BYRD AND SCHNABEL, 1985, EXAMPLE 2'/
      DATA
     &   NDATA(2), MDATA(2), NPDATA(2), NQDATA(2)
     &   /50, 2, 3, 1/
      DATA
     &   (BDATA(K,2),K=1,3)
     &   /-1.0E+0_R8, 1.0E+0_R8, 1.0E+0_R8/
      DATA
     &   YDATA( 1,1,2), XDATA( 1,1,2), XDATA( 1,2,2)
     &   /0.680832777217942900E+0_R8,
     &   0.625474598833994800E-1_R8, 0.110179064209783100E+0_R8/
      DATA
     &   YDATA( 2,1,2), XDATA( 2,1,2), XDATA( 2,2,2)
     &   /0.122183594595302200E+1_R8,
     &   0.202500343620642400E+0_R8, -0.196140862891327600E-1_R8/
      DATA
     &   YDATA( 3,1,2), XDATA( 3,1,2), XDATA( 3,2,2)
     &   /0.118958678734608200E+1_R8,
     &   0.164943738599876500E+0_R8, 0.166514874750996600E+0_R8/
      DATA
     &   YDATA( 4,1,2), XDATA( 4,1,2), XDATA( 4,2,2)
     &   /0.146982623764094600E+1_R8,
     &   0.304874137610506100E+0_R8, 0.612908688041490500E-2_R8/
      DATA
     &   YDATA( 5,1,2), XDATA( 5,1,2), XDATA( 5,2,2)
     &   /0.167775338189355300E+1_R8,
     &   0.532727445580665100E+0_R8, 0.938248787552444600E-1_R8/
      DATA
     &   YDATA( 6,1,2), XDATA( 6,1,2), XDATA( 6,2,2)
     &   /0.202485721906026200E+1_R8,
     &   0.508823707598910200E+0_R8, 0.499605775020505400E-2_R8/
      DATA
     &   YDATA( 7,1,2), XDATA( 7,1,2), XDATA( 7,2,2)
     &   /0.258912851935938800E+1_R8,
     &   0.704227041878554000E+0_R8, 0.819354849092326200E-1_R8/
      DATA
     &   YDATA( 8,1,2), XDATA( 8,1,2), XDATA( 8,2,2)
     &   /0.366894203254154800E+1_R8,
     &   0.592077736111512000E+0_R8, 0.127113960672389100E-1_R8/
      DATA
     &   YDATA( 9,1,2), XDATA( 9,1,2), XDATA( 9,2,2)
     &   /0.574609583351347300E+1_R8,
     &   0.104940945646421600E+1_R8, 0.258095243658316100E-1_R8/
      DATA
     &   YDATA(10,1,2), XDATA(10,1,2), XDATA(10,2,2)
     &   /0.127676424026489300E+2_R8,
     &    0.979382517558619200E+0_R8, 0.124280755181027900E+0_R8/
      DATA
     &   YDATA(11,1,2), XDATA(11,1,2), XDATA(11,2,2)
     &   /0.123473079693623100E+1_R8,
     &    0.637870453165538700E-1_R8, 0.304856401137196400E+0_R8/
      DATA
     &   YDATA(12,1,2), XDATA(12,1,2), XDATA(12,2,2)
     &   /0.142256120864082800E+1_R8,
     &    0.176123312906025700E+0_R8, 0.262387028078896900E+0_R8/
      DATA
     &   YDATA(13,1,2), XDATA(13,1,2), XDATA(13,2,2)
     &   /0.169889534013024700E+1_R8,
     &    0.310965082300263000E+0_R8, 0.226430765474758800E+0_R8/
      DATA
     &   YDATA(14,1,2), XDATA(14,1,2), XDATA(14,2,2)
     &   /0.173485577901204400E+1_R8,
     &    0.311394269116782100E+0_R8, 0.271375840410281800E+0_R8/
      DATA
     &   YDATA(15,1,2), XDATA(15,1,2), XDATA(15,2,2)
     &   /0.277761263972834600E+1_R8,
     &    0.447076126190612500E+0_R8, 0.255000858902618300E+0_R8/
      DATA
     &   YDATA(16,1,2), XDATA(16,1,2), XDATA(16,2,2)
     &   /0.339163324662617300E+1_R8,
     &    0.384786230998211100E+0_R8, 0.154958003178364000E+0_R8/
      DATA
     &   YDATA(17,1,2), XDATA(17,1,2), XDATA(17,2,2)
     &   /0.589615137312147500E+1_R8,
     &    0.649093176450780500E+0_R8, 0.258301685463773200E+0_R8/
      DATA
     &   YDATA(18,1,2), XDATA(18,1,2), XDATA(18,2,2)
     &   /0.124415625214576800E+2_R8,
     &    0.685612005372525500E+0_R8, 0.107391260603228600E+0_R8/
      DATA
     &   YDATA(19,1,2), XDATA(19,1,2), XDATA(19,2,2)
     &   /-0.498491739153861600E+2_R8,
     &    0.968747139425088400E+0_R8, 0.151932526135740700E+0_R8/
      DATA
     &   YDATA(20,1,2), XDATA(20,1,2), XDATA(20,2,2)
     &   /-0.832795509000618600E+1_R8,
     &    0.869789367989532900E+0_R8, 0.625507500586400000E-1_R8/
      DATA
     &   YDATA(21,1,2), XDATA(21,1,2), XDATA(21,2,2)
     &   /0.184934617774239900E+1_R8,
     &    -0.465309930332736600E-2_R8, 0.546795662595375200E+0_R8/
      DATA
     &   YDATA(22,1,2), XDATA(22,1,2), XDATA(22,2,2)
     &   /0.175192979176839200E+1_R8,
     &    0.604753397196646000E-2_R8, 0.230905749473922700E+0_R8/
      DATA
     &   YDATA(23,1,2), XDATA(23,1,2), XDATA(23,2,2)
     &   /0.253949381238535800E+1_R8,
     &    0.239418809621756000E+0_R8, 0.190752069681170700E+0_R8/
      DATA
     &   YDATA(24,1,2), XDATA(24,1,2), XDATA(24,2,2)
     &   /0.373500774928501700E+1_R8,
     &    0.456662468911699800E+0_R8, 0.328870615170984400E+0_R8/
      DATA
     &   YDATA(25,1,2), XDATA(25,1,2), XDATA(25,2,2)
     &   /0.548408128950331000E+1_R8,
     &    0.371115320522079500E+0_R8, 0.439978556640660500E+0_R8/
      DATA
     &   YDATA(26,1,2), XDATA(26,1,2), XDATA(26,2,2)
     &   /0.125256880521774300E+2_R8,
     &    0.586442107042503000E+0_R8, 0.490689043752286700E+0_R8/
      DATA
     &   YDATA(27,1,2), XDATA(27,1,2), XDATA(27,2,2)
     &   /-0.493587797164916600E+2_R8,
     &    0.579796274973298000E+0_R8, 0.521860998203383100E+0_R8/
      DATA
     &   YDATA(28,1,2), XDATA(28,1,2), XDATA(28,2,2)
     &   /-0.801158974965412700E+1_R8,
     &    0.805008094903899900E+0_R8, 0.292283538955391600E+0_R8/
      DATA
     &   YDATA(29,1,2), XDATA(29,1,2), XDATA(29,2,2)
     &   /-0.437399487061934100E+1_R8,
     &    0.637242340835710000E+0_R8, 0.402261740352486000E+0_R8/
      DATA
     &   YDATA(30,1,2), XDATA(30,1,2), XDATA(30,2,2)
     &   /-0.297800103425979600E+1_R8,
     &    0.982132817936118700E+0_R8, 0.392546836419047000E+0_R8/
      DATA
     &   YDATA(31,1,2), XDATA(31,1,2), XDATA(31,2,2)
     &   /0.271811057454661300E+1_R8,
     &    -0.223515657121262700E-1_R8, 0.650479019708978800E+0_R8/
      DATA
     &   YDATA(32,1,2), XDATA(32,1,2), XDATA(32,2,2)
     &   /0.377035865613392400E+1_R8,
     &    0.136081427545033600E+0_R8, 0.753020101897661800E+0_R8/
      DATA
     &   YDATA(33,1,2), XDATA(33,1,2), XDATA(33,2,2)
     &   /0.560111053917143100E+1_R8,
     &    0.145367053019870600E+0_R8, 0.611153532003093100E+0_R8/
      DATA
     &   YDATA(34,1,2), XDATA(34,1,2), XDATA(34,2,2)
     &   /0.128152376174926800E+2_R8,
     &    0.308221919576435500E+0_R8, 0.455217283290423900E+0_R8/
      DATA
     &   YDATA(35,1,2), XDATA(35,1,2), XDATA(35,2,2)
     &   /-0.498709177732467200E+2_R8,
     &    0.432658769133528300E+0_R8, 0.678607663414113000E+0_R8/
      DATA
     &   YDATA(36,1,2), XDATA(36,1,2), XDATA(36,2,2)
     &   /-0.815797696908314300E+1_R8,
     &    0.477785501079980300E+0_R8, 0.536178207572157000E+0_R8/
      DATA
     &   YDATA(37,1,2), XDATA(37,1,2), XDATA(37,2,2)
     &   /-0.440240491195158600E+1_R8,
     &    0.727986827616619000E+0_R8, 0.668497920573493900E+0_R8/
      DATA
     &   YDATA(38,1,2), XDATA(38,1,2), XDATA(38,2,2)
     &   /-0.276723957061767500E+1_R8,
     &    0.745950385588265100E+0_R8, 0.786077589007263700E+0_R8/
      DATA
     &   YDATA(39,1,2), XDATA(39,1,2), XDATA(39,2,2)
     &   /-0.223203667288734800E+1_R8,
     &    0.732537503527113500E+0_R8, 0.582625164046828400E+0_R8/
      DATA
     &   YDATA(40,1,2), XDATA(40,1,2), XDATA(40,2,2)
     &   /-0.169728270310622000E+1_R8,
     &    0.967352361433846300E+0_R8, 0.460779396016832800E+0_R8/
      DATA
     &   YDATA(41,1,2), XDATA(41,1,2), XDATA(41,2,2)
     &   /0.551015652153227000E+1_R8,
     &    0.129761784310891100E-1_R8, 0.700009537931860000E+0_R8/
      DATA
     &   YDATA(42,1,2), XDATA(42,1,2), XDATA(42,2,2)
     &   /0.128036180496215800E+2_R8,
     &    0.170163243950629700E+0_R8, 0.853131830764348700E+0_R8/
      DATA
     &   YDATA(43,1,2), XDATA(43,1,2), XDATA(43,2,2)
     &   /-0.498257683396339000E+2_R8,
     &    0.162768461906274000E+0_R8, 0.865315129048175000E+0_R8/
      DATA
     &   YDATA(44,1,2), XDATA(44,1,2), XDATA(44,2,2)
     &   /-0.877334550221761900E+1_R8,
     &    0.222914807946165800E+0_R8, 0.797511758502094500E+0_R8/
      DATA
     &   YDATA(45,1,2), XDATA(45,1,2), XDATA(45,2,2)
     &   /-0.453820192156867600E+1_R8,
     &    0.402910095604624900E+0_R8, 0.761492958727023100E+0_R8/
      DATA
     &   YDATA(46,1,2), XDATA(46,1,2), XDATA(46,2,2)
     &   /-0.297499315738677900E+1_R8,
     &    0.233770812593443200E+0_R8, 0.896000095844223500E+0_R8/
      DATA
     &   YDATA(47,1,2), XDATA(47,1,2), XDATA(47,2,2)
     &   /-0.212743255978538900E+1_R8,
     &    0.646528693486914700E+0_R8, 0.968574333700755700E+0_R8/
      DATA
     &   YDATA(48,1,2), XDATA(48,1,2), XDATA(48,2,2)
     &   /-0.209703205365401000E+1_R8,
     &    0.802811658568969400E+0_R8, 0.904866450476711600E+0_R8/
      DATA
     &   YDATA(49,1,2), XDATA(49,1,2), XDATA(49,2,2)
     &   /-0.155287292042086200E+1_R8,
     &    0.837137859891222900E+0_R8, 0.835684424990021900E+0_R8/
      DATA
     &   YDATA(50,1,2), XDATA(50,1,2), XDATA(50,2,2)
     &   /-0.161356673770480700E+1_R8,
     &    0.103165980756526600E+1_R8, 0.793902191912346100E+0_R8/

      DATA
     &   TDATA(3)
     &   /' BOGGS, BYRD AND SCHNABEL, 1985, EXAMPLE 3'/
      DATA
     &   NDATA(3), MDATA(3), NPDATA(3), NQDATA(3)
     &   /44, 1, 9, 1/
      DATA
     &   (BDATA(K,3),K=1,9)
     &   /0.281887509408440189E-5_R8,
     &   -0.231290549212363845E-2_R8, 0.583035555572801965E+1_R8,
     &    0.000000000000000000E+0_R8, 0.406910776203121026E+8_R8,
     &    0.138001105225000000E-2_R8, 0.596038513209999999E-1_R8,
     &    0.670582099359999998E+1_R8, 0.106994410000000000E+10_R8/
      DATA
     &   YDATA( 1,1,3), XDATA( 1,1,3)
     &   /0.988227696721327788E+0_R8, 0.25E-8_R8/
      DATA
     &   YDATA( 2,1,3), XDATA( 2,1,3)
     &   /0.988268083998559958E+0_R8, 0.64E-8_R8/
      DATA
     &   YDATA( 3,1,3), XDATA( 3,1,3)
     &   /0.988341022958438831E+0_R8, 1.0E-8_R8/
      DATA
     &   YDATA( 4,1,3), XDATA( 4,1,3)
     &   /0.988380557606306446E+0_R8, 0.9E-7_R8/
      DATA
     &   YDATA( 5,1,3), XDATA( 5,1,3)
     &   /0.988275062411751338E+0_R8, 1.0E-6_R8/
      DATA
     &   YDATA( 6,1,3), XDATA( 6,1,3)
     &   /0.988326680176446987E+0_R8, 0.4E-5_R8/
      DATA
     &   YDATA( 7,1,3), XDATA( 7,1,3)
     &   /0.988306058860433439E+0_R8, 0.9E-5_R8/
      DATA
     &   YDATA( 8,1,3), XDATA( 8,1,3)
     &   /0.988292880079125555E+0_R8, 0.16E-4_R8/
      DATA
     &   YDATA( 9,1,3), XDATA( 9,1,3)
     &   /0.988305279259496905E+0_R8, 0.36E-4_R8/
      DATA
     &   YDATA(10,1,3), XDATA(10,1,3)
     &   /0.988278142019574202E+0_R8, 0.64E-4_R8/
      DATA
     &   YDATA(11,1,3), XDATA(11,1,3)
     &   /0.988224953369819946E+0_R8, 1.0E-4_R8/
      DATA
     &   YDATA(12,1,3), XDATA(12,1,3)
     &   /0.988111989169778223E+0_R8, 0.144E-3_R8/
      DATA
     &   YDATA(13,1,3), XDATA(13,1,3)
     &   /0.988045627103840613E+0_R8, 0.225E-3_R8/
      DATA
     &   YDATA(14,1,3), XDATA(14,1,3)
     &   /0.987913715667047655E+0_R8, 0.400E-3_R8/
      DATA
     &   YDATA(15,1,3), XDATA(15,1,3)
     &   /0.987841994238525678E+0_R8, 0.625E-3_R8/
      DATA
     &   YDATA(16,1,3), XDATA(16,1,3)
     &   /0.987638450432434270E+0_R8, 0.900E-3_R8/
      DATA
     &   YDATA(17,1,3), XDATA(17,1,3)
     &   /0.987587364331771395E+0_R8, 0.1225E-2_R8/
      DATA
     &   YDATA(18,1,3), XDATA(18,1,3)
     &   /0.987576264149633684E+0_R8, 0.1600E-2_R8/
      DATA
     &   YDATA(19,1,3), XDATA(19,1,3)
     &   /0.987539209110983643E+0_R8, 0.2025E-2_R8/
      DATA
     &   YDATA(20,1,3), XDATA(20,1,3)
     &   /0.987621143807705698E+0_R8, 0.25E-2_R8/
      DATA
     &   YDATA(21,1,3), XDATA(21,1,3)
     &   /0.988023229785526217E+0_R8, 0.36E-2_R8/
      DATA
     &   YDATA(22,1,3), XDATA(22,1,3)
     &   /0.988558376710994197E+0_R8, 0.49E-2_R8/
      DATA
     &   YDATA(23,1,3), XDATA(23,1,3)
     &   /0.989304775352439885E+0_R8, 0.64E-2_R8/
      DATA
     &   YDATA(24,1,3), XDATA(24,1,3)
     &   /0.990210452265710472E+0_R8, 0.81E-2_R8/
      DATA
     &   YDATA(25,1,3), XDATA(25,1,3)
     &   /0.991095950592263900E+0_R8, 1.00E-2_R8/
      DATA
     &   YDATA(26,1,3), XDATA(26,1,3)
     &   /0.991475677297119272E+0_R8, 0.11025E-1_R8/
      DATA
     &   YDATA(27,1,3), XDATA(27,1,3)
     &   /0.991901306250746771E+0_R8, 0.12100E-1_R8/
      DATA
     &   YDATA(28,1,3), XDATA(28,1,3)
     &   /0.992619222425303263E+0_R8, 0.14400E-1_R8/
      DATA
     &   YDATA(29,1,3), XDATA(29,1,3)
     &   /0.993617037631973475E+0_R8, 0.16900E-1_R8/
      DATA
     &   YDATA(30,1,3), XDATA(30,1,3)
     &   /0.994727321698030676E+0_R8, 0.19600E-1_R8/
      DATA
     &   YDATA(31,1,3), XDATA(31,1,3)
     &   /0.996523114720326189E+0_R8, 0.25600E-1_R8/
      DATA
     &   YDATA(32,1,3), XDATA(32,1,3)
     &   /0.998036909563764020E+0_R8, 0.32400E-1_R8/
      DATA
     &   YDATA(33,1,3), XDATA(33,1,3)
     &   /0.999151968626971372E+0_R8, 0.40000E-1_R8/
      DATA
     &   YDATA(34,1,3), XDATA(34,1,3)
     &   /0.100017083706131769E+1_R8, 0.50625E-1_R8/
      DATA
     &   YDATA(35,1,3), XDATA(35,1,3)
     &   /0.100110046382923523E+1_R8, 0.75625E-1_R8/
      DATA
     &   YDATA(36,1,3), XDATA(36,1,3)
     &   /0.100059103180404652E+1_R8, 0.12250E+0_R8/
      DATA
     &   YDATA(37,1,3), XDATA(37,1,3)
     &   /0.999211829791257561E+0_R8, 0.16000E+0_R8/
      DATA
     &   YDATA(38,1,3), XDATA(38,1,3)
     &   /0.994711451526761862E+0_R8, 0.25000E+0_R8/
      DATA
     &   YDATA(39,1,3), XDATA(39,1,3)
     &   /0.989844132928847109E+0_R8, 0.33640E+0_R8/
      DATA
     &   YDATA(40,1,3), XDATA(40,1,3)
     &   /0.987234104554490439E+0_R8, 0.38440E+0_R8/
      DATA
     &   YDATA(41,1,3), XDATA(41,1,3)
     &   /0.980928240178404887E+0_R8, 0.49E+0_R8/
      DATA
     &   YDATA(42,1,3), XDATA(42,1,3)
     &   /0.970888680366055576E+0_R8, 0.64E+0_R8/
      DATA
     &   YDATA(43,1,3), XDATA(43,1,3)
     &   /0.960043769857327398E+0_R8, 0.81E+0_R8/
      DATA
     &   YDATA(44,1,3), XDATA(44,1,3)
     &   /0.947277159259551068E+0_R8, 1.00E+0_R8/

      DATA
     &   TDATA(4)
     &   /' HIMMELBLAU, 1970, EXAMPLE 6.2-4, PAGE 188'/
      DATA
     &   NDATA(4), MDATA(4), NPDATA(4), NQDATA(4)
     &   /13, 2, 3, 1/
      DATA
     &   (BDATA(K,4),K=1,3)
     &   /3.0E+0_R8, 3.0E+0_R8, -0.5E+0_R8/
      DATA
     &   YDATA( 1,1,4), XDATA( 1,1,4), XDATA( 1,2,4)
     &   /2.93E+0_R8, 0.0E+0_R8, 0.0E+0_R8/
      DATA
     &   YDATA( 2,1,4), XDATA( 2,1,4), XDATA( 2,2,4)
     &   /1.95E+0_R8, 0.0E+0_R8, 1.0E+0_R8/
      DATA
     &   YDATA( 3,1,4), XDATA( 3,1,4), XDATA( 3,2,4)
     &   /0.81E+0_R8, 0.0E+0_R8, 2.0E+0_R8/
      DATA
     &   YDATA( 4,1,4), XDATA( 4,1,4), XDATA( 4,2,4)
     &   /0.58E+0_R8, 0.0E+0_R8, 3.0E+0_R8/
      DATA
     &   YDATA( 5,1,4), XDATA( 5,1,4), XDATA( 5,2,4)
     &   /5.90E+0_R8, 1.0E+0_R8, 0.0E+0_R8/
      DATA
     &   YDATA( 6,1,4), XDATA( 6,1,4), XDATA( 6,2,4)
     &   /4.74E+0_R8, 1.0E+0_R8, 1.0E+0_R8/
      DATA
     &   YDATA( 7,1,4), XDATA( 7,1,4), XDATA( 7,2,4)
     &   /4.18E+0_R8, 1.0E+0_R8, 2.0E+0_R8/
      DATA
     &   YDATA( 8,1,4), XDATA( 8,1,4), XDATA( 8,2,4)
     &   /4.05E+0_R8, 1.0E+0_R8, 2.0E+0_R8/
      DATA
     &   YDATA( 9,1,4), XDATA( 9,1,4), XDATA( 9,2,4)
     &   /9.03E+0_R8, 2.0E+0_R8, 0.0E+0_R8/
      DATA
     &   YDATA(10,1,4), XDATA(10,1,4), XDATA(10,2,4)
     &   /7.85E+0_R8, 2.0E+0_R8, 1.0E+0_R8/
      DATA
     &   YDATA(11,1,4), XDATA(11,1,4), XDATA(11,2,4)
     &   /7.22E+0_R8, 2.0E+0_R8, 2.0E+0_R8/
      DATA
     &   YDATA(12,1,4), XDATA(12,1,4), XDATA(12,2,4)
     &   /8.50E+0_R8, 2.5E+0_R8, 2.0E+0_R8/
      DATA
     &   YDATA(13,1,4), XDATA(13,1,4), XDATA(13,2,4)
     &   /9.81E+0_R8, 2.9E+0_R8, 1.8E+0_R8/

      DATA
     &   TDATA(5)
     &   /' DRAPER AND SMITH, 1981, EXERCISE I, PAGE 521-522'/
      DATA
     &   NDATA(5), MDATA(5), NPDATA(5), NQDATA(5)
     &   /8, 2, 2, 1/
      DATA
     &   (BDATA(K,5),K=1,2)
     &   /0.01155E+0_R8, 5000.0E+0_R8/
      DATA
     &   YDATA(1,1,5), XDATA(1,1,5), XDATA(1,2,5)
     &   /0.912E+0_R8,  109.0E+0_R8, 600.0E+0_R8/
      DATA
     &   YDATA(2,1,5), XDATA(2,1,5), XDATA(2,2,5)
     &   /0.382E+0_R8,   65.0E+0_R8, 640.0E+0_R8/
      DATA
     &   YDATA(3,1,5), XDATA(3,1,5), XDATA(3,2,5)
     &   /0.397E+0_R8, 1180.0E+0_R8, 600.0E+0_R8/
      DATA
     &   YDATA(4,1,5), XDATA(4,1,5), XDATA(4,2,5)
     &   /0.376E+0_R8,   66.0E+0_R8, 640.0E+0_R8/
      DATA
     &   YDATA(5,1,5), XDATA(5,1,5), XDATA(5,2,5)
     &   /0.342E+0_R8, 1270.0E+0_R8, 600.0E+0_R8/
      DATA
     &   YDATA(6,1,5), XDATA(6,1,5), XDATA(6,2,5)
     &   /0.358E+0_R8,   69.0E+0_R8, 640.0E+0_R8/
      DATA
     &   YDATA(7,1,5), XDATA(7,1,5), XDATA(7,2,5)
     &   /0.348E+0_R8, 1230.0E+0_R8, 600.0E+0_R8/
      DATA
     &   YDATA(8,1,5), XDATA(8,1,5), XDATA(8,2,5)
     &   /0.376E+0_R8,   68.0E+0_R8, 640.0E+0_R8/

      DATA
     &   TDATA(6)
     &   /' POWELL AND MACDONALD, 1972, TABLES 7 AND 8, PAGES 153-154'/
      DATA
     &   NDATA(6), MDATA(6), NPDATA(6), NQDATA(6)
     &   /14, 1, 3, 1/
      DATA
     &   (BDATA(K,6),K=1,3)
     &   /25.0E+0_R8, 30.0E+0_R8, 6.0E+0_R8/
      DATA
     &   YDATA( 1,1,6), XDATA( 1,1,6)
     &   /26.38E+0_R8,  1.0E+0_R8/
      DATA
     &   YDATA( 2,1,6), XDATA( 2,1,6)
     &   /25.79E+0_R8,  2.0E+0_R8/
      DATA
     &   YDATA( 3,1,6), XDATA( 3,1,6)
     &   /25.29E+0_R8,  3.0E+0_R8/
      DATA
     &   YDATA( 4,1,6), XDATA( 4,1,6)
     &   /24.86E+0_R8,  4.0E+0_R8/
      DATA
     &   YDATA( 5,1,6), XDATA( 5,1,6)
     &   /24.46E+0_R8,  5.0E+0_R8/
      DATA
     &   YDATA( 6,1,6), XDATA( 6,1,6)
     &   /24.10E+0_R8,  6.0E+0_R8/
      DATA
     &   YDATA( 7,1,6), XDATA( 7,1,6)
     &   /23.78E+0_R8,  7.0E+0_R8/
      DATA
     &   YDATA( 8,1,6), XDATA( 8,1,6)
     &   /23.50E+0_R8,  8.0E+0_R8/
      DATA
     &   YDATA( 9,1,6), XDATA( 9,1,6)
     &   /23.24E+0_R8,  9.0E+0_R8/
      DATA
     &   YDATA(10,1,6), XDATA(10,1,6)
     &   /23.00E+0_R8, 10.0E+0_R8/
      DATA
     &   YDATA(11,1,6), XDATA(11,1,6)
     &   /22.78E+0_R8, 11.0E+0_R8/
      DATA
     &   YDATA(12,1,6), XDATA(12,1,6)
     &   /22.58E+0_R8, 12.0E+0_R8/
      DATA
     &   YDATA(13,1,6), XDATA(13,1,6)
     &   /22.39E+0_R8, 13.0E+0_R8/
      DATA
     &   YDATA(14,1,6), XDATA(14,1,6)
     &   /22.22E+0_R8, 14.0E+0_R8/

      DATA
     &   TDATA(7)
     &   /' FULLER, 1987, TABLE 3.2.10, PAGES 244-245'/
      DATA
     &   NDATA(7), MDATA(7), NPDATA(7), NQDATA(7)
     &   /20, 2, 5, 1/
      DATA
     &   (BDATA(K,7),K=1,5)
     &   /-1.0E+0_R8, -3.0E+0_R8, 0.09E+0_R8, 0.02E+0_R8, 0.08E+0_R8/
      DATA
     &   YDATA( 1,1,7), XDATA( 1,1,7), XDATA( 1,2,7)
     &   /0.0E+0_R8,  0.50E+0_R8, -0.12E+0_R8/
      DATA
     &   YDATA( 2,1,7), XDATA( 2,1,7), XDATA( 2,2,7)
     &   /0.0E+0_R8,  1.20E+0_R8, -0.60E+0_R8/
      DATA
     &   YDATA( 3,1,7), XDATA( 3,1,7), XDATA( 3,2,7)
     &   /0.0E+0_R8,  1.60E+0_R8, -1.00E+0_R8/
      DATA
     &   YDATA( 4,1,7), XDATA( 4,1,7), XDATA( 4,2,7)
     &   /0.0E+0_R8,  1.86E+0_R8, -1.40E+0_R8/
      DATA
     &   YDATA( 5,1,7), XDATA( 5,1,7), XDATA( 5,2,7)
     &   /0.0E+0_R8,  2.12E+0_R8, -2.54E+0_R8/
      DATA
     &   YDATA( 6,1,7), XDATA( 6,1,7), XDATA( 6,2,7)
     &   /0.0E+0_R8,  2.36E+0_R8, -3.36E+0_R8/
      DATA
     &   YDATA( 7,1,7), XDATA( 7,1,7), XDATA( 7,2,7)
     &   /0.0E+0_R8,  2.44E+0_R8, -4.00E+0_R8/
      DATA
     &   YDATA( 8,1,7), XDATA( 8,1,7), XDATA( 8,2,7)
     &   /0.0E+0_R8,  2.36E+0_R8, -4.75E+0_R8/
      DATA
     &   YDATA( 9,1,7), XDATA( 9,1,7), XDATA( 9,2,7)
     &   /0.0E+0_R8,  2.06E+0_R8, -5.25E+0_R8/
      DATA
     &   YDATA(10,1,7), XDATA(10,1,7), XDATA(10,2,7)
     &   /0.0E+0_R8,  1.74E+0_R8, -5.64E+0_R8/
      DATA
     &   YDATA(11,1,7), XDATA(11,1,7), XDATA(11,2,7)
     &   /0.0E+0_R8,  1.34E+0_R8, -5.97E+0_R8/
      DATA
     &   YDATA(12,1,7), XDATA(12,1,7), XDATA(12,2,7)
     &   /0.0E+0_R8,  0.90E+0_R8, -6.32E+0_R8/
      DATA
     &   YDATA(13,1,7), XDATA(13,1,7), XDATA(13,2,7)
     &   /0.0E+0_R8, -0.28E+0_R8, -6.44E+0_R8/
      DATA
     &   YDATA(14,1,7), XDATA(14,1,7), XDATA(14,2,7)
     &   /0.0E+0_R8, -0.78E+0_R8, -6.44E+0_R8/
      DATA
     &   YDATA(15,1,7), XDATA(15,1,7), XDATA(15,2,7)
     &   /0.0E+0_R8, -1.36E+0_R8, -6.41E+0_R8/
      DATA
     &   YDATA(16,1,7), XDATA(16,1,7), XDATA(16,2,7)
     &   /0.0E+0_R8, -1.90E+0_R8, -6.25E+0_R8/
      DATA
     &   YDATA(17,1,7), XDATA(17,1,7), XDATA(17,2,7)
     &   /0.0E+0_R8, -2.50E+0_R8, -5.88E+0_R8/
      DATA
     &   YDATA(18,1,7), XDATA(18,1,7), XDATA(18,2,7)
     &   /0.0E+0_R8, -2.88E+0_R8, -5.50E+0_R8/
      DATA
     &   YDATA(19,1,7), XDATA(19,1,7), XDATA(19,2,7)
     &   /0.0E+0_R8, -3.18E+0_R8, -5.24E+0_R8/
      DATA
     &   YDATA(20,1,7), XDATA(20,1,7), XDATA(20,2,7)
     &   /0.0E+0_R8, -3.44E+0_R8, -4.86E+0_R8/

      DATA
     &   TDATA(8)
     &   /' BATES AND WATTS, 1988, TABLE A1.13, PAGES 280-281'/
      DATA
     &   NDATA(8), MDATA(8), NPDATA(8), NQDATA(8)
     &   /23, 1, 5, 2/
      DATA
     &   (BDATA(K,8),K=1,5)
     &   /4.0E+0_R8, 2.0E+0_R8, 7.0E+0_R8, 0.40E+0_R8, 0.50E+0_R8/
      DATA
     &   YDATA( 1,1,8), YDATA( 1,2,8), XDATA( 1,1,8)
     &   /4.220E+0_R8, 0.136E+0_R8,     30.0E+0_R8/
      DATA
     &   YDATA( 2,1,8), YDATA( 2,2,8), XDATA( 2,1,8)
     &   /4.167E+0_R8, 0.167E+0_R8,     50.0E+0_R8/
      DATA
     &   YDATA( 3,1,8), YDATA( 3,2,8), XDATA( 3,1,8)
     &   /4.132E+0_R8, 0.188E+0_R8,     70.0E+0_R8/
      DATA
     &   YDATA( 4,1,8), YDATA( 4,2,8), XDATA( 4,1,8)
     &   /4.038E+0_R8, 0.212E+0_R8,    100.0E+0_R8/
      DATA
     &   YDATA( 5,1,8), YDATA( 5,2,8), XDATA( 5,1,8)
     &   /4.019E+0_R8, 0.236E+0_R8,    150.0E+0_R8/
      DATA
     &   YDATA( 6,1,8), YDATA( 6,2,8), XDATA( 6,1,8)
     &   /3.956E+0_R8, 0.257E+0_R8,    200.0E+0_R8/
      DATA
     &   YDATA( 7,1,8), YDATA( 7,2,8), XDATA( 7,1,8)
     &   /3.884E+0_R8, 0.276E+0_R8,    300.0E+0_R8/
      DATA
     &   YDATA( 8,1,8), YDATA( 8,2,8), XDATA( 8,1,8)
     &   /3.784E+0_R8, 0.297E+0_R8,    500.0E+0_R8/
      DATA
     &   YDATA( 9,1,8), YDATA( 9,2,8), XDATA( 9,1,8)
     &   /3.713E+0_R8, 0.309E+0_R8,    700.0E+0_R8/
      DATA
     &   YDATA(10,1,8), YDATA(10,2,8), XDATA(10,1,8)
     &   /3.633E+0_R8, 0.311E+0_R8,   1000.0E+0_R8/
      DATA
     &   YDATA(11,1,8), YDATA(11,2,8), XDATA(11,1,8)
     &   /3.540E+0_R8, 0.314E+0_R8,   1500.0E+0_R8/
      DATA
     &   YDATA(12,1,8), YDATA(12,2,8), XDATA(12,1,8)
     &   /3.433E+0_R8, 0.311E+0_R8,   2000.0E+0_R8/
      DATA
     &   YDATA(13,1,8), YDATA(13,2,8), XDATA(13,1,8)
     &   /3.358E+0_R8, 0.305E+0_R8,   3000.0E+0_R8/
      DATA
     &   YDATA(14,1,8), YDATA(14,2,8), XDATA(14,1,8)
     &   /3.258E+0_R8, 0.289E+0_R8,   5000.0E+0_R8/
      DATA
     &   YDATA(15,1,8), YDATA(15,2,8), XDATA(15,1,8)
     &   /3.193E+0_R8, 0.277E+0_R8,   7000.0E+0_R8/
      DATA
     &   YDATA(16,1,8), YDATA(16,2,8), XDATA(16,1,8)
     &   /3.128E+0_R8, 0.255E+0_R8,  10000.0E+0_R8/
      DATA
     &   YDATA(17,1,8), YDATA(17,2,8), XDATA(17,1,8)
     &   /3.059E+0_R8, 0.240E+0_R8,  15000.0E+0_R8/
      DATA
     &   YDATA(18,1,8), YDATA(18,2,8), XDATA(18,1,8)
     &   /2.984E+0_R8, 0.218E+0_R8,  20000.0E+0_R8/
      DATA
     &   YDATA(19,1,8), YDATA(19,2,8), XDATA(19,1,8)
     &   /2.934E+0_R8, 0.202E+0_R8,  30000.0E+0_R8/
      DATA
     &   YDATA(20,1,8), YDATA(20,2,8), XDATA(20,1,8)
     &   /2.876E+0_R8, 0.182E+0_R8,  50000.0E+0_R8/
      DATA
     &   YDATA(21,1,8), YDATA(21,2,8), XDATA(21,1,8)
     &   /2.838E+0_R8, 0.168E+0_R8,  70000.0E+0_R8/
      DATA
     &   YDATA(22,1,8), YDATA(22,2,8), XDATA(22,1,8)
     &   /2.798E+0_R8, 0.153E+0_R8, 100000.0E+0_R8/
      DATA
     &   YDATA(23,1,8), YDATA(23,2,8), XDATA(23,1,8)
     &   /2.759E+0_R8, 0.139E+0_R8, 150000.0E+0_R8/

      DATA
     &   TDATA(9)
     &   /' ZWOLAK, WATSON, AND TYSON, 2004.'/
      DATA
     &   NDATA(9), MDATA(9), NPDATA(9), NQDATA(9)
     &   /4, 1, 2, 1/
      DATA
     &   (BDATA(K,9),K=1,2)
     &   /200.0_R8, 5.0_R8/
      DATA
     &   YDATA( 1,1,9), XDATA( 1,1,9)
     &   /2.718281828459045_R8, 1.0_R8/
      DATA
     &   YDATA( 2,1,9), XDATA( 2,1,9)
     &   /7.389056098930650_R8, 2.0_R8/
      DATA
     &   YDATA( 3,1,9), XDATA( 3,1,9)
     &   /148.4131591025766_R8, 5.0_R8/
      DATA
     &   YDATA( 4,1,9), XDATA( 4,1,9)
     &   /403.4287934927353_R8, 6.0_R8/

      DATA
     &   TDATA(10)
     &   /' ZWOLAK, WATSON, AND TYSON, 2005.'/
      DATA
     &   NDATA(10), MDATA(10), NPDATA(10), NQDATA(10)
     &   /4, 1, 2, 1/
      DATA
     &   (BDATA(K,10),K=1,2)
     &   /200.0_R8, 5.0_R8/
      DATA
     &   YDATA( 1,1,10), XDATA( 1,1,10)
     &   /2.718281828459045_R8, 1.0_R8/
      DATA
     &   YDATA( 2,1,10), XDATA( 2,1,10)
     &   /7.389056098930650_R8, 2.0_R8/
      DATA
     &   YDATA( 3,1,10), XDATA( 3,1,10)
     &   /148.4131591025766_R8, 5.0_R8/
      DATA
     &   YDATA( 4,1,10), XDATA( 4,1,10)
     &   /403.4287934927353_R8, 6.0_R8/

C...Variable definitions (alphabetically)
C   BDATA:   The function parameter for each data set.
C   BETA:    The function parameters.
C   I:       An indexing variable.
C   J:       An indexing variable.
C   L:       An indexing variable.
C   LDX:     The leading dimension of array X.
C   M:       The number of columns of data in the explanatory variable.
C   MDATA:   The number of columns of data in the explanatory variable
C            in each data set.
C   N:       The number of observations.
C   NDATA:   The number of observations per data set.
C   NP:      The number of function parameters.
C   NPDATA:  The number of function parameters in each data set.
C   NQDATA:  The number of responses per observation in each data set.
C   SETNO:   The number of the data set being analyzed.
C   TDATA:   The reference for the each of the data sets.
C   TITLE:   The reference for the data set being analyzed.
C   X:       The explanatory variables.
C   XDATA:   The explanatory variables for each data set.
C   Y:       The response variable.
C   YDATA:   The response variables for each data set.


C***First executable statement  DODRXD


      TITLE = TDATA(SETNO)

      N = NDATA(SETNO)
      M = MDATA(SETNO)
      NP = NPDATA(SETNO)
      NQ = NQDATA(SETNO)

      DO 20 L=1,NQ
         DO 10 I=1,N
            Y(I,L) = YDATA(I,L,SETNO)
   10    CONTINUE
   20 CONTINUE

      DO 40 J=1,M
         DO 30 I=1,N
            X(I,J) = XDATA(I,J,SETNO)
   30    CONTINUE
   40 CONTINUE

      DO 50 K=1,NP
         BETA(K) = BDATA(K,SETNO)
   50 CONTINUE

      RETURN

      END
*DODRXF
      SUBROUTINE DODRXF
     &   (N,M,NP,NQ, 
     &    LDN,LDM,LDNP,
     &    BETA,XPLUSD,
     &    IFIXB,IFIXX,LDIFX,
     &    IDEVAL,F,FJACB,FJACD,
     &    ISTOP)
C***Begin Prologue  DODRXF
C***Refer to  ODR
C***Routines Called  (NONE)
C***Date Written   860529   (YYMMDD)
C***Revision Date  920619   (YYMMDD)
C***Purpose  Compute jacobian matricies for ODRPACK95 exerciser
C***End Prologue  DODRXF

C...Used modules
      USE REAL_PRECISION

C...Scalar arguments
      INTEGER
     &   IDEVAL,ISTOP,LDIFX,LDM,LDN,LDNP,M,N,NP,NQ

C...Array arguments
      REAL (KIND=R8)
     &   BETA(NP),F(LDN,NQ),FJACB(LDN,LDNP,NQ),FJACD(LDN,LDM,NQ),
     &   XPLUSD(LDN,M)
      INTEGER
     &   IFIXB(NP),IFIXX(LDIFX,M)

C...Scalar parameters
      INTEGER, PARAMETER :: MAXNP=10

C...Scalars in common
      INTEGER
     &   SETNO

C...Arrays in common
      REAL (KIND=R8)
     &   LOWER(MAXNP),UPPER(MAXNP)

C...Local scalars
      REAL (KIND=R8)
     &   CTHETA,FAC1,FAC2,FAC3,FAC4,FREQ,
     &   OMEGA,ONE,PHI,PI,R,STHETA,THETA,ZERO
      INTEGER
     &   I,J,K

C...Intrinsic functions
      INTRINSIC
     &   ATAN2,COS,EXP,SIN,SQRT

C...Common blocks
      COMMON /SETID/SETNO
      COMMON /BOUNDS/ LOWER,UPPER

C...Data statements
      DATA
     &   ZERO,ONE
     &   /0.0E0_R8,1.0E0_R8/

C...Variable definitions (alphabetically)
C   BETA:    Current values of parameters
C   F:       Predicted function values
C   FAC1:    A factors or terms used in computing the jacobians.
C   FAC2:    A factors or terms used in computing the jacobians.
C   FAC3:    A factors or terms used in computing the jacobians.
C   FAC4:    A factors or terms used in computing the jacobians.
C   FJACB:   Jacobian with respect to BETA
C   FJACD:   Jacobian with respect to errors DELTA
C   IDEVAL:  Indicator for selecting computation to be performed
C   IFIXB:   Indicators for "fixing" parameters (BETA)
C   IFIXX:   Indicators for "fixing" explanatory variable (X)
C   LDIFX:   Leading dimension of array IFIXX
C   ISTOP:   Stopping condition, where
C                     0 means current BETA and X+DELTA were
C                       acceptable and values were computed successfully
C                     1 means current BETA and X+DELTA are
C                       not acceptable;  ODRPACK95 should select
C                       values closer to most recently used values
C                    -1 means current BETA and X+DELTA are
C                       not acceptable;  ODRPACK95 should stop
C   LDN:     Leading dimension declarator equal or exceeding N
C   LDM:     Leading dimension declarator equal or exceeding M
C   LDNP:    Leading dimension declarator equal or exceeding NP
C   M:       The number of columns of data in the explanatory variable.
C   N:       The number of observations.
C   NP:      The number of function parameters.
C   NQ:      The number of responses per observation.
C   ONE:     The value 1.0E0_R8.
C   SETNO:   The number of the data set being analyzed.
C   XPLUSD:  Current value of explanatory variable, i.e., X + DELTA
C   ZERO:    The value 0.0E0_R8.


C***First executable statement  DODRXF


C  Do something with FJACD, FJACB, IFIXB and IFIXX to avoid warnings that they 
C  are not being used.  This is simply not to worry users that the example 
C  program is failing.
      IF (IFIXB(1) .GT. 0 .AND. IFIXX(1,1) .GT. 0 
     &    .AND. FJACB(1,1,1) .GT. 0 .AND. FJACD(1,1,1) .GT. 0 ) THEN
C        Do nothing.
      END IF


C  Check for BETA outside bounds.  Return with error if BETA outside bounds.

      IF (ANY(LOWER(1:NP).GT.BETA(1:NP))) THEN
         ISTOP = -1
         RETURN
      END IF

      IF (ANY(UPPER(1:NP).LT.BETA(1:NP))) THEN
         ISTOP = -2
         RETURN
      END IF


      IF (SETNO.EQ.1) THEN

C  Setno. 1:  Boggs, Byrd and Schnabel, 1985, example 1

         IF (BETA(1).LE.1.01E0_R8) THEN
            ISTOP = 0

            IF (MOD(IDEVAL,10).NE.0) THEN
               DO 100 I=1,N
                  F(I,1) = BETA(1)/(XPLUSD(I,1)-BETA(2))
  100          CONTINUE
            END IF

            IF (MOD(IDEVAL/10,10).NE.0) THEN
               DO 110 I=1,N
                  FJACB(I,1,1) = ONE/(XPLUSD(I,1)-BETA(2))
                  FJACB(I,2,1) = BETA(1)*(XPLUSD(I,1)-BETA(2))**(-2)
  110          CONTINUE
            END IF

            IF (MOD(IDEVAL/100,10).NE.0) THEN
               DO 120 I=1,N
                  FJACD(I,1,1) = -BETA(1)*(XPLUSD(I,1)-BETA(2))**(-2)
  120          CONTINUE
            END IF

         ELSE
            ISTOP = 1
         END IF

      ELSE IF (SETNO.EQ.2) THEN

C  Setno. 2:  Boggs, Byrd and Schnabel, 1985, example 2

         ISTOP = 0

         DO 200 I=1,N
            FAC1 = (BETA(2)*XPLUSD(I,1)+BETA(3)*XPLUSD(I,2)-ONE)

            IF (MOD(IDEVAL,10).NE.0) THEN
               F(I,1) = BETA(1)/FAC1
            END IF

            IF (MOD(IDEVAL/10,10).NE.0) THEN
               FJACB(I,1,1) = ONE/FAC1
               FJACB(I,2,1) = -BETA(1)*(FAC1**(-2))*XPLUSD(I,1)
               FJACB(I,3,1) = -BETA(1)*(FAC1**(-2))*XPLUSD(I,2)
            END IF

            IF (MOD(IDEVAL/100,10).NE.0) THEN
               FJACD(I,1,1) = -BETA(1)*(FAC1**(-2))*BETA(2)
               FJACD(I,2,1) = -BETA(1)*(FAC1**(-2))*BETA(3)
            END IF
  200    CONTINUE

      ELSE IF (SETNO.EQ.3) THEN

C  Setno. 3:  Boggs, Byrd and Schnabel, 1985, example 3

         ISTOP = 0

         IF (MOD(IDEVAL,10).NE.0) THEN
            DO 310 I=1,N
               F(I,1) = ZERO
               DO 300 J=1,4
                  F(I,1) = F(I,1) + BETA(J)/(XPLUSD(I,1)+BETA(J+5))
  300          CONTINUE
               F(I,1) = F(I,1) + BETA(5)
  310       CONTINUE
         END IF

         IF (MOD(IDEVAL/10,10).NE.0) THEN
            DO 330 I=1,N
               FJACB(I,5,1) = ONE
               DO 320 K=1,4
                  FJACB(I,K,1) = ONE/(XPLUSD(I,1)+BETA(K+5))
                  FJACB(I,K+5,1) = -BETA(K)*
     &                              (XPLUSD(I,1)+BETA(K+5))**(-2)
  320          CONTINUE
  330       CONTINUE
         END IF

         IF (MOD(IDEVAL/100,10).NE.0) THEN
            DO 350 I=1,N
               FJACD(I,1,1) = ZERO
               DO 340 K=4,1,-1
                  FJACD(I,1,1) = FJACD(I,1,1) -
     &                         BETA(K)*(XPLUSD(I,1)+BETA(K+5))**(-2)
  340          CONTINUE
  350       CONTINUE
         END IF

      ELSE IF (SETNO.EQ.4) THEN

C  Setno. 4:  Himmelblau, 1970, example 6.2-4, page 188

         ISTOP = 0

         IF (MOD(IDEVAL,10).NE.0) THEN
            DO 400 I = 1, N
               F(I,1) = BETA(1)*XPLUSD(I,1) +
     &                BETA(2)*EXP(BETA(3)*XPLUSD(I,2))
  400       CONTINUE
         END IF

         IF (MOD(IDEVAL/10,10).NE.0) THEN
            DO 410 I=1,N
               FJACB(I,1,1) = XPLUSD(I,1)
               FJACB(I,2,1) = EXP(BETA(3)*XPLUSD(I,2))
               FJACB(I,3,1) = BETA(2)*
     &                        EXP(BETA(3)*XPLUSD(I,2))*XPLUSD(I,2)
  410       CONTINUE
         END IF

         IF (MOD(IDEVAL/100,10).NE.0) THEN
            DO 420 I=1,N
               FJACD(I,1,1) = BETA(1)
               FJACD(I,2,1) = BETA(2)*EXP(BETA(3)*XPLUSD(I,2))*BETA(3)
  420       CONTINUE
         END IF

      ELSE IF (SETNO.EQ.5) THEN

C  Setno. 5:  Draper and Smith, 1981, exercise i, page 521-522

         ISTOP = 0

         IF (MOD(IDEVAL,10).NE.0) THEN
            DO 500 I=1,N
               F(I,1) = EXP(-BETA(1)*XPLUSD(I,1)*
     &                EXP(-BETA(2)*(ONE/XPLUSD(I,2) - ONE/620.0E0_R8)))
  500       CONTINUE
         END IF

         IF (MOD(IDEVAL/10,10).NE.0) THEN
            DO 510 I=1,N
               FAC1 = ONE/XPLUSD(I,2) - ONE/620.0E0_R8
               FAC2 = EXP(-BETA(2)*FAC1)
               FAC3 = BETA(1)*XPLUSD(I,1)
               FAC4 = EXP(-FAC3*FAC2)

               FJACB(I,1,1) = -FAC4*XPLUSD(I,1)*FAC2
               FJACB(I,2,1) = FAC4*FAC3*FAC2*FAC1

               IF (MOD(IDEVAL/100,10).NE.0) THEN
                  FJACD(I,1,1) = -FAC4*BETA(1)*FAC2
                  FJACD(I,2,1) = -FAC4*FAC3*FAC2*
     &                           BETA(2)/XPLUSD(I,2)**2
               END IF
  510       CONTINUE
         END IF

      ELSE IF (SETNO.EQ.6) THEN

C  Setno. 6:  Powell and Macdonald, 1972, tables 7 and 8, page 153-154
C             N.B.  this derivative is intentionally coded incorrectly

         ISTOP = 0

         IF (MOD(IDEVAL,10).NE.0) THEN
            DO 600 I=1,N
               F(I,1) = BETA(1)*
     &                (ONE+BETA(3)*XPLUSD(I,1)/BETA(2))**(-ONE/BETA(3))
  600       CONTINUE
         END IF

         IF (MOD(IDEVAL/10,10).NE.0) THEN
            DO 610 I=1,N
               FJACB(I,1,1) = ZERO
               FJACB(I,2,1) = ZERO
               FJACB(I,3,1) = ZERO

               IF (MOD(IDEVAL/100,10).NE.0) THEN
                  FJACD(I,1,1) = XPLUSD(I,1)
               END IF
  610       CONTINUE
         END IF

      ELSE IF (SETNO.EQ.7) THEN

C  Setno. 7:  Fuller, 1987, table 3.2.10, pages 244-245
C             N.B.  this derivative is intentionally coded incorrectly

         ISTOP = 0

         IF (MOD(IDEVAL,10).NE.0) THEN
            DO 700 I=1,N
               F(I,1) = BETA(3)*(XPLUSD(I,1)-BETA(1))**2 +
     &                  2*BETA(4)*(XPLUSD(I,1)-BETA(1))*
     &                            (XPLUSD(I,2)-BETA(2)) +
     &                  BETA(5)*(XPLUSD(I,2)-BETA(2))**2 - 1.0E0_R8
  700       CONTINUE
         END IF

         IF (MOD(IDEVAL/10,10).NE.0) THEN
            DO 710 I=1,N
               FJACB(I,1,1) = ZERO
               FJACB(I,2,1) = ZERO
               FJACB(I,3,1) = ZERO
               FJACB(I,4,1) = ZERO
               FJACB(I,5,1) = ZERO

               IF (MOD(IDEVAL/100,10).NE.0) THEN
                  FJACD(I,1,1) = ZERO
                  FJACD(I,2,1) = ZERO
               END IF
  710       CONTINUE
         END IF

      ELSE IF (SETNO.EQ.8) THEN

C  Setno. 8:  Bates and Watts, 1988, table A1.13, pages 280-281
C             N.B.  This derivative is intentionally coded incorrectly

         DO 800 I=1,N
            IF (XPLUSD(I,1).LT.0.0E0_R8) THEN
               ISTOP = 1
               RETURN
            END IF
  800    CONTINUE
         ISTOP = 0
 
         IF (MOD(IDEVAL,10).NE.0) THEN
            PI = 3.141592653589793238462643383279E0_R8
            THETA = PI*BETA(4)*0.5E0_R8
            CTHETA = COS(THETA)
            STHETA = SIN(THETA)
            DO 810 I=1,N
               FREQ   = XPLUSD(I,1)
               OMEGA  = (2.0E0_R8*PI*FREQ*EXP(-BETA(3)))**BETA(4)
               PHI    = ATAN2((OMEGA*STHETA),(1+OMEGA*CTHETA))
               R      = (BETA(1)-BETA(2)) *
     &                  SQRT((1+OMEGA*CTHETA)**2+(OMEGA*STHETA)**2)**
     &                  (-BETA(5))
               F(I,1) = BETA(2) + R*COS(BETA(5)*PHI)
               F(I,2) =           R*SIN(BETA(5)*PHI)
  810       CONTINUE
         END IF

         IF (MOD(IDEVAL/10,10).NE.0) THEN
            DO 820 I=1,N
               FJACB(I,1,1) = ZERO
               FJACB(I,2,1) = ZERO
               FJACB(I,3,1) = ZERO
               FJACB(I,4,1) = ZERO
               FJACB(I,5,1) = ZERO
   
               FJACB(I,1,2) = ZERO
               FJACB(I,2,2) = ZERO
               FJACB(I,3,2) = ZERO
               FJACB(I,4,2) = ZERO
               FJACB(I,5,2) = ZERO
   
               IF (MOD(IDEVAL/100,10).NE.0) THEN
                  FJACD(I,1,1) = ZERO
                  FJACD(I,1,2) = ZERO
               END IF
  820       CONTINUE
         END IF

      ELSE IF (SETNO.EQ.9) THEN

C  Setno. 9:  Zwolak, Watson, and Tyson, 2004.

          ISTOP = 0

          IF (MOD(IDEVAL,10).NE.0) THEN
             DO I=1,N
                F(I,1) = BETA(1)*EXP(BETA(2)*XPLUSD(I,1))
             END DO
          END IF

          IF (MOD(IDEVAL/10,10).NE.0) THEN
             DO I=1,N
                FJACB(I,1,1) = EXP(BETA(2)*XPLUSD(I,1))
                FJACB(I,2,1) = BETA(1)*XPLUSD(I,1)*EXP(BETA(2)*
     &             XPLUSD(I,1))
             END DO
          END IF

          IF (MOD(IDEVAL/100,10).NE.0) THEN
             DO I=1,N
                FJACD(I,1,1) = BETA(1)*BETA(2)*EXP(BETA(2)*XPLUSD(I,1))
             END DO
          END IF

      ELSE IF (SETNO.EQ.10) THEN

C  Setno. 10:  Zwolak, Watson, and Tyson, 2005.

          ISTOP = 0

          IF (MOD(IDEVAL,10).NE.0) THEN
             DO I=1,N
                F(I,1) = BETA(1)/2.0_R8*EXP(BETA(2)*XPLUSD(I,1))
             END DO
          END IF

          IF (MOD(IDEVAL/10,10).NE.0) THEN
             DO I=1,N
                FJACB(I,1,1) = EXP(BETA(2)*XPLUSD(I,1))/2.0_R8
                FJACB(I,2,1) = BETA(1)/2.0_R8*XPLUSD(I,1)*EXP(BETA(2)*
     &             XPLUSD(I,1))
             END DO
          END IF

          IF (MOD(IDEVAL/100,10).NE.0) THEN
             DO I=1,N
                FJACD(I,1,1) = BETA(1)/2.0_R8*BETA(2)*EXP(BETA(2)*
     &             XPLUSD(I,1))
             END DO
          END IF

      END IF

      RETURN

      END
*DODRXW
      SUBROUTINE DODRXW
     &   (MAXN,MAXM,MAXNP,MAXNQ,LDWE,LD2WE,ISODR,LIWMIN,LWMIN)
C***Begin Prologue  DODRXW
C***Refer to  ODR
C***Routines Called  (NONE)
C***Date Written   890205   (YYMMDD)
C***Revision Date  920619   (YYMMDD)
C***Purpose  Compute minimum lengths for work vectors
C***Routines Called  NONE
C***End Prologue  DODRXW

C...Used modules
      USE REAL_PRECISION

C...Scalar arguments
      INTEGER
     &   LDWE,LD2WE,LIWMIN,LWMIN,MAXN,MAXM,MAXNP,MAXNQ
      LOGICAL
     &   ISODR

C...Variable definitions (alphabetically)
C   ISODR:   The variable designating whether the solution is by odr
C            (ISODR=TRUE) or by ols (ISODR=FALSE).
C   LDWE:    The leading dimension of array WE.
C   LD2WE:   The second dimension of array WE.
C   LIWMIN:  The minimum length of vector IWORK for a given problem.
C   LWMIN:   The minimum length of vector WORK for a given problem.
C   MAXM:    The number of columns in the explanatory variable.
C   MAXN:    The number of observations.
C   MAXNP:   The number of function parameters.
C   MAXNQ:   The number of responses per observation.


C***First executable statement  DODRXW


      LIWMIN = 20+MAXNP+MAXNQ*(MAXNP+MAXM) 
      IF (ISODR) THEN
         LWMIN = 18 + 11*MAXNP + MAXNP**2 + MAXM + MAXM**2 + 
     &           4*MAXN*MAXNQ + 6*MAXN*MAXM + 2*MAXN*MAXNQ*MAXNP +  
     &           2*MAXN*MAXNQ*MAXM + MAXNQ**2 + 
     &           5*MAXNQ + MAXNQ*(MAXNP+MAXM) + LDWE*LD2WE*MAXNQ 
      ELSE
         LWMIN = 18 + 11*MAXNP + MAXNP**2 + MAXM + MAXM**2 + 
     &           4*MAXN*MAXNQ + 2*MAXN*MAXM + 2*MAXN*MAXNQ*MAXNP + 
     &           5*MAXNQ + MAXNQ*(MAXNP+MAXM) + LDWE*LD2WE*MAXNQ
      END IF

      RETURN
      END
