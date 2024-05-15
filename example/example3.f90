PROGRAM example3
   USE odrpack
   use odrpack_kinds, only: wp
   implicit none

   ! ODRPACK95 Argument Definitions
   !      ==> FCN      Name of the user supplied function subroutine
   !      ==> N        Number of observations
   !      ==> M        Columns of data in the explanatory variable
   !      ==> NP       Number of parameters
   !      ==> NQ       Number of responses per observation
   !     <==> BETA     Function parameters
   !      ==> Y        Response variable
   !      ==> X        Explanatory variable
   !      ==> WE       "Epsilon" weights
   !      ==> WD       "Delta" weights
   !      ==> IFIXB    Indicators for "fixing" parameters (BETA)
   !      ==> IFIXX    Indicators for "fixing" explanatory variable (X)
   !      ==> JOB      Task to be performed
   !      ==> NDIGIT   Good digits in subroutine function results
   !      ==> TAUFAC   Trust region initialization factor
   !      ==> SSTOL    Sum of squares convergence criterion
   !      ==> PARTOL   Parameter convergence criterion
   !      ==> MAXIT    Maximum number of iterations
   !      ==> IPRINT   Print control
   !      ==> LUNERR   Logical unit for error reports
   !      ==> LUNRPT   Logical unit for computation reports
   !      ==> STPB     Step sizes for finite difference derivatives wrt BETA
   !      ==> STPD     Step sizes for finite difference derivatives wrt DELTA
   !      ==> SCLB     Scale values for parameters BETA
   !      ==> SCLD     Scale values for errors delta in explanatory variable
   !     <==> WORK     REAL (KIND=wp) work vector
   !     <==  IWORK    Integer work vector
   !     <==  INFO     Stopping condition

   ! Parameters specifying maximum problem sizes handled by this driver
   !     MAXN          Maximum number of observations
   !     MAXM          Maximum number of columns in explanatory variable
   !     MAXNP         Maximum number of function parameters
   !     MAXNQ         Maximum number of responses per observation

! Parameter declarations and specifications
   INTEGER :: LDIFX, LDSCLD, LDSTPD, LDWD, LDWE, LDX, LDY, LD2WD, LD2WE, &
              LIWORK, LWORK, MAXM, MAXN, MAXNP, MAXNQ
   PARAMETER(MAXM=5, MAXN=100, MAXNP=25, MAXNQ=5, &
             LDY=MAXN, LDX=MAXN, &
             LDWE=MAXN, LD2WE=MAXNQ, LDWD=MAXN, LD2WD=1, &
             LDIFX=MAXN, LDSCLD=1, LDSTPD=1, &
             LWORK=18 + 11*MAXNP + MAXNP**2 + MAXM + MAXM**2 + &
             4*MAXN*MAXNQ + 6*MAXN*MAXM + 2*MAXN*MAXNQ*MAXNP + &
             2*MAXN*MAXNQ*MAXM + MAXNQ**2 + &
             5*MAXNQ + MAXNQ*(MAXNP + MAXM) + LDWE*LD2WE*MAXNQ, &
             LIWORK=20 + MAXNP + MAXNQ*(MAXNP + MAXM))

! Variable declarations
   INTEGER :: I, INFO, IPRINT, J, JOB, L, LUNERR, LUNRPT, M, MAXIT, N, NDIGIT, NP, NQ
   INTEGER :: IFIXB(MAXNP), IFIXX(LDIFX, MAXM), IWORK(:)
   REAL(KIND=wp) :: PARTOL, SSTOL, TAUFAC
   REAL(KIND=wp) :: BETA(MAXNP), DELTA(:, :), &
                    SCLB(MAXNP), SCLD(LDSCLD, MAXM), &
                    STPB(MAXNP), STPD(LDSTPD, MAXM), &
                    WD(LDWD, LD2WD, MAXM), WE(LDWE, LD2WE, MAXNQ), &
                    WORK(:), X(LDX, MAXM), Y(LDY, MAXNQ)
   EXTERNAL :: FCN
   POINTER  :: DELTA, IWORK, WORK

! Specify default values for DODRC arguments
   WE(1, 1, 1) = -1.0E0_wp
   WD(1, 1, 1) = -1.0E0_wp
   IFIXB(1) = -1
   IFIXX(1, 1) = -1
   JOB = -1
   NDIGIT = -1
   TAUFAC = -1.0E0_wp
   SSTOL = -1.0E0_wp
   PARTOL = -1.0E0_wp
   MAXIT = -1
   IPRINT = -1
   LUNERR = -1
   LUNRPT = -1
   STPB(1) = -1.0E0_wp
   STPD(1, 1) = -1.0E0_wp
   SCLB(1) = -1.0E0_wp
   SCLD(1, 1) = -1.0E0_wp

! Set up ODRPACK95 report files
   LUNERR = 9
   LUNRPT = 9
   OPEN (UNIT=9, FILE='./example/report3.dat')

! Read problem data
   OPEN (UNIT=5, FILE='./example/data3.dat')
   READ (5, FMT=*) N, M, NP, NQ
   READ (5, FMT=*) (BETA(I), I=1, NP)
   DO I = 1, N
      READ (5, FMT=*) (X(I, J), J=1, M), (Y(I, L), L=1, NQ)
   END DO

! Allocate work arrays
   ALLOCATE (DELTA(N, M), IWORK(LIWORK), WORK(LWORK))

! Specify task as explicit orthogonal distance regression
!         With central difference derivatives
!         Covariance matrix constructed with recomputed derivatives
!         DELTA initialized by user
!         Not a restart
! And indicate long initial report
!              No iteration reports
!              Long final report
   JOB = 01010
   IPRINT = 2002

! Initialize DELTA, and specify first decade of frequencies as fixed
   DO I = 1, N
      IF (X(I, 1) .LT. 100.0E0_wp) THEN
         DELTA(I, 1) = 0.0E0_wp
         IFIXX(I, 1) = 0
      ELSE IF (X(I, 1) .LE. 150.0E0_wp) THEN
         DELTA(I, 1) = 0.0E0_wp
         IFIXX(I, 1) = 1
      ELSE IF (X(I, 1) .LE. 1000.0E0_wp) THEN
         DELTA(I, 1) = 25.0E0_wp
         IFIXX(I, 1) = 1
      ELSE IF (X(I, 1) .LE. 10000.0E0_wp) THEN
         DELTA(I, 1) = 560.0E0_wp
         IFIXX(I, 1) = 1
      ELSE IF (X(I, 1) .LE. 100000.0E0_wp) THEN
         DELTA(I, 1) = 9500.0E0_wp
         IFIXX(I, 1) = 1
      ELSE
         DELTA(I, 1) = 144000.0E0_wp
         IFIXX(I, 1) = 1
      END IF
   END DO

! Set weights
   DO I = 1, N
      IF (X(I, 1) .EQ. 100.0E0_wp .OR. X(I, 1) .EQ. 150.0E0_wp) THEN
         WE(I, 1, 1) = 0.0E0_wp
         WE(I, 1, 2) = 0.0E0_wp
         WE(I, 2, 1) = 0.0E0_wp
         WE(I, 2, 2) = 0.0E0_wp
      ELSE
         WE(I, 1, 1) = 559.6E0_wp
         WE(I, 1, 2) = -1634.0E0_wp
         WE(I, 2, 1) = -1634.0E0_wp
         WE(I, 2, 2) = 8397.0E0_wp
      END IF
      WD(I, 1, 1) = (1.0E-4_wp)/(X(I, 1)**2)
   END DO

! Compute solution
   CALL ODR(FCN=FCN, &
            N=N, M=M, NP=NP, NQ=NQ, &
            BETA=BETA, &
            Y=Y, X=X, &
            DELTA=DELTA, &
            WE=WE, WD=WD, &
            IFIXB=IFIXB, IFIXX=IFIXX, &
            JOB=JOB, NDIGIT=NDIGIT, TAUFAC=TAUFAC, &
            SSTOL=SSTOL, PARTOL=PARTOL, MAXIT=MAXIT, &
            IPRINT=IPRINT, LUNERR=LUNERR, LUNRPT=LUNRPT, &
            STPB=STPB, STPD=STPD, &
            SCLB=SCLB, SCLD=SCLD, &
            WORK=WORK, IWORK=IWORK, &
            INFO=INFO)

END PROGRAM example3

SUBROUTINE FCN(N, M, NP, NQ, &
               LDN, LDM, LDNP, &
               BETA, XPLUSD, &
               IFIXB, IFIXX, LDIFX, &
               IDEVAL, F, FJACB, FJACD, &
               ISTOP)

! Subroutine arguments
!      ==> N        Number of observations
!      ==> M        Number of columns in explanatory variable
!      ==> NP       Number of parameters
!      ==> NQ       Number of responses per observation
!      ==> LDN      Leading dimension declarator equal or exceeding N
!      ==> LDM      Leading dimension declarator equal or exceeding M
!      ==> LDNP     Leading dimension declarator equal or exceeding NP
!      ==> BETA     Current values of parameters
!      ==> XPLUSD   Current value of explanatory variable, i.e., X + DELTA
!      ==> IFIXB    Indicators for "fixing" parameters (BETA)
!      ==> IFIXX    Indicators for "fixing" explanatory variable (X)
!      ==> LDIFX    Leading dimension of array IFIXX
!      ==> IDEVAL   Indicator for selecting computation to be performed
!     <==  F        Predicted function values
!     <==  FJACB    Jacobian with respect to BETA
!     <==  FJACD    Jacobian with respect to errors DELTA
!     <==  ISTOP    Stopping condition, where
!                     0 means current BETA and X+DELTA were
!                       acceptable and values were computed successfully
!                     1 means current BETA and X+DELTA are
!                       not acceptable;  ODRPACK95 should select values
!                       closer to most recently used values if possible
!                    -1 means current BETA and X+DELTA are
!                       not acceptable; ODRPACK95 should stop

!Used modules
   use odrpack_kinds, only: wp, ZERO, ONE
   implicit none

! Subroutine arguments:
   INTEGER, INTENT(IN) :: IDEVAL, LDIFX, LDM, LDN, LDNP, M, N, NP, NQ
   INTEGER, INTENT(IN) :: IFIXB(NP), IFIXX(LDIFX, M)
   REAL(KIND=wp), INTENT(IN) :: BETA(NP), XPLUSD(LDN, M)
   REAL(KIND=wp), INTENT(OUT) :: F(LDN, NQ), FJACB(LDN, LDNP, NQ), FJACD(LDN, LDM, NQ)
   INTEGER, INTENT(OUT) :: ISTOP
! Local variables
   REAL(KIND=wp) :: FREQ, OMEGA, CTHETA, STHETA, THETA, PHI, R
   REAL(KIND=wp), PARAMETER :: PI = 4*ATAN(ONE)
   INTEGER :: I

! Do something with FJACD, FJACB, IFIXB and IFIXX to avoid warnings that they
! are not being used.  This is simply not to worry users that the example
! program is failing.
   IF (IFIXB(1) .GT. 0 .AND. IFIXX(1, 1) .GT. 0 &
       .AND. FJACB(1, 1, 1) .GT. 0 .AND. FJACD(1, 1, 1) .GT. 0) THEN
      ! Do nothing.
   END IF

! Check for unacceptable values for this problem
   DO I = 1, N
      IF (XPLUSD(I, 1) .LT. ZERO) THEN
         ISTOP = 1
         RETURN
      END IF
   END DO
   ISTOP = 0

   THETA = PI*BETA(4)*0.5E0_wp
   CTHETA = COS(THETA)
   STHETA = SIN(THETA)

! Compute predicted values
   IF (MOD(IDEVAL, 10) .GE. 1) THEN
      DO I = 1, N
         FREQ = XPLUSD(I, 1)
         OMEGA = (2.0E0_wp*PI*FREQ*EXP(-BETA(3)))**BETA(4)
         PHI = ATAN2((OMEGA*STHETA), (1 + OMEGA*CTHETA))
         R = (BETA(1) - BETA(2))*SQRT((1 + OMEGA*CTHETA)**2 + (OMEGA*STHETA)**2)**(-BETA(5))
         F(I, 1) = BETA(2) + R*COS(BETA(5)*PHI)
         F(I, 2) = R*SIN(BETA(5)*PHI)
      END DO
   END IF

END SUBROUTINE FCN
