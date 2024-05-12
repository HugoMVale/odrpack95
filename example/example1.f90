
PROGRAM EXAMPLE1
USE ODRPACK95
use odrpack95_kinds, only: wp
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

! Parameter Declarations and Specifications
INTEGER :: LDIFX,LDSCLD,LDSTPD,LDWD,LDWE,LDX,LDY,LD2WD,LD2WE,    &
           LIWORK,LWORK,MAXM,MAXN,MAXNP,MAXNQ
PARAMETER (MAXM=5,MAXN=25,MAXNP=5,MAXNQ=1,                       &
             LDY=MAXN,LDX=MAXN,                                  &  
             LDWE=1,LD2WE=1,LDWD=1,LD2WD=1,                      & 
             LDIFX=MAXN,LDSTPD=1,LDSCLD=1,                       &
             LWORK=18 + 11*MAXNP + MAXNP**2 + MAXM + MAXM**2 +   &
             4*MAXN*MAXNQ + 6*MAXN*MAXM + 2*MAXN*MAXNQ*MAXNP +   &
             2*MAXN*MAXNQ*MAXM + MAXNQ**2 +                      &
             5*MAXNQ + MAXNQ*(MAXNP+MAXM) + LDWE*LD2WE*MAXNQ,    &
             LIWORK=20+MAXNP+MAXNQ*(MAXNP+MAXM))                 

! Variable Declarations 
INTEGER :: I,INFO,IPRINT,J,JOB,L,LUNERR,LUNRPT,M,MAXIT,N,NDIGIT,NP,NQ
INTEGER :: IFIXB(MAXNP),IFIXX(LDIFX,MAXM),IWORK(:)
REAL(KIND=wp) :: PARTOL,SSTOL,TAUFAC
REAL(KIND=wp) :: BETA(MAXNP),SCLB(MAXNP),SCLD(LDSCLD,MAXM),      &
                  STPB(MAXNP),STPD(LDSTPD,MAXM),                 &  
                  WD(LDWD,LD2WD,MAXM),WE(LDWE,LD2WE,MAXNQ),      &
                  WORK(:),X(LDX,MAXM),Y(LDY,MAXNQ)
EXTERNAL :: FCN
POINTER  :: IWORK, WORK

! Allocate work arrays
ALLOCATE(IWORK(LIWORK),WORK(LWORK))

! Specify default values for ODR arguments
WE(1,1,1)  = -1.0E0_wp
WD(1,1,1)  = -1.0E0_wp
IFIXB(1)   = -1
IFIXX(1,1) = -1
JOB        = -1
NDIGIT     = -1
TAUFAC     = -1.0E0_wp
SSTOL      = -1.0E0_wp
PARTOL     = -1.0E0_wp
MAXIT      = -1
IPRINT     = -1
LUNERR     = -1
LUNRPT     = -1
STPB(1)    = -1.0E0_wp
STPD(1,1)  = -1.0E0_wp
SCLB(1)    = -1.0E0_wp
SCLD(1,1)  = -1.0E0_wp

! Set up ODRPACK95 report files
LUNERR = 9
LUNRPT = 9
OPEN (UNIT=9, FILE='./example/report1.dat')

! Read problem data, and set nondefault value for argument IFIXX
OPEN (UNIT=5, FILE='./example/data1.dat')
READ (5,*) N,M,NP,NQ
READ (5,*) (BETA(I),I=1,NP)
DO I=1,N
    READ (5,*) (X(I,J),J=1,M),(Y(I,L),L=1,NQ)
    IF (X(I,1).EQ.0.0E0_wp .OR. X(I,1).EQ.100.0E0_wp) THEN
        IFIXX(I,1) = 0
    ELSE
        IFIXX(I,1) = 1
    END IF
END DO

! Specify task: Explicit orthogonal distance regression
!               With user supplied derivatives (checked)
!               Covariance matrix constructed with recomputed derivatives
!               Delta initialized to zero
!               Not a restart
! And indicate short initial report
!               Short iteration reports every iteration, and
!               Long final report

JOB = 00020
IPRINT = 1112

! Compute solution
CALL ODR(FCN=FCN,                                          &
         N=N,M=M,NP=NP,NQ=NQ,                              &
         BETA=BETA,                                        &
         Y=Y,X=X,                                          &
         WE=WE,WD=WD,                                      &
         IFIXB=IFIXB,IFIXX=IFIXX,                          &
         JOB=JOB,NDIGIT=NDIGIT,TAUFAC=TAUFAC,              &
         SSTOL=SSTOL,PARTOL=PARTOL,MAXIT=MAXIT,            &
         IPRINT=IPRINT,LUNERR=LUNERR,LUNRPT=LUNRPT,        &
         STPB=STPB,STPD=STPD,                              &
         SCLB=SCLB,SCLD=SCLD,                              &
         WORK=WORK,IWORK=IWORK,                            &
         INFO=INFO)

END PROGRAM


SUBROUTINE FCN(N,M,NP,NQ,                                  &
               LDN,LDM,LDNP,                               &
               BETA,XPLUSD,                                & 
               IFIXB,IFIXX,LDIFX,                          &
               IDEVAL,F,FJACB,FJACD,                       &
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

! Used modules
use odrpack95_kinds, only: wp, zero

! Input arguments, not to be changed by this routine:
INTEGER, intent(in) :: IDEVAL, LDIFX, LDM, LDN, LDNP, M, N, NP, NQ
INTEGER, intent(in) :: IFIXB(NP), IFIXX(LDIFX,M)
REAL(KIND=wp), intent(in) :: BETA(NP), XPLUSD(LDN,M)
! Output arguments:
REAL(KIND=wp), intent(out) :: F(LDN,NQ), FJACB(LDN,LDNP,NQ), FJACD(LDN,LDM,NQ)
integer, intent (out) :: ISTOP
! Local variables
INTRINSIC :: EXP
INTEGER :: I, L

! Do something with IFIXB and IFIXX to avoid warnings that they are not being
! used.  This is simply not to worry users that the example program is failing.
IF (IFIXB(1) .GT. 0 .AND. IFIXX(1,1) .GT. 0) THEN
!   Do nothing.
END IF

! Check for unacceptable values for this problem
IF (BETA(1) .LT. zero) THEN
    ISTOP = 1
    RETURN
ELSE
    ISTOP = 0
END IF

! Compute predicted values
IF (MOD(IDEVAL,10) .GE. 1) THEN
    DO L = 1,NQ
        DO I = 1,N
            F(I,L) = BETA(1) + BETA(2)*(EXP(BETA(3)*XPLUSD(I,1)) - 1.0E0_wp)**2
        END DO
    END DO
END IF

! Compute derivatives with respect to BETA
IF (MOD(IDEVAL/10,10) .GE. 1) THEN
    DO L = 1,NQ
        DO I = 1,N
            FJACB(I,1,L) = 1.0E0_wp
            FJACB(I,2,L) = (EXP(BETA(3)*XPLUSD(I,1)) - 1.0E0_wp)**2
            FJACB(I,3,L) = BETA(2)*2*                                &
                           (EXP(BETA(3)*XPLUSD(I,1)) - 1.0E0_wp)*    &
                           EXP(BETA(3)*XPLUSD(I,1))*XPLUSD(I,1)
        END DO
    END DO
END IF

! Compute derivatives with respect to DELTA
IF (MOD(IDEVAL/100,10) .GE. 1) THEN
    DO L = 1,NQ
        DO I = 1,N
            FJACD(I,1,L) = BETA(2)*2*                               & 
                           (EXP(BETA(3)*XPLUSD(I,1)) - 1.0E0_wp)*   &
                           EXP(BETA(3)*XPLUSD(I,1))*BETA(3)     
        END DO
    END DO
END IF

END SUBROUTINE
