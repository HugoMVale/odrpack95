
      PROGRAM SAMPLE
      USE ODRPACK95
      USE REAL_PRECISION

C  ODRPACK95 Argument Definitions
C      ==> FCN      Name of the user supplied function subroutine
C      ==> N        Number of observations 
C      ==> M        Columns of data in the explanatory variable
C      ==> NP       Number of parameters
C      ==> NQ       Number of responses per observation
C     <==> BETA     Function parameters
C      ==> Y        Response variable
C      ==> X        Explanatory variable
C      ==> WE       "Epsilon" weights
C      ==> WD       "Delta" weights
C      ==> IFIXB    Indicators for "fixing" parameters (BETA)
C      ==> IFIXX    Indicators for "fixing" explanatory variable (X)
C      ==> JOB      Task to be performed 
C      ==> NDIGIT   Good digits in subroutine function results
C      ==> TAUFAC   Trust region initialization factor
C      ==> SSTOL    Sum of squares convergence criterion
C      ==> PARTOL   Parameter convergence criterion
C      ==> MAXIT    Maximum number of iterations
C      ==> IPRINT   Print control 
c      ==> LUNERR   Logical unit for error reports 
C      ==> LUNRPT   Logical unit for computation reports 
C      ==> STPB     Step sizes for finite difference derivatives wrt BETA
C      ==> STPD     Step sizes for finite difference derivatives wrt DELTA
C      ==> SCLB     Scale values for parameters BETA
C      ==> SCLD     Scale values for errors delta in explanatory variable 
C     <==> WORK     REAL (KIND=R8) work vector
C     <==  IWORK    Integer work vector
C     <==  INFO     Stopping condition 
 
C  Parameters specifying maximum problem sizes handled by this driver
C     MAXN          Maximum number of observations 
C     MAXM          Maximum number of columns in explanatory variable
C     MAXNP         Maximum number of function parameters
C     MAXNQ         Maximum number of responses per observation

C  Parameter Declarations and Specifications
      INTEGER    LDIFX,LDSCLD,LDSTPD,LDWD,LDWE,LDX,LDY,LD2WD,LD2WE,
     &           LIWORK,LWORK,MAXM,MAXN,MAXNP,MAXNQ
      PARAMETER (MAXM=5,MAXN=25,MAXNP=5,MAXNQ=1,
     &           LDY=MAXN,LDX=MAXN,
     &           LDWE=1,LD2WE=1,LDWD=1,LD2WD=1,
     &           LDIFX=MAXN,LDSTPD=1,LDSCLD=1,
     &           LWORK=18 + 11*MAXNP + MAXNP**2 + MAXM + MAXM**2 + 
     &                 4*MAXN*MAXNQ + 6*MAXN*MAXM + 2*MAXN*MAXNQ*MAXNP +  
     &                 2*MAXN*MAXNQ*MAXM + MAXNQ**2 + 
     &                 5*MAXNQ + MAXNQ*(MAXNP+MAXM) + LDWE*LD2WE*MAXNQ,
     &           LIWORK=20+MAXNP+MAXNQ*(MAXNP+MAXM))

C  Variable Declarations 
      INTEGER        I,INFO,IPRINT,J,JOB,L,LUNERR,LUNRPT,M,MAXIT,N,
     &               NDIGIT,NP,NQ
      INTEGER        IFIXB(MAXNP),IFIXX(LDIFX,MAXM),IWORK(:)
      REAL (KIND=R8) PARTOL,SSTOL,TAUFAC
      REAL (KIND=R8) BETA(MAXNP),SCLB(MAXNP),SCLD(LDSCLD,MAXM),
     &               STPB(MAXNP),STPD(LDSTPD,MAXM),
     &               WD(LDWD,LD2WD,MAXM),WE(LDWE,LD2WE,MAXNQ),
     &               WORK(:),X(LDX,MAXM),Y(LDY,MAXNQ)
      EXTERNAL       FCN
      POINTER        IWORK,WORK


C  Allocate work arrays
      ALLOCATE(IWORK(LIWORK),WORK(LWORK))

C  Specify default values for ODR arguments
      WE(1,1,1)  = -1.0E0_R8
      WD(1,1,1)  = -1.0E0_R8
      IFIXB(1)   = -1
      IFIXX(1,1) = -1
      JOB        = -1
      NDIGIT     = -1
      TAUFAC     = -1.0E0_R8
      SSTOL      = -1.0E0_R8
      PARTOL     = -1.0E0_R8
      MAXIT      = -1
      IPRINT     = -1
      LUNERR     = -1
      LUNRPT     = -1
      STPB(1)    = -1.0E0_R8
      STPD(1,1)  = -1.0E0_R8
      SCLB(1)    = -1.0E0_R8
      SCLD(1,1)  = -1.0E0_R8

C  Set up ODRPACK95 report files
      LUNERR  =   9
      LUNRPT  =   9
      OPEN (UNIT=9,FILE='REPORT1')

C  Read problem data, and set nondefault value for argument IFIXX
      OPEN (UNIT=5,FILE='DATA1')
      READ (5,FMT=*) N,M,NP,NQ
      READ (5,FMT=*) (BETA(I),I=1,NP)
      DO 10 I=1,N
         READ (5,FMT=*) (X(I,J),J=1,M),(Y(I,L),L=1,NQ)
         IF (X(I,1).EQ.0.0E0_R8 .OR. X(I,1).EQ.100.0E0_R8) THEN
            IFIXX(I,1) = 0
         ELSE
            IFIXX(I,1) = 1
         END IF
   10 CONTINUE

C  Specify task: Explicit orthogonal distance regression
C                With user supplied derivatives (checked)
C                Covariance matrix constructed with recomputed derivatives
C                Delta initialized to zero
C                Not a restart
C  And indicate short initial report
C               Short iteration reports every iteration, and
C               Long final report
      JOB     = 00020
      IPRINT  = 1112

C  Compute solution
      CALL ODR(FCN=FCN,
     &         N=N,M=M,NP=NP,NQ=NQ,
     &         BETA=BETA,
     &         Y=Y,X=X,
     &         WE=WE,WD=WD,
     &         IFIXB=IFIXB,IFIXX=IFIXX,
     &         JOB=JOB,NDIGIT=NDIGIT,TAUFAC=TAUFAC,
     &         SSTOL=SSTOL,PARTOL=PARTOL,MAXIT=MAXIT,
     &         IPRINT=IPRINT,LUNERR=LUNERR,LUNRPT=LUNRPT,
     &         STPB=STPB,STPD=STPD,
     &         SCLB=SCLB,SCLD=SCLD,
     &         WORK=WORK,IWORK=IWORK,
     &         INFO=INFO)
      END


      SUBROUTINE FCN(N,M,NP,NQ,
     &               LDN,LDM,LDNP,
     &               BETA,XPLUSD,
     &               IFIXB,IFIXX,LDIFX,
     &               IDEVAL,F,FJACB,FJACD,
     &               ISTOP)

C  Subroutine arguments
C      ==> N        Number of observations
C      ==> M        Number of columns in explanatory variable
C      ==> NP       Number of parameters
C      ==> NQ       Number of responses per observation
C      ==> LDN      Leading dimension declarator equal or exceeding N
C      ==> LDM      Leading dimension declarator equal or exceeding M
C      ==> LDNP     Leading dimension declarator equal or exceeding NP
C      ==> BETA     Current values of parameters
C      ==> XPLUSD   Current value of explanatory variable, i.e., X + DELTA
C      ==> IFIXB    Indicators for "fixing" parameters (BETA)
C      ==> IFIXX    Indicators for "fixing" explanatory variable (X)
C      ==> LDIFX    Leading dimension of array IFIXX
C      ==> IDEVAL   Indicator for selecting computation to be performed
C     <==  F        Predicted function values
C     <==  FJACB    Jacobian with respect to BETA
C     <==  FJACD    Jacobian with respect to errors DELTA
C     <==  ISTOP    Stopping condition, where
C                     0 means current BETA and X+DELTA were
C                       acceptable and values were computed successfully
C                     1 means current BETA and X+DELTA are
C                       not acceptable;  ODRPACK95 should select values 
C                       closer to most recently used values if possible
C                    -1 means current BETA and X+DELTA are 
C                       not acceptable; ODRPACK95 should stop

C  Used modules
      USE REAL_PRECISION

C  Input arguments, not to be changed by this routine:
      INTEGER          I,IDEVAL,ISTOP,L,LDIFX,LDM,LDN,LDNP,M,N,NP,NQ
      REAL (KIND=R8) BETA(NP),XPLUSD(LDN,M)
      INTEGER          IFIXB(NP),IFIXX(LDIFX,M)
C  Output arguments:
      REAL (KIND=R8) F(LDN,NQ),FJACB(LDN,LDNP,NQ),FJACD(LDN,LDM,NQ)
C  Local variables
      INTRINSIC        EXP


C  Do something with IFIXB and IFIXX to avoid warnings that they are not being
C  used.  This is simply not to worry users that the example program is failing.
      IF (IFIXB(1) .GT. 0 .AND. IFIXX(1,1) .GT. 0 ) THEN
C        Do nothing.
      END IF


C  Check for unacceptable values for this problem
      IF (BETA(1) .LT. 0.0E0_R8) THEN
         ISTOP = 1
         RETURN
      ELSE
         ISTOP = 0
      END IF

C  Compute predicted values
      IF (MOD(IDEVAL,10).GE.1) THEN
         DO 110 L = 1,NQ
            DO 100 I = 1,N
               F(I,L) = BETA(1) + 
     &                  BETA(2)*(EXP(BETA(3)*XPLUSD(I,1)) - 1.0E0_R8)**2
  100       CONTINUE
  110    CONTINUE
      END IF

C  Compute derivatives with respect to BETA
      IF (MOD(IDEVAL/10,10).GE.1) THEN
         DO 210 L = 1,NQ
            DO 200 I = 1,N
               FJACB(I,1,L) = 1.0E0_R8
               FJACB(I,2,L) = (EXP(BETA(3)*XPLUSD(I,1)) - 1.0E0_R8)**2
               FJACB(I,3,L) = BETA(2)*2*
     &                        (EXP(BETA(3)*XPLUSD(I,1)) - 1.0E0_R8)*
     &                        EXP(BETA(3)*XPLUSD(I,1))*XPLUSD(I,1)
  200       CONTINUE
  210    CONTINUE
      END IF

C  Compute derivatives with respect to DELTA
      IF (MOD(IDEVAL/100,10).GE.1) THEN
         DO 310 L = 1,NQ
            DO 300 I = 1,N
               FJACD(I,1,L) = BETA(2)*2*
     &                        (EXP(BETA(3)*XPLUSD(I,1)) - 1.0E0_R8)*
     &                        EXP(BETA(3)*XPLUSD(I,1))*BETA(3)
  300       CONTINUE
  310    CONTINUE
      END IF

      RETURN
      END




