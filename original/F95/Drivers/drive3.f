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
C      ==> NDIGIT   Good digits in subroutine fcn results
C      ==> TAUFAC   Trust region initialization factor
C      ==> SSTOL    Sum of squares convergence criterion
C      ==> PARTOL   Parameter convergence criterion
C      ==> MAXIT    Maximum number of iterations
C      ==> IPRINT   Print control
C      ==> LUNERR   Logical unit for error reports
C      ==> LUNRPT   Logical unit for computation reports
C      ==> STPB     Step sizes for finite difference derivatives wrt BETA
C      ==> STPD     Step sizes for finite difference derivatives wrt DELTA
C      ==> SCLB     Scale values for parameters BETA 
C      ==> SCLD     Scale values for errors DELTA in explanatory variable
C     <==> WORK     REAL (KIND=R8) work vector
C     <==  IWORK    Integer work vector
C     <==  INFO     Stopping condition 
 
C  Parameters specifying maximum problem sizes handled by this driver
C     MAXN          Maximum number of observations 
C     MAXM          Maximum number of columns in explanatory variable
C     MAXNP         Maximum number of function parameters
C     MAXNQ         Maximum number of responses per observation

C  Parameter declarations and specifications
      INTEGER    LDIFX,LDSCLD,LDSTPD,LDWD,LDWE,LDX,LDY,LD2WD,LD2WE,
     &           LIWORK,LWORK,MAXM,MAXN,MAXNP,MAXNQ
      PARAMETER (MAXM=5,MAXN=100,MAXNP=25,MAXNQ=5,
     &           LDY=MAXN,LDX=MAXN,
     &           LDWE=MAXN,LD2WE=MAXNQ,LDWD=MAXN,LD2WD=1,
     &           LDIFX=MAXN,LDSCLD=1,LDSTPD=1,
     &           LWORK=18 + 11*MAXNP + MAXNP**2 + MAXM + MAXM**2 + 
     &                 4*MAXN*MAXNQ + 6*MAXN*MAXM + 2*MAXN*MAXNQ*MAXNP +  
     &                 2*MAXN*MAXNQ*MAXM + MAXNQ**2 + 
     &                 5*MAXNQ + MAXNQ*(MAXNP+MAXM) + LDWE*LD2WE*MAXNQ,
     &           LIWORK=20+MAXNP+MAXNQ*(MAXNP+MAXM))

C  Variable declarations 
      INTEGER        I,INFO,IPRINT,J,JOB,L,LUNERR,LUNRPT,M,MAXIT,N,
     &               NDIGIT,NP,NQ
      INTEGER        IFIXB(MAXNP),IFIXX(LDIFX,MAXM),IWORK(:)   
      REAL (KIND=R8) PARTOL,SSTOL,TAUFAC
      REAL (KIND=R8) BETA(MAXNP),DELTA(:,:),
     &               SCLB(MAXNP),SCLD(LDSCLD,MAXM),
     &               STPB(MAXNP),STPD(LDSTPD,MAXM),
     &               WD(LDWD,LD2WD,MAXM),WE(LDWE,LD2WE,MAXNQ),
     &               WORK(:),X(LDX,MAXM),Y(LDY,MAXNQ)
      EXTERNAL       FCN
      POINTER        DELTA,IWORK,WORK


C  Specify default values for DODRC arguments
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
      OPEN (UNIT=9,FILE='REPORT3')

C  Read problem data
      OPEN (UNIT=5,FILE='DATA3')
      READ (5,FMT=*) N,M,NP,NQ
      READ (5,FMT=*) (BETA(I),I=1,NP)
      DO 10 I=1,N
         READ (5,FMT=*) (X(I,J),J=1,M),(Y(I,L),L=1,NQ)
   10 CONTINUE

C  Allocate work arrays
      ALLOCATE(DELTA(N,M),IWORK(LIWORK),WORK(LWORK))

C  Specify task as explicit orthogonal distance regression
C                  With central difference derivatives
C                  Covariance matrix constructed with recomputed derivatives
C                  DELTA initialized by user
C                  Not a restart
C  And indicate long initial report
C               No iteration reports
C               Long final report
      JOB     = 01010
      IPRINT  = 2002

C  Initialize DELTA, and specify first decade of frequencies as fixed
      DO 20 I=1,N
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
   20 CONTINUE

C  Set weights
      DO 30 I=1,N
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
   30 CONTINUE

C  Compute solution
      CALL ODR(FCN=FCN,
     &         N=N,M=M,NP=NP,NQ=NQ,
     &         BETA=BETA,
     &         Y=Y,X=X,
     &         DELTA=DELTA,
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
C                     0 Means current BETA and X+DELTA were
C                       acceptable and values were computed successfully
C                     1 Means current BETA and X+DELTA are
C                       not acceptable;  ODRPACK95 should select values  
C                       closer to most recently used values if possible 
C                    -1 Means current BETA and X+DELTA are
C                       not acceptable;  ODRPACK95 should stop

C  Used modules
      USE REAL_PRECISION

C  Input arguments, not to be changed by this routine:
      INTEGER          I,IDEVAL,ISTOP,LDIFX,LDM,LDN,LDNP,M,N,NP,NQ
      REAL (KIND=R8) BETA(NP),XPLUSD(LDN,M)
      INTEGER          IFIXB(NP),IFIXX(LDIFX,M)
C  Output arguments:
      REAL (KIND=R8) F(LDN,NQ),FJACB(LDN,LDNP,NQ),FJACD(LDN,LDM,NQ)
C  Local variables
      REAL (KIND=R8) FREQ,PI,OMEGA,CTHETA,STHETA,THETA,PHI,R
      INTRINSIC        ATAN2,EXP,SQRT


C  Do something with FJACD, FJACB, IFIXB and IFIXX to avoid warnings that they 
C  are not being used.  This is simply not to worry users that the example 
C  program is failing.
      IF (IFIXB(1) .GT. 0 .AND. IFIXX(1,1) .GT. 0 
     &    .AND. FJACB(1,1,1) .GT. 0 .AND. FJACD(1,1,1) .GT. 0 ) THEN
C        Do nothing.
      END IF


C  Check for unacceptable values for this problem
      DO 10 I=1,N
         IF (XPLUSD(I,1).LT.0.0E0_R8) THEN
            ISTOP = 1
            RETURN
         END IF
   10 CONTINUE
      ISTOP = 0

      PI = 3.141592653589793238462643383279E0_R8

      THETA = PI*BETA(4)*0.5E0_R8
      CTHETA = COS(THETA)
      STHETA = SIN(THETA)

C  Compute predicted values
      IF (MOD(IDEVAL,10).GE.1) THEN
         DO 100 I = 1,N
            FREQ  = XPLUSD(I,1)
            OMEGA = (2.0E0_R8*PI*FREQ*EXP(-BETA(3)))**BETA(4)
            PHI   = ATAN2((OMEGA*STHETA),(1+OMEGA*CTHETA))
            R     = (BETA(1)-BETA(2)) * 
     &              SQRT((1+OMEGA*CTHETA)**2+
     &                       (OMEGA*STHETA)**2)**(-BETA(5))
            F(I,1) = BETA(2) + R*COS(BETA(5)*PHI)
            F(I,2) =           R*SIN(BETA(5)*PHI)
  100    CONTINUE
      END IF

      RETURN
      END












































