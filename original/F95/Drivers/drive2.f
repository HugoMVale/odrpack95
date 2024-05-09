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
C      ==> Y        Response variable (unused when model is implicit)
C      ==> X        Explanatory variable
C      ==> WE       Initial penalty parameter for implicit model
C      ==> WD       "Delta" weights
C      ==> JOB      Task to be performed
C      ==> IPRINT   Print control
C      ==> LUNERR   Logical unit for error reports
C      ==> LUNRPT   Logical unit for computation reports
C     <==> WORK     REAL (KIND=R8) work vector
C     <==  IWORK    Integer work vector
C     <==  INFO     Stopping condition 
 
C  Parameters specifying maximum problem sizes handled by this driver
C     MAXN          Maximum number of observations 
C     MAXM          Maximum number of columns in explanatory variable
C     MAXNP         Maximum number of function parameters
C     MAXNQ         Maximum number of responses per observation

C  Parameter declarations and specifications
      INTEGER    LDWD,LDWE,LDX,LDY,LD2WD,LD2WE,
     &           LIWORK,LWORK,MAXM,MAXN,MAXNP,MAXNQ
      PARAMETER (MAXM=5,MAXN=25,MAXNP=5,MAXNQ=2,
     &           LDY=MAXN,LDX=MAXN,
     &           LDWE=1,LD2WE=1,LDWD=1,LD2WD=1,
     &           LWORK=18 + 11*MAXNP + MAXNP**2 + MAXM + MAXM**2 + 
     &                 4*MAXN*MAXNQ + 6*MAXN*MAXM + 2*MAXN*MAXNQ*MAXNP +  
     &                 2*MAXN*MAXNQ*MAXM + MAXNQ**2 + 
     &                 5*MAXNQ + MAXNQ*(MAXNP+MAXM) + LDWE*LD2WE*MAXNQ,
     &           LIWORK=20+MAXNP+MAXNQ*(MAXNP+MAXM))

C  Variable declarations 
      INTEGER        I,INFO,IPRINT,J,JOB,LUNERR,LUNRPT,M,N,NP,NQ
      INTEGER        IWORK(:)   
      REAL (KIND=R8) BETA(MAXNP),
     &               WD(LDWD,LD2WD,MAXM),WE(LDWE,LD2WE,MAXNQ),
     &               WORK(:),X(LDX,MAXM),Y(LDY,MAXNQ)
      EXTERNAL       FCN
      POINTER        IWORK,WORK


C  Allocate work arrays
      ALLOCATE(IWORK(LIWORK),WORK(LWORK))

C  Specify default values for DODR arguments
      WE(1,1,1)  = -1.0E0_R8
      WD(1,1,1)  = -1.0E0_R8
      JOB        = -1
      IPRINT     = -1
      LUNERR     = -1
      LUNRPT     = -1
 
C  Set up ODRPACK95 report files
      LUNERR  =   9
      LUNRPT  =   9
      OPEN (UNIT=9,FILE='REPORT2')

C  Read problem data
      OPEN (UNIT=5,FILE='DATA2')
      READ (5,FMT=*) N,M,NP,NQ
      READ (5,FMT=*) (BETA(I),I=1,NP)
      DO 10 I=1,N
         READ (5,FMT=*) (X(I,J),J=1,M)
   10 CONTINUE

C  Specify task: Implicit orthogonal distance regression
C                With forward finite difference derivatives
C                Covariance matrix constructed with recomputed derivatives
C                DELTA initialized to zero
C                Not a restart
      JOB     = 00001

C  Compute solution
      CALL ODR(FCN=FCN,
     &         N=N,M=M,NP=NP,NQ=NQ,
     &         BETA=BETA,
     &         Y=Y,X=X,
     &         WE=WE,WD=WD,
     &         JOB=JOB,
     &         IPRINT=IPRINT,LUNERR=LUNERR,LUNRPT=LUNRPT,
     &         WORK=WORK,IWORK=IWORK,
     &         INFO=INFO)
      END


      SUBROUTINE FCN(N,M,NP,NQ,
     &               LDN,LDM,LDNP,
     &               BETA,XPLUSD,
     &               IFIXB,IFIXX,LDIFX,
     &               IDEVAL,F,FJACB,FJACD,
     &               ISTOP)

C  Subroutine Arguments
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

C  Used Modules
      USE REAL_PRECISION

C  Input arguments, not to be changed by this routine:
      INTEGER          I,IDEVAL,ISTOP,L,LDIFX,LDM,LDN,LDNP,M,N,NP,NQ
      REAL (KIND=R8) BETA(NP),XPLUSD(LDN,M)
      INTEGER          IFIXB(NP),IFIXX(LDIFX,M)
C  Output arguments:
      REAL (KIND=R8) F(LDN,NQ),FJACB(LDN,LDNP,NQ),FJACD(LDN,LDM,NQ)


C  Do something with FJACD, FJACB, IFIXB and IFIXX to avoid warnings that they 
C  are not being used.  This is simply not to worry users that the example 
C  program is failing.
      IF (IFIXB(1) .GT. 0 .AND. IFIXX(1,1) .GT. 0 
     &    .AND. FJACB(1,1,1) .GT. 0 .AND. FJACD(1,1,1) .GT. 0 ) THEN
C        Do nothing.
      END IF


C  Check for unacceptable values for this problem
      IF (BETA(1) .GT. 0.0E0_R8) THEN
         ISTOP = 1
         RETURN
      ELSE
         ISTOP = 0
      END IF

C  Compute predicted values
      IF (MOD(IDEVAL,10).GE.1) THEN
         DO 110 L = 1,NQ
            DO 100 I = 1,N
               F(I,L) = BETA(3)*(XPLUSD(I,1)-BETA(1))**2 +
     &                  2*BETA(4)*(XPLUSD(I,1)-BETA(1))*
     &                            (XPLUSD(I,2)-BETA(2)) +
     &                  BETA(5)*(XPLUSD(I,2)-BETA(2))**2 - 1.0E0_R8
  100       CONTINUE
  110    CONTINUE
      END IF

      RETURN
      END
