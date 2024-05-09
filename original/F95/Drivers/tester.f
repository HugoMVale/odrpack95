*TESTER
      PROGRAM TESTER
C***BEGIN PROLOGUE  TESTER
C***REFER TO ODR
C***ROUTINES CALLED  ODR
C***DATE WRITTEN   20040322   (YYYYMMDD)
C***REVISION DATE  20040322   (YYYYMMDD)
C***PURPOSE  EXCERCISE ERROR REPORTING OF THE F90 VERSION OF ODRPACK95
C***END PROLOGUE  TESTER

C...USED MODULES
      USE REAL_PRECISION
      USE ODRPACK95

C...LOCAL SCALARS
      INTEGER N, M, NQ, NP, INFO, LUN
C     STAT

C...LOCAL ARRAYS
      REAL (KIND=R8) 
     &   BETA(:),Y(:,:),X(:,:),UPPER(2),LOWER(2)

C...ALLOCATABLE ARRAYS
      ALLOCATABLE BETA,Y,X

C...EXTERNAL SUBPROGRAMS
      EXTERNAL FCN

      COMMON /BOUNDS/ UPPER,LOWER

C***FIRST EXECUTABLE STATEMENT  TESTER

      OPEN(UNIT=8,FILE="SUMMARY")
      WRITE(8,*) "NO SUMMARY AVAILABLE"
      CLOSE(8)

      LUN = 9
      OPEN(UNIT=LUN,FILE="REPORT")

C  ERROR IN PROBLEM SIZE

      N  = 0
      M  = 0
      NQ = 0
      NP = 0
      ALLOCATE(BETA(NP),Y(N,NQ),X(N,M))
      Y(:,:) = 0.0_R8
      X(:,:) = 0.0_R8
      BETA(:) = 0.0_R8

      CALL ODR(FCN,N,M,NP,NQ,BETA,Y,X,IPRINT=1,INFO=INFO,
     &   LUNRPT=LUN,LUNERR=LUN)

      WRITE(LUN,*) "INFO = ", INFO

C  ERROR IN JOB SPECIFICATION WITH WORK AND IWORK

      N  = 1
      M  = 1
      NQ = 1
      NP = 1
      DEALLOCATE(BETA,Y,X)
      ALLOCATE(BETA(NP),Y(N,NQ),X(N,M))
      Y(:,:) = 0.0_R8
      X(:,:) = 0.0_R8
      BETA(:) = 0.0_R8

      CALL ODR(FCN,N,M,NP,NQ,BETA,Y,X,IPRINT=1,INFO=INFO,JOB=10000,
     &   LUNRPT=LUN,LUNERR=LUN)

      WRITE(LUN,*) "INFO = ", INFO

C  ERROR IN JOB SPECIFICATION WITH DELTA

      N  = 1
      M  = 1
      NQ = 1
      NP = 1
      DEALLOCATE(BETA,Y,X)
      ALLOCATE(BETA(NP),Y(N,NQ),X(N,M))
      Y(:,:) = 0.0_R8
      X(:,:) = 0.0_R8
      BETA(:) = 0.0_R8

      CALL ODR(FCN,N,M,NP,NQ,BETA,Y,X,IPRINT=1,INFO=INFO,JOB=1000,
     &   LUNRPT=LUN,LUNERR=LUN)

      WRITE(LUN,*) "INFO = ", INFO

C  BOUNDS TOO SMALL FOR DERIVATIVE CHECKER WHEN DERIVATIVES DON'T AGREE.
      
      N  = 4
      M  = 1
      NQ = 1
      NP = 2
      DEALLOCATE(BETA,Y,X)
      ALLOCATE(BETA(NP),Y(N,NQ),X(N,M))
      BETA(:) = (/ -200.0_R8, -5.0_R8 /)
      UPPER(1:2) = (/ -200.0_R8, 0.0_R8 /)
      LOWER(1:2) = (/ -200.000029802322_R8, -5.0_R8 /)
      Y(:,1) = (/ 2.718281828459045_R8, 7.389056098930650_R8,
     &148.4131591025766_R8, 403.4287934927353_R8 /)
      X(:,1) = (/ 1.0_R8, 2.0_R8, 5.0_R8, 6.0_R8 /)

      CALL ODR(FCN,N,M,NP,NQ,BETA,Y,X,IPRINT=1,INFO=INFO,JOB=0020,
     &   LUNRPT=LUN,LUNERR=LUN,LOWER=LOWER,UPPER=UPPER)

      WRITE(LUN,*) "INFO = ", INFO

C  ERROR IN ARRAY ALLOCATION
C  The following code is intended to force memory allocation failure.  An
C  appropriate N for your machine must be chosen to ensure memory allocation
C  will fail within ODRPACK95.  A value of about 1/4 the total memory available 
C  to a process should do the trick.  However, most modern operating systems and
C  Fortran compilers will not likely deny ODRPACK95 memory before they fail for
C  another reason.  Therefore, the memory allocation checks in ODRPACK95 are not
C  easy to provoke.  An operating system may return successfull memory
C  allocation but fail to guarantee the memory causing a segfault when some
C  memory locations are accessed.  A Fortran compiler or operating system may
C  allow limited sized stacks during subroutine invocation causing the ODRPACK95
C  call to fail before ODRPACK95 executes its first line.
C
C     N  = 032000000
C     M  = 1
C     NQ = 1
C     NP = 1
C     DEALLOCATE(BETA,Y,X)
C     ALLOCATE(BETA(NP),Y(N,NQ),X(N,M),STAT=STAT)
C     IF (STAT.NE.0) THEN
C         WRITE(0,*) 
C    &       "SYSTEM ERROR: COULD NOT ALLOCATE MEMORY, TESTER ",
C    &       "FAILED TO RUN."
C         STOP
C     END IF
C     Y(:,:) = 0.0_R8
C     X(:,:) = 0.0_R8
C     BETA(:) = 0.0_R8
C
C     CALL ODR(FCN,N,M,NP,NQ,BETA,Y,X,IPRINT=1,INFO=INFO,
C    &   LUNRPT=LUN,LUNERR=LUN)
C
C     WRITE(LUN,*) "INFO = ", INFO

      CLOSE(LUN)

      END PROGRAM
*FCN
      SUBROUTINE FCN
     &   (N,M,NP,NQ, 
     &    LDN,LDM,LDNP,
     &    BETA,XPLUSD,
     &    IFIXB,IFIXX,LDIFX,
     &    IDEVAL,F,FJACB,FJACD,
     &    ISTOP)
C***BEGIN PROLOGUE  FCN
C***REFER TO  ODR
C***ROUTINES CALLED  (NONE)
C***DATE WRITTEN   20040322 (YYYYMMDD)
C***REVISION DATE  20040322 (YYYYMMDD)
C***PURPOSE  DUMMY ROUTINE FOR ODRPACK95 ERROR EXERCISER
C***END PROLOGUE  FCN

C...USED MODULES
      USE REAL_PRECISION

C...SCALAR ARGUMENTS
      INTEGER
     &   IDEVAL,ISTOP,LDIFX,LDM,LDN,LDNP,M,N,NP,NQ

C...ARRAY ARGUMENTS
      REAL (KIND=R8)
     &   BETA(NP),F(LDN,NQ),FJACB(LDN,LDNP,NQ),FJACD(LDN,LDM,NQ),
     &   XPLUSD(LDN,M)
      INTEGER
     &   IFIXB(NP),IFIXX(LDIFX,M)

C...ARRAYS IN COMMON
      REAL (KIND=R8)
     &   LOWER(2),UPPER(2)

C...LOCAL SCALARS
      INTEGER
     &   I

      COMMON /BOUNDS/ UPPER,LOWER


C***FIRST EXECUTABLE STATEMENT

C  Do something with FJACD, FJACB, IFIXB and IFIXX to avoid warnings that they 
C  are not being used.  This is simply not to worry users that the example 
C  program is failing.
      IF (IFIXB(1) .GT. 0 .AND. IFIXX(1,1) .GT. 0 
     &    .AND. FJACB(1,1,1) .GT. 0 .AND. FJACD(1,1,1) .GT. 0 ) THEN
C        Do nothing.
      END IF


      IF (ANY(LOWER(1:NP).GT.BETA(1:NP))) THEN
         WRITE(0,*) "LOWER BOUNDS VIOLATED"
         DO I=1,NP
            IF (LOWER(I).GT.BETA(I)) THEN
               WRITE(0,*) "   IN THE ", I, " POSITION WITH ", BETA(I), 
     &            "<", LOWER(I)
            END IF
         END DO
      END IF

      IF (ANY(UPPER(1:NP).LT.BETA(1:NP))) THEN
         WRITE(0,*) "UPPER BOUNDS VIOLATED"
         DO I=1,NP
            IF (UPPER(I).LT.BETA(I)) THEN
               WRITE(0,*) "   IN THE ", I, " POSITION WITH ", BETA(I), 
     &            ">", UPPER(I)
            END IF
         END DO
      END IF

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
     &         XPLUSD(I,1))
         END DO
      END IF

      IF (MOD(IDEVAL/100,10).NE.0) THEN
         DO I=1,N
            FJACD(I,1,1) = BETA(1)*BETA(2)*EXP(BETA(2)*XPLUSD(I,1))
         END DO
      END IF

      END SUBROUTINE
