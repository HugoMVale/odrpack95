C  This sample problem comes from Zwolak et al. 2001 (High Performance Computing
C  Symposium, "Estimating rate constants in cell cycle models").  The call to
C  ODRPACK95 is modified from the call the authors make to ODRPACK.  This is
C  done to illustrate the need for bounds.  The authors could just have easily
C  used the call statement here to solve their problem.
C
C  Curious users are encouraged to remove the bounds in the call statement,
C  run the code, and compare the results to the current call statement.
      PROGRAM EXAMPLE4
      use odrpack95_kinds, only: wp
      USE ODRPACK95
      IMPLICIT NONE
C      INTEGER :: I
C      REAL (KIND=wp) :: C, M, TOUT
      INTERFACE
          SUBROUTINE FCN(N,M,NP,NQ,LDN,LDM,LDNP,BETA,XPLUSD,IFIXB,
     +        IFIXX,LDIFX,IDEVAL,F,FJACB,FJACD,ISTOP)
          use odrpack95_kinds, only: wp
          INTEGER, INTENT(IN) :: IDEVAL,LDIFX,LDM,LDN,LDNP,M,N,NP,NQ
          INTEGER, INTENT(IN) :: IFIXB(NP),IFIXX(LDIFX,M)
          REAL(KIND=wp), INTENT(IN) :: BETA(NP),XPLUSD(LDN,M)
          INTEGER, INTENT(OUT) :: ISTOP
          REAL(KIND=wp), INTENT(OUT) :: F(LDN,NQ),FJACB(LDN,LDNP,NQ),
     +        FJACD(LDN,LDM,NQ)
          END SUBROUTINE FCN
      END INTERFACE
      REAL(KIND=wp) :: BETA(3) = (/ 1.1E-0_wp, 3.3E+0_wp, 8.7_wp /)
      OPEN(9,FILE="./example/report4.dat")
      CALL ODR(
     +    FCN,
     +    N = 5, M = 1, NP = 3, NQ = 1,
     +    BETA = BETA,
     +    Y = RESHAPE((/ 55.0_wp, 45.0_wp, 40.0_wp, 30.0_wp, 20.0_wp /),
     +        (/5,1/)),
     +    X = RESHAPE((/ 0.15_wp, 0.20_wp, 0.25_wp, 0.30_wp, 0.50_wp /),
     +        (/5,1/)),
     +    LOWER = (/ 0.0_wp, 0.0_wp, 0.0_wp /),
     +    UPPER = (/ 1000.0_wp, 1000.0_wp, 1000.0_wp /),
     +    IPRINT = 2122,
     +    LUNRPT = 9,
     +    MAXIT  = 20
     +)
      CLOSE(9)
C  The following code will reproduce the plot in Figure 2 of Zwolak et
C  al. 2001.
C      DO I = 0, 100
C          C = 0.05+(0.7-0.05)*I/100
C          TOUT = 1440.0D0
C          !CALL MPF(M,C,1.1D-10,3.3D-3,8.7D0,0.0D0,TOUT,C/2)
C          CALL MPF(M,C,1.15395968E-02_wp, 2.61676386E-03_wp, 
C     +             9.23138811E+00_wp,0.0D0,TOUT,C/2)
C          WRITE(*,*) C, TOUT
C      END DO
      END PROGRAM

      SUBROUTINE FCN(N,M,NP,NQ,LDN,LDM,LDNP,BETA,XPLUSD,IFIXB,
     +    IFIXX,LDIFX,IDEVAL,F,FJACB,FJACD,ISTOP)
      use odrpack95_kinds, only: wp
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IDEVAL,LDIFX,LDM,LDN,LDNP,M,N,NP,NQ
      INTEGER, INTENT(IN) :: IFIXB(NP),IFIXX(LDIFX,M)
      REAL(KIND=wp), INTENT(IN) :: BETA(NP),XPLUSD(LDN,M)
      INTEGER, INTENT(OUT) :: ISTOP
      REAL(KIND=wp), INTENT(OUT) :: F(LDN,NQ),FJACB(LDN,LDNP,NQ),
     +    FJACD(LDN,LDM,NQ)
      ! Local variables
      REAL(KIND=wp) :: MOUT
      INTEGER :: I
      ISTOP = 0
      FJACB(:,:,:) = 0.0E0_wp
      FJACD(:,:,:) = 0.0E0_wp
      IF ( MOD(IDEVAL,10).GE.1 ) THEN
          DO I = 1, N
              F(I,1) = 1440.0_wp
              CALL MPF(MOUT,XPLUSD(I,1),BETA(1),BETA(2),BETA(3),0.0_wp,
     +            F(I,1),XPLUSD(I,1)/2)
          END DO
      END IF
      END SUBROUTINE FCN

C-------------------------------------------------------------------------------
C
C  MPF
C
C  If ROOT is not zero then returns value of time when M==ROOT in TOUT.  Else,
C  runs until TOUT and returns value in M.  If PRINT_EVERY is non-zero then
C  the solution is printed every PRINT_EVERY time units or every H (which ever
C  is greater).
C
C  This routine is not meant to be precise, it is only intended to be good
C  enough for providing a working example of ODRPACK95 with bounds.  4th order 
C  Runge Kutta and linear interpolation are used for numerical integration and
C  root finding, respectively.
C
C  M - MPF
C  C - Total Cyclin
C  KWEE, K25, K25P - Model parameters (BETA(1:3))
C
      SUBROUTINE MPF(M,C,KWEE,K25,K25P,PRINT_EVERY,TOUT,ROOT)
          use odrpack95_kinds, only: wp
          REAL (KIND=wp), INTENT(OUT) :: M
          REAL (KIND=wp), INTENT(IN)  :: C, KWEE, K25, K25P,
     +        PRINT_EVERY, ROOT
          REAL (KIND=wp), INTENT(INOUT) :: TOUT
          !  Local variables
          REAL (KIND=wp), PARAMETER :: H = 1.0D-1
          REAL (KIND=wp) :: LAST_PRINT, LAST_M, LAST_T, T
          REAL (KIND=wp) :: K1, K2, K3, K4
          INTERFACE
              FUNCTION DMDT(M,C,KWEE,K25,K25P) RESULT(RES)
                  use odrpack95_kinds, only: wp
                  REAL (KIND=wp) :: M, C, KWEE, K25, K25P, RES
              END FUNCTION
          END INTERFACE
          M = 0.0D0
          T = 0.0D0
          LAST_PRINT = 0.0D0
          IF ( PRINT_EVERY .GT. 0.0D0 ) THEN
              WRITE(*,*) T, M
          END IF
          DO WHILE ( T .LT. TOUT )
              LAST_T = T
              LAST_M = M
              K1 = H*DMDT(M,C,KWEE,K25,K25P)
              K2 = H*DMDT(M+K1/2,C,KWEE,K25,K25P)
              K3 = H*DMDT(M+K2/2,C,KWEE,K25,K25P)
              K4 = H*DMDT(M+K3,C,KWEE,K25,K25P)
              M = M+(K1+2*K2+2*K3+K4)/6
              T = T + H
              IF ( T .GE. PRINT_EVERY+LAST_PRINT .AND. 
     +            PRINT_EVERY .GT. 0.0D0 )             
     +        THEN
                  WRITE(*,*) T, M
                  LAST_PRINT = T
              END IF
              IF ( ROOT .GT. 0.0D0 ) THEN
                  IF ( LAST_M .LE. ROOT .AND. ROOT .LT. M ) THEN
                      TOUT = (T-LAST_T)/(M-LAST_M)*(ROOT-LAST_M)+LAST_T
                      RETURN
                  END IF
              END IF
          END DO
      END SUBROUTINE MPF


C  Equation from Zwolak et al. 2001.
      FUNCTION DMDT(M,C,KWEE,K25,K25P) RESULT(RES)
          use odrpack95_kinds, only: wp
          REAL (KIND=wp) :: M, C, KWEE, K25, K25P, RES
          RES = KWEE*M+(K25+K25P*M**2)*(C-M)
      END FUNCTION DMDT
