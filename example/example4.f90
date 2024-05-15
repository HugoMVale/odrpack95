!   This sample problem comes from Zwolak et al. 2001 (High Performance Computing
! Symposium, "Estimating rate constants in cell cycle models").  The call to
! ODRPACK95 is modified from the call the authors make to ODRPACK.  This is
! done to illustrate the need for bounds.  The authors could just have easily
! used the call statement here to solve their problem.
!   Curious users are encouraged to remove the bounds in the call statement,
! run the code, and compare the results to the current call statement.
PROGRAM example4
   use odrpack_kinds, only: wp
   USE odrpack
   IMPLICIT NONE

   REAL(KIND=wp) :: BETA(3) = (/1.1E-0_wp, 3.3E+0_wp, 8.7_wp/)
   EXTERNAL :: FCN
!     INTEGER :: I
!     REAL (KIND=wp) :: C, M, TOUT

   OPEN (9, FILE="./example/report4.dat")

   CALL ODR(FCN, &
            N=5, M=1, NP=3, NQ=1, &
            BETA=BETA, &
            Y=RESHAPE((/55.0_wp, 45.0_wp, 40.0_wp, 30.0_wp, 20.0_wp/), (/5, 1/)), &
            X=RESHAPE((/0.15_wp, 0.20_wp, 0.25_wp, 0.30_wp, 0.50_wp/), (/5, 1/)), &
            LOWER=(/0.0_wp, 0.0_wp, 0.0_wp/), &
            UPPER=(/1000.0_wp, 1000.0_wp, 1000.0_wp/), &
            IPRINT=2122, &
            LUNRPT=9, &
            MAXIT=20)

   CLOSE (9)

   ! The following code will reproduce the plot in Figure 2 of Zwolak et al. 2001.
   !    DO I = 0, 100
   !       C = 0.05 + (0.7 - 0.05)*I/100
   !       TOUT = 1440.0D0
   !       !CALL MPF(M,C,1.1D-10,3.3D-3,8.7D0,0.0D0,TOUT,C/2)
   !       CALL MPF(M, C, 1.15395968E-02_wp, 2.61676386E-03_wp, &
   !                9.23138811E+00_wp, 0.0D0, TOUT, C/2)
   !       WRITE (*, *) C, TOUT
   !    END DO

END PROGRAM example4

SUBROUTINE FCN(N, M, NP, NQ, LDN, LDM, LDNP, BETA, XPLUSD, IFIXB, IFIXX, LDIFX, &
               IDEVAL, F, FJACB, FJACD, ISTOP)
   use odrpack_kinds, only: wp, ZERO
   IMPLICIT NONE
! Subroutine arguments:
   INTEGER, INTENT(IN) :: IDEVAL, LDIFX, LDM, LDN, LDNP, M, N, NP, NQ
   INTEGER, INTENT(IN) :: IFIXB(NP), IFIXX(LDIFX, M)
   REAL(KIND=wp), INTENT(IN) :: BETA(NP), XPLUSD(LDN, M)
   REAL(KIND=wp), INTENT(OUT) :: F(LDN, NQ), FJACB(LDN, LDNP, NQ), FJACD(LDN, LDM, NQ)
   INTEGER, INTENT(OUT) :: ISTOP
! Local variables
   REAL(KIND=wp) :: MOUT
   INTEGER :: I

   ISTOP = 0
   FJACB(:, :, :) = ZERO
   FJACD(:, :, :) = ZERO
   IF (MOD(IDEVAL, 10) .GE. 1) THEN
      DO I = 1, N
         F(I, 1) = 1440.0_wp
         CALL MPF(MOUT, XPLUSD(I, 1), &
                  BETA(1), BETA(2), BETA(3), ZERO, F(I, 1), XPLUSD(I, 1)/2)
      END DO
   END IF
END SUBROUTINE FCN

! -------------------------------------------------------------------------------
! MPF
!
! If ROOT is not zero then returns value of time when M==ROOT in TOUT.  Else,
! runs until TOUT and returns value in M.  If PRINT_EVERY is non-zero then
! the solution is printed every PRINT_EVERY time units or every H (which ever
! is greater).
!
! This routine is not meant to be precise, it is only intended to be good
! enough for providing a working example of ODRPACK95 with bounds.  4th order
! Runge Kutta and linear interpolation are used for numerical integration and
! root finding, respectively.
!
! M - MPF
! C - Total Cyclin
! KWEE, K25, K25P - Model parameters (BETA(1:3))
SUBROUTINE MPF(M, C, KWEE, K25, K25P, PRINT_EVERY, TOUT, ROOT)
   use odrpack_kinds, only: wp, ZERO
   implicit none

   REAL(KIND=wp), INTENT(OUT) :: M
   REAL(KIND=wp), INTENT(IN)  :: C, KWEE, K25, K25P, PRINT_EVERY, ROOT
   REAL(KIND=wp), INTENT(INOUT) :: TOUT
   !  Local variables
   REAL(KIND=wp), PARAMETER :: H = 1.0E-1_wp
   REAL(KIND=wp) :: LAST_PRINT, LAST_M, LAST_T, T
   REAL(KIND=wp) :: K1, K2, K3, K4

   M = ZERO
   T = ZERO
   LAST_PRINT = ZERO
   IF (PRINT_EVERY .GT. ZERO) THEN
      WRITE (*, *) T, M
   END IF
   DO WHILE (T .LT. TOUT)
      LAST_T = T
      LAST_M = M
      K1 = H*DMDT(M, C, KWEE, K25, K25P)
      K2 = H*DMDT(M + K1/2, C, KWEE, K25, K25P)
      K3 = H*DMDT(M + K2/2, C, KWEE, K25, K25P)
      K4 = H*DMDT(M + K3, C, KWEE, K25, K25P)
      M = M + (K1 + 2*K2 + 2*K3 + K4)/6
      T = T + H
      IF (T .GE. PRINT_EVERY + LAST_PRINT .AND. PRINT_EVERY .GT. ZERO) THEN
         WRITE (*, *) T, M
         LAST_PRINT = T
      END IF
      IF (ROOT .GT. ZERO) THEN
         IF (LAST_M .LE. ROOT .AND. ROOT .LT. M) THEN
            TOUT = (T - LAST_T)/(M - LAST_M)*(ROOT - LAST_M) + LAST_T
            RETURN
         END IF
      END IF
   END DO

CONTAINS

   ! Equation from Zwolak et al. 2001.
   FUNCTION DMDT(M_, C_, KWEE_, K25_, K25P_) RESULT(RES)
      REAL(KIND=wp) :: M_, C_, KWEE_, K25_, K25P_, RES
      RES = KWEE_*M_ + (K25_ + K25P_*M_**2)*(C_ - M_)
   END FUNCTION DMDT

END SUBROUTINE MPF
