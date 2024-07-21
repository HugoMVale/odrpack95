module example4_model
!! Model for example4.

   use odrpack_kinds, only: wp, one, zero
   implicit none

contains

   pure subroutine fcn(n, m, np, nq, ldn, ldm, ldnp, beta, xplusd, ifixb, ifixx, ldifx, &
                       ideval, f, fjacb, fjacd, istop)

      integer, intent(in) :: ideval, ldifx, ldm, ldn, ldnp, m, n, np, nq
      integer, intent(in) :: ifixb(np), ifixx(ldifx, m)
      real(kind=wp), intent(in) :: beta(np), xplusd(ldn, m)
      real(kind=wp), intent(out) :: f(ldn, nq), fjacb(ldn, ldnp, nq), fjacd(ldn, ldm, nq)
      integer, intent(out) :: istop

      real(kind=wp) :: mout
      integer :: i

      istop = 0
      fjacb(:, :, :) = zero
      fjacd(:, :, :) = zero
      if (mod(ideval, 10) .ge. 1) then
         do i = 1, n
            f(i, 1) = 1440.0_wp
            call mpf(mout, xplusd(i, 1), &
                     beta(1), beta(2), beta(3), zero, f(i, 1), xplusd(i, 1)/2)
         end do
      end if
   end subroutine fcn

   pure subroutine mpf(m, c, kwee, k25, k25p, print_every, tout, root)
   !! If ROOT is not zero then returns value of time when M==ROOT in TOUT.  Else,
   !! runs until TOUT and returns value in M.  If PRINT_EVERY is non-zero then
   !! the solution is printed every PRINT_EVERY time units or every H (which ever
   !! is greater).
   !!
   !! This routine is not meant to be precise, it is only intended to be good
   !! enough for providing a working example of ODRPACK95 with bounds.  4th order
   !! Runge Kutta and linear interpolation are used for numerical integration and
   !! root finding, respectively.
   !!
   !! M - MPF
   !! C - Total Cyclin
   !! KWEE, K25, K25P - Model parameters (BETA(1:3))

      real(kind=wp), intent(out) :: m
      real(kind=wp), intent(in) :: c, kwee, k25, k25p, print_every, root
      real(kind=wp), intent(inout) :: tout
      real(kind=wp), parameter :: h = 1.0E-1_wp
      real(kind=wp) :: last_print, last_m, last_t, t
      real(kind=wp) :: k1, k2, k3, k4

      m = zero
      t = zero

      last_print = zero
      if (print_every .gt. zero) then
         !write (*, *) t, m
      end if

      do while (t .lt. tout)
         last_t = t
         last_m = m
         k1 = h*dmdt(m, c, kwee, k25, k25p)
         k2 = h*dmdt(m + k1/2, c, kwee, k25, k25p)
         k3 = h*dmdt(m + k2/2, c, kwee, k25, k25p)
         k4 = h*dmdt(m + k3, c, kwee, k25, k25p)
         m = m + (k1 + 2*k2 + 2*k3 + k4)/6
         t = t + h
         if (t .ge. print_every + last_print .and. print_every .gt. zero) &
            then
            !write (*, *) t, m
            last_print = t
         end if
         if (root .gt. zero) then
            if (last_m .le. root .and. root .lt. m) then
               tout = (t - last_t)/(m - last_m)*(root - last_m) + last_t
               return
            end if
         end if
      end do

   contains

      ! Equation from Zwolak et al. 2001.
      real(kind=wp) pure function dmdt(m_, c_, kwee_, k25_, k25p_) result(res)
         real(kind=wp), intent(in) :: m_, c_, kwee_, k25_, k25p_
         res = kwee_*m_ + (k25_ + k25p_*m_**2)*(c_ - m_)
      end function dmdt

   end subroutine mpf

end module example4_model

program example4
!! Default ODR job, with parameter bounds.
!!   This sample problem comes from Zwolak et al. 2001 (High Performance Computing
!! Symposium, "Estimating rate constants in cell cycle models"). The call to
!! ODRPACK95 is modified from the call the authors make to ODRPACK. This is
!! done to illustrate the need for bounds. The authors could just have easily
!! used the call statement here to solve their problem.
!!   Curious users are encouraged to remove the bounds in the call statement,
!! run the code, and compare the results to the current call statement.
   use odrpack_kinds, only: wp
   use odrpack, only: odr
   use example4_model, only: fcn
   implicit none

   real(kind=wp) :: beta(3) = [1.1E-0_wp, 3.3E+0_wp, 8.7_wp]
   integer :: lunrpt
   ! INTEGER :: I
   ! REAL (KIND=wp) :: C, M, TOUT

   open (newunit=lunrpt, file="./example/report4.dat")

   call odr(fcn, &
            n=5, m=1, np=3, nq=1, &
            beta=beta, &
            y=reshape([55.0_wp, 45.0_wp, 40.0_wp, 30.0_wp, 20.0_wp], [5, 1]), &
            x=reshape([0.15_wp, 0.20_wp, 0.25_wp, 0.30_wp, 0.50_wp], [5, 1]), &
            lower=[0.0_wp, 0.0_wp, 0.0_wp], &
            upper=[1000.0_wp, 1000.0_wp, 1000.0_wp], &
            iprint=2122, &
            lunrpt=lunrpt, &
            maxit=20)

   close (9)

   ! The following code will reproduce the plot in Figure 2 of Zwolak et al. 2001.
   !    DO I = 0, 100
   !       C = 0.05 + (0.7 - 0.05)*I/100
   !       TOUT = 1440.0D0
   !       !CALL MPF(M,C,1.1D-10,3.3D-3,8.7D0,0.0D0,TOUT,C/2)
   !       CALL MPF(M, C, 1.15395968E-02_wp, 2.61676386E-03_wp, &
   !                9.23138811E+00_wp, 0.0D0, TOUT, C/2)
   !       WRITE (*, *) C, TOUT
   !    END DO

end program example4
