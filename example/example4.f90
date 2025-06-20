module example4_model
!! Model for example4.

   use odrpack_kinds, only: wp, one, zero
   implicit none

contains

   pure subroutine fcn( &
      n, m, q, np, beta, xplusd, ifixb, ifixx, ldifx, ideval, f, fjacb, fjacd, istop)
   !! User-supplied subroutine for evaluating the model.

      integer, intent(in) :: ideval, ldifx, m, n, np, q
      integer, intent(in) :: ifixb(np), ifixx(ldifx, m)
      real(kind=wp), intent(in) :: beta(np), xplusd(n, m)
      real(kind=wp), intent(out) :: f(n, q), fjacb(n, np, q), fjacd(n, m, q)
      integer, intent(out) :: istop

      ! Local variables
      real(kind=wp) :: uout
      integer :: i

      istop = 0
      fjacb = zero
      fjacd = zero
      if (mod(ideval, 10) >= 1) then
         do i = 1, n
            f(i, 1) = 1440.0_wp
            call mpf(uout, xplusd(i, 1), &
                     beta(1), beta(2), beta(3), zero, f(i, 1), xplusd(i, 1)/2)
         end do
      end if

   end subroutine fcn

   pure subroutine mpf(u, c, kwee, k25, k25p, print_every, tout, root)
   !! If ROOT is not zero then returns value of time when U==ROOT in TOUT.  Else,
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

      real(kind=wp), intent(out) :: u
      real(kind=wp), intent(in) :: c, kwee, k25, k25p, print_every, root
      real(kind=wp), intent(inout) :: tout
      real(kind=wp), parameter :: h = 1.0E-1_wp
      real(kind=wp) :: last_print, last_u, last_t, t
      real(kind=wp) :: k1, k2, k3, k4

      u = zero
      t = zero

      last_print = zero
      if (print_every > zero) then
         !write (*, *) t, u
      end if

      do while (t < tout)
         last_t = t
         last_u = u
         k1 = h*dudt(u, c, kwee, k25, k25p)
         k2 = h*dudt(u + k1/2, c, kwee, k25, k25p)
         k3 = h*dudt(u + k2/2, c, kwee, k25, k25p)
         k4 = h*dudt(u + k3, c, kwee, k25, k25p)
         u = u + (k1 + 2*k2 + 2*k3 + k4)/6
         t = t + h
         if (t >= print_every + last_print .and. print_every > zero) &
            then
            !write (*, *) t, m
            last_print = t
         end if
         if (root > zero) then
            if (last_u <= root .and. root < u) then
               tout = (t - last_t)/(u - last_u)*(root - last_u) + last_t
               return
            end if
         end if
      end do

   contains

      ! Equation from Zwolak et al. 2001.
      real(kind=wp) pure function dudt(u_, c_, kwee_, k25_, k25p_) result(res)
         real(kind=wp), intent(in) :: u_, c_, kwee_, k25_, k25p_
         res = kwee_*u_ + (k25_ + k25p_*u_**2)*(c_ - u_)
      end function dudt

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
   use example4_model, only: fcn, mpf
   implicit none

   real(kind=wp) :: beta(3)
   integer :: n, m, np, q, lunrpt
   ! integer :: i
   ! real (kind=wp) :: u, c, tout

   open (newunit=lunrpt, file="./example/report4.dat")

   n = 5
   m = 1
   np = 3
   q = 1

   beta = [1.1E-0_wp, 3.3E+0_wp, 8.7_wp]

   call odr(fcn, n, m, q, np, &
            beta=beta, &
            y=reshape([55.0_wp, 45.0_wp, 40.0_wp, 30.0_wp, 20.0_wp], [n, q]), &
            x=reshape([0.15_wp, 0.20_wp, 0.25_wp, 0.30_wp, 0.50_wp], [n, m]), &
            lower=[0.0_wp, 0.0_wp, 0.0_wp], &
            upper=[1000.0_wp, 1000.0_wp, 1000.0_wp], &
            iprint=2122, &
            lunrpt=lunrpt, &
            maxit=20)

   close (lunrpt)

   ! The following code will reproduce the plot in Figure 2 of Zwolak et al. 2001.
   ! do i = 0, 100
   !    c = 0.05_wp + (0.7_wp - 0.05_wp)*i/100
   !    tout = 1440.0_wp
   !    !call mpf(u, c, 1.1e-10_wp, 3.3e-3_wp, 8.7_wp, 0.0_wp, tout, c/2)
   !    call mpf(u, c, 1.15395968E-02_wp, 2.61676386E-03_wp, &
   !             9.23138811_wp, 0.0_wp, tout, c/2)
   !    write (*, *) c, tout
   ! end do

end program example4
