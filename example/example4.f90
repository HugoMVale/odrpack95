module example4_model
!! Model for example4.

   use iso_c_binding, only: c_ptr
   use odrpack_kinds, only: dp, one, zero
   implicit none

contains

   pure subroutine fcn( &
      n, m, q, np, ldifx, beta, xplusd, ifixb, ifixx, ideval, f, fjacb, fjacd, istop, data)
   !! User-supplied subroutine for evaluating the model.

      integer, intent(in) :: n, m, q, np, ldifx, ideval, ifixb(np), ifixx(ldifx, m)
      real(dp), intent(in) :: beta(np), xplusd(n, m)
      real(dp), intent(out) :: f(n, q), fjacb(n, np, q), fjacd(n, m, q)
      integer, intent(out) :: istop
      type(c_ptr), intent(in), value :: data

      ! Local variables
      real(dp) :: uout
      integer :: i

      istop = 0
      fjacb = zero
      fjacd = zero
      if (mod(ideval, 10) > 0) then
         do i = 1, ubound(f, 1)
            f(i, 1) = 1440.0_dp
            call mpf(uout, xplusd(i, 1), beta(1), beta(2), beta(3), zero, f(i, 1), xplusd(i, 1)/2)
         end do
      end if

   end subroutine fcn

   pure subroutine mpf(u, c, kwee, k25, k25p, print_every, tout, root)
   !! If `root` is not zero, then it returns value of time when `u=root` in `tout`. Else, runs
   !! until `tout` and returns the value in `u`. If `print_every` is non-zero then the solution
   !! is printed every `print_every` time units or every `h` (which ever is greater).
   !!
   !! This routine is not meant to be precise, it is only intended to be good enough for
   !! providing a working example of ODRPACK95 with bounds. 4th order Runge Kutta and linear
   !! interpolation are used for numerical integration and root finding, respectively.
   !!
   !! `u`: MPF
   !! `c`: Total Cyclin
   !! `kwee`, `k25`, `k25p`: Model parameters corresponding to `beta(1:3)`

      real(dp), intent(out) :: u
      real(dp), intent(in) :: c, kwee, k25, k25p, print_every, root
      real(dp), intent(inout) :: tout
      real(dp), parameter :: h = 1.0E-1_dp
      real(dp) :: last_print, last_u, last_t, t
      real(dp) :: k1, k2, k3, k4

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
      real(dp) pure function dudt(u_, c_, kwee_, k25_, k25p_) result(res)
         real(dp), intent(in) :: u_, c_, kwee_, k25_, k25p_
         res = kwee_*u_ + (k25_ + k25p_*u_**2)*(c_ - u_)
      end function dudt

   end subroutine mpf

end module example4_model

program example4
!! Default ODR job, with parameter bounds.
!!   This sample problem comes from Zwolak et al. 2001 (High Performance Computing
!! Symposium, "Estimating rate constants in cell cycle models"). The call to
!! [[odr]] is modified from the call the authors make to ODRPACK. This is
!! done to illustrate the need for bounds. The authors could just have easily
!! used the call statement here to solve their problem.
!!   Curious users are encouraged to remove the bounds in the call statement,
!! run the code, and compare the results to the current call statement.

   use odrpack_kinds, only: dp
   use odrpack, only: odr, odrpack_model
   use example4_model, only: fcn, mpf
   implicit none

   type(odrpack_model) :: model
   real(dp) :: beta(3)
   integer :: n, m, np, q, lunrpt
   ! integer :: i
   ! real (dp) :: u, c, tout

   model%fcn => fcn

   open (newunit=lunrpt, file="./example/report4.dat")

   n = 5
   m = 1
   np = 3
   q = 1

   beta = [1.1E-0_dp, 3.3E+0_dp, 8.7_dp]

   call odr(model, n, m, q, np, beta, &
            y=reshape([55.0_dp, 45.0_dp, 40.0_dp, 30.0_dp, 20.0_dp], [n, q]), &
            x=reshape([0.15_dp, 0.20_dp, 0.25_dp, 0.30_dp, 0.50_dp], [n, m]), &
            lower=[0.0_dp, 0.0_dp, 0.0_dp], &
            upper=[1000.0_dp, 1000.0_dp, 1000.0_dp], &
            iprint=2122, &
            lunrpt=lunrpt, &
            maxit=20)

   close (lunrpt)

   ! The following code will reproduce the plot in Figure 2 of Zwolak et al. 2001.
   ! do i = 0, 100
   !    c = 0.05_dp + (0.7_dp - 0.05_dp)*i/100
   !    tout = 1440.0_dp
   !    !call mpf(u, c, 1.1e-10_dp, 3.3e-3_dp, 8.7_dp, 0.0_dp, tout, c/2)
   !    call mpf(u, c, 1.15395968E-02_dp, 2.61676386E-03_dp, &
   !             9.23138811_dp, 0.0_dp, tout, c/2)
   !    write (*, *) c, tout
   ! end do

end program example4
