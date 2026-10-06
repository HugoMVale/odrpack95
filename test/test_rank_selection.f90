module test_rank_selection_m

   use iso_c_binding, only: c_ptr
   use odrpack_kinds, only: dp
   implicit none

   real(dp) :: signs(3)

contains

   pure subroutine fcn( &
      n, m, q, np, ldifx, beta, xplusd, ifixb, ifixx, ideval, f, fjacb, fjacd, istop, data)

      integer, intent(in) :: n, m, q, np, ldifx, ideval, ifixb(np), ifixx(ldifx, m)
      real(dp), intent(in) :: beta(np), xplusd(n, m)
      real(dp), intent(out) :: f(n, q), fjacb(n, np, q), fjacd(n, m, q)
      integer, intent(out) :: istop
      type(c_ptr), intent(in), value :: data

      real(dp) :: a, b, c, e(n)

      a = signs(1)*beta(1)
      b = signs(2)*beta(2)
      c = signs(3)*beta(3)
      e = exp(-xplusd(:, 1)/b)
      istop = 0

      if (mod(ideval, 10) /= 0) then
         f(:, 1) = a*e + c
      end if
      if (mod(ideval/10, 10) /= 0) then
         fjacb(:, 1, 1) = signs(1)*e
         fjacb(:, 2, 1) = signs(2)*a*e*xplusd(:, 1)/b**2
         fjacb(:, 3, 1) = signs(3)
      end if
      if (mod(ideval/100, 10) /= 0) then
         fjacd(:, 1, 1) = -a*e/b
      end if

   end subroutine fcn

end module test_rank_selection_m

program test_rank_selection
   !! Recover from a nearly rank-deficient start independently of parameter signs.

   use odrpack_kinds, only: dp
   use odrpack, only: odr, odrpack_model
   use test_rank_selection_m, only: fcn, signs
   implicit none

   integer, parameter :: n = 10, m = 1, q = 1, np = 3
   ! OLS with central differences and unchecked analytic derivatives.
   integer, parameter :: jobs(2) = [12, 32]
   real(dp), parameter :: truth(np) = [5000.0_dp, 0.1_dp, 200.0_dp]
   real(dp), parameter :: start(np) = [16000.0_dp, 0.002_dp, 0.0_dp]
   type(odrpack_model) :: model
   integer :: i, j, mask, info, istop
   real(dp) :: x(n, m), y(n, q), beta(np), fitted(np), ss
   real(dp) :: yest(n, q), fjacb(n, np, q), fjacd(n, m, q)

   model%fcn => fcn
   ! Distinct abscissae: 0.01, 0.1, 0.2, ..., 0.9.
   x(1, 1) = 0.01_dp
   do i = 2, n
      x(i, 1) = real(i - 1, dp)/10.0_dp
   end do
   y(:, 1) = truth(1)*exp(-x(:, 1)/truth(2)) + truth(3)

   ! At this start the amplitude and decay-length columns are nearly dependent.
   ! lcstep must select by absolute magnitude, as IDAMAX did, when removing a
   ! direction. Reversing parameter signs leaves the least-squares problem intact.
   do i = 1, size(jobs)
      do mask = 0, 7
         do j = 1, np
            signs(j) = merge(-1.0_dp, 1.0_dp, btest(mask, j - 1))
         end do
         beta = signs*start
         call odr(model, n, m, q, np, beta, y, x, &
                  job=jobs(i), sstol=1.0e-14_dp, maxit=1000, iprint=0, info=info)
         fitted = signs*beta
         call model%fcn(n, m, q, np, 1, beta, x, [1, 1, 1], reshape([1], [1, 1]), &
                        1, yest, fjacb, fjacd, istop, model%data)
         ss = sum((yest - y)**2)

         if (istop /= 0 .or. info < 1 .or. info > 3 .or. &
             .not. all(abs(fitted - truth) <= 1.0e-8_dp*abs(truth)) .or. &
             .not. (ss <= 1.0e-12_dp)) then
            write (*, *) 'job, sign mask, info:', jobs(i), mask, info
            write (*, *) 'parameters:', fitted
            write (*, *) 'sum of squares:', ss
            error stop 'Rank-selection regression failed'
         end if
      end do
   end do

end program test_rank_selection
