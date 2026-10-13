module test_checked_jacobian_m

   use iso_c_binding, only: c_ptr
   use odrpack_kinds, only: dp
   implicit none

contains

   pure subroutine fcn( &
      n, m, q, np, ldifx, beta, xplusd, ifixb, ifixx, ideval, f, fjacb, fjacd, istop, data)

      integer, intent(in) :: n, m, q, np, ldifx, ideval, ifixb(np), ifixx(ldifx, m)
      real(dp), intent(in) :: beta(np), xplusd(n, m)
      real(dp), intent(out) :: f(n, q), fjacb(n, np, q), fjacd(n, m, q)
      integer, intent(out) :: istop
      type(c_ptr), intent(in), value :: data

      istop = 0
      if (mod(ideval, 10) /= 0) then
         f(:, 1) = beta(1) + beta(2)*xplusd(:, 1) + beta(3)*xplusd(:, 1)**2
      end if
      if (mod(ideval/10, 10) /= 0) then
         fjacb(:, 1, 1) = 1.0_dp
         fjacb(:, 2, 1) = xplusd(:, 1)
         fjacb(:, 3, 1) = xplusd(:, 1)**2
      end if
      if (mod(ideval/100, 10) /= 0) then
         fjacd(:, 1, 1) = beta(2) + 2.0_dp*beta(3)*xplusd(:, 1)
      end if

   end subroutine fcn

end module test_checked_jacobian_m

program test_checked_jacobian
   !! Derivative checking must not change the weighted first step or free columns.

   use odrpack_kinds, only: dp
   use odrpack, only: odr, odrpack_model
   use test_checked_jacobian_m, only: fcn
   implicit none

   integer, parameter :: n = 6, m = 1, q = 1, np = 3
   integer, parameter :: limits(2) = [1, 50]
   real(dp), parameter :: truth(np) = [1.0_dp, 2.0_dp, 0.25_dp]
   real(dp), parameter :: start(np) = [1.1_dp, 1.8_dp, 0.3_dp]
   real(dp), parameter :: tol = 1.0e-10_dp
   type(odrpack_model) :: model
   integer :: i, task, weighted, fixed, limit, mode, ifixb(np), info(2), istop, failures
   real(dp) :: x(n, m), y(n, q), beta(np, 2), delta(n, m, 2), we(n, 1, q), wd(1, 1, m)
   real(dp) :: yest(n, q), fjacb(n, np, q), fjacd(n, m, q), ss(2)

   model%fcn => fcn
   x(:, 1) = [-1.0_dp, -0.6_dp, -0.2_dp, 0.2_dp, 0.6_dp, 1.0_dp]
   call model%fcn(n, m, q, np, 1, truth, x, [1, 1, 1], reshape([1], [1, 1]), &
                  1, y, fjacb, fjacd, istop, model%data)
   failures = 0
   wd = 25.0_dp

   do task = 0, 2, 2
      do weighted = 0, 1
         we = 1.0_dp
         if (weighted /= 0) then
            we(:, 1, 1) = [(100.0_dp*real(i, dp), i=1, n)]
         end if
         do fixed = 0, 1
            ifixb = 1
            ! Fix the first column so the remaining columns must be packed.
            if (fixed /= 0) ifixb(1) = 0
            do limit = 1, size(limits)
               do mode = 1, 2
                  beta(:, mode) = start
                  if (fixed /= 0) beta(1, mode) = truth(1)
                  delta(:, :, mode) = 0.0_dp
                  call odr(model, n, m, q, np, beta(:, mode), y, x, &
                           delta=delta(:, :, mode), we=we, wd=wd, ifixb=ifixb, &
                           job=10*(mode + 1) + task, maxit=limits(limit), iprint=0, info=info(mode))
                  call model%fcn(n, m, q, np, 1, beta(:, mode), x + delta(:, :, mode), &
                                 ifixb, reshape([1], [1, 1]), 1, yest, fjacb, fjacd, istop, model%data)
                  ss(mode) = sum(we(:, 1, 1)*(yest(:, 1) - y(:, 1))**2) &
                             + wd(1, 1, 1)*sum(delta(:, :, mode)**2)
               end do

               if (any(info < 1) .or. any(info > 4) .or. &
                   .not. all(abs(beta(:, 1) - beta(:, 2)) <= tol) .or. &
                   .not. all(abs(delta(:, :, 1) - delta(:, :, 2)) <= tol) .or. &
                   .not. (abs(ss(1) - ss(2)) <= tol*max(1.0_dp, ss(2)))) then
                  write (*, *) 'task, weighted, fixed, maxit, info:', task, weighted, fixed, limits(limit), info
                  write (*, *) 'checked beta:', beta(:, 1)
                  write (*, *) 'unchecked beta:', beta(:, 2)
                  write (*, *) 'sum of squares:', ss
                  failures = failures + 1
               end if
               if (fixed /= 0 .and. any(beta(1, :) /= truth(1))) error stop 'Fixed parameter changed'
               if (limit == 2) then
                  do mode = 1, 2
                     if (info(mode) > 3 .or. .not. all(abs(beta(:, mode) - truth) <= 1.0e-8_dp) .or. &
                         .not. (ss(mode) <= tol)) then
                        write (*, *) 'Solution failed:', task, weighted, fixed, mode, info(mode)
                        failures = failures + 1
                     end if
                  end do
               end if
            end do
         end do
      end do
   end do
   if (failures /= 0) error stop 'Checked-Jacobian regression failed'

end program test_checked_jacobian
