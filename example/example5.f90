module example5_model
!! Model for example5.

   use iso_c_binding, only: c_ptr
   use odrpack_kinds, only: dp
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

      istop = 0

      ! Calculate model
      if (mod(ideval, 10) > 0) then
         f(:, 1) = beta(1)*exp(beta(2)*xplusd(:, 1))
      end if

      ! Calculate model partials with respect to `beta`
      if (mod(ideval/10, 10) > 0) then
         fjacb(:, 1, 1) = exp(beta(2)*xplusd(:, 1))
         fjacb(:, 2, 1) = beta(1)*xplusd(:, 1)*exp(beta(2)*xplusd(:, 1))
      end if

      ! Calculate model partials with respect to `delta`
      if (mod(ideval/100, 10) > 0) then
         fjacd(:, 1, 1) = beta(1)*beta(2)*exp(beta(2)*xplusd(:, 1))
      end if

   end subroutine fcn

end module example5_model

program example5
!! Explicit ODR job, with parameter bounds, user-supplied derivatives, and output of work
!! arrays.

   use odrpack_kinds, only: dp
   use odrpack, only: odr, odrpack_model, workspace_dimensions
   use example5_model, only: fcn
   implicit none

   type(odrpack_model) :: model
   real(dp), allocatable :: beta(:), lower(:), upper(:), x(:, :), y(:, :), rwork(:)
   integer, allocatable :: iwork(:)
   integer :: np, n, m, q, job
   integer :: lrwork, liwork !, i

   model%fcn => fcn

   job = 20

   np = 2
   n = 4
   m = 1
   q = 1

   allocate (beta(np), lower(np), upper(np), x(n, m), y(n, q))

   beta(1:2) = [2.0_dp, 0.5_dp]
   lower(1:2) = [0.0_dp, 0.0_dp]
   upper(1:2) = [10.0_dp, 0.9_dp]
   x(1:4, 1) = [0.982_dp, 1.998_dp, 4.978_dp, 6.01_dp]
   y(1:4, 1) = [2.7_dp, 7.4_dp, 148.0_dp, 403.0_dp]

   ! Manual allocation of work arrays
   ! Not required! Just to show it can be done if so desired
   call workspace_dimensions(n, m, q, np, .true., lrwork, liwork)
   allocate (iwork(liwork))
   allocate (rwork(lrwork))

   call odr(model, n, m, q, np, beta, y, x, &
            lower=lower, upper=upper, &
            job=job, iwork=iwork, rwork=rwork)

   ! Remove the comments to print out 'iwork'
   ! print *, "iwork"
   ! do i=1, size(iwork)
   !    print *, i, iwork(i)
   ! end do

end program example5
