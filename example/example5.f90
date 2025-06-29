module example5_model
!! Model for example5.

   use odrpack_kinds, only: wp
   implicit none

contains

   pure subroutine fcn(beta, xplusd, ifixb, ifixx, ideval, f, fjacb, fjacd, istop)
   !! User-supplied subroutine for evaluating the model.

      integer, intent(in) :: ideval, ifixb(:), ifixx(:, :)
      real(kind=wp), intent(in) :: beta(:), xplusd(:, :)
      real(kind=wp), intent(out) :: f(:, :), fjacb(:, :, :), fjacd(:, :, :)
      integer, intent(out) :: istop

      istop = 0

      ! Calculate model
      if (mod(ideval, 10) /= 0) then
         f(:, 1) = beta(1)*exp(beta(2)*xplusd(:, 1))
      end if

      ! Calculate model partials with respect to `beta`
      if (mod(ideval/10, 10) /= 0) then
         fjacb(:, 1, 1) = exp(beta(2)*xplusd(:, 1))
         fjacb(:, 2, 1) = beta(1)*xplusd(:, 1)*exp(beta(2)*xplusd(:, 1))
      end if

      ! Calculate model partials with respect to `delta`
      if (mod(ideval/100, 10) /= 0) then
         fjacd(:, 1, 1) = beta(1)*beta(2)*exp(beta(2)*xplusd(:, 1))
      end if

   end subroutine fcn

end module example5_model

program example5
!! Explicit ODR job, with parameter bounds, user-supplied derivatives, and output of work
!! arrays.

   use odrpack_kinds, only: wp
   use odrpack, only: odr, workspace_dimensions
   use example5_model, only: fcn
   implicit none

   real(kind=wp), allocatable :: beta(:), lower(:), upper(:), x(:, :), y(:, :), rwork(:)
   integer, allocatable :: iwork(:)
   integer :: np, n, m, q, job
   integer :: lrwork, liwork !, i

   job = 20

   np = 2
   n = 4
   m = 1
   q = 1

   allocate (beta(np), lower(np), upper(np), x(n, m), y(n, q))

   beta(1:2) = [2.0_wp, 0.5_wp]
   lower(1:2) = [0.0_wp, 0.0_wp]
   upper(1:2) = [10.0_wp, 0.9_wp]
   x(1:4, 1) = [0.982_wp, 1.998_wp, 4.978_wp, 6.01_wp]
   y(1:4, 1) = [2.7_wp, 7.4_wp, 148.0_wp, 403.0_wp]

   ! Manual allocation of work arrays
   ! Not required! Just to show it can be done if so desired
   call workspace_dimensions(n, m, q, np, .true., lrwork, liwork)
   allocate (iwork(liwork))
   allocate (rwork(lrwork))

   call odr(fcn, n, m, q, np, beta, y, x, &
            lower=lower, upper=upper, &
            job=job, iwork=iwork, rwork=rwork)

   ! Remove the comments to print out 'iwork'
   ! print *, "iwork"
   ! do i=1, size(iwork)
   !    print *, i, iwork(i)
   ! end do

end program example5
