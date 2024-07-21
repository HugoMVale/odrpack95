module example5_model
!! Model for example5.

   use odrpack_kinds, only: wp
   implicit none

contains

   pure subroutine fcn(n, m, np, nq, ldn, ldm, ldnp, beta, xplusd, ifixb, ifixx, &
                       ldifx, ideval, f, fjacb, fjacd, istop)
   !! User-supplied subroutine for evaluating the model.
      integer, intent(in) :: ideval, ldifx, ldm, ldn, ldnp, m, n, np, nq
      integer, intent(in) :: ifixb(np), ifixx(ldifx, m)
      real(kind=wp), intent(in) :: beta(np), xplusd(ldn, m)
      real(kind=wp), intent(out) :: f(ldn, nq), fjacb(ldn, ldnp, nq), fjacd(ldn, ldm, nq)
      integer, intent(out) :: istop

      istop = 0

      ! Calculate model
      if (mod(ideval, 10) .ne. 0) then
         f(:, 1) = beta(1)*exp(beta(2)*xplusd(:, 1))
      end if

      ! Calculate model partials with respect to `beta`
      if (mod(ideval/10, 10) .ne. 0) then
         fjacb(:, 1, 1) = exp(beta(2)*xplusd(:, 1))
         fjacb(:, 2, 1) = beta(1)*xplusd(:, 1)*exp(beta(2)*xplusd(:, 1))
      end if

      ! Calculate model partials with respect to `delta`
      if (mod(ideval/100, 10) .ne. 0) then
         fjacd(:, 1, 1) = beta(1)*beta(2)*exp(beta(2)*xplusd(:, 1))
      end if

   end subroutine fcn

end module example5_model

program example5
!! Explicit ODR job, with parameter bounds, user-supplied derivatives, and output of work
!! arrays.
   use odrpack_kinds, only: wp
   use odrpack, only: odr
   use example5_model, only: fcn
   implicit none

   real(kind=wp), allocatable :: beta(:), lower(:), upper(:), x(:, :), y(:, :)
   real(kind=wp), pointer :: work(:)
   integer, pointer :: iwork(:)
   integer :: np, n, m, nq, job
   !integer :: i,

   np = 2
   n = 4
   m = 1
   nq = 1

   nullify (iwork, work)

   allocate (beta(np), lower(np), upper(np), x(n, m), y(n, nq))

   beta(1:2) = [2.0_wp, 0.5_wp]
   lower(1:2) = [0.0_wp, 0.0_wp]
   upper(1:2) = [10.0_wp, 0.9_wp]
   x(1:4, 1) = [0.982_wp, 1.998_wp, 4.978_wp, 6.01_wp]
   y(1:4, 1) = [2.7_wp, 7.4_wp, 148.0_wp, 403.0_wp]

   job = 20

   call odr(fcn, n, m, np, nq, beta, y, x, job=job, iwork=iwork, work=work, &
            lower=lower, upper=upper)

   ! print *
   ! print *, "iwork"
   ! do i=1, size(iwork)
   !    print *, i, iwork(i)
   ! end do

end program example5
