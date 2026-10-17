module bench_models
!! Explicit model functions used by the BLAS benchmark.

   use odrpack_kinds, only: dp
   use iso_c_binding, only: c_ptr
   implicit none

contains

   pure function linear_fn(beta, x) result(f)
      real(dp), intent(in) :: beta(:), x(:)
      real(dp) :: f(size(x))
      f = beta(1) + beta(2)*x
   end function linear_fn

   pure function expo_fn(beta, x) result(f)
      real(dp), intent(in) :: beta(:), x(:)
      real(dp) :: f(size(x))
      f = beta(1)*exp(beta(2)*x) + beta(3)
   end function expo_fn

   pure subroutine fcn_linear( &
      n, m, q, np, ldifx, beta, xplusd, ifixb, ifixx, ideval, f, fjacb, fjacd, istop, data)
   !! User-supplied model evaluation. Derivatives are left to ODRPACK finite differences.
      integer, intent(in) :: n, m, q, np, ldifx, ideval, ifixb(np), ifixx(ldifx, m)
      real(dp), intent(in) :: beta(np), xplusd(n, m)
      real(dp), intent(out) :: f(n, q), fjacb(n, np, q), fjacd(n, m, q)
      integer, intent(out) :: istop
      type(c_ptr), intent(in), value :: data

      istop = 0
      f(:, 1) = linear_fn(beta, xplusd(:, 1))
   end subroutine fcn_linear

   pure subroutine fcn_expo( &
      n, m, q, np, ldifx, beta, xplusd, ifixb, ifixx, ideval, f, fjacb, fjacd, istop, data)
   !! User-supplied model evaluation. Derivatives are left to ODRPACK finite differences.
      integer, intent(in) :: n, m, q, np, ldifx, ideval, ifixb(np), ifixx(ldifx, m)
      real(dp), intent(in) :: beta(np), xplusd(n, m)
      real(dp), intent(out) :: f(n, q), fjacb(n, np, q), fjacd(n, m, q)
      integer, intent(out) :: istop
      type(c_ptr), intent(in), value :: data

      istop = 0
      f(:, 1) = expo_fn(beta, xplusd(:, 1))
   end subroutine fcn_expo

end module bench_models

program bench_odr
!! Timing benchmark for ODR fits. Runs every case once and prints one CSV row per case.
!! Usage: bench_odr <rep>

   use iso_fortran_env, only: int64, output_unit
   use odrpack_kinds, only: dp
   use odrpack, only: odr, odrpack_model
   use bench_models, only: fcn_linear, fcn_expo, linear_fn, expo_fn
   implicit none

   integer, parameter :: sizes(3) = [1000, 10000, 100000]
   integer :: rep, im, is, ib
   character(len=32) :: arg

   call get_command_argument(1, arg)
   read (arg, *) rep

   do im = 1, 2
      do is = 1, size(sizes)
         do ib = 0, 1
            call run_case(im, sizes(is), ib == 1, rep)
         end do
      end do
   end do

contains

   subroutine run_case(model_id, n, bounded, rep)
   !! Generates seeded synthetic data, times one ODR fit and prints one CSV row.
      integer, intent(in) :: model_id, n, rep
      logical, intent(in) :: bounded

      type(odrpack_model) :: model
      integer :: np, q, m, info, nseed
      integer, allocatable :: seed(:)
      real(dp), allocatable :: r(:, :), xt(:), x(:, :), y(:, :)
      real(dp), allocatable :: beta(:), beta0(:), beta_true(:), lower(:), upper(:)
      integer(int64) :: c0, c1, rate
      real(dp) :: elapsed
      character(len=8) :: name
      character(len=128) :: beta_str

      if (model_id == 1) then
         np = 2
         name = 'linear'
         beta_true = [2.0_dp, 0.5_dp]
         beta0 = [1.0_dp, 1.0_dp]
         model%fcn => fcn_linear
      else
         np = 3
         name = 'expo'
         beta_true = [2.0_dp, -0.3_dp, 1.0_dp]
         beta0 = [1.0_dp, -0.1_dp, 0.5_dp]
         model%fcn => fcn_expo
      end if

      m = 1
      q = 1
      allocate (r(n, 3), xt(n), x(n, m), y(n, q), beta(np), lower(np), upper(np))
      lower = -10.0_dp
      upper = 10.0_dp

      ! Seeded data so that every build fits identical problems
      call random_seed(size=nseed)
      allocate (seed(nseed))
      seed = 42
      call random_seed(put=seed)
      call random_number(r)
      xt = 10.0_dp*r(:, 1)
      x(:, 1) = xt + 0.2_dp*(r(:, 2) - 0.5_dp)
      if (model_id == 1) then
         y(:, 1) = linear_fn(beta_true, xt) + 0.2_dp*(r(:, 3) - 0.5_dp)
      else
         y(:, 1) = expo_fn(beta_true, xt) + 0.2_dp*(r(:, 3) - 0.5_dp)
      end if

      beta = beta0
      call system_clock(c0, rate)
      if (bounded) then
         call odr(model, n, m, q, np, beta, y, x, lower=lower, upper=upper, &
                  iprint=0, lunerr=0, lunrpt=0, info=info)
      else
         call odr(model, n, m, q, np, beta, y, x, &
                  iprint=0, lunerr=0, lunrpt=0, info=info)
      end if
      call system_clock(c1)
      elapsed = real(c1 - c0, dp)/real(rate, dp)

      write (beta_str, '(*(es23.15,:,";"))') beta
      write (output_unit, '(a,",",i0,",",a,",",i0,",",es16.8,",",i0,",",a)') &
         trim(name), n, merge('bounded', 'free   ', bounded), rep, elapsed, info, &
         trim(adjustl(beta_str))

   end subroutine run_case

end program bench_odr
