module odrpack_c
   use, intrinsic :: iso_c_binding, only: c_int, c_double, c_char
   use odrpack, only: odr
   use odrpack_kinds, only: wp
   implicit none
   private

   public :: odr_short_c, dwinf_c, open_file, close_file

   abstract interface
      subroutine fcn_tc(n, m, np, nq, ldn, ldm, ldnp, beta, xplusd, ifixb, ifixx, ldifx, &
                        ideval, f, fjacb, fjacd, istop) bind(C)
         !! User-supplied subroutine for evaluating the model.
         import :: c_int, c_double
         integer(c_int), intent(in) :: n
         integer(c_int), intent(in) :: m
         integer(c_int), intent(in) :: np
         integer(c_int), intent(in) :: nq
         integer(c_int), intent(in) :: ldn
         integer(c_int), intent(in) :: ldm
         integer(c_int), intent(in) :: ldnp
         real(c_double), intent(in) :: beta(np)
         real(c_double), intent(in) :: xplusd(ldn, m)
         integer(c_int), intent(in) :: ifixb(np)
         integer(c_int), intent(in) :: ifixx(ldifx, m)
         integer(c_int), intent(in) :: ldifx
         integer(c_int), intent(in) :: ideval
         real(c_double), intent(out) :: f(ldn, nq)
         real(c_double), intent(out) :: fjacb(ldn, ldnp, nq)
         real(c_double), intent(out) :: fjacd(ldn, ldm, nq)
         integer(c_int), intent(out) :: istop
      end subroutine
   end interface

contains

   subroutine odr_short_c(fcn, n, m, np, nq, beta, y, x, job, lower, upper) bind(C)
      !! "Short" wrapper for the ODR routine including only principal arguments.
      procedure(fcn_tc) :: fcn
         !! User-supplied subroutine for evaluating the model.
      integer(c_int), intent(in) :: n
         !! Number of observations.
      integer(c_int), intent(in) :: m
         !! Number of columns of data in the independent variable.
      integer(c_int), intent(in) :: np
         !! Number of function parameters.
      integer(c_int), intent(in) :: nq
         !! Number of responses per observation.
      real(c_double), intent(inout) :: beta(np)
         !! Function parameters.
      real(c_double), intent(in) :: y(n, nq)
         !! Dependent variable. Unused when the model is implicit.
      real(c_double), intent(in) :: x(n, m)
         !! Explanatory variable.
      integer(c_int), intent(in), optional :: job
         !! Variable controlling problem initialization and computational method.
      real(c_double), intent(in), optional :: lower(np)
         !! Lower bound on `beta`.
      real(kind=wp), intent(in), optional :: upper(np)
         !! Upper bound on `beta`.

      call odr(fcn, n, m, np, nq, beta, y, x, job=job, lower=lower, upper=upper)

   end subroutine odr_short_c

   subroutine odr_c( &
      fcn, &
      n, m, np, nq, &
      beta, &
      y, x, &
      ! delta, &
      ! we, wd, &
      ifixb, & !ifixx, &
      job, ndigit, taufac, &
      sstol, partol, maxit, &
      iprint, lunerr, lunrpt, &
      stpb, & !stpd, &
      sclb, & !scld, &
      info, &
      lower, upper) bind(C)
      !! Driver routine for finding the weighted explicit or implicit orthogonal distance
      !! regression (ODR) or ordinary linear or nonlinear least squares (OLS) solution (long call
      !! statement).
      procedure(fcn_tc) :: fcn
         !! User-supplied subroutine for evaluating the model.
      integer(c_int), intent(in) :: n
         !! Number of observations.
      integer(c_int), intent(in) :: m
         !! Number of columns of data in the independent variable.
      integer(c_int), intent(in) :: np
         !! Number of function parameters.
      integer(c_int), intent(in) :: nq
         !! Number of responses per observation.
      real(c_double), intent(inout) :: beta(np)
         !! Function parameters.
      real(c_double), intent(in) :: y(n, nq)
         !! Dependent variable. Unused when the model is implicit.
      real(c_double), intent(in) :: x(n, m)
         !! Explanatory variable.
      ! real(c_double), intent(inout), optional :: delta(n, m)
      !    !! Initial error in the `x` data.
      !    real(c_double), intent(in), optional :: we(ldwe, ld2we, nq)
      !       !! `epsilon` weights.
      !    real(c_double), intent(in), optional :: wd(ldwd, ld2wd, m)
      !       !! `delta` weights.
      integer(c_int), intent(in), optional :: ifixb(np)
         !! Values designating whether the elements of `beta` are fixed at their input values or not.
      !    integer(c_int), intent(in), optional :: ifixx(ldifx, m)
      !    !! Values designating whether the elements of `x` are fixed at their input values or not.
      integer(c_int), intent(in), optional :: job
         !! Variable controlling problem initialization and computational method.
      integer(c_int), intent(in), optional :: ndigit
         !! Number of accurate digits in the function results, as supplied by the user.
      real(c_double), intent(in), optional :: taufac
         !! Factor used to compute the initial trust region diameter.
      real(c_double), intent(in), optional :: sstol
         !! Sum-of-squares convergence stopping tolerance.
      real(c_double), intent(in), optional :: partol
         !! Parameter convergence stopping tolerance.
      integer(c_int), intent(in), optional :: maxit
         !! Maximum number of iterations allowed.
      integer(c_int), intent(in), optional :: iprint
         !! Print control variable.
      integer(c_int), intent(in), optional :: lunerr
         !! Logical unit number for error messages.
      integer(c_int), intent(in), optional :: lunrpt
         !! Logical unit number for computation reports.
      real(c_double), intent(in), optional :: stpb(np)
         !! Relative step for computing finite difference derivatives with respect to `beta`.
      !    integer(c_int), intent(in) :: ldstpd
      !    real(c_double), intent(in), optional :: stpd(ldstpd, m)
      !       !! Relative step for computing finite difference derivatives with respect to `delta`.
      real(c_double), intent(in), optional :: sclb(np)
         !! Scaling values for `beta`.
      !    real(c_double), intent(in), optional :: scld(ldscld, m)
      !       !! Scaling values for `delta`.
      integer(c_int), intent(out), optional :: info
         !! Variable designating why the computations were stopped.
      real(c_double), intent(in), optional :: lower(np)
         !! Lower bound on `beta`.
      real(c_double), intent(in), optional :: upper(np)
         !! Upper bound on `beta`.

      call odr(fcn, n, m, np, nq, beta, y, x, &
               ! delta=delta, &
               ! we=we, wd=wd, &
               ifixb=ifixb, &
               ! ifixx=ifixx, &
               job=job, ndigit=ndigit, taufac=taufac, &
               sstol=sstol, partol=partol, maxit=maxit, &
               iprint=iprint, lunerr=lunerr, lunrpt=lunrpt, &
               stpb=stpb, &
               !stpd=stpd, &
               sclb=sclb, &
               !scld=scld, &
               info=info, &
               lower=lower, upper=upper)

   end subroutine odr_c

   subroutine open_file(lun, fn, fnlen, ierr) bind(C)
   !! Open a new file associated with a specified logical unit number.
      integer(c_int), intent(inout) :: lun
         !! Logical unit number. If `lun>0`, the user-supplied logical unit number is used.
         !! Otherwise, a new logical unit number is assigned.
      integer(c_int), intent(in) :: fnlen
         !! Length of the string containing the file name, i.e. `strlen(fn)`.
      character(kind=c_char), intent(in) :: fn(fnlen)
         !! String containing the file name.
      integer(c_int), intent(out) :: ierr
         !! Error code returned by `iostat`.

      character(len=fnlen) :: fn_
      integer :: i

      do i = 1, fnlen
         fn_(i:i) = fn(i)
      end do

      if (lun > 0) then
         open (file=fn_, unit=lun, status='new', iostat=ierr)
      else
         open (file=fn_, newunit=lun, status='new', iostat=ierr)
      end if

   end subroutine open_file

   subroutine close_file(lun) bind(C)
   !! Close a file associated with a specified logical unit number.
      integer(c_int), intent(in) :: lun
         !! Logical unit number.
      close (unit=lun)
   end subroutine close_file

   subroutine dwinf_c( &
      n, m, np, nq, ldwe, ld2we, isodr, &
      deltai, epsi, xplusi, fni, sdi, vcvi, &
      rvari, wssi, wssdei, wssepi, rcondi, etai, &
      olmavi, taui, alphai, actrsi, pnormi, rnorsi, prersi, &
      partli, sstoli, taufci, epsmai, &
      beta0i, betaci, betasi, betani, si, ssi, ssfi, qrauxi, ui, &
      fsi, fjacbi, we1i, diffi, &
      deltsi, deltni, ti, tti, omegai, fjacdi, &
      wrk1i, wrk2i, wrk3i, wrk4i, wrk5i, wrk6i, wrk7i, &
      loweri, upperi, &
      lwkmn) bind(C)
   !! Set storage locations within real work space
      integer(c_int), intent(in) :: n
         !! Number of observations.
      integer(c_int), intent(in) :: m
         !! Number of columns of data in the explanatory variable.
      integer(c_int), intent(in) :: np
         !! Number of function parameters.
      integer(c_int), intent(in) :: nq
         !! Number of responses per observation.
      integer(c_int), intent(in) :: ldwe
         !! Leading dimension of array `we`.
      integer(c_int), intent(in) :: ld2we
         !! Second dimension of array `we`.
      integer(c_int), intent(out) :: isodr
         !! Variable designating whether the solution is by ODR (isodr=.true.) or by OLS (isodr=.false.).
      integer(c_int), intent(out) :: deltai
         !! Starting location in array `work` of array `delta.
      integer(c_int), intent(out) :: epsi
         !! Starting location in array `work` of array `eps`.
      integer(c_int), intent(out) :: xplusi
         !! Starting location in array `work` of array `xplusd`.
      integer(c_int), intent(out) :: fni
         !! Starting location in array `work` of array `fn`.
      integer(c_int), intent(out) :: sdi
         !! Starting location in array `work` of array `sd`.
      integer(c_int), intent(out) :: vcvi
         !! Starting location in array `work` of array `vcv`.
      integer(c_int), intent(out) :: rvari
         !! Location in array `work` of variable `rvar`.
      integer(c_int), intent(out) :: wssi
         !! Location in array `work` of variable `wss`.
      integer(c_int), intent(out) :: wssdei
         !! Location in array `work` of variable `wssdel`.
      integer(c_int), intent(out) :: wssepi
         !! Location in array `work` of variable `wsep`.
      integer(c_int), intent(out) :: rcondi
         !! Location in array `work` of variable `rcond`.
      integer(c_int), intent(out) :: etai
         !! Location in array `work` of variable `eta`.
      integer(c_int), intent(out) :: olmavi
         !! Location in array `work` of variable `olmavg`.
      integer(c_int), intent(out) :: taui
         !! Location in array `work` of variable `tau`.
      integer(c_int), intent(out) :: alphai
         !! Location in array `work` of variable `alpha`.
      integer(c_int), intent(out) :: actrsi
         !! Location in array `work` of variable `actrs`.
      integer(c_int), intent(out) :: pnormi
         !! Location in array `work` of variable `pnorm`.
      integer(c_int), intent(out) :: rnorsi
         !! Location in array `work` of variable `rnorms`.
      integer(c_int), intent(out) :: prersi
         !! Location in array `work` of variable `prers`.
      integer(c_int), intent(out) :: partli
         !! Location in array `work` of variable `partol`.
      integer(c_int), intent(out) :: sstoli
         !! Location in array `work` of variable `sstol`.
      integer(c_int), intent(out) :: taufci
         !! Location in array `work` of variable `taufac`.
      integer(c_int), intent(out) :: epsmai
         !! Location in array `work` of variable `epsmac`.
      integer(c_int), intent(out) :: beta0i
         !! Starting location in array `work` of array `beta0`.
      integer(c_int), intent(out) :: betaci
         !! Starting location in array `work` of array `betac`.
      integer(c_int), intent(out) :: betasi
         !! Starting location in array `work` of array `betas`.
      integer(c_int), intent(out) :: betani
         !! Starting location in array `work` of array `betan`.
      integer(c_int), intent(out) :: si
         !! Starting location in array `work` of array `s`.
      integer(c_int), intent(out) :: ssi
         !! Starting location in array `work` of array `ss`.
      integer(c_int), intent(out) :: ssfi
         !! Starting location in array `work` of array `ssf`.
      integer(c_int), intent(out) :: qrauxi
         !! Starting location in array `work` of array `qraux`.
      integer(c_int), intent(out) :: ui
         !! Starting location in array `work` of array `u`.
      integer(c_int), intent(out) :: fsi
         !! Starting location in array `work` of array `fs`.
      integer(c_int), intent(out) :: fjacbi
         !! Starting location in array `work` of array `fjacb`.
      integer(c_int), intent(out) :: we1i
         !! Location in array `work` of variable `we1`.
      integer(c_int), intent(out) :: diffi
         !! Starting location in array `work` of array `diff`.
      integer(c_int), intent(out) :: deltsi
         !! Starting location in array `work` of array `deltas`.
      integer(c_int), intent(out) :: deltni
         !! Starting location in array `work` of array `deltan`.
      integer(c_int), intent(out) :: ti
         !! Starting location in array `work` of array `t`.
      integer(c_int), intent(out) :: tti
         !! Starting location in array `work` of array `tt`.
      integer(c_int), intent(out) :: omegai
         !! Starting location in array `work` of array `omega`.
      integer(c_int), intent(out) :: fjacdi
         !! Starting location in array `work` of array `fjacd`.
      integer(c_int), intent(out) :: wrk1i
         !! Starting location in array `work` of array `wrk1`.
      integer(c_int), intent(out) :: wrk2i
         !! Starting location in array `work` of array `wrk2`.
      integer(c_int), intent(out) :: wrk3i
         !! Starting location in array `work` of array `wrk3`.
      integer(c_int), intent(out) :: wrk4i
         !! Starting location in array `work` of array `wrk4`.
      integer(c_int), intent(out) :: wrk5i
         !! Starting location in array `work` of array `wrk5`.
      integer(c_int), intent(out) :: wrk6i
         !! Starting location in array `work` of array `wrk6`.
      integer(c_int), intent(out) :: wrk7i
         !! Starting location in array `work` of array `wrk7`.
      integer(c_int), intent(out) :: loweri
         !! Starting location in array `work` of array `lower`.
      integer(c_int), intent(out) :: upperi
         !! Starting location in array `work` of array `upper`.
      integer(c_int), intent(out) :: lwkmn
         !! Minimum acceptable length of vector work.

      external :: dwinf

      call dwinf(n, m, np, nq, ldwe, ld2we, isodr, &
                 deltai, epsi, xplusi, fni, sdi, vcvi, &
                 rvari, wssi, wssdei, wssepi, rcondi, etai, &
                 olmavi, taui, alphai, actrsi, pnormi, rnorsi, prersi, &
                 partli, sstoli, taufci, epsmai, &
                 beta0i, betaci, betasi, betani, si, ssi, ssfi, qrauxi, ui, &
                 fsi, fjacbi, we1i, diffi, &
                 deltsi, deltni, ti, tti, omegai, fjacdi, &
                 wrk1i, wrk2i, wrk3i, wrk4i, wrk5i, wrk6i, wrk7i, &
                 loweri, upperi, &
                 lwkmn)

   end subroutine dwinf_c

end module odrpack_c
