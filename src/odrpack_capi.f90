module odrpack_capi
   !! C-bindings for 'odrpack'.

   use, intrinsic :: iso_c_binding, only: c_bool, c_char, c_double, c_f_pointer, c_int, c_ptr
   use odrpack, only: odr, workspace_dimensions
   use odrpack_core, only: loc_iwork, loc_rwork
   implicit none

   abstract interface
      subroutine fcn_tc( &
         n, m, q, np, ldifx, beta, xplusd, ifixb, ifixx, ideval, f, fjacb, fjacd, istop) bind(C)
      !! User-supplied subroutine for evaluating the model.
         import :: c_int, c_double
         implicit none
         integer(c_int), intent(in) :: n
            !! Number of observations.
         integer(c_int), intent(in) :: m
            !! Number of columns of data in the independent variable.
         integer(c_int), intent(in) :: q
            !! Number of responses per observation.
         integer(c_int), intent(in) :: np
            !! Number of function parameters.
         integer(c_int), intent(in) :: ldifx
            !! Leading dimension of array `ifixx`.
         real(c_double), intent(in) :: beta(np)
            !! Current values of parameters.
         real(c_double), intent(in) :: xplusd(n, m)
            !! Current value of explanatory variable, i.e., `x + delta`.
         integer(c_int), intent(in) :: ifixb(np)
            !! Indicators for "fixing" parameters (`beta`).
         integer(c_int), intent(in) :: ifixx(ldifx, m)
            !! Indicators for "fixing" explanatory variable (`x`).
         integer(c_int), intent(in) :: ideval
            !! Indicator for selecting computation to be performed.
         real(c_double), intent(out) :: f(n, q)
            !! Predicted function values.
         real(c_double), intent(out) :: fjacb(n, np, q)
            !! Jacobian with respect to `beta`.
         real(c_double), intent(out) :: fjacd(n, m, q)
            !! Jacobian with respect to errors `delta`.
         integer(c_int), intent(out) :: istop
            !! Stopping condition, with meaning as follows.
            !!  `0`: Current `beta` and `x + delta` were acceptable and values were computed
            !!       successfully.
            !!  `1`: Current `beta` and `x + delta` are not acceptable; 'odrpack' should select
            !!       values closer to most recently used values if possible.
            !! `-1`: Current `beta` and `x + delta` are not acceptable; 'odrpack' should stop.
      end subroutine
   end interface

   interface
      integer(c_int) function strlen(string) bind(C)
      !! Length of a C-string.
         import :: c_int, c_ptr
         type(c_ptr), intent(in), value :: string
      end function
   end interface

   type, bind(C) :: rworkidx_t
   !! 0-based indices of the variables stored in the real work array.
      integer(c_int) :: delta
      integer(c_int) :: eps
      integer(c_int) :: xplusd
      integer(c_int) :: fn
      integer(c_int) :: sd
      integer(c_int) :: vcv
      integer(c_int) :: rvar
      integer(c_int) :: wss
      integer(c_int) :: wssdel
      integer(c_int) :: wsseps
      integer(c_int) :: rcond
      integer(c_int) :: eta
      integer(c_int) :: olmavg
      integer(c_int) :: tau
      integer(c_int) :: alpha
      integer(c_int) :: actrs
      integer(c_int) :: pnorm
      integer(c_int) :: rnorms
      integer(c_int) :: prers
      integer(c_int) :: partol
      integer(c_int) :: sstol
      integer(c_int) :: taufac
      integer(c_int) :: epsmac
      integer(c_int) :: beta0
      integer(c_int) :: betac
      integer(c_int) :: betas
      integer(c_int) :: betan
      integer(c_int) :: s
      integer(c_int) :: ss
      integer(c_int) :: ssf
      integer(c_int) :: qraux
      integer(c_int) :: u
      integer(c_int) :: fs
      integer(c_int) :: fjacb
      integer(c_int) :: we1
      integer(c_int) :: diff
      integer(c_int) :: deltas
      integer(c_int) :: deltan
      integer(c_int) :: t
      integer(c_int) :: tt
      integer(c_int) :: omega
      integer(c_int) :: fjacd
      integer(c_int) :: wrk1
      integer(c_int) :: wrk2
      integer(c_int) :: wrk3
      integer(c_int) :: wrk4
      integer(c_int) :: wrk5
      integer(c_int) :: wrk6
      integer(c_int) :: wrk7
      integer(c_int) :: lower
      integer(c_int) :: upper
      integer(c_int) :: lrwkmin
   end type rworkidx_t

   type, bind(C) :: iworkidx_t
   !! 0-based indices of the variables stored in the integer work array.
      integer(c_int) :: msgb
      integer(c_int) :: msgd
      integer(c_int) :: ifix2
      integer(c_int) :: istop
      integer(c_int) :: nnzw
      integer(c_int) :: npp
      integer(c_int) :: idf
      integer(c_int) :: job
      integer(c_int) :: iprint
      integer(c_int) :: lunerr
      integer(c_int) :: lunrpt
      integer(c_int) :: nrow
      integer(c_int) :: ntol
      integer(c_int) :: neta
      integer(c_int) :: maxit
      integer(c_int) :: niter
      integer(c_int) :: nfev
      integer(c_int) :: njev
      integer(c_int) :: int2
      integer(c_int) :: irank
      integer(c_int) :: ldtt
      integer(c_int) :: bound
      integer(c_int) :: liwkmin
   end type iworkidx_t

contains

   subroutine odr_short_c( &
      fcn, &
      n, m, q, np, &
      beta, y, x, &
      delta, &
      lower, upper, &
      job) bind(C)
   !! "Short-call" wrapper for the `odr` routine including very few optional arguments.

      procedure(fcn_tc) :: fcn
         !! User-supplied subroutine for evaluating the model.
      integer(c_int), intent(in) :: n
         !! Number of observations.
      integer(c_int), intent(in) :: m
         !! Number of columns of data in the independent variable.
      integer(c_int), intent(in) :: q
         !! Number of responses per observation.
      integer(c_int), intent(in) :: np
         !! Number of function parameters.
      real(c_double), intent(inout) :: beta(np)
         !! Function parameters.
      real(c_double), intent(in) :: y(n, q)
         !! Dependent variable. Unused when the model is implicit.
      real(c_double), intent(in) :: x(n, m)
         !! Explanatory variable.
      real(c_double), intent(inout), optional :: delta(n, m)
         !! Error in the `x` data. Initial guess on input and estimated value on output.
      real(c_double), intent(in), optional :: lower(np)
         !! Lower bound on `beta`.
      real(c_double), intent(in), optional :: upper(np)
         !! Upper bound on `beta`.
      integer(c_int), intent(in), optional :: job
         !! Variable controlling initialization and computational method.

      call odr(fcn_, n, m, q, np, beta, y, x, &
               delta=delta, &
               lower=lower, upper=upper, &
               job=job)

   contains

      subroutine fcn_(beta, xplusd, ifixb, ifixx, ideval, f, fjacb, fjacd, istop)
         real(c_double), intent(in) :: beta(:)
         real(c_double), intent(in) :: xplusd(:, :)
         integer(c_int), intent(in) :: ifixb(:)
         integer(c_int), intent(in) :: ifixx(:, :)
         integer(c_int), intent(in) :: ideval
         real(c_double), intent(out) :: f(:, :)
         real(c_double), intent(out) :: fjacb(:, :, :)
         real(c_double), intent(out) :: fjacd(:, :, :)
         integer(c_int), intent(out) :: istop

         call fcn(size(xplusd, 1), size(xplusd, 2), size(f, 2), size(beta), size(ifixx, 1), &
                  beta, xplusd, ifixb, ifixx, ideval, f, fjacb, fjacd, istop)

      end subroutine fcn_

   end subroutine odr_short_c

   subroutine odr_medium_c( &
      fcn, &
      n, m, q, np, &
      ldwe, ld2we, &
      ldwd, ld2wd, &
      ldifx, &
      beta, y, x, &
      we, wd, &
      ifixb, ifixx, &
      delta, &
      lower, upper, &
      job, iprint, lunerr, lunrpt, info) bind(C)
   !! "Medium-call" wrapper for the `odr` routine including most commonly used optional arguments.

      procedure(fcn_tc) :: fcn
         !! User-supplied subroutine for evaluating the model.
      integer(c_int), intent(in) :: n
         !! Number of observations.
      integer(c_int), intent(in) :: m
         !! Number of columns of data in the independent variable.
      integer(c_int), intent(in) :: q
         !! Number of responses per observation.
      integer(c_int), intent(in) :: np
         !! Number of function parameters.
      integer(c_int), intent(in) :: ldwe
         !! Leading dimension of array `we`, `ldwe ∈ {1, n}`.
      integer(c_int), intent(in) :: ld2we
         !! Second dimension of array `we`, `ld2we ∈ {1, q}`.
      integer(c_int), intent(in) :: ldwd
         !! Leading dimension of array `wd`, `ldwd ∈ {1, n}`.
      integer(c_int), intent(in) :: ld2wd
         !! Second dimension of array `wd`, `ld2wd ∈ {1, m}`.
      integer(c_int), intent(in) :: ldifx
         !! Leading dimension of array `ifixx`, `ldifx ∈ {1, n}`.
      real(c_double), intent(inout) :: beta(np)
         !! Function parameters.
      real(c_double), intent(in) :: y(n, q)
         !! Dependent variable. Unused when the model is implicit.
      real(c_double), intent(in) :: x(n, m)
         !! Explanatory variable.
      real(c_double), intent(in), optional :: we(ldwe, ld2we, q)
         !! `epsilon` weights.
      real(c_double), intent(in), optional :: wd(ldwd, ld2wd, m)
         !! `delta` weights.
      integer(c_int), intent(in), optional :: ifixb(np)
         !! Values designating whether the elements of `beta` are fixed at their input values or not.
      integer(c_int), intent(in), optional :: ifixx(ldifx, m)
         !! Values designating whether the elements of `x` are fixed at their input values or not.
      real(c_double), intent(inout), optional :: delta(n, m)
         !! Error in the `x` data. Initial guess on input and estimated value on output.
      real(c_double), intent(in), optional :: lower(np)
         !! Lower bound on `beta`.
      real(c_double), intent(in), optional :: upper(np)
         !! Upper bound on `beta`.
      integer(c_int), intent(in), optional :: job
         !! Variable controlling initialization and computational method.
      integer(c_int), intent(in), optional :: iprint
         !! Print control variable.
      integer(c_int), intent(in), optional :: lunerr
         !! Logical unit number for error messages.
         !!  `0`: No output.
         !!  `6`: Output to standard output.
         !!  `k /= 0,6`: Output to logical unit `k`.
      integer(c_int), intent(in), optional :: lunrpt
         !! Logical unit number for computation reports.
         !!  `0`: No output.
         !!  `6`: Output to standard output.
         !!  `k /= 0,6`: Output to logical unit `k`.
      integer(c_int), intent(out), optional :: info
         !! Logical unit number for computation reports.

      call odr(fcn_, n, m, q, np, beta, y, x, &
               we=we, wd=wd, &
               ifixb=ifixb, ifixx=ifixx, &
               delta=delta, &
               lower=lower, upper=upper, &
               job=job, &
               iprint=iprint, lunerr=lunerr, lunrpt=lunrpt, &
               info=info)
   contains

      subroutine fcn_(beta, xplusd, ifixb, ifixx, ideval, f, fjacb, fjacd, istop)
         real(c_double), intent(in) :: beta(:)
         real(c_double), intent(in) :: xplusd(:, :)
         integer(c_int), intent(in) :: ifixb(:)
         integer(c_int), intent(in) :: ifixx(:, :)
         integer(c_int), intent(in) :: ideval
         real(c_double), intent(out) :: f(:, :)
         real(c_double), intent(out) :: fjacb(:, :, :)
         real(c_double), intent(out) :: fjacd(:, :, :)
         integer(c_int), intent(out) :: istop

         call fcn(size(xplusd, 1), size(xplusd, 2), size(f, 2), size(beta), size(ifixx, 1), &
                  beta, xplusd, ifixb, ifixx, ideval, f, fjacb, fjacd, istop)

      end subroutine fcn_

   end subroutine odr_medium_c

   subroutine odr_long_c( &
      fcn, &
      n, m, q, np, &
      ldwe, ld2we, &
      ldwd, ld2wd, &
      ldifx, &
      ldstpd, ldscld, &
      lrwork, liwork, &
      beta, y, x, &
      we, wd, &
      ifixb, ifixx, &
      stpb, stpd, &
      sclb, scld, &
      delta, &
      lower, upper, &
      rwork, iwork, &
      job, ndigit, taufac, &
      sstol, partol, maxit, &
      iprint, lunerr, lunrpt, &
      info) bind(C)
   !! "Long-call" wrapper for the `odr` routine including all optional arguments.

      procedure(fcn_tc) :: fcn
         !! User-supplied subroutine for evaluating the model.
      integer(c_int), intent(in) :: n
         !! Number of observations.
      integer(c_int), intent(in) :: m
         !! Number of columns of data in the independent variable.
      integer(c_int), intent(in) :: q
         !! Number of responses per observation.
      integer(c_int), intent(in) :: np
         !! Number of function parameters.
      integer(c_int), intent(in) :: ldwe
         !! Leading dimension of array `we`, `ldwe ∈ {1, n}`.
      integer(c_int), intent(in) :: ld2we
         !! Second dimension of array `we`, `ld2we ∈ {1, q}`.
      integer(c_int), intent(in) :: ldwd
         !! Leading dimension of array `wd`, `ldwd ∈ {1, n}`.
      integer(c_int), intent(in) :: ld2wd
         !! Second dimension of array `wd`, `ld2wd ∈ {1, m}`.
      integer(c_int), intent(in) :: ldifx
         !! Leading dimension of array `ifixx`, `ldifx ∈ {1, n}`.
      integer(c_int), intent(in) :: ldstpd
         !! Leading dimension of array `stpd`, `ldstpd ∈ {1, n}`.
      integer(c_int), intent(in) :: ldscld
         !! Leading dimension of array `scld`, `ldscld ∈ {1, n}`.
      integer(c_int), intent(in) :: lrwork
         !! Length of array `rwork`.
      integer(c_int), intent(in) :: liwork
         !! Length of array `iwork`.
      real(c_double), intent(inout) :: beta(np)
         !! Function parameters.
      real(c_double), intent(in) :: y(n, q)
         !! Dependent variable. Unused when the model is implicit.
      real(c_double), intent(in) :: x(n, m)
         !! Explanatory variable.
      real(c_double), intent(in), optional :: we(ldwe, ld2we, q)
         !! `epsilon` weights.
      real(c_double), intent(in), optional :: wd(ldwd, ld2wd, m)
         !! `delta` weights.
      integer(c_int), intent(in), optional :: ifixb(np)
         !! Values designating whether the elements of `beta` are fixed at their input values or not.
      integer(c_int), intent(in), optional :: ifixx(ldifx, m)
         !! Values designating whether the elements of `x` are fixed at their input values or not.
      real(c_double), intent(in), optional :: stpb(np)
         !! Relative step for computing finite difference derivatives with respect to `beta`.
      real(c_double), intent(in), optional :: stpd(ldstpd, m)
         !! Relative step for computing finite difference derivatives with respect to `delta`.
      real(c_double), intent(in), optional :: sclb(np)
         !! Scaling values for `beta`.
      real(c_double), intent(in), optional :: scld(ldscld, m)
         !! Scaling values for `delta`.
      real(c_double), intent(inout), optional :: delta(n, m)
         !! Error in the `x` data. `Shape: (n, m)`. Initial guess on input and estimated value
         !! on output.
      real(c_double), intent(in), optional :: lower(np)
         !! Lower bound on `beta`.
      real(c_double), intent(in), optional :: upper(np)
         !! Upper bound on `beta`.
      real(c_double), intent(inout), optional :: rwork(lrwork)
         !! Real work space.
      integer(c_int), intent(inout), optional :: iwork(liwork)
         !! Integer work space.
      integer(c_int), intent(in), optional :: job
         !! Variable controlling initialization and computational method.
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
         !!  `0`: No output.
         !!  `6`: Output to standard output.
         !!  `k /= 0,6`: Output to logical unit `k`.
      integer(c_int), intent(in), optional :: lunrpt
         !! Logical unit number for computation reports.
         !!  `0`: No output.
         !!  `6`: Output to standard output.
         !!  `k /= 0,6`: Output to logical unit `k`.
      integer(c_int), intent(out), optional :: info
         !! Variable designating why the computations were stopped.

      call odr(fcn_, n, m, q, np, beta, y, x, &
               we=we, wd=wd, &
               ifixb=ifixb, ifixx=ifixx, &
               delta=delta, &
               lower=lower, upper=upper, &
               job=job, &
               ndigit=ndigit, taufac=taufac, sstol=sstol, partol=partol, maxit=maxit, &
               iprint=iprint, lunerr=lunerr, lunrpt=lunrpt, &
               stpb=stpb, stpd=stpd, &
               sclb=sclb, scld=scld, &
               info=info, &
               rwork=rwork, iwork=iwork)

   contains

      subroutine fcn_(beta, xplusd, ifixb, ifixx, ideval, f, fjacb, fjacd, istop)
         real(c_double), intent(in) :: beta(:)
         real(c_double), intent(in) :: xplusd(:, :)
         integer(c_int), intent(in) :: ifixb(:)
         integer(c_int), intent(in) :: ifixx(:, :)
         integer(c_int), intent(in) :: ideval
         real(c_double), intent(out) :: f(:, :)
         real(c_double), intent(out) :: fjacb(:, :, :)
         real(c_double), intent(out) :: fjacd(:, :, :)
         integer(c_int), intent(out) :: istop

         call fcn(size(xplusd, 1), size(xplusd, 2), size(f, 2), size(beta), size(ifixx, 1), &
                  beta, xplusd, ifixb, ifixx, ideval, f, fjacb, fjacd, istop)

      end subroutine fcn_

   end subroutine odr_long_c

   subroutine open_file(filename_cptr, lun, ierr) bind(C)
   !! Open a new file associated with a specified logical unit number.

      type(c_ptr), intent(in), value :: filename_cptr
         !! C-string containing the file name.
      integer(c_int), intent(inout) :: lun
         !! Logical unit number. If `lun > 0`, the user-supplied logical unit number is used.
         !! Otherwise, a new logical unit number is assigned.
      integer(c_int), intent(out) :: ierr
         !! Error code returned by `iostat`.

      character(kind=c_char), pointer :: filename_fptr(:)
      character(len=:), allocatable :: filename
      character(len=256) :: errmsg
      integer :: length, i

      length = strlen(filename_cptr)
      allocate (filename_fptr(length))
      call c_f_pointer(cptr=filename_cptr, fptr=filename_fptr, shape=[length])

      allocate (character(len=length) :: filename)
      do i = 1, length
         filename(i:i) = filename_fptr(i)
      end do

      if (lun > 0) then
         open (file=filename, unit=lun, status='replace', iostat=ierr, iomsg=errmsg)
      else
         open (file=filename, newunit=lun, status='replace', iostat=ierr, iomsg=errmsg)
      end if

      if (ierr /= 0) then
         print *, "I/O error: ", trim(errmsg)
      end if

   end subroutine open_file

   subroutine close_file(lun, ierr) bind(C)
   !! Close a file associated with a specified logical unit number.

      integer(c_int), intent(in) :: lun
         !! Logical unit number.
      integer(c_int), intent(out) :: ierr
         !! Error code returned by `iostat`.

      close (unit=lun, iostat=ierr)

   end subroutine close_file

   pure subroutine loc_iwork_c(m, q, np, iwi) bind(C)
   !! Get storage locations within integer work space.

      integer(c_int), intent(in) :: m
         !! Number of columns of data in the independent variable.
      integer(c_int), intent(in) :: q
         !! Number of responses per observation.
      integer(c_int), intent(in) :: np
         !! Number of function parameters.
      type(iworkidx_t), intent(out) :: iwi
         !! 0-based indexes of integer work array.

      integer :: msgbi, msgdi, ifix2i, istopi, nnzwi, nppi, idfi, jobi, iprinti, lunerri, &
                 lunrpti, nrowi, ntoli, netai, maxiti, niteri, nfevi, njevi, int2i, iranki, &
                 ldtti, boundi, liwkmin

      call loc_iwork(m, q, np, &
                     msgbi, msgdi, ifix2i, istopi, &
                     nnzwi, nppi, idfi, &
                     jobi, iprinti, lunerri, lunrpti, &
                     nrowi, ntoli, netai, &
                     maxiti, niteri, nfevi, njevi, int2i, iranki, ldtti, &
                     boundi, &
                     liwkmin)

      iwi%msgb = msgbi - 1
      iwi%msgd = msgdi - 1
      iwi%ifix2 = ifix2i - 1
      iwi%istop = istopi - 1
      iwi%nnzw = nnzwi - 1
      iwi%npp = nppi - 1
      iwi%idf = idfi - 1
      iwi%job = jobi - 1
      iwi%iprint = iprinti - 1
      iwi%lunerr = lunerri - 1
      iwi%lunrpt = lunrpti - 1
      iwi%nrow = nrowi - 1
      iwi%ntol = ntoli - 1
      iwi%neta = netai - 1
      iwi%maxit = maxiti - 1
      iwi%niter = niteri - 1
      iwi%nfev = nfevi - 1
      iwi%njev = njevi - 1
      iwi%int2 = int2i - 1
      iwi%irank = iranki - 1
      iwi%ldtt = ldtti - 1
      iwi%bound = boundi - 1
      iwi%liwkmin = liwkmin

   end subroutine loc_iwork_c

   pure subroutine loc_rwork_c(n, m, q, np, ldwe, ld2we, isodr, rwi) bind(C)
   !! Get storage locations within real work space.

      integer(c_int), intent(in) :: n
         !! Number of observations.
      integer(c_int), intent(in) :: m
         !! Number of columns of data in the explanatory variable.
      integer(c_int), intent(in) :: q
         !! Number of responses per observation.
      integer(c_int), intent(in) :: np
         !! Number of function parameters.
      integer(c_int), intent(in) :: ldwe
         !! Leading dimension of array `we`.
      integer(c_int), intent(in) :: ld2we
         !! Second dimension of array `we`.
      logical(c_bool), intent(in) :: isodr
         !! Variable designating whether the solution is by ODR (`.true.`) or by OLS (`.false.`).
      type(rworkidx_t), intent(out) :: rwi
         !! 0-based indexes of real work array.

      integer :: deltai, epsi, xplusdi, fni, sdi, vcvi, rvari, wssi, wssdeli, wssepsi, &
                 rcondi, etai, olmavgi, taui, alphai, actrsi, pnormi, rnormsi, prersi, partoli, &
                 sstoli, taufaci, epsmaci, beta0i, betaci, betasi, betani, si, ssi, ssfi, &
                 qrauxi, ui, fsi, fjacbi, we1i, diffi, deltasi, deltani, ti, tti, omegai, &
                 fjacdi, wrk1i, wrk2i, wrk3i, wrk4i, wrk5i, wrk6i, wrk7i, loweri, upperi, &
                 lrwkmin

      call loc_rwork(n, m, q, np, ldwe, ld2we, logical(isodr, kind=kind(.true.)), &
                     deltai, epsi, xplusdi, fni, sdi, vcvi, &
                     rvari, wssi, wssdeli, wssepsi, rcondi, etai, &
                     olmavgi, taui, alphai, actrsi, pnormi, rnormsi, prersi, &
                     partoli, sstoli, taufaci, epsmaci, &
                     beta0i, betaci, betasi, betani, si, ssi, ssfi, qrauxi, ui, &
                     fsi, fjacbi, we1i, diffi, &
                     deltasi, deltani, ti, tti, omegai, fjacdi, &
                     wrk1i, wrk2i, wrk3i, wrk4i, wrk5i, wrk6i, wrk7i, &
                     loweri, upperi, &
                     lrwkmin)

      rwi%delta = deltai - 1
      rwi%eps = epsi - 1
      rwi%xplusd = xplusdi - 1
      rwi%fn = fni - 1
      rwi%sd = sdi - 1
      rwi%vcv = vcvi - 1
      rwi%rvar = rvari - 1
      rwi%wss = wssi - 1
      rwi%wssdel = wssdeli - 1
      rwi%wsseps = wssepsi - 1
      rwi%rcond = rcondi - 1
      rwi%eta = etai - 1
      rwi%olmavg = olmavgi - 1
      rwi%tau = taui - 1
      rwi%alpha = alphai - 1
      rwi%actrs = actrsi - 1
      rwi%pnorm = pnormi - 1
      rwi%rnorms = rnormsi - 1
      rwi%prers = prersi - 1
      rwi%partol = partoli - 1
      rwi%sstol = sstoli - 1
      rwi%taufac = taufaci - 1
      rwi%epsmac = epsmaci - 1
      rwi%beta0 = beta0i - 1
      rwi%betac = betaci - 1
      rwi%betas = betasi - 1
      rwi%betan = betani - 1
      rwi%s = si - 1
      rwi%ss = ssi - 1
      rwi%ssf = ssfi - 1
      rwi%qraux = qrauxi - 1
      rwi%u = ui - 1
      rwi%fs = fsi - 1
      rwi%fjacb = fjacbi - 1
      rwi%we1 = we1i - 1
      rwi%diff = diffi - 1
      rwi%deltas = deltasi - 1
      rwi%deltan = deltani - 1
      rwi%t = ti - 1
      rwi%tt = tti - 1
      rwi%omega = omegai - 1
      rwi%fjacd = fjacdi - 1
      rwi%wrk1 = wrk1i - 1
      rwi%wrk2 = wrk2i - 1
      rwi%wrk3 = wrk3i - 1
      rwi%wrk4 = wrk4i - 1
      rwi%wrk5 = wrk5i - 1
      rwi%wrk6 = wrk6i - 1
      rwi%wrk7 = wrk7i - 1
      rwi%lower = loweri - 1
      rwi%upper = upperi - 1
      rwi%lrwkmin = lrwkmin

   end subroutine loc_rwork_c

   pure subroutine workspace_dimensions_c(n, m, q, np, isodr, lrwork, liwork) bind(C)
   !! Calculate the dimensions of the workspace arrays.

      integer(c_int), intent(in) :: n
         !! Number of observations.
      integer(c_int), intent(in) :: m
         !! Number of columns of data in the independent variable.
      integer(c_int), intent(in) :: q
         !! Number of responses per observation.
      integer(c_int), intent(in) :: np
         !! Number of function parameters.
      logical(c_bool), intent(in) :: isodr
         !! Variable designating whether the solution is by ODR (`.true.`) or by OLS (`.false.`).
      integer(c_int), intent(out) :: lrwork
         !! Length of real `rwork` array.
      integer(c_int), intent(out) :: liwork
         !! Length of integer `iwork` array.

      call workspace_dimensions(n, m, q, np, logical(isodr, kind=kind(.true.)), &
                                lrwork, liwork)

   end subroutine workspace_dimensions_c

end module odrpack_capi
