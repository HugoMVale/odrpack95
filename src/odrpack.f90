module odrpack

   use odrpack_kinds, only: wp
   implicit none
   private

   public :: odr, tempret

   ! A temporary work array for holding return values before copying to a lower rank array.
   real(kind=wp), allocatable :: tempret(:, :)

   abstract interface
      subroutine fcn_t(n, m, np, nq, ldn, ldm, ldnp, beta, xplusd, ifixb, ifixx, ldifx, &
                       ideval, f, fjacb, fjacd, istop)
      !! User-supplied subroutine for evaluating the model.
         import :: wp
         integer, intent(in) :: n
            !! Number of observations.
         integer, intent(in) :: m
            !! Number of columns of data in the independent variable.
         integer, intent(in) :: np
            !! Number of function parameters.
         integer, intent(in) :: nq
            !! Number of responses per observation.
         integer, intent(in) :: ldn
            !! Leading dimension declarator equal or exceeding `n`.
         integer, intent(in) :: ldm
            !! Leading dimension declarator equal or exceeding `m`.
         integer, intent(in) :: ldnp
            !! Leading dimension declarator equal or exceeding `np`.
         real(kind=wp), intent(in) :: beta(np)
            !! Current values of parameters.
         real(kind=wp), intent(in) :: xplusd(ldn, m)
            !! Current value of explanatory variable, i.e., `x + delta`.
         integer, intent(in) :: ifixb(np)
            !! Indicators for "fixing" parameters (`beta`).
         integer, intent(in) :: ifixx(ldifx, m)
            !!  Indicators for "fixing" explanatory variable (`x`).	
         integer, intent(in) :: ldifx
            !! Leading dimension of array `ifixx`.
         integer, intent(in) :: ideval
            !! Indicator for selecting computation to be performed.
         real(kind=wp), intent(out) :: f(ldn, nq)
            !! Predicted function values.
         real(kind=wp), intent(out) :: fjacb(ldn, ldnp, nq)
            !! Jacobian with respect to `beta`.
         real(kind=wp), intent(out) :: fjacd(ldn, ldm, nq)
            !! Jacobian with respect to errors `delta`.
         integer, intent(out) :: istop
            !! Stopping condition, with meaning as follows. 0 means current `beta` and
            !! `x+delta` were acceptable and values were computed successfully. 1 means current
            !! `beta` and `x+delta` are not acceptable;  ODRPACK95 should select values closer
            !! to most recently used values if possible. -1 means current `beta` and `x+delta`
            !! are not acceptable; ODRPACK95 should stop.
      end subroutine fcn_t
   end interface

contains

   subroutine odr &
      (fcn, &
       n, m, np, nq, &
       beta, &
       y, x, &
       delta, &
       we, wd, &
       ifixb, ifixx, &
       job, ndigit, taufac, &
       sstol, partol, maxit, &
       iprint, lunerr, lunrpt, &
       stpb, stpd, &
       sclb, scld, &
       work, iwork, &
       info, &
       lower, upper)
   !! Driver routine for finding the weighted explicit or implicit orthogonal distance
   !! regression (ODR) or ordinary linear or nonlinear least squares (OLS) solution (long call
   !! statement).
      ! Date Written:   860529   (YYMMDD)
      ! Revision Date:  20040301 (YYYYMMDD)
      ! Category No.:  G2E,I1B1
      ! Keywords:  Orthogonal distance regression,
      !   Nonlinear least squares,
      !   Measurement error models,
      !   Errors in variables
      ! Authors:
      !   Boggs, Paul T.
      !     Applied and Computational Mathematics Division
      !     National Institute of Standards and Technology
      !     Gaithersburg, MD 20899
      !   Byrd, Richard H.
      !     Department of Computer Science
      !     University of Colorado, Boulder, CO 80309
      !   Rogers, Janet E.
      !     Applied and Computational Mathematics Division
      !     National Institute of Standards and Technology
      !     Boulder, CO 80303-3328
      !   Schnabel, Robert B.
      !     Department of Computer Science
      !     University of Colorado, Boulder, CO 80309
      !     and
      !     Applied and Computational Mathematics Division
      !     National Institute of Standards and Technology
      !     Boulder, CO 80303-3328
      ! Purpose:  REAL (KIND=wp) driver routine for finding the weighted explicit or implicit
      !   orthogonal distance regression (ODR) or ordinary linear or nonlinear least squares (OLS)
      !   solution (long call statement)
      ! Description:
      !   For details, see ODRPACK95 User's Reference Guide.
      ! References:
      !   Boggs, P. T., R. H. Byrd, J. R. Donaldson, and R. B. Schnabel (1989),
      !     "Algorithm 676 --- ODRPACK: Software for Weighted Orthogonal Distance Regression,"
      !     ACM Trans. Math. Software., 15(4):348-364.
      !       Boggs, P. T., R. H. Byrd, J. E. Rogers, and
      !   R. B. Schnabel (1992),
      !     "User's Reference Guide for ODRPACK Version 2.01,
      !     Software for Weighted Orthogonal Distance Regression,"
      !     National Institute of Standards and Technology
      !     Internal Report Number 92-4834.
      !  Boggs, P. T., R. H. Byrd, and R. B. Schnabel (1987),
      !    "A Stable and Efficient Algorithm for Nonlinear Orthogonal Distance Regression,"
      !    SIAM J. Sci. Stat. Comput., 8(6):1052-1078.

      use odrpack_kinds, only: wp, negone, ZERO

      procedure(fcn_t) :: fcn
         !! User-supplied subroutine for evaluating the model.
      integer, intent(in) :: n
         !! Number of observations.
      integer, intent(in) :: m
         !! Number of columns of data in the independent variable.
      integer, intent(in) :: np
         !! Number of function parameters.
      integer, intent(in) :: nq
         !! Number of responses per observation.
      real(kind=wp), intent(inout) :: beta(:)
         !! Function parameters. `Shape: (np)`.
      real(kind=wp), intent(in) :: y(:, :)
         !! Dependent variable. `Shape: (n, nq)`. Unused when the model is implicit.
      real(kind=wp), intent(in) :: x(:, :)
         !! Explanatory variable. `Shape: (n, m)`.
      real(kind=wp), intent(inout), allocatable, optional :: delta(:, :)
         !! Initial error in the `x` data. `Shape: (n, m)`.
      real(kind=wp), intent(in), optional :: we(:, :, :)
         !! `epsilon` weights. `Shape: (1<=ldwe<=n, 1<=ld2we<=nq, nq)`. See p. 25.
      real(kind=wp), intent(in), optional :: wd(:, :, :)
         !! `delta` weights. `Shape: (1<=ldwd<=n, 1<=ld2wd<=m, m)`. See p. 26.
      integer, intent(in), optional :: ifixb(:)
         !! Values designating whether the elements of `beta` are fixed at their input values or not. `Shape: (np)`.
      integer, intent(in), optional :: ifixx(:, :)
         !! Values designating whether the elements of `x` are fixed at their input values or not. `Shape: (1<=ldifx<=n, m)`. See p. 27.
      integer, intent(in), optional :: job
         !! Variable controlling problem initialization and computational method.
      integer, intent(in), optional :: ndigit
         !! Number of accurate digits in the function results, as supplied by the user.
      real(kind=wp), intent(in), optional :: taufac
         !! Factor used to compute the initial trust region diameter.
      real(kind=wp), intent(in), optional :: sstol
         !! Sum-of-squares convergence stopping tolerance.
      real(kind=wp), intent(in), optional :: partol
      !! Parameter convergence stopping tolerance.
      integer, intent(in), optional :: maxit
         !! Maximum number of iterations allowed.
      integer, intent(in), optional :: iprint
         !! Print control variable.
      integer, intent(in), optional :: lunerr
         !! Logical unit number for error messages.
      integer, intent(in), optional :: lunrpt
         !! Logical unit number for computation reports.
      real(kind=wp), intent(in), optional :: stpb(:)
         !! Relative step for computing finite difference derivatives with respect to `beta`. `Shape: (np)`.
      real(kind=wp), intent(in), optional :: stpd(:, :)
         !! Relative step for computing finite difference derivatives with respect to `delta`. `Shape: (1<=ldstpd<=n, m)`. See p. 31.
      real(kind=wp), intent(in), optional :: sclb(:)
      !! Scaling values for `beta`. `Shape: (np)`.
      real(kind=wp), intent(in), optional :: scld(:, :)
         !! Scaling values for `delta`. `Shape: (1<=ldscld<=n, m)`. See p. 32.
      real(kind=wp), intent(inout), pointer, optional :: work(:)
         !! Real work space.
      integer, intent(out), pointer, optional :: iwork(:)
         !! Integer work space.
      integer, intent(out), optional :: info
         !! Variable designating why the computations were stopped.
      real(kind=wp), intent(in), optional :: lower(:)
         !! Lower bound on `beta`. `Shape: (np)`.
      real(kind=wp), intent(in), optional :: upper(:)
         !! Upper bound on `beta`. `Shape: (np)`.

      ! Local variables
      ! TODO: remove save? replace pointer with allocatable?
      integer :: ldwe, ld2we, ldwd, ld2wd, ldifx, ldscld, ldstpd, ljob, lndigit, lmaxit, &
                 liprint, llunerr, llunrpt, linfo, lenwork, leniwork, linfo1, linfo2, linfo3, &
                 linfo4, linfo5
      integer :: lifixb(np), lifixx(n, m)
      real(kind=wp) :: ltaufac, lsstol, lpartol
      real(kind=wp) :: llower(np), lwe(n, nq, nq), lwd(n, m, m), lstpb(np), lstpd(n, m), &
                       lsclb(np), lscld(n, m), lupper(np), wd1(1, 1, 1)

      real(kind=wp), pointer, save :: lwork(:)
      integer, pointer, save :: liwork(:)  
      logical :: head

      ! TODO: interface or inside module?
      external :: dodcnt

      ! Set LINFO to zero indicating no errors have been found thus far
      linfo = 0
      linfo1 = 0
      linfo2 = 0
      linfo3 = 0
      linfo4 = 0
      linfo5 = 0

      ! Set all scalar variable defaults except JOB
      ldwe = 1
      ld2we = 1
      ldwd = 1
      ld2wd = 1
      ldifx = 1
      ldscld = 1
      ldstpd = 1
      liprint = -1
      llunerr = -1
      llunrpt = -1
      lmaxit = -1
      lndigit = -1
      lpartol = negone
      lsstol = negone
      ltaufac = negone
      head = .true.

      !  Check for the option arguments for printing (so error messages can be
      !  printed appropriately from here on out
      if (present(iprint)) then
         liprint = iprint
      end if

      if (present(lunrpt)) then
         llunrpt = lunrpt
      end if
      if (llunrpt .lt. 0) then
         llunrpt = 6
      end if

      if (present(lunerr)) then
         llunerr = lunerr
      end if
      if (llunerr .lt. 0) then
         llunerr = 6
      end if

      ! Ensure the problem size is valid
      if (n .le. 0) then
         linfo5 = 1
         linfo4 = 1
      end if

      if (m .le. 0) then
         linfo5 = 1
         linfo3 = 1
      end if

      if (np .le. 0) then
         linfo5 = 1
         linfo2 = 1
      end if

      if (nq .le. 0) then
         linfo5 = 1
         linfo1 = 1
      end if

      if (linfo5 .ne. 0) then
         linfo = 10000*linfo5 + 1000*linfo4 + 100*linfo3 + 10*linfo2 + linfo1
         if (llunerr .gt. 0 .and. liprint .ne. 0) then
            call dodphd(head, llunrpt)
            call dodpe1( &
               llunerr, linfo, linfo5, linfo4, linfo3, linfo2, linfo1, &
               n, m, nq, &
               ldscld, ldstpd, ldwe, ld2we, ldwd, ld2wd, &
               lenwork, leniwork &
               )
         end if
         if (present(info)) then
            info = linfo
         end if
         return
      end if

      ! Define LJOB and check that necessary arguments are passed for JOB
      if (present(job)) then
         ljob = job
         if (mod(job, 10000)/1000 .ge. 1) then
            if (.not. present(delta)) then
               linfo5 = 7
               linfo4 = 1
            elseif (.not. allocated(delta)) then
               linfo5 = 7
               linfo4 = 1
            end if
         end if
         if (job .ge. 10000) then
            if (.not. present(iwork)) then
               linfo5 = 7
               linfo2 = 1
            elseif (.not. associated(iwork)) then
               linfo5 = 7
               linfo2 = 1
            end if
         end if
         if (job .ge. 10000) then
            if (.not. present(work)) then
               linfo5 = 7
               linfo3 = 1
            elseif (.not. associated(work)) then
               linfo5 = 7
               linfo3 = 1
            end if
         end if
      else
         ljob = -1
      end if

      if (linfo5 .ne. 0) then
         linfo = 10000*linfo5 + 1000*linfo4 + 100*linfo3 + 10*linfo2 + linfo1
         if (llunerr .gt. 0 .and. liprint .ne. 0) then
            call dodphd(head, llunrpt)
            call dodpe1( &
               llunerr, linfo, linfo5, linfo4, linfo3, linfo2, linfo1, &
               n, m, nq, &
               ldscld, ldstpd, ldwe, ld2we, ldwd, ld2wd, &
               lenwork, leniwork &
               )
         end if
         if (present(info)) then
            info = linfo
         end if
         return
      end if

      ! Determine the size of WORK
      if (ljob .lt. 0 .or. mod(ljob, 10) .le. 1) then
         lenwork = 18 + 13*np + np**2 + m + m**2 + 4*n*nq + 6*n*m + 2*n*nq*np + &
                   2*n*nq*m + nq**2 + 5*nq + nq*(np + m) + n*nq*nq
      else
         lenwork = 18 + 13*np + np**2 + m + m**2 + 4*n*nq + 2*n*m + 2*n*nq*np + &
                   5*nq + nq*(np + m) + n*nq*nq
      end if

      ! Determine the size of IWORK
      leniwork = 20 + 2*np + nq*(np + m)

      ! Allocate the work arrays
      allocate (lwork(lenwork), tempret(max(n, np), max(nq, m)), stat=linfo3)
      allocate (liwork(leniwork), stat=linfo2)
      lwork = ZERO
      liwork = 0
      if (present(delta)) then
         if (.not. allocated(delta)) then
            allocate (delta(n, m), stat=linfo4)
         end if
      end if
      if (linfo4 .ne. 0 .or. linfo3 .ne. 0 .or. linfo2 .ne. 0) then
         linfo5 = 8
      end if

      if (linfo5 .ne. 0) then
         linfo = 10000*mod(linfo5, 10) + 1000*mod(linfo4, 10) + &
                 100*mod(linfo3, 10) + 10*mod(linfo2, 10) + mod(linfo1, 10)
         if (llunerr .gt. 0 .and. liprint .ne. 0) then
            call dodphd(head, llunrpt)
            call dodpe1( &
               llunerr, linfo, linfo5, linfo4, linfo3, linfo2, linfo1, &
               n, m, nq, &
               ldscld, ldstpd, ldwe, ld2we, ldwd, ld2wd, &
               lenwork, leniwork)
         end if
         if (present(info)) then
            info = linfo
         end if
         return
      end if

      ! Set array variable defaults except IWORK
      lwork(1:n*m) = ZERO
      lifixb(1) = -1
      lifixx(1, 1) = -1
      llower(1:np) = -huge(ZERO)
      lsclb(1) = negone
      lscld(1, 1) = negone
      lstpb(1) = negone
      lstpd(1, 1) = negone
      lupper(1:np) = huge(ZERO)
      lwe(1, 1, 1) = negone
      lwd(1, 1, 1) = negone

      ! Check the size of required arguments and return errors if they are too small
      if (size(beta) .lt. np) then
         linfo1 = linfo1 + 1
      end if

      if (any(size(y) .lt. (/n, nq/))) then
         linfo1 = linfo1 + 2
      end if

      if (any(size(x) .lt. (/n, m/))) then
         linfo1 = linfo1 + 4
      end if

      ! Check the presence of optional arguments and copy their values internally or
      ! report errors as necessary
      if (present(ifixb)) then
         if (size(ifixb) .lt. np) then
            linfo1 = linfo1 + 64
         end if
         if (ifixb(1) .lt. 0) then
            lifixb(1) = ifixb(1)
         else
            lifixb(1:np) = ifixb(1:np)
         end if
      end if

      if (present(ifixx)) then
         ldifx = size(ifixx, 1)
         if (any(size(ifixx) .le. (/0, 0/))) then
            linfo1 = linfo1 + 128
         end if
         if (.not. (ifixx(1, 1) .lt. 0 .or. ldifx .eq. 1 .or. ldifx .ge. n) &
             .or. size(ifixx, 2) .lt. m) then
            linfo1 = linfo1 + 128
         end if
         if (ldifx .gt. n) then
            ldifx = n
         end if
         if (ifixx(1, 1) .lt. 0) then
            lifixx(1, 1) = ifixx(1, 1)
         else
            lifixx(1:ldifx, 1:m) = ifixx(1:ldifx, 1:m)
         end if
      end if

      if (present(iwork)) then
         if (associated(iwork)) then
            if (size(iwork) .lt. leniwork) then
               linfo1 = linfo1 + 8192
            end if
            !  This is a restart, copy IWORK.
            if (mod(ljob/10000, 10) .ge. 1) then
               liwork(1:leniwork) = iwork(1:leniwork)
            end if
         end if
      end if

      if (present(maxit)) then
         lmaxit = maxit
      end if

      if (present(ndigit)) then
         lndigit = ndigit
      end if

      if (present(partol)) then
         lpartol = partol
      end if

      if (present(sclb)) then
         if (size(sclb) .lt. np) then
            linfo1 = linfo1 + 1024
         end if
         if (sclb(1) .le. ZERO) then
            lsclb(1) = sclb(1)
         else
            lsclb(1:np) = sclb(1:np)
         end if
      end if

      if (present(scld)) then
         ldscld = size(scld, 1)
         if (any(size(scld) .le. (/0, 0/))) then
            linfo1 = linfo1 + 2048
         end if
         if (.not. (scld(1, 1) .le. ZERO .or. ldscld .eq. 1 .or. ldscld .ge. n) &
             .or. size(scld, 2) .lt. m) then
            linfo1 = linfo1 + 2048
         end if
         if (ldscld .gt. n) then
            ldscld = n
         end if
         if (scld(1, 1) .le. ZERO) then
            lscld(1, 1) = scld(1, 1)
         else
            lscld(1:ldscld, 1:m) = scld(1:ldscld, 1:m)
         end if
      end if

      if (present(sstol)) then
         lsstol = sstol
      end if

      if (present(stpb)) then
         if (size(stpb) .lt. np) then
            linfo1 = linfo1 + 256
         end if
         if (stpb(1) .le. ZERO) then
            lstpb(1) = stpb(1)
         else
            lstpb(1:np) = stpb(1:np)
         end if
      end if

      if (present(stpd)) then
         ldstpd = size(stpd, 1)
         if (any(size(stpd) .le. (/0, 0/))) then
            linfo1 = linfo1 + 512
         end if
         if (.not. (stpd(1, 1) .le. ZERO .or. ldstpd .eq. 1 .or. ldstpd .ge. n) &
             .or. size(stpd, 2) .lt. m) then
            linfo1 = linfo1 + 512
         end if
         if (ldstpd .gt. n) then
            ldstpd = n
         end if
         if (stpd(1, 1) .le. ZERO) then
            lstpd(1, 1) = stpd(1, 1)
         else
            lstpd(1:ldstpd, 1:m) = stpd(1:ldstpd, 1:m)
         end if
      end if

      if (present(taufac)) then
         ltaufac = taufac
      end if

      if (present(we)) then
         ldwe = size(we, 1)
         ld2we = size(we, 2)
         if (any(size(we) .le. (/0, 0, 0/))) then
            linfo1 = linfo1 + 16
         end if
         if (.not. (we(1, 1, 1) .lt. ZERO .or. &
                    ((ldwe .eq. 1 .or. ldwe .ge. n) .and. &
                     (ld2we .eq. 1 .or. ld2we .ge. nq))) .or. size(we, 3) .lt. nq) then
            linfo1 = linfo1 + 16
         end if
         if (ldwe .gt. n) then
            ldwe = n
         end if
         if (ld2we .gt. nq) then
            ld2we = nq
         end if
         if (we(1, 1, 1) .lt. ZERO) then
            lwe(1, 1, 1) = we(1, 1, 1)
         else
            lwe(1:ldwe, 1:ld2we, 1:nq) = we(1:ldwe, 1:ld2we, 1:nq)
         end if
      end if

      if (present(wd)) then
         ldwd = size(wd, 1)
         ld2wd = size(wd, 2)
         if (any(size(wd) .le. (/0, 0, 0/))) then
            linfo1 = linfo1 + 32
         end if
         if (.not. (wd(1, 1, 1) .lt. ZERO .or. &
                    ((ldwd .eq. 1 .or. ldwd .ge. n) .and. &
                     (ld2wd .eq. 1 .or. ld2wd .ge. m))) .or. size(wd, 3) .lt. m) then
            linfo1 = linfo1 + 32
         end if
         if (ldwd .gt. n) then
            ldwd = n
         end if
         if (ld2wd .gt. m) then
            ld2wd = m
         end if
         if (wd(1, 1, 1) .le. 0.0_wp) then
            lwd(1, 1, 1) = wd(1, 1, 1)
         else
            lwd(1:ldwd, 1:ld2wd, 1:m) = wd(1:ldwd, 1:ld2wd, 1:m)
         end if
      end if

      if (present(work)) then
         if (associated(work)) then
            if (size(work) .lt. lenwork) then
               linfo1 = linfo1 + 4096
            end if
            !  Deltas are in WORK, copy them.
            if (mod(ljob/1000, 10) .ge. 1 .and. .not. present(delta)) then
               lwork(1:n*m) = work(1:n*m)
            end if
            !  This is a restart, copy WORK.
            if (mod(ljob/10000, 10) .ge. 1) then
               lwork(1:lenwork) = work(1:lenwork)
            end if
         end if
      end if

      if (present(delta)) then
         if (allocated(delta)) then
            if (any(shape(delta) .lt. (/n, m/))) then
               linfo1 = linfo1 + 8
            end if
            lwork(1:n*m) = reshape(delta(1:n, 1:m), (/n*m/))
         end if
      end if

      if (present(lower)) then
         if (size(lower) .lt. np) then
            linfo1 = linfo1 + 32768
         end if
         llower(1:np) = lower(1:np)
      end if

      if (present(upper)) then
         if (size(upper) .lt. np) then
            linfo1 = linfo1 + 16384
         end if
         lupper(1:np) = upper(1:np)
      end if

      ! Report an error if any of the array sizes didn't match.
      if (linfo1 .ne. 0) then
         linfo = 100000 + linfo1
         linfo1 = 0
         if (llunerr .gt. 0 .and. liprint .ne. 0) then
            call dodphd(head, llunrpt)
            call dodpe1( &
               llunerr, linfo, linfo5, linfo4, linfo3, linfo2, linfo1, &
               n, m, nq, &
               ldscld, ldstpd, ldwe, ld2we, ldwd, ld2wd, &
               lenwork, leniwork &
               )
         end if
         if (present(info)) then
            info = linfo
         end if
         return
      end if

      if (lwd(1, 1, 1) .ne. 0) then
         call dodcnt &
            (fcn, &
             n, m, np, nq, &
             beta(1:np), &
             y(1:n, 1:nq), n, x(1:n, 1:m), n, &
             lwe(1:ldwe, 1:ld2we, 1:nq), ldwe, ld2we, &
             lwd(1:ldwd, 1:ld2wd, 1:m), ldwd, ld2wd, &
             lifixb, lifixx(1:ldifx, 1:m), ldifx, &
             ljob, lndigit, ltaufac, &
             lsstol, lpartol, lmaxit, &
             liprint, llunerr, llunrpt, &
             lstpb, lstpd(1:ldstpd, 1:m), ldstpd, &
             lsclb, lscld(1:ldscld, 1:m), ldscld, &
             lwork, lenwork, liwork, leniwork, &
             linfo, &
             llower, lupper)
      else
         wd1(1, 1, 1) = negone
         call dodcnt &
            (fcn, &
             n, m, np, nq, &
             beta(1:np), &
             y(1:n, 1:nq), n, x(1:n, 1:m), n, &
             lwe(1:ldwe, 1:ld2we, 1:nq), ldwe, ld2we, &
             wd1, 1, 1, &
             lifixb, lifixx(1:ldifx, 1:m), ldifx, &
             ljob, lndigit, ltaufac, &
             lsstol, lpartol, lmaxit, &
             liprint, llunerr, llunrpt, &
             lstpb, lstpd(1:ldstpd, 1:m), ldstpd, &
             lsclb, lscld(1:ldscld, 1:m), ldscld, &
             lwork, lenwork, liwork, leniwork, &
             linfo, &
             llower, lupper)
      end if

      if (present(delta)) then
         if (allocated(delta)) then
            delta(1:n, 1:m) = reshape(lwork(1:n*m), (/n, m/))
         end if
      end if

      if (present(info)) then
         info = linfo
      end if

      if (present(iwork)) then
         if (.not. associated(iwork)) then
            iwork => liwork
         else
            iwork(1:leniwork) = liwork(1:leniwork)
            deallocate (liwork)
         end if
      else
         deallocate (liwork)
      end if

      if (present(work)) then
         if (.not. associated(work)) then
            work => lwork
         else
            work(1:lenwork) = lwork(1:lenwork)
            deallocate (lwork)
         end if
      else
         deallocate (lwork)
      end if

      deallocate (tempret)

   end subroutine odr

end module odrpack

subroutine dacces                                                             &
         ( n, m, np, nq, ldwe, ld2we,                                         &
         work, lwork, iwork, liwork,                                          &
         access, isodr,                                                       &
         jpvt, omega, u, qraux, sd, vcv, wrk1, wrk2, wrk3, wrk4, wrk5, wrk6,  &
         nnzw, npp,                                                           &
         job, partol, sstol, maxit, taufac, eta, neta,                        &
         lunrpt, ipr1, ipr2, ipr2f, ipr3,                                     &
         wss, rvar, idf,                                                      &
         tau, alpha, niter, nfev, njev, int2, olmavg,                         &
         rcond, irank, actrs, pnorm, prers, rnorms, istop)
!***Begin Prologue  DACCES
!***Refer to  ODR
!***Routines Called  DIWINF,DWINF
!***Date Written   860529   (YYMMDD)
!***Revision Date  920619   (YYMMDD)
!***Purpose  Access or store values in the work arrays
!***End Prologue  DACESS
!
!...Used modules
        use odrpack_kinds,only: wp
!
!...Scalar arguments
        real(kind = wp)                                                       &
         actrs, alpha, eta, olmavg, partol, pnorm, prers, rcond,              &
         rnorms, rvar, sstol, tau, taufac
        integer                                                               &
         idf, int2, ipr1, ipr2, ipr2f, ipr3, irank, istop, istopi, job, jpvt, &
!-----------------^------------------------------------------------------------
!!! FPT - 1271 Non-standard Fortran intrinsic(s) used as local identifier(s)
!------------------------------------------------------------------------------
         ldwe, ld2we, liwork, lunrpt, lwork, m, maxit, n, neta, nfev, niter,  &
         njev, nnzw, np, npp, nq, omega, qraux, sd, u, vcv, wrk1, wrk2, wrk3, &
         wrk4, wrk5, wrk6
        logical                                                               &
         access, isodr
!--------------^---------------------------------------------------------------
!!! FPT - 1273 Fortran auxiliary keyword used as identifier name.
!------------------------------------------------------------------------------
!
!...Array arguments
        real(kind = wp)                                                       &
         work( lwork), wss(3)
        integer                                                               &
         iwork( liwork)
!
!...Local scalars
        integer                                                               &
         actrsi, alphai, betaci, betani, betasi, beta0i, boundi,              &
         deltai, deltni, deltsi, diffi, epsi,                                 &
         epsmai, etai, fjacbi, fjacdi, fni, fsi, idfi, int2i, iprini, iprint, &
         iranki, jobi, jpvti, ldtti, liwkmn, loweri, luneri, lunrpi, lwkmn,   &
         maxiti,                                                              &
         msgb, msgd, netai, nfevi, niteri, njevi, nnzwi, nppi, nrowi,         &
         ntoli, olmavi, omegai, partli, pnormi, prersi, qrauxi, rcondi,       &
         rnorsi, rvari, sdi, si, ssfi, ssi, sstoli, taufci, taui, ti, tti, ui &
         , upperi, vcvi, we1i, wrk1i, wrk2i, wrk3i, wrk4i, wrk5i, wrk6i,      &
         wrk7i, wssi, wssdei, wssepi, xplusi
!...External subroutines
        external                                                              &
         diwinf, dwinf
!
!...Variable Definitions (alphabetically)
!       ACCESS:  The variable designating whether information is to be
!       accessed from the work arrays (ACCESS=TRUE) or stored in
!       them (ACCESS=FALSE).
!       ACTRS:   The saved actual relative reduction in the sum-of-squares.
!       ACTRSI:  The location in array WORK of variable ACTRS.
!       ALPHA:   The Levenberg-Marquardt parameter.
!       ALPHAI:  The location in array WORK of variable ALPHA.
!       BETACI:  The starting location in array WORK of array BETAC.
!       BETANI:  The starting location in array WORK of array BETAN.
!       BETASI:  The starting location in array WORK of array BETAS.
!       BETA0I:  The starting location in array WORK of array BETA0.
!       DELTAI:  The starting location in array WORK of array DELTA.
!       DELTNI:  The starting location in array WORK of array DELTAN.
!       DELTSI:  The starting location in array WORK of array DELTAS.
!       DIFFI:   The starting location in array WORK of array DIFF.
!       EPSI:    The starting location in array WORK of array EPS.
!       EPSMAI:  The location in array WORK of variable EPSMAC.
!       ETA:     The relative noise in the function results.
!       ETAI:    The location in array WORK of variable ETA.
!       FJACBI:  The starting location in array WORK of array FJACB.
!       FJACDI:  The starting location in array WORK of array FJACD.
!       FNI:     The starting location in array WORK of array FN.
!       FSI:     The starting location in array WORK of array FS.
!       IDF:     The degrees of freedom of the fit, equal to the number of
!       observations with nonzero weighted derivatives minus the
!       number of parameters being estimated.
!       IDFI:    The starting location in array IWORK of variable IDF.
!       INT2:    The number of internal doubling steps.
!       INT2I:   The location in array IWORK of variable INT2.
!       IPR1:    The value of the fourth digit (from the right) of IPRINT,
!       which controls the initial summary report.
!       IPR2:    The value of the third digit (from the right) of IPRINT,
!       which controls the iteration reports.
!       IPR2F:   The value of the second digit (from the right) of IPRINT,
!       which controls the frequency of the iteration reports.
!       IPR3:    The value of the first digit (from the right) of IPRINT,
!       which controls the final summary report.
!       IPRINI:  The location in array IWORK of variable IPRINT.
!       IPRINT:  The print control variable.
!       IRANK:   The rank deficiency of the Jacobian wrt BETA.
!       IRANKI:  The location in array IWORK of variable IRANK.
!       ISODR:   The variable designating whether the solution is to be
!       found by ODR (ISODR=TRUE) or by OLS (ISODR=FALSE).
!       ISTOP:   The variable designating whether there are problems
!       computing the function at the current BETA and DELTA.
!       ISTOPI:  The location in array IWORK of variable ISTOP.
!       IWORK:   The integer work space.
!       JOB:     The variable controling problem initialization and
!       computational method.
!       JOBI:    The location in array IWORK of variable JOB.
!       JPVT:    The pivot vector.
!       JPVTI:   The starting location in array IWORK of variable JPVT.
!       LDTTI:   The starting location in array IWORK of variable LDTT.
!       LDWE:    The leading dimension of array WE.
!       LD2WE:   The second dimension of array WE.
!       LIWORK:  The length of vector IWORK.
!       LUNERI:  The location in array IWORK of variable LUNERR.
!       LUNERR:  The logical unit number used for error messages.
!       LUNRPI:  The location in array IWORK of variable LUNRPT.
!       LUNRPT:  The logical unit number used for computation reports.
!       LWKMN:   The minimum acceptable length of array WORK.
!       LWORK:   The length of vector WORK.
!       M:       The number of columns of data in the explanatory variable.
!       MAXIT:   The maximum number of iterations allowed.
!       MAXITI:  The location in array IWORK of variable MAXIT.
!       MSGB:    The starting location in array IWORK of array MSGB.
!       MSGD:    The starting location in array IWORK of array MSGD.
!       N:       The number of observations.
!       NETA:    The number of accurate digits in the function results.
!       NETAI:   The location in array IWORK of variable NETA.
!       NFEV:    The number of function evaluations.
!       NFEVI:   The location in array IWORK of variable NFEV.
!       NITER:   The number of iterations taken.
!       NITERI:  The location in array IWORK of variable NITER.
!       NJEV:    The number of Jacobian evaluations.
!       NJEVI:   The location in array IWORK of variable NJEV.
!       NNZW:    The number of nonzero weighted observations.
!       NNZWI:   The location in array IWORK of variable NNZW.
!       NP:      The number of function parameters.
!       NPP:     The number of function parameters actually estimated.
!       NPPI:    The location in array IWORK of variable NPP.
!       NQ:      The number of responses per observation.
!       NROWI:   The location in array IWORK of variable NROW.
!       NTOLI:   The location in array IWORK of variable NTOL.
!       OLMAVG:  The average number of Levenberg-Marquardt steps per
!       iteration.
!       OLMAVI:  The location in array WORK of variable OLMAVG.
!       OMEGA:   The starting location in array WORK of array OMEGA.
!       OMEGAI:  The starting location in array WORK of array OMEGA.
!       PARTLI:  The location in array work of variable PARTOL.
!       PARTOL:  The parameter convergence stopping tolerance.
!       PNORM:   The norm of the scaled estimated parameters.
!       PNORMI:  The location in array WORK of variable PNORM.
!       PRERS:   The saved predicted relative reduction in the
!       sum-of-squares.
!       PRERSI:  The location in array WORK of variable PRERS.
!       QRAUX:   The starting location in array WORK of array QRAUX.
!       QRAUXI:  The starting location in array WORK of array QRAUX.
!       RCOND:   The approximate reciprocal condition of FJACB.
!       RCONDI:  The location in array WORK of variable RCOND.
!       RESTRT:  The variable designating whether the call is a restart
!       (RESTRT=TRUE) or not (RESTRT=FALSE).
!       RNORMS:  The norm of the saved weighted EPSILONS and DELTAS.
!       RNORSI:  The location in array WORK of variable RNORMS.
!       RVAR:    The residual variance, i.e. standard deviation squared.
!       RVARI:   The location in array WORK of variable RVAR.
!       SCLB:    The scaling values used for BETA.
!       SCLD:    The scaling values used for DELTA.
!       SD:      The starting location in array WORK of array SD.
!       SDI:     The starting location in array WORK of array SD.
!       SI:      The starting location in array WORK of array S.
!       SSFI:    The starting location in array WORK of array SSF.
!       SSI:     The starting location in array WORK of array SS.
!       SSTOL:   The sum-of-squares convergence stopping tolerance.
!       SSTOLI:  The location in array WORK of variable SSTOL.
!       TAU:     The trust region diameter.
!       TAUFAC:  The factor used to compute the initial trust region
!       diameter.
!       TAUFCI:  The location in array WORK of variable TAUFAC.
!       TAUI:    the location in array WORK of variable TAU.
!       TI:      The starting location in array WORK of array T.
!       TTI:     The starting location in array WORK of array TT.
!       U:       The starting location in array WORK of array U.
!       UI:      The starting location in array WORK of array U.
!       VCV:     The starting location in array WORK of array VCV.
!       VCVI:    The starting location in array WORK of array VCV.
!       WE1I:    The starting location in array WORK of array WE1.
!       WORK:    The REAL (KIND=wp) work space.
!       WRK1:    The starting location in array WORK of array WRK1.
!       WRK1I:   The starting location in array WORK of array WRK1.
!       WRK2:    The starting location in array WORK of array WRK2.
!       WRK2I:   The starting location in array WORK of array WRK2.
!       WRK3:    The starting location in array WORK of array wrk3.
!       WRK3I:   The starting location in array WORK of array wrk3.
!       WRK4:    The starting location in array WORK of array wrk4.
!       WRK4I:   The starting location in array WORK of array wrk4.
!       WRK5:    The starting location in array WORK of array wrk5.
!       WRK5I:   The starting location in array WORK of array wrk5.
!       WRK6:    The starting location in array WORK of array wrk6.
!       WRK6I:   The starting location in array WORK of array wrk6.
!       WRK7I:   The starting location in array WORK of array wrk7.
!       WSS:     The sum of the squares of the weighted EPSILONS and DELTAS,
!       the sum of the squares of the weighted DELTAS, and
!       the sum of the squares of the weighted EPSILONS.
!       WSSI:    The starting location in array WORK of variable WSS(1).
!       WSSDEI:  The starting location in array WORK of variable WSS(2).
!       WSSEPI:  The starting location in array WORK of variable WSS(3).
!       XPLUSI:  The starting location in array WORK of array XPLUSD.
!
!
!***First executable statement  DACCES
!
!
!  Find starting locations within integer workspace
!
        call diwinf( m, np, nq,                                               &
         msgb, msgd, jpvti, istopi,                                           &
         nnzwi, nppi, idfi,                                                   &
         jobi, iprini, luneri, lunrpi,                                        &
         nrowi, ntoli, netai,                                                 &
         maxiti, niteri, nfevi, njevi, int2i, iranki, ldtti,                  &
         boundi,                                                              &
         liwkmn)
!
!  Find starting locations within REAL (KIND=wp) work space
!
        call dwinf( n, m, np, nq, ldwe, ld2we, isodr,                         &
         deltai, epsi, xplusi, fni, sdi, vcvi,                                &
         rvari, wssi, wssdei, wssepi, rcondi, etai,                           &
         olmavi, taui, alphai, actrsi, pnormi, rnorsi, prersi,                &
         partli, sstoli, taufci, epsmai,                                      &
         beta0i, betaci, betasi, betani, si, ssi, ssfi, qrauxi, ui,           &
         fsi, fjacbi, we1i, diffi,                                            &
         deltsi, deltni, ti, tti, omegai, fjacdi,                             &
         wrk1i, wrk2i, wrk3i, wrk4i, wrk5i, wrk6i, wrk7i,                     &
         loweri, upperi,                                                      &
         lwkmn)
!
        if ( access) then
!
!  Set starting locations for work vectors
!
           jpvt = jpvti
           omega = omegai
           qraux = qrauxi
           sd = sdi
           vcv = vcvi
           u = ui
           wrk1 = wrk1i
           wrk2 = wrk2i
           wrk3 = wrk3i
           wrk4 = wrk4i
           wrk5 = wrk5i
           wrk6 = wrk6i
!
!  Access values from the work vectors
!
           actrs = work( actrsi)
           alpha = work( alphai)
           eta = work( etai)
           olmavg = work( olmavi)
           partol = work( partli)
           pnorm = work( pnormi)
           prers = work( prersi)
           rcond = work( rcondi)
           wss(1) = work( wssi)
           wss(2) = work( wssdei)
           wss(3) = work( wssepi)
           rvar = work( rvari)
           rnorms = work( rnorsi)
           sstol = work( sstoli)
           tau = work( taui)
           taufac = work( taufci)
!
           neta = iwork( netai)
           irank = iwork( iranki)
           job = iwork( jobi)
           lunrpt = iwork( lunrpi)
           maxit = iwork( maxiti)
           nfev = iwork( nfevi)
           niter = iwork( niteri)
           njev = iwork( njevi)
           nnzw = iwork( nnzwi)
           npp = iwork( nppi)
           idf = iwork( idfi)
           int2 = iwork( int2i)
!
!  Set up print control variables
!
           iprint = iwork( iprini)
!
           ipr1 = mod( iprint,10000)/1000
           ipr2 = mod( iprint,1000)/100
           ipr2f = mod( iprint,100)/10
           ipr3 = mod( iprint,10)
!
        else
!
!  Store values into the work vectors
!
           work( actrsi) = actrs
           work( alphai) = alpha
           work( olmavi) = olmavg
           work( partli) = partol
           work( pnormi) = pnorm
           work( prersi) = prers
           work( rcondi) = rcond
           work( wssi) = wss(1)
           work( wssdei) = wss(2)
           work( wssepi) = wss(3)
           work( rvari) = rvar
           work( rnorsi) = rnorms
           work( sstoli) = sstol
           work( taui) = tau
!
           iwork( iranki) = irank
           iwork( istopi) = istop
           iwork( nfevi) = nfev
           iwork( niteri) = niter
           iwork( njevi) = njev
           iwork( idfi) = idf
           iwork( int2i) = int2
        endif
!
        return
end subroutine

subroutine desubi( n, m, wd, ldwd, ld2wd, alpha, tt, ldtt, i, e)
!***Begin Prologue  DESUBI
!***Refer to  ODR
!***Routines Called 
!***Date Written   860529   (YYMMDD)
!***Revision Date  920304   (YYMMDD)
!***Purpose  Compute E = WD + ALPHA*TT**2
!***End Prologue  DESUBI

!...Used modules
        use odrpack_kinds, only: wp, zero

!...Scalar arguments
        real(kind = wp) :: alpha
        integer :: ldtt, ldwd, ld2wd, m, n

!...Array arguments
        real(kind = wp) :: e( m, m), tt( ldtt, m), wd( ldwd, ld2wd, m)

!...Local scalars
        integer :: i, j, j1, j2

!...Variable Definitions (alphabetically)
!       ALPHA:  The Levenberg-Marquardt parameter.
!       E:      The value of the array E = WD + ALPHA*TT**2
!       I:      An indexing variable.
!       J:      An indexing variable.
!       J1:     An indexing variable.
!       J2:     An indexing variable.
!       LDWD:   The leading dimension of array WD.
!       LD2WD:  The second dimension of array WD.
!       M:      The number of columns of data in the independent variable.
!       N:      The number of observations.
!       NP:     The number of responses per observation.
!       TT:     The scaling values used for DELTA.
!       WD:     The squared DELTA weights, D**2.
!       ZERO:   The value 0.0E0_wp.

!***First executable statement  DESUBI

!       N.B. the locations of WD and TT accessed depend on the value
!       of the first element of each array and the leading dimensions
!       of the multiply subscripted arrays.

        if ( n .eq. 0 .or. m .eq. 0) return

        if ( wd(1,1,1) .ge. zero) then
           if ( ldwd .ge. n) then
!  The elements of WD have been individually specified
!
              if ( ld2wd .eq. 1) then
!  The arrays stored in WD are diagonal
                 e = zero 
                 do 10 j = 1, m
                    e( j, j) = wd( i,1, j)
10               continue
              else
!  The arrays stored in WD are full positive semidefinite matrices
                 do 30 j1 = 1, m
                    do 20 j2 = 1, m
                       e( j1, j2) = wd( i, j1, j2)
20                  continue
30               continue
              endif
!
              if ( tt(1,1) .gt. zero) then
                 if ( ldtt .ge. n) then
                    do 110 j = 1, m
                       e( j, j) = e( j, j)+ alpha* tt( i, j)**2
110                 continue
                 else
                    do 120 j = 1, m
                       e( j, j) = e( j, j)+ alpha* tt(1, j)**2
120                 continue
                 endif
              else
                 do 130 j = 1, m
                    e( j, j) = e( j, j)+ alpha* tt(1,1)**2
130              continue
              endif
           else
!  WD is an M by M matrix
!
              if ( ld2wd .eq. 1) then
!  The array stored in WD is diagonal
                 e = zero
                 do 140 j = 1, m
                    e( j, j) = wd(1,1, j)
140              continue
              else
!  The array stored in WD is a full positive semidefinite matrices
                 do 160 j1 = 1, m
                    do 150 j2 = 1, m
                       e( j1, j2) = wd(1, j1, j2)
150                 continue
160              continue
              endif
!
              if ( tt(1,1) .gt. zero) then
                 if ( ldtt .ge. n) then
                    do 210 j = 1, m
                       e( j, j) = e( j, j)+ alpha* tt( i, j)**2
210                 continue
                 else
                    do 220 j = 1, m
                       e( j, j) = e( j, j)+ alpha* tt(1, j)**2
220                 continue
                 endif
              else
                 do 230 j = 1, m
                    e( j, j) = e( j, j)+ alpha* tt(1,1)**2
230              continue
              endif
           endif
        else
!  WD is a diagonal matrix with elements ABS(WD(1,1,1))
           e = zero
           if ( tt(1,1) .gt. zero) then
              if ( ldtt .ge. n) then
                 do 310 j = 1, m
                    e( j, j) = abs( wd(1,1,1))+ alpha* tt( i, j)**2
310              continue
              else
                 do 320 j = 1, m
                    e( j, j) = abs( wd(1,1,1))+ alpha* tt(1, j)**2
320              continue
              endif
           else
              do 330 j = 1, m
                 e( j, j) = abs( wd(1,1,1))+ alpha* tt(1,1)**2
330           continue
           endif
        endif

end subroutine

subroutine detaf                                                              &
         ( fcn,                                                               &
         n, m, np, nq,                                                        &
         xplusd, beta, epsmac, nrow,                                          &
         partmp, pv0,                                                         &
         ifixb, ifixx, ldifx,                                                 &
         istop, nfev, eta, neta,                                              &
         wrk1, wrk2, wrk6, wrk7,                                              &
         info,                                                                &
         lower, upper)
!***Begin Prologue  DETAF
!***Refer to  ODR
!***Routines Called  FCN
!***Date Written   860529   (YYMMDD)
!***Revision Date  920619   (YYMMDD)
!***Purpose  Compute noise and number of good digits in function results
!       (Adapted from STARPAC subroutine ETAFUN)
!***End Prologue  DETAF
!
!...Used modules
        use odrpack_kinds,only: wp
!
!...Scalar arguments
        real(kind = wp)                                                       &
         epsmac, eta
        integer                                                               &
         info, istop, ldifx, m, n, neta, nfev, np, nq, nrow
!
!...Array arguments
        real(kind = wp)                                                       &
         beta( np), lower( np), partmp( np), pv0( n, nq), upper( np),         &
         wrk1( n, m, nq), wrk2( n, nq), wrk6( n, np, nq), wrk7(-2:2, nq),     &
         xplusd( n, m)
        integer                                                               &
         ifixb( np), ifixx( ldifx, m)
!
!...Subroutine arguments
        external                                                              &
         fcn
!
!...Local scalars
        real(kind = wp)                                                       &
         a, b, fac, hundrd, one, p1, p2, p5, shift, stp, two, zero
        integer                                                               &
         j, k, l, sbk
!
!...Local arrays
        real(kind = wp)                                                       &
         parpts(-2:2, np)
!
!...Data statements
        data                                                                  &
         zero, p1, p2, p5, one, two, hundrd                                   &
         /0.0E0_wp,0.1E0_wp,0.2E0_wp,0.5E0_wp,1.0E0_wp,2.0E0_wp,              &
         1.0E2_wp/
!
!...Routine names used as subprogram arguments
!       FCN:      The user supplied subroutine for evaluating the model.
!
!...Variable Definitions (ALPHABETICALLY)
!       A:       Parameters of the local fit.
!       B:       Parameters of the local fit.
!       BETA:    The function parameters.
!       EPSMAC:  The value of machine precision.
!       ETA:     The noise in the model results.
!       FAC:     A factor used in the computations.
!       HUNDRD:  The value 1.0E2_wp.
!       IFIXB:   The values designating whether the elements of BETA are
!       fixed at their input values or not.
!       IFIXX:   The values designating whether the elements of X are
!       fixed at their input values or not.
!       ISTOP:   The variable designating whether there are problems
!       Computing the function at the current BETA and DELTA.
!       J:       An index variable.
!       K:       An index variable.
!       L:       AN INDEX VARIABLE.
!       LDIFX:   The leading dimension of array IFIXX.
!       LOWER:   The lower bound of BETA.
!       M:       The number of columns of data in the explanatory variable.
!       N:       The number of observations.
!       NETA:    The number of accurate digits in the model results.
!       NFEV:    The number of function evaluations.
!       NP:      The number of function parameters.
!       NQ:      The number of responses per observation.
!       NROW:    The row number at which the derivative is to be checked.
!       ONE:     The value 1.0E0_wp.
!       P1:      The value 0.1E0_wp.
!       P2:      The value 0.2E0_wp.
!       P5:      The value 0.5E0_wp.
!       PARPTS:  The points that PARTMP will take on during FCN evaluations.
!       PARTMP:  The model parameters.
!       PV0:     The original predicted values.
!       SHIFT:   When PARPTS cross the parameter bounds they are shifted by SHIFT.
!       SBK:     The sign of BETA(K).
!       STP:     A small value used to perturb the parameters.
!       UPPER:   The upper bound of BETA.
!       WRK1:    A work array of (N BY M BY NQ) elements.
!       WRK2:    A work array of (N BY NQ) elements.
!       WRK6:    A work array of (N BY NP BY NQ) elements.
!       WRK7:    A work array of (5 BY NQ) elements.
!       XPLUSD:  The values of X + DELTA.
!       ZERO:    The value 0.0E0_wp.
!
!
!***First executable statement  DETAF
!
!
        stp = hundrd* epsmac
        eta = epsmac
!
!       Create points to use in calculating FCN for ETA and NETA.
        do j = -2,2
           if ( j .eq. 0) then
              parpts(0,:) = beta(:)
           else
              do k = 1, np
                 if ( ifixb(1) .lt. 0) then
                    parpts( j, k) = beta( k)+ j* stp* beta( k)
                 elseif ( ifixb( k) .ne. 0) then
                    parpts( j, k) = beta( k)+ j* stp* beta( k)
                 else
                    parpts( j, k) = beta( k)
                 endif
              enddo
           endif
        enddo
!
!       Adjust the points used in calculating FCN to uphold the boundary
!       constraints.
        do k = 1, np
           sbk = sign( one, parpts(2, k)- parpts(-2, k))
           if ( parpts( sbk*2, k) .gt. upper( k)) then
              shift = upper( k)- parpts( sbk*2, k)
              parpts( sbk*2, k) = upper( k)
              do j = - sbk*2, sbk*1, sbk
                 parpts( j, k) = parpts( j, k)+ shift
              enddo
              if ( parpts(- sbk*2, k) .lt. lower( k)) then
                 info = 90010
                 return
              endif
           endif
           if ( parpts(- sbk*2, k) .lt. lower( k)) then
              shift = lower( k)- parpts(- sbk*2, k)
              parpts(- sbk*2, k) = lower( k)
              do j = - sbk*1, sbk*2, sbk
                 parpts( j, k) = parpts( j, k)+ shift
              enddo
              if ( parpts( sbk*2, k) .gt. upper( k)) then
                 info = 90010
                 return
              endif
           endif
        enddo
!
!       Evaluate FCN for all points in PARPTS.
        do j = -2,2
           if (all( parpts( j,:) .eq. beta(:))) then
              do l = 1, nq
                 wrk7( j, l) = pv0( nrow, l)
              enddo
           else
              partmp(:) = parpts( j,:)
              istop = 0
              call fcn( n, m, np, nq,                                         &
               n, m, np,                                                      &
               partmp(:), xplusd,                                             &
               ifixb, ifixx, ldifx,                                           &
               003, wrk2, wrk6, wrk1, istop)
              if ( istop .ne. 0) then
                 return
              else
                 nfev = nfev+1
              endif
              do l = 1, nq
                 wrk7( j, l) = wrk2( nrow, l)
              enddo
           endif
        enddo
!
!       Calculate ETA and NETA.
        do 100 l = 1, nq
           a = zero
           b = zero
           do 50 j = -2,2
              a = a+ wrk7( j, l)
              b = b+ j* wrk7( j, l)
50         continue
           a = p2* a
           b = p1* b
           if (( wrk7(0, l) .ne. zero) .and.                                  &
!--------------------------------^---------------------------------------------
!!! FPT - 3087 REAL or COMPLEX quantity tested for exact equality/inequality
!------------------------------------------------------------------------------
               (abs( wrk7(1, l)+ wrk7(-1, l)) .gt. hundrd* epsmac)) then
              fac = one/abs( wrk7(0, l))
           else
              fac = one
           endif
           do 60 j = -2,2
              wrk7( j, l) = abs(( wrk7( j, l)-( a+ j* b))* fac)
              eta = max( wrk7( j, l), eta)
60         continue
100     continue
        neta = max( two, p5-log10( eta))
!
        return
end subroutine
!DEVJAC
subroutine devjac                                                             &
         ( fcn,                                                               &
         anajac, cdjac,                                                       &
         n, m, np, nq,                                                        &
         betac, beta, stpb,                                                   &
         ifixb, ifixx, ldifx,                                                 &
         x, ldx, delta, xplusd, stpd, ldstpd,                                 &
         ssf, tt, ldtt, neta, fn,                                             &
         stp, wrk1, wrk2, wrk3, wrk6,                                         &
         fjacb, isodr, fjacd, we1, ldwe, ld2we,                               &
         njev, nfev, istop, info,                                             &
         lower, upper)
!***Begin Prologue  DEVJAC
!***Refer to  ODR
!***Routines Called  FCN,DDOT,DIFIX,DJACCD,DJACFD,DWGHT,DUNPAC
!***Date Written   860529   (YYMMDD)
!***Revision Date  920304   (YYMMDD)
!***Purpose  Compute the weighted Jacobians wrt BETA and DELTA
!***End Prologue  DEVJAC
!
!...Used modules
        use odrpack_kinds, only: wp, zero
        use odrpack,only: tempret
!
!...Scalar arguments
        integer                                                               &
         info, istop, ldifx, ldstpd, ldtt, ldwe, ldx, ld2we,                  &
         m, n, neta, nfev, njev, np, nq
        logical                                                               &
         anajac, cdjac, isodr
!
!...Array arguments
        real(kind = wp)                                                       &
         beta( np), betac( np), delta( n, m), fjacb( n, np, nq), fjacd( n, m, &
         nq), fn( n, nq), lower( np), ssf( np), stp( n), stpb( np), stpd(     &
         ldstpd, m), tt( ldtt, m), upper( np), we1( ldwe, ld2we, nq), wrk1( n &
         , m, nq), wrk2( n, nq), wrk3( np), wrk6( n, np, nq), x( ldx, m),     &
         xplusd( n, m)
        integer                                                               &
         ifixb( np), ifixx( ldifx, m)

!...Subroutine arguments
        external                                                              &
         fcn

!...Local scalars
        integer :: ideval, j, k, k1, l
        logical :: error
!-------------^----------------------------------------------------------------
!!! FPT - 1269 Fortran keyword used as local identifier name.
!------------------------------------------------------------------------------

!...External subroutines
        external :: difix, djaccd, djacfd, dunpac

!...External functions
        real(kind = wp), external :: ddot

!...Interface blocks
        interface
           subroutine dwght                                                   &
            ( n, m, wt, ldwt, ld2wt, t, wtt)
           use odrpack_kinds,only: wp
           integer                                                            &
            ldwt, ld2wt, m, n
           real(kind = wp)                                                    &
            t(:,:), wt(:,:,:), wtt(:,:)
           end subroutine
        end interface
!
!...Routine names used as subprogram arguments
!       FCN:     The user-supplied subroutine for evaluating the model.
!
!...Variable Definitions (alphabetically)
!       ANAJAC:  The variable designating whether the Jacobians are
!       computed by finite differences (ANAJAC=FALSE) or not
!       (ANAJAC=TRUE).
!       BETA:    The function parameters.
!       BETAC:   The current estimated values of the unfixed BETA's.
!       CDJAC:   The variable designating whether the Jacobians are
!       computed by central differences (CDJAC=TRUE) or by forward
!       differences (CDJAC=FALSE).
!       DELTA:   The estimated values of DELTA.
!       ERROR:   The variable designating whether ODRPACK95 detected nonzero
!       values in array DELTA in the OLS case, and thus whether
!       the user may have overwritten important information
!       by computing FJACD in the OLS case.
!       FJACB:   The Jacobian with respect to BETA.
!       FJACD:   The Jacobian with respect to DELTA.
!       FN:      The predicted values of the function at the current point.
!       IDEVAL:  The variable designating what computations are to be
!       performed by user-supplied subroutine FCN.
!       IFIXB:   The values designating whether the elements of BETA are
!       fixed at their input values or not.
!       IFIXX:   The values designating whether the elements of DELTA are
!       fixed at their input values or not.
!       INFO:    The variable designating why the computations were stopped.
!       ISTOP:   The variable designating that the user wishes the
!       computations stopped.
!       ISODR:   The variable designating whether the solution is by ODR
!       (ISODR=TRUE) or OLS (ISODR=FALSE).
!       J:       An indexing variable.
!       K:       An indexing variable.
!       K1:      An indexing variable.
!       L:       An indexing variable.
!       LDIFX:   The leading dimension of array IFIXX.
!       LDSTPD:  The leading dimension of array STPD.
!       LDTT:    The leading dimension of array TT.
!       LDWE:    The leading dimension of arrays WE and WE1.
!       LDX:     The leading dimension of array X.
!       LD2WE:   The second dimension of arrays WE and WE1.
!       M:       The number of columns of data in the independent variable.
!       N:       The number of observations.
!       NETA:    The number of accurate digits in the function results.
!       NFEV:    The number of function evaluations.
!       NJEV:    The number of Jacobian evaluations.
!       NP:      The number of function parameters.
!       NQ:      The number of responses per observation.
!       SSF:     The scale used for the BETA's.
!       STP:     The step used for computing finite difference
!       derivatives with respect to DELTA.
!       STPB:    The relative step used for computing finite difference
!       derivatives with respect to BETA.
!       STPD:    The relative step used for computing finite difference
!       derivatives with respect to DELTA.
!       TT:      The scaling values used for DELTA.
!       WE1:     The square roots of the EPSILON weights in array WE.
!       WRK1:    A work array of (N by M by NQ) elements.
!       WRK2:    A work array of (N by NQ) elements.
!       WRK3:    A work array of (NP) elements.
!       WRK6:    A work array of (N BY NP BY NQ) elements.
!       X:       The independent variable.
!       XPLUSD:  The values of X + DELTA.
!       ZERO:    The value 0.0E0_wp.
!
!
!***First executable statement  DEVJAC
!
!
!  Insert current unfixed BETA estimates into BETA
      call dunpac( np, betac, beta, ifixb)

!  Compute XPLUSD = X + DELTA
      xplusd = x(1:n,:) + delta

!  Compute the Jacobian wrt the estimated BETAS (FJACB) and
!       the Jacobian wrt DELTA (FJACD)
        istop = 0
        if ( isodr) then
           ideval = 110
        else
           ideval = 010
        endif
        if ( anajac) then
           call fcn( n, m, np, nq,                                            &
            n, m, np,                                                         &
            beta, xplusd,                                                     &
            ifixb, ifixx, ldifx,                                              &
            ideval, wrk2, fjacb, fjacd,                                       &
            istop)
           if ( istop .ne. 0) then
              return
           else
              njev = njev+1
           endif
!  Make sure fixed elements of FJACD are zero
           if ( isodr) then
              do 10 l = 1, nq
                 call difix( n, m, ifixx, ldifx, fjacd(1,1, l), n, fjacd(1,1, &
                  l), n)
10            continue
           endif
        elseif ( cdjac) then
           call djaccd( fcn,                                                  &
            n, m, np, nq,                                                     &
            beta, x, ldx, delta, xplusd, ifixb, ifixx, ldifx,                 &
            stpb, stpd, ldstpd,                                               &
            ssf, tt, ldtt, neta, fn, stp, wrk1, wrk2, wrk3, wrk6,             &
            fjacb, isodr, fjacd, nfev, istop, info,                           &
            lower, upper)
        else
           call djacfd( fcn,                                                  &
            n, m, np, nq,                                                     &
            beta, x, ldx, delta, xplusd, ifixb, ifixx, ldifx,                 &
            stpb, stpd, ldstpd,                                               &
            ssf, tt, ldtt, neta, fn, stp, wrk1, wrk2, wrk3, wrk6,             &
            fjacb, isodr, fjacd, nfev, istop, info,                           &
            lower, upper)
        endif
        if ( istop .lt. 0 .or. info .ge. 10000) then
           return
        elseif (.not. isodr) then
!  Try to detect whether the user has computed JFACD
!  Within FCN in the OLS case
           error = ddot( n* m, delta,1, delta,1) .ne. zero
!-----------------------------------------------------^------------------------
!!! FPT - 3087 REAL or COMPLEX quantity tested for exact equality/inequality
!------------------------------------------------------------------------------
           if ( error) then
              info = 50300
              return
           endif
        endif
!
!  Weight the Jacobian wrt the estimated BETAS
!
        if ( ifixb(1) .lt. 0) then
           do 20 k = 1, np
              call dwght( n, nq, we1, ldwe, ld2we,                            &
               fjacb(1: n, k,1: nq), tempret(1: n,1: nq))
              fjacb(1: n, k,1: nq) = tempret(1: n,1: nq)
20         continue
        else
           k1 = 0
           do 30 k = 1, np
              if ( ifixb( k) .ge. 1) then
                 k1 = k1+1
                 call dwght( n, nq, we1, ldwe, ld2we,                         &
                  fjacb(1: n, k,1: nq), tempret(1: n,1: nq))
                 fjacb(1: n, k1,1: nq) = tempret(1: n,1: nq)
              endif
30         continue
        endif
!
!  Weight the Jacobian's wrt DELTA as appropriate
!
        if ( isodr) then
           do 40 j = 1, m
              call dwght( n, nq, we1, ldwe, ld2we,                            &
               fjacd(1: n, j,1: nq), tempret(1: n,1: nq))
              fjacd(1: n, j,1: nq) = tempret(1: n,1: nq)
40         continue
        endif
!
        return
end subroutine
!DFCTR
subroutine dfctr( oksemi, a, lda, n, info)
!***Begin Prologue  DFCTR
!***Refer to  ODR
!***Routines Called  DDOT
!***Date Written   910706   (YYMMDD)
!***Revision Date  920619   (YYMMDD)
!***Purpose  Factor the positive (semi)definite matrix A using a
!       modified Cholesky factorization
!       (adapted from LINPACK subroutine DPOFA)
!***References  Dongarra J.J., Bunch J.R., Moler C.B., Stewart G.W.,
!       *LINPACK Users Guide*, SIAM, 1979.
!***End PROLOGUE  DFCTR
!
!...Used modules
        use odrpack_kinds,only: wp
!
!...Scalar arguments
        integer info, lda, n
        logical oksemi
!
!...Array arguments
        real(kind = wp) a( lda, n)
!
!...Local scalars
        real(kind = wp) xi, s, t, ten, zero
        integer j, k
!
!...External functions
        external ddot
        real(kind = wp) ddot
!
        data                                                                  &
         zero, ten                                                            &
         /0.0E0_wp,10.0E0_wp/
!
!...Variable Definitions (alphabetically)
!       A:       The array to be factored.  Upon return, A contains the
!       upper triangular matrix  R  so that  A = trans(R)*R
!       where the strict lower triangle is set to zero
!       if  INFO .NE. 0 , the factorization is not complete.
!       I:       An indexing variable.
!       INFO:    An idicator variable, where if
!       INFO = 0  then factorization was completed
!       INFO = K  signals an error condition.  The leading minor
!       of order  K  is not positive (semi)definite.
!       J:       An indexing variable.
!       LDA:     The leading dimension of array A.
!       N:       The number of rows and columns of data in array A.
!       OKSEMI:  The indicating whether the factored array can be positive
!       semidefinite (OKSEMI=TRUE) or whether it must be found to
!       be positive definite (OKSEMI=FALSE).
!       TEN:     The value 10.0E0_wp.
!       XI:      A value used to test for non positive semidefiniteness.
!       ZERO:    The value 0.0E0_wp.
!
!
!***First executable statement  DFCTR
!
!
!  Set relative tolerance for detecting non positive semidefiniteness.
        xi = - ten*epsilon( zero)
!------------------^-----------------------------------------------------------
!!! FPT - 3437 Mixed real or complex sizes in expression - loss of precision
!------------------------------------------------------------------------------
!
!  Compute factorization, storing in upper triangular portion of A
        do 20 j = 1, n
           info = j
           s = zero
           do 10 k = 1, j-1
              if ( a( k, k) .eq. zero) then
!--------------------------------^---------------------------------------------
!!! FPT - 3087 REAL or COMPLEX quantity tested for exact equality/inequality
!------------------------------------------------------------------------------
                 t = zero
              else
                 t = a( k, j)- ddot( k-1, a(1, k),1, a(1, j),1)
                 t = t/ a( k, k)
              endif
              a( k, j) = t
              s = s+ t* t
10         continue
           s = a( j, j)- s
!          ......Exit
           if ( a( j, j) .lt. zero .or. s .lt. xi*abs( a( j, j))) then
              return
           elseif (.not. oksemi .and. s .le. zero) then
              return
           elseif ( s .le. zero) then
              a( j, j) = zero
           else
              a( j, j) = sqrt( s)
           endif
20      continue
        info = 0
!
!  Zero out lower portion of A
        do 40 j = 2, n
           do 30 k = 1, j-1
              a( j, k) = zero
30         continue
40      continue
!
        return
end subroutine
!DFCTRW
subroutine dfctrw                                                             &
         ( n, m, nq, npp,                                                     &
         isodr,                                                               &
         we, ldwe, ld2we, wd, ldwd, ld2wd,                                    &
         wrk0, wrk4,                                                          &
         we1, nnzw, info)
!***Begin Prologue  DFCTRW
!***Refer to  ODR
!***Routines Called  DFCTR
!***Date Written   860529   (YYMMDD)
!***Revision Date  920619   (YYMMDD)
!***Purpose  Check input parameters, indicating errors found using
!       nonzero values of argument INFO as described in the
!       ODRPACK95 reference guide
!***End Prologue  DFCTRW
!
!...Used modules
        use odrpack_kinds,only: wp
!
!...Scalar arguments
        integer                                                               &
         info, ldwd, ldwe, ld2wd, ld2we,                                      &
         m, n, nnzw, npp, nq
        logical                                                               &
         isodr
!
!...Array arguments
        real(kind = wp)                                                       &
         we( ldwe, ld2we, nq), we1( ldwe, ld2we, nq), wd( ldwd, ld2wd, m),    &
         wrk0( nq, nq), wrk4( m, m)
!
!...Local scalars
        real(kind = wp)                                                       &
         zero
        integer                                                               &
         i, inf, j, j1, j2, l, l1, l2
        logical                                                               &
         notzro
!
!...External subroutines
        external                                                              &
         dfctr
!
!...Data statements
        data                                                                  &
         zero                                                                 &
         /0.0E0_wp/
!
!...Variable Definitions (alphabetically)
!       I:       An indexing variable.
!       INFO:    The variable designating why the computations were stopped.
!       ISODR:   The variable designating whether the solution is by ODR
!       (ISODR=TRUE) or by OLS (ISODR=FALSE).
!       J:       An indexing variable.
!       J1:      An indexing variable.
!       J2:      An indexing variable.
!       L:       An indexing variable.
!       L1:      An indexing variable.
!       L2:      An indexing variable.
!       LAST:    The last row of the array to be accessed.
!       LDWD:    The leading dimension of array WD.
!       LDWE:    The leading dimension of array WE.
!       LD2WD:   The second dimension of array WD.
!       LD2WE:   The second dimension of array WE.
!       M:       The number of columns of data in the explanatory variable.
!       N:       The number of observations.
!       NNZW:    The number of nonzero weighted observations.
!       NOTZRO:  The variable designating whether a given component of the
!       weight array WE contains a nonzero element (NOTZRO=FALSE)
!       or not (NOTZRO=TRUE).
!       NPP:     The number of function parameters being estimated.
!       NQ:      The number of responses per observations.
!       WE:      The (squared) EPSILON weights.
!       WE1:     The factored EPSILON weights, S.T. trans(WE1)*WE1 = WE.
!       WD:      The (squared) DELTA weights.
!       WRK0:    A work array of (NQ BY NQ) elements.
!       WRK4:    A work array of (M BY M) elements.
!       ZERO:    The value 0.0E0_wp.
!
!
!***First executable statement  DFCTRW
!
!
!  Check EPSILON weights, and store factorization in WE1
!
        if ( we(1,1,1) .lt. zero) then
!  WE contains a scalar
           we1(1,1,1) = -sqrt(abs( we(1,1,1)))
           nnzw = n
!
        else
           nnzw = 0
!
           if ( ldwe .eq. 1) then
!
              if ( ld2we .eq. 1) then
!  WE contains a diagonal matrix
                 do 110 l = 1, nq
                    if ( we(1,1, l) .gt. zero) then
                       nnzw = n
                       we1(1,1, l) = sqrt( we(1,1, l))
                    elseif ( we(1,1, l) .lt. zero) then
                       info = 30010
                       goto 300
                    endif
110              continue
              else
!
!  WE contains a full NQ by NQ semidefinite matrix
                 do 130 l1 = 1, nq
                    do 120 l2 = l1, nq
                       wrk0( l1, l2) = we(1, l1, l2)
120                 continue
130              continue
                 call dfctr(.true., wrk0, nq, nq, inf)
                 if ( inf .ne. 0) then
                    info = 30010
                    goto 300
                 else
                    do 150 l1 = 1, nq
                       do 140 l2 = 1, nq
                          we1(1, l1, l2) = wrk0( l1, l2)
140                    continue
                       if ( we1(1, l1, l1) .ne. zero) then
!-----------------------------------------------^------------------------------
!!! FPT - 3087 REAL or COMPLEX quantity tested for exact equality/inequality
!------------------------------------------------------------------------------
                          nnzw = n
                       endif
150                 continue
                 endif
              endif
!
           else
!
              if ( ld2we .eq. 1) then
!  WE contains an array of  diagonal matrix
                 do 220 i = 1, n
                    notzro = .false.
                    do 210 l = 1, nq
                       if ( we( i,1, l) .gt. zero) then
                          notzro = .true.
                          we1( i,1, l) = sqrt( we( i,1, l))
                       elseif ( we( i,1, l) .lt. zero) then
                          info = 30010
                          goto 300
                       endif
210                 continue
                    if ( notzro) then
                       nnzw = nnzw+1
                    endif
220              continue
              else
!
!  WE contains an array of full NQ by NQ semidefinite matrices
                 do 270 i = 1, n
                    do 240 l1 = 1, nq
                       do 230 l2 = l1, nq
                          wrk0( l1, l2) = we( i, l1, l2)
230                    continue
240                 continue
                    call dfctr(.true., wrk0, nq, nq, inf)
                    if ( inf .ne. 0) then
                       info = 30010
                       goto 300
                    else
                       notzro = .false.
                       do 260 l1 = 1, nq
                          do 250 l2 = 1, nq
                             we1( i, l1, l2) = wrk0( l1, l2)
250                       continue
                          if ( we1( i, l1, l1) .ne. zero) then
!---------------------------------------------------^--------------------------
!!! FPT - 3087 REAL or COMPLEX quantity tested for exact equality/inequality
!------------------------------------------------------------------------------
                             notzro = .true.
                          endif
260                    continue
                    endif
                    if ( notzro) then
                       nnzw = nnzw+1
                    endif
270              continue
              endif
           endif
        endif
!
!  Check for a sufficient number of nonzero EPSILON weights
!
        if ( nnzw .lt. npp) then
           info = 30020
        endif
!
!
!  Check DELTA weights
!
300     continue
        if (.not. isodr .or. wd(1,1,1) .lt. zero) then
!  Problem is not ODR, or WD contains a scalar
           return
!
        else
!
           if ( ldwd .eq. 1) then
!
              if ( ld2wd .eq. 1) then
!  WD contains a diagonal matrix
                 do 310 j = 1, m
                    if ( wd(1,1, j) .le. zero) then
                       info = max(30001, info+1)
                       return
                    endif
310              continue
              else
!
!  WD contains a full M by M positive definite matrix
                 do 330 j1 = 1, m
                    do 320 j2 = j1, m
                       wrk4( j1, j2) = wd(1, j1, j2)
320                 continue
330              continue
                 call dfctr(.false., wrk4, m, m, inf)
                 if ( inf .ne. 0) then
                    info = max(30001, info+1)
                    return
                 endif
              endif
!
           else
!
              if ( ld2wd .eq. 1) then
!  WD contains an array of diagonal matrices
                 do 420 i = 1, n
                    do 410 j = 1, m
                       if ( wd( i,1, j) .le. zero) then
                          info = max(30001, info+1)
                          return
                       endif
410                 continue
420              continue
              else
!
!  WD contains an array of full M by M positive definite matrices
                 do 470 i = 1, n
                    do 440 j1 = 1, m
                       do 430 j2 = j1, m
                          wrk4( j1, j2) = wd( i, j1, j2)
430                    continue
440                 continue
                    call dfctr(.false., wrk4, m, m, inf)
                    if ( inf .ne. 0) then
                       info = max(30001, info+1)
                       return
                    endif
470              continue
              endif
           endif
        endif
!
        return
end subroutine
!DFLAGS
subroutine dflags                                                             &
         ( job, restrt, initd, dovcv, redoj, anajac, cdjac, chkjac, isodr,    &
         implct)
!***Begin Prologue  DFLAGS
!***Refer to  ODR
!***Routines Called  (None)
!***Date Written   860529   (YYMMDD)
!***Revision Date  920304   (YYMMDD)
!***Purpose  Set flags indicating conditions specified by JOB
!***End Prologue  DFLAGS
!
!...Scalar arguments
        integer                                                               &
         job
        logical                                                               &
         anajac, cdjac, chkjac, dovcv, implct, initd, isodr, redoj, restrt
!
!...Local scalars
        integer                                                               &
         j
!
!...Variable Definitions (alphabetically)
!       ANAJAC:  The variable designating whether the Jacobians are computed
!       by finite differences (ANAJAC=FALSE) or not (ANAJAC=TRUE).
!       CDJAC:   The variable designating whether the Jacobians are computed
!       by central differences (CDJAC=TRUE) or by forward
!       differences (CDJAC=FALSE).
!       CHKJAC:  The variable designating whether the user-supplied
!       Jacobians are to be checked (CHKJAC=TRUE) or not
!       (CHKJAC=FALSE).
!       DOVCV:   The variable designating whether the covariance matrix is
!       to be computed (DOVCV=TRUE) or not (DOVCV=FALSE).
!       IMPLCT:  The variable designating whether the solution is by
!       implicit ODR (IMPLCT=TRUE) or explicit ODR (IMPLCT=FALSE).
!       INITD:   The variable designating whether DELTA is to be initialized
!       to zero (INITD=TRUE) or to the first N by M elements of
!       array WORK (INITD=FALSE).
!       ISODR:   The variable designating whether the solution is by ODR
!       (ISODR=TRUE) or by OLS (ISODR=FALSE).
!       J:       The value of a specific digit of JOB.
!       JOB:     The variable controling problem initialization and
!       computational method.
!       REDOJ:   The variable designating whether the Jacobian matrix is to
!       be recomputed for the computation of the covariance matrix
!       (REDOJ=TRUE) or not (REDOJ=FALSE).
!       RESTRT:  The variable designating whether the call is a restart
!       (RESTRT=TRUE) or not (RESTRT=FALSE).
!
!
!***First executable statement  DFLAGS
!
!
        if ( job .ge. 0) then
!
           restrt = job .ge. 10000
!
           initd = mod( job,10000)/1000 .eq. 0
!
           j = mod( job,1000)/100
           if ( j .eq. 0) then
              dovcv = .true.
              redoj = .true.
           elseif ( j .eq. 1) then
              dovcv = .true.
              redoj = .false.
           else
              dovcv = .false.
              redoj = .false.
           endif
!
           j = mod( job,100)/10
           if ( j .eq. 0) then
              anajac = .false.
              cdjac = .false.
              chkjac = .false.
           elseif ( j .eq. 1) then
              anajac = .false.
              cdjac = .true.
              chkjac = .false.
           elseif ( j .eq. 2) then
              anajac = .true.
              cdjac = .false.
              chkjac = .true.
           else
              anajac = .true.
              cdjac = .false.
              chkjac = .false.
           endif
!
           j = mod( job,10)
           if ( j .eq. 0) then
              isodr = .true.
              implct = .false.
           elseif ( j .eq. 1) then
              isodr = .true.
              implct = .true.
           else
              isodr = .false.
              implct = .false.
           endif
!
        else
!
           restrt = .false.
           initd = .true.
           dovcv = .true.
           redoj = .true.
           anajac = .false.
           cdjac = .false.
           chkjac = .false.
           isodr = .true.
           implct = .false.
!
        endif
!
        return
end subroutine
!DHSTEP
function dhstep                                                               &
         ( itype, neta, i, j, stp, ldstp)                                     &
         result( dhstepr)
!***Begin Prologue  DHSTEP
!***Refer to  ODR
!***Routines Called  (NONE)
!***Date Written   860529   (YYMMDD)
!***Revision Date  920304   (YYMMDD)
!***Purpose  Set relative step size for finite difference derivatives
!***End Prologue  DHSTEP
!
!...Used modules
        use odrpack_kinds,only: wp
!
!...Scalar arguments
        integer                                                               &
         i, itype, j, ldstp, neta
!
!...Array arguments
        real(kind = wp)                                                       &
         stp( ldstp, j)
!
!...Result
        real(kind = wp)                                                       &
         dhstepr
!
!...Local scalars
        real(kind = wp)                                                       &
         ten, three, two, zero
!
!...Data statements
        data                                                                  &
         zero, two, three, ten                                                &
         /0.0E0_wp,2.0E0_wp,3.0E0_wp,10.0E0_wp/
!
!...Variable Definitions (alphabetically)
!       I:       An identifier for selecting user supplied step sizes.
!       ITYPE:   The finite difference method being used, where
!       ITYPE = 0 indicates forward finite differences, and
!       ITYPE = 1 indicates central finite differences.
!       J:       An identifier for selecting user supplied step sizes.
!       LDSTP:   The leading dimension of array STP.
!       NETA:    The number of good digits in the function results.
!       STP:     The step size for the finite difference derivative.
!       TEN:     The value 10.0E0_wp.
!       THREE:   The value 3.0E0_wp.
!       TWO:     The value 2.0E0_wp.
!       ZERO:    The value 0.0E0_wp.
!
!
!
!***First executable statement  DHSTEP
!
!
!  Set DHSTEP to relative finite difference step size
!
        if ( stp(1,1) .le. zero) then
!
           if ( itype .eq. 0) then
!  Use default forward finite difference step size
              dhstepr = ten**(-abs( neta)/ two- two)
!
           else
!  Use default central finite difference step size
              dhstepr = ten**(-abs( neta)/ three)
           endif
!
        elseif ( ldstp .eq. 1) then
           dhstepr = stp(1, j)
!
        else
           dhstepr = stp( i, j)
        endif
!
        return
end function
!DIFIX
subroutine difix                                                              &
         ( n, m, ifix, ldifix, t, ldt, tfix, ldtfix)
!***Begin Prologue  DIFIX
!***Refer to  ODR
!***Routines Called  (None)
!***Date Written   910612   (YYMMDD)
!***Revision Date  920304   (YYMMDD)
!***Purpose  Set elements of T to zero according to IFIX
!***End Prologue  DIFIX
!
!...Used modules
        use odrpack_kinds,only: wp
!
!...Scalar arguments
        integer                                                               &
         ldifix, ldt, ldtfix, m, n
!
!...Array arguments
        real(kind = wp)                                                       &
         t( ldt, m), tfix( ldtfix, m)
        integer                                                               &
         ifix( ldifix, m)
!------------^-----------------------------------------------------------------
!!! FPT - 2517 ANSI FORTRAN 77 intrinsic used as a local identifier
!------------------------------------------------------------------------------
!
!...Local scalars
        real(kind = wp)                                                       &
         zero
        integer                                                               &
         i, j
!
!...Data statements
        data                                                                  &
         zero                                                                 &
         /0.0E0_wp/
!
!...Variable Definitions (alphabetically)
!       I:       An indexing variable.
!       IFIX:    The array designating whether an element of T is to be
!       set to zero.
!       J:       an indexing variable.
!       LDT:     The leading dimension of array T.
!       LDIFIX:  The leading dimension of array IFIX.
!       LDTFIX:  The leading dimension of array TFIX.
!       M:       The number of columns of data in the array.
!       N:       The number of rows of data in the array.
!       T:       The array being set to zero according to the elements
!       of IFIX.
!       TFIX:    The resulting array.
!       ZERO:    The value 0.0E0_wp.
!
!
!***First executable statement  DIFIX
!
!
        if ( n .eq. 0 .or. m .eq. 0) return
!
        if ( ifix(1,1) .ge. zero) then
!---------------------------^--------------------------------------------------
!!! FPT - 3085 Objects of .EQ., .NE. .GT. etc. are of different data types
!------------------------------------------------------------------------------
           if ( ldifix .ge. n) then
              do 20 j = 1, m
                 do 10 i = 1, n
                    if ( ifix( i, j) .eq. 0) then
                       tfix( i, j) = zero
                    else
                       tfix( i, j) = t( i, j)
                    endif
10               continue
20            continue
           else
              do 100 j = 1, m
                 if ( ifix(1, j) .eq. 0) then
                    do 30 i = 1, n
                       tfix( i, j) = zero
30                  continue
                 else
                    do 90 i = 1, n
                       tfix( i, j) = t( i, j)
90                  continue
                 endif
100           continue
           endif
        endif
!
        return
end subroutine
!DINIWK
subroutine diniwk                                                             &
         ( n, m, np, work, lwork, iwork, liwork,                              &
         x, ldx, ifixx, ldifx, scld, ldscld,                                  &
         beta, sclb,                                                          &
         sstol, partol, maxit, taufac,                                        &
         job, iprint, lunerr, lunrpt,                                         &
         lower, upper,                                                        &
         epsmai, sstoli, partli, maxiti, taufci,                              &
         jobi, iprini, luneri, lunrpi,                                        &
         ssfi, tti, ldtti, deltai,                                            &
         loweri, upperi, boundi)
!***Begin Prologue  DINIWK
!***Refer to  ODR
!***Routines Called  DFLAGS,DSCLB,DSCLD
!***Date Written   860529   (YYMMDD)
!***Revision Date  920304   (YYMMDD)
!***Purpose  Initialize work vectors as necessary
!***End Prologue  DINIWK
!
!...Used modules
        use odrpack_kinds, only: wp, zero, one, two, three

!...Scalar arguments
        real(kind = wp)                                                       &
         partol, sstol, taufac
        integer                                                               &
         boundi, deltai, epsmai, iprini, iprint, job, jobi, ldifx,            &
         ldscld, ldtti, ldx, liwork, loweri, luneri, lunerr, lunrpi, lunrpt,  &
         lwork, m, maxit, maxiti, n, np, partli, ssfi, sstoli, taufci, tti,   &
         upperi

!...Array arguments
        real(kind = wp)                                                       &
         beta( np), lower( np), sclb( np), scld( ldscld, m), upper( np),      &
         work( lwork), x( ldx, m)
        integer                                                               &
         ifixx( ldifx, m), iwork( liwork)

!...Local scalars
        integer :: i, j, istart
        logical :: anajac, cdjac, chkjac, dovcv, implct, initd, isodr, redoj, restrt

!...External subroutines
        external :: dcopy, dflags, dsclb, dscld


!...Variable Definitions (alphabetically)
!       ANAJAC:  The variable designating whether the Jacobians are
!       computed by finite differences (ANAJAC=FALSE) or not
!       (ANAJAC=TRUE).
!       BETA:    The function parameters.
!       CDJAC:   The variable designating whether the Jacobians are
!       computed by central differences (CDJAC=TRUE) or by forward
!       differences (CDJAC=FALSE).
!       CHKJAC:  The variable designating whether the user-supplied
!       Jacobians are to be checked (CHKJAC=TRUE) or not
!       (CHKJAC=FALSE).
!       DELTAI:  The starting location in array WORK of array DELTA.
!       DOVCV:   The variable designating whether the covariance matrix is
!       to be computed (DOVCV=TRUE) or not (DOVCV=FALSE).
!       EPSMAI:  The location in array WORK of variable EPSMAC.
!       I:       An indexing variable.
!       IFIXX:   The values designating whether the elements of X are fixed
!       at their input values or not.
!       IMPLCT:  The variable designating whether the solution is by
!       implicit ODR (IMPLCT=TRUE) or explicit ODR (IMPLCT=FALSE).
!       INITD:   The variable designating whether DELTA is to be initialized
!       to zero (INITD=TRUE) or to the values in the first N by M
!       elements of array WORK (INITD=FALSE).
!       IPRINI:  The location in array IWORK of variable IPRINT.
!       IPRINT:  The print control variable.
!       ISODR:   The variable designating whether the solution is by ODR
!       (ISODR=TRUE) or by OLS (ISODR=FALSE).
!       IWORK:   The integer work space.
!       J:       An indexing variable.
!       JOB:     The variable controling problem initialization and
!       computational method.
!       JOBI:    The location in array IWORK of variable JOB.
!       LDIFX:   The leading dimension of array IFIXX.
!       LDSCLD:  The leading dimension of array SCLD.
!       LDTTI:   The leading dimension of array TT.
!       LDX:     The leading dimension of array X.
!       LIWORK:  The length of vector IWORK.
!       LUNERI:  The location in array IWORK of variable LUNERR.
!       LUNERR:  The logical unit number used for error messages.
!       LUNRPI:  The location in array iwork of variable LUNRPT.
!       LUNRPT:  The logical unit number used for computation reports.
!       LWORK:   The length of vector WORK.
!       M:       The number of columns of data in the independent variable.
!       MAXIT:   The maximum number of iterations allowed.
!       MAXITI:  The location in array IWORK of variable MAXIT.
!       N:       The number of observations.
!       NP:      The number of function parameters.
!       ONE:     The value 1.0E0_wp.
!       PARTLI:  The location in array work of variable partol.
!       PARTOL:  The parameter convergence stopping criteria.
!       REDOJ:   The variable designating whether the Jacobian matrix is to
!       be recomputed for the computation of the covariance matrix
!       (REDOJ=TRUE) or not (REDOJ=FALSE).
!       RESTRT:  The variable designating whether the call is a restart
!       (RESTRT=TRUE) or not (RESTRT=FALSE).
!       SCLB:    The scaling values for BETA.
!       SCLD:    The scaling values for DELTA.
!       SSFI:    The starting location in array WORK of array SSF.
!       SSTOL:   The sum-of-squares convergence stopping criteria.
!       SSTOLI:  The location in array WORK of variable SSTOL.
!       TAUFAC:  The factor used to compute the initial trust region
!       diameter.
!       TAUFCI:  The location in array WORK of variable TAUFAC.
!       THREE:   The value 3.0E0_wp.
!       TTI:     The starting location in array WORK of the ARRAY TT.
!       TWO:     The value 2.0E0_wp.
!       WORK:    The REAL (KIND=wp) work space.
!       X:       The independent variable.
!       ZERO:    The value 0.0E0_wp.
!
!
!***First executable statement  DINIWK
!
!
        call dflags( job, restrt, initd, dovcv, redoj,                        &
         anajac, cdjac, chkjac, isodr, implct)
!
!  Store value of machine precision in work vector
!
        work( epsmai) = epsilon( zero)
!
!  Set tolerance for stopping criteria based on the change in the
!  parameters  (see also subprogram DODCNT)
!
        if ( partol .lt. zero) then
           work( partli) = work( epsmai)**( two/ three)
        else
           work( partli) = min( partol, one)
        endif
!
!  Set tolerance for stopping criteria based on the change in the
!  sum of squares of the weighted observational errors
!
        if ( sstol .lt. zero) then
           work( sstoli) = sqrt( work( epsmai))
        else
           work( sstoli) = min( sstol, one)
        endif
!
!  Set factor for computing trust region diameter at first iteration
!
        if ( taufac .le. zero) then
           work( taufci) = one
        else
           work( taufci) = min( taufac, one)
        endif
!
!  Set maximum number of iterations
!
        if ( maxit .lt. 0) then
           iwork( maxiti) = 50
        else
           iwork( maxiti) = maxit
        endif
!
!  Store problem initialization and computational method control
!  variable
!
        if ( job .le. 0) then
           iwork( jobi) = 0
        else
           iwork( jobi) = job
        endif
!
!  Set print control
!
        if ( iprint .lt. 0) then
           iwork( iprini) = 2001
        else
           iwork( iprini) = iprint
        endif
!
!  Set logical unit number for error messages
!
        if ( lunerr .lt. 0) then
           iwork( luneri) = 6
        else
           iwork( luneri) = lunerr
        endif
!
!  Set logical unit number for computation reports
!
        if ( lunrpt .lt. 0) then
           iwork( lunrpi) = 6
        else
           iwork( lunrpi) = lunrpt
        endif
!
!  Compute scaling for BETA's and DELTA's
!
        if ( sclb(1) .le. zero) then
           call dsclb( np, beta, work( ssfi))
        else
           call dcopy( np, sclb,1, work( ssfi),1)
        endif
        if ( isodr) then
           if ( scld(1,1) .le. zero) then
              iwork( ldtti) = n
              call dscld( n, m, x, ldx, work( tti), iwork( ldtti))
           else
              if ( ldscld .eq. 1) then
                 iwork( ldtti) = 1
                 call dcopy( m, scld(1,1),1, work( tti),1)
              else
                 iwork( ldtti) = n
                 do 10 j = 1, m
                    call dcopy( n, scld(1, j),1,                              &
                     work( tti+( j-1)* iwork( ldtti)),1)
10               continue
              endif
           endif
        endif
!
!  Initialize DELTA's as necessary
!
        if ( isodr) then
           if (initd) then
              !call dzero( n, m, work( deltai), n)
              work(deltai:deltai+(n*m-1)) = zero
           else
              if ( ifixx(1,1) .ge. 0) then
                 if ( ldifx .eq. 1) then
                    do 20 j = 1, m
                       if ( ifixx(1, j) .eq. 0) then
                          istart = deltai+(j-1)*n
                          work(istart:istart+(n-1)) = zero
                          !call dzero( n,1, work( deltai+( j-1)* n), n)
                       endif
20                  continue
                 else
                    do 40 j = 1, m
                       do 30 i = 1, n
                          if ( ifixx( i, j) .eq. 0) then
                             work( deltai-1+ i+( j-1)* n) = zero
                          endif
30                     continue
40                  continue
                 endif
              endif
           endif
        else
           !call dzero( n, m, work( deltai), n)
           work(deltai:deltai+(n*m-1)) = zero
        endif
!
!  Copy bounds into WORK
!
        work( loweri: loweri+ np-1) = lower(1: np)
        work( upperi: upperi+ np-1) = upper(1: np)
!
!  Initialize parameters on bounds in IWORK
!
        iwork( boundi: boundi+ np-1) = 0

end subroutine

subroutine diwinf                                                             &
         ( m, np, nq,                                                         &
         msgbi, msgdi, ifix2i, istopi,                                        &
         nnzwi, nppi, idfi,                                                   &
         jobi, iprini, luneri, lunrpi,                                        &
         nrowi, ntoli, netai,                                                 &
         maxiti, niteri, nfevi, njevi, int2i, iranki, ldtti,                  &
         boundi,                                                              &
         liwkmn)
!***Begin Prologue  DIWINF
!***Refer to  ODR
!***Routines Called  (NONE)
!***Date Written   860529   (YYMMDD)
!***Revision Date  920304   (YYMMDD)
!***Purpose  Set storage locations within integer work space
!***End Prologue  DIWINF
!
!...Scalar arguments
        integer                                                               &
         boundi, idfi, int2i, iprini, iranki, istopi, jobi, ifix2i, ldtti,    &
         liwkmn, luneri, lunrpi, m, maxiti, msgbi, msgdi, netai, nfevi,       &
         niteri, njevi, nnzwi, np, nppi, nq, nrowi, ntoli
!
!...Variable Definitions (alphabetically)
!       IDFI:    The location in array IWORK of variable IDF.
!       IFIX2I:  The starting location in array IWORK of array IFIX2.
!       INT2I:   The location in array IWORK of variable INT2.
!       IPRINI:  The location in array IWORK of variable IPRINT.
!       IRANKI:  The location in array IWORK of variable IRANK.
!       ISTOPI:  The location in array IWORK of variable ISTOP.
!       JOBI:    The location in array IWORK of variable JOB.
!       LDTTI:   The location in array IWORK of variable LDTT.
!       LIWKMN:  The minimum acceptable length of array IWORK.
!       LUNERI:  The location in array IWORK of variable LUNERR.
!       LUNRPI:  The location in array IWORK of variable LUNRPT.
!       M:       The number of columns of data in the independent variable.
!       MAXITI:  The location in array iwork of variable MAXIT.
!       MSGBI:   The starting location in array IWORK of array MSGB.
!       MSGDI:   The starting location in array IWORK of array MSGD.
!       NETAI:   The location in array IWORK of variable NETA.
!       NFEVI:   The location in array IWORK of variable NFEV.
!       NITERI:  The location in array IWORK of variabel NITER.
!       NJEVI:   The location in array IWORK of variable NJEV.
!       NNZWI:   The location in array IWORK of variable NNZW.
!       NP:      The number of function parameters.
!       NPPI:    The location in array IWORK of variable NPP.
!       NQ:      The number of responses per observation.
!       NROWI:   The location in array IWORK of variable NROW.
!       NTOLI:   The location in array IWORK of variable NTOL.
!
!
!***First executable statement  DIWINF
!
!
        if ( np .ge. 1 .and. m .ge. 1) then
           msgbi = 1
           msgdi = msgbi+ nq* np+1
           ifix2i = msgdi+ nq* m+1
           istopi = ifix2i+ np
           nnzwi = istopi+1
           nppi = nnzwi+1
           idfi = nppi+1
           jobi = idfi+1
           iprini = jobi+1
           luneri = iprini+1
           lunrpi = luneri+1
           nrowi = lunrpi+1
           ntoli = nrowi+1
           netai = ntoli+1
           maxiti = netai+1
           niteri = maxiti+1
           nfevi = niteri+1
           njevi = nfevi+1
           int2i = njevi+1
           iranki = int2i+1
           ldtti = iranki+1
           boundi = ldtti+1
           liwkmn = boundi+ np-1
        else
           msgbi = 1
           msgdi = 1
           ifix2i = 1
           istopi = 1
           nnzwi = 1
           nppi = 1
           idfi = 1
           jobi = 1
           iprini = 1
           luneri = 1
           lunrpi = 1
           nrowi = 1
           ntoli = 1
           netai = 1
           maxiti = 1
           niteri = 1
           nfevi = 1
           njevi = 1
           int2i = 1
           iranki = 1
           ldtti = 1
           boundi = 1
           liwkmn = 1
        endif
!
        return
end subroutine
!DJACCD
subroutine djaccd                                                             &
         ( fcn,                                                               &
         n, m, np, nq,                                                        &
         beta, x, ldx, delta, xplusd, ifixb, ifixx, ldifx,                    &
         stpb, stpd, ldstpd,                                                  &
         ssf, tt, ldtt, neta, fn, stp, wrk1, wrk2, wrk3, wrk6,                &
         fjacb, isodr, fjacd, nfev, istop, info,                              &
         lower, upper)
!***Begin Prologue  DJACCD
!***Refer to  ODR
!***Routines Called  FCN,DHSTEP
!***Date Written   860529   (YYMMDD)
!***Revision Date  920619   (YYMMDD)
!***Purpose  Compute central difference approximations to the
!       Jacobian wrt the estimated BETAS and wrt the DELTAS
!***End Prologue  DJACCD
!
!...Used modules
        use odrpack_kinds,only: wp, zero, one
!
!...Scalar arguments
        integer                                                               &
         info, istop, ldifx, ldstpd, ldtt, ldx, m, n, neta, nfev, np, nq
        logical                                                               &
         isodr
!
!...Array arguments
        real(kind = wp)                                                       &
         beta( np), delta( n, m), fjacb( n, np, nq), fjacd( n, m, nq), fn( n, &
         nq), lower( np), ssf( np), stp( n), stpb( np), stpd( ldstpd, m), tt( &
         ldtt, m), upper( np), wrk1( n, m, nq), wrk2( n, nq), wrk3( np), wrk6 &
         ( n, np, nq), x( ldx, m), xplusd( n, m)
        integer                                                               &
         ifixb( np), ifixx( ldifx, m)
!
!...Subroutine arguments
        external                                                              &
         fcn

!...Local scalars
        real(kind = wp) :: betak, typj
        integer :: i, j, k, l
        logical :: doit, setzro

!...External functions
        real(kind = wp), external :: dhstep, derstep

!...Routine names used as subprogram arguments
!       FCN:     The user supplied subroutine for evaluating the model.
!
!...Variable Definitions (alphabetically)
!       BETA:    The function parameters.
!       BETAK:   The K-th function parameter.
!       DELTA:   The estimated errors in the explanatory variables.
!       DOIT:    The variable designating whether the derivative wrt a given
!       BETA or DELTA needs to be computed (DOIT=TRUE) or not
!       (DOIT=FALSE).
!       FJACB:   The Jacobian with respect to BETA.
!       FJACD:   The Jacobian with respect to DELTA.
!       FN:      The new predicted values from the function.  Used when parameter is
!       on a boundary.
!       I:       An indexing variable.
!       IFIXB:   The values designating whether the elements of BETA are
!       fixed at their input values or not.
!       IFIXX:   The values designating whether the elements of X are fixed
!       at their input values or not.
!       INFO:    The variable designating why the computations were stopped.
!       ISODR:   The variable designating whether the solution is by ODR
!       (ISODR=TRUE) or by OLS (ISODR=FALSE).
!       ISTOP:   The variable designating whether there are problems
!       computing the function at the current BETA and DELTA.
!       J:       An indexing variable.
!       K:       An indexing variable.
!       L:       An indexing variable.
!       LDIFX:   The leading dimension of array IFIXX.
!       LDSTPD:  The leading dimension of array STPD.
!       LDTT:    The leading dimension of array TT.
!       LDX:     The leading dimension of array X.
!       LOWER:   The lower bound on BETA.
!       M:       The number of columns of data in the explanatory variable.
!       N:       The number of observations.
!       NETA:    The number of good digits in the function results.
!       NFEV:    The number of function evaluations.
!       NP:      The number of function parameters.
!       SETZRO:  The variable designating whether the derivative wrt some
!       DELTA needs to be set to zero (SETZRO=TRUE) or not
!       (SETZRO=FALSE).
!       SSF:     The scaling values used for BETA.
!       STP:     The step used for computing finite difference
!       derivatives with respect to each DELTA.
!       STPB:    the relative step used for computing finite difference
!       derivatives with respect to each BETA.
!       STPD:    The relative step used for computing finite difference
!       derivatives with respect to each DELTA.
!       TT:      The scaling values used for DELTA.
!       TYPJ:    The typical size of the J-th unknown BETA or DELTA.
!       UPPER:   The upper bound on BETA.
!       X:       The explanatory variable.
!       XPLUSD:  The values of X + DELTA.
!       WRK1:    A work array of (N BY M BY NQ) elements.
!       WRK2:    A work array of (N BY NQ) elements.
!       WRK3:    A work array of (NP) elements.
!       WRK6:    A WORK ARRAY OF (N BY NP BY NQ) elements.
!
!
!***First executable statement  DJACCD
!
!
!  Compute the Jacobian wrt the estimated BETAS
!
        do 60 k = 1, np
           if ( ifixb(1) .ge. 0) then
              if ( ifixb( k) .eq. 0) then
                 doit = .false.
              else
                 doit = .true.
              endif
           else
              doit = .true.
           endif
           if (.not. doit) then
              do 10 l = 1, nq
                 fjacb(1:n, k, l) = zero
10            continue
           else
              betak = beta( k)
              wrk3( k) = betak                                                &
               + derstep(1, k, betak, ssf, stpb, neta)
              wrk3( k) = wrk3( k)- betak
!
              beta( k) = betak+ wrk3( k)
              if ( beta( k) .gt. upper( k)) then
                 beta( k) = upper( k)
              elseif ( beta( k) .lt. lower( k)) then
                 beta( k) = lower( k)
              endif
              if ( beta( k)-2* wrk3( k) .lt. lower( k)) then
                 beta( k) = lower( k)+2* wrk3( k)
              elseif ( beta( k)-2* wrk3( k) .gt. upper( k)) then
                 beta( k) = upper( k)+2* wrk3( k)
              endif
              if ( beta( k) .gt. upper( k) .or. beta( k) .lt. lower( k)) then
                 info = 60001
                 return
              endif
              istop = 0
              if ( beta( k) .eq. betak) then
!--------------------------------^---------------------------------------------
!!! FPT - 3087 REAL or COMPLEX quantity tested for exact equality/inequality
!------------------------------------------------------------------------------
                 wrk2(1: n,1: nq) = fn(1: n,1: nq)
              else
                 call fcn( n, m, np, nq,                                      &
                  n, m, np,                                                   &
                  beta, xplusd,                                               &
                  ifixb, ifixx, ldifx,                                        &
                  001, wrk2, wrk6, wrk1,                                      &
                  istop)
                 if ( istop .ne. 0) then
                    return
                 else
                    nfev = nfev+1
                 endif
              endif
              do 30 l = 1, nq
                 do 20 i = 1, n
                    fjacb( i, k, l) = wrk2( i, l)
20               continue
30            continue
!
              beta( k) = beta( k)-2* wrk3( k)
              if ( beta( k) .gt. upper( k)) then
                 info = 60001
                 return
              endif
              if ( beta( k) .lt. lower( k)) then
                 info = 60001
                 return
              endif
              istop = 0
              if ( beta( k) .eq. betak) then
!--------------------------------^---------------------------------------------
!!! FPT - 3087 REAL or COMPLEX quantity tested for exact equality/inequality
!------------------------------------------------------------------------------
                 wrk2(1: n,1: nq) = fn(1: n,1: nq)
              else
                 call fcn( n, m, np, nq,                                      &
                  n, m, np,                                                   &
                  beta, xplusd,                                               &
                  ifixb, ifixx, ldifx,                                        &
                  001, wrk2, wrk6, wrk1,                                      &
                  istop)
                 if ( istop .ne. 0) then
                    return
                 else
                    nfev = nfev+1
                 endif
              endif
!
              do 50 l = 1, nq
                 do 40 i = 1, n
                    fjacb( i, k, l) = ( fjacb( i, k, l)- wrk2( i, l))/(2*     &
                     wrk3( k))
40               continue
50            continue
              beta( k) = betak
           endif
60      continue
!
!  Compute the Jacobian wrt the X'S
!
        if ( isodr) then
           do 220 j = 1, m
              if ( ifixx(1,1) .lt. 0) then
                 doit = .true.
                 setzro = .false.
              elseif ( ldifx .eq. 1) then
                 if ( ifixx(1, j) .eq. 0) then
                    doit = .false.
                 else
                    doit = .true.
                 endif
                 setzro = .false.
              else
                 doit = .false.
                 setzro = .false.
                 do 100 i = 1, n
                    if ( ifixx( i, j) .ne. 0) then
                       doit = .true.
                    else
                       setzro = .true.
                    endif
100              continue
              endif
              if (.not. doit) then
                 do 110 l = 1, nq
                    fjacd(1:n, j, l) = zero
110              continue
              else
                 do 120 i = 1, n
                    if ( xplusd( i, j) .eq. zero) then
!-------------------------------------------^----------------------------------
!!! FPT - 3087 REAL or COMPLEX quantity tested for exact equality/inequality
!------------------------------------------------------------------------------
                       if ( tt(1,1) .lt. zero) then
                          typj = one/abs( tt(1,1))
                       elseif ( ldtt .eq. 1) then
                          typj = one/ tt(1, j)
                       else
                          typj = one/ tt( i, j)
                       endif
                    else
                       typj = abs( xplusd( i, j))
                    endif
                    stp( i) = xplusd( i, j)                                   &
                     +sign( one, xplusd( i, j))                               &
                     * typj* dhstep(1, neta, i, j, stpd, ldstpd)
                    stp( i) = stp( i)- xplusd( i, j)
                    xplusd( i, j) = xplusd( i, j)+ stp( i)
120              continue
                 istop = 0
                 call fcn( n, m, np, nq,                                      &
                  n, m, np,                                                   &
                  beta, xplusd,                                               &
                  ifixb, ifixx, ldifx,                                        &
                  001, wrk2, wrk6, wrk1,                                      &
                  istop)
                 if ( istop .ne. 0) then
                    return
                 else
                    nfev = nfev+1
                    do 140 l = 1, nq
                       do 130 i = 1, n
                          fjacd( i, j, l) = wrk2( i, l)
130                    continue
140                 continue
                 endif
!
                 do 150 i = 1, n
                    xplusd( i, j) = x( i, j)+ delta( i, j)- stp( i)
150              continue
                 istop = 0
                 call fcn( n, m, np, nq,                                      &
                  n, m, np,                                                   &
                  beta, xplusd,                                               &
                  ifixb, ifixx, ldifx,                                        &
                  001, wrk2, wrk6, wrk1,                                      &
                  istop)
                 if ( istop .ne. 0) then
                    return
                 else
                    nfev = nfev+1
                 endif
!
                 if ( setzro) then
                    do 180 i = 1, n
                       if ( ifixx( i, j) .eq. 0) then
                          do 160 l = 1, nq
                             fjacd( i, j, l) = zero
160                       continue
                       else
                          do 170 l = 1, nq
                             fjacd( i, j, l) = ( fjacd( i, j, l)- wrk2( i, l) &
                              )/(2* stp( i))
170                       continue
                       endif
180                 continue
                 else
                    do 200 l = 1, nq
                       do 190 i = 1, n
                          fjacd( i, j, l) = ( fjacd( i, j, l)- wrk2( i, l))/  &
                           (2* stp( i))
190                    continue
200                 continue
                 endif
                 do 210 i = 1, n
                    xplusd( i, j) = x( i, j)+ delta( i, j)
210              continue
              endif
220        continue
        endif
!
        return
end subroutine
!MBFB
subroutine mbfb                                                               &
         ( np, beta, lower, upper, ssf, stpb, neta, eta, interval)
!***BEGIN PROLOGUE  MBFB
!***REFER TO  ODR
!***ROUTINES CALLED  DHSTEP
!***DATE WRITTEN   20040624   (YYYYMMDD)
!***REVISION DATE  20040624   (YYYYMMDD)
!***PURPOSE  ENSURE RANGE OF BOUNDS IS LARGE ENOUGH FOR DERIVATIVE CHECKING.
!***         MOVE BETA AWAY FROM BOUNDS SO THAT DERIVATIVES CAN BE CALCULATED.
!***END PROLOGUE  MBFB
!
!...USED MODULES
        use odrpack_kinds,only: wp
!
!...SCALAR ARGUMENTS
        integer                                                               &
         neta, np
        real(kind = wp)                                                       &
         eta
!
!...ARRAY ARGUMENTS
        integer                                                               &
         interval( np)
        real(kind = wp)                                                       &
         beta( np), lower( np), ssf( np), stpb( np), upper( np)
!
!...LOCAL SCALARS
        integer                                                               &
         k
        real(kind = wp)                                                       &
         h, h0, h1, hc, hc0, hc1, hundred, one, stpr, stpl, ten, three, typj, &
         zero
!
!...EXTERNAL FUNCTIONS
        real(kind = wp)                                                       &
         dhstep
        external                                                              &
         dhstep
!
!...DATA STATEMENTS
        data                                                                  &
         zero, one, ten, hundred, three                                       &
         /0.0E0_wp,1.0E0_wp,10.0E0_wp,100.0E0_wp,3.0E0_wp/
!
!...VARIABLE DEFINITIONS (ALPHABETICALLY)
!       BETA:    BETA for the jacobian checker.  BETA will be moved far enough from
!       the bounds so that the derivative checker may proceed.
!       H:       Relative step size for forward differences.
!       H0:      Initial relative step size for forward differences.
!       H1:      Default relative step size for forward differences.
!       HC:      Relative step size for center differences.
!       HC0:     Initial relative step size for center differences.
!       HC1:     Default relative step size for center differences.
!       HUNDRED: 100.0E0_wp
!       INTERVAL: Specifies which difference methods and step sizes are supported by
!       the current intervale UPPER-LOWER.
!       K:       Index variable for BETA.
!       NETA:    Number of good digits in the function results.
!       ONE:     The value 1.0E0_wp.
!       SSF:     The scale used for the BETA'S.
!       STPB:    The relative step used for computing finite difference derivatives
!       with respect to BETA.
!       STPL:    Maximum step to the left of BETA (-) the derivative checker will
!       use.
!       STPR:    Maximum step to the right of BETA (+) the derivative checker will
!       use.
!       TEN:     10.0E0_wp
!       THREE:   3.0E0_wp
!       TYPJ:    The typical size of the J-th unkonwn BETA.
!       ZERO:    The value 0.0E0_wp.
!
        interval(:) = 111
        do k = 1, np
           h0 = dhstep(0, neta,1, k, stpb,1)
           hc0 = h0
           h1 = sqrt( eta)
           hc1 = eta**( one/ three)
           h = max( ten* h1,min( hundred* h0, one))
           hc = max( ten* hc1,min( hundred* hc0, one))
           if ( beta( k) .eq. zero) then
!-----------------------------^------------------------------------------------
!!! FPT - 3087 REAL or COMPLEX quantity tested for exact equality/inequality
!------------------------------------------------------------------------------
              if ( ssf(1) .lt. zero) then
                 typj = one/abs( ssf(1))
              else
                 typj = one/ ssf( k)
              endif
           else
              typj = abs( beta( k))
           endif
           stpr = ( h* typj*sign( one, beta( k))+ beta( k))- beta( k)
           stpl = ( hc* typj*sign( one, beta( k))+ beta( k))- beta( k)
!          Check outer interval.
           if ( lower( k)+2*abs( stpl) .gt. upper( k)) then
              if ( interval( k) .ge. 100) then
                 interval( k) = interval( k)-100
              endif
           elseif ( beta( k)+ stpl .gt. upper( k) .or. beta( k)- stpl .gt.    &
            upper( k))                                                        &
            then
              beta( k) = upper( k)-abs( stpl)
           elseif ( beta( k)+ stpl .lt. lower( k) .or. beta( k)- stpl .lt.    &
            lower( k))                                                        &
            then
              beta( k) = lower( k)+abs( stpl)
           endif
!          Check middle interval.
           if ( lower( k)+2*abs( stpr) .gt. upper( k)) then
              if (mod( interval( k),100) .ge. 10) then
                 interval( k) = interval( k)-10
              endif
           elseif ( beta( k)+ stpr .gt. upper( k) .or. beta( k)- stpr .gt.    &
            upper( k))                                                        &
            then
              beta( k) = upper( k)-abs( stpr)
           elseif ( beta( k)+ stpr .lt. lower( k) .or. beta( k)- stpr .lt.    &
            lower( k))                                                        &
            then
              beta( k) = lower( k)+abs( stpr)
           endif
!          Check inner interval
           if ( lower( k)+abs( stpr) .gt. upper( k)) then
              interval( k) = 0
           elseif ( beta( k)+ stpr .gt. upper( k)) then
              beta( k) = upper( k)- stpr
           elseif ( beta( k)+ stpr .lt. lower( k)) then
              beta( k) = lower( k)- stpr
           endif
        enddo
!
end subroutine
!DERSTEP
function derstep                                                              &
         ( itype, k, betak, ssf, stpb, neta)                                  &
         result( derstepr)
!***Begin Prologue  DERSTEP
!***Refer to  ODR
!***Routines Called  DHSTEP
!***Date Written   20040616   (YYYYMMDD)
!***Revision Date  20040616   (YYYYMMDD)
!***Purpose  Compute step size for center and forward difference calculations
!***End Prologue  DERSTEP
!
!...Used modules
        use odrpack_kinds,only: wp
!
!...Scalar arguments
        integer                                                               &
         itype, k, neta
        real(kind = wp)                                                       &
         betak
!
!...Array arguments
        real(kind = wp)                                                       &
         ssf( k), stpb( k)
!
!...Result
        real(kind = wp)                                                       &
         derstepr
!
!...Local scalars
        real(kind = wp)                                                       &
         one, typj, zero
!
!...External functions
        real(kind = wp)                                                       &
         dhstep
        external                                                              &
         dhstep
!
!...Data statements
        data                                                                  &
         zero, one                                                            &
         /0.0E0_wp,1.0E0_wp/
!
!...Variable definitions (alphabetically)
!       BETAK:   The K-th function parameter.
!       ITYPE:   0 - calc foward difference step, 1 - calc center difference step.
!       K:       Index into beta where BETAK resides.
!       NETA:    Number of good digits in the function results.
!       ONE:     The value 1.0E0_wp.
!       SSF:     The scale used for the BETA'S.
!       STPB:    The relative step used for computing finite difference derivatives
!       with respect to BETA.
!       TYPJ:    The typical size of the J-th unkonwn BETA.
!       ZERO:    The value 0.0E0_wp.
!
!
!***First executable statement  DERSTEP
!
!
        if ( betak .eq. zero) then
!-----------------------^------------------------------------------------------
!!! FPT - 3087 REAL or COMPLEX quantity tested for exact equality/inequality
!------------------------------------------------------------------------------
           if ( ssf(1) .lt. zero) then
              typj = one/abs( ssf(1))
           else
              typj = one/ ssf( k)
           endif
        else
           typj = abs( betak)
        endif
        derstepr = sign( one, betak)* typj* dhstep( itype, neta,1, k, stpb,1)
!
        return
end function
!DJACFD
subroutine djacfd                                                             &
         ( fcn,                                                               &
         n, m, np, nq,                                                        &
         beta, x, ldx, delta, xplusd, ifixb, ifixx, ldifx,                    &
         stpb, stpd, ldstpd,                                                  &
         ssf, tt, ldtt, neta, fn, stp, wrk1, wrk2, wrk3, wrk6,                &
         fjacb, isodr, fjacd, nfev, istop, info,                              &
         lower, upper)
!***Begin Prologue  DJACFD
!***Refer to  ODR
!***Routines Called  FCN,DHSTEP,DERSTEP
!***Date Written   860529   (YYMMDD)
!***Revision Date  920619   (YYMMDD)
!***Purpose  Compute forward difference approximations to the
!       Jacobian wrt the estimated BETAS and wrt the DELTAS
!***End Prologue  DJACFD
!
!...Used modules
        use odrpack_kinds, only: wp, zero, one
!
!...Scalar arguments
        integer                                                               &
         info, istop, ldifx, ldstpd, ldtt, ldx, m, n, neta, nfev, np, nq
        logical                                                               &
         isodr
!
!...Array arguments
        real(kind = wp)                                                       &
         beta( np), delta( n, m), fjacb( n, np, nq), fjacd( n, m, nq), fn( n, &
         nq), lower( np), ssf( np), stp( n), stpb( np), stpd( ldstpd, m), tt( &
         ldtt, m), upper( np), wrk1( n, m, nq), wrk2( n, nq), wrk3( np), wrk6 &
         ( n, np, nq), x( ldx, m), xplusd( n, m)
        integer                                                               &
         ifixb( np), ifixx( ldifx, m)

!...Subroutine arguments
        external                                                              &
         fcn

!...Local scalars
        real(kind = wp) :: betak, step, typj
        integer :: i, j, k, l
        logical :: doit, setzro

!...External functions
        real(kind = wp), external :: dhstep, derstep

!...Routine names used as subprogram arguments
!       FCN:     The user supplied subroutine for evaluating the model.
!
!...Variable Definitions (alphabetically)
!       BETA:    The function parameters.
!       BETAK:   The K-th function parameter.
!       DELTA:   The estimated errors in the explanatory variables.
!       DOIT:    The variable designating whether the derivative wrt a
!       given BETA or DELTA needs to be computed (DOIT=TRUE)
!       or not (DOIT=FALSE).
!       FJACB:   The Jacobian with respect to BETA.
!       FJACD:   The Jacobian with respect to DELTA.
!       FN:      The new predicted values from the function.
!       I:       An indexing variable.
!       IFIXB:   The values designating whether the elements of BETA are
!       fixed at their input values or not.
!       IFIXX:   The values designating whether the elements of X are
!       fixed at their input values or not.
!       ISODR:   The variable designating whether the solution is by ODR
!       (ISODR=TRUE) or by OLS (ISODR=FALSE).
!       ISTOP:   The variable designating whether there are problems
!       computing the function at the current BETA and DELTA.
!       J:       An indexing variable.
!       K:       An indexing variable.
!       L:       An indexing variable.
!       LDIFX:   The leading dimension of array IFIXX.
!       LDSTPD:  The leading dimension of array STPD.
!       LDTT:    The leading dimension of array TT.
!       LDX:     The leading dimension of array X.
!       M:       The number of columns of data in the explanatory variable.
!       N:       The number of observations.
!       NETA:    The number of good digits in the function results.
!       NFEV:    The number of function evaluations.
!       NP:      The number of function parameters.
!       SETZRO:  The variable designating whether the derivative wrt some
!       DELTA needs to be set to zero (SETZRO=TRUE) or not
!       (SETZRO=FALSE).
!       SSF:     The scale used for the BETA'S.
!       STP:     The step used for computing finite difference
!       derivatives with respect to DELTA.
!       STPB:    The relative step used for computing finite difference
!       derivatives with respect to BETA.
!       STPD:    The relative step used for computing finite difference
!       derivatives with respect to DELTA.
!       TT:      The scaling values used for DELTA.
!       TYPJ:    The typical size of the J-th unknown BETA or DELTA.
!       X:       The explanatory variable.
!       XPLUSD:  The values of X + DELTA.
!       WRK1:    A work array of (N by M by NQ) elements.
!       WRK2:    A work array of (N BY NQ) elements.
!       WRK3:    A work array of (NP) elements.
!       WRK6:    A work array of (N BY NP BY NQ) elements.
!
!
!***First executable statement  DJACFD
!
!
!  Compute the Jacobian wrt the estimated BETAS
!
        do 40 k = 1, np
           if ( ifixb(1) .ge. 0) then
              if ( ifixb( k) .eq. 0) then
                 doit = .false.
              else
                 doit = .true.
              endif
           else
              doit = .true.
           endif
           if (.not. doit) then
              do 10 l = 1, nq
                 fjacb(1:n, k, l) = zero
10            continue
           else
              betak = beta( k)
              step = derstep(0, k, betak, ssf, stpb, neta)
              wrk3( k) = betak+ step
              wrk3( k) = wrk3( k)- betak
              beta( k) = betak+ wrk3( k)
              if ( beta( k) .gt. upper( k)) then
                 step = - step
                 wrk3( k) = betak+ step
                 wrk3( k) = wrk3( k)- betak
                 beta( k) = betak+ wrk3( k)
              endif
              if ( beta( k) .lt. lower( k)) then
                 step = - step
                 wrk3( k) = betak+ step
                 wrk3( k) = wrk3( k)- betak
                 beta( k) = betak+ wrk3( k)
                 if ( beta( k) .gt. upper( k)) then
                    info = 60001
                    return
                 endif
              endif
              istop = 0
              call fcn( n, m, np, nq,                                         &
               n, m, np,                                                      &
               beta, xplusd,                                                  &
               ifixb, ifixx, ldifx,                                           &
               001, wrk2, wrk6, wrk1,                                         &
               istop)
              if ( istop .ne. 0) then
                 return
              else
                 nfev = nfev+1
              endif
              do 30 l = 1, nq
                 do 20 i = 1, n
                    fjacb( i, k, l) = ( wrk2( i, l)- fn( i, l))/ wrk3( k)
20               continue
30            continue
              beta( k) = betak
           endif
40      continue
!
!  Compute the Jacobian wrt the X'S
!
        if ( isodr) then
           do 220 j = 1, m
              if ( ifixx(1,1) .lt. 0) then
                 doit = .true.
                 setzro = .false.
              elseif ( ldifx .eq. 1) then
                 if ( ifixx(1, j) .eq. 0) then
                    doit = .false.
                 else
                    doit = .true.
                 endif
                 setzro = .false.
              else
                 doit = .false.
                 setzro = .false.
                 do 100 i = 1, n
                    if ( ifixx( i, j) .ne. 0) then
                       doit = .true.
                    else
                       setzro = .true.
                    endif
100              continue
              endif
              if (.not. doit) then
                 do 110 l = 1, nq
                    fjacd(1:n, j, l) = zero
110              continue
              else
                 do 120 i = 1, n
                    if ( xplusd( i, j) .eq. zero) then
!-------------------------------------------^----------------------------------
!!! FPT - 3087 REAL or COMPLEX quantity tested for exact equality/inequality
!------------------------------------------------------------------------------
                       if ( tt(1,1) .lt. zero) then
                          typj = one/abs( tt(1,1))
                       elseif ( ldtt .eq. 1) then
                          typj = one/ tt(1, j)
                       else
                          typj = one/ tt( i, j)
                       endif
                    else
                       typj = abs( xplusd( i, j))
                    endif
!
                    stp( i) = xplusd( i, j)                                   &
                     +sign( one, xplusd( i, j))                               &
                     * typj* dhstep(0, neta, i, j, stpd, ldstpd)
                    stp( i) = stp( i)- xplusd( i, j)
                    xplusd( i, j) = xplusd( i, j)+ stp( i)
120              continue
!
                 istop = 0
                 call fcn( n, m, np, nq,                                      &
                  n, m, np,                                                   &
                  beta, xplusd,                                               &
                  ifixb, ifixx, ldifx,                                        &
                  001, wrk2, wrk6, wrk1,                                      &
                  istop)
                 if ( istop .ne. 0) then
                    return
                 else
                    nfev = nfev+1
                    do 140 l = 1, nq
                       do 130 i = 1, n
                          fjacd( i, j, l) = wrk2( i, l)
130                    continue
140                 continue
!
                 endif
!
                 if ( setzro) then
                    do 180 i = 1, n
                       if ( ifixx( i, j) .eq. 0) then
                          do 160 l = 1, nq
                             fjacd( i, j, l) = zero
160                       continue
                       else
                          do 170 l = 1, nq
                             fjacd( i, j, l) = ( fjacd( i, j, l)- fn( i, l))/ &
                              stp( i)
170                       continue
                       endif
180                 continue
                 else
                    do 200 l = 1, nq
                       do 190 i = 1, n
                          fjacd( i, j, l) = ( fjacd( i, j, l)- fn( i, l))/    &
                           stp( i)
190                    continue
200                 continue
                 endif
                 do 210 i = 1, n
                    xplusd( i, j) = x( i, j)+ delta( i, j)
210              continue
              endif
220        continue
        endif
!
        return
end subroutine
!DJCK
subroutine djck                                                               &
         ( fcn,                                                               &
         n, m, np, nq,                                                        &
         beta, betaj, xplusd,                                                 &
         ifixb, ifixx, ldifx, stpb, stpd, ldstpd,                             &
         ssf, tt, ldtt,                                                       &
         eta, neta, ntol, nrow, isodr, epsmac,                                &
         pv0i, fjacb, fjacd,                                                  &
         msgb, msgd, diff, istop, nfev, njev,                                 &
         wrk1, wrk2, wrk6,                                                    &
         interval)
!***Begin Prologue  DJCK
!***Refer to  ODR
!***Routines Called  FCN,DHSTEP,DJCKM
!***Date Written   860529   (YYMMDD)
!***Revision Date  920619   (YYMMDD)
!***Purpose  Driver routine for the derivative checking process
!       (adapted from STARPAC subroutine DCKCNT)
!***End Prologue  DJCK
!
!...Used modules
        use odrpack_kinds,only: wp
!
!...Scalar arguments
        real(kind = wp)                                                       &
         epsmac, eta
        integer                                                               &
         istop, ldifx, ldstpd, ldtt,                                          &
         m, n, neta, nfev, njev, np, nq, nrow, ntol
        logical                                                               &
         isodr
!
!...Array arguments
        real(kind = wp)                                                       &
         beta( np), betaj( np), diff( nq, np+ m), fjacb( n, np, nq), fjacd( n &
         , m, nq), pv0i( n, nq), ssf( np), stpb( np), stpd( ldstpd, m), tt(   &
         ldtt, m), wrk1( n, m, nq), wrk2( n, nq), wrk6( n, np, nq), xplusd( n &
         , m)
        integer                                                               &
         ifixb( np), ifixx( ldifx, m), interval( np), msgb(1+ nq* np),        &
         msgd(1+ nq* m)
!
!...Subroutine arguments
        external                                                              &
         fcn
!
!...Local scalars
        real(kind = wp)                                                       &
         diffj, h0, hc0, one, p5, pv, tol, typj, zero
        integer                                                               &
         ideval, j, lq, msgb1, msgd1
        logical                                                               &
         isfixd, iswrtb
!
!...Local arrays
        real(kind = wp)                                                       &
         pv0( n, nq)
!
!...External subroutines
        external                                                              &
         djckm
!
!...External functions
        real(kind = wp)                                                       &
         dhstep
        external                                                              &
         dhstep
!
!...Data statements
        data                                                                  &
         zero, p5, one                                                        &
         /0.0E0_wp,0.5E0_wp,1.0E0_wp/
!
!...Routine names used as subprogram arguments
!       FCN:     The user supplied subroutine for evaluating the model.
!
!...Variable Definitions (alphabetically)
!       BETA:    The function parameters.
!       BETAJ:   The function parameters offset such that steps don't cross
!       bounds.
!       DIFF:    The relative differences between the user supplied and
!       finite difference derivatives for each derivative checked.
!       DIFFJ:   The relative differences between the user supplied and
!       finite difference derivatives for the derivative being
!       checked.
!       EPSMAC:  The value of machine precision.
!       ETA:     The relative noise in the function results.
!       FJACB:   The Jacobian with respect to BETA.
!       FJACD:   The Jacobian with respect to DELTA.
!       H0:      The initial relative step size for forward differences.
!       HC0:     The initial relative step size for central differences.
!       IDEVAL:  The variable designating what computations are to be
!       performed by user supplied subroutine FCN.
!       IFIXB:   The values designating whether the elements of BETA are
!       fixed at their input values or not.
!       IFIXX:   The values designating whether the elements of X are
!       fixed at their input values or not.
!       INTERVAL: Specifies which checks can be performed when checking derivatives
!       based on the interval of the bound constraints.
!       ISFIXD:  The variable designating whether the parameter is fixed
!       (ISFIXD=TRUE) or not (ISFIXD=FALSE).
!       ISTOP:   The variable designating whether there are problems
!       computing the function at the current BETA and DELTA.
!       ISODR:   The variable designating whether the solution is by ODR
!       (ISODR=.TRUE.) or by OLS (ISODR=.FALSE.).
!       ISWRTB:  The variable designating whether the derivatives wrt BETA
!       (ISWRTB=TRUE) or DELTA (ISWRTB=FALSE) are being checked.
!       J:       An index variable.
!       LDIFX:   The leading dimension of array IFIXX.
!       LDSTPD:  The leading dimension of array STPD.
!       LDTT:    The leading dimension of array TT.
!       LQ:      The response currently being examined.
!       M:       The number of columns of data in the explanatory variable.
!       MSGB:    The error checking results for the Jacobian wrt BETA.
!       MSGB1:   The error checking results for the Jacobian wrt BETA.
!       MSGD:    The error checking results for the Jacobian wrt DELTA.
!       MSGD1:   The error checking results for the Jacobian wrt DELTA.
!       N:       The number of observations.
!       NETA:    The number of reliable digits in the model results, either
!       set by the user or computed by DETAF.
!       NFEV:    The number of function evaluations.
!       NJEV:    The number of Jacobian evaluations.
!       NP:      The number of function parameters.
!       NQ:      The number of responses per observation.
!       NROW:    The row number of the explanatory variable array at which
!       the derivative is checked.
!       NTOL:    The number of digits of agreement required between the
!       numerical derivatives and the user supplied derivatives.
!       ONE:     The value 1.0E0_wp.
!       P5:      The value 0.5E0_wp.
!       PV:      The scalar in which the predicted value from the model for
!       row   NROW   is stored.
!       PV0:     The predicted values using the current parameter estimates
!       (possibly offset from the user supplied estimates to create
!       distance between parameters and the bounds on the parameters).
!       PV0I:    The predicted values using the user supplied parameter estimates.
!       SSF:     The scaling values used for BETA.
!       STPB:    The step size for finite difference derivatives wrt BETA.
!       STPD:    The step size for finite difference derivatives wrt DELTA.
!       TOL:     The agreement tolerance.
!       TT:      The scaling values used for DELTA.
!       TYPJ:    The typical size of the J-th unknown BETA or DELTA.
!       WRK1:    A work array of (N BY M BY NQ) elements.
!       WRK2:    A work array of (N BY NQ) elements.
!       WRK6:    A work array of (N BY NP BY NQ) elements.
!       XPLUSD:  The values of X + DELTA.
!       ZERO:    The value 0.0E0_wp.
!
!
!***First executable statement  DJCK
!
!
!  Set tolerance for checking derivatives
!
        tol = eta**(0.25E0_wp)
        ntol = max( one, p5-log10( tol))
!
!
!  Compute, if necessary, PV0
!
        pv0 = pv0i
        if (any( beta(:) .ne. betaj(:))) then
           istop = 0
           ideval = 001
           call fcn( n, m, np, nq,                                            &
            n, m, np,                                                         &
            betaj, xplusd,                                                    &
            ifixb, ifixx, ldifx,                                              &
            ideval, pv0, fjacb, fjacd,                                        &
            istop)
           if ( istop .ne. 0) then
              return
           else
              njev = njev+1
           endif
        endif
!
!
!  Compute user supplied derivative values
!
        istop = 0
        if ( isodr) then
           ideval = 110
        else
           ideval = 010
        endif
        call fcn( n, m, np, nq,                                               &
         n, m, np,                                                            &
         betaj, xplusd,                                                       &
         ifixb, ifixx, ldifx,                                                 &
         ideval, wrk2, fjacb, fjacd,                                          &
         istop)
        if ( istop .ne. 0) then
           return
        else
           njev = njev+1
        endif
!
!  Check derivatives wrt BETA for each response of observation NROW
!
        msgb1 = 0
        msgd1 = 0
!
        do 30 lq = 1, nq
!
!  Set predicted value of model at current parameter estimates
           pv = pv0( nrow, lq)
!
           iswrtb = .true.
           do 10 j = 1, np
!
              if ( ifixb(1) .lt. 0) then
                 isfixd = .false.
              elseif ( ifixb( j) .eq. 0) then
                 isfixd = .true.
              else
                 isfixd = .false.
              endif
!
              if ( isfixd) then
                 msgb(1+ lq+( j-1)* nq) = -1
              else
                 if ( beta( j) .eq. zero) then
!-----------------------------------^------------------------------------------
!!! FPT - 3087 REAL or COMPLEX quantity tested for exact equality/inequality
!------------------------------------------------------------------------------
                    if ( ssf(1) .lt. zero) then
                       typj = one/abs( ssf(1))
                    else
                       typj = one/ ssf( j)
                    endif
                 else
                    typj = abs( beta( j))
                 endif
!
                 h0 = dhstep(0, neta,1, j, stpb,1)
                 hc0 = h0
!
!  Check derivative wrt the J-th parameter at the NROW-th row
!
                 if ( interval( j) .ge. 1) then
                    call djckm( fcn,                                          &
                     n, m, np, nq,                                            &
                     betaj, xplusd,                                           &
                     ifixb, ifixx, ldifx,                                     &
                     eta, tol, nrow, epsmac, j, lq, typj, h0, hc0,            &
                     iswrtb, pv, fjacb( nrow, j, lq),                         &
                     diffj, msgb1, msgb(2), istop, nfev,                      &
                     wrk1, wrk2, wrk6, interval)
                    if ( istop .ne. 0) then
                       msgb(1) = -1
                       return
                    else
                       diff( lq, j) = diffj
                    endif
                 else
                    msgb(1+ j) = 9
                 endif
              endif
!
10         continue
!
!  Check derivatives wrt X for each response of observation NROW
!
           if ( isodr) then
              iswrtb = .false.
              do 20 j = 1, m
!
                 if ( ifixx(1,1) .lt. 0) then
                    isfixd = .false.
                 elseif ( ldifx .eq. 1) then
                    if ( ifixx(1, j) .eq. 0) then
                       isfixd = .true.
                    else
                       isfixd = .false.
                    endif
                 else
                    isfixd = .false.
                 endif
!
                 if ( isfixd) then
                    msgd(1+ lq+( j-1)* nq) = -1
                 else
!
                    if ( xplusd( nrow, j) .eq. zero) then
!----------------------------------------------^-------------------------------
!!! FPT - 3087 REAL or COMPLEX quantity tested for exact equality/inequality
!------------------------------------------------------------------------------
                       if ( tt(1,1) .lt. zero) then
                          typj = one/abs( tt(1,1))
                       elseif ( ldtt .eq. 1) then
                          typj = one/ tt(1, j)
                       else
                          typj = one/ tt( nrow, j)
                       endif
                    else
                       typj = abs( xplusd( nrow, j))
                    endif
!
                    h0 = dhstep(0, neta, nrow, j, stpd, ldstpd)
                    hc0 = dhstep(1, neta, nrow, j, stpd, ldstpd)
!
!  Check derivative wrt the J-th column of DELTA at row NROW
!
                    call djckm( fcn,                                          &
                     n, m, np, nq,                                            &
                     betaj, xplusd,                                           &
                     ifixb, ifixx, ldifx,                                     &
                     eta, tol, nrow, epsmac, j, lq, typj, h0, hc0,            &
                     iswrtb, pv, fjacd( nrow, j, lq),                         &
                     diffj, msgd1, msgd(2), istop, nfev,                      &
                     wrk1, wrk2, wrk6, interval)
                    if ( istop .ne. 0) then
                       msgd(1) = -1
                       return
                    else
                       diff( lq, np+ j) = diffj
                    endif
                 endif
!
20            continue
           endif
30      continue
        msgb(1) = msgb1
        msgd(1) = msgd1
!
        return
end subroutine
!DJCKC
subroutine djckc                                                              &
         ( fcn,                                                               &
         n, m, np, nq,                                                        &
         beta, xplusd, ifixb, ifixx, ldifx,                                   &
         eta, tol, nrow, epsmac, j, lq, hc, iswrtb,                           &
         fd, typj, pvpstp, stp0,                                              &
         pv, d,                                                               &
         diffj, msg, istop, nfev,                                             &
         wrk1, wrk2, wrk6)
!***Begin Prologue  DJCKC
!***Refer to  ODR
!***Routines Called  DJCKF,DPVB,DPVD
!***Date Written   860529   (YYMMDD)
!***Revision Date  920619   (YYMMDD)
!***Purpose  Check whether high curvature could be the cause of the
!       disagreement between the numerical and analytic derviatives
!       (adapted from STARPAC subroutine DCKCRV)
!***End prologue  DJCKC
!
!...Used modules
        use odrpack_kinds,only: wp
!
!...Scalar arguments
        real(kind = wp)                                                       &
         d, diffj, epsmac, eta, fd, hc, pv, pvpstp, stp0, tol, typj
        integer                                                               &
         istop, j, ldifx, lq, m, n, nfev, np, nq, nrow
        logical                                                               &
         iswrtb
!
!...Array arguments
        real(kind = wp)                                                       &
         beta( np), wrk1( n, m, nq), wrk2( n, nq), wrk6( n, np, nq), xplusd(  &
         n, m)
        integer                                                               &
         ifixb( np), ifixx( ldifx, m), msg( nq, j)
!
!...Subroutine arguments
        external                                                              &
         fcn
!
!...Local scalars
        real(kind = wp)                                                       &
         curve, one, pvmcrv, pvpcrv, p01, stp, stpcrv, ten, two
!
!...External subroutines
        external                                                              &
         djckf, dpvb, dpvd
!
!...Data statements
        data                                                                  &
         p01, one, two, ten                                                   &
         /0.01E0_wp,1.0E0_wp,2.0E0_wp,10.0E0_wp/
!
!...Routine names used as subprogram arguments
!       FCN:     The user supplied subroutine for evaluating the model.
!
!...Variable Definitions (alphabetically)
!       BETA:    The function parameters.
!       CURVE:   A measure of the curvature in the model.
!       D:       The derivative with respect to the Jth unknown parameter.
!       DIFFJ:   The relative differences between the user supplied and
!       finite difference derivatives for the derivative being
!       checked.
!       EPSMAC:  The value of machine precision.
!       ETA:     The relative noise in the model
!       FD:      The forward difference derivative wrt the Jth parameter.
!       HC:      The relative step size for central finite differences.
!       IFIXB:   The values designating whether the elements of BETA are
!       fixed at their input values or not.
!       IFIXX:   The values designating whether the elements of X are
!       fixed at their input values or not.
!       ISTOP:   The variable designating whether there are problems
!       computing the function at the current BETA and DELTA.
!       ISWRTB:  The variable designating whether the derivatives wrt BETA
!       (ISWRTB=TRUE) or DELTA(ISWRTB=FALSE) are being checked.
!       J:       The index of the partial derivative being examined.
!       LDIFX:   The leading dimension of array IFIXX.
!       LQ:      The response currently being examined.
!       M:       The number of columns of data in the explanatory variable.
!       MSG:     The error checking results.
!       N:       The number of observations.
!       NFEV:    The number of function evaluations.
!       NP:      The number of function parameters.
!       NQ:      The number of responses per observation.
!       NROW:    The row number of the explanatory variable array at which
!       the derivative is to be checked.
!       ONE:     The value 1.0E0_wp.
!       PV:      The predicted value of the model for row   NROW   .
!       PVMCRV:  The predicted value for row    NROW   of the model
!       based on the current parameter estimates for all but the
!       Jth parameter value, which is BETA(J)-STPCRV.
!       PVPCRV:  The predicted value for row    NROW   of the model
!       based on the current parameter estimates for all but the
!       Jth parameter value, which is BETA(J)+STPCRV.
!       PVPSTP:  The predicted value for row    NROW   of the model
!       based on the current parameter estimates for all but the
!       Jth parameter value, which is BETA(J) + STP0.
!       P01:     The value 0.01E0_wp.
!       STP0:    The initial step size for the finite difference derivative.
!       STP:     A step size for the finite difference derivative.
!       STPCRV:  The step size selected to check for curvature in the model.
!       TEN:     The value 10.0E0_wp.
!       TOL:     The agreement tolerance.
!       TWO:     The value 2.0E0_wp.
!       TYPJ:    The typical size of the J-th unknown BETA or DELTA.
!       WRK1:    A work array of (N BY M BY NQ) elements.
!       WRK2:    A work array of (N BY NQ) elements.
!       WRK6:    A work array of (N BY NP BY NQ) elements.
!       XPLUSD:  The values of X + DELTA.
!
!
!***First executable statement  DJCKC
!
!
        if ( iswrtb) then
!
!  Perform central difference computations for derivatives wrt BETA
!
           stpcrv = ( hc* typj*sign( one, beta( j))+ beta( j))- beta( j)
           call dpvb( fcn,                                                    &
            n, m, np, nq,                                                     &
            beta, xplusd, ifixb, ifixx, ldifx,                                &
            nrow, j, lq, stpcrv,                                              &
            istop, nfev, pvpcrv,                                              &
            wrk1, wrk2, wrk6)
           if ( istop .ne. 0) then
              return
           endif
           call dpvb( fcn,                                                    &
            n, m, np, nq,                                                     &
            beta, xplusd, ifixb, ifixx, ldifx,                                &
            nrow, j, lq,- stpcrv,                                             &
            istop, nfev, pvmcrv,                                              &
            wrk1, wrk2, wrk6)
           if ( istop .ne. 0) then
              return
           endif
        else
!
!  Perform central difference computations for derivatives wrt DELTA
!
           stpcrv = ( hc* typj*sign( one, xplusd( nrow, j))+ xplusd( nrow, j) &
            )- xplusd( nrow, j)
           call dpvd( fcn,                                                    &
            n, m, np, nq,                                                     &
            beta, xplusd, ifixb, ifixx, ldifx,                                &
            nrow, j, lq, stpcrv,                                              &
            istop, nfev, pvpcrv,                                              &
            wrk1, wrk2, wrk6)
           if ( istop .ne. 0) then
              return
           endif
           call dpvd( fcn,                                                    &
            n, m, np, nq,                                                     &
            beta, xplusd, ifixb, ifixx, ldifx,                                &
            nrow, j, lq,- stpcrv,                                             &
            istop, nfev, pvmcrv,                                              &
            wrk1, wrk2, wrk6)
           if ( istop .ne. 0) then
              return
           endif
        endif
!
!  Estimate curvature by second derivative of model
!
        curve = abs(( pvpcrv- pv)+( pvmcrv- pv))/( stpcrv* stpcrv)
        curve = curve+                                                        &
         eta*(abs( pvpcrv)+abs( pvmcrv)+ two*abs( pv))/( stpcrv**2)
!
!
!  Check if finite precision arithmetic could be the culprit.
        call djckf( fcn,                                                      &
         n, m, np, nq,                                                        &
         beta, xplusd, ifixb, ifixx, ldifx,                                   &
         eta, tol, nrow, j, lq, iswrtb,                                       &
         fd, typj, pvpstp, stp0, curve, pv, d,                                &
         diffj, msg, istop, nfev,                                             &
         wrk1, wrk2, wrk6)
        if ( istop .ne. 0) then
           return
        endif
        if ( msg( lq, j) .eq. 0) then
           return
        endif
!
!  Check if high curvature could be the problem.
!
        stp = two*max( tol*abs( d)/ curve, epsmac)
        if ( stp .lt. abs( ten* stp0)) then
           stp = min( stp, p01*abs( stp0))
        endif
!
!
        if ( iswrtb) then
!
!  Perform computations for derivatives wrt BETA
           stp = ( stp*sign( one, beta( j))+ beta( j))- beta( j)
           call dpvb( fcn,                                                    &
            n, m, np, nq,                                                     &
            beta, xplusd, ifixb, ifixx, ldifx,                                &
            nrow, j, lq, stp,                                                 &
            istop, nfev, pvpstp,                                              &
            wrk1, wrk2, wrk6)
           if ( istop .ne. 0) then
              return
           endif
        else
!
!  Perform computations for derivatives wrt DELTA
           stp = ( stp*sign( one, xplusd( nrow, j))+ xplusd( nrow, j))-       &
            xplusd( nrow, j)
           call dpvd( fcn,                                                    &
            n, m, np, nq,                                                     &
            beta, xplusd, ifixb, ifixx, ldifx,                                &
            nrow, j, lq, stp,                                                 &
            istop, nfev, pvpstp,                                              &
            wrk1, wrk2, wrk6)
           if ( istop .ne. 0) then
              return
           endif
        endif
!
!  Compute the new numerical derivative
!
        fd = ( pvpstp- pv)/ stp
        diffj = min( diffj,abs( fd- d)/abs( d))
!
!  Check whether the new numerical derivative is ok
        if (abs( fd- d) .le. tol*abs( d)) then
           msg( lq, j) = 0
!
!  Check if finite precision may be the culprit (fudge factor = 2)
        elseif (abs( stp*( fd- d)) .lt. two* eta*(abs( pv)+abs( pvpstp))+     &
         curve*( epsmac* typj)**2) then
           msg( lq, j) = 5
        endif
!
        return
end subroutine
!DJCKF
subroutine djckf                                                              &
         ( fcn,                                                               &
         n, m, np, nq,                                                        &
         beta, xplusd, ifixb, ifixx, ldifx,                                   &
         eta, tol, nrow, j, lq, iswrtb,                                       &
         fd, typj, pvpstp, stp0, curve, pv, d,                                &
         diffj, msg, istop, nfev,                                             &
         wrk1, wrk2, wrk6)
!***Begin Prologue  DJCKF
!***Refer to  ODR
!***Routines Called  DPVB,DPVD
!***Date Written   860529   (YYMMDD)
!***Revision Date  920619   (YYMMDD)
!***Purpose  Check whether finite precision arithmetic could be the
!       cause of the disagreement between the derivatives
!       (adapted from STARPAC subroutine DCKFPA)
!***End Prologue  DJCKF
!
!...Used modules
        use odrpack_kinds,only: wp
!
!...Scalar arguments
        real(kind = wp)                                                       &
         curve, d, diffj, eta, fd, pv, pvpstp, stp0, tol, typj
        integer                                                               &
         istop, j, ldifx, lq, m, n, nfev, np, nq, nrow
        logical                                                               &
         iswrtb
!
!...Array arguments
        real(kind = wp)                                                       &
         beta( np), wrk1( n, m, nq), wrk2( n, nq), wrk6( n, np, nq), xplusd(  &
         n, m)
        integer                                                               &
         ifixb( np), ifixx( ldifx, m), msg( nq, j)
!
!...Subroutine arguments
        external                                                              &
         fcn
!
!...Local scalars
        real(kind = wp)                                                       &
         hundrd, one, p1, stp, two
        logical                                                               &
         large
!
!...External subroutines
        external                                                              &
         dpvb, dpvd
!
!...Data statements
        data                                                                  &
         p1, one, two, hundrd                                                 &
         /0.1E0_wp,1.0E0_wp,2.0E0_wp,100.0E0_wp/
!
!...Routine names used as subprogram arguments
!       FCN:     The user supplied subroutine for evaluating the model.
!
!...Variable Definitions (alphabetically)
!       BETA:    The function parameters.
!       CURVE:   A measure of the curvature in the model.
!       D:       The derivative with respect to the Jth unknown parameter.
!       DIFFJ:   The relative differences between the user supplied and
!       finite difference derivatives for the derivative being
!       checked.
!       ETA:     The relative noise in the model
!       FD:      The forward difference derivative wrt the Jth parameter.
!       HUNDRD:  The value 100.0E0_wp.
!       IFIXB:   The values designating whether the elements of BETA are
!       fixed at their input values or not.
!       IFIXX:   The values designating whether the elements of X are
!       fixed at their input values or not.
!       ISTOP:   The variable designating whether there are problems
!       computing the function at the current BETA and DELTA.
!       ISWRTB:  The variable designating whether the derivatives wrt BETA
!       (ISWRTB=TRUE) or DELTA(ISWRTB=FALSE) are being checked.
!       J:       The index of the partial derivative being examined.
!       LARGE:   The value designating whether the recommended increase in
!       the step size would be greater than TYPJ.
!       LDIFX:   The leading dimension of array IFIXX.
!       LQ:      The response currently being examined.
!       M:       The number of columns of data in the explanatory variable.
!       MSG:     The error checking results.
!       N:       The number of observations.
!       NFEV:    The number of function evaluations.
!       NP:      The number of function parameters.
!       NQ:      The number of responses per observation.
!       NROW:    The row number of the explanatory variable array at which
!       the derivative is to be checked.
!       ONE:     The value 1.0E0_wp.
!       PV:      The predicted value for row   NROW   .
!       PVPSTP:  The predicted value for row    NROW   of the model
!       based on the current parameter estimates for all but the
!       Jth parameter value, which is BETA(J) + STP0.
!       P1:      The value 0.1E0_wp.
!       STP0:    The step size for the finite difference derivative.
!       TOL:     The agreement tolerance.
!       TWO:     The value 2.0E0_wp.
!       TYPJ:    The typical size of the J-th unknown BETA or DELTA.
!       WRK1:    A work array of (N BY M BY NQ) elements.
!       WRK2:    A work array of (N BY NQ) elements.
!       WRK6:    A work array of (N BY NP BY NQ) elements.
!       XPLUSD:  The values of X + DELTA.
!
!
!***First executable statement  DJCKF
!
!
!  Finite precision arithmetic could be the problem.
!  Try a larger step size based on estimate of condition error
!
        stp = eta*(abs( pv)+abs( pvpstp))/( tol*abs( d))
        if ( stp .gt. abs( p1* stp0)) then
           stp = max( stp, hundrd*abs( stp0))
        endif
        if ( stp .gt. typj) then
           stp = typj
           large = .true.
        else
           large = .false.
        endif
!
        if ( iswrtb) then
!
!  Perform computations for derivatives wrt BETA
           stp = ( stp*sign( one, beta( j))+ beta( j))- beta( j)
           call dpvb( fcn,                                                    &
            n, m, np, nq,                                                     &
            beta, xplusd, ifixb, ifixx, ldifx,                                &
            nrow, j, lq, stp,                                                 &
            istop, nfev, pvpstp,                                              &
            wrk1, wrk2, wrk6)
        else
!
!  Perform computations for derivatives wrt DELTA
           stp = ( stp*sign( one, xplusd( nrow, j))+ xplusd( nrow, j))-       &
            xplusd( nrow, j)
           call dpvd( fcn,                                                    &
            n, m, np, nq,                                                     &
            beta, xplusd, ifixb, ifixx, ldifx,                                &
            nrow, j, lq, stp,                                                 &
            istop, nfev, pvpstp,                                              &
            wrk1, wrk2, wrk6)
        endif
        if ( istop .ne. 0) then
           return
        endif
!
        fd = ( pvpstp- pv)/ stp
        diffj = min( diffj,abs( fd- d)/abs( d))
!
!  Check for agreement
!
        if ((abs( fd- d)) .le. tol*abs( d)) then
!  Forward difference quotient and analytic derivatives agree.
           msg( lq, j) = 0
!
        elseif ((abs( fd- d) .le. abs( two* curve* stp)) .or. large) then
!  Curvature may be the culprit (fudge factor = 2)
           if ( large) then
              msg( lq, j) = 4
           else
              msg( lq, j) = 5
           endif
        endif
!
        return
end subroutine
!DJCKM
subroutine djckm                                                              &
         ( fcn,                                                               &
         n, m, np, nq,                                                        &
         beta, xplusd, ifixb, ifixx, ldifx,                                   &
         eta, tol, nrow, epsmac, j, lq, typj, h0, hc0,                        &
         iswrtb, pv, d,                                                       &
         diffj, msg1, msg, istop, nfev,                                       &
         wrk1, wrk2, wrk6, interval)
!***Begin Prologue  DJCKM
!***Refer to  ODR
!***Routines Called  DJCKC,DJCKZ,DPVB,DPVD
!***Date Written   860529   (YYMMDD)
!***Revision Date  920619   (YYMMDD)
!***Purpose  Check user supplied analytic derivatives against numerical
!       derivatives
!       (adapted from STARPAC subroutine DCKMN)
!***End prologue  DJCKM
!
!...Used modules
        use odrpack_kinds,only: wp
!
!...Scalar arguments
        real(kind = wp)                                                       &
         d, diffj, epsmac, eta, h0, hc0, pv, tol, typj
        integer                                                               &
         istop, j, ldifx, lq, m, msg1, n, nfev, np, nq, nrow
        logical                                                               &
         iswrtb
!
!...Array arguments
        real(kind = wp)                                                       &
         beta( np), wrk1( n, m, nq), wrk2( n, nq), wrk6( n, np, nq), xplusd(  &
         n, m)
        integer                                                               &
         ifixb( np), ifixx( ldifx, m), interval( np), msg( nq, j)
!
!...Subroutine arguments
        external                                                              &
         fcn
!
!...Local scalars
        real(kind = wp)                                                       &
         big, fd, h, hc, h1, hc1, hundrd, one, pvpstp, p01, p1, stp0,         &
         ten, three, tol2, two, zero
        integer                                                               &
         i
!
!...External subroutines
        external                                                              &
         djckc, djckz, dpvb, dpvd
!
!...Data statements
        data                                                                  &
         zero, p01, p1, one, two, three, ten, hundrd                          &
         /0.0E0_wp,0.01E0_wp,0.1E0_wp,1.0E0_wp,2.0E0_wp,3.0E0_wp,             &
         1.0E1_wp,1.0E2_wp/
        data                                                                  &
         big, tol2                                                            &
         /1.0E19_wp,5.0E-2_wp/
!
!...Routine names used as subprogram arguments
!       FCN:     The user supplied subroutine for evaluating the model.
!
!...Variable Definitions (alphabetically)
!       BETA:    The function parameters.
!       BIG:     A big value, used to initialize DIFFJ.
!       D:       The derivative with respect to the Jth unknown parameter.
!       DIFFJ:   The relative differences between the user supplied and
!       finite difference derivatives for the derivative being
!       checked.
!       EPSMAC:  The value of machine precision.
!       ETA:     The relative noise in the function results.
!       FD:      The forward difference derivative wrt the Jth parameter.
!       H:       The relative step size for forward differences.
!       H0:      The initial relative step size for forward differences.
!       H1:      The default relative step size for forward differences.
!       HC:      The relative step size for central differences.
!       HC0:     The initial relative step size for central differences.
!       HC1:     The default relative step size for central differences.
!       HUNDRD:  The value 100.0E0_wp.
!       IFIXB:   The values designating whether the elements of BETA are
!       fixed at their input values or not.
!       IFIXX:   The values designating whether the elements of X are
!       fixed at their input values or not.
!       INTERVAL: Specifies which checks can be performed when checking derivatives
!       based on the interval of the bound constraints.
!       ISTOP:   The variable designating whether there are problems
!       computing the function at the current BETA and DELTA.
!       ISWRTB:  The variable designating whether the derivatives wrt BETA
!       (ISWRTB=TRUE) or DELTAS (ISWRTB=FALSE) are being checked.
!       J:       The index of the partial derivative being examined.
!       LDIFX:   The leading dimension of array IFIXX.
!       LQ:      The response currently being examined.
!       M:       The number of columns of data in the explanatory variable.
!       MSG:     The error checking results.
!       MSG1:    The error checking results summary.
!       N:       The number of observations.
!       NFEV:    The number of function evaluations.
!       NP:      The number of function parameters.
!       NQ:      The number of responses per observation.
!       NROW:    The row number of the explanatory variable array at which
!       the derivative is to be checked.
!       ONE:     The value 1.0E0_wp.
!       PV:      The predicted value from the model for row   NROW   .
!       PVPSTP:  The predicted value for row    NROW   of the model
!       Using the current parameter estimates for all but the Jth
!       parameter value, which is BETA(J) + STP0.
!       P01:     The value 0.01E0_wp.
!       P1:      The value 0.1E0_wp.
!       STP0:    The initial step size for the finite difference derivative.
!       TEN:     The value 10.0E0_wp.
!       THREE:   The value 3.0E0_wp.
!       TWO:     The value 2.0E0_wp.
!       TOL:     The agreement tolerance.
!       TOL2:    A minimum agreement tolerance.
!       TYPJ:    The typical size of the J-th unknown BETA or DELTA.
!       WRK1:    A work array of (N BY M BY NQ) elements.
!       WRK2:    A work array of (N BY NQ) elements.
!       WRK6:    A work array of (N BY NP BY NQ) elements.
!       XPLUSD:  The values of X + DELTA.
!       ZERO:    The value 0.0E0_wp.
!
!
!***First executable statement  DJCKM
!
!
!  Calculate the Jth partial derivative using forward difference
!  quotients and decide if it agrees with user supplied values
!
        h1 = sqrt( eta)
        hc1 = eta**( one/ three)
!
        msg( lq, j) = 7
        diffj = big
!
        do 10 i = 1,3
!
           if ( i .eq. 1) then
!  Try initial relative step size
              h = h0
              hc = hc0
!
           elseif ( i .eq. 2) then
!  Try larger relative step size
              h = max( ten* h1,min( hundrd* h0, one))
              hc = max( ten* hc1,min( hundrd* hc0, one))
!
           elseif ( i .eq. 3) then
!  Try smaller relative step size
              h = min( p1* h1,max( p01* h, two* epsmac))
              hc = min( p1* hc1,max( p01* hc, two* epsmac))
           endif
!
           if ( iswrtb) then
!
!  Perform computations for derivatives wrt BETA
!
              stp0 = ( h* typj*sign( one, beta( j))+ beta( j))- beta( j)
              call dpvb( fcn,                                                 &
               n, m, np, nq,                                                  &
               beta, xplusd, ifixb, ifixx, ldifx,                             &
               nrow, j, lq, stp0,                                             &
               istop, nfev, pvpstp,                                           &
               wrk1, wrk2, wrk6)
           else
!
!  Perform computations for derivatives wrt DELTA
!
              stp0 = ( h* typj*sign( one, xplusd( nrow, j))+ xplusd( nrow, j) &
               )- xplusd( nrow, j)
              call dpvd( fcn,                                                 &
               n, m, np, nq,                                                  &
               beta, xplusd, ifixb, ifixx, ldifx,                             &
               nrow, j, lq, stp0,                                             &
               istop, nfev, pvpstp,                                           &
               wrk1, wrk2, wrk6)
           endif
           if ( istop .ne. 0) then
              return
           endif
!
           fd = ( pvpstp- pv)/ stp0
!
!  Check for agreement
!
           if (abs( fd- d) .le. tol*abs( d)) then
!  Numerical and analytic derivatives agree
!
!  Set relative difference for derivative checking report
              if (( d .eq. zero) .or.                                         &
!--------------------------^---------------------------------------------------
!!! FPT - 3087 REAL or COMPLEX quantity tested for exact equality/inequality
!------------------------------------------------------------------------------
                  ( fd .eq. zero)) then
!---------------------------^--------------------------------------------------
!!! FPT - 3087 REAL or COMPLEX quantity tested for exact equality/inequality
!------------------------------------------------------------------------------
                 diffj = abs( fd- d)
              else
                 diffj = abs( fd- d)/abs( d)
              endif
!
!  Set MSG flag.
              if ( d .eq. zero) then
!-------------------------^----------------------------------------------------
!!! FPT - 3087 REAL or COMPLEX quantity tested for exact equality/inequality
!------------------------------------------------------------------------------
!
!  JTH analytic and numerical derivatives are both zero.
                 msg( lq, j) = 1
!
              else
!  JTH analytic and numerical derivatives are both nonzero.
                 msg( lq, j) = 0
              endif
!
           else
!
!  Numerical and analytic derivatives disagree.  Check why
              if (( d .eq. zero) .or.                                         &
!--------------------------^---------------------------------------------------
!!! FPT - 3087 REAL or COMPLEX quantity tested for exact equality/inequality
!------------------------------------------------------------------------------
                  ( fd .eq. zero)) then
!---------------------------^--------------------------------------------------
!!! FPT - 3087 REAL or COMPLEX quantity tested for exact equality/inequality
!------------------------------------------------------------------------------
                 if ( interval( j) .ge. 10 .or. .not. iswrtb) then
                    call djckz( fcn,                                          &
                     n, m, np, nq,                                            &
                     beta, xplusd, ifixb, ifixx, ldifx,                       &
                     nrow, epsmac, j, lq, iswrtb,                             &
                     tol, d, fd, typj, pvpstp, stp0, pv,                      &
                     diffj, msg, istop, nfev,                                 &
                     wrk1, wrk2, wrk6)
                 else
                    msg( lq, j) = 8
                 endif
              else
                 if ( interval( j) .ge. 100 .or. .not. iswrtb) then
                    call djckc( fcn,                                          &
                     n, m, np, nq,                                            &
                     beta, xplusd, ifixb, ifixx, ldifx,                       &
                     eta, tol, nrow, epsmac, j, lq, hc, iswrtb,               &
                     fd, typj, pvpstp, stp0, pv, d,                           &
                     diffj, msg, istop, nfev,                                 &
                     wrk1, wrk2, wrk6)
                 else
                    msg( lq, j) = 8
                 endif
              endif
              if ( msg( lq, j) .le. 2) then
                 goto 20
              endif
           endif
10      continue
!
!  Set summary flag to indicate questionable results
20      continue
        if (( msg( lq, j) .ge. 7) .and.                                       &
            ( diffj .le. tol2)) msg( lq, j) = 6
        if (( msg( lq, j) .ge. 1) .and.                                       &
            ( msg( lq, j) .le. 6)) then
           msg1 = max( msg1,1)
        elseif ( msg( lq, j) .ge. 7) then
           msg1 = 2
        endif
!
        return
end subroutine
!DJCKZ
subroutine djckz                                                              &
         ( fcn,                                                               &
         n, m, np, nq,                                                        &
         beta, xplusd, ifixb, ifixx, ldifx,                                   &
         nrow, epsmac, j, lq, iswrtb,                                         &
         tol, d, fd, typj, pvpstp, stp0, pv,                                  &
         diffj, msg, istop, nfev,                                             &
         wrk1, wrk2, wrk6)
!***Begin Prologue  DJCKZ
!***Refer to  ODR
!***Routines Called  DPVB,DPVD
!***Date Written   860529   (YYMMDD)
!***Revision Date  920619   (YYMMDD)
!***Purpose  Recheck the derivatives in the case where the finite
!       difference derivative disagrees with the analytic
!       derivative and the analytic derivative is zero
!       (adapted from STARPAC subroutine DCKZRO)
!***End Prologue  DJCKZ
!
!...Used modules
        use odrpack_kinds,only: wp
!
!...Scalar arguments
        real(kind = wp)                                                       &
         d, diffj, epsmac, fd, pv, pvpstp, stp0, tol, typj
        integer                                                               &
         istop, j, ldifx, lq, m, n, nfev, np, nq, nrow
        logical                                                               &
         iswrtb
!
!...Array arguments
        real(kind = wp)                                                       &
         beta( np), wrk1( n, m, nq), wrk2( n, nq), wrk6( n, np, nq), xplusd(  &
         n, m)
        integer                                                               &
         ifixb( np), ifixx( ldifx, m), msg( nq, j)
!
!...Subroutine arguments
        external                                                              &
         fcn
!
!...Local scalars
        real(kind = wp)                                                       &
         cd, one, pvmstp, three, two, zero
!
!...External subroutines
        external                                                              &
         dpvb, dpvd
!
!...Data statements
        data                                                                  &
         zero, one, two, three                                                &
         /0.0E0_wp,1.0E0_wp,2.0E0_wp,3.0E0_wp/
!
!...Routine names used as subprogram arguments
!       FCN:     THE USER SUPPLIED SUBROUTINE FOR EVALUATING THE MODEL.
!
!...Variable Definitions (alphabetically)
!       BETA:    The function parameters.
!       CD:      The central difference derivative wrt the Jth parameter.
!       D:       The derivative with respect to the Jth unknown parameter.
!       DIFFJ:   The relative differences between the user supplied and
!       finite difference derivatives for the derivative being
!       checked.
!       EPSMAC:  The value of machine precision.
!       FD:      The forward difference derivative wrt the Jth parameter.
!       IFIXB:   The values designating whether the elements of BETA are
!       Fixed at their input values or not.
!       IFIXX:   The values designating whether the elements of X are
!       fixed at their input values or not.
!       ISTOP:   The variable designating whether there are problems
!       computing the function at the current BETA and DELTA.
!       ISWRTB:  The variable designating whether the derivatives wrt BETA
!       (ISWRTB=TRUE) or X (ISWRTB=FALSE) are being checked.
!       J:       The index of the partial derivative being examined.
!       LDIFX:   The leading dimension of array IFIXX.
!       LQ:      The response currently being examined.
!       M:       The number of columns of data in the explanatory variable.
!       MSG:     The error checking results.
!       N:       The number of observations.
!       NFEV:    The number of function evaluations.
!       NP:      The number of function parameters.
!       NQ:      The number of responses per observation.
!       NROW:    The row number of the explanatory variable array at which
!       The derivative is to be checked.
!       ONE:     The value 1.0E0_wp.
!       PV:      The predicted value from the model for row   NROW   .
!       PVMSTP:  The predicted value for row    NROW   of the model
!       using the current parameter estimates for all but the
!       Jth parameter value, which is BETA(J) - STP0.
!       PVPSTP:  The predicted value for row    NROW   of the model
!       using the current parameter estimates for all but the
!       JTH parameter value, which is BETA(J) + STP0.
!       STP0:    The initial step size for the finite difference derivative.
!       THREE:   The value 3.0E0_wp.
!       TWO:     The value 2.0E0_wp.
!       TOL:     The agreement tolerance.
!       TYPJ:    The typical size of the J-th unknown BETA or DELTA.
!       WRK1:    A work array of (N BY M BY NQ) elements.
!       WRK2:    A work array of (N BY NQ) elements.
!       WRK6:    A work array of (N BY NP BY NQ) elements.
!       XPLUSD:  The values of X + DELTA.
!       ZERO:    The value 0.0E0_wp.
!
!
!***First executable statement  DJCKZ
!
!
!  Recalculate numerical derivative using central difference and step
!  size of 2*STP0
!
        if ( iswrtb) then
!
!  Perform computations for derivatives wrt BETA
!
           call dpvb( fcn,                                                    &
            n, m, np, nq,                                                     &
            beta, xplusd, ifixb, ifixx, ldifx,                                &
            nrow, j, lq,- stp0,                                               &
            istop, nfev, pvmstp,                                              &
            wrk1, wrk2, wrk6)
        else
!
!  Perform computations for derivatives wrt DELTA
!
           call dpvd( fcn,                                                    &
            n, m, np, nq,                                                     &
            beta, xplusd, ifixb, ifixx, ldifx,                                &
            nrow, j, lq,- stp0,                                               &
            istop, nfev, pvmstp,                                              &
            wrk1, wrk2, wrk6)
        endif
        if ( istop .ne. 0) then
           return
        endif
!
        cd = ( pvpstp- pvmstp)/( two* stp0)
        diffj = min(abs( cd- d),abs( fd- d))
!
!  Check for agreement
!
        if ( diffj .le. tol*abs( d)) then
!
!  Finite difference and analytic derivatives now agree.
           if ( d .eq. zero) then
!----------------------^-------------------------------------------------------
!!! FPT - 3087 REAL or COMPLEX quantity tested for exact equality/inequality
!------------------------------------------------------------------------------
              msg( lq, j) = 1
           else
              msg( lq, j) = 0
           endif
!
        elseif ( diffj* typj .le. abs( pv* epsmac**( one/ three))) then
!  Derivatives are both close to zero
           msg( lq, j) = 2
!
        else
!  Derivatives are not both close to zero
           msg( lq, j) = 3
        endif
!
        return
end subroutine
!DODCHK
subroutine dodchk                                                             &
         ( n, m, np, nq,                                                      &
         isodr, anajac, implct,                                               &
         beta, ifixb,                                                         &
         ldx, ldifx, ldscld, ldstpd, ldwe, ld2we, ldwd, ld2wd,                &
         ldy,                                                                 &
         lwork, lwkmn, liwork, liwkmn,                                        &
         sclb, scld, stpb, stpd,                                              &
         info,                                                                &
         lower, upper)
!***Begin Prologue  DODCHK
!***Refer to  ODR
!***Routines Called  (None)
!***Date Written   860529   (YYMMDD)
!***Revision Date  920619   (YYMMDD)
!***Purpose  Check input parameters, indicating errors found using
!       nonzero values of argument INFO
!***End Prologue  DODCHK
!
!...Used modules
        use odrpack_kinds,only: wp
!
!...Scalar arguments
        integer                                                               &
         info, ldifx, ldscld, ldstpd, ldwd, ldwe, ldx, ldy, ld2wd, ld2we,     &
         liwkmn, liwork, lwkmn, lwork, m, n, np, nq
        logical                                                               &
         anajac, implct, isodr
!
!...Array arguments
        real(kind = wp)                                                       &
         beta( np), lower( np), sclb( np), scld( ldscld, m), stpb( np),       &
         stpd( ldstpd, m), upper( np)
        integer                                                               &
         ifixb( np)
!
!...Local scalars
        integer                                                               &
         i, j, k, last, npp
!
!...Variable Definitions (alphabetically)
!       ANAJAC:  The variable designating whether the Jacobians are
!       computed by finite differences (ANAJAC=FALSE) or not
!       (ANAJAC=TRUE).
!       I:       An indexing variable.
!       IFIXB:   The values designating whether the elements of BETA are
!       fixed at their input values or not.
!       IMPLCT:  The variable designating whether the solution is by
!       implicit ODR (IMPLCT=TRUE) or explicit ODR (IMPLCT=FALSE).
!       INFO:    The variable designating why the computations were stopped.
!       ISODR:   The variable designating whether the solution is by ODR
!       (ISODR=TRUE) or by OLS (ISODR=FALSE).
!       J:       An indexing variable.
!       K:       An indexing variable.
!       LAST:    The last row of the array to be accessed.
!       LDIFX:   The leading dimension of array IFIXX.
!       LDSCLD:  The leading dimension of array SCLD.
!       LDSTPD:  The leading dimension of array STPD.
!       LDWD:    The leading dimension of array WD.
!       LDWE:    The leading dimension of array WE.
!       LDX:     The leading dimension of array X.
!       LDY:     The leading dimension of array X.
!       LD2WD:   The second dimension of array WD.
!       LD2WE:   The second dimension of array WE.
!       LIWKMN:  The minimum acceptable length of array IWORK.
!       LIWORK:  The length of vector IWORK.
!       LWKMN:   The minimum acceptable length of array WORK.
!       LWORK:   The length of vector WORK.
!       M:       The number of columns of data in the explanatory variable.
!       N:       The number of observations.
!       NP:      The number of function parameters.
!       NPP:     The number of function parameters being estimated.
!       NQ:      The number of responses per observations.
!       SCLB:    The scaling values for BETA.
!       SCLD:    The scaling value for DELTA.
!       STPB:    The step for the finite difference derivitive wrt BETA.
!       STPD:    The step for the finite difference derivitive wrt DELTA.
!
!
!***First executable statement  DODCHK
!
!
!  Find actual number of parameters being estimated
!
        if ( np .le. 0 .or. ifixb(1) .lt. 0) then
           npp = np
        else
           npp = 0
           do 10 k = 1, np
              if ( ifixb( k) .ne. 0) then
                 npp = npp+1
              endif
10         continue
        endif
!
!  Check problem specification parameters
!
        if ( n .le. 0 .or. m .le. 0 .or.                                      &
            ( npp .le. 0 .or. npp .gt. n) .or.                                &
            ( nq .le. 0)) then
!
           info = 10000
           if ( n .le. 0) then
              info = info+1000
           endif
           if ( m .le. 0) then
              info = info+100
           endif
           if ( npp .le. 0 .or. npp .gt. n) then
              info = info+10
           endif
           if ( nq .le. 0) then
              info = info+1
           endif
!
           return
!
        endif
!
!  Check dimension specification parameters
!
        if ((.not. implct .and. ldy .lt. n) .or.                              &
            ( ldx .lt. n) .or.                                                &
            ( ldwe .ne. 1 .and. ldwe .lt. n) .or.                             &
            ( ld2we .ne. 1 .and. ld2we .lt. nq) .or.                          &
            ( isodr .and.                                                     &
             ( ldwd .ne. 1 .and. ldwd .lt. n)) .or.                           &
            ( isodr .and.                                                     &
             ( ld2wd .ne. 1 .and. ld2wd .lt. m)) .or.                         &
            ( isodr .and.                                                     &
             ( ldifx .ne. 1 .and. ldifx .lt. n)) .or.                         &
            ( isodr .and.                                                     &
             ( ldstpd .ne. 1 .and. ldstpd .lt. n)) .or.                       &
            ( isodr .and.                                                     &
             ( ldscld .ne. 1 .and. ldscld .lt. n)) .or.                       &
            ( lwork .lt. lwkmn) .or.                                          &
            ( liwork .lt. liwkmn)) then
!
           info = 20000
           if (.not. implct .and. ldy .lt. n) then
              info = info+1000
           endif
           if ( ldx .lt. n) then
              info = info+2000
           endif
!
           if (( ldwe .ne. 1 .and. ldwe .lt. n) .or.                          &
               ( ld2we .ne. 1 .and. ld2we .lt. nq)) then
              info = info+100
           endif
           if ( isodr .and.                                                   &
               (( ldwd .ne. 1 .and. ldwd .lt. n) .or.                         &
                ( ld2wd .ne. 1 .and. ld2wd .lt. m))) then
              info = info+200
           endif
!
           if ( isodr .and.                                                   &
               ( ldifx .ne. 1 .and. ldifx .lt. n)) then
              info = info+10
           endif
           if ( isodr .and.                                                   &
               ( ldstpd .ne. 1 .and. ldstpd .lt. n)) then
              info = info+20
           endif
           if ( isodr .and.                                                   &
               ( ldscld .ne. 1 .and. ldscld .lt. n)) then
              info = info+40
           endif
!
           if ( lwork .lt. lwkmn) then
              info = info+1
           endif
           if ( liwork .lt. liwkmn) then
              info = info+2
           endif
           return
!
        endif
!
!  Check DELTA scaling
!
        if ( isodr .and. scld(1,1) .gt. 0) then
!---------------------------------------^--------------------------------------
!!! FPT - 3085 Objects of .EQ., .NE. .GT. etc. are of different data types
!------------------------------------------------------------------------------
           if ( ldscld .ge. n) then
              last = n
           else
              last = 1
           endif
           do 120 j = 1, m
              do 110 i = 1, last
                 if ( scld( i, j) .le. 0) then
!--------------------------------------^---------------------------------------
!!! FPT - 3085 Objects of .EQ., .NE. .GT. etc. are of different data types
!------------------------------------------------------------------------------
                    info = 30200
                    goto 130
                 endif
110           continue
120        continue
        endif
130     continue
!
!  Check BETA scaling
!
        if ( sclb(1) .gt. 0) then
!-------------------------^----------------------------------------------------
!!! FPT - 3085 Objects of .EQ., .NE. .GT. etc. are of different data types
!------------------------------------------------------------------------------
           do 210 k = 1, np
              if ( sclb( k) .le. 0) then
!--------------------------------^---------------------------------------------
!!! FPT - 3085 Objects of .EQ., .NE. .GT. etc. are of different data types
!------------------------------------------------------------------------------
                 if ( info .eq. 0) then
                    info = 30100
                 else
                    info = info+100
                 endif
                 goto 220
              endif
210        continue
        endif
220     continue
!
!  Check DELTA finite difference step sizes
!
        if ( anajac .and. isodr .and. stpd(1,1) .gt. 0) then
!----------------------------------------------------^-------------------------
!!! FPT - 3085 Objects of .EQ., .NE. .GT. etc. are of different data types
!------------------------------------------------------------------------------
           if ( ldstpd .ge. n) then
              last = n
           else
              last = 1
           endif
           do 320 j = 1, m
              do 310 i = 1, last
                 if ( stpd( i, j) .le. 0) then
!--------------------------------------^---------------------------------------
!!! FPT - 3085 Objects of .EQ., .NE. .GT. etc. are of different data types
!------------------------------------------------------------------------------
                    if ( info .eq. 0) then
                       info = 32000
                    else
                       info = info+2000
                    endif
                    goto 330
                 endif
310           continue
320        continue
        endif
330     continue
!
!  Check BETA finite difference step sizes
!
        if ( anajac .and. stpb(1) .gt. 0) then
!--------------------------------------^---------------------------------------
!!! FPT - 3085 Objects of .EQ., .NE. .GT. etc. are of different data types
!------------------------------------------------------------------------------
           do 410 k = 1, np
              if ( stpb( k) .le. 0) then
!--------------------------------^---------------------------------------------
!!! FPT - 3085 Objects of .EQ., .NE. .GT. etc. are of different data types
!------------------------------------------------------------------------------
                 if ( info .eq. 0) then
                    info = 31000
                 else
                    info = info+1000
                 endif
                 goto 420
              endif
410        continue
        endif
420     continue
!
!  Check bounds
!
        if (any( upper(1: np) .lt. lower(1: np))) then
           if ( info .eq. 0) then
              info = 91000
           endif
        endif
!
        if (any(( upper(1: np) .lt. beta(1: np) .or. lower(1: np) .gt. beta(1 &
         : np)) .and. .not. upper(1: np) .lt. lower(1: np))) then
           if ( info .ge. 90000) then
              info = info+100
           else
              info = 90100
           endif
        endif
!
        return
end subroutine
!DODCNT
subroutine dodcnt                                                             &
         ( fcn, n, m, np, nq, beta, y, ldy, x, ldx,                           &
         we, ldwe, ld2we, wd, ldwd, ld2wd, ifixb, ifixx, ldifx,               &
         job, ndigit, taufac, sstol, partol, maxit, iprint, lunerr, lunrpt,   &
         stpb, stpd, ldstpd, sclb, scld, ldscld,                              &
         work, lwork, iwork, liwork,                                          &
         info,                                                                &
         lower, upper)
!***Begin Prologue  DODCNT
!***Refer to  ODR
!***Routines Called  DODDRV
!***Date Written   860529   (YYMMDD)
!***Revision Date  920304   (YYMMDD)
!***Purpose  REAL (KIND=wp) driver routine for finding
!       the weighted explicit or implicit orthogonal distance
!       regression (ODR) or ordinary linear or nonlinear least
!       squares (OLS) solution
!***End Prologue  DODCNT
!
!...Used modules
        use odrpack_kinds,only: wp
!
!...Scalar arguments
        real(kind = wp)                                                       &
         partol, sstol, taufac
        integer                                                               &
         info, iprint, job, ldifx, ldscld, ldstpd, ldwd, ldwe, ldx, ldy,      &
         ld2wd, ld2we, liwork, lunerr, lunrpt, lwork, m, maxit, n, ndigit, np &
         , nq
!
!...Array arguments
        real(kind = wp)                                                       &
         beta( np), lower( np), sclb( np), scld( ldscld, m), stpb( np),       &
         stpd( ldstpd, m), upper( np), wd( ldwd, ld2wd, m),                   &
         we( ldwe, ld2we, nq), work( lwork), x( ldx, m), y( ldy, nq)
        integer                                                               &
         ifixb( np), ifixx( ldifx, m), iwork( liwork)
!
!...Subroutine arguments
        external                                                              &
         fcn
!
!...Local scalars
        real(kind = wp)                                                       &
         cnvtol, one, pcheck, pfac, pstart, three, tstimp, zero
        integer                                                               &
         iprnti, ipr1, ipr2, ipr2f, ipr3, jobi, job1, job2, job3, job4, job5, &
         maxiti, maxit1
        logical                                                               &
         done, fstitr, head, implct, prtpen
!
!...Local arrays
        real(kind = wp)                                                       &
         pnlty(1,1,1)
!
!...External subroutines
        external                                                              &
         doddrv
!
!...External functions
!
!...Data statements
        data                                                                  &
         pcheck, pstart, pfac, zero, one, three                               &
         /1.0E3_wp,1.0E1_wp,1.0E1_wp,0.0E0_wp,1.0E0_wp,3.0E0_wp/
!
!...Routine names used as subprogram arguments
!       FCN:     The user-supplied subroutine for evaluating the model.
!
!...Variable Definitions (alphabetically)
!       BETA:    The function parameters.
!       CNVTOL:  The convergence tolerance for implicit models.
!       DONE:    The variable designating whether the inplicit solution has
!       been found (DONE=TRUE) or not (DONE=FALSE).
!       FSTITR:  The variable designating whether this is the first
!       iteration (FSTITR=TRUE) or not (FSTITR=FALSE).
!       HEAD:    The variable designating whether the heading is to be
!       printed (HEAD=TRUE) or not (HEAD=FALSE).
!       IFIXB:   The values designating whether the elements of BETA are
!       fixed at their input values or not.
!       IFIXX:   The values designating whether the elements of X are
!       fixed at their input values or not.
!       IMPLCT:  The variable designating whether the solution is by
!       implicit ODR (IMPLCT=TRUE) or explicit ODR (IMPLCT=FALSE).
!       INFO:    The variable designating why the computations were stopped.
!       IPRINT:  The print control variables.
!       IPRNTI:  The print control variables.
!       IPR1:    The 1st digit of the print control variable.
!       IPR2:    The 2nd digit of the print control variable.
!       IPR3:    The 3rd digit of the print control variable.
!       IPR4:    The 4th digit of the print control variable.
!       IWORK:   The integer work space.
!       JOB:     The variable controling problem initialization and
!       computational method.
!       JOBI:    The variable controling problem initialization and
!       computational method.
!       JOB1:    The 1st digit of the variable controling problem
!       initialization and computational method.
!       JOB2:    The 2nd digit of the variable controling problem
!       initialization and computational method.
!       JOB3:    The 3rd digit of the variable controling problem
!       initialization and computational method.
!       JOB4:    The 4th digit of the variable controling problem
!       initialization and computational method.
!       JOB5:    The 5th digit of the variable controling problem
!       initialization and computational method.
!       LDIFX:   The leading dimension of array IFIXX.
!       LDSCLD:  The leading dimension of array SCLD.
!       LDSTPD:  The leading dimension of array STPD.
!       LDWD:    The leading dimension of array WD.
!       LDWE:    The leading dimension of array WE.
!       LDX:     The leading dimension of array X.
!       LDY:     The leading dimension of array Y.
!       LD2WD:   The second dimension of array WD.
!       LD2WE:   The second dimension of array WE.
!       LIWORK:  The length of vector IWORK.
!       LOWER:   The lower bound for BETA.
!       LUNERR:  The logical unit number used for error messages.
!       LUNRPT:  The logical unit number used for computation reports.
!       LWORK:   The length of vector work.
!       M:       The number of columns of data in the independent variable.
!       MAXIT:   The maximum number of iterations allowed.
!       MAXITI:  For implicit models, the number of iterations allowed for
!       The current penalty parameter value.
!       MAXIT1:  For implicit models, the number of iterations allowed for
!       the next penalty parameter value.
!       N:       The number of observations.
!       NDIGIT:  The number of accurate digits in the function results, as
!       supplied by the user.
!       NP:      The number of function parameters.
!       NQ:      The number of responses per observation.
!       ONE:     The value 1.0E0_wp.
!       PARTOL:  The user supplied parameter convergence stopping tolerance.
!       PCHECK:  The value designating the minimum penalty parameter allowed
!       before the implicit problem can be considered solved.
!       PFAC:    The factor for increasing the penalty parameter.
!       PNLTY:   The penalty parameter for an implicit model.
!       PRTPEN:  The value designating whether the penalty parameter is to be
!       printed in the iteration report (PRTPEN=TRUE) or not
!       (PRTPEN=FALSE).
!       PSTART:  The factor for increasing the penalty parameter.
!       SCLB:    The scaling values for BETA.
!       SCLD:    The scaling values for DELTA.
!       STPB:    The relative step for computing finite difference
!       Derivatives with respect to BETA.
!       STPD:    The relative step for computing finite difference
!       Derivatives with respect to DELTA.
!       SSTOL:   The sum-of-squares convergence stopping tolerance.
!       TAUFAC:  The factor used to compute the initial trust region
!       diameter.
!       THREE:   The value 3.0E0_wp.
!       TSTIMP:  The relative change in the parameters between the initial
!       values and the solution.
!       UPPER:   The upper bound for BETA.
!       WD:      The DELTA weights.
!       WE:      The EPSILON weights.
!       WORK:    The REAL (KIND=wp) work space.
!       X:       The independent variable.
!       Y:       The dependent variable.  Unused when the model is implicit.
!       ZERO:    The value 0.0E0_wp.
!
!
!***First executable statement  DODCNT
!
!
        implct = mod( job,10) .eq. 1
        fstitr = .true.
        head = .true.
        prtpen = .false.
!
        if ( implct) then
!
!  Set up for implicit problem
!
           if ( iprint .ge. 0) then
              ipr1 = mod( iprint,10000)/1000
              ipr2 = mod( iprint,1000)/100
              ipr2f = mod( iprint,100)/10
              ipr3 = mod( iprint,10)
           else
              ipr1 = 2
              ipr2 = 0
              ipr2f = 0
              ipr3 = 1
           endif
           iprnti = ipr1*1000+ ipr2*100+ ipr2f*10
!
           job5 = mod( job,100000)/10000
           job4 = mod( job,10000)/1000
           job3 = mod( job,1000)/100
           job2 = mod( job,100)/10
           job1 = mod( job,10)
           jobi = job5*10000+ job4*1000+ job3*100+ job2*10+ job1
!
           if ( we(1,1,1) .le. zero) then
              pnlty(1,1,1) = - pstart
           else
              pnlty(1,1,1) = - we(1,1,1)
           endif
!
           if ( partol .lt. zero) then
              cnvtol = epsilon( zero)**( one/ three)
!--------------------------------------^---------------------------------------
!!! FPT - 3437 Mixed real or complex sizes in expression - loss of precision
!------------------------------------------------------------------------------
           else
              cnvtol = min( partol, one)
           endif
!
           if ( maxit .ge. 1) then
              maxiti = maxit
           else
              maxiti = 100
           endif
!
           done = maxiti .eq. 0
           prtpen = .true.
!
10         continue
           call doddrv                                                        &
            ( head, fstitr, prtpen,                                           &
            fcn, n, m, np, nq, beta, y, ldy, x, ldx,                          &
            pnlty,1,1, wd, ldwd, ld2wd, ifixb, ifixx, ldifx,                  &
            jobi, ndigit, taufac, sstol, cnvtol, maxiti,                      &
            iprnti, lunerr, lunrpt,                                           &
            stpb, stpd, ldstpd, sclb, scld, ldscld,                           &
            work, lwork, iwork, liwork,                                       &
            maxit1, tstimp, info, lower, upper)
!
           if ( done) then
              return
           else
              done = maxit1 .le. 0 .or.                                       &
               (abs( pnlty(1,1,1)) .ge. pcheck .and.                          &
               tstimp .le. cnvtol)
           endif
!
           if ( done) then
              if ( tstimp .le. cnvtol) then
                 info = ( info/10)*10+2
              else
                 info = ( info/10)*10+4
              endif
              jobi = 10000+1000+ job3*100+ job2*10+ job1
              maxiti = 0
              iprnti = ipr3
           else
              prtpen = .true.
              pnlty(1,1,1) = pfac* pnlty(1,1,1)
              jobi = 10000+1000+000+ job2*10+ job1
              maxiti = maxit1
              iprnti = 0000+ ipr2*100+ ipr2f*10
           endif
           goto 10
        else
           call doddrv                                                        &
            ( head, fstitr, prtpen,                                           &
            fcn, n, m, np, nq, beta, y, ldy, x, ldx,                          &
            we, ldwe, ld2we, wd, ldwd, ld2wd, ifixb, ifixx, ldifx,            &
            job, ndigit, taufac, sstol, partol, maxit,                        &
            iprint, lunerr, lunrpt,                                           &
            stpb, stpd, ldstpd, sclb, scld, ldscld,                           &
            work, lwork, iwork, liwork,                                       &
            maxit1, tstimp, info, lower, upper)
        endif

end subroutine

subroutine doddrv                                                             &
         ( head, fstitr, prtpen,                                              &
         fcn, n, m, np, nq, beta, y, ldy, x, ldx,                             &
         we, ldwe, ld2we, wd, ldwd, ld2wd, ifixb, ifixx, ldifx,               &
         job, ndigit, taufac, sstol, partol, maxit,                           &
         iprint, lunerr, lunrpt,                                              &
         stpb, stpd, ldstpd, sclb, scld, ldscld,                              &
         work, lwork, iwork, liwork,                                          &
         maxit1, tstimp, info, lower, upper)
!***Begin Prologue  DODDRV
!***Refer to  ODR
!***Routines Called  FCN,DCOPY,DDOT,DETAF,DFCTRW,DFLAGS,
!       DINIWK,DIWINF,DJCK,DNRM2,DODCHK,DODMN,
!       DODPER,DPACK,DSETN,DUNPAC,DWGHT,DWINF,DXMY
!       DERSTEP
!***Date Written   860529   (YYMMDD)
!***Revision Date  920619   (YYMMDD)
!***Purpose  Perform error checking and initialization, and begin
!       procedure for performing orthogonal distance regression
!       (ODR) or ordinary linear or nonlinear least squares (OLS)
!***End Prologue  DODDRV
!
!...Used modules
        use odrpack_kinds,only: wp
        use odrpack,only: tempret
!
!...Scalar arguments
        real(kind = wp)                                                       &
         partol, sstol, taufac, tstimp
        integer                                                               &
         info, iprint, job, ldifx, ldscld, ldstpd, ldwd, ldwe, ldx, ldy,      &
         ld2wd, ld2we, liwork, lunerr, lunrpt, lwork, m, maxit, maxit1,       &
         n, ndigit, np, nq
        logical                                                               &
         fstitr, head, prtpen
!
!...Array arguments
        real(kind = wp)                                                       &
         beta( np), lower( np), sclb( np), scld( ldscld, m), stpb( np),       &
         stpd( ldstpd, m), upper( np), we( ldwe, ld2we, nq),                  &
         wd( ldwd, ld2wd, m), work( lwork), x( ldx, m), y( ldy, nq)
        integer                                                               &
         ifixb( np), ifixx( ldifx, m), iwork( liwork)
!
!...Subroutine arguments
        external                                                              &
         fcn
!
!...Local scalars
        real(kind = wp)                                                       &
         epsmac, eta, p5, one, ten, zero
        integer                                                               &
         actrsi, alphai, betaci, betani, betasi, beta0i, boundi, deltai,      &
         deltni, deltsi, diffi, epsmai, etai, fi, fjacbi, fjacdi, fni, fsi, i &
         , idfi, int2i, iprini, iranki, istop, istopi, jobi, jpvti, k, ldtt,  &
         ldtti, liwkmn, loweri, luneri, lunrpi, lwkmn, lwrk, maxiti, msgb,    &
         msgd, neta, netai, nfev, nfevi, niteri, njev, njevi, nnzw, nnzwi,    &
         npp, nppi, nrow, nrowi, ntol, ntoli, olmavi, omegai, partli, pnormi, &
         prersi, qrauxi, rcondi, rnorsi, rvari, sdi, si, ssfi, ssi, sstoli,   &
         taufci, taui, ti, tti, ui, upperi, vcvi, we1i, wrk1i, wrk2i, wrk3i,  &
         wrk4i, wrk5i, wrk6i, wrk7i, wrk, wssi, wssdei, wssepi, xplusi
        logical                                                               &
         anajac, cdjac, chkjac, dovcv, implct, initd, isodr, redoj, restrt
!
!...Local arrays
        real(kind = wp)                                                       &
         betaj( np)
        integer                                                               &
         interval( np)
!
!...External functions
        real(kind = wp)                                                       &
         ddot, dnrm2, derstep
        external                                                              &
         ddot, dnrm2, derstep
!
!...External subroutines
        external                                                              &
         dcopy, detaf, dfctrw, dflags, diniwk, diwinf, djck, dodchk,          &
         dodmn, dodper, dpack, dsetn, dunpac, dwinf, dxmy
!
!...Data statements
        data                                                                  &
         zero, p5, one, ten                                                   &
         /0.0E0_wp,0.5E0_wp,1.0E0_wp,10.0E0_wp/
!
!...Interface blocks
        interface
           subroutine dwght                                                   &
            ( n, m, wt, ldwt, ld2wt, t, wtt)
           use odrpack_kinds,only: wp
           integer                                                            &
            ldwt, ld2wt, m, n
           real(kind = wp)                                                    &
            t(:,:), wt(:,:,:), wtt(:,:)
           end subroutine
        end interface
!
!...Routine names used as subprogram arguments
!       FCN:     THE USER SUPPLIED SUBROUTINE FOR EVALUATING THE MODEL.
!
!...Variable Definitions (alphabetically)
!       ACTRSI:  The location in array work of variable ACTRS.
!       ALPHAI:  The location in array work of variable ALPHA.
!       ANAJAC:  The variable designating whether the Jacobians are
!       computed by finite differences (ANAJAC=FALSE) or not
!       (ANAJAC=TRUE).
!       BETA:    The function parameters.
!       BETACI:  The starting location in array WORK of array BETAC.
!       BETAJ:   The parameters to use in the derivative checking algorithm.
!       BETANI:  The starting location in array WORK of array BETAN.
!       BETASI:  The starting location in array WORK of array BETAS.
!       BETA0I:  The starting location in array WORK of array BETA0.
!       CDJAC:   The variable designating whether the Jacobians are
!       Computed by central differences (CDJAC=TRUE) or forward
!       differences (CDJAC=FALSE).
!       CHKJAC:  The variable designating whether the user supplied
!       Jacobians are to be checked (CHKJAC=TRUE) or not
!       (CHKJAC=FALSE).
!       DELTAI:  The starting location in array WORK of array DELTA.
!       DELTNI:  The starting location in array WORK of array DELTAN.
!       DELTSI:  The starting location in array WORK of array DELTAS.
!       DIFFI:   The starting location in array WORK of array DIFF.
!       DOVCV:   The variable designating whether the covariance matrix is
!       to be computed (DOVCV=TRUE) or not (DOVCV=FALSE).
!       EPSMAI:  The location in array WORK of variable EPSMAC.
!       ETA:     The relative noise in the function results.
!       ETAI:    The location in array WORK of variable ETA.
!       FI:      The starting location in array WORK of array F.
!       FJACBI:  The starting location in array WORK of array FJACB.
!       FJACDI:  The starting location in array WORK of array FJACD.
!       FNI:     The starting location in array WORK of array FN.
!       FSI:     The starting location in array WORK of array FS.
!       FSTITR:  The variable designating whether this is the first
!       iteration (FSTITR=TRUE) or not (FSTITR=FALSE).
!       HEAD:    The variable designating whether the heading is to be
!       printed (HEAD=TRUE) or not (HEAD=FALSE).
!       I:       An index variable.
!       IDFI:    The location in array iwork of variable IDF.
!       IFIXB:   The values designating whether the elements of BETA are
!       fixed at their input values or not.
!       IFIXX:   The values designating whether the elements of X are
!       fixed at their input values or not.
!       IMPLCT:  The variable designating whether the solution is by
!       implicit ODR (IMPLCT=TRUE) or explicit ODR (IMPLCT=FALSE).
!       INFO:    The variable designating why the computations were stopped.
!       INITD:   The variable designating whether DELTA is to be initialized
!       to zero (INITD=TRUE) or to the values in the first N by M
!       elements of array WORK (INITD=FALSE).
!       INT2I:   The location in array IWORK of variable INT2.
!       INTERVAL: Specifies which checks can be performed when checking derivatives
!       based on the interval of the bound constraints.
!       IPRINI:  The location in array iwork of variable IPRINT.
!       IPRINT:  The print control variable.
!       IRANKI:  The location in array IWORK of variable IRANK.
!       ISODR:   The variable designating whether the solution is by ODR
!       (ISODR=TRUE) or by OLS (ISODR=FALSE).
!       ISTOP:   The variable designating whether there are problems
!       computing the function at the current BETA and DELTA.
!       ISTOPI:  The location in array IWORK of variable ISTOP.
!       IWORK:   The integer work space.
!       JOB:     The variable controling problem initialization and
!       computational method.
!       JOBI:    The location in array IWORK of variable JOB.
!       JPVTI:   The starting location in array IWORK of array JPVT.
!       K:       An index variable.
!       LDIFX:   The leading dimension of array IFIXX.
!       LDSCLD:  The leading dimension of array SCLD.
!       LDSTPD:  The leading dimension of array STPD.
!       LDTT:    The leading dimension of array TT.
!       LDTTI:   The location in array IWORK of variable LDTT.
!       LDWD:    The leading dimension of array WD.
!       LDWE:    The leading dimension of array WE.
!       LDX:     The leading dimension of array X.
!       LDY:     The leading dimension of array Y.
!       LD2WD:   The second dimension of array WD.
!       LD2WE:   The second dimension of array WE.
!       LIWKMN:  The minimum acceptable length of array IWORK.
!       LIWORK:  The length of vector IWORK.
!       LOWER:   The lower bound for BETA.
!       LUNERI:  The location in array IWORK of variable LUNERR.
!       LUNERR:  The logical unit number used for error messages.
!       LUNRPI:  The location in array IWORK of variable LUNRPT.
!       LUNRPT:  The logical unit number used for computation reports.
!       LWKMN:   The minimum acceptable length of array WORK.
!       LWORK:   The length of vector WORK.
!       LWRK:    The length of vector WRK.
!       M:       The number of columns of data in the explanatory variable.
!       MAXIT:   The maximum number of iterations allowed.
!       MAXIT1:  For implicit models, the iterations allowed for the next
!       penalty parameter value.
!       MAXITI:  The location in array IWORK of variable MAXIT.
!       MSGB:    The starting location in array IWORK of array MSGB.
!       MSGD:    The starting location in ARRAY IWORK of array MSGD.
!       N:       The number of observations.
!       NDIGIT:  The number of accurate digits in the function results, as
!       supplied by the user.
!       NETA:    The number of accurate digits in the function results.
!       NETAI:   The location in array IWORK of variable NETA.
!       NFEV:    The number of function evaluations.
!       NFEVI:   The location in array IWORK of variable NFEV.
!       NITERI:  The location in array IWORK of variable NITER.
!       NJEV:    The number of Jacobian evaluations.
!       NJEVI:   The location in array IWORK of variable NJEV.
!       NNZW:    The number of nonzero observational error weights.
!       NNZWI:   The location in array IWORK of variable NNZW.
!       NP:      The number of function parameters.
!       NPP:     The number of function parameters being estimated.
!       NPPI:    The location in array IWORK of variable NPP.
!       NQ:      The number of responses per observation.
!       NROW:    The row number at which the derivative is to be checked.
!       NROWI:   The location in array IWORK of variable NROW.
!       NTOL:    The number of digits of agreement required between the
!       numerical derivatives and the user supplied derivatives,
!       set by DJCK.
!       NTOLI:   The location in array IWORK of variable NTOL.
!       OLMAVI:  The location in array WORK of variable OLMAVG.
!       OMEGAI:  The starting location in array WORK of array OMEGA.
!       ONE:     The value 1.0E0_wp.
!       PARTLI:  The location in array WORK of variable PARTOL.
!       PARTOL:  The parameter convergence stopping tolerance.
!       PNORM:   The norm of the scaled estimated parameters.
!       PNORMI:  The location in array WORK of variable PNORM.
!       PRERSI:  The location in array WORK of variable PRERS.
!       PRTPEN:  The variable designating whether the penalty parameter is
!       to be printed in the iteration report (PRTPEN=TRUE) or not
!       (PRTPEN=FALSE).
!       P5:      The value 0.5E0_wp.
!       QRAUXI:  The starting location in array WORK of array QRAUX.
!       RCONDI:  The location in array WORK of variable RCOND.
!       REDOJ:   The variable designating whether the Jacobian matrix is to
!       be recomputed for the computation of the covariance matrix
!       (REDOJ=TRUE) or not (REDOJ=FALSE).
!       RESTRT:  The variable designating whether the call is a restart
!       (RESTRT=TRUE) or not (RESTRT=FALSE).
!       RNORSI:  The location in array WORK of variable RNORMS.
!       RVARI:   The location in array WORK of variable RVAR.
!       SCLB:    The scaling values for BETA.
!       SCLD:    The scaling values for DELTA.
!       SDI:     The starting location in array WORK of array SD.
!       SI:      The starting location in array WORK of array S.
!       SSFI:    The starting location in array WORK of array SSF.
!       SSI:     The starting location in array WORK of array SS.
!       SSTOL:   The sum-of-squares convergence stopping tolerance.
!       SSTOLI:  The location in array WORK of variable SSTOL.
!       STPB:    The step size for finite difference derivatives wrt BETA.
!       STPD:    The step size for finite difference derivatives wrt DELTA.
!       TAUFAC:  The factor used to compute the initial trust region
!       diameter.
!       TAUFCI:  The location in array WORK of variable TAUFAC.
!       TAUI:    The location in array WORK of variable TAU.
!       TEN:     The value 10.0E0_wp.
!       TI:      The starting location in array WORK of array T.
!       TSTIMP:  The relative change in the parameters between the initial
!       values and the solution.
!       TTI:     The starting location in array WORK of array TT.
!       UI:      The starting location in array WORK of array U.
!       UPPER:   The upper bound for BETA.
!       VCVI:    The starting location in array WORK of array VCV.
!       WD:      The DELTA weights.
!       WE:      The EPSILON weights.
!       WE1I:    The starting location in array WORK of array WE1.
!       WORK:    The REAL (KIND=wp) work space.
!       WRK:     The starting location in array WORK of array WRK,
!       equivalenced to WRK1 and WRK2.
!       WRK1I:   The starting location in array WORK of array WRK1.
!       WRK2I:   The starting location in array WORK of array WRK2.
!       WRK3I:   The starting location in array WORK of array WRK3.
!       WRK4I:   The starting location in array WORK of array WRK4.
!       WRK5I:   The starting location in array WORK of array WRK5.
!       WRK6I:   The starting location in array WORK of array WRK6.
!       WRK7I:   The starting location in array WORK of array WRK7.
!       WSSI:    The location in array WORK of variable wss.
!       WSSDEI:  The location in array WORK of variable WSSDEL.
!       WSSEPI:  The location in array WORK of variable WSSEPS.
!       X:       The explanatory variable.
!       XPLUSI:  The starting location in array WORK of array XPLUSD.
!       Y:       The dependent variable.  Unused when the model is implicit.
!       ZERO:    The value 0.0E0_wp.
!
!
!***First executable statement  DODDRV
!
!
!  Initialize necessary variables
!
        call dflags( job, restrt, initd, dovcv, redoj,                        &
         anajac, cdjac, chkjac, isodr, implct)
!
!  Set starting locations within integer workspace
!  (invalid values of M, NP and/or NQ are handled reasonably by DIWINF)
!
        call diwinf( m, np, nq,                                               &
         msgb, msgd, jpvti, istopi,                                           &
         nnzwi, nppi, idfi,                                                   &
         jobi, iprini, luneri, lunrpi,                                        &
         nrowi, ntoli, netai,                                                 &
         maxiti, niteri, nfevi, njevi, int2i, iranki, ldtti,                  &
         boundi,                                                              &
         liwkmn)
!
!  Set starting locations within REAL (KIND=wp) work space
!  (invalid values of N, M, NP, NQ, LDWE and/or LD2WE
!  are handled reasonably by DWINF)
!
        call dwinf( n, m, np, nq, ldwe, ld2we, isodr,                         &
         deltai, fi, xplusi, fni, sdi, vcvi,                                  &
         rvari, wssi, wssdei, wssepi, rcondi, etai,                           &
         olmavi, taui, alphai, actrsi, pnormi, rnorsi, prersi,                &
         partli, sstoli, taufci, epsmai,                                      &
         beta0i, betaci, betasi, betani, si, ssi, ssfi, qrauxi, ui,           &
         fsi, fjacbi, we1i, diffi,                                            &
         deltsi, deltni, ti, tti, omegai, fjacdi,                             &
         wrk1i, wrk2i, wrk3i, wrk4i, wrk5i, wrk6i, wrk7i,                     &
         loweri, upperi,                                                      &
         lwkmn)
        if ( isodr) then
           wrk = wrk1i
           lwrk = n* m* nq+ n* nq
        else
           wrk = wrk2i
           lwrk = n* nq
        endif
!
!  Update the penalty parameters
!  (WE(1,1,1) is not a user supplied array in this case)
        if ( restrt .and. implct) then
           we(1,1,1) = max( work( we1i)**2,abs( we(1,1,1)))
           work( we1i) = -sqrt(abs( we(1,1,1)))
        endif
!
        if ( restrt) then
!
!  Reset maximum number of iterations
!
           if ( maxit .ge. 0) then
              iwork( maxiti) = iwork( niteri)+ maxit
           else
              iwork( maxiti) = iwork( niteri)+10
           endif
!
           if ( iwork( niteri) .lt. iwork( maxiti)) then
              info = 0
           endif
!
           if ( job .ge. 0) iwork( jobi) = job
           if ( iprint .ge. 0) iwork( iprini) = iprint
           if ( partol .ge. zero .and. partol .lt. one) work( partli) =       &
            partol
           if ( sstol .ge. zero .and. sstol .lt. one) work( sstoli) = sstol
!
           work( olmavi) = work( olmavi)* iwork( niteri)
!
           if ( implct) then
              call dcopy( n* nq, work( fni),1, work( fi),1)
           else
              call dxmy( n, nq, work( fni), n, y, ldy, work( fi), n)
           endif
           call dwght( n, nq,                                                 &
            reshape( work( we1i: we1i+ ldwe* ld2we* nq-1),(/ ldwe, ld2we, nq  &
            /)), ldwe, ld2we,reshape( work( fi: fi+ n* nq-1),(/ n, nq/)),     &
            tempret(1: n,1: nq))
           work( fi: fi+ n* nq-1) = reshape( tempret(1: n,1: nq),(/ n* nq/))
           work( wssepi) = ddot( n* nq, work( fi),1, work( fi),1)
           work( wssi) = work( wssepi)+ work( wssdei)
!
        else
!
!  Perform error checking
!
           info = 0
!
           call dodchk( n, m, np, nq,                                         &
            isodr, anajac, implct,                                            &
            beta, ifixb,                                                      &
            ldx, ldifx, ldscld, ldstpd, ldwe, ld2we, ldwd, ld2wd,             &
            ldy,                                                              &
            lwork, lwkmn, liwork, liwkmn,                                     &
            sclb, scld, stpb, stpd,                                           &
            info,                                                             &
            lower, upper)
           if ( info .gt. 0) then
              goto 50
           endif
!
!  Initialize work vectors as necessary
!
           do 10 i = n* m+ n* nq+1, lwork
              work( i) = zero
10         continue
           do 20 i = 1, liwork
              iwork( i) = 0
20         continue
!
           call diniwk( n, m, np,                                             &
            work, lwork, iwork, liwork,                                       &
            x, ldx, ifixx, ldifx, scld, ldscld,                               &
            beta, sclb,                                                       &
            sstol, partol, maxit, taufac,                                     &
            job, iprint, lunerr, lunrpt,                                      &
            lower, upper,                                                     &
            epsmai, sstoli, partli, maxiti, taufci,                           &
            jobi, iprini, luneri, lunrpi,                                     &
            ssfi, tti, ldtti, deltai,                                         &
            loweri, upperi, boundi)
!
           iwork( msgb) = -1
           iwork( msgd) = -1
           work( taui) = - work( taufci)
!
!  Set up for parameter estimation -
!  Pull BETA's to be estimated and corresponding scale values
!  and store in WORK(BETACI) and WORK(SSI), respectively
!
           call dpack( np, iwork( nppi), work( betaci), beta, ifixb)
           call dpack( np, iwork( nppi), work( ssi), work( ssfi), ifixb)
           npp = iwork( nppi)
!
!  Check that WD is positive definite and WE is positive semidefinite,
!  saving factorization of WE, and counting number of nonzero weights
!
           call dfctrw( n, m, nq, npp,                                        &
            isodr,                                                            &
            we, ldwe, ld2we, wd, ldwd, ld2wd,                                 &
            work( wrk2i), work( wrk4i),                                       &
            work( we1i), nnzw, info)
           iwork( nnzwi) = nnzw
!
           if ( info .ne. 0) then
              goto 50
           endif
!
!  Evaluate the predicted values and
!          weighted EPSILONS at the starting point
!
           call dunpac( np, work( betaci), beta, ifixb)
           work(xplusi:xplusi+(n*m-1)) = work(deltai:deltai+(n*m-1)) + reshape(x(1:n,:), shape=[n*m])
           istop = 0
           call fcn( n, m, np, nq,                                            &
            n, m, np,                                                         &
            beta, work(xplusi),                                              &
            ifixb, ifixx, ldifx,                                              &
            002, work( fni), work( wrk6i), work( wrk1i),                      &
            istop)
           iwork( istopi) = istop
           if ( istop .eq. 0) then
              iwork( nfevi) = iwork( nfevi)+1
              if ( implct) then
                 call dcopy( n* nq, work( fni),1, work( fi),1)
              else
                 call dxmy( n, nq, work( fni), n, y, ldy, work( fi), n)
              endif
              call dwght( n, nq,                                              &
               reshape( work( we1i: we1i+ ldwe* ld2we* nq-1),                 &
               (/ ldwe, ld2we, nq/)), ldwe, ld2we,                            &
               reshape( work( fi: fi+ n* nq-1),(/ n, nq/)),                   &
               tempret(1: n,1: nq))
              work( fi: fi+ n* nq-1) = reshape( tempret(1: n,1: nq),(/ n* nq  &
               /))
           else
              info = 52000
              goto 50
           endif
!
!  Compute norm of the initial estimates
!
           call dwght( npp,1,reshape( work( ssi: ssi+ npp-1),(/ npp,1,1/)),   &
            npp,1,reshape( work( betaci: betaci+ npp-1),(/ npp,1/)),          &
            tempret(1: npp,1:1))
           work( wrk: wrk+ npp-1) = tempret(1: npp,1)
           if ( isodr) then
              call dwght( n, m,reshape( work( tti: tti+ iwork( ldtti)*1* m-1) &
               ,(/ iwork( ldtti),1, m/)), iwork( ldtti),1,reshape( work(      &
               deltai: deltai+ n* m-1),(/ n, m/)), tempret(1: n,1: m))
              work( wrk+ npp: wrk+ npp+ n* m-1) =                             &
               reshape( tempret(1: n,1: m),(/ n* m/))
              work( pnormi) = dnrm2( npp+ n* m, work( wrk),1)
           else
              work( pnormi) = dnrm2( npp, work( wrk),1)
           endif
!
!  Compute sum of squares of the weighted EPSILONS and weighted DELTAS
!
           work( wssepi) = ddot( n* nq, work( fi),1, work( fi),1)
           if ( isodr) then
              call dwght( n, m, wd, ldwd, ld2wd,                              &
               reshape( work( deltai: deltai+ n* m),(/ n, m/)),               &
               tempret(1: n,1: m))
              work( wrk: wrk+ n* m-1) = reshape( tempret(1: n,1: m),(/ n* m/) &
               )
              work( wssdei) = ddot( n* m, work( deltai),1, work( wrk),1)
           else
              work( wssdei) = zero
           endif
           work( wssi) = work( wssepi)+ work( wssdei)
!
!  Select first row of X + DELTA that contains no zeros
!
           nrow = -1
           call dsetn( n, m, work( xplusi), n, nrow)
           iwork( nrowi) = nrow
!
!  Set number of good digits in function results
!
           epsmac = work( epsmai)
           if ( ndigit .lt. 2) then
              iwork( netai) = -1
              nfev = iwork( nfevi)
              call detaf( fcn,                                                &
               n, m, np, nq,                                                  &
               work( xplusi), beta, epsmac, nrow,                             &
               work( betani), work( fni),                                     &
               ifixb, ifixx, ldifx,                                           &
               istop, nfev, eta, neta,                                        &
               work( wrk1i), work( wrk2i), work( wrk6i), work( wrk7i),        &
               info,                                                          &
               lower, upper)
              iwork( istopi) = istop
              iwork( nfevi) = nfev
              if ( istop .ne. 0 .or. info .ne. 0) then
                 if ( info .eq. 0) then
                    info = 53000
                 endif
                 iwork( netai) = 0
                 work( etai) = zero
                 goto 50
              else
                 iwork( netai) = - neta
                 work( etai) = eta
              endif
           else
              iwork( netai) = min( ndigit,int( p5-log10( epsmac)))
              work( etai) = max( epsmac, ten**(- ndigit))
           endif
!
!  Check bounds are large enough for derivative calculations.
!
           if (.not. anajac .or. chkjac) then
              if ( cdjac) then
                 do k = 1, np
                    if ( upper( k)-abs(2* derstep(1, k, upper( k), work( ssfi &
                     ), stpb, neta)) .lt. lower( k)) then
                       info = 90020
                       goto 50
                    endif
                 enddo
              else
                 do k = 1, np
                    if ( upper( k)-abs(2* derstep(0, k, upper( k), work( ssfi &
                     ), stpb, neta)) .lt. lower( k)) then
                       info = 90020
                       goto 50
                    endif
                 enddo
              endif
           endif
!
!  CHECK DERIVATIVES IF NECESSARY
!
           if ( chkjac .and. anajac) then
              ntol = -1
              nfev = iwork( nfevi)
              njev = iwork( njevi)
              neta = iwork( netai)
              ldtt = iwork( ldtti)
              eta = work( etai)
              epsmac = work( epsmai)
!  ENSURE BETA IS NOT TOO CLOSE TO BOUNDS FOR THE DERIVATIVE CHECK.
              betaj(:) = beta(:)
              call mbfb( np, betaj, lower, upper, work( ssfi), stpb, neta,    &
               eta, interval)
!  CHECK THE DERIVATIVES.
              call djck( fcn,                                                 &
               n, m, np, nq,                                                  &
               beta, betaj, work( xplusi),                                    &
               ifixb, ifixx, ldifx, stpb, stpd, ldstpd,                       &
               work( ssfi), work( tti), ldtt,                                 &
               eta, neta, ntol, nrow, isodr, epsmac,                          &
               work( fni), work( fjacbi), work( fjacdi),                      &
               iwork( msgb), iwork( msgd), work( diffi),                      &
               istop, nfev, njev,                                             &
               work( wrk1i), work( wrk2i), work( wrk6i),                      &
               interval)
              iwork( istopi) = istop
              iwork( nfevi) = nfev
              iwork( njevi) = njev
              iwork( ntoli) = ntol
              if ( istop .ne. 0) then
                 info = 54000
              elseif ( iwork( msgb) .ne. 0 .or. iwork( msgd) .ne. 0) then
                 info = 40000
              endif
           else
!
!  Indicate user supplied derivatives were not checked
              iwork( msgb) = -1
              iwork( msgd) = -1
           endif
!
!  Print appropriate error messages
!
50         if (( info .ne. 0) .or.                                            &
               ( iwork( msgb) .ne. -1)) then
              if ( lunerr .ne. 0 .and. iprint .ne. 0) then
                 call dodper                                                  &
                  ( info, lunerr,                                             &
                  n, m, np, nq,                                               &
                  ldscld, ldstpd, ldwe, ld2we, ldwd, ld2wd,                   &
                  lwkmn, liwkmn,                                              &
                  work( fjacbi), work( fjacdi),                               &
                  work( diffi), iwork( msgb), isodr, iwork( msgd),            &
                  work( xplusi), iwork( nrowi), iwork( netai), iwork( ntoli))
              endif
!
!  Set INFO to reflect errors in the user supplied Jacobians
!
              if ( info .eq. 40000) then
                 if ( iwork( msgb) .eq. 2 .or. iwork( msgd) .eq. 2) then
                    if ( iwork( msgb) .eq. 2) then
                       info = info+1000
                    endif
                    if ( iwork( msgd) .eq. 2) then
                       info = info+100
                    endif
                 else
                    info = 0
                 endif
              endif
              if ( info .ne. 0) then
                 return
              endif
           endif
        endif
!
!  Save the initial values of BETA
        call dcopy( np, beta,1, work( beta0i),1)
!
!  Find least squares solution
!
        call dcopy( n* nq, work( fni),1, work( fsi),1)
        ldtt = iwork( ldtti)
        call dodmn( head, fstitr, prtpen,                                     &
         fcn, n, m, np, nq, job, beta, y, ldy, x, ldx,                        &
         we, work( we1i), ldwe, ld2we, wd, ldwd, ld2wd,                       &
         ifixb, ifixx, ldifx,                                                 &
         work( betaci), work( betani), work( betasi), work( si),              &
         work( deltai), work( deltni), work( deltsi),                         &
         work( loweri), work( upperi),                                        &
         work( ti), work( fi), work( fni), work( fsi),                        &
         work( fjacbi), iwork( msgb), work( fjacdi), iwork( msgd),            &
         work( ssfi), work( ssi), work( tti), ldtt,                           &
         stpb, stpd, ldstpd,                                                  &
         work( xplusi), work( wrk), lwrk,                                     &
         work, lwork, iwork, liwork, info,                                    &
         iwork( boundi))
        maxit1 = iwork( maxiti)- iwork( niteri)
        tstimp = zero
        do 100 k = 1, np
           if ( beta( k) .eq. zero) then
!-----------------------------^------------------------------------------------
!!! FPT - 3087 REAL or COMPLEX quantity tested for exact equality/inequality
!------------------------------------------------------------------------------
              tstimp = max( tstimp,                                           &
               abs( beta( k)- work( beta0i-1+ k))/ work( ssfi-1+ k))
           else
              tstimp = max( tstimp,                                           &
               abs( beta( k)- work( beta0i-1+ k))/abs( beta( k)))
           endif
100     continue
!
        return
!
end subroutine
!DODLM
subroutine dodlm                                                              &
         ( n, m, np, nq, npp,                                                 &
         f, fjacb, fjacd,                                                     &
         wd, ldwd, ld2wd, ss, tt, ldtt, delta,                                &
         alpha2, tau, epsfcn, isodr,                                          &
         tfjacb, omega, u, qraux, jpvt,                                       &
         s, t, nlms, rcond, irank,                                            &
         wrk1, wrk2, wrk3, wrk4, wrk5, wrk, lwrk, istopc)
!***Begin Prologue  DODLM
!***Refer to  ODR
!***Routines Called  DDOT,DNRM2,DODSTP,DSCALE,DWGHT
!***Date Written   860529   (YYMMDD)
!***Revision Date  920619   (YYMMDD)
!***Purpose  Compute Levenberg-Marquardt parameter and steps S AND T
!       using analog of the trust-region Levenberg-Marquardt
!       algorithm
!***End Prologue  DODLM
!
!...Used modules
        use odrpack_kinds,only: wp
        use odrpack,only: tempret
!
!...Scalar arguments
        real(kind = wp)                                                       &
         alpha2, epsfcn, rcond, tau
        integer                                                               &
         irank, istopc, ldtt, ldwd, ld2wd, lwrk, m, n, nlms, np, npp, nq
        logical                                                               &
         isodr
!
!...Array arguments
        real(kind = wp)                                                       &
         delta( n, m), f( n, nq), fjacb( n, np, nq), fjacd( n, m, nq),        &
         omega( nq, nq), qraux( np), s( np), ss( np),                         &
         t( n, m), tfjacb( n, nq, np), tt( ldtt, m), u( np), wd( ldwd, ld2wd, &
         m), wrk( lwrk), wrk1( n, nq, m), wrk2( n, nq), wrk3( np), wrk4( m, m &
         ), wrk5( m)
        integer                                                               &
         jpvt( np)
!
!...Local scalars
        real(kind = wp)                                                       &
         alpha1, alphan, bot, p001, p1, phi1, phi2, sa, top, zero
        integer                                                               &
         i, iwrk, j, k, l
        logical                                                               &
         forvcv
!
!...External functions
        real(kind = wp)                                                       &
         ddot, dnrm2
        external                                                              &
         ddot, dnrm2
!
!...External subroutines
        external                                                              &
         dodstp, dscale
!
!...Data statements
        data                                                                  &
         zero, p001, p1                                                       &
         /0.0E0_wp,0.001E0_wp,0.1E0_wp/
!
!...Interface blocks
        interface
           subroutine dwght                                                   &
            ( n, m, wt, ldwt, ld2wt, t, wtt)
           use odrpack_kinds,only: wp
           integer                                                            &
            ldwt, ld2wt, m, n
           real(kind = wp)                                                    &
            t(:,:), wt(:,:,:), wtt(:,:)
           end subroutine
        end interface
!
!...Variable Definitions (alphabetically)
!       ALPHAN:  The new Levenberg-Marquardt parameter.
!       ALPHA1:  The previous Levenberg-Marquardt parameter.
!       ALPHA2:  The current Levenberg-Marquardt parameter.
!       BOT:     The lower limit for setting ALPHA.
!       DELTA:   The estimated errors in the explanatory variables.
!       EPSFCN:  The function's precision.
!       F:       The (weighted) estimated values of EPSILON.
!       FJACB:   The Jacobian with respect to BETA.
!       FJACD:   The Jacobian with respect to DELTA.
!       FORVCV:  The variable designating whether this subroutine was
!       called to set up for the covariance matrix computations
!       (FORVCV=TRUE) or not (FORVCV=FALSE).
!       I:       An indexing variable.
!       IRANK:   The rank deficiency of the Jacobian wrt BETA.
!       ISODR:   The variable designating whether the solution is by ODR
!       (ISODR=TRUE) or by OLS (ISODR=FALSE).
!       ISTOPC:  The variable designating whether the computations were
!       stoped due to some numerical error detected within
!       subroutine DODSTP.
!       IWRK:    An indexing variable.
!       J:       An indexing variable.
!       K:       An indexing variable.
!       L:       An indexing variable.
!       JPVT:    The pivot vector.
!       LDTT:    The leading dimension of array TT.
!       LDWD:    The leading dimension of array WD.
!       LD2WD:   The second dimension of array WD.
!       LWRK:    The length of vector WRK.
!       M:       The number of columns of data in the explanatory variable.
!       N:       The number of observations.
!       NLMS:    The number of Levenberg-Marquardt steps taken.
!       NP:      The number of function parameters.
!       NPP:     The number of function parameters being estimated.
!       NQ:      The number of responses per observation.
!       OMEGA:   The array (I-FJACD*INV(P)*trans(FJACD))**(-1/2)  where
!       P = trans(FJACD)*FJACD + D**2 + ALPHA*TT**2
!       P001:    The value 0.001E0_wp
!       P1:      The value 0.1E0_wp
!       PHI1:    The previous difference between the norm of the scaled step
!       and the trust region diameter.
!       PHI2:    The current difference between the norm of the scaled step
!       and the trust region diameter.
!       QRAUX:   The array required to recover the orthogonal part of the
!       Q-R decomposition.
!       RCOND:   The approximate reciprocal condition of TFJACB.
!       S:       The step for BETA.
!       SA:      The scalar PHI2*(ALPHA1-ALPHA2)/(PHI1-PHI2).
!       SS:      The scaling values used for the unfixed BETAS.
!       T:       The step for DELTA.
!       TAU:     The trust region diameter.
!       TFJACB:  The array OMEGA*FJACB.
!       TOP:     The upper limit for setting ALPHA.
!       TT:      The scale used for the DELTA'S.
!       U:       The approximate null vector for TFJACB.
!       WD:      The DELTA weights.
!       WRK:     A work array of (LWRK) elements,
!       equivalenced to WRK1 and WRK2.
!       WRK1:    A work array of (N by NQ by M) elements.
!       WRK2:    A work array of (N by NQ) elements.
!       WRK3:    A work array of (NP) elements.
!       WRK4:    A work array of (M by M) elements.
!       WRK5:    A work array of (M) elements.
!       ZERO:    The value 0.0E0_wp.
!
!
!***First executable statement  DODLM
!
        forvcv = .false.
        istopc = 0
!
!  Compute full Gauss-Newton step (ALPHA=0)
!
        alpha1 = zero
        call dodstp( n, m, np, nq, npp,                                       &
         f, fjacb, fjacd,                                                     &
         wd, ldwd, ld2wd, ss, tt, ldtt, delta,                                &
         alpha1, epsfcn, isodr,                                               &
         tfjacb, omega, u, qraux, jpvt,                                       &
         s, t, phi1, irank, rcond, forvcv,                                    &
         wrk1, wrk2, wrk3, wrk4, wrk5, wrk, lwrk, istopc)
        if ( istopc .ne. 0) then
           return
        endif
!
!  Initialize TAU if necessary
!
        if ( tau .lt. zero) then
           tau = abs( tau)* phi1
        endif
!
!  Check if full Gauss-Newton step is optimal
!
        if (( phi1- tau) .le. p1* tau) then
           nlms = 1
           alpha2 = zero
           return
        endif
!
!  Full Gauss-Newton step is outside trust region -
!  find locally constrained optimal step
!
        phi1 = phi1- tau
!
!  Initialize upper and lower bounds for ALPHA
!
        bot = zero
!
        do 30 k = 1, npp
           do 20 l = 1, nq
              do 10 i = 1, n
                 tfjacb( i, l, k) = fjacb( i, k, l)
10            continue
20         continue
           wrk( k) = ddot( n* nq, tfjacb(1,1, k),1, f(1,1),1)
30      continue
        call dscale( npp,1, ss, npp, wrk, npp, wrk, npp)
!
        if ( isodr) then
           call dwght( n, m, wd, ldwd, ld2wd, delta, tempret(1: n,1: m))
           wrk( npp+1: npp+1+ n* m-1) = reshape( tempret(1: n,1: m),(/ n* m/) &
            )
           iwrk = npp
           do 50 j = 1, m
              do 40 i = 1, n
                 iwrk = iwrk+1
                 wrk( iwrk) = wrk( iwrk)+                                     &
                  ddot( nq, fjacd( i, j,1), n* m, f( i,1), n)
40            continue
50         continue
           call dscale( n, m, tt, ldtt, wrk( npp+1), n, wrk( npp+1), n)
           top = dnrm2( npp+ n* m, wrk,1)/ tau
        else
           top = dnrm2( npp, wrk,1)/ tau
        endif
!
        if ( alpha2 .gt. top .or. alpha2 .eq. zero) then
!---------------------------------------------^--------------------------------
!!! FPT - 3087 REAL or COMPLEX quantity tested for exact equality/inequality
!------------------------------------------------------------------------------
           alpha2 = p001* top
        endif
!
!  Main loop
!
        do 60 i = 1,10
!
!  Compute locally constrained steps S and T and PHI(ALPHA) for
!  current value of ALPHA
!
           call dodstp( n, m, np, nq, npp,                                    &
            f, fjacb, fjacd,                                                  &
            wd, ldwd, ld2wd, ss, tt, ldtt, delta,                             &
            alpha2, epsfcn, isodr,                                            &
            tfjacb, omega, u, qraux, jpvt,                                    &
            s, t, phi2, irank, rcond, forvcv,                                 &
            wrk1, wrk2, wrk3, wrk4, wrk5, wrk, lwrk, istopc)
           if ( istopc .ne. 0) then
              return
           endif
           phi2 = phi2- tau
!
!  Check whether current step is optimal
!
           if (abs( phi2) .le. p1* tau .or.                                   &
               ( alpha2 .eq. bot .and. phi2 .lt. zero)) then
!----------------------------^-------------------------------------------------
!!! FPT - 3087 REAL or COMPLEX quantity tested for exact equality/inequality
!------------------------------------------------------------------------------
              nlms = i+1
              return
           endif
!
!  Current step is not optimaL
!
!  Update bounds for ALPHA and compute new ALPHA
!
           if ( phi1- phi2 .eq. zero) then
!-------------------------------^----------------------------------------------
!!! FPT - 3087 REAL or COMPLEX quantity tested for exact equality/inequality
!------------------------------------------------------------------------------
              nlms = 12
              return
           endif
           sa = phi2*( alpha1- alpha2)/( phi1- phi2)
           if ( phi2 .lt. zero) then
              top = min( top, alpha2)
           else
              bot = max( bot, alpha2)
           endif
           if ( phi1* phi2 .gt. zero) then
              bot = max( bot, alpha2- sa)
           else
              top = min( top, alpha2- sa)
           endif
!
           alphan = alpha2- sa*( phi1+ tau)/ tau
           if ( alphan .ge. top .or. alphan .le. bot) then
              alphan = max( p001* top,sqrt( top* bot))
           endif
!
!  Get ready for next iteration
!
           alpha1 = alpha2
           alpha2 = alphan
           phi1 = phi2
60      continue
!
!  Set NLMS to indicate an optimal step could not be found in 10 trys
!
        nlms = 12
!
        return
end subroutine

subroutine dodmn                                                              &
         ( head, fstitr, prtpen,                                              &
         fcn, n, m, np, nq, job, beta, y, ldy, x, ldx,                        &
         we, we1, ldwe, ld2we, wd, ldwd, ld2wd,                               &
         ifixb, ifixx, ldifx,                                                 &
         betac, betan, betas, s, delta, deltan, deltas,                       &
         lower, upper,                                                        &
         t, f, fn, fs, fjacb, msgb, fjacd, msgd,                              &
         ssf, ss, tt, ldtt, stpb, stpd, ldstpd,                               &
         xplusd, wrk, lwrk, work, lwork, iwork, liwork, info,                 &
         bound)
!***Begin Prologue  DODMN
!***Refer to  ODR
!***Routines Called  FCN,DACCES,DCOPY,DDOT,DEVJAC,DFLAGS,DNRM2,DODLM,
!       DODPCR,DODVCV,DUNPAC,DWGHT,DXMY
!***Date Written   860529   (YYMMDD)
!***Revision Date  920619   (YYMMDD)
!***Purpose  Iteratively compute least squares solution
!***End Prologue  DODMN
!
!...Used modules
        use odrpack_kinds, only: wp, zero, one
        use odrpack,only: tempret
!
!...Scalar arguments
        integer                                                               &
         info, job, ldifx, ldstpd, ldtt, ldwd, ldwe, ldx, ldy, ld2wd, ld2we,  &
         liwork, lwork, lwrk, m, n, np, nq
!
!...Array arguments
        real(kind = wp)                                                       &
         beta( np), betac( np), betan( np), betas( np),                       &
         delta( n, m), deltan( n, m), deltas( n, m),                          &
         f( n, nq), fjacb( n, np, nq), fjacd( n, m, nq), fn( n, nq), fs( n,   &
         nq), lower( np), s( np), ss( np), ssf( np), stpb( np), stpd( ldstpd, &
         m), t( n, m), tt( ldtt, m), upper( np), wd( ldwd, ld2wd, m), we(     &
         ldwe, ld2we, nq), we1( ldwe, ld2we, nq), work( lwork), x( ldx, m),   &
         xplusd( n, m), wrk( lwrk), y( ldy, nq)
        integer                                                               &
         bound( np), ifixb( np), ifixx( ldifx, m), iwork( liwork),            &
         msgb( nq* np+1), msgd( nq* m+1)
        logical                                                               &
         fstitr, head, prtpen

!...Subroutine arguments
        external                                                              &
         fcn

!...Local scalars
        real(kind = wp) :: actred, actrs, alpha, dirder, eta, olmavg,         &
         p0001, p1, p25, p5, p75, partol, pnorm, prered, prers,               &
         ratio, rcond, rnorm, rnormn, rnorms, rss, rvar, sstol, tau, taufac,  &
         temp, temp1, temp2, tsnorm
        integer                                                               &
         i, idf, iflag, int2, ipr, ipr1, ipr2, ipr2f, ipr3, irank,            &
!---------------------------^--------------------------------------------------
!!! FPT - 1271 Non-standard Fortran intrinsic(s) used as local identifier(s)
!------------------------------------------------------------------------------
         istop, istopc, iwrk, j, jpvt, l, looped, ludflt, lunr, lunrpt,       &
         maxit, neta, nfev, niter, njev, nlms, nnzw, npp, npr, npu, omega,    &
         qraux, sd, u, vcv, wrk1, wrk2, wrk3, wrk4, wrk5, wrk6
        logical                                                               &
         access, anajac, cdjac, chkjac, cnvpar, cnvss, didvcv, dovcv,         &
!--------------^---------------------------------------------------------------
!!! FPT - 1273 Fortran auxiliary keyword used as identifier name.
!------------------------------------------------------------------------------
         implct, initd, intdbl, isodr, lstep, redoj, restrt

!...Local arrays
        real(kind = wp) :: loweru( np), upperu( np), wss(3)

!...External functions
        real(kind = wp), external :: ddot, dnrm2

!...External subroutines
        external :: dacces, dcopy, devjac, dflags, dodlm, dodpcr, dodvcv, dunpac, dxmy

!...Data statements
        data                                                                  &
         p0001, p1, p25, p5, p75                                              &
         /0.00010E0_wp,0.10E0_wp,0.250E0_wp,0.50E0_wp,0.750E0_wp/
        data                                                                  &
         ludflt                                                               &
         /6/
!
!...Interface blocks
        interface
           subroutine dwght                                                   &
            ( n, m, wt, ldwt, ld2wt, t, wtt)
           use odrpack_kinds,only: wp
           integer                                                            &
            ldwt, ld2wt, m, n
           real(kind = wp)                                                    &
            t(:,:), wt(:,:,:), wtt(:,:)
           end subroutine
        end interface
!
!...Routine names used as subprogram arguments
!       FCN:     The user supplied subroutine for evaluating the model.
!
!...Variable Definitions (alphabetically)
!       ACCESS:  The variable designating whether information is to be
!       accessed from the work arrays (ACCESS=TRUE) or stored in
!       them (ACCESS=FALSE).
!       ACTRED:  The actual relative reduction in the sum-of-squares.
!       ACTRS:   The saved actual relative reduction in the sum-of-squares.
!       ALPHA:   The Levenberg-Marquardt parameter.
!       ANAJAC:  The variable designating whether the Jacobians are computed
!       by finite differences (ANAJAC=FALSE) or not (ANAJAC=TRUE).
!       BETA:    The function parameters.
!       BETAC:   The current estimated values of the unfixed BETA'S.
!       BETAN:   The new estimated values of the unfixed BETA'S.
!       BETAS:   The saved estimated values of the unfixed BETA'S.
!       CDJAC:   The variable designating whether the Jacobians are computed
!       by central differences (cdjac=true) or by forward
!       differences (CDJAC=FALSE).
!       CHKJAC:  The variable designating whether the user supplied
!       Jacobians are to be checked (CHKJAC=TRUE) or not
!       (CHKJAC=FALSE).
!       CNVPAR:  The variable designating whether parameter convergence was
!       attained (CNVPAR=TRUE) or not (CNVPAR=FALSE).
!       CNVSS:   The variable designating whether sum-of-squares convergence
!       was attained (CNVSS=TRUE) or not (CNVSS=FALSE).
!       DELTA:   The estimated errors in the explanatory variables.
!       DELTAN:  The new estimated errors in the explanatory variables.
!       DELTAS:  The saved estimated errors in the explanatory variables.
!       DIDVCV:  The variable designating whether the covariance matrix was
!       computed (DIDVCV=TRUE) or not (DIDVCV=FALSE).
!       DIRDER:  The directional derivative.
!       DOVCV:   The variable designating whether the covariance matrix
!       should to be computed (DOVCV=TRUE) or not (DOVCV=FALSE).
!       ETA:     The relative noise in the function results.
!       F:       The (weighted) estimated values of EPSILON.
!       FJACB:   The Jacobian with respect to BETA.
!       FJACD:   The Jacobian with respect to DELTA.
!       FN:      The new predicted values from the function.
!       FS:      The saved predicted values from the function.
!       FSTITR:  The variable designating whether this is the first
!       iteration (FSTITR=TRUE) or not (FSTITR=FALSE).
!       HEAD:    The variable designating whether the heading is to be
!       printed (HEAD=TRUE) or not (HEAD=FALSE).
!       I:       An indexing variable.
!       IDF:     The degrees of freedom of the fit, equal to the number of
!       observations with nonzero weighted derivatives minus the
!       number of parameters being estimated.
!       IFIXB:   The values designating whether the elements of BETA are
!       fixed at their input values or not.
!       IFIXX:   The values designating whether the elements of X are
!       fixed at their input values or not.
!       IFLAG:   The variable designating which report is to be printed.
!       IMPLCT:  The variable designating whether the solution is by
!       implicit ODR (IMPLCT=TRUE) or explicit ODR (IMPLCT=FALSE).
!       INFO:    The variable designating why the computations were stopped.
!       INITD:   The variable designating whether delta is initialized to
!       zero (INITD=TRUE) or to the values in the first N by M
!       elements of array work (INITD=FALSE).
!       INT2:    The number of internal doubling steps taken.
!       INTDBL:  The variable designating whether internal doubling is to be
!       used (INTDBL=TRUE) or NOT (INTDBL=FALSE).
!       IPR:     The values designating the length of the printed report.
!       IPR1:    The value of the 4th digit (from the right) of iprint,
!       which controls the initial summary report.
!       IPR2:    The value of the 3rd digit (from the right) of iprint,
!       which controls the iteration report.
!       IPR2F:   The value of the 2nd digit (from the right) of iprint,
!       which controls the frequency of the iteration reports.
!       IPR3:    The value of the 1st digit (from the right) of iprint,
!       which controls the final summary report.
!       IRANK:   The rank deficiency of the Jacobian wrt BETA.
!       ISODR:   The variable designating whether the solution is by ODR
!       (ISODR=TRUE) or OLS (ISODR=FALSE).
!       ISTOP:   The variable designating whether there are problems
!       computing the function at the current BETA and DELTA.
!       ISTOPC:  The variable designating whether the computations were
!       stoped due to some numerical error within routine  DODSTP.
!       IWORK:   The integer work space.
!       IWRK:    An index variable.
!       J:       An index variable.
!       JOB:     The variable controling problem initialization and
!       computational method.
!       JPVT:    The starting location in IWORK of array JPVT.
!       L:       An index variable.
!       LDIFX:   The leading dimension of array IFIXX.
!       LDTT:    The leading dimension of array TT.
!       LDWD:    The leading dimension of array WD.
!       LDWE:    The leading dimension of array WE and WE1.
!       LDX:     The leading dimension of array X.
!       LDY:     The leading dimension of array Y.
!       LD2WD:   The second dimension of array WD.
!       LD2WE:   The second dimension of array WE and WE1.
!       LIWORK:  The length of vector IWORK.
!       LOOPED:  A counter used to determine how many times the subloop
!       has been executed, where if the count becomes large
!       enough the computations will be stopped.
!       LOWERU:  The lower bound for unfixed BETAs.
!       LSTEP:   The variable designating whether a successful step has
!       been found (LSTEP=TRUE) or not (LSTEP=FALSE).
!       LUDFLT:  The default logical unit number, used for computation
!       reports to the screen.
!       LUNR:    The logical unit number used for computation reports.
!       LUNRPT:  The logical unit number used for computation reports.
!       LWORK:   The length of vector WORK.
!       LWRK:    The length of vector WRK.
!       M:       The number of columns of data in the explanatory variable.
!       MAXIT:   The maximum number of iterations allowed.
!       MSGB:    The error checking results for the Jacobian wrt BETA.
!       MSGD:    The error checking results for the Jacobian wrt DELTA.
!       N:       The number of observations.
!       NETA:    The number of accurate digits in the function results.
!       NFEV:    The number of function evaluations.
!       NITER:   The number of iterations taken.
!       NJEV:    The number of Jacobian evaluations.
!       NLMS:    The number of Levenberg-Marquardt steps taken.
!       NNZW:    The number of nonzero weighted observations.
!       NP:      The number of function parameters.
!       NPP:     The number of function parameters being estimated.
!       NPR:     The number of times the report is to be written.
!       NPU:     The number of unfixed parameters.
!       NQ:      The number of responses per observation.
!       OLMAVG:  The average number of Levenberg-Marquardt steps per
!       iteration.
!       OMEGA:   The starting location in WORK of array OMEGA.
!       P0001:   The value 0.0001E0_wp.
!       P1:      The value 0.1E0_wp.
!       P25:     The value 0.25E0_wp.
!       P5:      The value 0.5E0_wp.
!       P75:     The value 0.75E0_wp.
!       PARTOL:  The parameter convergence stopping tolerance.
!       PNORM:   The norm of the scaled estimated parameters.
!       PRERED:  The predicted relative reduction in the sum-of-squares.
!       PRERS:   The old predicted relative reduction in the sum-of-squares.
!       PRTPEN:  The value designating whether the penalty parameter is to
!       be printed in the iteration report (PRTPEN=TRUE) or not
!       (PRTPEN=FALSE).
!       QRAUX:   The starting location in array WORK of array QRAUX.
!       RATIO:   The ratio of the actual relative reduction to the predicted
!       relative reduction in the sum-of-squares.
!       RCOND:   The approximate reciprocal condition of FJACB.
!       REDOJ:   The variable designating whether the Jacobian matrix is to
!       be recomputed for the computation of the covariance matrix
!       (REDOJ=TRUE) or not (REDOJ=FALSE).
!       RESTRT:  The variable designating whether the call is a restart
!       (RESTRT=TRUE) or not (RESTRT=FALSE).
!       RNORM:   The norm of the weighted errors.
!       RNORMN:  The new norm of the weighted errors.
!       RNORMS:  The saved norm of the weighted errors.
!       RSS:     The residual sum of squares.
!       RVAR:    The residual variance.
!       S:       The step for BETA.
!       SD:      The starting location in array work of array SD.
!       SS:      The scaling values used for the unfixed BETAS.
!       SSF:     The scaling values used for BETA.
!       SSTOL:   The sum-of-squares convergence stopping tolerance.
!       STPB:    The relative step used for computing finite difference
!       derivatives with respect to each BETA.
!       STPD:    The relative step used for computing finite difference
!       derivatives with respect to DELTA.
!       T:       The step for DELTA.
!       TAU:     The trust region diameter.
!       TAUFAC:  The factor used to compute the initial trust region
!       diameter.
!       TEMP:    A temporary storage location.
!       TEMP1:   A temporary storage location.
!       TEMP2:   A temporary storage location.
!       TSNORM:  The norm of the scaled step.
!       TT:      The scaling values used for DELTA.
!       U:       The starting location in array WORK of array U.
!       UPPERU:  The upper bound for unfixed BETAs.
!       VCV:     The starting location in array WORK of array VCV.
!       WE:      The EPSILON weights.
!       WE1:     The square root of the EPSILON weights.
!       WD:      The DELTA weights.
!       WORK:    The REAL (KIND=wp) work space.
!       WSS:     The sum-of-squares of the weighted EPSILONS and DELTAS,
!       the sum-of-squares of the weighted DELTAS, and
!       the sum-of-squares of the weighted EPSILONS.
!       WRK:     A work array, equivalenced to WRK1 and WRK2
!       WRK1:    The starting location in array WORK of array WRK1.
!       WRK2:    The starting location in array WORK of array WRK2.
!       WRK3:    The starting location in array WORK of array WRK3.
!       WRK4:    The starting location in array WORK of array WRK4.
!       WRK5:    The starting location in array WORK of array WRK5.
!       WRK6:    The starting location in array WORK of array WRK6.
!       X:       The explanatory variable.
!       XPLUSD:  The values of X + DELTA.
!       Y:       The dependent variable.  Unused when the model is implicit.
!
!
!***First executable statement  DODMN
!
!
!  Initialize necessary variables
!
        call dpack( np, npu, loweru, lower, ifixb)
        call dpack( np, npu, upperu, upper, ifixb)
        call dflags( job, restrt, initd, dovcv, redoj,                        &
         anajac, cdjac, chkjac, isodr, implct)
        access = .true.
        call dacces( n, m, np, nq, ldwe, ld2we,                               &
         work, lwork, iwork, liwork,                                          &
         access, isodr,                                                       &
         jpvt, omega, u, qraux, sd, vcv,                                      &
         wrk1, wrk2, wrk3, wrk4, wrk5, wrk6,                                  &
         nnzw, npp,                                                           &
         job, partol, sstol, maxit, taufac, eta, neta,                        &
         lunrpt, ipr1, ipr2, ipr2f, ipr3,                                     &
         wss, rvar, idf,                                                      &
         tau, alpha, niter, nfev, njev, int2, olmavg,                         &
         rcond, irank, actrs, pnorm, prers, rnorms, istop)
        rnorm = sqrt( wss(1))
!
        didvcv = .false.
        intdbl = .false.
        lstep = .true.
!
!  Print initial summary if desired
!
        if ( ipr1 .ne. 0 .and. lunrpt .ne. 0) then
           iflag = 1
           if ( ipr1 .ge. 3 .and. lunrpt .ne. ludflt) then
              npr = 2
           else
              npr = 1
           endif
           if ( ipr1 .ge. 6) then
              ipr = 2
           else
              ipr = 2-mod( ipr1,2)
           endif
           lunr = lunrpt
           do 10 i = 1, npr
              call dodpcr( ipr, lunr,                                         &
               head, prtpen, fstitr, didvcv, iflag,                           &
               n, m, np, nq, npp, nnzw,                                       &
               msgb, msgd, beta, y, ldy, x, ldx, delta,                       &
               we, ldwe, ld2we, wd, ldwd, ld2wd,                              &
               ifixb, ifixx, ldifx,                                           &
               lower, upper,                                                  &
               ssf, tt, ldtt, stpb, stpd, ldstpd,                             &
               job, neta, taufac, sstol, partol, maxit,                       &
               wss, rvar, idf, work( sd),                                     &
               niter, nfev, njev, actred, prered,                             &
               tau, pnorm, alpha, f, rcond, irank, info, istop)
              if ( ipr1 .ge. 5) then
                 ipr = 2
              else
                 ipr = 1
              endif
              lunr = ludflt
10         continue
!
        endif
!
!  Stop if initial estimates are exact solution
!
        if ( rnorm .eq. zero) then
!-----------------------^------------------------------------------------------
!!! FPT - 3087 REAL or COMPLEX quantity tested for exact equality/inequality
!------------------------------------------------------------------------------
           info = 1
           olmavg = zero
           istop = 0
           goto 150
        endif
!
!  Stop if number of iterations already equals maximum permitted
!
        if ( restrt .and.                                                     &
            ( niter .ge. maxit)) then
           istop = 0
           goto 150
        elseif ( niter .ge. maxit) then
           info = 4
           istop = 0
           goto 150
        endif
!
!  Main loop
!
100     continue
!
        niter = niter+1
        rnorms = rnorm
        looped = 0
!
!  Evaluate jacobian using best estimate of function (FS)
!
        if (( niter .eq. 1) .and.                                             &
            ( anajac .and. chkjac)) then
           istop = 0
        else
           call devjac( fcn,                                                  &
            anajac, cdjac,                                                    &
            n, m, np, nq,                                                     &
            betac, beta, stpb,                                                &
            ifixb, ifixx, ldifx,                                              &
            x, ldx, delta, xplusd, stpd, ldstpd,                              &
            ssf, tt, ldtt, neta, fs,                                          &
            t, work( wrk1), work( wrk2), work( wrk3), work( wrk6),            &
            fjacb, isodr, fjacd, we1, ldwe, ld2we,                            &
            njev, nfev, istop, info,                                          &
            lower, upper)
        endif
        if ( istop .ne. 0) then
           info = 51000
           goto 200
        elseif ( info .eq. 50300) then
           goto 200
        endif
!
!  Sub loop for
!       internal doubling or
!       computing new step when old failed
!
110     continue
!
!  Compute steps S and T
!
        if ( looped .gt. 100) then
           info = 60000
           goto 200
        else
           looped = looped+1
           call dodlm( n, m, np, nq, npp,                                     &
            f, fjacb, fjacd,                                                  &
            wd, ldwd, ld2wd, ss, tt, ldtt, delta,                             &
            alpha, tau, eta, isodr,                                           &
            work( wrk6), work( omega),                                        &
            work( u), work( qraux), iwork( jpvt),                             &
            s, t, nlms, rcond, irank,                                         &
            work( wrk1), work( wrk2), work( wrk3), work( wrk4),               &
            work( wrk5), wrk, lwrk, istopc)
        endif
        if ( istopc .ne. 0) then
           info = istopc
           goto 200
        endif
        olmavg = olmavg+ nlms

!  Compute BETAN = BETAC + S
!       DELTAN = DELTA + T
        betan = betac + s
        if (isodr) deltan = delta + t
!
!  Project the step wrt the bounds
        do i = 1, npu
           if ( loweru( i) .eq. upperu( i)) then
!-------------------------------^----------------------------------------------
!!! FPT - 3087 REAL or COMPLEX quantity tested for exact equality/inequality
!------------------------------------------------------------------------------
              betan( i) = upperu( i)
              s( i) = upperu( i)- betac( i)
              bound( i) = 3
           elseif ( betan( i) .le. loweru( i)) then
              betan( i) = loweru( i)
              s( i) = loweru( i)- betac( i)
              bound( i) = 2
           elseif ( betan( i) .ge. upperu( i)) then
              betan( i) = upperu( i)
              s( i) = upperu( i)- betac( i)
              bound( i) = 1
           else
              bound( i) = 0
           endif
        enddo
!
!  Compute norm of scaled steps S and T (TSNORM)
!
        call dwght( npp,1,reshape( ss,(/ npp,1,1/)), npp,1,                   &
         reshape( s,(/ npp,1/)), tempret(1: npp,1:1))
        wrk(1: npp) = tempret(1: npp,1)
        if ( isodr) then
           call dwght( n, m,reshape( tt,(/ ldtt,1, m/)), ldtt,1,              &
            t, tempret(1: n,1: m))
           wrk( npp+1: npp+1+ n* m-1) = reshape( tempret(1: n,1: m),(/ n* m/) &
            )
           tsnorm = dnrm2( npp+ n* m, wrk,1)
        else
           tsnorm = dnrm2( npp, wrk,1)
        endif
!
!  Compute scaled predicted reduction
!
        iwrk = 0
        do 130 l = 1, nq
           do 120 i = 1, n
              iwrk = iwrk+1
              wrk( iwrk) = ddot( npp, fjacb( i,1, l), n, s,1)
              if ( isodr) wrk( iwrk) = wrk( iwrk)+                            &
               ddot( m, fjacd( i,1, l), n, t( i,1), n)
120        continue
130     continue
        if ( isodr) then
           call dwght( n, m, wd, ldwd, ld2wd, t, tempret(1: n,1: m))
           wrk( n* nq+1: n* nq+1+ n* m-1) = reshape( tempret(1: n,1: m),(/ n* &
            m/))
           temp1 = ddot( n* nq, wrk,1, wrk,1)+ ddot( n* m, t,1, wrk( n* nq+1) &
            ,1)
           temp1 = sqrt( temp1)/ rnorm
        else
           temp1 = dnrm2( n* nq, wrk,1)/ rnorm
        endif
        temp2 = sqrt( alpha)* tsnorm/ rnorm
        prered = temp1**2+ temp2**2/ p5
!
        dirder = -( temp1**2+ temp2**2)
!
!  Evaluate predicted values at new point
!
        call dunpac( np, betan, beta, ifixb)
        xplusd = x(1:n,:) + deltan
        istop = 0
        call fcn( n, m, np, nq,                                               &
         n, m, np,                                                            &
         beta, xplusd,                                                        &
         ifixb, ifixx, ldifx,                                                 &
         002, fn, work( wrk6), work( wrk1),                                   &
         istop)
        if ( istop .eq. 0) then
           nfev = nfev+1
        endif
!
        if ( istop .lt. 0) then
!
!  Set INFO to indicate user has stopped the computations in FCN
!
           info = 51000
           goto 200
        elseif ( istop .gt. 0) then
!
!  Set norm to indicate step should be rejected
!
           rnormn = rnorm/( p1* p75)
        else
!
!  Compute norm of new weighted EPSILONS and weighted DELTAS (RNORMN)
!
           if ( implct) then
              call dcopy( n* nq, fn,1, wrk,1)
           else
              call dxmy( n, nq, fn, n, y, ldy, wrk, n)
           endif
           call dwght( n, nq, we1, ldwe, ld2we,reshape( wrk,(/ n, nq/)),      &
            tempret(1: n,1: nq))
           wrk(1: n* nq) = reshape( tempret(1: n,1: nq),(/ n* nq/))
           if ( isodr) then
              call dwght( n, m, wd, ldwd, ld2wd, deltan, tempret(1: n,1: m))
              wrk( n* nq+1: n* nq+1+ n* m-1) = reshape( tempret(1: n,1: m),(/ &
               n* m/))
              rnormn = sqrt( ddot( n* nq, wrk,1, wrk,1)+                      &
               ddot( n* m, deltan,1, wrk( n* nq+1),1))
           else
              rnormn = dnrm2( n* nq, wrk,1)
           endif
        endif
!
!  Compute scaled actual reduction
!
        if ( p1* rnormn .lt. rnorm) then
           actred = one-( rnormn/ rnorm)**2
        else
           actred = - one
        endif
!
!  Compute ratio of actual reduction to predicted reduction
!
        if ( prered .eq. zero) then
!------------------------^-----------------------------------------------------
!!! FPT - 3087 REAL or COMPLEX quantity tested for exact equality/inequality
!------------------------------------------------------------------------------
           ratio = zero
        else
           ratio = actred/ prered
        endif
!
!  Check on lack of reduction in internal doubling case
!
        if ( intdbl .and.                                                     &
            ( ratio .lt. p0001 .or. rnormn .gt. rnorms)) then
           istop = 0
           tau = tau* p5
           alpha = alpha/ p5
           call dcopy( npp, betas,1, betan,1)
           call dcopy( n* m, deltas,1, deltan,1)
           call dcopy( n* nq, fs,1, fn,1)
           actred = actrs
           prered = prers
           rnormn = rnorms
           ratio = p5
        endif
!
!  Update step bound
!
        intdbl = .false.
        if ( ratio .lt. p25) then
           if ( actred .ge. zero) then
              temp = p5
           else
              temp = p5* dirder/( dirder+ p5* actred)
           endif
           if ( p1* rnormn .ge. rnorm .or. temp .lt. p1) then
              temp = p1
           endif
           tau = temp*min( tau, tsnorm/ p1)
           alpha = alpha/ temp
!
        elseif ( alpha .eq. zero) then
!---------------------------^--------------------------------------------------
!!! FPT - 3087 REAL or COMPLEX quantity tested for exact equality/inequality
!------------------------------------------------------------------------------
           tau = tsnorm/ p5
!
        elseif ( ratio .ge. p75 .and. nlms .le. 11) then
!
!  Step qualifies for internal doubling
!          - Update TAU and ALPHA
!          - Save information for current point
!
           intdbl = .true.
!
           tau = tsnorm/ p5
           alpha = alpha* p5
!
           call dcopy( npp, betan,1, betas,1)
           call dcopy( n* m, deltan,1, deltas,1)
           call dcopy( n* nq, fn,1, fs,1)
           actrs = actred
           prers = prered
           rnorms = rnormn
        endif
!
!  If internal doubling, skip convergence checks
!
        if ( intdbl .and. tau .gt. zero) then
           int2 = int2+1
           goto 110
        endif
!
!  Check acceptance
!
        if ( ratio .ge. p0001) then
           call dcopy( n* nq, fn,1, fs,1)
           if ( implct) then
              call dcopy( n* nq, fs,1, f,1)
           else
              call dxmy( n, nq, fs, n, y, ldy, f, n)
           endif
           call dwght( n, nq, we1, ldwe, ld2we, f, tempret(1: n,1: nq))
           f(1: n,1: nq) = tempret(1: n,1: nq)
           call dcopy( npp, betan,1, betac,1)
           call dcopy( n* m, deltan,1, delta,1)
           rnorm = rnormn
           call dwght( npp,1,reshape( ss,(/ npp,1,1/)), npp,1,                &
            reshape( betac,(/ npp,1/)), tempret(1: npp,1:1))
           wrk(1: npp) = tempret(1: npp,1)
           if ( isodr) then
              call dwght( n, m,reshape( tt,(/ ldtt,1, m/)), ldtt,1,           &
               delta, tempret(1: n,1: m))
              wrk( npp+1: npp+1+ n* m-1) = reshape( tempret(1: n,1: m),(/ n*  &
               m/))
              pnorm = dnrm2( npp+ n* m, wrk,1)
           else
              pnorm = dnrm2( npp, wrk,1)
           endif
           lstep = .true.
        else
           lstep = .false.
        endif
!
!  TEST CONVERGENCE
!
        info = 0
        cnvss = rnorm .eq. zero                                               &
!--------------------------^---------------------------------------------------
!!! FPT - 3087 REAL or COMPLEX quantity tested for exact equality/inequality
!------------------------------------------------------------------------------
         .or.                                                                 &
         (abs( actred) .le. sstol .and.                                       &
         prered .le. sstol .and.                                              &
         p5* ratio .le. one)
        cnvpar = ( tau .le. partol* pnorm) .and. (.not. implct)
        if ( cnvss) info = 1
        if ( cnvpar) info = 2
        if ( cnvss .and. cnvpar) info = 3
!
!  Print iteration report
!
        if ( info .ne. 0 .or. lstep) then
           if ( ipr2 .ne. 0 .and. ipr2f .ne. 0 .and. lunrpt .ne. 0) then
              if ( ipr2f .eq. 1 .or. mod( niter, ipr2f) .eq. 1) then
                 iflag = 2
                 call dunpac( np, betac, beta, ifixb)
                 wss(1) = rnorm* rnorm
                 if ( ipr2 .ge. 3 .and. lunrpt .ne. ludflt) then
                    npr = 2
                 else
                    npr = 1
                 endif
                 if ( ipr2 .ge. 6) then
                    ipr = 2
                 else
                    ipr = 2-mod( ipr2,2)
                 endif
                 lunr = lunrpt
                 do 140 i = 1, npr
                    call dodpcr( ipr, lunr,                                   &
                     head, prtpen, fstitr, didvcv, iflag,                     &
                     n, m, np, nq, npp, nnzw,                                 &
                     msgb, msgd, beta, y, ldy, x, ldx, delta,                 &
                     we, ldwe, ld2we, wd, ldwd, ld2wd,                        &
                     ifixb, ifixx, ldifx,                                     &
                     lower, upper,                                            &
                     ssf, tt, ldtt, stpb, stpd, ldstpd,                       &
                     job, neta, taufac, sstol, partol, maxit,                 &
                     wss, rvar, idf, work( sd),                               &
                     niter, nfev, njev, actred, prered,                       &
                     tau, pnorm, alpha, f, rcond, irank, info, istop)
                    if ( ipr2 .ge. 5) then
                       ipr = 2
                    else
                       ipr = 1
                    endif
                    lunr = ludflt
140              continue
                 fstitr = .false.
                 prtpen = .false.
              endif
           endif
        endif
!
!  Check if finished
!
        if ( info .eq. 0) then
           if ( lstep) then
!
!  Begin next interation unless a stopping criteria has been met
!
              if ( niter .ge. maxit) then
                 info = 4
              else
                 goto 100
              endif
           else
!
!  Step failed - recompute unless a stopping criteria has been met
!
              goto 110
           endif
        endif
!
150     continue
!
        if ( istop .gt. 0) info = info+100
!
!  Store unweighted EPSILONS and X+DELTA to return to user
!
        if ( implct) then
           call dcopy( n* nq, fs,1, f,1)
        else
           call dxmy( n, nq, fs, n, y, ldy, f, n)
        endif
        call dunpac( np, betac, beta, ifixb)
        xplusd = x(1:n, :) + delta
!
!  Compute covariance matrix of estimated parameters
!  in upper NP by NP portion of WORK(VCV) if requested
!
        if ( dovcv .and. istop .eq. 0) then
!
!  Re-evaluate Jacobian at final solution, if requested
!  Otherwise, Jacobian from beginning of last iteration will be used
!  to compute covariance matrix
!
           if ( redoj) then
              call devjac( fcn,                                               &
               anajac, cdjac,                                                 &
               n, m, np, nq,                                                  &
               betac, beta, stpb,                                             &
               ifixb, ifixx, ldifx,                                           &
               x, ldx, delta, xplusd, stpd, ldstpd,                           &
               ssf, tt, ldtt, neta, fs,                                       &
               t, work( wrk1), work( wrk2), work( wrk3), work( wrk6),         &
               fjacb, isodr, fjacd, we1, ldwe, ld2we,                         &
               njev, nfev, istop, info,                                       &
               lower, upper)
!
!
              if ( istop .ne. 0) then
                 info = 51000
                 goto 200
              elseif ( info .eq. 50300) then
                 goto 200
              endif
           endif
!
           if ( implct) then
              call dwght( n, m, wd, ldwd, ld2wd, delta, tempret(1: n,1: m))
              wrk( n* nq+1: n* nq+1+ n* m-1) = reshape( tempret(1: n,1: m),(/ &
               n* m/))
              rss = ddot( n* m, delta,1, wrk( n* nq+1),1)
           else
              rss = rnorm* rnorm
           endif
           if ( redoj .or. niter .ge. 1) then
              call dodvcv( n, m, np, nq, npp,                                 &
               f, fjacb, fjacd,                                               &
               wd, ldwd, ld2wd, ssf, ss, tt, ldtt, delta,                     &
               eta, isodr,                                                    &
               work( vcv), work( sd),                                         &
               work( wrk6), work( omega),                                     &
               work( u), work( qraux), iwork( jpvt),                          &
               s, t, irank, rcond, rss, idf, rvar, ifixb,                     &
               work( wrk1), work( wrk2), work( wrk3), work( wrk4),            &
               work( wrk5), wrk, lwrk, istopc)
              if ( istopc .ne. 0) then
                 info = istopc
                 goto 200
              endif
              didvcv = .true.
           endif
!
        endif
!
!  Set JPVT to indicate dropped, fixed and estimated parameters
!
200     do 210 i = 0, np-1
           work( wrk3+ i) = iwork( jpvt+ i)
           iwork( jpvt+ i) = -2
210     continue
        if ( redoj .or. niter .ge. 1) then
           do 220 i = 0, npp-1
              j = int( work( wrk3+ i))-1
              if ( i .le. npp- irank-1) then
                 iwork( jpvt+ j) = 1
              else
                 iwork( jpvt+ j) = -1
              endif
220        continue
           if ( npp .lt. np) then
              j = npp-1
              do 230 i = np-1,0,-1
                 if ( ifixb( i+1) .eq. 0) then
                    iwork( jpvt+ i) = 0
                 else
                    iwork( jpvt+ i) = iwork( jpvt+ j)
                    j = j-1
                 endif
230           continue
           endif
        endif
!
!  Store various scalars in work arrays for return to user
!
        if ( niter .ge. 1) then
           olmavg = olmavg/ niter
        else
           olmavg = zero
        endif
!
!  Compute weighted sums of squares for return to user
!
        call dwght( n, nq, we1, ldwe, ld2we, f, tempret(1: n,1: nq))
        wrk(1: n* nq) = reshape( tempret(1: n,1: nq),(/ n* nq/))
        wss(3) = ddot( n* nq, wrk,1, wrk,1)
        if ( isodr) then
           call dwght( n, m, wd, ldwd, ld2wd, delta, tempret(1: n,1: m))
           wrk( n* nq+1: n* nq+1+ n* m-1) = reshape( tempret(1: n,1: m),(/ n* &
            m/))
           wss(2) = ddot( n* m, delta,1, wrk( n* nq+1),1)
        else
           wss(2) = zero
        endif
        wss(1) = wss(2)+ wss(3)
!
        access = .false.
        call dacces( n, m, np, nq, ldwe, ld2we,                               &
         work, lwork, iwork, liwork,                                          &
         access, isodr,                                                       &
         jpvt, omega, u, qraux, sd, vcv,                                      &
         wrk1, wrk2, wrk3, wrk4, wrk5, wrk6,                                  &
         nnzw, npp,                                                           &
         job, partol, sstol, maxit, taufac, eta, neta,                        &
         lunrpt, ipr1, ipr2, ipr2f, ipr3,                                     &
         wss, rvar, idf,                                                      &
         tau, alpha, niter, nfev, njev, int2, olmavg,                         &
         rcond, irank, actrs, pnorm, prers, rnorms, istop)
!
!  Encode existance of questionable results into info
!
        if ( info .le. 9 .or. info .ge. 60000) then
           if ( msgb(1) .eq. 1 .or. msgd(1) .eq. 1) then
              info = info+1000
           endif
           if ( istop .ne. 0) then
              info = info+100
           endif
           if ( irank .ge. 1) then
              if ( npp .gt. irank) then
                 info = info+10
              else
                 info = info+20
              endif
           endif
        endif
!
!  Print final summary
!
        if ( ipr3 .ne. 0 .and. lunrpt .ne. 0) then
           iflag = 3
!
           if ( ipr3 .ge. 3 .and. lunrpt .ne. ludflt) then
              npr = 2
           else
              npr = 1
           endif
           if ( ipr3 .ge. 6) then
              ipr = 2
           else
              ipr = 2-mod( ipr3,2)
           endif
           lunr = lunrpt
           do 240 i = 1, npr
              call dodpcr( ipr, lunr,                                         &
               head, prtpen, fstitr, didvcv, iflag,                           &
               n, m, np, nq, npp, nnzw,                                       &
               msgb, msgd, beta, y, ldy, x, ldx, delta,                       &
               we, ldwe, ld2we, wd, ldwd, ld2wd,                              &
               iwork( jpvt), ifixx, ldifx,                                    &
               lower, upper,                                                  &
               ssf, tt, ldtt, stpb, stpd, ldstpd,                             &
               job, neta, taufac, sstol, partol, maxit,                       &
               wss, rvar, idf, work( sd),                                     &
               niter, nfev, njev, actred, prered,                             &
               tau, pnorm, alpha, f, rcond, irank, info, istop)
              if ( ipr3 .ge. 5) then
                 ipr = 2
              else
                 ipr = 1
              endif
              lunr = ludflt
240        continue
        endif
!
        return
!
end subroutine
!DODPC1
subroutine dodpc1                                                             &
         ( ipr, lunrpt,                                                       &
         anajac, cdjac, chkjac, initd, restrt, isodr, implct, dovcv, redoj,   &
         msgb1, msgb, msgd1, msgd,                                            &
         n, m, np, nq, npp, nnzw,                                             &
         x, ldx, ifixx, ldifx, delta, wd, ldwd, ld2wd, tt, ldtt, stpd, ldstpd &
         , y, ldy, we, ldwe, ld2we, pnlty, beta, ifixb, ssf, stpb, lower,     &
         upper, job, neta, taufac, sstol, partol, maxit, wss, wssdel, wsseps)
!***Begin Prologue  DODPC1
!***Refer to  ODR
!***Routines Called  DHSTEP
!***Date Written   860529   (YYMMDD)
!***Revision Date  920619   (YYMMDD)
!***Purpose  Generate initial summary report
!***End Prologue  DODPC1
!
!...Used modules
        use odrpack_kinds,only: wp
!
!...Scalar arguments
        real(kind = wp)                                                       &
         partol, pnlty, sstol, taufac, wss, wssdel, wsseps
        integer                                                               &
         ipr, job, ldifx, ldstpd, ldtt, ldwd, ldwe, ldx, ldy, ld2wd, ld2we,   &
         lunrpt, m, maxit, msgb1, msgd1, n, neta, nnzw, np, npp, nq
        logical                                                               &
         anajac, cdjac, chkjac, dovcv, implct, initd, isodr, redoj, restrt
!
!...Array arguments
        real(kind = wp)                                                       &
         beta( np), delta( n, m), lower( np), ssf( np), stpb( np), stpd(      &
         ldstpd, m), tt( ldtt, m), upper( np), wd( ldwd, ld2wd, m), we( ldwe, &
         ld2we, nq), x( ldx, m), y( ldy, nq)
        integer                                                               &
         ifixb( np), ifixx( ldifx, m), msgb( nq, np), msgd( nq, m)
!
!...Local scalars
        real(kind = wp)                                                       &
         temp1, temp2, temp3, zero
        integer                                                               &
         i, itemp, j, job1, job2, job3, job4, job5, l
!
!...Local arrays
        character tempc0*2, tempc1*5, tempc2*13
!
!...External functions
        real(kind = wp)                                                       &
         dhstep
        external                                                              &
         dhstep
!
!...Data statements
        data                                                                  &
         zero                                                                 &
         /0.0E0_wp/
!
!...Variable Definitions (alphabetically)
!       ANAJAC:  The variable designating whether the Jacobians are computed
!       by finite differences (ANAJAC=FALSE) or not (ANAJAC=TRUE).
!       BETA:    The function parameters.
!       CDJAC:   The variable designating whether the Jacobians are computed
!       by central differences (CDJAC=TRUE) or forward differences
!       (CDJAC=FALSE).
!       CHKJAC:  The variable designating whether the user supplied
!       Jacobians are to be checked (CHKJAC=TRUE) or not
!       (CHKJAC=FALSE).
!       DELTA:   The estimated errors in the explanatory variables.
!       DOVCV:   The variable designating whether the covariance matrix is
!       to be computed (DOVCV=TRUE) or not (DOVCV=FALSE).
!       I:       An indexing variable.
!       IFIXB:   The values designating whether the elements of BETA are
!       fixed at their input values or not.
!       IFIXX:   The values designating whether the elements of X are
!       fixed at their input values or not.
!       IMPLCT:  The variable designating whether the solution is by
!       implicit ODR (IMPLCT=TRUE) or explicit ODR (IMPLCT=FALSE).
!       INITD:   The variable designating whether DELTA is initialized to
!       zero (INITD=TRUE) or to the values in the first N by M
!       elements of array WORK (INITD=FALSE).
!       IPR:     The value indicating the report to be printed.
!       ISODR:   The variable designating whether the solution is by ODR
!       (ISODR=TRUE) or by OLS (ISODR=FALSE).
!       ITEMP:   A temporary integer value.
!       J:       An indexing variable.
!       JOB:     The variable controling problem initialization and
!       computational method.
!       JOB1:    The 1st digit (from the left) of variable JOB.
!       JOB2:    The 2nd digit (from the left) of variable JOB.
!       JOB3:    The 3rd digit (from the left) of variable JOB.
!       JOB4:    The 4th digit (from the left) of variable JOB.
!       JOB5:    The 5th digit (from the left) of variable JOB.
!       L:       An indexing variable.
!       LDIFX:   The leading dimension of array IFIXX.
!       LDTT:    The leading dimension of array TT.
!       LDWD:    The leading dimension of array WD.
!       LDWE:    The leading dimension of array WE.
!       LDX:     The leading dimension of array X.
!       LDY:     The leading dimension of array Y.
!       LD2WD:   The second dimension of array WD.
!       LD2WE:   The second dimension of array WE.
!       LUNRPT:  The logical unit number for the computation reports.
!       M:       The number of columns of data in the explanatory variable.
!       MAXIT:   The maximum number of iterations allowed.
!       MSGB:    The error checking results for the Jacobian wrt beta.
!       MSGB1:   The error checking results for the Jacobian wrt BETA.
!       MSGD:    The error checking results for the Jacobian wrt DELTA.
!       MSGD1:   The error checking results for the Jacobian wrt DELTA.
!       N:       The number of observations.
!       NETA:    The number of accurate digits in the function results.
!       A negative value indicates that NETA was estimated by
!       ODRPACK95. A positive value indictes the value was supplied
!       by the user.
!       NNZW:    The number of nonzero observational error weights.
!       NP:      The number of function parameters.
!       NPP:     The number of function parameters being estimated.
!       NQ:      The number of responses per observation.
!       PARTOL:  The parameter convergence stopping tolerance.
!       PNLTY:   The penalty parameter for an implicit model.
!       REDOJ:   The variable designating whether the Jacobian matrix is to
!       be recomputed for the computation of the covariance matrix
!       (REDOJ=TRUE) or not (REDOJ=FALSE).
!       RESTRT:  The variable designating whether the call is a restart
!       (RESTRT=TRUE) or not (RESTRT=FALSE).
!       SSF:     The scaling values for BETA.
!       SSTOL:   The sum-of-squares convergence stopping tolerance.
!       STPB:    The relative step used for computing finite difference
!       derivatives with respect to BETA.
!       STPD:    The relative step used for computing finite difference
!       derivatives with respect to DELTA.
!       TAUFAC:  The factor used to compute the initial trust region
!       diameter.
!       TEMPC0:  A temporary CHARACTER*2 value.
!       TEMPC1:  A temporary CHARACTER*5 value.
!       TEMPC2:  A temporary CHARACTER*13 value.
!       TEMP1:   A temporary REAL (KIND=wp) value.
!       TEMP2:   A temporary REAL (KIND=wp) value.
!       TEMP3:   A temporary REAL (KIND=wp) value.
!       TT:      The scaling values for DELTA.
!       WD:      The DELTA weights.
!       WE:      The EPSILON weights.
!       WSS:     The sum-of-squares of the weighted EPSILONS and DELTAS.
!       WSSDEL:  The sum-of-squares of the weighted DELTAS.
!       WSSEPS:  The sum-of-squares of the weighted EPSILONS.
!       X:       The explanatory variable.
!       Y:       The response variable.  Unused when the model is implicit.
!       ZERO:    The value 0.0E0_wp.
!
!
!***First executable statement  DODPC1
!
!
!  Print problem size specification
!
        write ( lunrpt,1000) n, nnzw, nq, m, np, npp
!
!
!  Print control values
!
        job1 = job/10000
        job2 = mod( job,10000)/1000
        job3 = mod( job,1000)/100
        job4 = mod( job,100)/10
        job5 = mod( job,10)
        write ( lunrpt,1100) job
        if ( restrt) then
           write ( lunrpt,1110) job1
        else
           write ( lunrpt,1111) job1
        endif
        if ( isodr) then
           if ( initd) then
              write ( lunrpt,1120) job2
           else
              write ( lunrpt,1121) job2
           endif
        else
           write ( lunrpt,1122) job2, job5
        endif
        if ( dovcv) then
           write ( lunrpt,1130) job3
           if ( redoj) then
              write ( lunrpt,1131)
           else
              write ( lunrpt,1132)
           endif
        else
           write ( lunrpt,1133) job3
        endif
        if ( anajac) then
           write ( lunrpt,1140) job4
           if ( chkjac) then
              if ( msgb1 .ge. 1 .or. msgd1 .ge. 1) then
                 write ( lunrpt,1141)
              else
                 write ( lunrpt,1142)
              endif
           else
              write ( lunrpt,1143)
           endif
        elseif ( cdjac) then
           write ( lunrpt,1144) job4
        else
           write ( lunrpt,1145) job4
        endif
        if ( isodr) then
           if ( implct) then
              write ( lunrpt,1150) job5
           else
              write ( lunrpt,1151) job5
           endif
        else
           write ( lunrpt,1152) job5
        endif
        if ( neta .lt. 0) then
           write ( lunrpt,1200)- neta
        else
           write ( lunrpt,1210) neta
        endif
        write ( lunrpt,1300) taufac
!
!
!  Print stopping criteria
!
        write ( lunrpt,1400) sstol, partol, maxit
!
!
!  Print initial sum of squares
!
        if ( implct) then
           write ( lunrpt,1500) wssdel
           if ( isodr) then
              write ( lunrpt,1510) wss, wsseps, pnlty
           endif
        else
           write ( lunrpt,1600) wss
           if ( isodr) then
              write ( lunrpt,1610) wssdel, wsseps
           endif
        endif
!
!
        if ( ipr .ge. 2) then
!
!
!  Print function parameter data
!
           write ( lunrpt,4000)
           if ( chkjac .and.                                                  &
               (( msgb1 .ge. 1) .or.                                          &
                ( msgd1 .ge. 1))) then
              write ( lunrpt,4110)
           elseif ( anajac) then
              write ( lunrpt,4120)
           else
              write ( lunrpt,4200)
           endif
           do 130 j = 1, np
              if ( ifixb(1) .lt. 0) then
                 tempc1 = '   NO'
              else
                 if ( ifixb( j) .ne. 0) then
                    tempc1 = '   NO'
                 else
                    tempc1 = '  YES'
                 endif
              endif
              if ( anajac) then
                 if ( chkjac .and.                                            &
                     (( msgb1 .ge. 1) .or.                                    &
                      ( msgd1 .ge. 1))) then
                    itemp = -1
                    do 110 l = 1, nq
                       itemp = max( itemp, msgb( l, j))
110                 continue
                    if ( itemp .le. -1) then
                       tempc2 = '    UNCHECKED'
                    elseif ( itemp .eq. 0) then
                       tempc2 = '     VERIFIED'
                    elseif ( itemp .ge. 1) then
                       tempc2 = ' QUESTIONABLE'
                    endif
                 else
                    tempc2 = '             '
                 endif
              else
                 tempc2 = '             '
              endif
              if ( ssf(1) .lt. zero) then
                 temp1 = abs( ssf(1))
              else
                 temp1 = ssf( j)
              endif
              if ( anajac) then
                 write ( lunrpt,4310) j, beta( j), tempc1, temp1, lower( j),  &
                  upper( j), tempc2
              else
                 if ( cdjac) then
                    temp2 = dhstep(1, neta,1, j, stpb,1)
                 else
                    temp2 = dhstep(0, neta,1, j, stpb,1)
                 endif
                 write ( lunrpt,4320) j, beta( j), tempc1, temp1,             &
                  lower( j), upper( j), temp2
              endif
130        continue
!
!  Print explanatory variable data
!
           if ( isodr) then
              write ( lunrpt,2010)
              if ( chkjac .and.                                               &
                  (( msgb1 .ge. 1) .or.                                       &
                   ( msgd1 .ge. 1))) then
                 write ( lunrpt,2110)
              elseif ( anajac) then
                 write ( lunrpt,2120)
              else
                 write ( lunrpt,2130)
              endif
           else
              write ( lunrpt,2020)
              write ( lunrpt,2140)
           endif
           if ( isodr) then
              do 240 j = 1, m
                 tempc0 = '1,'
                 do 230 i = 1, n, n-1
!
                    if ( ifixx(1,1) .lt. 0) then
                       tempc1 = '   NO'
                    else
                       if ( ldifx .eq. 1) then
                          if ( ifixx(1, j) .eq. 0) then
                             tempc1 = '  YES'
                          else
                             tempc1 = '   NO'
                          endif
                       else
                          if ( ifixx( i, j) .eq. 0) then
                             tempc1 = '  YES'
                          else
                             tempc1 = '   NO'
                          endif
                       endif
                    endif
!
                    if ( tt(1,1) .lt. zero) then
                       temp1 = abs( tt(1,1))
                    else
                       if ( ldtt .eq. 1) then
                          temp1 = tt(1, j)
                       else
                          temp1 = tt( i, j)
                       endif
                    endif
!
                    if ( wd(1,1,1) .lt. zero) then
                       temp2 = abs( wd(1,1,1))
                    else
                       if ( ldwd .eq. 1) then
                          if ( ld2wd .eq. 1) then
                             temp2 = wd(1,1, j)
                          else
                             temp2 = wd(1, j, j)
                          endif
                       else
                          if ( ld2wd .eq. 1) then
                             temp2 = wd( i,1, j)
                          else
                             temp2 = wd( i, j, j)
                          endif
                       endif
                    endif
!
                    if ( anajac) then
                       if ( chkjac .and.                                      &
                           ((( msgb1 .ge. 1) .or.                             &
                             ( msgd1 .ge. 1)) .and.                           &
                            ( i .eq. 1))) then
                          itemp = -1
                          do 210 l = 1, nq
                             itemp = max( itemp, msgd( l, j))
210                       continue
                          if ( itemp .le. -1) then
                             tempc2 = '    UNCHECKED'
                          elseif ( itemp .eq. 0) then
                             tempc2 = '     VERIFIED'
                          elseif ( itemp .ge. 1) then
                             tempc2 = ' QUESTIONABLE'
                          endif
                       else
                          tempc2 = '             '
                       endif
                       if ( m .le. 9) then
                          write ( lunrpt,5110)                                &
                           tempc0, j, x( i, j),                               &
                           delta( i, j), tempc1, temp1, temp2, tempc2
                       else
                          write ( lunrpt,5120)                                &
                           tempc0, j, x( i, j),                               &
                           delta( i, j), tempc1, temp1, temp2, tempc2
                       endif
                    else
                       tempc2 = '             '
                       if ( cdjac) then
                          temp3 = dhstep(1, neta, i, j, stpd, ldstpd)
                       else
                          temp3 = dhstep(0, neta, i, j, stpd, ldstpd)
                       endif
                       if ( m .le. 9) then
                          write ( lunrpt,5210)                                &
                           tempc0, j, x( i, j),                               &
                           delta( i, j), tempc1, temp1, temp2, temp3
                       else
                          write ( lunrpt,5220)                                &
                           tempc0, j, x( i, j),                               &
                           delta( i, j), tempc1, temp1, temp2, temp3
                       endif
                    endif
!
                    tempc0 = 'N,'
!
230              continue
                 if ( j .lt. m)write ( lunrpt,6000)
240           continue
           else
!
              do 260 j = 1, m
                 tempc0 = '1,'
                 do 250 i = 1, n, n-1
                    if ( m .le. 9) then
                       write ( lunrpt,5110)                                   &
                        tempc0, j, x( i, j)
                    else
                       write ( lunrpt,5120)                                   &
                        tempc0, j, x( i, j)
                    endif
                    tempc0 = 'N,'
250              continue
                 if ( j .lt. m)write ( lunrpt,6000)
260           continue
           endif
!
!  Print response variable data and observation error weights
!
           if (.not. implct) then
              write ( lunrpt,3000)
              write ( lunrpt,3100)
              do 310 l = 1, nq
                 tempc0 = '1,'
                 do 300 i = 1, n, n-1
                    if ( we(1,1,1) .lt. zero) then
                       temp1 = abs( we(1,1,1))
                    elseif ( ldwe .eq. 1) then
                       if ( ld2we .eq. 1) then
                          temp1 = we(1,1, l)
                       else
                          temp1 = we(1, l, l)
                       endif
                    else
                       if ( ld2we .eq. 1) then
                          temp1 = we( i,1, l)
                       else
                          temp1 = we( i, l, l)
                       endif
                    endif
                    if ( nq .le. 9) then
                       write ( lunrpt,5110)                                   &
                        tempc0, l, y( i, l), temp1
                    else
                       write ( lunrpt,5120)                                   &
                        tempc0, l, y( i, l), temp1
                    endif
                    tempc0 = 'N,'
300              continue
                 if ( l .lt. nq)write ( lunrpt,6000)
310           continue
           endif
        endif
!
        return
!
!  Format statements
!
1000    format                                                                &
         (/' --- Problem Size:'/                                              &
         '            N = ',I5,                                               &
         '          (number with nonzero weight = ',I5,')'/                   &
         '           NQ = ',I5/                                               &
         '            M = ',I5/                                               &
         '           NP = ',I5,                                               &
         '          (number unfixed = ',I5,')')
1100    format                                                                &
         (/' --- Control Values:'/                                            &
         '          JOB = ',I5.5/                                             &
         '              = ABCDE, where')
1110    format                                                                &
         ('                       A=',I1,' ==> fit is a restart.')
1111    format                                                                &
         ('                       A=',I1,' ==> fit is not a restart.')
1120    format                                                                &
         ('                       B=',I1,' ==> deltas are initialized',       &
         ' to zero.')
1121    format                                                                &
         ('                       B=',I1,' ==> deltas are initialized',       &
         ' by user.')
1122    format                                                                &
         ('                       B=',I1,' ==> deltas are fixed at',          &
         ' zero since E=',I1,'.')
1130    format                                                                &
         ('                       C=',I1,' ==> covariance matrix will',       &
         ' be computed using')
1131    format                                                                &
         ('                               derivatives re-',                   &
         'evaluated at the solution.')
1132    format                                                                &
         ('                               derivatives from the',              &
         ' last iteration.')
1133    format                                                                &
         ('                       C=',I1,' ==> covariance matrix will',       &
         ' not be computed.')
1140    format                                                                &
         ('                       D=',I1,' ==> derivatives are',              &
         ' supplied by user.')
1141    format                                                                &
         ('                               derivatives were checked.'/         &
         '                               results appear questionable.')
1142    format                                                                &
         ('                               derivatives were checked.'/         &
         '                               results appear correct.')
1143    format                                                                &
         ('                               derivatives were not',              &
         ' checked.')
1144    format                                                                &
         ('                       D=',I1,' ==> derivatives are',              &
         ' estimated by central',                                             &
         ' differences.')
1145    format                                                                &
         ('                       D=',I1,' ==> derivatives are',              &
         ' estimated by forward',                                             &
         ' differences.')
1150    format                                                                &
         ('                       E=',I1,' ==> method is implicit ODR.')
1151    format                                                                &
         ('                       E=',I1,' ==> method is explicit ODR.')
1152    format                                                                &
         ('                       E=',I1,' ==> method is explicit OLS.')
1200    format                                                                &
         ('       NDIGIT = ',I5,'          (estimated by ODRPACK95)')
1210    format                                                                &
         ('       NDIGIT = ',I5,'          (supplied by user)')
1300    format                                                                &
         ('       TAUFAC = ',1P,E12.2)
1400    format                                                                &
         (/' --- Stopping Criteria:'/                                         &
         '        SSTOL = ',1P,E12.2,                                         &
         '   (sum of squares stopping tolerance)'/                            &
         '       PARTOL = ',1P,E12.2,                                         &
         '   (parameter stopping tolerance)'/                                 &
         '        MAXIT = ',I5,                                               &
         '          (maximum number of iterations)')
1500    format                                                                &
         (/' --- Initial Sum of Squared Weighted Deltas =',                   &
         17X,1P,E17.8)
1510    format                                                                &
         ('         Initial Penalty Function Value     =',1P,E17.8/           &
         '                 Penalty Term               =',1P,E17.8/            &
         '                 Penalty Parameter          =',1P,E10.1)
1600    format                                                                &
         (/' --- Initial Weighted Sum of Squares        =',                   &
         17X,1P,E17.8)
1610    format                                                                &
         ('         Sum of Squared Weighted Deltas     =',1P,E17.8/           &
         '         Sum of Squared Weighted Epsilons   =',1P,E17.8)
2010    format                                                                &
         (/' --- Explanatory Variable and Delta Weight Summary:')
2020    format                                                                &
         (/' --- Explanatory Variable Summary:')
2110    format                                                                &
         (/'       Index      X(I,J)  DELTA(I,J)    Fixed',                   &
         '     Scale    Weight    Derivative'/                                &
         '                                             ',                     &
         '                        Assessment'/,                               &
         '       (I,J)                          (IFIXX)',                     &
         '    (SCLD)      (WD)              '/)
2120    format                                                                &
         (/'       Index      X(I,J)  DELTA(I,J)    Fixed',                   &
         '     Scale    Weight              '/                                &
         '                                             ',                     &
         '                                  '/,                               &
         '       (I,J)                          (IFIXX)',                     &
         '    (SCLD)      (WD)              '/)
2130    format                                                                &
         (/'       Index      X(I,J)  DELTA(I,J)    Fixed',                   &
         '     Scale    Weight    Derivative'/                                &
         '                                             ',                     &
         '                         Step Size'/,                               &
         '       (I,J)                          (IFIXX)',                     &
         '    (SCLD)      (WD)        (STPD)'/)
2140    format                                                                &
         (/'       Index      X(I,J)'/                                        &
         '       (I,J)            '/)
3000    format                                                                &
         (/' --- Response Variable and Epsilon Error Weight',                 &
         ' Summary:')
3100    format                                                                &
         (/'       Index      Y(I,L)      Weight'/                            &
         '       (I,L)                    (WE)'/)
4000    format                                                                &
         (/' --- Function Parameter Summary:')
4110    format                                                                &
         (/'       Index   BETA(K)    Fixed     Scale   LOWER(K)',            &
         '   UPPER(K)    Derivative'/                                         &
         '                                                    ',              &
         '               Assessment'/,                                        &
         '         (K)            (IFIXB)    (SCLB)           ',              &
         '                         '/)
4120    format                                                                &
         (/'       Index   BETA(K)    Fixed     Scale   LOWER(K)',            &
         '   UPPER(K)              '/                                         &
         '                                                    ',              &
         '                         '/,                                        &
         '         (K)            (IFIXB)    (SCLB)           ',              &
         '                         '/)
4200    format                                                                &
         (/'       Index   BETA(K)    Fixed     Scale   LOWER(K)',            &
         '   UPPER(K)    Derivative'/                                         &
         '                                                    ',              &
         '                Step Size'/,                                        &
         '         (K)            (IFIXB)    (SCLB)           ',              &
         '                   (STPB)'/)
4310    format                                                                &
         (7X,I5,1P,E10.2,4X,A5,E10.2,E11.2E3,E11.2E3,1X,A13)
4320    format                                                                &
         (7X,I5,1P,E10.2,4X,A5,E10.2,E11.2E3,E11.2E3,1X,E13.5)
5110    format                                                                &
         (9X,A2,I1,1P,2E12.3,4X,A5,2E10.2,1X,A13)
5120    format                                                                &
         (8X,A2,I2,1P,2E12.3,4X,A5,2E10.2,1X,A13)
5210    format                                                                &
         (9X,A2,I1,1P,2E12.3,4X,A5,2E10.2,1X,E13.5)
5220    format                                                                &
         (8X,A2,I2,1P,2E12.3,4X,A5,2E10.2,1X,E13.5)
6000    format                                                                &
         (' ')
end subroutine
!DODPC2
subroutine dodpc2                                                             &
         ( ipr, lunrpt, fstitr, implct, prtpen,                               &
         pnlty,                                                               &
         niter, nfev, wss, actred, prered, alpha, tau, pnorm, np, beta)
!***Begin Prologue  DODPC2
!***Refer to  ODR
!***Routines Called  (NONE)
!***Date Written   860529   (YYMMDD)
!***Revision Date  920304   (YYMMDD)
!***Purpose  Generate iteration reports
!***End Prologue  DODPC2
!
!...Used modules
        use odrpack_kinds,only: wp
!
!...Scalar arguments
        real(kind = wp)                                                       &
         actred, alpha, pnlty, pnorm, prered, tau, wss
        integer                                                               &
         ipr, lunrpt, nfev, niter, np
        logical                                                               &
         fstitr, implct, prtpen
!
!...Array arguments
        real(kind = wp)                                                       &
         beta( np)
!
!...Local scalars
        real(kind = wp)                                                       &
         ratio, zero
        integer                                                               &
         j, k, l
        character gn*3
!
!...Data statements
        data                                                                  &
         zero                                                                 &
         /0.0E0_wp/
!
!...Variable Definitions (alphabetically)
!       ACTRED:  The actual relative reduction in the sum-of-squares.
!       ALPHA:   The Levenberg-Marquardt parameter.
!       BETA:    The function parameters.
!       FSTITR:  The variable designating whether this is the first
!       iteration (FSTITR=.TRUE.) or not (FSTITR=.FALSE.).
!       GN:      The CHARACTER*3 variable indicating whether a Gauss-Newton
!       step was taken.
!       IMPLCT:  The variable designating whether the solution is by
!       implicit ODR (IMPLCT=TRUE) or explicit ODR (IMPLCT=FALSE).
!       IPR:     The value indicating the report to be printed.
!       J:       An indexing variable.
!       K:       An indexing variable.
!       L:       An indexing variable.
!       LUNRPT:  The logical unit number used for computation reports.
!       NFEV:    The number of function evaluations.
!       NITER:   The number of iterations.
!       NP:      The number of function parameters.
!       PNLTY:   The penalty parameter for an implicit model.
!       PNORM:   The norm of the scaled estimated parameters.
!       PRERED:  The predicted relative reduction in the sum-of-squares.
!       PRTPEN:  The variable designating whether the penalty parameter is
!       to be printed in the iteration report (PRTPEN=TRUE) or not
!       (PRTPEN=FALSE).
!       RATIO:   The ratio of TAU to PNORM.
!       TAU:     The trust region diameter.
!       WSS:     The sum-of-squares of the weighted EPSILONS and DELTAS.
!       ZERO:    The value 0.0E0_wp.
!
!
!***First executable statement  DODPC2
!
!
        if ( fstitr) then
           if ( ipr .eq. 1) then
              if ( implct) then
                 write ( lunrpt,1121)
              else
                 write ( lunrpt,1122)
              endif
           else
              if ( implct) then
                 write ( lunrpt,1131)
              else
                 write ( lunrpt,1132)
              endif
           endif
        endif
        if ( prtpen) then
           write ( lunrpt,1133) pnlty
        endif
!
        if ( alpha .eq. zero) then
!-----------------------^------------------------------------------------------
!!! FPT - 3087 REAL or COMPLEX quantity tested for exact equality/inequality
!------------------------------------------------------------------------------
           gn = 'YES'
        else
           gn = ' NO'
        endif
        if ( pnorm .ne. zero) then
!-----------------------^------------------------------------------------------
!!! FPT - 3087 REAL or COMPLEX quantity tested for exact equality/inequality
!------------------------------------------------------------------------------
           ratio = tau/ pnorm
        else
           ratio = zero
        endif
        if ( ipr .eq. 1) then
           write ( lunrpt,1141) niter, nfev, wss, actred, prered,             &
            ratio, gn
        else
           j = 1
           k = min(3, np)
           if ( j .eq. k) then
              write ( lunrpt,1141) niter, nfev, wss, actred, prered,          &
               ratio, gn, j, beta( j)
           else
              write ( lunrpt,1142) niter, nfev, wss, actred, prered,          &
               ratio, gn, j, k,( beta( l), l = j, k)
           endif
           if ( np .gt. 3) then
              do 10 j = 4, np,3
                 k = min( j+2, np)
                 if ( j .eq. k) then
                    write ( lunrpt,1151) j, beta( j)
                 else
                    write ( lunrpt,1152) j, k,( beta( l), l = j, k)
                 endif
10            continue
           endif
        endif
!
        return
!
!  Format statements
!
1121    format                                                                &
         (//                                                                  &
         '         Cum.      Penalty    Act. Rel.   Pred. Rel.'/              &
         '  It.  No. FN     Function   Sum-of-Sqs   Sum-of-Sqs',              &
         '              G-N'/                                                 &
         ' Num.   Evals        Value    Reduction    Reduction',              &
         '  TAU/PNORM  Step'/                                                 &
         ' ----  ------  -----------  -----------  -----------',              &
         '  ---------  ----')
1122    format                                                                &
         (//                                                                  &
         '         Cum.                 Act. Rel.   Pred. Rel.'/              &
         '  It.  No. FN     Weighted   Sum-of-Sqs   Sum-of-Sqs',              &
         '              G-N'/                                                 &
         ' Num.   Evals   Sum-of-Sqs    Reduction    Reduction',              &
         '  TAU/PNORM  Step'/                                                 &
         ' ----  ------  -----------  -----------  -----------',              &
         '  ---------  ----'/)
1131    format                                                                &
         (//                                                                  &
         '         Cum.      Penalty    Act. Rel.   Pred. Rel.'/              &
         '  It.  No. FN     Function   Sum-of-Sqs   Sum-of-Sqs',              &
         '              G-N      BETA -------------->'/                       &
         ' Num.   Evals        Value    Reduction    Reduction',              &
         '  TAU/PNORM  Step     Index           Value'/                       &
         ' ----  ------  -----------  -----------  -----------',              &
         '  ---------  ----     -----           -----')
1132    format                                                                &
         (//                                                                  &
         '         Cum.                 Act. Rel.   Pred. Rel.'/              &
         '  It.  No. FN     Weighted   Sum-of-Sqs   Sum-of-Sqs',              &
         '              G-N      BETA -------------->'/                       &
         ' Num.   Evals   Sum-of-Sqs    Reduction    Reduction',              &
         '  TAU/PNORM  Step     Index           Value'/                       &
         ' ----  ------  -----------  -----------  -----------',              &
         '  ---------  ----     -----           -----'/)
1133    format                                                                &
         (/' Penalty Parameter Value = ',1P,E10.1)
1141    format                                                                &
         (1X,I4,I8,1X,1P,E12.5,2E13.4,E11.3,3X,A3,7X,I3,3E16.8)
1142    format                                                                &
         (1X,I4,I8,1X,1P,E12.5,2E13.4,E11.3,3X,A3,1X,I3,' To',I3,3E16.8)
1151    format                                                                &
         (76X,I3,1P,E16.8)
1152    format                                                                &
         (70X,I3,' To',I3,1P,3E16.8)
end subroutine
!DODPC3
subroutine dodpc3                                                             &
         ( ipr, lunrpt,                                                       &
         isodr, implct, didvcv, dovcv, redoj, anajac,                         &
         n, m, np, nq, npp,                                                   &
         info, niter, nfev, njev, irank, rcond, istop,                        &
         wss, wssdel, wsseps, pnlty, rvar, idf,                               &
         beta, sdbeta, ifixb2, f, delta,                                      &
         lower, upper)
!***Begin Prologue  DODPC3
!***Refer to  ODR
!***Routines Called  DPPT
!***Date Written   860529   (YYMMDD)
!***REvision Date  920619   (YYMMDD)
!***Purpose  Generate final summary report
!***End Prologue  DODPC3
!
!...Used modules
        use odrpack_kinds,only: wp
!
!...Scalar arguments
        real(kind = wp)                                                       &
         pnlty, rcond, rvar, wss, wssdel, wsseps
        integer                                                               &
         idf, info, ipr, irank, istop, lunrpt, m,                             &
         n, nfev, niter, njev, np, npp, nq
        logical                                                               &
         anajac, didvcv, dovcv, implct, isodr, redoj
!
!...Array arguments
        real(kind = wp)                                                       &
         beta( np), delta( n, m), f( n, nq), lower( np), upper( np), sdbeta(  &
         np)
        integer                                                               &
         ifixb2( np)
!
!...Local scalars
        real(kind = wp)                                                       &
         tval
        integer                                                               &
         d1, d2, d3, d4, d5, i, j, k, l, nplm1
        character fmt1*90
!
!...External functions
        real(kind = wp)                                                       &
         dppt
        external                                                              &
         dppt
!
!...Variable Definitions (alphabetically)
!       ANAJAC:  The variable designating whether the JACOBIANS are computed
!       by finite differences (ANAJAC=FALSE) or not (ANAJAC=TRUE).
!       BETA:    The function parameters.
!       D1:      The first digit of INFO.
!       D2:      The second digit of INFO.
!       D3:      The third digit of INFO.
!       D4:      The fourth digit of INFO.
!       D5:      The fifth digit of INFO.
!       DELTA:   The estimated errors in the explanatory variables.
!       DIDVCV:  The variable designating whether the covariance matrix was
!       computed (DIDVCV=TRUE) or not (DIDVCV=FALSE).
!       DOVCV:   The variable designating whether the covariance matrix was
!       to be computed (DOVCV=TRUE) or not (DOVCV=FALSE).
!       F:       The estimated values of EPSILON.
!       FMT1:    A CHARACTER*90 variable used for formats.
!       I:       An indexing variable.
!       IDF:     The degrees of freedom of the fit, equal to the number of
!       observations with nonzero weighted derivatives minus the
!       number of parameters being estimated.
!       IFIXB2:  The values designating whether the elements of BETA were
!       estimated, fixed, or dropped because they caused rank
!       deficiency, corresponding to values of IFIXB2 equaling 1,
!       0, and -1, respectively.  If IFIXB2 is -2, then no attempt
!       was made to estimate the parameters because MAXIT = 0.
!       IMPLCT:  The variable designating whether the solution is by
!       implicit ODR (IMPLCT=TRUE) or explicit ODR (IMPLCT=FALSE).
!       INFO:    The variable designating why the computations were stopped.
!       IPR:     The variable indicating what is to be printed.
!       IRANK:   The rank deficiency of the Jacobian wrt BETA.
!       ISODR:   The variable designating whether the solution is by ODR
!       (ISODR=TRUE) or by OLS (ISODR=FALSE).
!       ISTOP:   The variable designating whether there are problems
!       computing the function at the current BETA and DELTA.
!       J:       An indexing variable.
!       K:       An indexing variable.
!       L:       An indexing variable.
!       LOWER:   Lower bound on BETA.
!       LUNRPT:  The logical unit number used for computation reports.
!       M:       The number of columns of data in the explanatory variable.
!       N:       The number of observations.
!       NFEV:    The number of function evaluations.
!       NITER:   The number of iterations.
!       NJEV:    The number of Jacobian evaluations.
!       NP:      The number of function parameters.
!       NPLM1:   The number of items to be printed per line, minus one.
!       NPP:     The number of function parameters being estimated.
!       NQ:      The number of responses per observation.
!       PNLTY:   The penalty parameter for an implicit model.
!       RCOND:   The approximate reciprocal condition of TFJACB.
!       REDOJ:   The variable designating whether the Jacobian matrix is
!       to be recomputed for the computation of the covariance
!       matrix (REDOJ=TRUE) or not (REDOJ=FALSE).
!       RVAR:    The residual variance.
!       SDBETA:  The standard errors of the estimated parameters.
!       TVAL:    The value of the 97.5 percent point function for the
!       T distribution.
!       UPPER:   Upper bound on BETA.
!       WSS:     The sum-of-squares of the weighted EPSILONS and DELTAS.
!       WSSDEL:  The sum-of-squares of the weighted DELTAS.
!       WSSEPS:  The sum-of-squares of the weighted EPSILONS.
!
!
!***First executable statement  DODPC3
!
!
        d1 = info/10000
        d2 = mod( info,10000)/1000
        d3 = mod( info,1000)/100
        d4 = mod( info,100)/10
        d5 = mod( info,10)
!
!  Print stopping conditions
!
        write ( lunrpt,1000)
        if ( info .le. 9) then
           if ( info .eq. 1) then
              write ( lunrpt,1011) info
           elseif ( info .eq. 2) then
              write ( lunrpt,1012) info
           elseif ( info .eq. 3) then
              write ( lunrpt,1013) info
           elseif ( info .eq. 4) then
              write ( lunrpt,1014) info
           elseif ( info .le. 9) then
              write ( lunrpt,1015) info
           endif
        elseif ( info .le. 9999) then
!
!  Print warning diagnostics
!
           write ( lunrpt,1020) info
           if ( d2 .eq. 1)write ( lunrpt,1021)
           if ( d3 .eq. 1)write ( lunrpt,1022)
           if ( d4 .eq. 1)write ( lunrpt,1023)
           if ( d4 .eq. 2)write ( lunrpt,1024)
           if ( d5 .eq. 1) then
              write ( lunrpt,1031)
           elseif ( d5 .eq. 2) then
              write ( lunrpt,1032)
           elseif ( d5 .eq. 3) then
              write ( lunrpt,1033)
           elseif ( d5 .eq. 4) then
              write ( lunrpt,1034)
           elseif ( d5 .le. 9) then
              write ( lunrpt,1035) d5
           endif
        else
!
!  Print error messages
!
           write ( lunrpt,1040) info
           if ( d1 .eq. 5) then
              write ( lunrpt,1042)
              if ( d2 .ne. 0)write ( lunrpt,1043) d2
              if ( d3 .eq. 3) then
                 write ( lunrpt,1044) d3
              elseif ( d3 .ne. 0) then
                 write ( lunrpt,1045) d3
              endif
           elseif ( d1 .eq. 6) then
              write ( lunrpt,1050)
           else
              write ( lunrpt,1060) d1
           endif
        endif
!
!  Print misc. stopping info
!
        write ( lunrpt,1300) niter
        write ( lunrpt,1310) nfev
        if ( anajac)write ( lunrpt,1320) njev
        write ( lunrpt,1330) irank
        write ( lunrpt,1340) rcond
        write ( lunrpt,1350) istop
!
!  Print final sum of squares
!
        if ( implct) then
           write ( lunrpt,2000) wssdel
           if ( isodr) then
              write ( lunrpt,2010) wss, wsseps, pnlty
           endif
        else
           write ( lunrpt,2100) wss
           if ( isodr) then
              write ( lunrpt,2110) wssdel, wsseps
           endif
        endif
        if ( didvcv) then
           write ( lunrpt,2200)sqrt( rvar), idf
        endif
!
        nplm1 = 3
!
!  Print estimated BETA's, and,
!  if, full rank, their standard errors
!
        write ( lunrpt,3000)
        if ( didvcv) then
           write ( lunrpt,7300)
           tval = dppt(0.975E0_wp, idf)
           do 10 j = 1, np
              if ( ifixb2( j) .ge. 1) then
                 write ( lunrpt,8400) j, beta( j),                            &
                  lower( j), upper( j),                                       &
                  sdbeta( j),                                                 &
                  beta( j)- tval* sdbeta( j),                                 &
                  beta( j)+ tval* sdbeta( j)
              elseif ( ifixb2( j) .eq. 0) then
                 write ( lunrpt,8600) j, beta( j), lower( j), upper( j)
              else
                 write ( lunrpt,8700) j, beta( j), lower( j), upper( j)
              endif
10         continue
           if (.not. redoj)write ( lunrpt,7310)
        else
           if ( dovcv) then
              if ( d1 .le. 5) then
                 write ( lunrpt,7410)
              else
                 write ( lunrpt,7420)
              endif
           endif
!
           if (( irank .eq. 0 .and. npp .eq. np) .or. niter .eq. 0) then
              if ( np .eq. 1) then
                 write ( lunrpt,7100)
              else
                 write ( lunrpt,7200)
              endif
              do 20 j = 1, np, nplm1+1
                 k = min( j+ nplm1, np)
                 if ( k .eq. j) then
                    write ( lunrpt,8100) j, beta( j)
                 else
                    write ( lunrpt,8200) j, k,( beta( l), l = j, k)
                 endif
20            continue
              if ( niter .ge. 1) then
                 write ( lunrpt,8800)
              else
                 write ( lunrpt,8900)
              endif
           else
              write ( lunrpt,7500)
              do 30 j = 1, np
                 if ( ifixb2( j) .ge. 1) then
                    write ( lunrpt,8500) j, beta( j), lower( j), upper( j)
                 elseif ( ifixb2( j) .eq. 0) then
                    write ( lunrpt,8600) j, beta( j), lower( j), upper( j)
                 else
                    write ( lunrpt,8700) j, beta( j), lower( j), upper( j)
                 endif
30            continue
           endif
        endif
!
        if ( ipr .eq. 1) return
!
!
!  Print EPSILON's and DELTA's together in a column if the number of
!  columns of data in EPSILON and DELTA is less than or equal to three.
!
        if ( implct .and.                                                     &
            ( m .le. 4)) then
           write ( lunrpt,4100)
           write ( fmt1,9110) m
           write ( lunrpt, fmt1)( j, j = 1, m)
           do 40 i = 1, n
              write ( lunrpt,4130) i,( delta( i, j), j = 1, m)
40         continue
!
        elseif ( isodr .and.                                                  &
                ( nq+ m .le. 4)) then
           write ( lunrpt,4110)
           write ( fmt1,9120) nq, m
           write ( lunrpt, fmt1)( l, l = 1, nq),( j, j = 1, m)
           do 50 i = 1, n
              write ( lunrpt,4130) i,( f( i, l), l = 1, nq),( delta( i, j), j &
               =1, m)
50         continue
!
        elseif (.not. isodr .and.                                             &
                (( nq .ge. 2) .and.                                           &
                 ( nq .le. 4))) then
           write ( lunrpt,4120)
           write ( fmt1,9130) nq
           write ( lunrpt, fmt1)( l, l = 1, nq)
           do 60 i = 1, n
              write ( lunrpt,4130) i,( f( i, l), l = 1, nq)
60         continue
        else
!
!  Print EPSILON's and DELTA's separately
!
           if (.not. implct) then
!
!  Print EPSILON'S
!
              do 80 j = 1, nq
                 write ( lunrpt,4200) j
                 if ( n .eq. 1) then
                    write ( lunrpt,7100)
                 else
                    write ( lunrpt,7200)
                 endif
                 do 70 i = 1, n, nplm1+1
                    k = min( i+ nplm1, n)
                    if ( i .eq. k) then
                       write ( lunrpt,8100) i, f( i, j)
                    else
                       write ( lunrpt,8200) i, k,( f( l, j), l = i, k)
                    endif
70               continue
80            continue
           endif
!
!  Print DELTA'S
!
           if ( isodr) then
              do 100 j = 1, m
                 write ( lunrpt,4300) j
                 if ( n .eq. 1) then
                    write ( lunrpt,7100)
                 else
                    write ( lunrpt,7200)
                 endif
                 do 90 i = 1, n, nplm1+1
                    k = min( i+ nplm1, n)
                    if ( i .eq. k) then
                       write ( lunrpt,8100) i, delta( i, j)
                    else
                       write ( lunrpt,8200) i, k,( delta( l, j), l = i, k)
                    endif
90               continue
100           continue
           endif
        endif
!
        return
!
!  Format statements
!
1000    format                                                                &
         (/' --- Stopping Conditions:')
1011    format                                                                &
         ('         INFO = ',I5,' ==> sum of squares convergence.')
1012    format                                                                &
         ('         INFO = ',I5,' ==> parameter convergence.')
1013    format                                                                &
         ('         INFO = ',I5,' ==> sum of squares convergence and',        &
         ' parameter convergence.')
1014    format                                                                &
         ('         INFO = ',I5,' ==> iteration limit reached.')
1015    format                                                                &
         ('         INFO = ',I5,' ==> unexpected value,',                     &
         ' probably indicating'/                                              &
         '                           incorrectly specified',                  &
         ' user input.')
1020    format                                                                &
         ('         INFO = ',I5.4/                                            &
         '              =  ABCD, where a nonzero value for digit A,',         &
         ' B, or C indicates why'/                                            &
         '                       the results might be questionable,',         &
         ' and digit D indicates'/                                            &
         '                       the actual stopping condition.')
1021    format                                                                &
         ('                       A=1 ==> derivatives are',                   &
         ' questionable.')
1022    format                                                                &
         ('                       B=1 ==> user set ISTOP to',                 &
         ' nonzero value during last'/                                        &
         '                               call to subroutine FCN.')
1023    format                                                                &
         ('                       C=1 ==> derivatives are not',               &
         ' full rank at the solution.')
1024    format                                                                &
         ('                       C=2 ==> derivatives are zero',              &
         ' rank at the solution.')
1031    format                                                                &
         ('                       D=1 ==> sum of squares convergence.')
1032    format                                                                &
         ('                       D=2 ==> parameter convergence.')
1033    format                                                                &
         ('                       D=3 ==> sum of squares convergence',        &
         ' and parameter convergence.')
1034    format                                                                &
         ('                       D=4 ==> iteration limit reached.')
1035    format                                                                &
         ('                       D=',I1,' ==> unexpected value,',            &
         ' probably indicating'/                                              &
         '                               incorrectly specified',              &
         ' user input.')
1040    format                                                                &
         ('         INFO = ',I5.5/                                            &
         '              = ABCDE, where a nonzero value for a given',          &
         ' digit indicates an'/                                               &
         '                       abnormal stopping condition.')
1042    format                                                                &
         ('                       A=5 ==> user stopped computations',         &
         ' in subroutine FCN.')
1043    format                                                                &
         ('                       B=',I1,' ==> computations were',            &
         ' stopped during the'/                                               &
         '                                    function evaluation.')
1044    format                                                                &
         ('                       C=',I1,' ==> computations were',            &
         ' stopped because'/                                                  &
         '                                    derivatives with',              &
         ' respect to delta were'/                                            &
         '                                    computed by',                   &
         ' subroutine FCN when'/                                              &
         '                                    fit is OLS.')
1045    format                                                                &
         ('                       C=',I1,' ==> computations were',            &
         ' stopped during the'/                                               &
         '                                    jacobian evaluation.')
1050    format                                                                &
         ('                       A=6 ==> numerical instabilities',           &
         ' have been detected,'/                                              &
         '                               possibly indicating',                &
         ' a discontinuity in the'/                                           &
         '                               derivatives or a poor',              &
         ' poor choice of problem'/                                           &
         '                               scale or weights.')
1060    format                                                                &
         ('                       A=',I1,' ==> unexpected value,',            &
         ' probably indicating'/                                              &
         '                               incorrectly specified',              &
         ' user input.')
1300    format                                                                &
         ('        NITER = ',I5,                                              &
         '          (number of iterations)')
1310    format                                                                &
         ('         NFEV = ',I5,                                              &
         '          (number of function evaluations)')
1320    format                                                                &
         ('         NJEV = ',I5,                                              &
         '          (number of jacobian evaluations)')
1330    format                                                                &
         ('        IRANK = ',I5,                                              &
         '          (rank deficiency)')
1340    format                                                                &
         ('        RCOND = ',1P,E12.2,                                        &
         '   (inverse condition number)')
!1341 FORMAT
!       +  ('                      ==> POSSIBLY FEWER THAN 2 SIGNIFICANT',
!       +                        ' DIGITS IN RESULTS;'/
!       +   '                          SEE ODRPACK95 REFERENCE',
!       +                        ' GUIDE, SECTION 4.C.')
1350    format                                                                &
         ('        ISTOP = ',I5,                                              &
         '          (returned by user from',                                  &
         ' subroutine FCN)')
2000    format                                                                &
         (/' --- Final Sum of Squared Weighted Deltas = ',                    &
         17X,1P,E17.8)
2010    format                                                                &
         ('         Final Penalty Function Value     = ',1P,E17.8/            &
         '               Penalty Term               = ',1P,E17.8/             &
         '               Penalty Parameter          = ',1P,E10.1)
2100    format                                                                &
         (/' --- Final Weighted Sums of Squares       = ',17X,1P,E17.8)
2110    format                                                                &
         ('         Sum of Squared Weighted Deltas   = ',1P,E17.8/            &
         '         Sum of Squared Weighted Epsilons = ',1P,E17.8)
2200    format                                                                &
         (/' --- Residual Standard Deviation          = ',                    &
         17X,1P,E17.8/                                                        &
         '         Degrees of Freedom               =',I5)
3000    format                                                                &
         (/' --- Estimated BETA(J), J = 1, ..., NP:')
4100    format                                                                &
         (/' --- Estimated DELTA(I,*), I = 1, ..., N:')
4110    format                                                                &
         (/' --- Estimated EPSILON(I) and DELTA(I,*), I = 1, ..., N:')
4120    format                                                                &
         (/' --- Estimated EPSILON(I), I = 1, ..., N:')
4130    format (5X,I5,1P,5E16.8)
4200    format                                                                &
         (/' --- Estimated EPSILON(I,',I3,'), I = 1, ..., N:')
4300    format                                                                &
         (/' --- Estimated DELTA(I,',I3,'), I = 1, ..., N:')
7100    format                                                                &
         (/'           Index           Value'/)
7200    format                                                                &
         (/'           Index           Value -------------->'/)
7300    format                                                                &
         (/'                     BETA      LOWER     UPPER      S.D. ',       &
         ' ___ 95% Confidence ___'/                                           &
         '                                                    BETA ',         &
         '        Interval'/)
7310    format                                                                &
         (/'     N.B. standard errors and confidence intervals are',          &
         ' computed using'/                                                   &
         '          derivatives calculated at the beginning',                 &
         ' of the last iteration,'/                                           &
         '          and not using derivatives re-evaluated at the',           &
         ' final solution.')
7410    format                                                                &
         (/'     N.B. the standard errors of the estimated betas were',       &
         ' not computed because'/                                             &
         '          the derivatives were not available.  Either MAXIT',       &
         ' is 0 and the third'/                                               &
         '          digit of JOB is greater than 1, or the most',             &
         ' recently tried values of'/                                         &
         '          BETA and/or X+DELTA were identified as',                  &
         ' unacceptable by user supplied'/                                    &
         '          subroutine FCN.')
7420    format                                                                &
         (/'     N.B. the standard errors of the estimated betas were',       &
         ' not computed.'/                                                    &
         '          (see info above.)')
7500    format                                                                &
         (/'                     BETA         Status')
8100    format                                                                &
         (11X,I5,1P,E16.8)
8200    format                                                                &
         (3X,I5,' to',I5,1P,7E16.8)
8400    format                                                                &
         (3X,I5,1X,1P,E16.8,1X,E10.2,E10.2,E10.2,1X,E10.2,1X,'to',E10.2)
8500    format                                                                &
         (3X,I5,1X,1P,E16.8,1X,E10.2,E10.2,4X,'Estimated')
8600    format                                                                &
         (3X,I5,1X,1P,E16.8,1X,E10.2,E10.2,4X,'    Fixed')
8700    format                                                                &
         (3X,I5,1X,1P,E16.8,1X,E10.2,E10.2,4X,'  Dropped')
8800    format                                                                &
         (/'     N.B. no parameters were fixed by the user or',               &
         ' dropped at the last'/                                              &
         '          iteration because they caused the model to be',           &
         ' rank deficient.')
8900    format                                                                &
         (/'     N.B. no change was made to the user supplied parameter',     &
         ' values because'/                                                   &
         '          MAXIT=0.')
9110    format                                                                &
         ('(/''         I'',',                                                &
         I2,'(''      DELTA(I,'',I1,'')'')/)')
9120    format                                                                &
         ('(/''         I'',',                                                &
         I2,'(''    EPSILON(I,'',I1,'')''),',                                 &
         I2,'(''      DELTA(I,'',I1,'')'')/)')
9130    format                                                                &
         ('(/''         I'',',                                                &
         I2,'(''    EPSILON(I,'',I1,'')'')/)')
!
end subroutine
!DODPCR
subroutine dodpcr                                                             &
         ( ipr, lunrpt,                                                       &
         head, prtpen, fstitr, didvcv, iflag,                                 &
         n, m, np, nq, npp, nnzw,                                             &
         msgb, msgd, beta, y, ldy, x, ldx, delta,                             &
         we, ldwe, ld2we, wd, ldwd, ld2wd,                                    &
         ifixb, ifixx, ldifx,                                                 &
         lower, upper,                                                        &
         ssf, tt, ldtt, stpb, stpd, ldstpd,                                   &
         job, neta, taufac, sstol, partol, maxit,                             &
         wss, rvar, idf, sdbeta,                                              &
         niter, nfev, njev, actred, prered,                                   &
         tau, pnorm, alpha, f, rcond, irank, info, istop)
!***Begin Prologue  DODPCR
!***Refer to  ODR
!***Routines Called  DFLAGS,DODPC1,DODPC2,DODPC3,DODPHD
!***Date Written   860529   (YYMMDD)
!***Revision Date  920619   (YYMMDD)
!***Purpose  Generate computation reports
!***End Prologue  DODPCR
!
!...Used modules
        use odrpack_kinds,only: wp
!
!...Scalar arguments
        real(kind = wp)                                                       &
         actred, alpha, partol, pnorm, prered, rcond, rvar,                   &
         sstol, tau, taufac
        integer                                                               &
         idf, iflag, info, ipr, irank, istop, job, ldifx, ldstpd, ldtt, ldwd, &
         ldwe, ldx, ldy, ld2wd, ld2we, lunrpt, m, maxit, n, neta, nfev, niter &
         , njev, nnzw, np, npp, nq
        logical                                                               &
         didvcv, fstitr, head, prtpen
!
!...Array arguments
        real(kind = wp)                                                       &
         beta( np), delta( n, m), f( n, nq), lower( np), sdbeta( np), ssf( np &
         ), stpb( np), stpd( ldstpd, m), tt( ldtt, m), upper( np), wd( ldwd,  &
         ld2wd, m), we( ldwe, ld2we, nq), wss(3), x( ldx, m), y( ldy, nq)
        integer                                                               &
         ifixb( np), ifixx( ldifx, m), msgb( nq* np+1), msgd( nq* m+1)
!
!...Local scalars
        real(kind = wp)                                                       &
         pnlty
        logical                                                               &
         anajac, cdjac, chkjac, dovcv, implct, initd, isodr, redoj, restrt
        character typ*3
!
!...External subroutines
        external                                                              &
         dflags, dodpc1, dodpc2, dodpc3, dodphd
!
!...Variable Definitions (alphabetically)
!       ACTRED:  The actual relative reduction in the sum-of-squares.
!       ALPHA:   The Levenberg-Marquardt parameter.
!       ANAJAC:  The variable designating whether the Jacobians are computed
!       by finite differences (ANAJAC=FALSE) or not (ANAJAC=TRUE).
!       BETA:    The function parameters.
!       CDJAC:   The variable designating whether the jacobians are computed
!       by central differences (CDJAC=TRUE) or by forward
!       differences (CDJAC=FALSE).
!       CHKJAC:  The variable designating whether the user supplied
!       Jacobians are to be checked (CHKJAC=TRUE) or not
!       (CHKJAC=FALSE).
!       DELTA:   The estimated errors in the explanatory variables.
!       DIDVCV:  The variable designating whether the covariance matrix was
!       computed (DIDVCV=TRUE) or not (DIDVCV=FALSE).
!       DOVCV:   The variable designating whether the covariance matrix is
!       to be computed (DOVCV=TRUE) or not (DOVCV=FALSE).
!       F:       The (weighted) estimated values of EPSILON.
!       FSTITR:  The variable designating whether this is the first
!       iteration (FSTITR=TRUE) or not (FSTITR=FALSE).
!       HEAD:    The variable designating whether the heading is to be
!       printed (HEAD=TRUE) or not (HEAD=FALSE).
!       IDF:     The degrees of freedom of the fit, equal to the number of
!       observations with nonzero weighted derivatives minus the
!       number of parameters being estimated.
!       IFIXB:   The values designating whether the elements of BETA are
!       fixed at their input values or not.
!       IFIXX:   The values designating whether the elements of X are
!       fixed at their input values or not.
!       IFLAG:   The variable designating what is to be printed.
!       IMPLCT:  The variable designating whether the solution is by
!       implicit ODR (IMPLCT=TRUE) or explicit ODR (IMPLCT=FALSE).
!       INFO:    The variable designating why the computations were stopped.
!       INITD:   The variable designating whether DELTA is initialized to
!       zero (INITD=TRUE) or to the values in the first N  by M
!       elements of array WORK (INITD=FALSE).
!       IPR:     The value indicating the report to be printed.
!       IRANK:   The rank deficiency of the Jacobian wrt BETA.
!       ISODR:   The variable designating whether the solution is by ODR
!       (ISODR=TRUE) or by OLS (ISODR=FALSE).
!       ISTOP:   The variable designating whether there are problems
!       computing the function at the current BETA and DELTA.
!       JOB:     The variable controling problem initialization and
!       computational method.
!       LDIFX:   The leading dimension of array IFIXX.
!       LDSTPD:  The leading dimension of array STPD.
!       LDTT:    The leading dimension of array TT.
!       LDWD:    The leading dimension of array WD.
!       LDWE:    The leading dimension of array WE.
!       LDX:     The leading dimension of array X.
!       LDY:     The leading dimension of array Y.
!       LD2WD:   The second dimension of array WD.
!       LD2WE:   The second dimension of array WE.
!       LOWER:   Lower bound on BETA.
!       LUNRPT:  The logical unit number for computation reports.
!       M:       The number of columns of data in the explanatory variable.
!       MAXIT:   The maximum number of iterations allowed.
!       MSGB:    The error checking results for the Jacobian wrt BETA.
!       MSGD:    The error checking results for the Jacobian wrt DELTA.
!       N:       The number of observations.
!       NETA:    The number of accurate digits in the function results.
!       NFEV:    The number of function evaluations.
!       NITER:   The number of iterations.
!       NJEV:    The number of Jacobian evaluations.
!       NNZW:    The number of nonzero weighted observations.
!       NP:      The number of function parameters.
!       NQ:      The number of responses per observation.
!       NPP:     The number of function parameters being estimated.
!       PARTOL:  The parameter convergence stopping tolerance.
!       PNLTY:   The penalty parameter for an implicit model.
!       PNORM:   The norm of the scaled estimated parameters.
!       PRERED:  The predicted relative reduction in the sum-of-squares.
!       PRTPEN:  The variable designating whether the penalty parameter is
!       to be printed in the iteration report (PRTPEN=TRUE) or not
!       (PRTPEN=FALSE).
!       RCOND:   The approximate reciprocal condition number of TFJACB.
!       REDOJ:   The variable designating whether the Jacobian matrix is to
!       be recomputed for the computation of the covariance matrix
!       (REDOJ=TRUE) or not (REDOJ=FALSE).
!       RESTRT:  The variable designating whether the call is a restart
!       (RESTRT=TRUE) OR NOT (RESTRT=FALSE).
!       RVAR:    The residual variance.
!       SDBETA:  The standard deviations of the estimated BETA'S.
!       SSF:     The scaling values for BETA.
!       SSTOL:   The sum-of-squares convergence stopping tolerance.
!       STPB:    The relative step for computing finite difference
!       derivatives with respect to BETA.
!       STPD:    The relative step for computing finite difference
!       derivatives with respect to DELTA.
!       TAU:     The trust region diameter.
!       TAUFAC:  The factor used to compute the initial trust region
!       diameter.
!       TT:      The scaling values for DELTA.
!       TYP:     The CHARACTER*3 string "ODR" or "OLS".
!       UPPER:   Upper bound on BETA.
!       WE:      The EPSILON weights.
!       WD:      The DELTA weights.
!       WSS:     The sum-of-squares of the weighted EPSILONS and DELTAS,
!       the sum-of-squares of the weighted DELTAS, and
!       the sum-of-squares of the weighted EPSILONS.
!       X:       The explanatory variable.
!       Y:       The dependent variable.  Unused when the model is implicit.
!
!
!***First executable statement  DODPCR
!
!
        call dflags( job, restrt, initd, dovcv, redoj,                        &
         anajac, cdjac, chkjac, isodr, implct)
        pnlty = abs( we(1,1,1))
!
        if ( head) then
           call dodphd( head, lunrpt)
        endif
        if ( isodr) then
           typ = 'ODR'
        else
           typ = 'OLS'
        endif
!
!  Print initial summary
!
        if ( iflag .eq. 1) then
           write ( lunrpt,1200) typ
           call dodpc1                                                        &
            ( ipr, lunrpt,                                                    &
            anajac, cdjac, chkjac, initd, restrt, isodr, implct, dovcv, redoj &
            , msgb(1), msgb(2), msgd(1), msgd(2), n, m, np, nq, npp, nnzw, x, &
            ldx, ifixx, ldifx, delta, wd, ldwd, ld2wd, tt, ldtt, stpd, ldstpd &
            , y, ldy, we, ldwe, ld2we, pnlty, beta, ifixb, ssf, stpb, lower,  &
            upper, job, neta, taufac, sstol, partol, maxit, wss(1), wss(2),   &
            wss(3))
!
!  Print iteration reports
!
        elseif ( iflag .eq. 2) then
!
           if ( fstitr) then
              write ( lunrpt,1300) typ
           endif
           call dodpc2                                                        &
            ( ipr, lunrpt, fstitr, implct, prtpen,                            &
            pnlty,                                                            &
            niter, nfev, wss(1), actred, prered, alpha, tau, pnorm, np, beta)
!
!  Print final summary
!
        elseif ( iflag .eq. 3) then
!
           write ( lunrpt,1400) typ
           call dodpc3                                                        &
            ( ipr, lunrpt,                                                    &
            isodr, implct, didvcv, dovcv, redoj, anajac,                      &
            n, m, np, nq, npp,                                                &
            info, niter, nfev, njev, irank, rcond, istop,                     &
            wss(1), wss(2), wss(3), pnlty, rvar, idf,                         &
            beta, sdbeta, ifixb, f, delta, lower, upper)
        endif
!
        return
!
!  Format statements
!
1200    format                                                                &
         (/' *** Initial summary for fit by method of ',A3,' ***')
1300    format                                                                &
         (/' *** Iteration reports for fit by method of ',A3,' ***')
1400    format                                                                &
         (/' *** Final summary for fit by method of ',A3,' ***')
!
end subroutine
!DODPE1
subroutine dodpe1                                                             &
         ( unit, info, d1, d2, d3, d4, d5,                                    &
         n, m, nq,                                                            &
         ldscld, ldstpd, ldwe, ld2we, ldwd, ld2wd,                            &
         lwkmn, liwkmn)
!***Begin Prologue  DODPE1
!***Refer to  ODR
!***Routines Called  (NONE)
!***Date Written   860529   (YYMMDD)
!***Revision Date  920619   (YYMMDD)
!***Purpose  Print error reports
!***End Prologue  DODPE1
!
!...Scalar arguments
        integer                                                               &
         d1, d2, d3, d4, d5, info, ldscld, ldstpd, ldwd, ldwe, ld2wd, ld2we,  &
         liwkmn, lwkmn, m, n, nq, unit
!-------------------------------------^----------------------------------------
!!! FPT - 1273 Fortran auxiliary keyword used as identifier name.
!------------------------------------------------------------------------------
!
!...Variable Definitions (alphabetically)
!       D1:      The 1st digit (from the left) of INFO.
!       D2:      The 2nd digit (from the left) of INFO.
!       D3:      The 3rd digit (from the left) of INFO.
!       D4:      The 4th digit (from the left) of INFO.
!       D5:      The 5th digit (from the left) of INFO.
!       INFO:    The variable designating why the computations were stopped.
!       LDSCLD:  The leading dimension of array SCLD.
!       LDSTPD:  The leading dimension of array STPD.
!       LDWD:    The leading dimension of array WD.
!       LDWE:    The leading dimension of array WE.
!       LIWKMN:  The minimum acceptable length of array IWORK.
!       LWKMN:   The minimum acceptable length of array WORK.
!       LD2WD:   The second dimension of array WD.
!       LD2WE:   The second dimension of array WE.
!       M:       The number of columns of data in the explanatory variable.
!       N:       The number of observations.
!       NQ:      The number of responses per observation.
!       UNIT:    The logical unit number used for error messages.
!
!
!***First executable statement  DODPE1
!
!
!  Print appropriate messages for errors in problem specification
!  parameters
!
        if ( d1 .eq. 1) then
           if ( d2 .ne. 0) then
              write ( unit,1100)
           endif
           if ( d3 .ne. 0) then
              write ( unit,1200)
           endif
           if ( d4 .ne. 0) then
              write ( unit,1300)
           endif
           if ( d5 .ne. 0) then
              write ( unit,1400)
           endif
!
!  Print appropriate messages for errors in dimension specification
!  parameters
!
        elseif ( d1 .eq. 2) then
!
           if ( d2 .ne. 0) then
              if ( d2 .eq. 1 .or. d2 .eq. 3) then
                 write ( unit,2110)
              endif
              if ( d2 .eq. 2 .or. d2 .eq. 3) then
                 write ( unit,2120)
              endif
           endif
!
           if ( d3 .ne. 0) then
              if ( d3 .eq. 1 .or. d3 .eq. 3 .or. d3 .eq. 5 .or. d3 .eq. 7)    &
               then
                 write ( unit,2210)
              endif
              if ( d3 .eq. 2 .or. d3 .eq. 3 .or. d3 .eq. 6 .or. d3 .eq. 7)    &
               then
                 write ( unit,2220)
              endif
              if ( d3 .eq. 4 .or. d3 .eq. 5 .or. d3 .eq. 6 .or. d3 .eq. 7)    &
               then
                 write ( unit,2230)
              endif
           endif
!
           if ( d4 .ne. 0) then
              if ( d4 .eq. 1 .or. d4 .eq. 3) then
                 write ( unit,2310)
              endif
              if ( d4 .eq. 2 .or. d4 .eq. 3) then
                 write ( unit,2320)
              endif
           endif
!
           if ( d5 .ne. 0) then
              if ( d5 .eq. 1 .or. d5 .eq. 3) then
                 write ( unit,2410) lwkmn
              endif
              if ( d5 .eq. 2 .or. d5 .eq. 3) then
                 write ( unit,2420) liwkmn
              endif
           endif
!
        elseif ( d1 .eq. 3) then
!
!  Print appropriate messages for errors in scale values
!
           if ( d3 .ne. 0) then
              if ( d3 .eq. 2 .or. d3 .eq. 3) then
                 if ( ldscld .ge. n) then
                    write ( unit,3110)
                 else
                    write ( unit,3120)
                 endif
              endif
              if ( d3 .eq. 1 .or. d3 .eq. 3) then
                 write ( unit,3130)
              endif
           endif
!
!  Print appropriate messages for errors in derivative step values
!
           if ( d2 .ne. 0) then
              if ( d2 .eq. 2 .or. d2 .eq. 3) then
                 if ( ldstpd .ge. n) then
                    write ( unit,3210)
                 else
                    write ( unit,3220)
                 endif
              endif
              if ( d2 .eq. 1 .or. d2 .eq. 3) then
                 write ( unit,3230)
              endif
           endif
!
!  Print appropriate messages for errors in observational error weights
!
           if ( d4 .ne. 0) then
              if ( d4 .eq. 1) then
                 if ( ldwe .ge. n) then
                    if ( ld2we .ge. nq) then
                       write ( unit,3310)
                    else
                       write ( unit,3320)
                    endif
                 else
                    if ( ld2we .ge. nq) then
                       write ( unit,3410)
                    else
                       write ( unit,3420)
                    endif
                 endif
              endif
              if ( d4 .eq. 2) then
                 write ( unit,3500)
              endif
           endif
!
!  Print appropriate messages for errors in DELTA weights
!
           if ( d5 .ne. 0) then
              if ( ldwd .ge. n) then
                 if ( ld2wd .ge. m) then
                    write ( unit,4310)
                 else
                    write ( unit,4320)
                 endif
              else
                 if ( ld2wd .ge. m) then
                    write ( unit,4410)
                 else
                    write ( unit,4420)
                 endif
              endif
           endif
!
        elseif ( d1 .eq. 7) then
!
!  Print the appropriate messages for errors in JOB
!
           if ( d2 .ne. 0) then
              write ( unit,5000)
           endif
!
           if ( d3 .ne. 0) then
              write ( unit,5100)
           endif
!
           if ( d4 .ne. 0) then
              write ( unit,5200)
           endif
!
        elseif ( d1 .eq. 8) then
!
!  Print the appropriate messages for errors in array allocation
!
           if ( d2 .ne. 0) then
              write ( unit,7200)
           endif
!
           if ( d3 .ne. 0) then
              write ( unit,7300)
           endif
!
           if ( d4 .ne. 0) then
              write ( unit,7400)
           endif
!
        elseif ( d1 .eq. 9) then
!
!  Print the appropriate messages for errors in bounds
!
           if ( d2 .ne. 0) then
              write ( unit,6000)
           endif
!
           if ( d3 .ne. 0) then
              write ( unit,6100)
           endif
!
           if ( d4 .eq. 1) then
              write ( unit,6210)
           endif
!
           if ( d4 .eq. 2) then
              write ( unit,6220)
           endif
!
        endif
!
!  Print error messages for array sizes incorrect
!
        if ( info/100000 .eq. 1) then
           info = info-100000
           if ( info .ge. 32768) then
              info = info-32768
              write ( unit,8015)
           endif
           if ( info .ge. 16384) then
              info = info-16384
              write ( unit,8014)
           endif
           if ( info .ge. 8192) then
              info = info-8192
              write ( unit,8013)
           endif
           if ( info .ge. 4096) then
              info = info-4096
              write ( unit,8012)
           endif
           if ( info .ge. 2048) then
              info = info-2048
              write ( unit,8011)
           endif
           if ( info .ge. 1024) then
              info = info-1024
              write ( unit,8010)
           endif
           if ( info .ge. 512) then
              info = info-512
              write ( unit,8009)
           endif
           if ( info .ge. 256) then
              info = info-256
              write ( unit,8008)
           endif
           if ( info .ge. 128) then
              info = info-128
              write ( unit,8007)
           endif
           if ( info .ge. 64) then
              info = info-64
              write ( unit,8006)
           endif
           if ( info .ge. 32) then
              info = info-32
              write ( unit,8005)
           endif
           if ( info .ge. 16) then
              info = info-16
              write ( unit,8004)
           endif
           if ( info .ge. 8) then
              info = info-8
              write ( unit,8003)
           endif
           if ( info .ge. 4) then
              info = info-4
              write ( unit,8002)
           endif
           if ( info .ge. 2) then
              info = info-2
              write ( unit,8001)
           endif
           if ( info .ge. 1) then
              info = info-1
              write ( unit,8000)
           endif
        endif
!
!  Format statements
!
1100    format                                                                &
         (/' ERROR :  N is less than one.')
1200    format                                                                &
         (/' ERROR :  M is less than one.')
1300    format                                                                &
         (/' ERROR :  NP is less than one'/                                   &
         '          or NP is greater than N.')
1400    format                                                                &
         (/' ERROR :  NQ is less than one.')
2110    format                                                                &
         (/' ERROR :  LDX is less than N.')
2120    format                                                                &
         (/' ERROR :  LDY is less than N.')
2210    format                                                                &
         (/' ERROR :  LDIFX is less than N'/                                  &
         '          and LDIFX is not equal to one.')
2220    format                                                                &
         (/' ERROR :  LDSCLD is less than N'/                                 &
         '          and LDSCLD is not equal to one.')
2230    format                                                                &
         (/' ERROR :  LDSTPD is less than N'/                                 &
         '          and LDSTPD is not equal to one.')
2310    format                                                                &
         (/' ERROR :  LDWE is less than N'/                                   &
         '          and LDWE is not equal to one or'/                         &
         '          or'/                                                      &
         '          LD2WE is less than NQ'/                                   &
         '          and LD2WE is not equal to one.')
2320    format                                                                &
         (/' ERROR :  LDWD is less than N'/                                   &
         '          and LDWD is not equal to one.')
2410    format                                                                &
         (/' ERROR :  LWORK is less than ',I7,','/                            &
         '          the smallest acceptable dimension of array WORK.')
2420    format                                                                &
         (/' ERROR :  LIWORK is less than ',I7,','/                           &
         '          the smallest acceptable dimension of array',              &
         ' IWORK.')
3110    format                                                                &
         (/' ERROR :  SCLD(I,J) is less than or equal to zero'/               &
         '          for some I = 1, ..., N and J = 1, ..., M.'//              &
         '          when SCLD(1,1) is greater than zero'/                     &
         '          and LDSCLD is greater than or equal to N then'/           &
         '          each of the N by M elements of'/                          &
         '          SCLD must be greater than zero.')
3120    format                                                                &
         (/' ERROR :  SCLD(1,J) is less than or equal to zero'/               &
         '          for some J = 1, ..., M.'//                                &
         '          when SCLD(1,1) is greater than zero'/                     &
         '          and LDSCLD is equal to one then'/                         &
         '          each of the 1 by M elements of'/                          &
         '          SCLD must be greater than zero.')
3130    format                                                                &
         (/' ERROR :  SCLB(K) is less than or equal to zero'/                 &
         '          for some K = 1, ..., NP.'//                               &
         '          all NP elements of',                                      &
         '          SCLB must be greater than zero.')
3210    format                                                                &
         (/' ERROR :  STPD(I,J) is less than or equal to zero'/               &
         '          for some I = 1, ..., N and J = 1, ..., M.'//              &
         '          when STPD(1,1) is greater than zero'/                     &
         '          and LDSTPD is greater than or equal to N then'/           &
         '          each of the N by M elements of'/                          &
         '          STPD must be greater than zero.')
3220    format                                                                &
         (/' ERROR :  STPD(1,J) is less than or equal to zero'/               &
         '          for some J = 1, ..., M.'//                                &
         '          when STPD(1,1) is greater than zero'/                     &
         '          and LDSTPD is equal to one then'/                         &
         '          each of the 1 by M elements of'/                          &
         '          STPD must be greater than zero.')
3230    format                                                                &
         (/' ERROR :  STPB(K) is less than or equal to zero'/                 &
         '          for some K = 1, ..., NP.'//                               &
         '          all NP elements of',                                      &
         ' STPB must be greater than zero.')
3310    format                                                                &
         (/' ERROR :  At least one of the (NQ by NQ) arrays starting'/        &
         '          in WE(I,1,1), I = 1, ..., N, is not positive'/            &
         '          semidefinite.  When WE(1,1,1) is greater than'/           &
         '          or equal to zero, and LDWE is greater than or'/           &
         '          equal to N, and LD2WE is greater than or equal'/          &
         '          to NQ, then each of the (NQ by NQ) arrays in WE'/         &
         '          must be positive semidefinite.')
3320    format                                                                &
         (/' ERROR :  At least one of the (1 by NQ) arrays starting'/         &
         '          in WE(I,1,1), I = 1, ..., N, has a negative'/             &
         '          element.  When WE(1,1,1) is greater than or'/             &
         '          equal to zero, and LDWE is greater than or equal'/        &
         '          to N, and LD2WE is equal to 1, then each of the'/         &
         '          (1 by NQ) arrays in WE must have only non-'/              &
         '          negative elements.')
3410    format                                                                &
         (/' ERROR :  The (NQ by NQ) array starting in WE(1,1,1) is'/         &
         '          not positive semidefinite.  When WE(1,1,1) is'/           &
         '          greater than or equal to zero, and LDWE is equal'/        &
         '          to 1, and LD2WE is greater than or equal to NQ,'/         &
         '          then the (NQ by NQ) array in WE must be positive'/        &
         '          semidefinite.')
3420    format                                                                &
         (/' ERROR :  The (1 by NQ) array starting in WE(1,1,1) has'/         &
         '          a negative element.  When WE(1,1,1) is greater'/          &
         '          than or equal to zero, and LDWE is equal to 1,'/          &
         '          and LD2WE is equal to 1, then the (1 by NQ)'/             &
         '          array in WE must have only nonnegative elements.')
3500    format                                                                &
         (/' ERROR :  The number of nonzero arrays in array WE is'/           &
         '          less than NP.')
4310    format                                                                &
         (/' ERROR :  At least one of the (M by M) arrays starting'/          &
         '          in WD(I,1,1), I = 1, ..., N, is not positive'/            &
         '          definite.  When WD(1,1,1) is greater than zero,'/         &
         '          and LDWD is greater than or equal to N, and'/             &
         '          LD2WD is greater than or equal to M, then each'/          &
         '          of the (M by M) arrays in WD must be positive'/           &
         '          definite.')
4320    format                                                                &
         (/' ERROR :  At least one of the (1 by M) arrays starting'/          &
         '          in WD(I,1,1), I = 1, ..., N, has a nonpositive'/          &
         '          element.  When WD(1,1,1) is greater than zero,'/          &
         '          and LDWD is greater than or equal to N, and'/             &
         '          LD2WD is equal to 1, then each of the (1 by M)'/          &
         '          arrays in WD must have only positive elements.')
4410    format                                                                &
         (/' ERROR :  The (M by M) array starting in WD(1,1,1) is'/           &
         '          not positive definite.  When WD(1,1,1) is'/               &
         '          greater than zero, and LDWD is equal to 1, and'/          &
         '          LD2WD is greater than or equal to M, then the'/           &
         '          (M by M) array in WD must be positive definite.')
4420    format                                                                &
         (/' ERROR :  The (1 by M) array starting in WD(1,1,1) has a'/        &
         '          nonpositive element.  When WD(1,1,1) is greater'/         &
         '          than zero, and LDWD is equal to 1, and LD2WD is'/         &
         '          equal to 1, then the (1 by M) array in WD must'/          &
         '          have only positive elements.')
5000    format                                                                &
         (/' ERROR :  JOB requires the optional argument DELTA and'/          &
         '          DELTA is not present or not associated.')
5100    format                                                                &
         (/' ERROR :  JOB requires the optional argument WORK and'/           &
         '          WORK is not present or not associated.')
5200    format                                                                &
         (/' ERROR :  JOB requires the optional argument IWORK and'/          &
         '          IWORK is not present or not associated.')
6000    format                                                                &
         (/' ERROR :  LOWER(K).GT.UPPER(K) for some K.  Adjust the'/          &
         '          the bounds so that LOWER(K).LE.UPPER(K) holds'/           &
         '          for all K.')
6100    format                                                                &
         (/' ERROR :  BETA(K).GT.UPPER(K) or BETA(K).LT.LOWER(K) '/           &
         '          for some K.  Adjust the bounds or BETA so '/              &
         '          that LOWER(K).LE.BETA(K).LE.UPPER(K) holds'/              &
         '          for all K.')
6210    format                                                                &
         (/' ERROR :  UPPER(K)-LOWER(K) .LT. 400*BETA(K)*EPSMAC  '/           &
         '          for some K and EPSMAC having the largest '/               &
         '          value such that 1+EPSMAC.NE.1.  This '/                   &
         '          constraint on UPPER and LOWER is necessary'/              &
         '          for the calculation of NDIGIT.  Increase the'/            &
         '          range of the bounds or specify NDIGIT '/                  &
         '          explicitly.')
6220    format                                                                &
         (/' ERROR :  UPPER(K)-LOWER(K) .LT. ABS(STEP) for some'/             &
         '          K where step is the step size for numeric'/               &
         '          derivatives.  Increase the bounds or supply'/             &
         '          an analytic jacobian.')
7200    format                                                                &
         (/' ERROR :  DELTA could not be allocated. ')
7300    format                                                                &
         (/' ERROR :  WORK could not be allocated. ')
7400    format                                                                &
         (/' ERROR :  IWORK could not be allocated. ')
8000    format                                                                &
         (/' ERROR :  BETA has incorrect size. ')
8001    format                                                                &
         (/' ERROR :  Y has incorrect size. ')
8002    format                                                                &
         (/' ERROR :  X has incorrect size. ')
8003    format                                                                &
         (/' ERROR :  DELTA has incorrect size. ')
8004    format                                                                &
         (/' ERROR :  WE has incorrect size. ')
8005    format                                                                &
         (/' ERROR :  WD has incorrect size. ')
8006    format                                                                &
         (/' ERROR :  IFIXB has incorrect size. ')
8007    format                                                                &
         (/' ERROR :  IFIXX has incorrect size. ')
8008    format                                                                &
         (/' ERROR :  STPB has incorrect size. ')
8009    format                                                                &
         (/' ERROR :  STPD has incorrect size. ')
8010    format                                                                &
         (/' ERROR :  SCLB has incorrect size. ')
8011    format                                                                &
         (/' ERROR :  SCLD has incorrect size. ')
8012    format                                                                &
         (/' ERROR :  WORK has incorrect size. ')
8013    format                                                                &
         (/' ERROR :  IWORK has incorrect size. ')
8014    format                                                                &
         (/' ERROR :  UPPER has incorrect size. ')
8015    format                                                                &
         (/' ERROR :  LOWER has incorrect size. ')
end subroutine
!DODPE2
subroutine dodpe2                                                             &
         ( unit,                                                              &
         n, m, np, nq,                                                        &
         fjacb, fjacd,                                                        &
         diff, msgb1, msgb, isodr, msgd1, msgd,                               &
         xplusd, nrow, neta, ntol)
!***Begin Prologue  DODPE2
!***Refer to  ODR
!***Routines Called  (NONE)
!***Date Written   860529   (YYMMDD)
!***Revision Date  920619   (YYMMDD)
!***Purpose  Generate the derivative checking report
!***End Prologue  DODPE2
!
!...Used modules
        use odrpack_kinds,only: wp
!
!...Scalar arguments
        integer                                                               &
         m, msgb1, msgd1, n, neta, np, nq, nrow, ntol, unit
!----------------------------------------------------------^-------------------
!!! FPT - 1273 Fortran auxiliary keyword used as identifier name.
!------------------------------------------------------------------------------
        logical                                                               &
         isodr
!
!...Array arguments
        real(kind = wp)                                                       &
         diff( nq, np+ m), fjacb( n, np, nq), fjacd( n, m, nq), xplusd( n, m)
        integer                                                               &
         msgb( nq, np), msgd( nq, m)
!
!...Local scalars
        integer                                                               &
         i, j, k, l
        character flag*1, typ*3
!
!...Local arrays
        logical                                                               &
         ftnote(0:9)
!
!...Variable Definitions (alphabetically)
!       DIFF:    The relative differences between the user supplied and
!       finite difference derivatives for each derivative checked.
!       FJACB:   The Jacobian with respect to BETA.
!       FJACD:   The Jacobian with respect to DELTA.
!       FLAG:    The character string indicating highly questionable results.
!       FTNOTE:  The array controling footnotes.
!       I:       An index variable.
!       ISODR:   The variable designating whether the solution is by ODR
!       (ISODR=.TRUE.) or by OLS (ISODR=.FALSE.).
!       J:       An index variable.
!       K:       An index variable.
!       L:       An index variable.
!       M:       The number of columns of data in the explanatory variable.
!       MSGB:    The error checking results for the Jacobian wrt BETA.
!       MSGB1:   The error checking results for the Jacobian wrt BETA.
!       MSGD:    The error checking results for the Jacobian wrt DELTA.
!       MSGD1:   The error checking results for the Jacobian wrt DELTA.
!       N:       The number of observations.
!       NETA:    The number of reliable digits in the model.
!       NP:      The number of function parameters.
!       NQ:      The number of responses per observation.
!       NROW:    The row number of the explanatory variable array at
!       which the derivative is to be checked.
!       NTOL:    The number of digits of agreement required between the
!       finite difference and the user supplied derivatives.
!       TYP:     The character string indicating solution type, ODR or OLS.
!       UNIT:    The logical unit number used for error messages.
!       XPLUSD:  The values of X + DELTA.
!
!
!***First executable statement  DODPE2
!
!
!  Set up for footnotes
!
        do 10 i = 0,9
           ftnote( i) = .false.
10      continue
!
        do 40 l = 1, nq
           if ( msgb1 .ge. 1) then
              do 20 i = 1, np
                 if ( msgb( l, i) .ge. 1) then
                    ftnote(0) = .true.
                    ftnote( msgb( l, i)) = .true.
                 endif
20            continue
           endif
!
           if ( msgd1 .ge. 1) then
              do 30 i = 1, m
                 if ( msgd( l, i) .ge. 1) then
                    ftnote(0) = .true.
                    ftnote( msgd( l, i)) = .true.
                 endif
30            continue
           endif
40      continue
!
!       Print report
!
        if ( isodr) then
           typ = 'ODR'
        else
           typ = 'OLS'
        endif
        write ( unit,1000) typ
!
        do 70 l = 1, nq
!
           write ( unit,2100) l, nrow
           write ( unit,2200)
!
           do 50 i = 1, np
              k = msgb( l, i)
              if ( k .eq. 7) then
                 flag = '*'
              else
                 flag = ' '
              endif
              if ( k .le. -1) then
                 write ( unit,3100) i
              elseif ( k .eq. 0) then
                 write ( unit,3200) i, fjacb( nrow, i, l), diff( l, i), flag
              elseif ( k .eq. 8) then
                 write ( unit,3400) i, fjacb( nrow, i, l), flag, k
              elseif ( k .eq. 9) then
                 write ( unit,3500) i, flag, k
              elseif ( k .ge. 1) then
                 write ( unit,3300) i, fjacb( nrow, i, l), diff( l, i), flag, &
                  k
              endif
50         continue
           if ( isodr) then
              do 60 i = 1, m
                 k = msgd( l, i)
                 if ( k .eq. 7) then
                    flag = '*'
                 else
                    flag = ' '
                 endif
                 if ( k .le. -1) then
                    write ( unit,4100) nrow, i
                 elseif ( k .eq. 0) then
                    write ( unit,4200) nrow, i,                               &
                     fjacd( nrow, i, l), diff( l, np+ i), flag
                 elseif ( k .ge. 1) then
                    write ( unit,4300) nrow, i,                               &
                     fjacd( nrow, i, l), diff( l, np+ i), flag, k
                 endif
60            continue
           endif
70      continue
!
!       Print footnotes
!
        if ( ftnote(0)) then
!
           write ( unit,5000)
           if ( ftnote(1))write ( unit,5100)
           if ( ftnote(2))write ( unit,5200)
           if ( ftnote(3))write ( unit,5300)
           if ( ftnote(4))write ( unit,5400)
           if ( ftnote(5))write ( unit,5500)
           if ( ftnote(6))write ( unit,5600)
           if ( ftnote(7))write ( unit,5700)
           if ( ftnote(8))write ( unit,5800)
           if ( ftnote(9))write ( unit,5900)
        endif
!
        if ( neta .lt. 0) then
           write ( unit,6000)- neta
        else
           write ( unit,6100) neta
        endif
        write ( unit,7000) ntol
!
!  Print out row of explanatory variable which was checked.
!
        write ( unit,8100) nrow
!
        do 80 j = 1, m
           write ( unit,8110) nrow, j, xplusd( nrow, j)
80      continue
!
        return
!
!       Format statements
!
1000    format                                                                &
         (//' *** Derivative checking report for fit by method of ',A3,       &
         ' ***'/)
2100    format (/'     For response ',I2,' of observation ',I5/)
2200    format ('                      ','         User',                     &
         '               ','                '/                                &
         '                      ','     Supplied',                            &
         '     Relative','    Derivative '/                                   &
         '        Derivative WRT','        Value',                            &
         '   Difference','    Assessment '/)
3100    format ('             BETA(',I3,')','       ---   ',                  &
         '       ---   ','    Unchecked')
3200    format ('             BETA(',I3,')',1P,2E13.2,3X,A1,                  &
         'Verified')
3300    format ('             BETA(',I3,')',1P,2E13.2,3X,A1,                  &
         'Questionable (see note ',I1,')')
3400    format ('             BETA(',I3,')',1P,1E13.2,13X,3X,A1,              &
         'Questionable (see note ',I1,')')
3500    format ('             BETA(',I3,')',1P,13X,13X,3X,A1,                 &
         'Small bounds (see note ',I1,')')
4100    format ('          DELTA(',I2,',',I2,')','       ---   ',             &
         '       ---   ','    Unchecked')
4200    format ('          DELTA(',I2,',',I2,')',1P,2E13.2,3X,A1,             &
         'Verified')
4300    format ('          DELTA(',I2,',',I2,')',1P,2E13.2,3X,A1,             &
         'Questionable (see note ',I1,')')
5000    format                                                                &
         (/'     NOTES:')
5100    format                                                                &
         (/'      (1) User supplied and finite difference derivatives',       &
         ' agree, but'/                                                       &
         '          results are questionable because both are zero.')
5200    format                                                                &
         (/'      (2) User supplied and finite difference derivatives',       &
         ' agree, but'/                                                       &
         '          results are questionable because one is',                 &
         ' identically zero'/                                                 &
         '          and the other is only approximately zero.')
5300    format                                                                &
         (/'      (3) User supplied and finite difference derivatives',       &
         ' disagree, but'/                                                    &
         '          results are questionable because one is',                 &
         ' identically zero'/                                                 &
         '          and the other is not.')
5400    format                                                                &
         (/'      (4) User supplied and finite difference derivatives',       &
         ' disagree, but'/                                                    &
         '          finite difference derivative is questionable',            &
         ' because either'/                                                   &
         '          the ratio of relative curvature to relative',             &
         ' slope is too high'/                                                &
         '          or the scale is wrong.')
5500    format                                                                &
         (/'      (5) User supplied and finite difference derivatives',       &
         ' disagree, but'/                                                    &
         '          finite difference derivative is questionable',            &
         ' because the'/                                                      &
         '          ratio of relative curvature to relative slope is',        &
         ' too high.')
5600    format                                                                &
         (/'      (6) User supplied and finite difference derivatives',       &
         ' disagree, but'/                                                    &
         '          have at least 2 digits in common.')
5700    format                                                                &
         (/'      (7) User supplied and finite difference derivatives',       &
         ' disagree, and'/                                                    &
         '          have fewer than 2 digits in common.  derivative',         &
         ' checking must'/                                                    &
         '          be turned off in order to proceed.')
5800    format                                                                &
         (/'      (8) User supplied and finite difference derivatives',       &
         ' disagree, and'/                                                    &
         '          bound constraints are too small to calculate',            &
         ' further'/                                                          &
         '          information.')
5900    format                                                                &
         (/'      (9) Bound constraints too small to check derivative.')
6000    format                                                                &
         (/'     Number of reliable digits in function results       ',       &
         I5/                                                                  &
         '        (estimated by ODRPACK95)')
6100    format                                                                &
         (/'     Number of reliable digits in function results       ',       &
         I5/                                                                  &
         '        (supplied by user)')
7000    format                                                                &
         (/'     Number of digits of agreement required between      '/       &
         '     user supplied and finite difference derivative for  '/         &
         '     user supplied derivative to be considered verified  ',         &
         I5)
8100    format                                                                &
         (/'     Row number at which derivatives were checked        ',       &
         I5//                                                                 &
         '       -values of the explanatory variables at this row'/)
8110    format                                                                &
         (10X,'X(',I2,',',I2,')',1X,1P,3E16.8)
end subroutine
!DODPE3
subroutine dodpe3                                                             &
         ( unit, d2, d3)
!***Begin Prologue  DODPE3
!***Refer to  ODR
!***Routines Called  (NONE)
!***Date Written   860529   (YYMMDD)
!***Revision Date  920619   (YYMMDD)
!***Purpose  Print error reports indicating that computations were
!       stopped in user supplied subroutines FCN
!***End Prologue  DODPE3
!
!...Scalar arguments
        integer                                                               &
         d2, d3, unit
!--------------------^---------------------------------------------------------
!!! FPT - 1273 Fortran auxiliary keyword used as identifier name.
!------------------------------------------------------------------------------
!
!...Variable Definitions (alphabetically)
!       D2:      The 2nd digit (from the left) of INFO.
!       D3:      The 3rd digit (from the left) of INFO.
!       UNIT:    The logical unit number used for error messages.
!
!
!***First executable statement  DODPE3
!
!
!  Print appropriate messages to indicate where computations were
!  stopped
!
        if ( d2 .eq. 2) then
           write ( unit,1100)
        elseif ( d2 .eq. 3) then
           write ( unit,1200)
        elseif ( d2 .eq. 4) then
           write ( unit,1300)
        endif
        if ( d3 .eq. 2) then
           write ( unit,1400)
        endif
!
!  Format statements
!
1100    format                                                                &
         (//' Variable ISTOP has been returned with a nonzero value  '/       &
         ' from user supplied subroutine FCN when invoked using the'/         &
         ' initial estimates of BETA and DELTA supplied by the     '/         &
         ' user.  The initial estimates must be adjusted to allow  '/         &
         ' proper evaluation of subroutine FCN before the          '/         &
         ' regression procedure can continue.')
1200    format                                                                &
         (//' Variable ISTOP has been returned with a nonzero value  '/       &
         ' from user supplied subroutine FCN.  This occurred during'/         &
         ' the computation of the number of reliable digits in the '/         &
         ' predicted values (F) returned from subroutine FCN, indi-'/         &
         ' cating that changes in the initial estimates of BETA(K),'/         &
         ' K=1,NP, as small as 2*BETA(K)*SQRT(MACHINE PRECISION),  '/         &
         ' where MACHINE PRECISION is defined as the smallest value'/         &
         ' E such that 1+E>1 on the computer being used, prevent   '/         &
         ' subroutine FCN from being properly evaluated.  The      '/         &
         ' initial estimates must be adjusted to allow proper      '/         &
         ' evaluation of subroutine FCN during these computations  '/         &
         ' before the regression procedure can continue.')
1300    format                                                                &
         (//' Variable ISTOP has been returned with a nonzero value  '/       &
         ' from user supplied subroutine FCN.  This occurred during'/         &
         ' the derivative checking procedure, indicating that      '/         &
         ' changes in the initial estimates of BETA(K), K=1,NP, as '/         &
         ' small as MAX[BETA(K),1/SCLB(K)]*10**(-NETA/2), and/or   '/         &
         ' of DELTA(I,J), I=1,N and J=1,M, as small as             '/         &
         ' MAX[DELTA(I,J),1/SCLD(I,J)]*10**(-NETA/2), where NETA   '/         &
         ' is defined to be the number of reliable digits in       '/         &
         ' predicted values (F) returned from subroutine FCN,      '/         &
         ' prevent subroutine FCN from being properly evaluated.   '/         &
         ' the initial estimates must be adjusted to allow proper  '/         &
         ' evaluation of subroutine FCN during these computations  '/         &
         ' before the regression procedure can continue.')
1400    format                                                                &
         (//' Variable ISTOP has been returned with a nonzero value  '/       &
         ' from user supplied subroutine FCN when invoked for '/              &
         ' derivative evaluations using the initial estimates of '/           &
         ' BETA and DELTA supplied by the user.  The initial '/               &
         ' estimates must be adjusted to allow proper evaluation '/           &
         ' of subroutine FCN before the regression procedure can '/           &
         ' continue.')
end subroutine
!DODPER
subroutine dodper                                                             &
         ( info, lunerr,                                                      &
         n, m, np, nq,                                                        &
         ldscld, ldstpd, ldwe, ld2we, ldwd, ld2wd,                            &
         lwkmn, liwkmn,                                                       &
         fjacb, fjacd,                                                        &
         diff, msgb, isodr, msgd,                                             &
         xplusd, nrow, neta, ntol)
!***Begin Prologue  DODPER
!***Refer to  ODR
!***Routines Called  DODPE1,DODPE2,DODPE3,DODPHD
!***Date Written   860529   (YYMMDD)
!***Revision Date  920619   (YYMMDD)
!***Purpose  Controlling routine for printing error reports
!***End Prologue  DODPER
!
!...Used modules
        use odrpack_kinds,only: wp
!
!...Scalar arguments
        integer                                                               &
         info, ldscld, ldstpd, ldwd, ldwe, ld2wd, ld2we, liwkmn, lunerr,      &
         lwkmn, m, n, neta, np, nq, nrow, ntol
        logical                                                               &
         isodr
!
!...Array arguments
        real(kind = wp)                                                       &
         diff( nq, np+ m), fjacb( n, np, nq), fjacd( n, m, nq), xplusd( n, m)
        integer                                                               &
         msgb( nq* np+1), msgd( nq* m+1)
!
!...Local scalars
        integer                                                               &
         d1, d2, d3, d4, d5, unit
!--------------------------------^---------------------------------------------
!!! FPT - 1273 Fortran auxiliary keyword used as identifier name.
!------------------------------------------------------------------------------
        logical                                                               &
         head
!
!...External subroutines
        external                                                              &
         dodpe1, dodpe2, dodpe3, dodphd
!
!...Variable Definitions (alphabetically)
!       D1:      The 1st digit (from the left) of INFO.
!       D2:      The 2nd digit (from the left) of INFO.
!       D3:      The 3rd digit (from the left) of INFO.
!       D4:      The 4th digit (from the left) of INFO.
!       D5:      The 5th digit (from the left) of INFO.
!       DIFF:    The relative differences between the user supplied and
!       finite difference derivatives for each derivative checked.
!       FJACB:   The Jacobian with respect to BETA.
!       FJACD:   The Jacobian with respect to DELTA.
!       HEAD:    The variable designating whether the heading is to be
!       printed (HEAD=.TRUE.) or not (HEAD=.FALSE.).
!       INFO:    The variable designating why the computations were stopped.
!       ISODR:   The variable designating whether the solution is by ODR
!       (ISODR=.TRUE.) or by OLS (ISODR=.FALSE.).
!       LDSCLD:  The leading dimension of array SCLD.
!       LDSTPD:  The leading dimension of array STPD.
!       LDWD:    The leading dimension of array WD.
!       LDWE:    The leading dimension of array WE.
!       LD2WD:   The second dimension of array WD.
!       LD2WE:   The second dimension of array WE.
!       LIWKMN:  The minimum acceptable length of array IWORK.
!       LUNERR:  The logical unit number used for error messages.
!       LWKMN:   The minimum acceptable length of array WORK.
!       M:       The number of columns of data in the explanatory variable.
!       MSGB:    The error checking results for the Jacobian wrt BETA.
!       MSGD:    The error checking results for the Jacobian wrt DELTA.
!       N:       The number of observations.
!       NETA:    The number of reliable digits in the model.
!       NP:      The number of function parameters.
!       NQ:      The number of responses per observation.
!       NROW:    The row number of the explanatory variable array at
!       which the derivative is to be checked.
!       NTOL:    The number of digits of agreement required between the
!       finite difference and the user supplied derivatives.
!       UNIT:    The logical unit number for error messages.
!       XPLUSD:  The values X + DELTA.
!
!
!***First executable statement  DODPER
!
!
!  Set logical unit number for error report
!
        if ( lunerr .eq. 0) then
           return
        elseif ( lunerr .lt. 0) then
           unit = 6
        else
           unit = lunerr
        endif
!
!  Print heading
!
        head = .true.
        call dodphd( head, unit)
!
!  Extract individual digits from variable INFO
!
        d1 = mod( info,100000)/10000
        d2 = mod( info,10000)/1000
        d3 = mod( info,1000)/100
        d4 = mod( info,100)/10
        d5 = mod( info,10)
!
!  Print appropriate error messages for ODRPACK95 invoked stop
!
        if (( d1 .ge. 1 .and. d1 .le. 3) .or.                                 &
            ( d1 .eq. 7 .or. d1 .eq. 9)) then
!
!  Print appropriate messages for errors in
!          problem specification parameters
!          dimension specification parameters
!          number of good digits in X
!          weights
!
           call dodpe1( unit, info, d1, d2, d3, d4, d5,                       &
            n, m, nq,                                                         &
            ldscld, ldstpd, ldwe, ld2we, ldwd, ld2wd,                         &
            lwkmn, liwkmn)
!
        elseif (( d1 .eq. 4) .or.                                             &
                ( msgb(1) .ge. 0)) then
!
!  Print appropriate messages for derivative checking
!
           call dodpe2( unit,                                                 &
            n, m, np, nq,                                                     &
            fjacb, fjacd,                                                     &
            diff, msgb(1), msgb(2), isodr, msgd(1), msgd(2),                  &
            xplusd, nrow, neta, ntol)
!
        elseif ( d1 .eq. 5) then
!
!  Print appropriate error message for user invoked stop from FCN
!
           call dodpe3( unit, d2, d3)
!
        endif
!
!  Print correct form of call statement
!
        if (( d1 .ge. 1 .and. d1 .le. 3) .or.                                 &
            ( d1 .eq. 4 .and.                                                 &
             ( d2 .eq. 2 .or. d3 .eq. 2)) .or.                                &
            ( d1 .eq. 5)) then
           write ( unit,1100)
        endif
!
        return
!
!  Format statements
!
1100    format                                                                &
         (//' The correct form of the call statement is '//                   &
         '       CALL ODR'/                                                   &
         '      +     (FCN,'/                                                 &
         '      +     N,M,NP,NQ,'/                                            &
         '      +     BETA,'/                                                 &
         '      +     Y,X,'/                                                  &
         '      +     DELTA*,'/                                               &
         '      +     WE*,WD*,'/                                              &
         '      +     IFIXB*,IFIXX*,'/                                        &
         '      +     JOB*,NDIGIT*,TAUFAC*,'/                                 &
         '      +     SSTOL*,PARTOL*,MAXIT*,'/                                &
         '      +     IPRINT*,LUNERR*,LUNRPT*,'/                              &
         '      +     STPB*,STPD*,'/                                          &
         '      +     SCLB*,SCLD*,'/                                          &
         '      +     WORK*,IWORK*,'/                                         &
         '      +     INFO*,'/                                                &
         '      +     LOWER*,UPPER*)'/                                        &
         ' * optional argument')
!
end subroutine
!DODPHD
subroutine dodphd                                                             &
         ( head, unit)
!***Begin Prologue  DODPHD
!***Refer to  ODR
!***Routines Called  (NONE)
!***Date Written   860529   (YYMMDD)
!***Revision Date  920619   (YYMMDD)
!***Purpose  Print ODRPACK95 heading
!***End Prologue  DODPHD
!
!...Scalar arguments
        integer                                                               &
         unit
!------------^-----------------------------------------------------------------
!!! FPT - 1273 Fortran auxiliary keyword used as identifier name.
!------------------------------------------------------------------------------
        logical                                                               &
         head
!
!...Variable Definitions (alphabetically)
!       HEAD:    The variable designating whether the heading is to be
!       printed (HEAD=.TRUE.) or not (HEAD=.FALSE.).
!       UNIT:    The logical unit number to which the heading is written.
!
!
!***First executable statement  DODPHD
!
!
        if ( head) then
           write ( unit,1000)
           head = .false.
        endif
!
        return
!
!       Format statements
!
1000    format (                                                              &
         ' ********************************************************* '/       &
         ' * ODRPACK95 version 1.00 of 12-27-2005 (REAL (KIND=wp)) * '/       &
         ' ********************************************************* '/)
end subroutine

subroutine dodstp                                                             &
         ( n, m, np, nq, npp,                                                 &
         f, fjacb, fjacd,                                                     &
         wd, ldwd, ld2wd, ss, tt, ldtt, delta,                                &
         alpha, epsfcn, isodr,                                                &
         tfjacb, omega, u, qraux, kpvt,                                       &
         s, t, phi, irank, rcond, forvcv,                                     &
         wrk1, wrk2, wrk3, wrk4, wrk5, wrk, lwrk, istopc)
!***Begin Prologue  DODSTP
!***Refer to  ODR
!***Routines Called  IDAMAX,DCHEX,DESUBI,DFCTR,DNRM2,DQRDC,DQRSL,DROT,
!       DROTG,DSOLVE,DTRCO,DTRSL,DVEVTR,DWGHT
!***Date Written   860529   (YYMMDD)
!***Revision Date  920619   (YYMMDD)
!***Purpose  Compute locally constrained steps S and T, and PHI(ALPHA)
!***End Prologue  DODSTP
!
!...Used modules
        use odrpack_kinds, only: wp, zero, one
        use odrpack, only: tempret
!
!...Scalar arguments
        real(kind = wp)                                                       &
         alpha, epsfcn, phi, rcond
        integer                                                               &
         irank, istopc, ldtt, ldwd, ld2wd, lwrk, m, n, np, npp, nq
        logical                                                               &
         isodr
!
!...Array arguments
        real(kind = wp)                                                       &
         delta( n, m), f( n, nq), fjacb( n, np, nq), fjacd( n, m, nq),        &
         omega( nq, nq), qraux( np), s( np), ss( np),                         &
         t( n, m), tfjacb( n, nq, np), tt( ldtt, m), u( np), wd( ldwd, ld2wd, &
         m), wrk1( n, nq, m), wrk2( n, nq), wrk3( np), wrk4( m, m), wrk5( m), &
         wrk( lwrk)
        integer                                                               &
         kpvt( np)

!...Local scalars
        real(kind = wp) :: co, si, temp
        integer :: i, imax, inf, ipvt, j, k, k1, k2, kp, l
        logical :: elim, forvcv

!...LOCAL ARRAYS
        real(kind = wp) :: dum(2)

!...External functions
        real(kind = wp), external :: dnrm2
        integer, external :: idamax

!...External subroutines
        external :: dchex, desubi, dfctr, dqrdc, dqrsl, drot, drotg, dsolve, dtrco, dtrsl, &
                    dvevtr

!...Interface blocks
        interface
           subroutine dwght                                                   &
            ( n, m, wt, ldwt, ld2wt, t, wtt)
           use odrpack_kinds,only: wp
           integer                                                            &
            ldwt, ld2wt, m, n
           real(kind = wp)                                                    &
            t(:,:), wt(:,:,:), wtt(:,:)
           end subroutine
        end interface

!...Variable definitions (alphabetically)
!       ALPHA:   The Levenberg-Marquardt parameter.
!       CO:      The cosine from the plane rotation.
!       DELTA:   The estimated errors in the explanatory variables.
!       DUM:     A dummy array.
!       ELIM:    The variable designating whether columns of the Jacobian
!       wrt BETA have been eliminated (ELIM=TRUE) or not
!       (ELIM=FALSE).
!       EPSFCN:  The function's precision.
!       F:       The (weighted) estimated values of EPSILON.
!       FJACB:   The Jacobian with respect to BETA.
!       FJACD:   The Jacobian with respect to DELTA.
!       FORVCV:  The variable designating whether this subroutine was
!       called to set up for the covariance matrix computations
!       (FORVCV=TRUE) or not (FORVCV=FALSE).
!       I:       An indexing variable.
!       IMAX:    The index of the element of U having the largest absolute
!       value.
!       INF:     The return code from LINPACK routines.
!       IPVT:    The variable designating whether pivoting is to be done.
!       IRANK:   The rank deficiency of the Jacobian wrt BETA.
!       ISODR:   The variable designating whether the solution is by ODR
!       (ISODR=TRUE) or by OLS (ISODR=FALSE).
!       ISTOPC:  The variable designating whether the computations were
!       stoped due to a numerical error within subroutine DODSTP.
!       J:       An indexing variable.
!       K:       An indexing variable.
!       K1:      An indexing variable.
!       K2:      An indexing variable.
!       KP:      The rank of the Jacobian wrt BETA.
!       KPVT:    The pivot vector.
!       L:       An indexing variable.
!       LDTT:    The leading dimension of array TT.
!       LDWD:    The leading dimension of array WD.
!       LD2WD:   The second dimension of array WD.
!       LWRK:    The length of vector WRK.
!       M:       The number of columns of data in the explanatory variable.
!       N:       The number of observations.
!       NP:      The number of function parameters.
!       NPP:     The number of function parameters being estimated.
!       OMEGA:   The array defined S.T.
!       OMEGA*trans(OMEGA) = inv(I+FJACD*inv(E)*trans(FJACD))
!       = (I-FJACD*inv(P)*trans(FJACD))
!       where E = D**2 + ALPHA*TT**2
!       P = trans(FJACD)*FJACD + D**2 + ALPHA*TT**2
!       ONE:     The value 1.0E0_wp.
!       PHI:     The difference between the norm of the scaled step
!       And the trust region diameter.
!       QRAUX:   The array required to recover the orthogonal part of the
!       Q-R decomposition.
!       RCOND:   The approximate reciprocal condition number of TFJACB.
!       S:       The step for BETA.
!       SI:      The sine from the plane rotation.
!       SS:      The scaling values for the unfixed BETAS.
!       T:       The step for DELTA.
!       TEMP:    A temporary storage LOCATION.
!       TFJACB:  The array OMEGA*FJACB.
!       TT:      The scaling values for DELTA.
!       U:       The approximate null vector for TFJACB.
!       WD:      The (squared) DELTA weights.
!       WRK:     A work array of (LWRK) elements,
!       equivalenced to WRK1 and WRK2.
!       WRK1:    A work array of (N by NQ by M) elements.
!       WRK2:    A work array of (N by NQ) elements.
!       WRK3:    A work array of (NP) elements.
!       WRK4:    A work array of (M by M) elements.
!       WRK5:    A work array of (M) elements.
!       ZERO:    The value 0.0E0_wp.
!
!
!***First executable statement  DODSTP
!
!
!  Compute loop parameters which depend on weight structure
!
!  Set up KPVT if ALPHA = 0
!
        if ( alpha .eq. zero) then
!-----------------------^------------------------------------------------------
!!! FPT - 3087 REAL or COMPLEX quantity tested for exact equality/inequality
!------------------------------------------------------------------------------
           kp = npp
           do 10 k = 1, np
              kpvt( k) = k
10         continue
        else
           if ( npp .ge. 1) then
              kp = npp- irank
           else
              kp = npp
           endif
        endif
!
        if ( isodr) then
!
!  T = WD * DELTA = D*G2
           call dwght( n, m, wd, ldwd, ld2wd, delta, t)
!
           do 300 i = 1, n
!
!  Compute WRK4, such that
!             TRANS(WRK4)*WRK4 = E = (D**2 + ALPHA*TT**2)
              call desubi( n, m, wd, ldwd, ld2wd, alpha, tt, ldtt, i, wrk4)
              call dfctr(.false., wrk4, m, m, inf)
              if ( inf .ne. 0) then
                 istopc = 60000
                 return
              endif
!
!  Compute OMEGA, such that
!             trans(OMEGA)*OMEGA = I+FJACD*inv(E)*trans(FJACD)
!             inv(trans(OMEGA)*OMEGA) = I-FJACD*inv(P)*trans(FJACD)
              call dvevtr( m, nq, i,                                          &
               fjacd, n, m, wrk4, m, wrk1, n, nq, omega, nq, wrk5)
              do 110 l = 1, nq
                 omega( l, l) = one+ omega( l, l)
110           continue
              call dfctr(.false., omega, nq, nq, inf)
              if ( inf .ne. 0) then
                 istopc = 60000
                 return
              endif
!
!  Compute WRK1 = trans(FJACD)*(I-FJACD*inv(P)*trans(JFACD))
!             = trans(FJACD)*inv(trans(OMEGA)*OMEGA)
              do 130 j = 1, m
                 do 120 l = 1, nq
                    wrk1( i, l, j) = fjacd( i, j, l)
120              continue
                 call dsolve( nq, omega, nq, wrk1( i,1: nq, j),4)
                 call dsolve( nq, omega, nq, wrk1( i,1: nq, j),2)
130           continue
!
!  Compute WRK5 = inv(E)*D*G2
              do 140 j = 1, m
                 wrk5( j) = t( i, j)
140           continue
              call dsolve( m, wrk4, m, wrk5,4)
              call dsolve( m, wrk4, m, wrk5,2)
!
!  Compute TFJACB = inv(trans(OMEGA))*FJACB
              do 170 k = 1, kp
                 do 150 l = 1, nq
                    tfjacb( i, l, k) = fjacb( i, kpvt( k), l)
150              continue
                 call dsolve( nq, omega, nq, tfjacb( i,1: nq, k),4)
                 do 160 l = 1, nq
                    if ( ss(1) .gt. zero) then
                       tfjacb( i, l, k) = tfjacb( i, l, k)/ ss( kpvt( k))
                    else
                       tfjacb( i, l, k) = tfjacb( i, l, k)/abs( ss(1))
                    endif
160              continue
170           continue
!
!  Compute WRK2 = (V*inv(E)*D**2*G2 - G1)
              do 190 l = 1, nq
                 wrk2( i, l) = zero
                 do 180 j = 1, m
                    wrk2( i, l) = wrk2( i, l)+ fjacd( i, j, l)* wrk5( j)
180              continue
                 wrk2( i, l) = wrk2( i, l)- f( i, l)
190           continue
!
!  Compute WRK2 = inv(trans(OMEGA))*(V*inv(E)*D**2*G2 - G1)
              call dsolve( nq, omega, nq, wrk2( i,1: nq),4)
300        continue
!
        else
           do 360 i = 1, n
              do 350 l = 1, nq
                 do 340 k = 1, kp
                    tfjacb( i, l, k) = fjacb( i, kpvt( k), l)
                    if ( ss(1) .gt. zero) then
                       tfjacb( i, l, k) = tfjacb( i, l, k)/ ss( kpvt( k))
                    else
                       tfjacb( i, l, k) = tfjacb( i, l, k)/abs( ss(1))
                    endif
340              continue
                 wrk2( i, l) = - f( i, l)
350           continue
360        continue
        endif
!
!  Compute S
!
!  Do QR factorization (with column pivoting of TFJACB if ALPHA = 0)
!
        if ( alpha .eq. zero) then
!-----------------------^------------------------------------------------------
!!! FPT - 3087 REAL or COMPLEX quantity tested for exact equality/inequality
!------------------------------------------------------------------------------
           ipvt = 1
           do 410 k = 1, np
              kpvt( k) = 0
410        continue
        else
           ipvt = 0
        endif
!
        call dqrdc( tfjacb, n* nq, n* nq, kp, qraux, kpvt, wrk3, ipvt)
        call dqrsl( tfjacb, n* nq, n* nq, kp,                                 &
         qraux, wrk2, dum, wrk2, dum, dum, dum,1000, inf)
        if ( inf .ne. 0) then
           istopc = 60000
           return
        endif
!
!  Eliminate alpha part using givens rotations
!
        if ( alpha .ne. zero) then
!-----------------------^------------------------------------------------------
!!! FPT - 3087 REAL or COMPLEX quantity tested for exact equality/inequality
!------------------------------------------------------------------------------
           s(1:npp) = zero
           do 430 k1 = 1, kp
              wrk3(1:kp) = zero
              wrk3( k1) = sqrt( alpha)
              do 420 k2 = k1, kp
                 call drotg( tfjacb( k2,1, k2), wrk3( k2), co, si)
                 if ( kp- k2 .ge. 1) then
                    call drot( kp- k2, tfjacb( k2,1, k2+1), n* nq,            &
                     wrk3( k2+1),1, co, si)
                 endif
                 temp = co* wrk2( k2,1)+ si* s( kpvt( k1))
                 s( kpvt( k1)) = - si* wrk2( k2,1)+ co* s( kpvt( k1))
                 wrk2( k2,1) = temp
420           continue
430        continue
        endif
!
!  Compute solution - eliminate variables if necessary
!
        if ( npp .ge. 1) then
           if ( alpha .eq. zero) then
!--------------------------^---------------------------------------------------
!!! FPT - 3087 REAL or COMPLEX quantity tested for exact equality/inequality
!------------------------------------------------------------------------------
              kp = npp
!
!  Estimate RCOND - U will contain approx null vector
!
440           call dtrco( tfjacb, n* nq, kp, rcond, u,1)
              if ( rcond .le. epsfcn) then
                 elim = .true.
                 imax = idamax( kp, u,1)
!
! IMAX is the column to remove - use DCHEX and fix KPVT
!
                 if ( imax .ne. kp) then
                    call dchex( tfjacb, n* nq, kp, imax, kp, wrk2, n* nq,1,   &
                     qraux, wrk3,2)
                    k = kpvt( imax)
                    do 450 i = imax, kp-1
                       kpvt( i) = kpvt( i+1)
450                 continue
                    kpvt( kp) = k
                 endif
                 kp = kp-1
              else
                 elim = .false.
              endif
              if ( elim .and. kp .ge. 1) then
                 goto 440
              else
                 irank = npp- kp
              endif
           endif
        endif
!
        if ( forvcv) return
!
!  Backsolve and unscramble
!
        if ( npp .ge. 1) then
           do 510 i = kp+1, npp
              wrk2( i,1) = zero
510        continue
           if ( kp .ge. 1) then
              call dtrsl( tfjacb, n* nq, kp, wrk2,01, inf)
              if ( inf .ne. 0) then
                 istopc = 60000
                 return
              endif
           endif
           do 520 i = 1, npp
              if ( ss(1) .gt. zero) then
                 s( kpvt( i)) = wrk2( i,1)/ ss( kpvt( i))
              else
                 s( kpvt( i)) = wrk2( i,1)/abs( ss(1))
              endif
520        continue
        endif
!
        if ( isodr) then
!
!  NOTE: T and WRK1 have been initialized above,
!          where T    = WD * DELTA = D*G2
!          WRK1 = trans(FJACD)*(I-FJACD*inv(P)*trans(JFACD))
!
           do 670 i = 1, n
!
!  Compute WRK4, such that
!             trans(WRK4)*WRK4 = E = (D**2 + ALPHA*TT**2)
              call desubi( n, m, wd, ldwd, ld2wd, alpha, tt, ldtt, i, wrk4)
              call dfctr(.false., wrk4, m, m, inf)
              if ( inf .ne. 0) then
                 istopc = 60000
                 return
              endif
!
!  Compute WRK5 = inv(E)*D*G2
              do 610 j = 1, m
                 wrk5( j) = t( i, j)
610           continue
              call dsolve( m, wrk4, m, wrk5,4)
              call dsolve( m, wrk4, m, wrk5,2)
!
              do 640 l = 1, nq
                 wrk2( i, l) = f( i, l)
                 do 620 k = 1, npp
                    wrk2( i, l) = wrk2( i, l)+ fjacb( i, k, l)* s( k)
620              continue
                 do 630 j = 1, m
                    wrk2( i, l) = wrk2( i, l)- fjacd( i, j, l)* wrk5( j)
630              continue
640           continue
!
              do 660 j = 1, m
                 wrk5( j) = zero
                 do 650 l = 1, nq
                    wrk5( j) = wrk5( j)+ wrk1( i, l, j)* wrk2( i, l)
650              continue
                 t( i, j) = -( wrk5( j)+ t( i, j))
660           continue
              call dsolve( m, wrk4, m, t( i,1: m),4)
              call dsolve( m, wrk4, m, t( i,1: m),2)
670        continue
!
        endif
!
!  Compute PHI(ALPHA) from scaled S and T
!
        call dwght( npp,1,reshape( ss,(/ npp,1,1/)), npp,1,                   &
         reshape( s,(/ npp,1/)), tempret(1: npp,1:1))
        wrk(1: npp) = tempret(1: npp,1)
        if ( isodr) then
           call dwght( n, m,reshape( tt,(/ ldtt,1, m/)), ldtt,1,              &
            t, tempret(1: n,1: m))
           wrk( npp+1: npp+1+ n* m-1) = reshape( tempret(1: n,1: m),(/ n* m/) &
            )
           phi = dnrm2( npp+ n* m, wrk,1)
        else
           phi = dnrm2( npp, wrk,1)
        endif
!
        return
end subroutine
!DODVCV
subroutine dodvcv                                                             &
         ( n, m, np, nq, npp,                                                 &
         f, fjacb, fjacd,                                                     &
         wd, ldwd, ld2wd, ssf, ss, tt, ldtt, delta,                           &
         epsfcn, isodr,                                                       &
         vcv, sd,                                                             &
         wrk6, omega, u, qraux, jpvt,                                         &
         s, t, irank, rcond, rss, idf, rvar, ifixb,                           &
         wrk1, wrk2, wrk3, wrk4, wrk5, wrk, lwrk, istopc)
!***Begin Prologue  DODVCV
!***Refer to  ODR
!***Routines Called  DPODI,DODSTP
!***Date Written   901207   (YYMMDD)
!***Revision Date  920619   (YYMMDD)
!***Purpose  Compute covariance matrix of estimated parameters
!***End Prologue  DODVCV
!
!...Used modules
        use odrpack_kinds,only: wp
!
!...Scalar arguments
        real(kind = wp)                                                       &
         epsfcn, rcond, rss, rvar
        integer                                                               &
         idf, irank, istopc, ldtt, ldwd, ld2wd, lwrk, m, n, np, npp, nq
        logical                                                               &
         isodr
!
!...Array arguments
        real(kind = wp)                                                       &
         delta( n, m), f( n, nq),                                             &
         fjacb( n, np, nq), fjacd( n, m, nq),                                 &
         omega( nq, nq), qraux( np), s( np), sd( np), ss( np), ssf( np),      &
         t( n, m), tt( ldtt, m), u( np), vcv( np, np), wd( ldwd, ld2wd, m),   &
         wrk1( n, nq, m), wrk2( n, nq), wrk3( np), wrk4( m, m), wrk5( m),     &
         wrk6( n* nq, np), wrk( lwrk)
        integer                                                               &
         ifixb( np), jpvt( np)
!
!...Local scalars
        real(kind = wp)                                                       &
         temp, zero
        integer                                                               &
         i, iunfix, j, junfix, kp, l
        logical                                                               &
         forvcv
!
!...External subroutines
        external                                                              &
         dpodi, dodstp
!
!...Data statements
        data                                                                  &
         zero                                                                 &
         /0.0E0_wp/
!
!...Variable definitions (alphabetically)
!       DELTA:   The estimated errors in the explanatory variables.
!       EPSFCN:  The function's precision.
!       F:       The (weighted) estimated values of EPSILON.
!       FJACB:   The Jacobian with respect to BETA.
!       FJACD:   The Jacobian with respect to DELTA.
!       FORVCV:  The variable designating whether subroutine DODSTP is
!       called to set up for the covariance matrix computations
!       (FORVCV=TRUE) or not (FORVCV=FALSE).
!       I:       An indexing variable.
!       IDF:     The degrees of freedom of the fit, equal to the number of
!       observations with nonzero weighted derivatives minus the
!       number of parameters being estimated.
!       IFIXB:   The values designating whether the elements of BETA are
!       fixed at their input values or not.
!       IMAX:    The index of the element of U having the largest absolute
!       value.
!       IRANK:   The rank deficiency of the Jacobian wrt BETA.
!       ISODR:   The variable designating whether the solution is by ODR
!       (ISODR=TRUE) or by OLS (ISODR=FALSE).
!       ISTOPC:  The variable designating whether the computations were
!       stoped due to a numerical error within subroutine DODSTP.
!       IUNFIX:  The index of the next unfixed parameter.
!       J:       An indexing variable.
!       JPVT:    The pivot vector.
!       JUNFIX:  The index of the next unfixed parameter.
!       KP:      The rank of the Jacobian wrt BETA.
!       L:       An indexing variable.
!       LDTT:    The leading dimension of array TT.
!       LDWD:    The leading dimension of array WD.
!       LD2WD:   The second dimension of array WD.
!       LWRK:    The length of vector WRK.
!       M:       The number of columns of data in the explanatory variable.
!       N:       The number of observations.
!       NP:      The number of function parameters.
!       NPP:     The number of function parameters being estimated.
!       NQ:      The number of responses per observation.
!       OMEGA:   The array defined S.T.
!       OMEGA*trans(OMEGA) = inv(I+FJACD*inv(E)*trans(FJACD))
!       = (I-FJACD*inv(P)*trans(FJACD))
!       where E = D**2 + ALPHA*TT**2
!       P = trans(FJACD)*FJACD + D**2 + ALPHA*TT**2
!       QRAUX:   The array required to recover the orthogonal part of the
!       Q-R decomposition.
!       RCOND:   The approximate reciprocal condition of FJACB.
!       RSS:     The residual sum of squares.
!       RVAR:    The residual variance.
!       S:       The step for BETA.
!       SD:      The standard deviations of the estimated BETAS.
!       SS:      The scaling values for the unfixed BETAS.
!       SSF:     The scaling values used for BETA.
!       T:       The step for DELTA.
!       TEMP:    A temporary storage location
!       TT:      The scaling values for DELTA.
!       U:       The approximate null vector for FJACB.
!       VCV:     The covariance matrix of the estimated BETAS.
!       WD:      The DELTA weights.
!       WRK:     A work array of (LWRK) elements,
!       equivalenced to WRK1 and WRK2.
!       WRK1:    A work array of (N by NQ by M) elements.
!       WRK2:    A work array of (N by NQ) elements.
!       WRK3:    A work array of (NP) elements.
!       WRK4:    A work array of (M by M) elements.
!       WRK5:    A work array of (M) elements.
!       WRK6:    A work array of (N*NQ by P) elements.
!       ZERO:    The value 0.0E0_wp.
!
!
!***First executable statement  DODVCV
!
!
        forvcv = .true.
        istopc = 0
!
        call dodstp( n, m, np, nq, npp,                                       &
         f, fjacb, fjacd,                                                     &
         wd, ldwd, ld2wd, ss, tt, ldtt, delta,                                &
         zero, epsfcn, isodr,                                                 &
         wrk6, omega, u, qraux, jpvt,                                         &
         s, t, temp, irank, rcond, forvcv,                                    &
         wrk1, wrk2, wrk3, wrk4, wrk5, wrk, lwrk, istopc)
        if ( istopc .ne. 0) then
           return
        endif
        kp = npp- irank
        call dpodi( wrk6, n* nq, kp, wrk3,1)
!
        idf = 0
        do 150 i = 1, n
           do 120 j = 1, npp
              do 110 l = 1, nq
                 if ( fjacb( i, j, l) .ne. zero) then
!------------------------------------------^-----------------------------------
!!! FPT - 3087 REAL or COMPLEX quantity tested for exact equality/inequality
!------------------------------------------------------------------------------
                    idf = idf+1
                    goto 150
                 endif
110           continue
120        continue
           if ( isodr) then
              do 140 j = 1, m
                 do 130 l = 1, nq
                    if ( fjacd( i, j, l) .ne. zero) then
!---------------------------------------------^--------------------------------
!!! FPT - 3087 REAL or COMPLEX quantity tested for exact equality/inequality
!------------------------------------------------------------------------------
                       idf = idf+1
                       goto 150
                    endif
130              continue
140           continue
           endif
150     continue
!
        if ( idf .gt. kp) then
           idf = idf- kp
           rvar = rss/ idf
        else
           idf = 0
           rvar = rss
        endif
!
!  Store variances in SD, restoring original order
!
        do 200 i = 1, np
           sd( i) = zero
200     continue
        do 210 i = 1, kp
           sd( jpvt( i)) = wrk6( i, i)
210     continue
        if ( np .gt. npp) then
           junfix = npp
           do 220 j = np,1,-1
              if ( ifixb( j) .eq. 0) then
                 sd( j) = zero
              else
                 sd( j) = sd( junfix)
                 junfix = junfix-1
              endif
220        continue
        endif
!
!  Store covariance matrix in VCV, restoring original order
!
        do 310 i = 1, np
           do 300 j = 1, i
              vcv( i, j) = zero
300        continue
310     continue
        do 330 i = 1, kp
           do 320 j = i+1, kp
              if ( jpvt( i) .gt. jpvt( j)) then
                 vcv( jpvt( i), jpvt( j)) = wrk6( i, j)
              else
                 vcv( jpvt( j), jpvt( i)) = wrk6( i, j)
              endif
320        continue
330     continue
        if ( np .gt. npp) then
           iunfix = npp
           do 360 i = np,1,-1
              if ( ifixb( i) .eq. 0) then
                 do 340 j = i,1,-1
                    vcv( i, j) = zero
340              continue
              else
                 junfix = npp
                 do 350 j = np,1,-1
                    if ( ifixb( j) .eq. 0) then
                       vcv( i, j) = zero
                    else
                       vcv( i, j) = vcv( iunfix, junfix)
                       junfix = junfix-1
                    endif
350              continue
                 iunfix = iunfix-1
              endif
360        continue
        endif
!
        do 380 i = 1, np
           vcv( i, i) = sd( i)
           sd( i) = sqrt( rvar* sd( i))
           do 370 j = 1, i
              vcv( j, i) = vcv( i, j)
370        continue
380     continue
!
!  Unscale standard errors and covariance matrix
        do 410 i = 1, np
           if ( ssf(1) .gt. zero) then
              sd( i) = sd( i)/ ssf( i)
           else
              sd( i) = sd( i)/abs( ssf(1))
           endif
           do 400 j = 1, np
              if ( ssf(1) .gt. zero) then
                 vcv( i, j) = vcv( i, j)/( ssf( i)* ssf( j))
              else
                 vcv( i, j) = vcv( i, j)/( ssf(1)* ssf(1))
              endif
400        continue
410     continue
!
        return
end subroutine
!DPACK
subroutine dpack                                                              &
         ( n2, n1, v1, v2, ifix)
!***Begin Prologue  DPACK
!***Refer to  ODR
!***Routines Called  DCOPY
!***Date Written   860529   (YYMMDD)
!***Revision Date  920304   (YYMMDD)
!***Purpose  Select the unfixed elements of V2 and return them in V1
!***End Prologue  DPACK
!
!...Used modules
        use odrpack_kinds,only: wp
!
!...Scalar arguments
        integer                                                               &
         n1, n2
!
!...Array arguments
        real(kind = wp)                                                       &
         v1( n2), v2( n2)
        integer                                                               &
         ifix( n2)
!------------^-----------------------------------------------------------------
!!! FPT - 2517 ANSI FORTRAN 77 intrinsic used as a local identifier
!------------------------------------------------------------------------------
!
!...Local scalars
        integer                                                               &
         i
!
!...External subroutines
        external                                                              &
         dcopy
!
!...Variable definitions (alphabetically)
!       I:       An indexing variable.
!       IFIX:    The values designating whether the elements of V2 are
!       fixed at their input values or not.
!       N1:      The number of items in V1.
!       N2:      The number of items in V2.
!       V1:      The vector of the unfixed items from V2.
!       V2:      The vector of the fixed and unfixed items from which the
!       unfixed elements are to be extracted.
!
!
!***First executable statement  DPACK
!
!
        n1 = 0
        if ( ifix(1) .ge. 0) then
           do 10 i = 1, n2
              if ( ifix( i) .ne. 0) then
                 n1 = n1+1
                 v1( n1) = v2( i)
              endif
10         continue
        else
           n1 = n2
           call dcopy( n2, v2,1, v1,1)
        endif
!
        return
end subroutine
!DPPNML
function dppnml                                                               &
         ( p)                                                                 &
         result( dppnmlr)
!***Begin Prologue  DPPNML
!***Refer to  ODR
!***Routines Called  (None)
!***Date Written   901207   (YYMMDD)
!***Revision Date  920304   (YYMMDD)
!***Author  Filliben, James J.,
!       Statistical Engineering Division
!       National Bureau of Standards
!       Washington, D. C. 20234
!       (Original Version--June      1972.
!       (Updated         --September 1975,
!       November  1975, AND
!       October   1976.
!***Purpose  Compute the percent point function value for the
!       normal (Gaussian) distribution with mean 0 and standard
!       deviation 1, and with probability density function
!       F(X) = (1/SQRT(2*PI))*EXP(-X*X/2).
!       (Adapted from DATAPAC subroutine TPPF, with modifications
!       to facilitate conversion to REAL (KIND=wp) automatically)
!***Description
!       --The coding as presented below is essentially
!       identical to that presented by Odeh and Evans
!       as Algortihm 70 of Applied Statistics.
!       --As pointed out by Odeh and Evans in Applied
!       Statistics, their algorithm representes a
!       substantial improvement over the previously employed
!       Hastings approximation for the normal percent point
!       function, with accuracy improving from 4.5*(10**-4)
!       to 1.5*(10**-8).
!***References  Odeh and Evans, the Percentage Points of the Normal
!       Distribution, Algortihm 70, Applied Statistics, 1974,
!       Pages 96-97.
!       Evans, Algorithms for Minimal Degree Polynomial and
!       Rational Approximation, M. Sc. Thesis, 1972,
!       University of Victoria, B. C., Canada.
!       Hastings, Approximations for Digital Computers, 1955,
!       Pages 113, 191, 192.
!       National Bureau of Standards Applied Mathematics
!       Series 55, 1964, Page 933, Formula 26.2.23.
!       Filliben, Simple and Robust Linear Estimation of the
!       Location Parameter of a Symmetric Distribution
!       (Unpublished Ph.D. Dissertation, Princeton
!       University), 1969, Pages 21-44, 229-231.
!       Filliben, "The Percent Point Function",
!       (Unpublished Manuscript), 1970, Pages 28-31.
!       Johnson and Kotz, Continuous Univariate Distributions,
!       Volume 1, 1970, Pages 40-111.
!       Kelley Statistical Tables, 1948.
!       Owen, Handbook of Statistical Tables, 1962, Pages 3-16.
!       Pearson and Hartley, Biometrika Tables for
!       Statisticians, Volume 1, 1954, Pages 104-113.
!***End Prologue  DPPNML
!
!...Used modules
        use odrpack_kinds,only: wp
!
!...Scalar arguments
        real(kind = wp)                                                       &
         p
!
!...Result
        real(kind = wp)                                                       &
         dppnmlr
!
!...Local scalars
        real(kind = wp)                                                       &
         aden, anum, half, one, p0, p1, p2, p3, p4, q0, q1, q2, q3, q4, r, t, &
         two, zero
!
!...Data statements
        data                                                                  &
         p0, p1, p2, p3, p4                                                   &
         /-0.322232431088E0_wp,-1.0E0_wp,-0.342242088547E0_wp,                &
         -0.204231210245E-1_wp,-0.453642210148E-4_wp/
        data                                                                  &
         q0, q1, q2, q3, q4                                                   &
         /0.993484626060E-1_wp,0.588581570495E0_wp,                           &
         0.531103462366E0_wp,0.103537752850E0_wp,0.38560700634E-2_wp/
        data                                                                  &
         zero, half, one, two                                                 &
         /0.0E0_wp,0.5E0_wp,1.0E0_wp,2.0E0_wp/
!
!...Variable Definitions (alphabetically)
!       ADEN:    A value used in the approximation.
!       ANUM:    A value used in the approximation.
!       HALF:    The value 0.5E0_wp.
!       ONE:     The value 1.0E0_wp.
!       P:       The probability at which the percent point is to be
!       evaluated.  P must be between 0.0E0_wp and 1.0E0_wp, exclusive.
!       P0:      A parameter used in the approximation.
!       P1:      A parameter used in the approximation.
!       P2:      A parameter used in the approximation.
!       P3:      A parameter used in the approximation.
!       P4:      A parameter used in the approximation.
!       Q0:      A parameter used in the approximation.
!       Q1:      A parameter used in the approximation.
!       Q2:      A parameter used in the approximation.
!       Q3:      A parameter used in the approximation.
!       Q4:      A parameter used in the approximation.
!       R:       The probability at which the percent point is evaluated.
!       T:       A value used in the approximation.
!       TWO:     The value 2.0E0_wp.
!       ZERO:    The value 0.0E0_wp.
!
!
!***First executable statement  DPPT
!
!
        if ( p .eq. half) then
!-------------------^----------------------------------------------------------
!!! FPT - 3087 REAL or COMPLEX quantity tested for exact equality/inequality
!------------------------------------------------------------------------------
           dppnmlr = zero
!
        else
           r = p
           if ( p .gt. half) r = one- r
           t = sqrt(- two*log( r))
           anum = (((( t* p4+ p3)* t+ p2)* t+ p1)* t+ p0)
           aden = (((( t* q4+ q3)* t+ q2)* t+ q1)* t+ q0)
           dppnmlr = t+( anum/ aden)
!
           if ( p .lt. half) dppnmlr = - dppnmlr
        endif
!
        return
!
end function
!DPPT
function dppt                                                                 &
         ( p, idf)                                                            &
         result( dpptr)
!***Begin Prologue  DPPT
!***Refer to  ODR
!***Routines Called  DPPNML
!***Date Written   901207   (YYMMDD)
!***Revision Date  920304   (YYMMDD)
!***Author  Filliben, James J.,
!       Statistical Engineering Division
!       National Bureau of Standards
!       Washington, D. C. 20234
!       (Original Version--October   1975.)
!       (Updated         --November  1975.)
!***Purpose  Compute the percent point function value for the
!       student's T distribution with IDF degrees of freedom.
!       (Adapted from DATAPAC subroutine TPPF, with modifications
!       to facilitate conversion to REAL (KIND=wp) automatically)
!***Description
!       --For IDF = 1 AND IDF = 2, the percent point function
!       for the T distribution exists in simple closed form
!       and so the computed percent points are exact.
!       --For IDF between 3 and 6, inclusively, the approximation
!       is augmented by 3 iterations of Newton's method to
!       improve the accuracy, especially for P near 0 or 1.
!***References  National Bureau of Standards Applied Mathmatics
!       Series 55, 1964, Page 949, Formula 26.7.5.
!       Johnson and Kotz, Continuous Univariate Distributions,
!       Volume 2, 1970, Page 102, Formula 11.
!       Federighi, "Extended Tables of the Percentage Points
!       of Student"S T Distribution, Journal of the American
!       Statistical Association, 1969, Pages 683-688.
!       Hastings and Peacock, Statistical Distributions, A
!       Handbook for Students and Practitioners, 1975,
!       Pages 120-123.
!***End Prologue  DPPT
!
!...Used modules
        use odrpack_kinds,only: wp
!
!...Scalar arguments
        real(kind = wp)                                                       &
         p
        integer                                                               &
         idf
!
!...Result
        real(kind = wp)                                                       &
         dpptr
!
!...Local scalars
        real(kind = wp)                                                       &
         arg, b21, b31, b32, b33, b34, b41, b42, b43, b44, b45,               &
         b51, b52, b53, b54, b55, b56, c, con, d1, d3, d5, d7, d9, df, eight, &
         fiftn, half, one, pi, ppfn, s, term1, term2, term3, term4, term5,    &
         three, two, z, zero
        integer                                                               &
         ipass, maxit
!
!...External functions
        real(kind = wp)                                                       &
         dppnml
        external                                                              &
         dppnml
!
!...Data statements
        data                                                                  &
         b21                                                                  &
         /4.0E0_wp/
        data                                                                  &
         b31, b32, b33, b34                                                   &
         /96.0E0_wp,5.0E0_wp,16.0E0_wp,3.0E0_wp/
        data                                                                  &
         b41, b42, b43, b44, b45                                              &
         /384.0E0_wp,3.0E0_wp,19.0E0_wp,17.0E0_wp,-15.0E0_wp/
        data                                                                  &
         b51, b52, b53, b54, b55, b56                                         &
         /9216.0E0_wp,79.0E0_wp,776.0E0_wp,1482.0E0_wp,-1920.0E0_wp,          &
         -945.0E0_wp/
        data                                                                  &
         zero, half, one, two, three, eight, fiftn                            &
         /0.0E0_wp,0.5E0_wp,1.0E0_wp,2.0E0_wp,3.0E0_wp,8.0E0_wp,              &
         15.0E0_wp/
!
!...Variable definitions (alphabetically)
!       ARG:    A value used in the approximation.
!       B21:    A parameter used in the approximation.
!       B31:    A parameter used in the approximation.
!       B32:    A parameter used in the approximation.
!       B33:    A parameter used in the approximation.
!       B34:    A parameter used in the approximation.
!       B41:    A parameter used in the approximation.
!       B42:    A parameter used in the approximation.
!       B43:    A parameter used in the approximation.
!       B44:    A parameter used in the approximation.
!       B45:    A parameter used in the approximation.
!       B51:    A parameter used in the approximation.
!       B52:    A parameter used in the approximation.
!       B53:    A parameter used in the approximation.
!       B54:    A parameter used in the approximation.
!       B55:    A parameter used in the approximation.
!       B56:    A parameter used in the approximation.
!       C:      A value used in the approximation.
!       CON:    A value used in the approximation.
!       DF:     The degrees of freedom.
!       D1:     A value used in the approximation.
!       D3:     A value used in the approximation.
!       D5:     A value used in the approximation.
!       D7:     A value used in the approximation.
!       D9:     A value used in the approximation.
!       EIGHT:  The value 8.0E0_wp.
!       FIFTN:  The value 15.0E0_wp.
!       HALF:   The value 0.5E0_wp.
!       IDF:    The (positive integer) degrees of freedom.
!       IPASS:  A value used in the approximation.
!       MAXIT:  The maximum number of iterations allowed for the approx.
!       ONE:    The value 1.0E0_wp.
!       P:      The probability at which the percent point is to be
!       evaluated.  P must lie between 0.0DO and 1.0E0_wp, exclusive.
!       PI:     The value of pi.
!       PPFN:   The normal percent point value.
!       S:      A value used in the approximation.
!       TERM1:  A value used in the approximation.
!       TERM2:  A value used in the approximation.
!       TERM3:  A value used in the approximation.
!       TERM4:  A value used in the approximation.
!       TERM5:  A value used in the approximation.
!       THREE:  The value 3.0E0_wp.
!       TWO:    The value 2.0E0_wp.
!       Z:      A value used in the approximation.
!       ZERO:   The value 0.0E0_wp.
!
!
!***First executable statement  DPPT
!
!
        pi = 3.141592653589793238462643383279E0_wp
        df = idf
        maxit = 5
!
        if ( idf .le. 0) then
!
!  Treat the IDF < 1 case
           dpptr = zero
!
        elseif ( idf .eq. 1) then
!
!  Treat the IDF = 1 (Cauchy) case
           arg = pi* p
           dpptr = -cos( arg)/sin( arg)
!
        elseif ( idf .eq. 2) then
!
!  Treat the IDF = 2 case
           term1 = sqrt( two)/ two
           term2 = two* p- one
           term3 = sqrt( p*( one- p))
           dpptr = term1* term2/ term3
!
        elseif ( idf .ge. 3) then
!
!  Treat the IDF greater than or equal to 3 case
           ppfn = dppnml( p)
           d1 = ppfn
           d3 = ppfn**3
           d5 = ppfn**5
           d7 = ppfn**7
           d9 = ppfn**9
           term1 = d1
           term2 = ( one/ b21)*( d3+ d1)/ df
           term3 = ( one/ b31)*( b32* d5+ b33* d3+ b34* d1)/( df**2)
           term4 = ( one/ b41)*( b42* d7+ b43* d5+ b44* d3+ b45* d1)/( df**3)
           term5 = ( one/ b51)*( b52* d9+ b53* d7+ b54* d5+ b55* d3+ b56* d1) &
            /( df**4)
           dpptr = term1+ term2+ term3+ term4+ term5
!
           if ( idf .eq. 3) then
!
!  Augment the results for the IDF = 3 case
              con = pi*( p- half)
              arg = dpptr/sqrt( df)
              z = atan( arg)
              do 70 ipass = 1, maxit
                 s = sin( z)
                 c = cos( z)
                 z = z-( z+ s* c- con)/( two* c**2)
70            continue
              dpptr = sqrt( df)* s/ c
!
           elseif ( idf .eq. 4) then
!
!  Augment the results for the IDF = 4 case
              con = two*( p- half)
              arg = dpptr/sqrt( df)
              z = atan( arg)
              do 90 ipass = 1, maxit
                 s = sin( z)
                 c = cos( z)
                 z = z-(( one+ half* c**2)* s- con)/(( one+ half)* c**3)
90            continue
              dpptr = sqrt( df)* s/ c
!
           elseif ( idf .eq. 5) then
!
!  Augment the results for the IDF = 5 case
!
              con = pi*( p- half)
              arg = dpptr/sqrt( df)
              z = atan( arg)
              do 110 ipass = 1, maxit
                 s = sin( z)
                 c = cos( z)
                 z = z-( z+( c+( two/ three)* c**3)* s- con)/                 &
                  (( eight/ three)* c**4)
110           continue
              dpptr = sqrt( df)* s/ c
!
           elseif ( idf .eq. 6) then
!
!  Augment the results for the IDF = 6 case
              con = two*( p- half)
              arg = dpptr/sqrt( df)
              z = atan( arg)
              do 130 ipass = 1, maxit
                 s = sin( z)
                 c = cos( z)
                 z = z-(( one+ half* c**2+( three/ eight)* c**4)* s- con)/    &
                  (( fiftn/ eight)* c**5)
130           continue
              dpptr = sqrt( df)* s/ c
           endif
        endif
!
        return
!
end function
!DPVB
subroutine dpvb                                                               &
         ( fcn,                                                               &
         n, m, np, nq,                                                        &
         beta, xplusd, ifixb, ifixx, ldifx,                                   &
         nrow, j, lq, stp,                                                    &
         istop, nfev, pvb,                                                    &
         wrk1, wrk2, wrk6)
!***Begin Prologue  DPVB
!***Refer to  ODR
!***Routines Called  FCN
!***Date Written   860529   (YYMMDD)
!***Revision Date  920304   (YYMMDD)
!***Purpose  Compute the NROW-th function value using BETA(J) + STP
!***End Prologue  DPVB
!
!...Used modules
        use odrpack_kinds,only: wp
!
!...Scalar arguments
        real(kind = wp)                                                       &
         pvb, stp
        integer                                                               &
         istop, j, ldifx, lq, m, n, nfev, np, nq, nrow
!
!...Array arguments
        real(kind = wp)                                                       &
         beta( np), wrk1( n, m, nq), wrk2( n, nq), wrk6( n, np, nq), xplusd(  &
         n, m)
        integer                                                               &
         ifixb( np), ifixx( ldifx, m)
!
!...Subroutine arguments
        external                                                              &
         fcn
!
!...Local scalars
        real(kind = wp)                                                       &
         betaj
!
!...Routine names used as subprogram arguments
!       FCN:     The user-supplied subroutine for evaluating the model.
!
!...Variable Definitions (alphabetically)
!       BETA:    The function parameters.
!       BETAJ:   The current estimate of the jth parameter.
!       IFIXB:   The values designating whether the elements of BETA are
!       fixed at their input values or not.
!       IFIXX:   The values designating whether the elements of X are
!       fixed at their input values or not.
!       ISTOP:   The variable designating whether there are problems
!       computing the function at the current BETA and DELTA.
!       J:       The index of the partial derivative being examined.
!       LDIFX:   The leading dimension of array IFIXX.
!       LQ:      The response currently being examined.
!       M:       The number of columns of data in the independent variable.
!       N:       The number of observations.
!       NFEV:    The number of function evaluations.
!       NP:      The number of function parameters.
!       NQ:      The number of responses per observation.
!       NROW:    The row number of the independent variable array at
!       which the derivative is to be checked.
!       PVB:     The function value for the selected observation & response.
!       STP:     The step size for the finite difference derivative.
!       XPLUSD:  The values of X + DELTA.
!
!
!***First executable statement  DPVB
!
!
!  Compute predicted values
!
        betaj = beta( j)
        beta( j) = beta( j)+ stp
        istop = 0
        call fcn( n, m, np, nq,                                               &
         n, m, np,                                                            &
         beta, xplusd,                                                        &
         ifixb, ifixx, ldifx,                                                 &
         003, wrk2, wrk6, wrk1,                                               &
         istop)
        if ( istop .eq. 0) then
           nfev = nfev+1
        else
           return
        endif
        beta( j) = betaj
!
        pvb = wrk2( nrow, lq)
!
        return
end subroutine
!DPVD
subroutine dpvd                                                               &
         ( fcn,                                                               &
         n, m, np, nq,                                                        &
         beta, xplusd, ifixb, ifixx, ldifx,                                   &
         nrow, j, lq, stp,                                                    &
         istop, nfev, pvd,                                                    &
         wrk1, wrk2, wrk6)
!***Begin Prologue  DPVD
!***Refer to  ODR
!***Routines Called  FCN
!***Date Written   860529   (YYMMDD)
!***Revision Date  920304   (YYMMDD)
!***Purpose  Compute NROW-th function value using
!       X(NROW,J) + DELTA(NROW,J) + STP
!***End Prologue  DPVD
!
!...Used modules
        use odrpack_kinds,only: wp
!
!...Scalar arguments
        real(kind = wp)                                                       &
         pvd, stp
        integer                                                               &
         istop, j, ldifx, lq, m, n, nfev, np, nq, nrow
!
!...Array arguments
        real(kind = wp)                                                       &
         beta( np), wrk1( n, m, nq), wrk2( n, nq), wrk6( n, np, nq), xplusd(  &
         n, m)
        integer                                                               &
         ifixb( np), ifixx( ldifx, m)
!
!...Subroutine arguments
        external                                                              &
         fcn
!
!...Local scalars
        real(kind = wp)                                                       &
         xpdj
!
!...Routine names used as subprogram arguments
!       FCN:     The user-supplied subroutine for evaluating the model.
!
!...Variable Definitions (alphabetically)
!       BETA:    The function parameters.
!       IFIXB:   The values designating whether the elements of BETA are
!       fixed at their input values or not.
!       IFIXX:   The values designating whether the elements of X are
!       fixed at their input values or not.
!       ISTOP:   The variable designating whether there are problems
!       computing the function at the current BETA and DELTA.
!       J:       The index of the partial derivative being examined.
!       LDIFX:   The leading dimension of array IFIXX.
!       LQ:      The response currently being examined.
!       M:       The number of columns of data in the independent variable.
!       N:       The number of observations.
!       NFEV:    The number of function evaluations.
!       NP:      The number of function parameters.
!       NQ:      The number of responses per observation.
!       NROW:    The row number of the independent variable array at
!       which the derivative is to be checked.
!       PVD:     The function value for the selected observation & response.
!       STP:     The step size for the finite difference derivative.
!       XPDJ:    The (NROW,J)th element of XPLUSD.
!       XPLUSD:  The values of X + DELTA.
!
!
!***First executable statement  DPVD
!
!
!  Compute predicted values
!
        xpdj = xplusd( nrow, j)
        xplusd( nrow, j) = xplusd( nrow, j)+ stp
        istop = 0
        call fcn( n, m, np, nq,                                               &
         n, m, np,                                                            &
         beta, xplusd,                                                        &
         ifixb, ifixx, ldifx,                                                 &
         003, wrk2, wrk6, wrk1,                                               &
         istop)
        if ( istop .eq. 0) then
           nfev = nfev+1
        else
           return
        endif
        xplusd( nrow, j) = xpdj
!
        pvd = wrk2( nrow, lq)
!
        return
end subroutine
!DSCALE
subroutine dscale                                                             &
         ( n, m, scl, ldscl, t, ldt, sclt, ldsclt)
!***Begin Prologue  DSCALE
!***Refer to  ODR
!***Routines Called  (NONE)
!***Date Written   860529   (YYMMDD)
!***Revision Date  920304   (YYMMDD)
!***Purpose  Scale T by the inverse of SCL, I.E., compute T/SCL
!***End Prologue  DSCALE
!
!...Used modules
        use odrpack_kinds,only: wp
!
!...Scalar arguments
        integer                                                               &
         ldt, ldscl, ldsclt, m, n
!
!...Array arguments
        real(kind = wp)                                                       &
         t( ldt, m), scl( ldscl, m), sclt( ldsclt, m)
!
!...Local scalars
        real(kind = wp)                                                       &
         one, temp, zero
        integer                                                               &
         i, j
!
!...Data statements
        data                                                                  &
         one, zero                                                            &
         /1.0E0_wp,0.0E0_wp/
!
!...Variable Definitions (alphabetically)
!       I:       An indexing variable.
!       J:       An indexing variable.
!       LDSCL:   The leading dimension of array SCL.
!       LDSCLT:  The leading dimension of array SCLT.
!       LDT:     The leading dimension of array T.
!       M:       The number of columns of data in T.
!       N:       The number of rows of data in T.
!       ONE:     The value 1.0E0_wp.
!       SCL:     The scale values.
!       SCLT:    The inversely scaled matrix.
!       T:       The array to be inversely scaled by SCL.
!       TEMP:    A temporary scalar.
!       ZERO:    The value 0.0E0_wp.
!
!
!***First executable statement  DSCALE
!
!
        if ( n .eq. 0 .or. m .eq. 0) return
!
        if ( scl(1,1) .ge. zero) then
           if ( ldscl .ge. n) then
              do 80 j = 1, m
                 do 70 i = 1, n
                    sclt( i, j) = t( i, j)/ scl( i, j)
70               continue
80            continue
           else
              do 100 j = 1, m
                 temp = one/ scl(1, j)
                 do 90 i = 1, n
                    sclt( i, j) = t( i, j)* temp
90               continue
100           continue
           endif
        else
           temp = one/abs( scl(1,1))
           do 120 j = 1, m
              do 110 i = 1, n
                 sclt( i, j) = t( i, j)* temp
110           continue
120        continue
        endif
!
        return
end subroutine
!DSCLB
subroutine dsclb                                                              &
         ( np, beta, ssf)
!***Begin Prologue  DSCLB
!***Refer to  ODR
!***Routines Called  (NONE)
!***Date Written   860529   (YYMMDD)
!***Revision Date  920304   (YYMMDD)
!***Purpose  Select scaling values for BETA according to the
!       algorithm given in the ODRPACK95 reference guide
!***End Prologue  DSCLB
!
!...Used modules
        use odrpack_kinds,only: wp
!
!...Scalar arguments
        integer                                                               &
         np
!
!...Array arguments
        real(kind = wp)                                                       &
         beta( np), ssf( np)
!
!...Local scalars
        real(kind = wp)                                                       &
         bmax, bmin, one, ten, zero
        integer                                                               &
         k
        logical                                                               &
         bigdif
!
!...Data statements
        data                                                                  &
         zero, one, ten                                                       &
         /0.0E0_wp,1.0E0_wp,10.0E0_wp/
!
!...Variable Definitions (alphabetically)
!       BETA:    The function parameters.
!       BIGDIF:  The variable designating whether there is a significant
!       difference in the magnitudes of the nonzero elements of
!       BETA (BIGDIF=.TRUE.) or not (BIGDIF=.FALSE.).
!       BMAX:    The largest nonzero magnitude.
!       BMIN:    The smallest nonzero magnitude.
!       K:       An indexing variable.
!       NP:      The number of function parameters.
!       ONE:     The value 1.0E0_wp.
!       SSF:     The scaling values for BETA.
!       TEN:     The value 10.0E0_wp.
!       ZERO:    The value 0.0E0_wp.
!
!
!***First executable statement  DSCLB
!
!
        bmax = abs( beta(1))
        do 10 k = 2, np
           bmax = max( bmax,abs( beta( k)))
10      continue
!
        if ( bmax .eq. zero) then
!----------------------^-------------------------------------------------------
!!! FPT - 3087 REAL or COMPLEX quantity tested for exact equality/inequality
!------------------------------------------------------------------------------
!
!  All input values of BETA are zero
!
           do 20 k = 1, np
              ssf( k) = one
20         continue
!
        else
!
!  Some of the input values are nonzero
!
           bmin = bmax
           do 30 k = 1, np
              if ( beta( k) .ne. zero) then
!--------------------------------^---------------------------------------------
!!! FPT - 3087 REAL or COMPLEX quantity tested for exact equality/inequality
!------------------------------------------------------------------------------
                 bmin = min( bmin,abs( beta( k)))
              endif
30         continue
           bigdif = log10( bmax)-log10( bmin) .ge. one
           do 40 k = 1, np
              if ( beta( k) .eq. zero) then
!--------------------------------^---------------------------------------------
!!! FPT - 3087 REAL or COMPLEX quantity tested for exact equality/inequality
!------------------------------------------------------------------------------
                 ssf( k) = ten/ bmin
              else
                 if ( bigdif) then
                    ssf( k) = one/abs( beta( k))
                 else
                    ssf( k) = one/ bmax
                 endif
              endif
40         continue
!
        endif
!
        return
end subroutine
!DSCLD
subroutine dscld                                                              &
         ( n, m, x, ldx, tt, ldtt)
!***Begin Prologue  DSCLD
!***Refer to  ODR
!***Routines Called  (None)
!***Date Written   860529   (YYMMDD)
!***Revision Date  920304   (YYMMDD)
!***Purpose  Select scaling values for DELTA according to the
!       algorithm given in the ODRPACK95 reference guide
!***End Prologue  DSCLD
!
!...Used modules
        use odrpack_kinds,only: wp
!
!...Scalar arguments
        integer                                                               &
         ldtt, ldx, m, n
!
!...Array arguments
        real(kind = wp)                                                       &
         tt( ldtt, m), x( ldx, m)
!
!...Local scalars
        real(kind = wp)                                                       &
         one, ten, xmax, xmin, zero
        integer                                                               &
         i, j
        logical                                                               &
         bigdif
!
!...Data statements
        data                                                                  &
         zero, one, ten                                                       &
         /0.0E0_wp,1.0E0_wp,10.0E0_wp/
!
!...Variable Definitions (alphabetically)
!       BIGDIF:  The variable designating whether there is a significant
!       difference in the magnitudes of the nonzero elements of
!       X (BIGDIF=.TRUE.) or not (BIGDIF=.FALSE.).
!       I:       An indexing variable.
!       J:       An indexing variable.
!       LDTT:    The leading dimension of array TT.
!       LDX:     The leading dimension of array X.
!       M:       The number of columns of data in the independent variable.
!       N:       The number of observations.
!       ONE:     The value 1.0E0_wp.
!       TT:      THE SCALING VALUES FOR DELTA.
!       X:       The independent variable.
!       XMAX:    The largest nonzero magnitude.
!       XMIN:    THE SMALLEST NONZERO MAGNITUDE.
!       ZERO:    The value 0.0E0_wp.
!
!
!***First executable statement  DSCLD
!
!
        do 50 j = 1, m
           xmax = abs( x(1, j))
           do 10 i = 2, n
              xmax = max( xmax,abs( x( i, j)))
10         continue
!
           if ( xmax .eq. zero) then
!-------------------------^----------------------------------------------------
!!! FPT - 3087 REAL or COMPLEX quantity tested for exact equality/inequality
!------------------------------------------------------------------------------
!
!  All input values of X(I,J), I=1,...,N, are zero
!
              do 20 i = 1, n
                 tt( i, j) = one
20            continue
!
           else
!
!  Some of the input values are nonzero
!
              xmin = xmax
              do 30 i = 1, n
                 if ( x( i, j) .ne. zero) then
!-----------------------------------^------------------------------------------
!!! FPT - 3087 REAL or COMPLEX quantity tested for exact equality/inequality
!------------------------------------------------------------------------------
                    xmin = min( xmin,abs( x( i, j)))
                 endif
30            continue
              bigdif = log10( xmax)-log10( xmin) .ge. one
              do 40 i = 1, n
                 if ( x( i, j) .ne. zero) then
!-----------------------------------^------------------------------------------
!!! FPT - 3087 REAL or COMPLEX quantity tested for exact equality/inequality
!------------------------------------------------------------------------------
                    if ( bigdif) then
                       tt( i, j) = one/abs( x( i, j))
                    else
                       tt( i, j) = one/ xmax
                    endif
                 else
                    tt( i, j) = ten/ xmin
                 endif
40            continue
           endif
50      continue
!
        return
end subroutine
!DSETN
subroutine dsetn                                                              &
         ( n, m, x, ldx, nrow)
!***Begin Prologue  DSETN
!***Refer to  ODR
!***Routines Called  (None)
!***Date Written   860529   (YYMMDD)
!***Revision Date  920304   (YYMMDD)
!***Purpose  Select the row at which the derivative will be checked
!***End Prologue  DSETN
!
!...Used modules
        use odrpack_kinds,only: wp
!
!...Scalar arguments
        integer                                                               &
         ldx, m, n, nrow
!
!...Array arguments
        real(kind = wp)                                                       &
         x( ldx, m)
!
!...Local scalars
        integer                                                               &
         i, j
!
!...Variable Definitions (alphabetically)
!       I:       An index variable.
!       J:       An index variable.
!       LDX:     The leading dimension of array X.
!       M:       The number of columns of data in the independent variable.
!       N:       The number of observations.
!       NROW:    The selected row number of the independent variable.
!       X:       The independent variable.
!
!
!***First executable statement  DSETN
!
!
        if (( nrow .ge. 1) .and.                                              &
            ( nrow .le. n)) return
!
!       Select first row of independent variables which contains no zeros
!       if there is one, otherwise first row is used.
!
        do 20 i = 1, n
           do 10 j = 1, m
              if ( x( i, j) .eq. 0.0) goto 20
!--------------------------------^---------------------------------------------
!!! FPT - 3087 REAL or COMPLEX quantity tested for exact equality/inequality
!!! FPT - 3089 Comparison of REAL or COMPLEX values which differ in precision
!------------------------------------------------------------------------------
10         continue
           nrow = i
           return
20      continue
!
        nrow = 1
!
        return
end subroutine
!DSOLVE
subroutine dsolve( n, t, ldt, b, job)
!***Begin Prologue  DSOLVE
!***Refer to  ODR
!***Routines Called  DAXPY,DDOT
!***Date Written   920220   (YYMMDD)
!***Revision Date  920619   (YYMMDD)
!***Purpose  Solve systems of the form
!       T * X = B  or  trans(T) * X = B
!       where T is an upper or lower triangular matrix of order N,
!       and the solution X overwrites the RHS B.
!       (adapted from LINPACK subroutine DTRSL)
!***References  Dongarra J.J., Bunch J.R., Moler C.B., Stewart G.W.,
!       *LINPACK Users Guide*, SIAM, 1979.
!***End Prologue  DSOLVE
!
!...Used modules
        use odrpack_kinds,only: wp
!
!...Scalar arguments
        integer                                                               &
         job, ldt, n
!
!...Array arguments
        real(kind = wp)                                                       &
         b( n), t( ldt, n)
!
!...Local scalars
        real(kind = wp)                                                       &
         temp, zero
        integer                                                               &
         j1, j, jn
!
!...External functions
        real(kind = wp)                                                       &
         ddot
        external                                                              &
         ddot
!
!...External subroutines
        external                                                              &
         daxpy
!
!...Data statements
        data                                                                  &
         zero                                                                 &
         /0.0E0_wp/
!
!...Variable Definitions (alphabetically)
!       B:       On input:  the right hand side;  On exit:  the solution
!       J1:      The first nonzero entry in T.
!       J:       An indexing variable.
!       JN:      The last nonzero entry in T.
!       JOB:     What kind of system is to be solved, where if JOB is
!       1   Solve T*X=B, T lower triangular,
!       2   Solve T*X=B, T upper triangular,
!       3   Solve trans(T)*X=B, T lower triangular,
!       4   Solve trans(T)*X=B, T upper triangular.
!       LDT:     The leading dimension of array T.
!       N:       The number of rows and columns of data in array T.
!       T:       The upper or lower tridiagonal system.
!       ZERO:    The value 0.0E0_wp.
!
!
!***First executable statement  DSOLVE
!
!
!  Find first nonzero diagonal entry in T
        j1 = 0
        do 10 j = 1, n
           if ( j1 .eq. 0 .and. t( j, j) .ne. zero) then
!---------------------------------------------^--------------------------------
!!! FPT - 3087 REAL or COMPLEX quantity tested for exact equality/inequality
!------------------------------------------------------------------------------
              j1 = j
           elseif ( t( j, j) .eq. zero) then
!---------------------------------^--------------------------------------------
!!! FPT - 3087 REAL or COMPLEX quantity tested for exact equality/inequality
!------------------------------------------------------------------------------
              b( j) = zero
           endif
10      continue
        if ( j1 .eq. 0) return
!
!  Find last nonzero diagonal entry in T
        jn = 0
        do 20 j = n, j1,-1
           if ( jn .eq. 0 .and. t( j, j) .ne. zero) then
!---------------------------------------------^--------------------------------
!!! FPT - 3087 REAL or COMPLEX quantity tested for exact equality/inequality
!------------------------------------------------------------------------------
              jn = j
           elseif ( t( j, j) .eq. zero) then
!---------------------------------^--------------------------------------------
!!! FPT - 3087 REAL or COMPLEX quantity tested for exact equality/inequality
!------------------------------------------------------------------------------
              b( j) = zero
           endif
20      continue
!
        if ( job .eq. 1) then
!
!  Solve T*X=B for T lower triangular
           b( j1) = b( j1)/ t( j1, j1)
           do 30 j = j1+1, jn
              temp = - b( j-1)
              call daxpy( jn- j+1, temp, t( j, j-1),1, b( j),1)
              if ( t( j, j) .ne. zero) then
!--------------------------------^---------------------------------------------
!!! FPT - 3087 REAL or COMPLEX quantity tested for exact equality/inequality
!------------------------------------------------------------------------------
                 b( j) = b( j)/ t( j, j)
              else
                 b( j) = zero
              endif
30         continue
!
        elseif ( job .eq. 2) then
!
!  Solve T*X=B for T upper triangular.
           b( jn) = b( jn)/ t( jn, jn)
           do 40 j = jn-1, j1,-1
              temp = - b( j+1)
              call daxpy( j, temp, t(1, j+1),1, b(1),1)
              if ( t( j, j) .ne. zero) then
!--------------------------------^---------------------------------------------
!!! FPT - 3087 REAL or COMPLEX quantity tested for exact equality/inequality
!------------------------------------------------------------------------------
                 b( j) = b( j)/ t( j, j)
              else
                 b( j) = zero
              endif
40         continue
!
        elseif ( job .eq. 3) then
!
!  Solve trans(T)*X=B for T lower triangular.
           b( jn) = b( jn)/ t( jn, jn)
           do 50 j = jn-1, j1,-1
              b( j) = b( j)- ddot( jn- j+1, t( j+1, j),1, b( j+1),1)
              if ( t( j, j) .ne. zero) then
!--------------------------------^---------------------------------------------
!!! FPT - 3087 REAL or COMPLEX quantity tested for exact equality/inequality
!------------------------------------------------------------------------------
                 b( j) = b( j)/ t( j, j)
              else
                 b( j) = zero
              endif
50         continue
!
        elseif ( job .eq. 4) then
!
!  Solve trans(T)*X=B for T upper triangular.
           b( j1) = b( j1)/ t( j1, j1)
           do 60 j = j1+1, jn
              b( j) = b( j)- ddot( j-1, t(1, j),1, b(1),1)
              if ( t( j, j) .ne. zero) then
!--------------------------------^---------------------------------------------
!!! FPT - 3087 REAL or COMPLEX quantity tested for exact equality/inequality
!------------------------------------------------------------------------------
                 b( j) = b( j)/ t( j, j)
              else
                 b( j) = zero
              endif
60         continue
        endif
!
        return
end subroutine
!DUNPAC
subroutine dunpac                                                             &
         ( n2, v1, v2, ifix)
!***Begin Prologue  DUNPAC
!***Refer to  ODR
!***Routines Called  DCOPY
!***Date Written   860529   (YYMMDD)
!***Revision Date  920304   (YYMMDD)
!***Purpose  Copy the elements of V1 into the locations of V2 which are
!       unfixed
!***End Prologue  DUNPAC
!
!...Used modules
        use odrpack_kinds,only: wp
!
!...Scalar arguments
        integer                                                               &
         n2
!
!...Array arguments
        real(kind = wp)                                                       &
         v1( n2), v2( n2)
        integer                                                               &
         ifix( n2)
!------------^-----------------------------------------------------------------
!!! FPT - 2517 ANSI FORTRAN 77 intrinsic used as a local identifier
!------------------------------------------------------------------------------
!
!...Local scalars
        integer                                                               &
         i, n1
!
!...External subroutines
        external                                                              &
         dcopy
!
!...Variable Definitions (alphabetically)
!       I:       An indexing variable.
!       IFIX:    The values designating whether the elements of V2 are
!       fixed at their input values or not.
!       ODRPACK95 reference guide.)
!       N1:      The number of items in V1.
!       N2:      The number of items in V2.
!       V1:      The vector of the unfixed items.
!       V2:      The vector of the fixed and unfixed items into which the
!       elements of V1 are to be inserted.
!
!
!***First executable statement  DUNPAC
!
!
        n1 = 0
        if ( ifix(1) .ge. 0) then
           do 10 i = 1, n2
              if ( ifix( i) .ne. 0) then
                 n1 = n1+1
                 v2( i) = v1( n1)
              endif
10         continue
        else
           n1 = n2
           call dcopy( n2, v1,1, v2,1)
        endif
!
        return
end subroutine
!DVEVTR
subroutine dvevtr                                                             &
         ( m, nq, indx,                                                       &
         v, ldv, ld2v, e, lde, ve, ldve, ld2ve, vev, ldvev,                   &
         wrk5)
!***Begin Prologue  DVEVTR
!***Refer to  ODR
!***Routines Called  DSOLVE
!***Date Written   910613   (YYMMDD)
!***Revision Date  920304   (YYMMDD)
!***Purpose  Compute  V*E*trans(V) for the (INDX)TH M by NQ array in V
!***End Prologue  DVEVTR
!
!...Used modules
        use odrpack_kinds,only: wp
!
!...Scalar arguments
        integer                                                               &
         indx, lde, ldv, ldve, ldvev, ld2v, ld2ve, m, nq
!
!...Array arguments
        real(kind = wp)                                                       &
         e( lde, m), v( ldv, ld2v, nq), ve( ldve, ld2ve, m), vev( ldvev, nq), &
         wrk5( m)
!
!...Local scalars
        real(kind = wp)                                                       &
         zero
        integer                                                               &
         j, l1, l2
!
!...External subroutines
        external                                                              &
         dsolve
!
!...Data statements
        data                                                                  &
         zero                                                                 &
         /0.0E0_wp/
!
!...Variable Definitions (alphabetically)
!       INDX:    The row in V in which the M by NQ array is stored.
!       J:       An indexing variable.
!       LDE:     The leading dimension of array E.
!       LDV:     The leading dimension of array V.
!       LDVE:    The leading dimension of array VE.
!       LDVEV:   The leading dimension of array VEV.
!       LD2V:    The second dimension of array V.
!       L1:      An indexing variable.
!       L2:      An indexing variable.
!       M:       The number of columns of data in the independent variable.
!       NQ:      The number of responses per observation.
!       E:       The M by M matrix of the factors so ETE = (D**2 + ALPHA*T**2).
!       V:       An array of NQ by M matrices.
!       VE:      The NQ by M array VE = V * inv(E)
!       VEV:     The NQ by NQ array VEV = V * inv(ETE) * trans(V).
!       WRK5:    An M work vector.
!       ZERO:    The value 0.0E0_wp.
!
!
!***First executable statement  DVEVTR
!
!
        if ( nq .eq. 0 .or. m .eq. 0) return
!
        do 140 l1 = 1, nq
           do 110 j = 1, m
              wrk5( j) = v( indx, j, l1)
110        continue
           call dsolve( m, e, lde, wrk5,4)
           do 120 j = 1, m
              ve( indx, l1, j) = wrk5( j)
120        continue
140     continue
!
        do 230 l1 = 1, nq
           do 220 l2 = 1, l1
              vev( l1, l2) = zero
              do 210 j = 1, m
                 vev( l1, l2) = vev( l1, l2)+ ve( indx, l1, j)* ve( indx, l2, &
                  j)
210           continue
              vev( l2, l1) = vev( l1, l2)
220        continue
230     continue
!
        return
end subroutine
!DWGHT
subroutine dwght                                                              &
         ( n, m, wt, ldwt, ld2wt, t, wtt)
!***Begin Prologue  DWGHT
!***Refer to  ODR
!***Routines Called  (NONE)
!***Date Written   860529   (YYMMDD)
!***Revision Date  920304   (YYMMDD)
!***Purpose  Scale matrix T using WT, i.e., compute WTT = WT*T
!***End Prologue  DWGHT
!
!...Used modules
        use odrpack_kinds,only: wp
!
!...Scalar arguments
        integer                                                               &
         ldwt, ld2wt, m, n
!
!...Array arguments
        real(kind = wp)                                                       &
         t(:,:), wt(:,:,:), wtt(:,:)
!
!...Local scalars
        real(kind = wp)                                                       &
         temp, zero
        integer                                                               &
         i, j, k
!
!...Data statements
        data                                                                  &
         zero                                                                 &
         /0.0E0_wp/
!
!...Variable Definitions (alphabetically)
!       I:       An indexing variable.
!       J:       An indexing variable.
!       K:       An indexing variable.
!       LDWT:    The leading dimension of array WT.
!       LD2WT:   The second dimension of array WT.
!       M:       The number of columns of data in T.
!       N:       The number of rows of data in T.
!       T:       The array being scaled by WT.
!       TEMP:    A temporary scalar.
!       WT:      The weights.
!       WTT:     The results of weighting array T by WT.
!       Array WTT can be the same as T only if the arrays in WT
!       are upper triangular with zeros below the diagonal.
!       ZERO:    The value 0.0E0_wp.
!
!
!***First executable statement  DWGHT
!
!
        if ( n .eq. 0 .or. m .eq. 0) return
!
        if ( wt(1,1,1) .ge. zero) then
           if ( ldwt .ge. n) then
              if ( ld2wt .ge. m) then
!  WT is an N-array of M by M matrices
                 do 130 i = 1, n
                    do 120 j = 1, m
                       temp = zero
                       do 110 k = 1, m
                          temp = temp+ wt( i, j, k)* t( i, k)
110                    continue
                       wtt( i, j) = temp
120                 continue
130              continue
              else
!  WT is an N-array of diagonal matrices
                 do 230 i = 1, n
                    do 220 j = 1, m
                       wtt( i, j) = wt( i,1, j)* t( i, j)
220                 continue
230              continue
              endif
           else
              if ( ld2wt .ge. m) then
!  WT is an M by M matrix
                 do 330 i = 1, n
                    do 320 j = 1, m
                       temp = zero
                       do 310 k = 1, m
                          temp = temp+ wt(1, j, k)* t( i, k)
310                    continue
                       wtt( i, j) = temp
320                 continue
330              continue
              else
!  WT is a diagonal matrice
                 do 430 i = 1, n
                    do 420 j = 1, m
                       wtt( i, j) = wt(1,1, j)* t( i, j)
420                 continue
430              continue
              endif
           endif
        else
!  WT is a scalar
           do 520 j = 1, m
              do 510 i = 1, n
                 wtt( i, j) = abs( wt(1,1,1))* t( i, j)
510           continue
520        continue
        endif
!
        return
end subroutine
!DWINF
subroutine dwinf                                                              &
         ( n, m, np, nq, ldwe, ld2we, isodr,                                  &
         deltai, epsi, xplusi, fni, sdi, vcvi,                                &
         rvari, wssi, wssdei, wssepi, rcondi, etai,                           &
         olmavi, taui, alphai, actrsi, pnormi, rnorsi, prersi,                &
         partli, sstoli, taufci, epsmai,                                      &
         beta0i, betaci, betasi, betani, si, ssi, ssfi, qrauxi, ui,           &
         fsi, fjacbi, we1i, diffi,                                            &
         deltsi, deltni, ti, tti, omegai, fjacdi,                             &
         wrk1i, wrk2i, wrk3i, wrk4i, wrk5i, wrk6i, wrk7i,                     &
         loweri, upperi,                                                      &
         lwkmn)
!***Begin Prologue  DWINF
!***Refer to  ODR
!***Routines Called  (NONE)
!***Date Written   860529   (YYMMDD)
!***Revision Date  920619   (YYMMDD)
!***Purpose  Set storage locations within REAL (KIND=wp) work space
!***End Prologue  DWINF
!
!...Scalar arguments
        integer                                                               &
         actrsi, alphai, betaci, betani, betasi, beta0i, deltai, deltni,      &
         deltsi, diffi, epsi, epsmai, etai, fjacbi, fjacdi, fni, fsi, ldwe,   &
         ld2we, loweri, lwkmn, m, n, np, nq, olmavi, omegai, partli, pnormi,  &
         prersi, qrauxi, rcondi, rnorsi, rvari, sdi, si, ssfi, ssi, sstoli,   &
         taufci, taui, ti, tti, ui, upperi, vcvi, we1i, wrk1i, wrk2i, wrk3i,  &
         wrk4i, wrk5i, wrk6i, wrk7i, wssi, wssdei, wssepi, xplusi
        logical                                                               &
         isodr
!
!...Local scalars
        integer                                                               &
         next
!------------^-----------------------------------------------------------------
!!! FPT - 1273 Fortran auxiliary keyword used as identifier name.
!------------------------------------------------------------------------------
!
!...Variable Definitions (alphabetically)
!       ACTRSI:  The location in array WORK of variable ACTRS.
!       ALPHAI:  The location in array WORK of variable ALPHA.
!       BETACI:  The starting location in array WORK of array BETAC.
!       BETANI:  The starting location in array WORK of array BETAN.
!       BETASI:  The starting location in array WORK of array BETAS.
!       BETA0I:  The starting location in array WORK of array BETA0.
!       DELTAI:  The starting location in array WORK of array DELTA.
!       DELTNI:  The starting location in array WORK of array DELTAN.
!       DELTSI:  The starting location in array WORK of array DELTAS.
!       DIFFI:   The starting location in array WORK of array DIFF.
!       EPSI:    The starting location in array WORK of array EPS.
!       EPSMAI:  The location in array WORK of variable EPSMAC.
!       ETAI:    The location in array WORK of variable ETA.
!       FJACBI:  The starting location in array WORK of array FJACB.
!       FJACDI:  The starting location in array WORK of array FJACD.
!       FNI:     The starting location in array WORK of array FN.
!       FSI:     The starting location in array WORK of array FS.
!       ISODR:   The variable designating whether the solution is by ODR
!       (ISODR=TRUE) or by OLS (ISODR=FALSE).
!       LDWE:    The leading dimension of array WE.
!       LD2WE:   The second dimension of array WE.
!       LWKMN:   The minimum acceptable length of vector work.
!       M:       The number of columns of data in the explanatory variable.
!       N:       The number of observations.
!       NEXT:    The next available location with WORK.
!       NP:      The number of function parameters.
!       NQ:      The number of responses per observation.
!       OLMAVI:  The location in array WORK of variable OLMAVG.
!       OMEGAI:  The starting location in array WORK of array OMEGA.
!       PARTLI:  The location in array WORK of variable PARTOL.
!       PNORMI:  The location in array WORK of variable PNORM.
!       PRERSI:  The location in array WORK of variable PRERS.
!       QRAUXI:  The starting location in array WORK of array QRAUX.
!       RCONDI:  The location in array WORK of variable RCONDI.
!       RNORSI:  The location in array WORK of variable RNORMS.
!       RVARI:   The location in array WORK of variable RVAR.
!       SDI:     The starting location in array WORK of array SD.
!       SI:      The starting location in array WORK of array S.
!       SSFI:    The starting location in array WORK of array SSF.
!       SSI:     The starting location in array WORK of array SS.
!       SSTOLI:  The location in array WORK of variable SSTOL.
!       TAUFCI:  The location in array WORK of variable TAUFAC.
!       TAUI:    The location in array WORK of variable TAU.
!       TI:      The starting location in array WORK of array T.
!       TTI:     The starting location in array WORK of array TT.
!       UI:      The starting location in array WORK of array U.
!       VCVI:    The starting location in array WORK of array VCV.
!       WE1I:    The starting location in array WORK of array WE1.
!       WRK1I:   The starting location in array WORK of array WRK1.
!       WRK2I:   The starting location in array WORK of array WRK2.
!       WRK3I:   The starting location in array WORK of array WRK3.
!       WRK4I:   The starting location in array WORK of array WRK4.
!       WRK5I:   The starting location in array WORK of array WRK5.
!       WRK6I:   The starting location in array WORK of array WRK6.
!       WRK7I:   The starting location in array WORK of array WRK7.
!       WSSI:    The location in array WORK of variable WSS.
!       WSSDEI:  The location in array WORK of variable WSSDEL.
!       WSSEPI:  The location in array work of variable WSSEPS.
!       XPLUSI:  The starting location in array WORK of array XPLUSD.
!
!
!***First executable statement  DWINF
!
!
        if ( n .ge. 1 .and. m .ge. 1 .and. np .ge. 1 .and. nq .ge. 1 .and.    &
         ldwe .ge. 1 .and. ld2we .ge. 1) then
!
           deltai = 1
           epsi = deltai+ n* m
           xplusi = epsi+ n* nq
           fni = xplusi+ n* m
           sdi = fni+ n* nq
           vcvi = sdi+ np
           rvari = vcvi+ np* np
!
           wssi = rvari+1
           wssdei = wssi+1
           wssepi = wssdei+1
           rcondi = wssepi+1
           etai = rcondi+1
           olmavi = etai+1
!
           taui = olmavi+1
           alphai = taui+1
           actrsi = alphai+1
           pnormi = actrsi+1
           rnorsi = pnormi+1
           prersi = rnorsi+1
           partli = prersi+1
           sstoli = partli+1
           taufci = sstoli+1
           epsmai = taufci+1
           beta0i = epsmai+1
!
           betaci = beta0i+ np
           betasi = betaci+ np
           betani = betasi+ np
           si = betani+ np
           ssi = si+ np
           ssfi = ssi+ np
           qrauxi = ssfi+ np
           ui = qrauxi+ np
           fsi = ui+ np
!
           fjacbi = fsi+ n* nq
!
           we1i = fjacbi+ n* np* nq
!
           diffi = we1i+ ldwe* ld2we* nq
!
           next = diffi+ nq*( np+ m)
!
           if ( isodr) then
              deltsi = next
              deltni = deltsi+ n* m
              ti = deltni+ n* m
              tti = ti+ n* m
              omegai = tti+ n* m
              fjacdi = omegai+ nq* nq
              wrk1i = fjacdi+ n* m* nq
              next = wrk1i+ n* m* nq
           else
              deltsi = deltai
              deltni = deltai
              ti = deltai
              tti = deltai
              omegai = deltai
              fjacdi = deltai
              wrk1i = deltai
           endif
!
           wrk2i = next
           wrk3i = wrk2i+ n* nq
           wrk4i = wrk3i+ np
           wrk5i = wrk4i+ m* m
           wrk6i = wrk5i+ m
           wrk7i = wrk6i+ n* nq* np
           loweri = wrk7i+5* nq
           upperi = loweri+ np
           next = upperi+ np
!
           lwkmn = next
        else
           deltai = 1
           epsi = 1
           xplusi = 1
           fni = 1
           sdi = 1
           vcvi = 1
           rvari = 1
           wssi = 1
           wssdei = 1
           wssepi = 1
           rcondi = 1
           etai = 1
           olmavi = 1
           taui = 1
           alphai = 1
           actrsi = 1
           pnormi = 1
           rnorsi = 1
           prersi = 1
           partli = 1
           sstoli = 1
           taufci = 1
           epsmai = 1
           beta0i = 1
           betaci = 1
           betasi = 1
           betani = 1
           si = 1
           ssi = 1
           ssfi = 1
           qrauxi = 1
           fsi = 1
           ui = 1
           fjacbi = 1
           we1i = 1
           diffi = 1
           deltsi = 1
           deltni = 1
           ti = 1
           tti = 1
           fjacdi = 1
           omegai = 1
           wrk1i = 1
           wrk2i = 1
           wrk3i = 1
           wrk4i = 1
           wrk5i = 1
           wrk6i = 1
           wrk7i = 1
           loweri = 1
           upperi = 1
           lwkmn = 1
        endif

end subroutine

subroutine dxmy                                                               &
         ( n, m, x, ldx, y, ldy, xmy, ldxmy)
!***Begin Prologue  DXMY
!***Refer to  ODR
!***Routines Called  (NONE)
!***Date Written   860529   (YYMMDD)
!***Revision Date  920304   (YYMMDD)
!***Purpose  Compute XMY = X - Y
!***End Prologue  DXMY
!
!...Used modules
        use odrpack_kinds,only: wp
!
!...Scalar arguments
        integer                                                               &
         ldx, ldxmy, ldy, m, n
!
!...Array arguments
        real(kind = wp)                                                       &
         x( ldx, m), xmy( ldxmy, m), y( ldy, m)
!
!...Local scalars
        integer                                                               &
         i, j
!
!...Variable Definitions (alphabetically)
!       I:       An indexing variable.
!       J:       An indexing variable.
!       LDX:     The leading dimension of array X.
!       LDXMY:   The leading dimension of array XMY.
!       LDY:     The leading dimension of array Y.
!       M:       The number of columns of data in arrays X and Y.
!       N:       The number of rows of data in arrays X and Y.
!       X:       The first of the two arrays.
!       XMY:     The values of X-Y.
!       Y:       The second of the two arrays.
!
!
!***First executable statement  DXMY
!
!
        do 20 j = 1, m
           do 10 i = 1, n
              xmy( i, j) = x( i, j)- y( i, j)
10         continue
20      continue

end subroutine