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
      real(kind=wp), intent(inout), optional :: delta(:, :)
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
            ! elseif (.not. allocated(delta)) then
            !    linfo5 = 7
            !    linfo4 = 1
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
      lwork = zero
      liwork = 0
      ! if (present(delta)) then
      !    if (.not. allocated(delta)) then
      !       allocate (delta(n, m), stat=linfo4)
      !    end if
      ! end if
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
         !if (allocated(delta)) then
            if (any(shape(delta) .lt. (/n, m/))) then
               linfo1 = linfo1 + 8
            end if
            lwork(1:n*m) = reshape(delta(1:n, 1:m), (/n*m/))
         !end if
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
         !if (allocated(delta)) then
            delta(1:n, 1:m) = reshape(lwork(1:n*m), (/n, m/))
         !end if
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
!       DODPER,DPACK,DSETN,DUNPAC,DWGHT,DWINF
!       DERSTEP
!***Date Written   860529   (YYMMDD)
!***Revision Date  920619   (YYMMDD)
!***Purpose  Perform error checking and initialization, and begin
!       procedure for performing orthogonal distance regression
!       (ODR) or ordinary linear or nonlinear least squares (OLS)
!***End Prologue  DODDRV
!
!...Used modules
        use odrpack_kinds, only: wp, zero, one
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

!...Local scalars
        real(kind = wp) :: epsmac, eta, p5, ten
        integer :: actrsi, alphai, betaci, betani, betasi, beta0i, boundi, deltai,  &
                   deltni, deltsi, diffi, epsmai, etai, fi, fjacbi, fjacdi, fni, fsi, i, &
                   idfi, int2i, iprini, iranki, istop, istopi, jobi, jpvti, k, ldtt,  &
                   ldtti, liwkmn, loweri, luneri, lunrpi, lwkmn, lwrk, maxiti, msgb,    &
                   msgd, neta, netai, nfev, nfevi, niteri, njev, njevi, nnzw, nnzwi,    &
                   npp, nppi, nrow, nrowi, ntol, ntoli, olmavi, omegai, partli, pnormi, &
                   prersi, qrauxi, rcondi, rnorsi, rvari, sdi, si, ssfi, ssi, sstoli,   &
                   taufci, taui, ti, tti, ui, upperi, vcvi, we1i, wrk1i, wrk2i, wrk3i,  &
                   wrk4i, wrk5i, wrk6i, wrk7i, wrk, wssi, wssdei, wssepi, xplusi
        logical :: anajac, cdjac, chkjac, dovcv, implct, initd, isodr, redoj, restrt

!...Local arrays
        real(kind = wp) :: betaj( np)
        integer :: interval( np)

!...External functions
        real(kind = wp), external :: ddot, dnrm2, derstep

!...External subroutines
        external :: dcopy, detaf, dfctrw, dflags, diniwk, diwinf, djck, dodchk, &
                    dodmn, dodper, dpack, dsetn, dunpac, dwinf
!
!...Data statements
        data                                                                  &
         p5, ten                                                              &
         /0.5E0_wp,10.0E0_wp/
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
              !call dxmy( n, nq, work( fni), n, y, ldy, work( fi), n)
              work(fi:fi+(n*nq-1)) = work(fni:fni+(n*nq-1)) - reshape(y(1:n,:), shape=[n*nq]) 
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
            beta, work(xplusi),                                               &
            ifixb, ifixx, ldifx,                                              &
            002, work( fni), work( wrk6i), work( wrk1i),                      &
            istop)
           iwork( istopi) = istop
           if ( istop .eq. 0) then
              iwork( nfevi) = iwork( nfevi)+1
              if ( implct) then
                 call dcopy( n* nq, work( fni),1, work( fi),1)
              else
                 !call dxmy( n, nq, work( fni), n, y, ldy, work( fi), n)
                 work(fi:fi+(n*nq-1)) = work(fni:fni+(n*nq-1)) - reshape(y(1:n,:), shape=[n*nq]) 
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
!       DODPCR,DODVCV,DUNPAC,DWGHT
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
        external :: dacces, dcopy, devjac, dflags, dodlm, dodpcr, dodvcv, dunpac

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
              !call dxmy( n, nq, fn, n, y, ldy, wrk, n)
              wrk(1:n*nq) = reshape(fn - y(1:n,:), [n*nq])
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
              !call dxmy( n, nq, fs, n, y, ldy, f, n)
              f = fs - y(1:n,:)
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
           !call dxmy( n, nq, fs, n, y, ldy, f, n)
           f = fs - y(1:n,:)
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

