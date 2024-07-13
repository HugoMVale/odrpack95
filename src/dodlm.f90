subroutine dodlm &
   (n, m, np, nq, npp, &
    f, fjacb, fjacd, &
    wd, ldwd, ld2wd, ss, tt, ldtt, delta, &
    alpha2, tau, epsfcn, isodr, &
    tfjacb, omega, u, qraux, jpvt, &
    s, t, nlms, rcond, irank, &
    wrk1, wrk2, wrk3, wrk4, wrk5, wrk, lwrk, istopc)
!! Compute Levenberg-Marquardt parameter and steps `s` and `t` using analog of the
!! trust-region Levenberg-Marquardt algorithm.
! Routines Called  DDOT, DNRM2, DODSTP, DSCALE, DWGHT
! Date Written   860529   (YYMMDD)
! Revision Date  920619   (YYMMDD)

   use odrpack_kinds, only: wp, zero
   use odrpack, only: tempret

   integer, intent(in) :: n
      !! The number of observations.
   integer, intent(in) :: m
      !! The number of columns of data in the explanatory variable.
   integer, intent(in) :: np
      !! The number of function parameters.
   integer, intent(in) :: nq
      !! The number of responses per observation.
   integer, intent(in) :: npp
      !! The number of function parameters being estimated.
   real(kind=wp), intent(in) :: f(n, nq)
      !! The (weighted) estimated values of `epsilon`.
   real(kind=wp), intent(in) :: fjacb(n, np, nq)
      !! The Jacobian with respect to `beta`.
   real(kind=wp), intent(in) :: fjacd(n, m, nq)
      !! The Jacobian with respect to `delta`.
   real(kind=wp), intent(in) :: wd(ldwd, ld2wd, m)
      !! The `delta` weights.
   integer, intent(in) :: ldwd
      !! The leading dimension of array `wd`.
   integer, intent(in) :: ld2wd
      !! The second dimension of array `wd`.
   real(kind=wp), intent(in) :: ss(np)
      !! The scaling values used for the unfixed `beta`s.
   real(kind=wp), intent(in) :: tt(ldtt, m)
      !! The scale used for the `delta`s.
   integer, intent(in) :: ldtt
      !! The leading dimension of array `tt`.
   real(kind=wp), intent(in) :: delta(n, m)
      !! The estimated errors in the explanatory variables.
   real(kind=wp), intent(inout) :: alpha2
      !! The current Levenberg-Marquardt parameter.
   real(kind=wp), intent(inout) :: tau
      !! The trust region diameter.
   real(kind=wp), intent(in) :: epsfcn
      !! The function's precision.
   logical, intent(in) :: isodr
      !! The variable designating whether the solution is by ODR (`isodr = .true.`)
      !! or by OLS (`isodr = .false.`).
   real(kind=wp), intent(out) :: tfjacb(n, nq, np)
      !! The array `omega*fjacb`.
   real(kind=wp), intent(out) :: omega(nq, nq)
      !! The array `(I-fjacd*inv(p)*trans(fjacd))**(-1/2)`.
   real(kind=wp), intent(out) :: u(np)
      !! The approximate null vector for tfjacb.
   real(kind=wp), intent(out) :: qraux(np)
      !! The array required to recover the orthogonal part of the Q-R decomposition.
   integer, intent(out) :: jpvt(np)
      !! The pivot vector.
   real(kind=wp), intent(out) :: s(np)
      !! The step for `beta`.
   real(kind=wp), intent(out) :: t(n, m)
      !! The step for `delta`.
   integer, intent(out) :: nlms
      !! The number of Levenberg-Marquardt steps taken.
   real(kind=wp), intent(out) :: rcond
      !! The approximate reciprocal condition of `tfjacb`.
   integer, intent(out) :: irank
      !! The rank deficiency of the Jacobian wrt `beta`.
   real(kind=wp), intent(out) :: wrk1(n, nq, m)
      !! A work array of `(n, nq, m)` elements.
   real(kind=wp), intent(out) :: wrk2(n, nq)
      !! A work array of (n, nq) elements.
   real(kind=wp), intent(out) :: wrk3(np)
      !! A work array of `(np)` elements.
   real(kind=wp), intent(out) :: wrk4(m, m)
      !! A work array of `(m, m)` elements.
   real(kind=wp), intent(out) :: wrk5(m)
      !! A work array of `(m)` elements.
   real(kind=wp), intent(out) :: wrk(lwrk)
      !! A work array of `(lwrk)` elements, _equivalenced_ to `wrk1` and `wrk2`.
   integer, intent(in) :: lwrk
      !! The length of vector `wrk`.
   integer, intent(out) :: istopc
      !! The variable designating whether the computations were stopped due to some other
      !! numerical error detected within subroutine `dodstp`.

   ! Local scalars
   real(kind=wp), parameter :: p001 = 0.001_wp, p1 = 0.1_wp
   real(kind=wp) :: alpha1, alphan, bot, phi1, phi2, sa, top
   integer :: i, iwrk, j, k
   logical :: forvcv

   ! External functions
   real(kind=wp), external :: ddot, dnrm2

   ! External subroutines
   external :: dodstp, dscale

   ! Interface blocks
   interface
      subroutine dwght(n, m, wt, ldwt, ld2wt, t, wtt)
         use odrpack_kinds, only: wp
         integer, intent(in) :: n, m, ldwt, ld2wt
         real(kind=wp), intent(in) :: t(:, :), wt(:, :, :)
         real(kind=wp), intent(out) :: wtt(:, :)
      end subroutine
   end interface

   ! Variable Definitions (alphabetically)
   !  ALPHAN:  The new Levenberg-Marquardt parameter.
   !  ALPHA1:  The previous Levenberg-Marquardt parameter.
   !  ALPHA2:  The current Levenberg-Marquardt parameter.
   !  BOT:     The lower limit for setting ALPHA.
   !  DELTA:   The estimated errors in the explanatory variables.
   !  EPSFCN:  The function's precision.
   !  F:       The (weighted) estimated values of EPSILON.
   !  FJACB:   The Jacobian with respect to BETA.
   !  FJACD:   The Jacobian with respect to DELTA.
   !  FORVCV:  The variable designating whether this subroutine was called to set up for the
   !           covariance matrix computations (FORVCV=TRUE) or not (FORVCV=FALSE).
   !  I:       An indexing variable.
   !  IRANK:   The rank deficiency of the Jacobian wrt BETA.
   !  ISODR:   The variable designating whether the solution is by ODR (ISODR=TRUE) or
   !           by OLS (ISODR=FALSE).
   !  ISTOPC:  The variable designating whether the computations were stoped due to some
   !           other numerical error detected within subroutine DODSTP.
   !  IWRK:    An indexing variable.
   !  J:       An indexing variable.
   !  K:       An indexing variable.
   !  L:       An indexing variable.
   !  JPVT:    The pivot vector.
   !  LDTT:    The leading dimension of array TT.
   !  LDWD:    The leading dimension of array WD.
   !  LD2WD:   The second dimension of array WD.
   !  LWRK:    The length of vector WRK.
   !  M:       The number of columns of data in the explanatory variable.
   !  N:       The number of observations.
   !  NLMS:    The number of Levenberg-Marquardt steps taken.
   !  NP:      The number of function parameters.
   !  NPP:     The number of function parameters being estimated.
   !  NQ:      The number of responses per observation.
   !  OMEGA:   The array (I-FJACD*INV(P)*trans(FJACD))**(-1/2)  where
   !            P = trans(FJACD)*FJACD + D**2 + ALPHA*TT**2
   !  P001:    The value 0.001E0_wp
   !  P1:      The value 0.1E0_wp
   !  PHI1:    The previous difference between the norm of the scaled step and the trust
   !           region diameter.
   !  PHI2:    The current difference between the norm of the scaled step and the trust region
   !           diameter.
   !  QRAUX:   The array required to recover the orthogonal part of the Q-R decomposition.
   !  RCOND:   The approximate reciprocal condition of TFJACB.
   !  S:       The step for BETA.
   !  SA:      The scalar PHI2*(ALPHA1-ALPHA2)/(PHI1-PHI2).
   !  SS:      The scaling values used for the unfixed BETAS.
   !  T:       The step for DELTA.
   !  TAU:     The trust region diameter.
   !  TFJACB:  The array OMEGA*FJACB.
   !  TOP:     The upper limit for setting ALPHA.
   !  TT:      The scale used for the DELTA'S.
   !  U:       The approximate null vector for TFJACB.
   !  WD:      The DELTA weights.
   !  WRK:     A work array of (LWRK) elements, equivalenced to WRK1 and WRK2.
   !  WRK1:    A work array of (N by NQ by M) elements.
   !  WRK2:    A work array of (N by NQ) elements.
   !  WRK3:    A work array of (NP) elements.
   !  WRK4:    A work array of (M by M) elements.
   !  WRK5:    A work array of (M) elements.

   forvcv = .false.
   istopc = 0

   ! Compute full Gauss-Newton step (ALPHA=0)
   alpha1 = zero
   call dodstp(n, m, np, nq, npp, &
               f, fjacb, fjacd, &
               wd, ldwd, ld2wd, ss, tt, ldtt, delta, &
               alpha1, epsfcn, isodr, &
               tfjacb, omega, u, qraux, jpvt, &
               s, t, phi1, irank, rcond, forvcv, &
               wrk1, wrk2, wrk3, wrk4, wrk5, wrk, lwrk, istopc)
   if (istopc .ne. 0) then
      return
   end if

   ! Initialize TAU if necessary
   if (tau .lt. zero) then
      tau = abs(tau)*phi1
   end if

   ! Check if full Gauss-Newton step is optimal
   if ((phi1 - tau) .le. p1*tau) then
      nlms = 1
      alpha2 = zero
      return
   end if

   ! Full Gauss-Newton step is outside trust region - find locally constrained optimal step
   phi1 = phi1 - tau

   ! Initialize upper and lower bounds for ALPHA
   bot = zero

   do k = 1, npp
      tfjacb(1:n, 1:nq, k) = fjacb(1:n, k, 1:nq)
      wrk(k) = ddot(n*nq, tfjacb(1, 1, k), 1, f(1, 1), 1)
   end do
   call dscale(npp, 1, ss, npp, wrk, npp, wrk, npp)

   if (isodr) then
      call dwght(n, m, wd, ldwd, ld2wd, delta, tempret(1:n, 1:m))
      wrk(npp + 1:npp + 1 + n*m - 1) = reshape(tempret(1:n, 1:m), (/n*m/))
      iwrk = npp
      do j = 1, m
         do i = 1, n
            iwrk = iwrk + 1
            wrk(iwrk) = wrk(iwrk) + ddot(nq, fjacd(i, j, 1), n*m, f(i, 1), n)
         end do
      end do
      call dscale(n, m, tt, ldtt, wrk(npp + 1), n, wrk(npp + 1), n)
      top = dnrm2(npp + n*m, wrk, 1)/tau
   else
      top = dnrm2(npp, wrk, 1)/tau
   end if

   if (alpha2 .gt. top .or. alpha2 .eq. zero) then
      alpha2 = p001*top
   end if

   ! Main loop

   do i = 1, 10

      ! Compute locally constrained steps S and T and PHI(ALPHA) for current value of ALPHA
      call dodstp(n, m, np, nq, npp, &
                  f, fjacb, fjacd, &
                  wd, ldwd, ld2wd, ss, tt, ldtt, delta, &
                  alpha2, epsfcn, isodr, &
                  tfjacb, omega, u, qraux, jpvt, &
                  s, t, phi2, irank, rcond, forvcv, &
                  wrk1, wrk2, wrk3, wrk4, wrk5, wrk, lwrk, istopc)
      if (istopc .ne. 0) then
         return
      end if
      phi2 = phi2 - tau

      ! Check whether current step is optimal
      if (abs(phi2) .le. p1*tau .or. (alpha2 .eq. bot .and. phi2 .lt. zero)) then
         nlms = i + 1
         return
      end if

      ! Current step is not optimaL

      ! Update bounds for ALPHA and compute new ALPHA
      if (phi1 - phi2 .eq. zero) then
         nlms = 12
         return
      end if
      sa = phi2*(alpha1 - alpha2)/(phi1 - phi2)
      if (phi2 .lt. zero) then
         top = min(top, alpha2)
      else
         bot = max(bot, alpha2)
      end if
      if (phi1*phi2 .gt. zero) then
         bot = max(bot, alpha2 - sa)
      else
         top = min(top, alpha2 - sa)
      end if

      alphan = alpha2 - sa*(phi1 + tau)/tau
      if (alphan .ge. top .or. alphan .le. bot) then
         alphan = max(p001*top, sqrt(top*bot))
      end if

      ! Get ready for next iteration
      alpha1 = alpha2
      alpha2 = alphan
      phi1 = phi2

   end do

   ! Set NLMS to indicate an optimal step could not be found in 10 trys
   nlms = 12

end subroutine dodlm
