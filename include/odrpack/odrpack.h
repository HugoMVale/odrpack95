#ifndef ODRPACK_H
#define ODRPACK_H

#ifdef __cplusplus
#define ODRPACK_EXTERN extern "C"
#else
#define ODRPACK_EXTERN extern
#endif

/**
 * @brief Open a new file associated with a specified logical unit number.
 *
 * @param lun    Logical unit number.
 * @param fn     String containing the file name.
 * @param fnlen  Length of the string containing the file name.
 * @param ierr   Error code.
 */
ODRPACK_EXTERN void open_file(
    int *lun,
    const char *fn,
    const int *fnlen,
    int *ierr);

/**
 * @brief Close a file associated with a specified logical unit number.
 *
 * @param lun    Logical unit number.
 * @param ierr   Error code.
 */
ODRPACK_EXTERN void close_file(
    const int *lun,
    int *ierr);

/**
 * @brief User-supplied function for evaluating the model, computing predicted values and their Jacobians.
 *
 * @param n       `==>` Number of observations.
 * @param m       `==>` Number of columns of data in the independent variable.
 * @param np      `==>` Number of function parameters.
 * @param nq      `==>` Number of responses per observation.
 * @param ldn     `==>` Leading dimension declarator for `n`, `ldn >= n`.
 * @param ldm     `==>` Leading dimension declarator for `m`, `ldm >= m`.
 * @param ldnp    `==>` Leading dimension declarator for `np`, `ldnp >= np`.
 * @param beta    `==>` Array [np] of current parameter values.
 * @param xplusd  `==>` Array [m][ldn] of current explanatory variable values, i.e., `x + delta`.
 * @param ifixb   `==>` Array [np] of indicators for fixing parameters `beta`.
 * @param ifixx   `==>` Array [m][ldifx] of indicators for fixing explanatory variable `x`.
 * @param ldifx   `==>` Leading dimension of array `ifixx`, `ifixx ∈ {1, n}`
 * @param ideval  `==>` Indicator for selecting computation to be performed.
 * @param f       `<==` Array [nq][ldn] for predicted function values.
 * @param fjacb   `<==` Array [nq][ldnp][ldn] for Jacobian with respect to `beta`.
 * @param fjacd   `<==` Array [nq][ldm][ldn] for Jacobian with respect to errors `delta`.
 * @param istop   `<==` Integer for stopping condition. Values:
 *                0 - current `beta` and `x + delta` were acceptable and values were computed successfully,
 *                1 - current `beta` and `x + delta` are not acceptable; ODRPACK95 should select values closer to most recently used values if possible,
 *               -1 - current `beta` and `x + delta` are not acceptable; ODRPACK95 should stop.
 */
typedef void (*odrpack_fcn)(
    int *n,
    int *m,
    int *np,
    int *nq,
    int *ldn,
    int *ldm,
    int *ldnp,
    double *beta,
    double *xplusd,
    int *ifixb,
    int *ifixx,
    int *ldifx,
    int *ideval,
    double *f,
    double *fjacb,
    double *fjacd,
    int *istop);

/**
 * @brief "Basic" wrapper for the ODR routine including mandatory arguments and very few
 * optional arguments.
 *
 * @param fcn    `==>` User-supplied subroutine for evaluating the model.
 * @param n      `==>` Number of observations.
 * @param m      `==>` Number of columns of data in the independent variable.
 * @param np     `==>` Number of function parameters.
 * @param nq     `==>` Number of responses per observation.
 * @param beta   `<=>` Array [np] of function parameters.
 * @param y      `==>` Array [nq][n] of dependent variable. Unused when the model is implicit.
 * @param x      `==>` Array [m][n] of explanatory variable.
 * @param job    `==>` Optional variable controlling initialization and computational method.
 * @param lower  `==>` Optional array [np] with lower bound on `beta`.
 * @param upper  `==>` Optional array [np] with upper bound on `beta`.
 */
ODRPACK_EXTERN void odr_basic_c(
    odrpack_fcn fcn,
    const int *n,
    const int *m,
    const int *np,
    const int *nq,
    double *beta,
    const double *y,
    const double *x,
    const double *lower,
    const double *upper,
    const int *job);

/**
 * @brief "Short" wrapper for the ODR routine including mandatory arguments and most commonly
 * used optional arguments. Similar to the short-call statement of the original ODRPACK `DODR`.
 *
 * @param fcn    `==>` User-supplied subroutine for evaluating the model.
 * @param n      `==>` Number of observations.
 * @param m      `==>` Number of columns of data in the independent variable.
 * @param np     `==>` Number of function parameters.
 * @param nq     `==>` Number of responses per observation.
 * @param beta   `<=>` Array [np] of function parameters.
 * @param y      `==>` Array [nq][n] of dependent variable. Unused when the model is implicit.
 * @param x      `==>` Array [m][n] of explanatory variable.
 * @param we     `==>` Array [nq][ld2we][ldwe] with `epsilon` weights.
 * @param ldwe   `==>` Leading dimension of array `we`, `ldwe ∈ {1, n}`.
 * @param ld2we  `==>` Second dimension of array `we`, `ld2we ∈ {1, n}`.
 * @param wd     `==>` Array [m][ld2wd][ldwd] with `delta` weights.
 * @param ldwd   `==>` Leading dimension of array `wd`, `ldwd ∈ {1, n}`.
 * @param ld2wd  `==>` Second dimension of array `wd`, `ld2wd ∈ {1, m}`.
 * @param lower  `==>` Optional array [np] with lower bound on `beta`.
 * @param upper  `==>` Optional array [np] with upper bound on `beta`.
 * @param job    `==>` Optional variable controlling initialization and computational method.
 * @param iprint `==>` Optional print control variable.
 * @param lunerr `==>` Optional logical unit number for error messages.
 * @param lunrpt `==>` Optional logical unit number for computation reports.
 * @param info   `<==` Optional variable designating why the computations were stopped.
 */
ODRPACK_EXTERN void odr_short_c(
    odrpack_fcn fcn,
    const int *n,
    const int *m,
    const int *np,
    const int *nq,
    double *beta,
    const double *y,
    const double *x,
    const double *we,
    const int *ldwe,
    const int *ld2we,
    const double *wd,
    const int *ldwd,
    const int *ld2wd,
    const double *lower,
    const double *upper,
    const int *job,
    const int *iprint,
    const int *lunerr,
    const int *lunrpt,
    int *info);

/**
 * @brief "Long" wrapper for the ODR routine including mandatory arguments and most commonly
 * used optional arguments. Similar to the long-call statement of the original ODRPACK `DODC`.
 *
 * @param fcn    `==>` User-supplied subroutine for evaluating the model.
 * @param n      `==>` Number of observations.
 * @param m      `==>` Number of columns of data in the independent variable.
 * @param np     `==>` Number of function parameters.
 * @param nq     `==>` Number of responses per observation.
 * @param beta   `<=>` Array [np] of function parameters.
 * @param y      `==>` Array [nq][n] of dependent variable. Unused when the model is implicit.
 * @param x      `==>` Array [m][n] of explanatory variable.
 * @param we     `==>` Array [nq][ld2we][ldwe] with `epsilon` weights.
 * @param ldwe   `==>` Leading dimension of array `we`, `ldwe ∈ {1, n}`.
 * @param ld2we  `==>` Second dimension of array `we`, `ld2we ∈ {1, n}`.
 * @param wd     `==>` Array [m][ld2wd][ldwd] with `delta` weights.
 * @param ldwd   `==>` Leading dimension of array `wd`, `ldwd ∈ {1, n}`.
 * @param ld2wd  `==>` Second dimension of array `wd`, `ld2wd ∈ {1, m}`.
 * @param ifixb  `==>` Array [np] with values designating whether the elements of `beta` are fixed at their input values or not.
 * @param ifixx  `==>` Array [m][ldifx] with values designating whether the elements of `x` are fixed at their input values or not.
 * @param ldifx  `==>` Leading dimension of array `ifixx`, `ldifx ∈ {1, n}`.
 * @param stpb   `==>` Input array [np] with relative step for computing finite difference derivatives with respect to `beta`.
 * @param stpd   `==>` Input array [m][ldstpd] with relative step for computing finite difference derivatives with respect to `delta`.
 * @param ldstpd `==>` Leading dimension of array `stpd`, `ldstpd ∈ {1, n}`.
 * @param sclb   `==>` Input array [np] with scaling values for `beta`.
 * @param scld   `==>` Input array [m][ldscld] with scaling values for `delta`.
 * @param ldscld `==>` Leading dimension of array `scld`, `ldscld ∈ {1, n}`.
 * @param lower  `==>` Optional array [np] with lower bound on `beta`.
 * @param upper  `==>` Optional array [np] with upper bound on `beta`.
 * @param delta  `<=>` Optional array [m][n] with initial error in the `x` data.
 * @param job    `==>` Optional variable controlling initialization and computational method.
 * @param ndigit `==>` Optional number of accurate digits in the function results, as supplied by the user.
 * @param taufac `==>` Optional factor used to compute the initial trust region diameter.
 * @param sstol  `==>` Optional sum-of-squares convergence stopping tolerance.
 * @param partol `==>` Optional parameter convergence stopping tolerance.
 * @param maxit  `==>` Optional maximum number of iterations allowed.
 * @param iprint `==>` Optional print control variable.
 * @param lunerr `==>` Optional logical unit number for error messages.
 * @param lunrpt `==>` Optional logical unit number for computation reports.
 * @param info   `<==` Optional variable designating why the computations were stopped.
 */
ODRPACK_EXTERN void odr_long_c(
    odrpack_fcn fcn,
    const int *n,
    const int *m,
    const int *np,
    const int *nq,
    double *beta,
    const double *y,
    const double *x,
    const double *we,
    const int *ldwe,
    const int *ld2we,
    const double *wd,
    const int *ldwd,
    const int *ld2wd,
    const int *ifixb,
    const int *ifixx,
    const int *ldifx,
    const double *stpb,
    const double *scpd,
    const int *ldstpd,
    const double *sclb,
    const double *scld,
    const int *ldscld,
    const double *lower,
    const double *upper,
    const double *delta,
    const int *job,
    const int *ndigit,
    const double *taufac,
    const double *sstol,
    const double *partol,
    const int *maxit,
    const int *iprint,
    const int *lunerr,
    const int *lunrpt,
    int *info);

/**
 * @brief Set storage locations within real work space.
 *
 * @param n      `==>` Number of observations.
 * @param m      `==>` Number of columns of data in the explanatory variable.
 * @param np     `==>` Number of function parameters.
 * @param nq     `==>` Number of responses per observation.
 * @param ldwe   `==>` Leading dimension of array `we`.
 * @param ld2we  `==>` Second dimension of array `we`.
 * @param isodr  `==>` Variable designating whether the solution is by ODR (`isodr=.true.`) or by OLS (`isodr=.false.`).
 * @param deltai `<==` Starting location in array `work` of array `delta`.
 * @param epsi   `<==` Starting location in array `work` of array `eps`.
 * @param xplusi `<==` Starting location in array `work` of array `xplusd`.
 * @param fni    `<==` Starting location in array `work` of array `fn`.
 * @param sdi    `<==` Starting location in array `work` of array `sd`.
 * @param vcvi   `<==` Starting location in array `work` of array `vcv`.
 * @param rvari  `<==` Location in array `work` of variable `rvar`.
 * @param wssi   `<==` Location in array `work` of variable `wss`.
 * @param wssdei `<==` Location in array `work` of variable `wssdel`.
 * @param wssepi `<==` Location in array `work` of variable `wsep`.
 * @param rcondi `<==` Location in array `work` of variable `rcond`.
 * @param etai   `<==` Location in array `work` of variable `eta`.
 * @param olmavi `<==` Location in array `work` of variable `olmavg`.
 * @param taui   `<==` Location in array `work` of variable `tau`.
 * @param alphai `<==` Location in array `work` of variable `alpha`.
 * @param actrsi `<==` Location in array `work` of variable `actrs`.
 * @param pnormi `<==` Location in array `work` of variable `pnorm`.
 * @param rnorsi `<==` Location in array `work` of variable `rnorms`.
 * @param prersi `<==` Location in array `work` of variable `prers`.
 * @param partli `<==` Location in array `work` of variable `partol`.
 * @param sstoli `<==` Location in array `work` of variable `sstol`.
 * @param taufci `<==` Location in array `work` of variable `taufac`.
 * @param epsmai `<==` Location in array `work` of variable `epsmac`.
 * @param beta0i `<==` Starting location in array `work` of array `beta0`.
 * @param betaci `<==` Starting location in array `work` of array `betac`.
 * @param betasi `<==` Starting location in array `work` of array `betas`.
 * @param betani `<==` Starting location in array `work` of array `betan`.
 * @param si     `<==` Starting location in array `work` of array `s`.
 * @param ssi    `<==` Starting location in array `work` of array `ss`.
 * @param ssfi   `<==` Starting location in array `work` of array `ssf`.
 * @param qrauxi `<==` Starting location in array `work` of array `qraux`.
 * @param ui     `<==` Starting location in array `work` of array `u`.
 * @param fsi    `<==` Starting location in array `work` of array `fs`.
 * @param fjacbi `<==` Starting location in array `work` of array `fjacb`.
 * @param we1i   `<==` Location in array `work` of variable `we1`.
 * @param diffi  `<==` Starting location in array `work` of array `diff`.
 * @param deltsi `<==` Starting location in array `work` of array `deltas`.
 * @param deltni `<==` Starting location in array `work` of array `deltan`.
 * @param ti     `<==` Starting location in array `work` of array `t`.
 * @param tti    `<==` Starting location in array `work` of array `tt`.
 * @param omegai `<==` Starting location in array `work` of array `omega`.
 * @param fjacdi `<==` Starting location in array `work` of array `fjacd`.
 * @param wrk1i  `<==` Starting location in array `work` of array `wrk1`.
 * @param wrk2i  `<==` Starting location in array `work` of array `wrk2`.
 * @param wrk3i  `<==` Starting location in array `work` of array `wrk3`.
 * @param wrk4i  `<==` Starting location in array `work` of array `wrk4`.
 * @param wrk5i  `<==` Starting location in array `work` of array `wrk5`.
 * @param wrk6i  `<==` Starting location in array `work` of array `wrk6`.
 * @param wrk7i  `<==` Starting location in array `work` of array `wrk7`.
 * @param loweri `<==` Starting location in array `work` of array `lower`.
 * @param upperi `<==` Starting location in array `work` of array `upper`.
 * @param lwkmn  `<==` Minimum acceptable length of vector `work`.
 */
ODRPACK_EXTERN void dwinf_c(
    const int *n,
    const int *m,
    const int *np,
    const int *nq,
    const int *ldwe,
    const int *ld2we,
    const int *isodr,
    int *deltai,
    int *epsi,
    int *xplusi,
    int *fni,
    int *sdi,
    int *vcvi,
    int *rvari,
    int *wssi,
    int *wssdei,
    int *wssepi,
    int *rcondi,
    int *etai,
    int *olmavi,
    int *taui,
    int *alphai,
    int *actrsi,
    int *pnormi,
    int *rnorsi,
    int *prersi,
    int *partli,
    int *sstoli,
    int *taufci,
    int *epsmai,
    int *beta0i,
    int *betaci,
    int *betasi,
    int *betani,
    int *si,
    int *ssi,
    int *ssfi,
    int *qrauxi,
    int *ui,
    int *fsi,
    int *fjacbi,
    int *we1i,
    int *diffi,
    int *deltsi,
    int *deltni,
    int *ti,
    int *tti,
    int *omegai,
    int *fjacdi,
    int *wrk1i,
    int *wrk2i,
    int *wrk3i,
    int *wrk4i,
    int *wrk5i,
    int *wrk6i,
    int *wrk7i,
    int *loweri,
    int *upperi,
    int *lwkmn);

#endif // ODRPACK_H
