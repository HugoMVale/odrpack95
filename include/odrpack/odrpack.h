#include <stdbool.h>

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
 * @param lun      `<=>` Logical unit number.
 * @param filename `==>` String containing the file name.
 * @param ierr     `<==` Error code.
 */
ODRPACK_EXTERN void open_file(
    int *lun,
    const char *filename,
    int *ierr);

/**
 * @brief Close a file associated with a specified logical unit number.
 *
 * @param lun    `==>` Logical unit number.
 * @param ierr   `==>` Error code.
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
 * @brief 0-based locations within real work array.
 */
typedef struct
{
    int delta; /**< Starting location of array `delta`. */
    int eps;   /**< Starting location of array `eps`. */
    int xplus; /**< Starting location of array `xplusd`. */
    int fn;    /**< Starting location of array `fn`. */
    int sd;    /**< Starting location of array `sd`. */
    int vcv;   /**< Starting location of array `vcv`. */
    int rvar;  /**< Location of variable `rvar`. */
    int wss;   /**< Location of variable `wss`. */
    int wssde; /**< Location of variable `wssdel`. */
    int wssep; /**< Location of variable `wsep`. */
    int rcond; /**< Location of variable `rcond`. */
    int eta;   /**< Location of variable `eta`. */
    int olmav; /**< Location of variable `olmavg`. */
    int tau;   /**< Location of variable `tau`. */
    int alpha; /**< Location of variable `alpha`. */
    int actrs; /**< Location of variable `actrs`. */
    int pnorm; /**< Location of variable `pnorm`. */
    int rnors; /**< Location of variable `rnorms`. */
    int prers; /**< Location of variable `prers`. */
    int partl; /**< Location of variable `partol`. */
    int sstol; /**< Location of variable `sstol`. */
    int taufc; /**< Location of variable `taufac`. */
    int epsma; /**< Location of variable `epsmac`. */
    int beta0; /**< Starting location of array `beta0`. */
    int betac; /**< Starting location of array `betac`. */
    int betas; /**< Starting location of array `betas`. */
    int betan; /**< Starting location of array `betan`. */
    int s;     /**< Starting location of array `s`. */
    int ss;    /**< Starting location of array `ss`. */
    int ssf;   /**< Starting location of array `ssf`. */
    int qraux; /**< Starting location of array `qraux`. */
    int u;     /**< Starting location of array `u`. */
    int fs;    /**< Starting location of array `fs`. */
    int fjacb; /**< Starting location of array `fjacb`. */
    int we1;   /**< Location of variable `we1`. */
    int diff;  /**< Starting location of array `diff`. */
    int delts; /**< Starting location of array `deltas`. */
    int deltn; /**< Starting location of array `deltan`. */
    int t;     /**< Starting location of array `t`. */
    int tt;    /**< Starting location of array `tt`. */
    int omega; /**< Starting location of array `omega`. */
    int fjacd; /**< Starting location of array `fjacd`. */
    int wrk1;  /**< Starting location of array `wrk1`. */
    int wrk2;  /**< Starting location of array `wrk2`. */
    int wrk3;  /**< Starting location of array `wrk3`. */
    int wrk4;  /**< Starting location of array `wrk4`. */
    int wrk5;  /**< Starting location of array `wrk5`. */
    int wrk6;  /**< Starting location of array `wrk6`. */
    int wrk7;  /**< Starting location of array `wrk7`. */
    int lower; /**< Starting location of array `lower`. */
    int upper; /**< Starting location of array `upper`. */
    int lwkmn; /**< Minimum acceptable length of vector `work`. */
} workidx_t;

/**
 * @brief Set storage locations within real work space.
 *
 * @param n       `==>` Number of observations.
 * @param m       `==>` Number of columns of data in the explanatory variable.
 * @param np      `==>` Number of function parameters.
 * @param nq      `==>` Number of responses per observation.
 * @param ldwe    `==>` Leading dimension of array `we`.
 * @param ld2we   `==>` Second dimension of array `we`.
 * @param isodr   `==>` Variable designating whether the solution is by ODR (`isodr=.true.`) or by OLS (`isodr=.false.`).
 * @param workidx `<==` 0-based indexes of real work array.
 */
ODRPACK_EXTERN void dwinf_c(
    const int *n,
    const int *m,
    const int *np,
    const int *nq,
    const int *ldwe,
    const int *ld2we,
    const _Bool *isodr,
    workidx_t *workidx);

#endif // ODRPACK_H
