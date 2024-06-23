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
 * @param lun Logical unit number.
 */
ODRPACK_EXTERN void close_file(int *lun);

/**
 * @brief User-supplied function for evaluating the model, computing predicted values and their Jacobians.
 *
 * @param n       Number of observations.
 * @param m       Number of columns of data in the independent variable.
 * @param np      Number of function parameters.
 * @param nq      Number of responses per observation.
 * @param ldn     Leading dimension declarator for `n`, `ldn >= n`.
 * @param ldm     Leading dimension declarator for `m`, `ldm >= m`.
 * @param ldnp    Leading dimension declarator for `np`, `ldnp >= np`.
 * @param beta    Input array [np] of current parameter values.
 * @param xplusd  Input array [m][ldn] of current explanatory variable values, i.e., `x + delta`.
 * @param ifixb   Input array [np] of indicators for fixing parameters `beta`.
 * @param ifixx   Input array [m][ldifx] of indicators for fixing explanatory variable `x`.
 * @param ldifx   Leading dimension of array `ifixx`.
 * @param ideval  Indicator for selecting computation to be performed.
 * @param f       Output array [nq][ldn] for predicted function values.
 * @param fjacb   Output array [nq][ldnp][ldn] for Jacobian with respect to `beta`.
 * @param fjacd   Output array [nq][ldm][ldn] for Jacobian with respect to errors `delta`.
 * @param istop   Output integer for stopping condition. Values:
 *                0 - current `beta` and `x+delta` were acceptable and values were computed successfully,
 *                1 - current `beta` and `x+delta` are not acceptable; ODRPACK95 should select values closer to most recently used values if possible,
 *               -1 - current `beta` and `x+delta` are not acceptable; ODRPACK95 should stop.
 */
typedef void (*odrpack_fcn_callback)(
    int *n,
    int *m,
    int *np,
    int *nq,
    int *ldn,
    int *ldm,
    int *ldnp,
    double *beta,   // const
    double *xplusd, // const
    int *ifixb,     // const
    int *ifixx,     // const
    int *ldifx,
    int *ideval,
    double *f,
    double *fjacb,
    double *fjacd,
    int *istop);

/**
 * @brief Wrapper for the ODR routine including only mandatory arguments.
 *
 * @param fcn     User-supplied subroutine for evaluating the model.
 * @param n       Number of observations.
 * @param m       Number of columns of data in the independent variable.
 * @param np      Number of function parameters.
 * @param nq      Number of responses per observation.
 * @param beta    Input/output array [np] of function parameters.
 * @param y       Input array [nq][n] of dependent variable. Unused when the model is implicit.
 * @param x       Input array [m][n] of explanatory variable.
 * @param ifixb   Input array [np] of indicators for fixing parameters `beta`.
 * @param job     Variable controlling problem initialization and computational method.
 * @param ndigit  Number of accurate digits in the function results, as supplied by the user.
 * @param taufac  Factor used to compute the initial trust region diameter.
 * @param sstol   Sum-of-squares convergence stopping tolerance.
 * @param partol  Parameter convergence stopping tolerance.
 * @param maxint  Maximum number of iterations allowed.
 * @param iprint  Print control variable.
 * @param lunerr  Logical unit number for error messages.
 * @param lunrpt  Logical unit number for computation reports.
 * @param stpb    Input array [np] with relative step for computing finite difference derivatives with respect to `beta`.
 * @param sclb    Input array [np] with scaling values for `beta`.
 * @param info    Variable designating why the computations were stopped.
 * @param lower   Input array [np] with lower bound on `beta`.
 * @param upper   Input array [np] with upper bound on `beta`.
 */
ODRPACK_EXTERN void odr_c(
    odrpack_fcn_callback fcn,
    const int *n,
    const int *m,
    const int *np,
    const int *nq,
    double *beta,
    const double *y,
    const double *x,
    const int *ifixb,

    const int *job,
    const int *ndigit,
    const double *taufac,
    const double *sstol,
    const double *partol,
    const int *maxint,
    const int *iprint,
    const int *lunerr,
    const int *lunrpt,

    const double *stpb,
    const double *sclb,

    int *info,
    const double *lower,
    const double *upper);

/**
 * @brief Short wrapper for the ODR routine including only mandatory arguments.
 *
 * @param fcn    User-supplied subroutine for evaluating the model.
 * @param n      Number of observations.
 * @param m      Number of columns of data in the independent variable.
 * @param np     Number of function parameters.
 * @param nq     Number of responses per observation.
 * @param beta   Input/output array [np] of function parameters.
 * @param y      Input array [nq][n] of dependent variable. Unused when the model is implicit.
 * @param x      Input array [m][n] of explanatory variable.
 * @param job    Variable controlling problem initialization and computational method.
 * @param lower  Input array [np] with lower bound on `beta`.
 * @param upper  Input array [np] with upper bound on `beta`.
 */
ODRPACK_EXTERN void odr_short_c(
    odrpack_fcn_callback fcn,
    const int *n,
    const int *m,
    const int *np,
    const int *nq,
    double *beta,
    const double *y,
    const double *x,
    int *job,
    const double *lower,
    const double *upper);

#endif // ODRPACK_H