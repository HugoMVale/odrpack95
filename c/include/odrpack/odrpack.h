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
 * @param filename `==>` String containing the file name.
 * @param lun      `<=>` Logical unit number.
 * @param ierr     `<==` Error code (compiler dependent).
 */
ODRPACK_EXTERN void open_file(
    const char *filename,
    int *lun,
    int *ierr);

/**
 * @brief Close a file associated with a specified logical unit number.
 *
 * @param lun    `==>` Logical unit number.
 * @param ierr   `<==` Error code (compiler dependent).
 */
ODRPACK_EXTERN void close_file(
    const int *lun,
    int *ierr);

/**
 * @brief User-supplied function for evaluating the model, computing predicted values and their Jacobians.
 *
 * @param n       `==>` Number of observations.
 * @param m       `==>` Number of columns of data in the independent variable.
 * @param q       `==>` Number of responses per observation.
 * @param np      `==>` Number of function parameters.
 * @param ldifx   `==>` Leading dimension of array `ifixx`, `ifixx ∈ {1, n}`
 * @param beta    `==>` Array [np] of current parameter values.
 * @param xplusd  `==>` Array [m][ldn] of current explanatory variable values, i.e., `x + delta`.
 * @param ifixb   `==>` Array [np] of indicators for fixing parameters `beta`.
 * @param ifixx   `==>` Array [m][ldifx] of indicators for fixing explanatory variable `x`.
 * @param ideval  `==>` Indicator for selecting computation to be performed.
 * @param f       `<==` Array [q][n] for predicted function values.
 * @param fjacb   `<==` Array [q][np][n] for Jacobian with respect to `beta`.
 * @param fjacd   `<==` Array [q][m][n] for Jacobian with respect to errors `delta`.
 * @param istop   `<==` Integer for stopping condition. Values:
 *                0 - current `beta` and `x + delta` were acceptable and values were computed successfully,
 *                1 - current `beta` and `x + delta` are not acceptable; ODRPACK95 should select values closer to most recently used values if possible,
 *               -1 - current `beta` and `x + delta` are not acceptable; ODRPACK95 should stop.
 */
typedef void (*odrpack_fcn_t)(
    const int *n,
    const int *m,
    const int *q,
    const int *np,
    const int *ldifx,
    const double beta[],
    const double xplusd[],
    const int ifixb[],
    const int ifixx[],
    const int *ideval,
    double f[],
    double fjacb[],
    double fjacd[],
    int *istop);

/**
 * @brief "Short-call" wrapper for the ODR routine including mandatory arguments and very few
 * optional arguments.
 *
 * @param fcn    `==>` User-supplied subroutine for evaluating the model.
 * @param n      `==>` Number of observations.
 * @param m      `==>` Number of columns of data in the independent variable.
 * @param q      `==>` Number of responses per observation.
 * @param np     `==>` Number of function parameters.
 * @param beta   `<=>` Array [np] of function parameters.
 * @param y      `==>` Array [q][n] of dependent variable. Unused when the model is implicit.
 * @param x      `==>` Array [m][n] of explanatory variable.
 * @param delta  `<=>` Optional array [m][n] with initial error in the `x` data.
 * @param lower  `==>` Optional array [np] with lower bound on `beta`.
 * @param upper  `==>` Optional array [np] with upper bound on `beta`.
 * @param job    `==>` Optional variable controlling initialization and computational method.
 */
ODRPACK_EXTERN void odr_short_c(
    odrpack_fcn_t fcn,
    const int *n,
    const int *m,
    const int *q,
    const int *np,
    double beta[],
    const double y[],
    const double x[],
    double delta[],
    const double lower[],
    const double upper[],
    const int *job);

/**
 * @brief "Medium-call" wrapper for the ODR routine including mandatory arguments and most
 * commonly used optional arguments.
 *
 * @param fcn    `==>` User-supplied subroutine for evaluating the model.
 * @param n      `==>` Number of observations.
 * @param m      `==>` Number of columns of data in the independent variable.
 * @param q     `==>` Number of responses per observation.
 * @param np     `==>` Number of function parameters.
 * @param ldwe   `==>` Leading dimension of array `we`, `ldwe ∈ {1, n}`.
 * @param ld2we  `==>` Second dimension of array `we`, `ld2we ∈ {1, q}`.
 * @param ldwd   `==>` Leading dimension of array `wd`, `ldwd ∈ {1, n}`.
 * @param ld2wd  `==>` Second dimension of array `wd`, `ld2wd ∈ {1, m}`.
 * @param ldifx  `==>` Leading dimension of array `ifixx`, `ldifx ∈ {1, n}`.
 * @param beta   `<=>` Array [np] of function parameters.
 * @param y      `==>` Array [q][n] of dependent variable. Unused when the model is implicit.
 * @param x      `==>` Array [m][n] of explanatory variable.
 * @param we     `==>` Optional array [q][ld2we][ldwe] with `epsilon` weights.
 * @param wd     `==>` Optional array [m][ld2wd][ldwd] with `delta` weights.
 * @param ifixb  `==>` Optional array [np] with values designating whether the elements of `beta` are fixed at their input values or not.
 * @param ifixx  `==>` Optional array [m][ldifx] with values designating whether the elements of `x` are fixed at their input values or not.
 * @param delta  `<=>` Optional array [m][n] with initial error in the `x` data.
 * @param lower  `==>` Optional array [np] with lower bound on `beta`.
 * @param upper  `==>` Optional array [np] with upper bound on `beta`.
 * @param job    `==>` Optional variable controlling initialization and computational method.
 * @param iprint `==>` Optional print control variable.
 * @param lunerr `==>` Optional logical unit number for error messages.
 * @param lunrpt `==>` Optional logical unit number for computation reports.
 * @param info   `<==` Optional variable designating why the computations were stopped.
 */
ODRPACK_EXTERN void odr_medium_c(
    odrpack_fcn_t fcn,
    const int *n,
    const int *m,
    const int *q,
    const int *np,
    const int *ldwe,
    const int *ld2we,
    const int *ldwd,
    const int *ld2wd,
    const int *ldifx,
    double beta[],
    const double y[],
    const double x[],
    const double we[],
    const double wd[],
    const int ifixb[],
    const int ifixx[],
    double delta[],
    const double lower[],
    const double upper[],
    const int *job,
    const int *iprint,
    const int *lunerr,
    const int *lunrpt,
    int *info);

/**
 * @brief "Long-call" wrapper for the ODR routine including mandatory arguments and all
 * optional arguments.
 *
 * @param fcn    `==>` User-supplied subroutine for evaluating the model.
 * @param n      `==>` Number of observations.
 * @param m      `==>` Number of columns of data in the independent variable.
 * @param q     `==>` Number of responses per observation.
 * @param np     `==>` Number of function parameters.
 * @param ldwe   `==>` Leading dimension of array `we`, `ldwe ∈ {1, n}`.
 * @param ld2we  `==>` Second dimension of array `we`, `ld2we ∈ {1, q}`.
 * @param ldwd   `==>` Leading dimension of array `wd`, `ldwd ∈ {1, n}`.
 * @param ld2wd  `==>` Second dimension of array `wd`, `ld2wd ∈ {1, m}`.
 * @param ldifx  `==>` Leading dimension of array `ifixx`, `ldifx ∈ {1, n}`.
 * @param ldstpd `==>` Leading dimension of array `stpd`, `ldstpd ∈ {1, n}`.
 * @param ldscld `==>` Leading dimension of array `scld`, `ldscld ∈ {1, n}`.
 * @param lrwork `==>` Length of array `rwork`.
 * @param liwork `==>` Length of array `iwork`.
 * @param beta   `<=>` Array [np] of function parameters.
 * @param y      `==>` Array [q][n] of dependent variable. Unused when the model is implicit.
 * @param x      `==>` Array [m][n] of explanatory variable.
 * @param we     `==>` Optional array [q][ld2we][ldwe] with `epsilon` weights.
 * @param wd     `==>` Optional array [m][ld2wd][ldwd] with `delta` weights.
 * @param ifixb  `==>` Optional array [np] with values designating whether the elements of `beta` are fixed at their input values or not.
 * @param ifixx  `==>` Optional array [m][ldifx] with values designating whether the elements of `x` are fixed at their input values or not.
 * @param stpb   `==>` Optional array [np] with relative step for computing finite difference derivatives with respect to `beta`.
 * @param stpd   `==>` Optional array [m][ldstpd] with relative step for computing finite difference derivatives with respect to `delta`.
 * @param sclb   `==>` Optional array [np] with scaling values for `beta`.
 * @param scld   `==>` Optional array [m][ldscld] with scaling values for `delta`.
 * @param delta  `<=>` Optional array [m][n] with initial error in the `x` data.
 * @param lower  `==>` Optional array [np] with lower bound on `beta`.
 * @param upper  `==>` Optional array [np] with upper bound on `beta`.
 * @param rwork  `<=>` Optional real work space.
 * @param iwork  `<=>` Optional integer work space.
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
    odrpack_fcn_t fcn,
    const int *n,
    const int *m,
    const int *q,
    const int *np,
    const int *ldwe,
    const int *ld2we,
    const int *ldwd,
    const int *ld2wd,
    const int *ldifx,
    const int *ldstpd,
    const int *ldscld,
    const int *lrwork,
    const int *liwork,
    double beta[],
    const double y[],
    const double x[],
    const double we[],
    const double wd[],
    const int ifixb[],
    const int ifixx[],
    const double stpb[],
    const double stpd[],
    const double sclb[],
    const double scld[],
    double delta[],
    const double lower[],
    const double upper[],
    double rwork[],
    int iwork[],
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
 * @brief 0-based locations within integer work array.
 */
typedef struct
{
    int msgb;   /**< The starting location in array `iwork` of array `msgb`. */
    int msgd;   /**< The starting location in array `iwork` of array `msgd`. */
    int ifix2;  /**< The starting location in array `iwork` of array `ifix2`. */
    int istop;  /**< The location in array `iwork` of variable `istop`. */
    int nnzw;   /**< The location in array `iwork` of variable `nnzw`. */
    int npp;    /**< The location in array `iwork` of variable `npp`. */
    int idf;    /**< The location in array `iwork` of variable `idf`. */
    int job;    /**< The location in array `iwork` of variable `job`. */
    int iprint; /**< The location in array `iwork` of variable `iprint`. */
    int lunerr; /**< The location in array `iwork` of variable `lunerr`. */
    int lunrpt; /**< The location in array `iwork` of variable `lunrpt`. */
    int nrow;   /**< The location in array `iwork` of variable `nrow`. */
    int ntol;   /**< The location in array `iwork` of variable `ntol`. */
    int neta;   /**< The location in array `iwork` of variable `neta`. */
    int maxit;  /**< The location in array `iwork` of variable `maxit`. */
    int niter;  /**< The location in array `iwork` of variable `niter`. */
    int nfev;   /**< The location in array `iwork` of variable `nfev`. */
    int njev;   /**< The location in array `iwork` of variable `njev`. */
    int int2;   /**< The location in array `iwork` of variable `int2`. */
    int irank;  /**< The location in array `iwork` of variable `irank`. */
    int ldtt;   /**< The location in array `iwork` of variable `ldtt`. */
    int bound;  /**< The location in array `iwork` of variable `bound`. */
    int liwkmn; /**< The minimum acceptable length of array `iwork`. */
} iworkidx_t;

/**
 * @brief 0-based locations within real work array.
 */
typedef struct
{
    int delta;   /**< Starting location of array `delta`. */
    int eps;     /**< Starting location of array `eps`. */
    int xplusd;  /**< Starting location of array `xplusd`. */
    int fn;      /**< Starting location of array `fn`. */
    int sd;      /**< Starting location of array `sd`. */
    int vcv;     /**< Starting location of array `vcv`. */
    int rvar;    /**< Location of variable `rvar`. */
    int wss;     /**< Location of variable `wss`. */
    int wssdel;  /**< Location of variable `wssdel`. */
    int wsseps;  /**< Location of variable `wsep`. */
    int rcond;   /**< Location of variable `rcond`. */
    int eta;     /**< Location of variable `eta`. */
    int olmavg;  /**< Location of variable `olmavg`. */
    int tau;     /**< Location of variable `tau`. */
    int alpha;   /**< Location of variable `alpha`. */
    int actrs;   /**< Location of variable `actrs`. */
    int pnorm;   /**< Location of variable `pnorm`. */
    int rnorms;  /**< Location of variable `rnorms`. */
    int prers;   /**< Location of variable `prers`. */
    int partol;  /**< Location of variable `partol`. */
    int sstol;   /**< Location of variable `sstol`. */
    int taufac;  /**< Location of variable `taufac`. */
    int epsmac;  /**< Location of variable `epsmac`. */
    int beta0;   /**< Starting location of array `beta0`. */
    int betac;   /**< Starting location of array `betac`. */
    int betas;   /**< Starting location of array `betas`. */
    int betan;   /**< Starting location of array `betan`. */
    int s;       /**< Starting location of array `s`. */
    int ss;      /**< Starting location of array `ss`. */
    int ssf;     /**< Starting location of array `ssf`. */
    int qraux;   /**< Starting location of array `qraux`. */
    int u;       /**< Starting location of array `u`. */
    int fs;      /**< Starting location of array `fs`. */
    int fjacb;   /**< Starting location of array `fjacb`. */
    int we1;     /**< Location of variable `we1`. */
    int diff;    /**< Starting location of array `diff`. */
    int deltas;  /**< Starting location of array `deltas`. */
    int deltan;  /**< Starting location of array `deltan`. */
    int t;       /**< Starting location of array `t`. */
    int tt;      /**< Starting location of array `tt`. */
    int omega;   /**< Starting location of array `omega`. */
    int fjacd;   /**< Starting location of array `fjacd`. */
    int wrk1;    /**< Starting location of array `wrk1`. */
    int wrk2;    /**< Starting location of array `wrk2`. */
    int wrk3;    /**< Starting location of array `wrk3`. */
    int wrk4;    /**< Starting location of array `wrk4`. */
    int wrk5;    /**< Starting location of array `wrk5`. */
    int wrk6;    /**< Starting location of array `wrk6`. */
    int wrk7;    /**< Starting location of array `wrk7`. */
    int lower;   /**< Starting location of array `lower`. */
    int upper;   /**< Starting location of array `upper`. */
    int lrwkmin; /**< Minimum acceptable length of vector `rwork`. */
} rworkidx_t;

/**
 * @brief Get storage locations within integer work space.
 *
 * @param m   `==>` Number of columns of data in the explanatory variable.
 * @param q   `==>` Number of responses per observation.
 * @param np  `==>` Number of function parameters.
 * @param iwi `<==` 0-based indexes of integer work array.
 */
ODRPACK_EXTERN void loc_iwork_c(
    const int *m,
    const int *q,
    const int *np,
    iworkidx_t *iwi);

/**
 * @brief Get storage locations within real work space.
 *
 * @param n     `==>` Number of observations.
 * @param m     `==>` Number of columns of data in the explanatory variable.
 * @param q     `==>` Number of responses per observation.
 * @param np    `==>` Number of function parameters.
 * @param ldwe  `==>` Leading dimension of array `we`.
 * @param ld2we `==>` Second dimension of array `we`.
 * @param isodr `==>` Variable designating whether the solution is by ODR (`isodr=.true.`) or by OLS (`isodr=.false.`).
 * @param rwi   `<==` 0-based indexes of real work array.
 */
ODRPACK_EXTERN void loc_rwork_c(
    const int *n,
    const int *m,
    const int *q,
    const int *np,
    const int *ldwe,
    const int *ld2we,
    const bool *isodr,
    rworkidx_t *rwi);

/**
 * @brief Calculate the dimensions of the workspace arrays.
 *
 * @param n      `==>` Number of observations.
 * @param m      `==>` Number of columns of data in the explanatory variable.
 * @param q      `==>` Number of responses per observation.
 * @param np     `==>` Number of function parameters.
 * @param isodr  `==>` Variable designating whether the solution is by ODR (`isodr=.true.`) or by OLS (`isodr=.false.`).
 * @param lrwork `<==` Length of real `rwork` array.
 * @param liwork `<==` Length of integer `iwork` array.
 */
ODRPACK_EXTERN void workspace_dimensions_c(
    const int *n,
    const int *m,
    const int *q,
    const int *np,
    const bool *isodr,
    int *lrwork,
    int *liwork);

#endif  // ODRPACK_H
