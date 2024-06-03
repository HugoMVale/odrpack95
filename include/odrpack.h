#pragma once

#ifdef __cplusplus
#define ODRPACK_EXTERN extern "C"
#else
#define ODRPACK_EXTERN extern
#endif

typedef void (*minpack_func)(
    int /* n */,
    const double * /* x */,
    double * /* fvec */,
    int * /* iflag */,
    void * /* udata */);

/**
 * @brief User-supplied function for evaluating the model, computing predicted values and their Jacobians.
 *
 * @param n       Number of observations.
 * @param m       Number of columns of data in the independent variable.
 * @param np      Number of function parameters.
 * @param nq      Number of responses per observation.
 * @param ldn     Leading dimension declarator for `n`, must be greater than or equal to `n`.
 * @param ldm     Leading dimension declarator for `m`, must be greater than or equal to `m`.
 * @param ldnp    Leading dimension declarator for `np`, must be greater than or equal to `np`.
 * @param beta    Array of current parameter values (size: np).
 * @param xplusd  Array of current explanatory variable values, i.e., `x + delta` (size: ldn by m).
 * @param ifixb   Array of indicators for fixing parameters `beta` (size: np).
 * @param ifixx   Array of indicators for fixing explanatory variable `x` (size: ldifx by m).
 * @param ldifx   Leading dimension of array `ifixx`, must be greater than or equal to `ldifx`.
 * @param ideval  Indicator for selecting computation to be performed.
 * @param f       Output array for predicted function values (size: ldn by nq).
 * @param fjacb   Output array for Jacobian with respect to `beta` (size: ldn by ldnp by nq).
 * @param fjacd   Output array for Jacobian with respect to errors `delta` (size: ldn by ldm by nq).
 * @param istop   Output integer for stopping condition. Values:
 *                0 - current `beta` and `x+delta` were acceptable and values were computed successfully,
 *                1 - current `beta` and `x+delta` are not acceptable; ODRPACK95 should select values closer to most recently used values if possible,
 *               -1 - current `beta` and `x+delta` are not acceptable; ODRPACK95 should stop.
 */
typedef void (*odrpack_fcn)(
    int /* n */,
    int /* m */,
    int /* np */,
    int /* nq */,
    int /* ldn */,
    int /* ldm */,
    int /* ldnp */,
    const double * /* beta */,
    const double * /* xplusd */,
    const int * /* ifixb */,
    const int * /* ifixx */,
    int /* ldifx */,
    int /* ideval */,
    double * /* f */,
    double * /* fjacb */,
    double * /* fjacd */,
    int * /* istop */

);