/*
This is a translation of example 5 from the ODRPACK95 documentation.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../../include/odrpack/odrpack.h"

// User-supplied function for evaluating the model and its partial derivatives
void fcn(int *n, int *m, int *np, int *nq, int *ldn, int *ldm, int *ldnp, double *beta, double *xplusd,
         int *ifixb, int *ifixx, int *ldifx, int *ideval, double *f, double *fjacb, double *fjacd, int *istop)
{
    int i;

    *istop = 0;

    // Calculate model
    if (*ideval % 10 != 0)
    {
        for (i = 0; i < *n; i++)
        {
            f[i] = beta[0] * exp(beta[1] * xplusd[i]);
        }
    }

    // Calculate model partials with respect to `beta`
    if ((*ideval / 10) % 10 != 0)
    {
        for (i = 0; i < *n; i++)
        {
            fjacb[i] = exp(beta[1] * xplusd[i]);
            fjacb[*n + i] = beta[0] * xplusd[i] * exp(beta[1] * xplusd[i]);
        }
    }

    // Calculate model partials with respect to `delta`
    if ((*ideval / 100) % 10 != 0)
    {
        for (i = 0; i < *n; i++)
        {
            fjacd[i] = beta[0] * beta[1] * exp(beta[1] * xplusd[i]);
        }
    }
}

int main()
{
#define np 2
#define n 4
#define m 1
#define nq 1

    double beta[np] = {2.0, 0.5};
    double lower[np] = {0.0, 0.0};
    double upper[np] = {10.0, 0.9};
    double x[n * m] = {0.982, 1.998, 4.978, 6.01};
    double y[n * nq] = {2.7, 7.4, 148.0, 403.0};

    odr_c(fcn, n, m, np, nq, beta, y, x, lower, upper);

    return 0;
}
