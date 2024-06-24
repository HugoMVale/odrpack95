/*
This is a translation of example 5 from the ODRPACK95 documentation.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../../include/odrpack/odrpack.h"

// User-supplied function for evaluating the model and its partial derivatives
void fcn(int *n, int *m, int *np, int *nq, int *ldn, int *ldm, int *ldnp, double *beta,
         double *xplusd, int *ifixb, int *ifixx, int *ldifx, int *ideval, double *f,
         double *fjacb, double *fjacd, int *istop)
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
#define NP 2
#define N 4
#define M 1
#define NQ 1

    int n = NP, m = M, np = NP, nq = NQ;
    int job = 20;
    double beta[NP] = {2.0, 0.5};
    double lower[NP] = {0.0, 0.0};
    double upper[NP] = {10.0, 0.9};
    double x[M][N] = {{0.982, 1.998, 4.978, 6.01}};
    double y[NQ][N] = {{2.7, 7.4, 148.0, 403.0}};

    // odr_basic_c(fcn, &n, &m, &np, &nq, beta, (double *)y, (double *)x, lower, upper, &job);

    int ldwe = 1;
    int ld2we = 1;
    int ldwd = 1;
    int ld2wd = 1;
    double we[NQ][1][1] = {-1.0};
    double wd[M][1][1] = {-1.0};
    int iprint = 1001;
    int lunerr = -1;
    int lunrpt = -1;
    int info = 0;

    // int ifixb[NP] = {1, 1};
    // int ndigit = -1;
    // double taufac = -1.0;
    // double sstol = -1.0;
    // double partol = -1.0;
    // int maxint = -1;
    // double *stpb = NULL;
    // double *sclb = NULL;

    odr_short_c(fcn, &n, &m, &np, &nq, beta, (double *)y, (double *)x,
                (double *)we, &ldwe, &ld2we,
                (double *)wd, &ldwd, &ld2wd,
                lower, upper,
                &job, &iprint, &lunerr, &lunrpt,
                &info);

    // odr_basic_c(fcn, &n, &m, &np, &nq, beta, (double *)y, (double *)x,
    //             ifixb,
    //             &job, &ndigit, &taufac, &sstol, &partol, &maxint, &iprint, &lunerr, &lunrpt,
    //             stpb, sclb,
    //             &info,
    //             lower, upper);

    return 0;
}
