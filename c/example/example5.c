/*
This is a translation of example 5 from the ODRPACK95 documentation.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
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

    int n = N, m = M, np = NP, nq = NQ;
    int job = 20;
    double beta[NP] = {2.0, 0.5};
    double lower[NP] = {0.0, 0.0};
    double upper[NP] = {10.0, 0.9};
    double x[M][N] = {{0.982, 1.998, 4.978, 6.01}};
    double y[NQ][N] = {{2.7, 7.4, 148.0, 403.0}};

    // odr_basic_c(fcn, &n, &m, &np, &nq, beta, (double *)y, (double *)x, lower, upper, &job);

    // extra arguments for use with odr_short_c and odr_long_c
    int ldwe = 1;
    int ld2we = 1;
    int ldwd = 1;
    int ld2wd = 1;
    double we[NQ][1][1] = {-1.0};
    double wd[M][1][1] = {-1.0};
    int iprint = 2002;
    int info = 0;

    char *fn = "report5.txt";
    int lunrpt = 0;
    int lunerr = 0;
    int ierr = 0;

    open_file(&lunrpt, fn, &ierr);
    // printf("Error code (ierr): %d\n", ierr);
    lunerr = lunrpt;

    // odr_short_c(fcn, &n, &m, &np, &nq, beta, (double *)y, (double *)x,
    //             (double *)we, &ldwe, &ld2we,
    //             (double *)wd, &ldwd, &ld2wd,
    //             lower, upper,
    //             &job, &iprint, &lunerr, &lunrpt,
    //             &info);

    // extra arguments for use with odr_long_c
    double delta[M][N] = {{0.1, 0.1, 0.1, 0.1}};
    int ifixb[NP] = {1, 1};
    int ifixx[M][1] = {-1};
    int ldifx = 1;
    double stpb[NP] = {-1.0, -1.0};
    double stpd[M][1] = {-1.0};
    int ldstpd = 1;
    double sclb[NP] = {-1.0, -1.0};
    double scld[M][1] = {-1.0};
    int ldscld = 1;
    int ndigit = -1;
    double taufac = -1.0;
    double sstol = -1.0;
    double partol = -1.0;
    int maxit = -1;
    job = 1020;

    odr_long_c(fcn, &n, &m, &np, &nq, beta, (double *)y, (double *)x,
               (double *)we, &ldwe, &ld2we,
               (double *)wd, &ldwd, &ld2wd,
               ifixb, (int *)ifixx, &ldifx,
               stpb, (double *)stpd, &ldstpd,
               sclb, (double *)scld, &ldscld,
               lower, upper,
               (double *)delta,
               &job, &ndigit, &taufac, &sstol, &partol,
               &maxit, &iprint, &lunerr, &lunrpt,
               &info);

    close_file(&lunrpt, &ierr);
    // printf("Error code (ierr): %d\n", ierr);

    workidx_t workidx;
    _Bool isodr = true;

    dwinf_c(&n, &m, &np, &nq, &ldwe, &ld2we, &isodr, &workidx);
    printf("length of `work`: %d\n", workidx.lwkmn);

    // for (int i = 0; i < M; i++)
    // {
    //     for (int j = 0; j < N; j++)
    //     {
    //         printf("%lf ", delta[i][j]);
    //     }
    //     printf("\n");
    // }

    return 0;
}
