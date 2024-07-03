/*
This is a translation of example 3 from the ODRPACK95 documentation.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../../include/odrpack/odrpack.h"

// User-supplied function for evaluating the model
void fcn(int *n, int *m, int *np, int *nq, int *ldn, int *ldm, int *ldnp, double *beta,
         double *xplusd, int *ifixb, int *ifixx, int *ldifx, int *ideval, double *f,
         double *fjacb, double *fjacd, int *istop)
{

    // Local variables
    double freq, omega, ctheta, stheta, theta, phi, r;
    const double pi = 4 * atan(1.0);
    int i;

    // Check for unacceptable values for this problem
    for (i = 0; i < *n; i++)
    {
        if (xplusd[i * *m] < 0.0)
        {
            *istop = 1;
            return;
        }
    }
    *istop = 0;

    theta = pi * beta[3] * 0.5;
    ctheta = cos(theta);
    stheta = sin(theta);

    // Compute predicted values
    if ((*ideval % 10) >= 1)
    {
        for (i = 0; i < *n; i++)
        {
            freq = xplusd[i];
            omega = pow(2.0 * pi * freq * exp(-beta[2]), beta[3]);
            phi = atan2(omega * stheta, 1 + omega * ctheta);
            r = (beta[0] - beta[1]) * pow(sqrt(pow(1 + omega * ctheta, 2) + pow(omega * stheta, 2)), -beta[4]);
            f[i] = beta[1] + r * cos(beta[4] * phi);
            f[i + *n] = r * sin(beta[4] * phi);
        }
    }
}

int main()
{

#define NP 5
#define N 23
#define M 1
#define NQ 2

    // Variable declarations
    int n = N, m = M, np = NP, nq = NQ;
    int i, info, iprint, j, job, lunerr, lunrpt;
    int ifixx[M][N];
    double beta[NP], x[M][N], y[NQ][N], delta[M][N], wd[M][M][N], we[NQ][NQ][N];

    // Read problem data
    FILE *data_file = fopen("./data3.dat", "r");
    for (i = 0; i < NP; i++)
    {
        fscanf(data_file, "%lf", &beta[i]);
    }
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < m; j++)
        {
            fscanf(data_file, "%lf", &x[j][i]);
        }
        for (j = 0; j < nq; j++)
        {
            fscanf(data_file, "%lf", &y[j][i]);
        }
    }

    // Close files
    fclose(data_file);

    // Specify task as explicit orthogonal distance regression with central difference derivatives
    job = 01010;
    iprint = 2002;

    // Initialize delta and specify first decade of frequencies as fixed
    for (i = 0; i < n; i++)
    {
        if (x[0][i] < 100.0)
        {
            delta[0][i] = 0.0;
            ifixx[0][i] = 0;
        }
        else if (x[0][i] <= 150.0)
        {
            delta[0][i] = 0.0;
            ifixx[0][i] = 1;
        }
        else if (x[0][i] <= 1000.0)
        {
            delta[0][i] = 25.0;
            ifixx[0][i] = 1;
        }
        else if (x[0][i] <= 10000.0)
        {
            delta[0][i] = 560.0;
            ifixx[0][i] = 1;
        }
        else if (x[0][i] <= 100000.0)
        {
            delta[0][i] = 9500.0;
            ifixx[0][i] = 1;
        }
        else
        {
            delta[0][i] = 144000.0;
            ifixx[0][i] = 1;
        }
    }

    // Set weights
    for (i = 0; i < n; i++)
    {
        if (x[0][i] == 100.0 || x[0][i] == 150.0)
        {
            we[0][0][i] = 0.0;
            we[1][0][i] = 0.0;
            we[0][1][i] = 0.0;
            we[1][1][i] = 0.0;
        }
        else
        {
            we[0][0][i] = 559.6;
            we[1][0][i] = -1634.0;
            we[0][1][i] = -1634.0;
            we[1][1][i] = 8397.0;
        }
        wd[0][0][i] = 1.0e-4 / pow(x[0][i], 2);
    }

    // Compute solution
    odr_c(fcn, n, m, NP, nq, &beta, (double *)y, (double *)x, delta, we, wd, ifixx, &job, &iprint, &lunerr, &lunrpt, &info);

    return 0;
}
