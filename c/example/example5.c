/*
This is an adaptation of example 5 from the ODRPACK95 documentation.
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../include/odrpack/odrpack.h"

// User-supplied function for evaluating the model and its partial derivatives
void fcn(const int *n, const int *m, const int *q, const int *np, const int *ldifx,
         const double beta[], const double xplusd[],
         const int ifixb[], const int ifixx[], const int *ideval,
         double f[], double fjacb[], double fjacd[], int *istop, void *rdata) {
    *istop = 0;

    // Model function
    if (*ideval % 10 > 0) {
        for (int i = 0; i < *n; i++) {
            f[i] = beta[0] * exp(beta[1] * xplusd[i]);
        }
    }

    // Model partial derivatives wrt `beta`
    if ((*ideval / 10) % 10 > 0) {
        for (int i = 0; i < *n; i++) {
            fjacb[i] = exp(beta[1] * xplusd[i]);
            fjacb[*n + i] = beta[0] * xplusd[i] * exp(beta[1] * xplusd[i]);
        }
    }

    // Model partial derivatives wrt `delta`
    if ((*ideval / 100) % 10 > 0) {
        for (int i = 0; i < *n; i++) {
            fjacd[i] = beta[0] * beta[1] * exp(beta[1] * xplusd[i]);
        }
    }
}

int main() {
#define NP 2
#define N 4
#define M 1
#define Q 1
#define LDWE 1
#define LD2WE 1
#define LDWD 1
#define LD2WD 1
#define LDIFX 1
#define LDSTPD 1
#define LDSCLD 1

    int n = N, m = M, np = NP, q = Q;
    double beta[NP] = {2.0, 0.5};
    double lower[NP] = {0.0, 0.0};
    double upper[NP] = {10.0, 0.9};
    double x[M][N] = {{0.982, 1.998, 4.978, 6.01}};
    double y[Q][N] = {{2.7, 7.4, 148.0, 403.0}};
    double delta[M][N] = {{0.0, 0.0, 0.0, 0.0}};

    int ldwe = LDWE, ld2we = LD2WE, ldwd = LDWD, ld2wd = LD2WD;
    double we[Q][LD2WE][LDWE] = {{{-1.0}}};
    double wd[M][LD2WE][LDWD] = {{{-1.0}}};

    int ldifx = LDIFX;
    int ifixb[NP] = {1, 1};
    int ifixx[M][LDIFX] = {{-1}};

    int ldstpd = LDSTPD;
    double stpb[NP] = {-1.0, -1.0};
    double stpd[M][LDSTPD] = {{-1.0}};
    int ldscld = LDSCLD;
    double sclb[NP] = {-1.0, -1.0};
    double scld[M][LDSCLD] = {{-1.0}};

    int ndigit = -1;
    double taufac = -1.0;
    double sstol = -1.0;
    double partol = -1.0;
    int maxit = 45;
    _Bool isodr = true;

    int job = 1020;
    int iprint = 2001;
    int info = 0;

    // Determine workspace requirements
    int lrwork, liwork;
    workspace_dimensions_c(&n, &m, &q, &np, &isodr, &lrwork, &liwork);
    // printf("lrwork: %d\n", lrwork);
    // printf("liwork: %d\n", liwork);

    // Allocate memory for array `rwork`
    double *rwork = (double *)malloc(lrwork * sizeof(double));
    if (rwork == NULL) {
        fprintf(stderr, "Failed to allocate memory for 'rwork' array\n");
    }

    // Allocate memory for array `iwork`
    int *iwork = (int *)malloc(liwork * sizeof(int));
    if (iwork == NULL) {
        fprintf(stderr, "Failed to allocate memory for 'iwork' array\n");
    }

    // Open report file
    char *fn = "report5.txt";
    int lunrpt = 0;
    int lunerr = 0;
    int ierr = 0;
    open_file(fn, &lunrpt, &ierr);
    // printf("Error code (ierr): %d\n", ierr);
    lunerr = lunrpt;

    // Call odr
    odr_long_c(fcn, NULL,
               &n, &m, &q, &np, &ldwe, &ld2we, &ldwd, &ld2wd, &ldifx, &ldstpd, &ldscld,
               &lrwork, &liwork,
               beta,
               (double *)y, (double *)x,
               (double *)we, (double *)wd,
               ifixb, (int *)ifixx,
               stpb, (double *)stpd,
               sclb, (double *)scld,
               (double *)delta,
               lower, upper,
               rwork, iwork,
               &job, &ndigit, &taufac, &sstol, &partol, &maxit,
               &iprint, &lunerr, &lunrpt, &info);

    close_file(&lunrpt, &ierr);
    // printf("Error code (ierr): %d\n", ierr);

    // Get the variable locations within the integer and real work space
    iworkidx_t iwi;
    rworkidx_t rwi;
    loc_iwork_c(&m, &q, &np, &iwi);
    loc_rwork_c(&n, &m, &q, &np, &ldwe, &ld2we, &isodr, &rwi);

    // Print some outputs

    char message[256] = {};
    stop_message_c(info, (char *)message, sizeof(message));
    printf("Stop reason (info = %d): %s\n", info, message);

    for (int i = 0; i < NP; i++) {
        printf("beta[%d] = %f\n", i, beta[i]);
    }

    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            printf("delta[%d][%d] = %f\n", i, j, delta[i][j]);
        }
    }

    printf("nfev: %d\n", iwork[iwi.nfev]);
    printf("rcond : %f\n", rwork[rwi.rcond]);

    printf("\n");
    for (int i = 0; i < liwork; i++) {
        printf("iwork[%d] = %d\n", i, iwork[i]);
    }
    printf("\n");
    for (int i = 0; i < lrwork; i++) {
        printf("rwork[%d] = %f\n", i, rwork[i]);
    }

    free(rwork);
    free(iwork);

    return 0;
}
