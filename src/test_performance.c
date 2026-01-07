#include "lib_poisson1D.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

/* Timer précis */
double get_time() {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

int main() {

    int kl = 1, ku = 1;          // Tridiagonale pure
    int NRHS = 1;
    int lab = 2*kl + ku + 1;     
    int tailles[] = {100, 1000, 5000, 10000, 20000, 50000, 100000};
    int n_tailles = sizeof(tailles)/sizeof(int);
    int NREP = 50;

    double T0 = 0.0, T1 = 1.0;

    FILE *f = fopen("perf_TP.txt", "w");
    if (!f) { perror("fopen"); return 1; }

    printf("n  DGBTRF+DGBTRS   LU_tridiag   DGBSV\n");

    for (int i = 0; i < n_tailles; i++) {

        int n = tailles[i];
        double t_lapack = 0.0, t_tridiag = 0.0, t_dgbsv = 0.0;

        double *AB  = malloc(lab * n * sizeof(double));
        double *RHS = malloc(n * sizeof(double));
        int    *ipiv = malloc(n * sizeof(int));

        if (!AB || !RHS || !ipiv) { perror("malloc"); return 1; }

        for (int r = 0; r < NREP; r++) {

            int info;

            /* ----------------------------- */
            /* 1. LAPACK : DGBTRF + DGBTRS */
            set_GB_operator_colMajor_poisson1D(AB, &lab, &n, &kl);
            set_dense_RHS_DBC_1D(RHS, &n, &T0, &T1);

            double t0 = get_time();
            dgbtrf_(&n, &n, &kl, &ku, AB, &lab, ipiv, &info);
            if (info != 0) { fprintf(stderr,"DGBTRF error n=%d info=%d\n", n, info); return 1; }

            dgbtrs_("N", &n, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &n, &info, 1);
            if (info != 0) { fprintf(stderr,"DGBTRS error n=%d info=%d\n", n, info); return 1; }

            double t1 = get_time();
            t_lapack += (t1 - t0);

            /* ----------------------------- */
            /* 2. LU TRIDIAGONALE */
            set_GB_operator_colMajor_poisson1D(AB, &lab, &n, &kl);
            set_dense_RHS_DBC_1D(RHS, &n, &T0, &T1);

            t0 = get_time();
            dgbtrftridiag(&n, &n, &kl, &ku, AB, &lab, ipiv, &info);
            if (info != 0) { fprintf(stderr,"LU tridiag error n=%d info=%d\n", n, info); return 1; }

            dgbtrs_("N", &n, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &n, &info, 1);
            if (info != 0) { fprintf(stderr,"LU tridiag solve error n=%d info=%d\n", n, info); return 1; }

            t1 = get_time();
            t_tridiag += (t1 - t0);

            /* ----------------------------- */
            /* 3. DGBSV */
            set_GB_operator_colMajor_poisson1D(AB, &lab, &n, &kl);
            set_dense_RHS_DBC_1D(RHS, &n, &T0, &T1);

            t0 = get_time();
            dgbsv_(&n, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &n, &info);
            if (info != 0) { fprintf(stderr,"DGBSV error n=%d info=%d\n", n, info); return 1; }

            t1 = get_time();
            t_dgbsv += (t1 - t0);
        }

        t_lapack /= NREP;
        t_tridiag /= NREP;
        t_dgbsv /= NREP;

        printf("%d  %.3e  %.3e  %.3e\n", n, t_lapack, t_tridiag, t_dgbsv);
        fprintf(f, "%d %.6e %.6e %.6e\n", n, t_lapack, t_tridiag, t_dgbsv);

        free(AB); free(RHS); free(ipiv);
    }

    fclose(f);
    return 0;
}