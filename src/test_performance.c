#include "lib_poisson1D.h"
#include <stdio.h>
#include <time.h>

int main() {
    int kl = 1, ku = 1, kv = 1;
    int lab = kv + kl + ku + 1;
    int NRHS = 1;
    
    /* Tailles à tester */
    int tailles[] = {100, 500, 1000, 2000, 5000};
    int n_tailles = 5;
    
    /* Tableau pour stocker les résultats */
    double temps_lapack[5], temps_manuel[5], temps_dgbsv[5];
        
    printf("Démarrage des tests...\n");
    
    for (int i = 0; i < n_tailles; i++) {
        int la = tailles[i];
        
        double *AB = malloc(lab * la * sizeof(double));
        double *RHS = malloc(la * sizeof(double));
        int *ipiv = malloc(la * sizeof(int));
        
        set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
        double T0 = 0.0, T1 = 1.0;
        set_dense_RHS_DBC_1D(RHS, &la, &T0, &T1);
        
        /* --- Méthode 1: LAPACK --- */
        clock_t t1 = clock();
        int info;
        dgbtrf_(&la, &la, &kl, &ku, AB, &lab, ipiv, &info);
        if (info == 0) {
            int trans_len = 1;
            dgbtrs_("N", &la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &la, &info, trans_len);
        }
        clock_t t2 = clock();
        temps_lapack[i] = (double)(t2 - t1) / CLOCKS_PER_SEC;
        
        /* --- Méthode 2: Mon LU --- */
        set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
        for (int j = 0; j < la; j++) RHS[j] = (j == 0) ? T0 : ((j == la-1) ? T1 : 0.0);
        
        t1 = clock();
        dgbtrftridiag(&la, &la, &kl, &ku, AB, &lab, ipiv, &info);
        if (info == 0) {
            int trans_len = 1;
            dgbtrs_("N", &la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &la, &info, trans_len);
        }
        t2 = clock();
        temps_manuel[i] = (double)(t2 - t1) / CLOCKS_PER_SEC;
        
        /* --- Méthode 3: dgbsv --- */
        set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
        for (int j = 0; j < la; j++) RHS[j] = (j == 0) ? T0 : ((j == la-1) ? T1 : 0.0);
        
        t1 = clock();
        dgbsv_(&la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &la, &info);
        t2 = clock();
        temps_dgbsv[i] = (double)(t2 - t1) / CLOCKS_PER_SEC;
        
        free(AB);
        free(RHS);
        free(ipiv);
    }
        
    /* SAUVEGARDE DANS FICHIERS */
    FILE *f1 = fopen("mesures_lapack.txt", "w");
    FILE *f2 = fopen("mesures_manuel.txt", "w");
    FILE *f3 = fopen("mesures_dgbsv.txt", "w");
    
    if (f1 && f2 && f3) {
        for (int i = 0; i < n_tailles; i++) {
            fprintf(f1, "%d %.6f\n", tailles[i], temps_lapack[i]);
            fprintf(f2, "%d %.6f\n", tailles[i], temps_manuel[i]);
            fprintf(f3, "%d %.6f\n", tailles[i], temps_dgbsv[i]);
        }
        fclose(f1);
        fclose(f2);
        fclose(f3);
    }
    
    return 0;
}