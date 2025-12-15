/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"
#include <cblas.h>
#include <math.h>

void set_GB_operator_colMajor_poisson1D(double* AB, int *lab, int *la, int *kv){
    int N = *la;
    int KU = *kv;
    int LDAB = *lab;

    for(int j = 0; j < N; j++)
        for(int i = 0; i < LDAB; i++)
            AB[i + j*LDAB] = 0.0;

    for(int j = 0; j < N; j++)
    {
        AB[KU + 1 + j*LDAB] = 2.0;
        if(j < N-1)
            AB[KU+2 + j*LDAB] = -1.0;
        if(j > 0)
            AB[KU + j*LDAB] = -1.0;
    }
}

void set_GB_operator_colMajor_poisson1D_Id(double* AB, int *lab, int *la, int *kv){
    int N = *la;
    int KU = *kv;
    int LDAB = *lab;

    for(int j = 0; j < N; j++)
        for(int i = 0; i < LDAB; i++)
            AB[i + j*LDAB] = 0.0;

    for(int j = 0; j < N; j++)
        AB[KU + 1 + j*LDAB] = 1.0;
}

void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1){
    int N = *la;

    for(int i = 0; i < N; i++)
        RHS[i] = 0.0;

    RHS[0] += *BC0;
    RHS[N - 1] += *BC1;
}  

void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1){
    int N = *la;
  
    for(int i = 0; i < N; i++){
        double x = X[i];
        EX_SOL[i] = -0.5*x*x + (0.5 + *BC1 - *BC0)*x + *BC0;
    }
}  

void set_grid_points_1D(double* x, int* la){
    int N = *la;
    double h = 1.0/(N+1);
  
    for(int i = 0; i < N; i++){
        x[i] = (i+1)*h;
    }
}

double relative_forward_error(double* x, double* y, int* la){
    int N = *la;
    double* diff = (double*)malloc(N * sizeof(double));
  
    cblas_dcopy(N, x, 1, diff, 1);
    cblas_daxpy(N, -1.0, y, 1, diff, 1);
  
    double norm_diff = cblas_dnrm2(N, diff, 1);
    double norm_x = cblas_dnrm2(N, x, 1);
  
    free(diff);
  
    if(norm_x < 1e-15)
        return norm_diff;
  
    return norm_diff / norm_x;
}

int indexABCol(int i, int j, int *lab){
    int KU = (*lab - 1)/2;
    int LDAB = *lab;
  
    if(i < 0 || j < 0) return -1;
  
    int row_in_AB = KU + i - j;
    if(row_in_AB < 0 || row_in_AB >= LDAB) 
        return -1;
  
    return row_in_AB + j*LDAB;
}

int dgbtrftridiag(int *la, int*n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info){
    int N = *la;
    int KU = *ku;
    int LDAB = *lab;
  
    *info = 0;
  
    for(int j = 0; j < N; j++){
        if(j < N-1){
            if(fabs(AB[1 + j*LDAB]) > 1e-15){
                double factor = AB[2 + j*LDAB] / AB[1 + j*LDAB];
                AB[2 + j*LDAB] = factor;
                AB[1 + (j+1)*LDAB] -= factor * AB[0 + j*LDAB];
            }
        }
        ipiv[j] = j+1;
    }
  
    return *info;
}
