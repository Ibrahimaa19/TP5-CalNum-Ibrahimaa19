/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"
#include <cblas.h>
#include <math.h>
#include <lapacke.h>


void set_GB_operator_colMajor_poisson1D(double *AB, int *lab, int *la, int *ku)
{
  int i, j;
  int kl = *ku;   /* tridiagonale */

  /* Initialisation à zéro */
  for (j = 0; j < *la; j++)
    for (i = 0; i < *lab; i++)
      AB[i + j*(*lab)] = 0.0;

  for (j = 0; j < *la; j++) {

    /* Diagonale */
    AB[kl + (*ku) + j*(*lab)] = 2.0;

    /* Sur-diagonale */
    if (j < *la - 1)
      AB[kl + (*ku) - 1 + j*(*lab)] = -1.0;

    /* Sous-diagonale */
    if (j > 0)
      AB[kl + (*ku) + 1 + j*(*lab)] = -1.0;
  }
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



/* Exercice 5: Fonctions LAPACK               */

int dgbsv_tridiag(int *n, int *kl, int *ku, double *AB, 
                  int *lab, double *b, double *x, int *info) {
    int nrhs = 1;
    int* ipiv = (int*)malloc(*n * sizeof(int));
    
    // Copie b dans x (sera écrasé par la solution)
    dcopy_(n, b, &(int){1}, x, &(int){1});
    
    // Appel à la fonction Fortran dgbsv
    dgbsv_(n, kl, ku, &nrhs, AB, lab, ipiv, x, n, info);
    
    free(ipiv);
    return *info;
}

double validate_tridiag_poisson1D(double *x, double *b, int *la)
{
    int N = *la;
    double *Ax = (double *) malloc(N * sizeof(double));
    double *r  = (double *) malloc(N * sizeof(double));

    /* Calcul Ax */
    for (int i = 0; i < N; i++) {
        Ax[i] = 2.0 * x[i];
        if (i > 0)     Ax[i] -= x[i-1];
        if (i < N-1)   Ax[i] -= x[i+1];
    }

    /* r = b - Ax */
    for (int i = 0; i < N; i++) {
        r[i] = b[i] - Ax[i];
    }

    /* Normes */
    double norm_r = cblas_dnrm2(N, r, 1);
    double norm_b = cblas_dnrm2(N, b, 1);

    free(Ax);
    free(r);

    if (norm_b < 1e-15) return norm_r;
    return norm_r / norm_b;
}



/* Exercice 6: LU personnalisée  */

void lu_tridiag_simple(double *AB, int *lab, int *la, int *kv,
                       double *L, double *U, int *info)
{
    int N = *la;
    int ku = *kv;
    int kl = *kv;
    int ldab = *lab;
    int diag = kl + ku;

    *info = 0;

    // Init
    for(int j = 0; j < N; j++)
        for(int i = 0; i < ldab; i++)
            L[i + j*ldab] = U[i + j*ldab] = 0.0;

    for(int j = 0; j < N; j++) {
        U[diag + j*ldab] = AB[diag + j*ldab];
        L[diag + j*ldab] = 1.0;

        if(j < N-1)
            U[diag - 1 + (j+1)*ldab] = AB[diag - 1 + (j+1)*ldab];
    }

    for(int i = 0; i < N-1; i++) {
        if (fabs(U[diag + i*ldab]) < 1e-14) {
            *info = i+1;
            return;
        }

        double l = AB[diag + 1 + i*ldab] / U[diag + i*ldab];
        L[diag + 1 + i*ldab] = l;
        U[diag + (i+1)*ldab] -= l * U[diag - 1 + (i+1)*ldab];
    }
}


void solve_lu_tridiag(double *L, double *U, int *lab, int *la, int *kv,
                      double *b, double *x, int *info)
{
    int N = *la;
    int diag = (*kv) * 2;
    int ldab = *lab;

    double *y = malloc(N*sizeof(double));

    y[0] = b[0];
    for(int i = 1; i < N; i++)
        y[i] = b[i] - L[diag + 1 + (i-1)*ldab] * y[i-1];

    x[N-1] = y[N-1] / U[diag + (N-1)*ldab];
    for(int i = N-2; i >= 0; i--)
        x[i] = (y[i] - U[diag - 1 + (i+1)*ldab] * x[i+1])
               / U[diag + i*ldab];

    free(y);
}
*/