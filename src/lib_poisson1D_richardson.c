/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"
#include <math.h>

void eig_poisson1D(double* eigval, int *la) {
    int N = *la;
    double h = 1.0 / (N + 1.0);
    double factor = 2.0 / (h * h);
    
    for (int j = 1; j <= N; j++) {
        // Valeurs propres: lambda_j = 2/h^2 * (1 - cos(j*π/(N+1)))
        eigval[j-1] = factor * (1.0 - cos(j * M_PI / (N + 1.0)));
    }
}

double eigmax_poisson1D(int *la) {
    int N = *la;
    double h = 1.0 / (N + 1.0);
    // Plus grande valeur propre: j = N
    return (2.0 / (h * h)) * (1.0 - cos(N * M_PI / (N + 1.0)));
}

double eigmin_poisson1D(int *la) {
    int N = *la;
    double h = 1.0 / (N + 1.0);
    // Plus petite valeur propre: j = 1
    return (2.0 / (h * h)) * (1.0 - cos(M_PI / (N + 1.0)));
}

double richardson_alpha_opt(int *la) {
    double lmax = eigmax_poisson1D(la);
    double lmin = eigmin_poisson1D(la);
    // alpha optimal = 2 / (lambda_max + lambda_min)
    return 2.0 / (lmax + lmin);
}

/**
 * Solve linear system Ax=b using Richardson iteration with fixed relaxation parameter alpha.
 * The iteration is: x^(k+1) = x^(k) + alpha*(b - A*x^(k))
 * Stops when ||b - A*x^(k)||_2  / ||b||_2 < tol or when reaching maxit iterations.
 */
void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
  // TODO: Implement Richardson iteration
  // 1. Compute residual r = b - A*x (use dgbmv for matrix-vector product)
  // 2. Update x = x + alpha*r (use daxpy)
  // 3. Check convergence: ||r||_2 < tol (use dnrm2)
  // 4. Store residual norm in resvec and repeat
  int i, iter;
  double *r = (double *)malloc((*la) * sizeof(double));  // résidu
  double *Ax = (double *)malloc((*la) * sizeof(double)); // A*x
  double norm_r, norm_b;
  char trans = 'N';
  int nrhs = 1;
  
  for (i = 0; i < *la; i++) {
      X[i] = 0.0;
  }
  
  // Calcul de la norme de b pour le critère relatif
  norm_b = cblas_dnrm2(*la, RHS, 1);
  if (norm_b < 1e-15) norm_b = 1.0;
  
  for (iter = 0; iter < *maxit; iter++) {
      // Ax = A*X
      cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku,
                  1.0, AB, *lab, X, 1, 0.0, Ax, 1);
      
      // r = RHS - Ax
      cblas_dcopy(*la, RHS, 1, r, 1);
      cblas_daxpy(*la, -1.0, Ax, 1, r, 1);
      
      norm_r = cblas_dnrm2(*la, r, 1);
      
      if (resvec != NULL) {
          resvec[iter] = norm_r;
      }
      
      if (norm_r / norm_b < *tol) {
          break;
      }
      
      cblas_daxpy(*la, *alpha_rich, r, 1, X, 1);
  }
  
  *nbite = iter + 1;
  
  free(r);
  free(Ax);
}

/**
 * Extract MB for Jacobi method from tridiagonal matrix.
 * Such as the Jacobi iterative process is: x^(k+1) = x^(k) + D^(-1)*(b - A*x^(k))
 */
void extract_MB_jacobi_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
  // TODO: Extract diagonal elements from AB and store in MB
  // MB should contain only the diagonal of A
  int i;
    
  for (i = 0; i < (*lab) * (*la); i++) {
      MB[i] = 0.0;
  }
  
  for (i = 0; i < *la; i++) {
      MB[(*kv) + i * (*lab)] = AB[(*kv) + i * (*lab)];
  }
}

/**
 * Extract MB for Gauss-Seidel method from tridiagonal matrix.
 * Such as the Gauss-Seidel iterative process is: x^(k+1) = x^(k) + (D-E)^(-1)*(b - A*x^(k))
 */
void extract_MB_gauss_seidel_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
  // TODO: Extract diagonal and lower diagonal from AB
  // MB should contain the lower triangular part (including diagonal) of A
  int i;  
  
  for (i = 0; i < (*lab) * (*la); i++) {
      MB[i] = 0.0;
  }
  
  
  for (i = 0; i < *la; i++) {
      MB[(*kv) + i * (*lab)] = AB[(*kv) + i * (*lab)];
  }
  
  for (i = 0; i < *la - 1; i++) {
      MB[(*kv + 1) + i * (*lab)] = AB[(*kv + 1) + i * (*lab)];
  }  
}

/**
 * Solve linear system Ax=b using preconditioned Richardson iteration.
 * The iteration is: x^(k+1) = x^(k) + M^(-1)*(b - A*x^(k))
 * where M is either D for Jacobi or (D-E) for Gauss-Seidel.
 * Stops when ||b - A*x^(k)||_2  / ||b||_2 < tol or when reaching maxit iterations.
 */
void richardson_MB(double *AB, double *RHS, double *X, double *MB, 
                   int *lab, int *la, int *ku, int *kl, 
                   double *tol, int *maxit, double *resvec, int *nbite) {
    
    int iter;
    double *res = (double *)malloc((*la) * sizeof(double));
    double *Ax = (double *)malloc((*la) * sizeof(double));
    double norm_r, norm_b, rel_res;
    
    int kv_local = (*lab) - (*kl) - (*ku) - 1;
    int diag_idx = kv_local + (*ku);
    int sub_idx = diag_idx + 1;
    
    int use_gs = 0;
    for (int i = 0; i < (*la) - 1; i++) {
        if (fabs(MB[sub_idx + i * (*lab)]) > 1e-15) {
            use_gs = 1;
            break;
        }
    }
    
    for (int i = 0; i < *la; i++) {
        X[i] = 0.0;
    }
    
    norm_b = cblas_dnrm2(*la, RHS, 1);
    if (norm_b < 1e-15) norm_b = 1.0;
    
    for (iter = 0; iter < *maxit; iter++) {
        cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku,
                    1.0, AB, *lab, X, 1, 0.0, Ax, 1);
        
        cblas_dcopy(*la, RHS, 1, res, 1);
        cblas_daxpy(*la, -1.0, Ax, 1, res, 1);
        
        norm_r = cblas_dnrm2(*la, res, 1);
        rel_res = norm_r / norm_b;
        
        if (resvec != NULL) {
            resvec[iter] = rel_res;
        }
        
        if (rel_res < *tol) {
            break;
        }
        
        if (use_gs == 0) {
            for (int i = 0; i < *la; i++) {
                double d = MB[diag_idx + i * (*lab)];
                if (fabs(d) > 1e-15) {
                    res[i] /= d;
                } else {
                    res[i] = 0.0;
                }
            }
        } else {
            for (int i = 0; i < *la; i++) {
                double sum = 0.0;
                double d = MB[diag_idx + i * (*lab)];
                
                if (i > 0) {
                    double s = MB[sub_idx + (i - 1) * (*lab)];
                    sum += s * res[i - 1];
                }
                
                res[i] = (res[i] - sum) / d;
            }
        }
        
        cblas_daxpy(*la, 1.0, res, 1, X, 1);
    }
    
    *nbite = iter + 1;
    
    free(res);
    free(Ax);
}