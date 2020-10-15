#include <stdlib.h>
#include <stdio.h>
#include "util_prox.h"
#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include <float.h>
#include <stdbool.h>

#ifndef Rpackage
#include "blas.h"
#include "lapack.h"
#else
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#endif

int     izero       =  0;
int     ione        =  1;
double  dzero       =  0.0;
double  done        =  1.0;
double  dminusone   = -1.0;

// for debugging //


void print_bmatrix(const bool* matrix,int m,int n){
    int i,j;
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            printf("%d, ", matrix[i*n+j]);
        }
        printf("\n");
    }
    printf("\n");
}

void print_dmatrix(const double* matrix,int m,int n){
    int i,j;
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            printf("%5.4f, ", matrix[i*n+j]);
        }
        printf("\n");
    }
    printf("\n");
}
// for debugging //


void dmat_vcopy(const int n, const double *src, double *dst)
{
    F77_CALL(dcopy)(&n, src, &ione, dst, &ione);
}

void dmat_vset(int n, const double val, double *dst)
{
    while (n-- != 0)
        *dst++ = val;
}

void dmat_waxpby(int n, double alpha, const double *x, double beta,
                 const double *y, double *w)

{
#if 1
    if (w != x && w != y)
    {
        dmat_vset(n, 0, w);
        F77_CALL(daxpy)(&n, &alpha, x, &ione, w, &ione);
        F77_CALL(daxpy)(&n, &beta , y, &ione, w, &ione);
    }
    else if (w == x && w == y)
    {
        double tmp;
        tmp = alpha+beta;
        F77_CALL(dscal)(&n, &tmp, w, &ione);
    }
    else if (w == x /*&& w != y */)
    {
        if (alpha != 1.0) F77_CALL(dscal)(&n, &alpha, w, &ione);
        if (beta  != 0.0) F77_CALL(daxpy)(&n, &beta , y, &ione, w, &ione);
    }
    else /* if (w == y && w != x ) */
    {
        if (beta  != 1.0) F77_CALL(dscal)(&n, &beta , w, &ione);
        if (alpha != 0.0) F77_CALL(daxpy)(&n, &alpha, x, &ione, w, &ione);
    }
#else
    int i;
    
    if (beta == 0.0)
    {
        if (alpha == -1.0)
        {
            for (i = 0; i < n; i++)
                w[i] = -x[i];
        }
        else
        {
            for (i = 0; i < n; i++)
                w[i] = alpha*x[i];
        }
    }
    else if (beta == 1.0)
    {
        if (alpha == -1.0)
        {
            for (i = 0; i < n; i++)
                w[i] = -x[i] + y[i];
        }
        else
        {
            for (i = 0; i < n; i++)
                w[i] = alpha*x[i] + y[i];
        }
    }
    else if (beta == -1.0)
    {
        if (alpha == -1.0)
        {
            for (i = 0; i < n; i++)
                w[i] = -x[i] - y[i];
        }
        else
        {
            for (i = 0; i < n; i++)
                w[i] = alpha*x[i] - y[i];
        }
    }
    else
    {
        for (i = 0; i < n; i++)
            w[i] = alpha*x[i] + beta*y[i];
    }
#endif
}


double dmat_norm2(const int n, const double *x)
{
    return F77_CALL(dnrm2)(&n, x, &ione);
}

void dmat_C_ATB(int m, int n1, int n2, double* beta,
                const double* A, const double* B, double* C)
{
    F77_CALL(dgemm)("N","T",&n2,&n1,&m,&done,B,&n2,A,&n1,beta,C,&n2);
}



void soft_threshold_gen(double *x, bool* thred, const int n, const double lam)
{
    int i;
    // #pragma omp parallel for
    for (i = 0; i < n; i++) {
        if (thred[i]) x[i] = fmax(0, x[i] - lam) - fmax(0, -x[i] - lam);
    }
}
void soft_threshold(double *x, const int n, const double lam)
{
    int i;
    // #pragma omp parallel for
    for (i = 0; i < n; i++) {
        x[i] = fmax(0, x[i] - lam) - fmax(0, -x[i] - lam);
    }
}



void eigen_decomp(int n, double* X, double *eigvec, double *eigval) {
    
    double *WORK;
    double abstol, WORKopt, vl, vu;
    int *IWORK;
    int numeig, sizeWORK, sizeIWORK, IWORKopt, il, iu,info;
    vl = 0.0;
    vu = 0.0;
    il = 0;
    iu = 0;
    /*  The support of the eigenvectors. We will not use this but the routine needs it  */
    int ISUPPZ[2*n];
    abstol = -1.0; // default tolerance
    
    /*  Query the Lapack routine for optimal sizes for workspace arrays  */
    sizeWORK = -1;
    sizeIWORK = -1;
    F77_CALL(dsyevr)("V", "A", "L", &n, X, &n, &vl,&vu,&il,&iu, &abstol, &numeig, eigval, eigvec, &n, ISUPPZ, &WORKopt, &sizeWORK, &IWORKopt, &sizeIWORK,&info);
    sizeWORK = (int)WORKopt;
    sizeIWORK = IWORKopt;
    WORK = (double*)malloc (sizeWORK*sizeof(double));
    IWORK = (int*)malloc (sizeIWORK*sizeof(int));
    /*  Now calculate the eigenvalues and vectors using optimal workspaces  */
    F77_CALL(dsyevr)("V", "A", "L", &n, X, &n, &vl,&vu,&il,&iu, &abstol, &numeig, eigval, eigvec, &n, ISUPPZ, WORK, &sizeWORK, IWORK, &sizeIWORK,&info);
    /*  Cleanup  */
    free((void*)(WORK)); free((void*)IWORK);
}



/* -------------------------------------------------------------------------- */
/* function solves minimize_x .5*sum c_i*(x_i-a_i)^2+sum b_i*|x_i|+lam*||x||_2*/
int eval_prox_gen_group(const double *a, const double *b, const double* c,
                        const double* lam, int* n, int* max_num_iter, double* x)
/* -------------------------------------------------------------------------- */
{
    int converge = 0;
    int i, j, max_iter = max_num_iter[0], p = n[0];
    double tmp_norm2,c_min=c[0], c_max=c[0], min, max, I,tmp1,tmp2;
    /* Allocate memory */
    double* tmp     = (double *) malloc(sizeof(double) * p);
    
    for (i = 0; i < p; i++){
        c_min = fmin(c[i], c_min);
        c_max = fmax(c[i], c_max);
    }
    for (i = 0; i < p; i++)
    { tmp[i] = c[i] * (fmax(0, a[i] - b[i]/c[i]) - fmax(0, -a[i] - b[i]/c[i])); }
    // print_dmatrix(tmp,1,p);
    tmp_norm2 = dmat_norm2(p, tmp);
    if (tmp_norm2 <= lam[0]) {
        for (i = 0; i < p; i++) { x[i] = .0; }
    } else {
        max = (tmp_norm2 - lam[0]) / c_min;
        min = (tmp_norm2 - lam[0]) / c_max;
        I = .5*(min + max);
        for (i = 0; i < max_iter; i++)
        {
            tmp1 = .0;
            for (j = 0; j < p; j++){
                tmp2 = tmp[j] / (c[j]*I + lam[0]);
                tmp1 += (tmp2 * tmp2);
            }
            if (tmp1 > 1.0) { min = I; } else { max = I; }
            if (abs(tmp1 - 1.0) < 1e-8) { converge = 1; break; }
        }
        I = .5*(min + max);
        for (i = 0; i < p; i++)
            x[i] = I * tmp[i] / (c[i]*I + lam[0]);
    }
    
    /* Deallocate memory */
    if (tmp != NULL) free(tmp);
    return converge;
}

void clime_sym(const double* Sigma_hat, const double* lam,
               const int* pp, const double* max_sig_val_inv,
               const int* num_iter, const double* eps,
               double* Rho, double* Omega)
{
    int p = pp[0];
    int max_iter = num_iter[0], p_square = p*p;
    int i,j;
    double rho = Rho[0], primal_res, dual_res;
    double* Gamma = (double *) malloc(sizeof(double) * p_square);
    double* Gamma_old = (double *) malloc(sizeof(double) * p_square);
    double* tmp1 = (double *) malloc(sizeof(double) * p_square);
    double* tmp2 = (double *) malloc(sizeof(double) * p_square);
    
    dmat_vset(p_square,.0,Gamma);
    dmat_vset(p_square,.0,Gamma_old);
    // dmat_vset(p_square,.0,tmp1);
    dmat_vset(p_square,.0,tmp2);
    for (i = 0; i < max_iter; i++) {
        
        // update Omega //
        dmat_waxpby(p_square,2.0,Gamma,-1.0,Gamma_old,tmp2);
        // make a copy of the previous Omega //
        dmat_vcopy(p_square,Omega,Gamma_old);
        dmat_C_ATB(p,p,p,&dzero,Sigma_hat,tmp2,tmp1);
        dmat_C_ATB(p,p,p,&done,tmp2,Sigma_hat,tmp1);
        dmat_waxpby(p_square,-.5*max_sig_val_inv[0],tmp1,1.0,Omega,Omega);
        for (j = 0; j < p; j++) { tmp1[j] = Omega[j*p+j]; }
        soft_threshold(Omega,p_square,max_sig_val_inv[0]/rho);
        for (j = 0; j < p; j++) { Omega[j*p+j] = tmp1[j]; }
        
        // calculate d(Omega,Omega_old) //
        dmat_waxpby(p_square,1.0,Omega,-1.0,Gamma_old,Gamma_old);
        dual_res = rho * dmat_norm2(p_square, Gamma_old) / p;
        
        // store Gamma to Gamma_old //
        dmat_vcopy(p_square,Gamma,Gamma_old);
        
        // update Gamma //
        dmat_C_ATB(p,p,p,&dzero,Sigma_hat,Omega,tmp1);
        dmat_C_ATB(p,p,p,&done,Omega,Sigma_hat,tmp1);
        dmat_waxpby(p_square,.5,tmp1,.0,tmp1,tmp1);
        for (j = 0; j < p; j++) { tmp1[j*p+j] -= 1.0; }
        dmat_waxpby(p_square,1.0,tmp1,1.0,Gamma,Gamma);
        soft_threshold(Gamma, p_square, lam[0]);
        
        // caluclate d(Gamma, Gamma_old)
        dmat_waxpby(p_square,1.0,Gamma,-1.0,Gamma_old,tmp2);
        primal_res = dmat_norm2(p_square, tmp2) / p;
        
        // printf("#iter: %d, primal residual: %f, dual residual: %f. \n", i, primal_res, dual_res);
        // stopping criterion //
        if (fmax(primal_res,dual_res) < eps[0]) {
            printf("converges in %d steps! \n", i);
            // printf("#iter: %d, primal residual: %f, dual residual: %f. \n", i, primal_res, dual_res);
            break;
        }
    }
    if (Gamma != NULL) free(Gamma);
    if (Gamma_old != NULL) free(Gamma_old);
    if (tmp1 != NULL) free(tmp1);
    if (tmp2 != NULL) free(tmp2);
}

void clime_sym_gen(const double* Sigma_hat, const double* lam,
                   const int* pp, const double* max_sig_val_inv,
                   const int* num_iter, const double* eps,
                   bool* thred,
                   double* Rho, double* Omega)
{
    int p = pp[0];
    int max_iter = num_iter[0], p_square = p*p;
    int i,j;
    double rho = Rho[0], primal_res, dual_res;
    double* Gamma = (double *) malloc(sizeof(double) * p_square);
    double* Gamma_old = (double *) malloc(sizeof(double) * p_square);
    double* tmp1 = (double *) malloc(sizeof(double) * p_square);
    double* tmp2 = (double *) malloc(sizeof(double) * p_square);
    
    dmat_vset(p_square,.0,Gamma);
    // update Gamma once if Omega is available for warm-start //
//    dmat_C_ATB(p,p,p,&dzero,Sigma_hat,Omega,tmp1);
//    dmat_C_ATB(p,p,p,&done,Omega,Sigma_hat,tmp1);
//    dmat_waxpby(p_square,.5,tmp1,.0,tmp1,tmp1);
//    for (j = 0; j < p; j++) { tmp1[j*p+j] -= 1.0; }
//    dmat_waxpby(p_square,1.0,tmp1,1.0,Gamma,Gamma);
//    soft_threshold(Gamma, p_square, lam[0]);
    
    dmat_vset(p_square,.0,Gamma_old);
    // dmat_vset(p_square,.0,tmp1);
    dmat_vset(p_square,.0,tmp2);
    for (i = 0; i < max_iter; i++) {
        
        // update Omega //
        dmat_waxpby(p_square,2.0,Gamma,-1.0,Gamma_old,tmp2);
        // make a copy of the previous Omega //
        dmat_vcopy(p_square,Omega,Gamma_old);
        dmat_C_ATB(p,p,p,&dzero,Sigma_hat,tmp2,tmp1);
        dmat_C_ATB(p,p,p,&done,tmp2,Sigma_hat,tmp1);
        dmat_waxpby(p_square,-.5*max_sig_val_inv[0],tmp1,1.0,Omega,Omega);
        for (j = 0; j < p; j++) { tmp1[j] = Omega[j*p+j]; }
        //soft_threshold(Omega,p_square,max_sig_val_inv[0]/rho);
        //for (j = 0; j < p; j++) { Omega[j*p+j] = tmp1[j]; }
        soft_threshold_gen(Omega,thred,p_square,max_sig_val_inv[0]/rho);
        // calculate d(Omega,Omega_old) //
        dmat_waxpby(p_square,1.0,Omega,-1.0,Gamma_old,Gamma_old);
        dual_res = rho * dmat_norm2(p_square, Gamma_old) / p;
        
        // store Gamma to Gamma_old //
        dmat_vcopy(p_square,Gamma,Gamma_old);
        
        // update Gamma //
        dmat_C_ATB(p,p,p,&dzero,Sigma_hat,Omega,tmp1);
        dmat_C_ATB(p,p,p,&done,Omega,Sigma_hat,tmp1);
        dmat_waxpby(p_square,.5,tmp1,.0,tmp1,tmp1);
        for (j = 0; j < p; j++) { tmp1[j*p+j] -= 1.0; }
        dmat_waxpby(p_square,1.0,tmp1,1.0,Gamma,Gamma);
        soft_threshold(Gamma, p_square, lam[0]);
        
        // caluclate d(Gamma, Gamma_old)
        dmat_waxpby(p_square,1.0,Gamma,-1.0,Gamma_old,tmp2);
        primal_res = dmat_norm2(p_square, tmp2) / p;
        
        // printf("#iter: %d, primal residual: %f, dual residual: %f. \n", i, primal_res, dual_res);
        // stopping criterion //
        if (fmax(primal_res,dual_res) < eps[0]) {
            printf("converges in %d steps! \n", i);
            print_dmatrix(Omega,p,p);
            print_dmatrix(Gamma,p,p);
            // printf("#iter: %d, primal residual: %f, dual residual: %f. \n", i, primal_res, dual_res);
            break;
        }
    }
    if (Gamma != NULL) free(Gamma);
    if (Gamma_old != NULL) free(Gamma_old);
    if (tmp1 != NULL) free(tmp1);
    if (tmp2 != NULL) free(tmp2);
}

void clime_sym_L1(const double* Sigma_hat, const double* lam,
                  const int* pp, const double* max_sig_val_inv,
                  const int* num_iter, const double* eps,
                  double* Rho, double* Omega)
{
    int i,p_square = pp[0]*pp[0];
    bool* thred = (bool *) malloc(sizeof(bool)*p_square);
    for (i = 0; i < p_square; i++) { thred[i] = true; }
    for (i = 0; i < pp[0]; i++) { thred[i*pp[0]+i] = false; }
    clime_sym_gen(Sigma_hat,lam,pp,max_sig_val_inv,num_iter,eps,thred,Rho,Omega);
    if (thred != NULL) free(thred);
}

// nonconvex clime does not work, dual variables always equal to zero! //
void clime_sym_trL1(const double* Sigma_hat, const double* lam,
                const double* tau,
                const int* pp, const double* max_sig_val_inv,
                const int* num_iter, const int* dc_max_iter,
                const double* eps,
                double* Rho, double* Omega)
{
    int i,j,p_square = pp[0]*pp[0];
    bool sparsity_pattern_changed;
    bool* thred = (bool *) malloc(sizeof(bool)*p_square);
    for (i = 0; i < p_square; i++) { thred[i] = true; }
    for (i = 0; i < pp[0]; i++) { thred[i*pp[0]+i] = false; }
    clime_sym_gen(Sigma_hat,lam,pp,max_sig_val_inv,num_iter,eps,thred,Rho,Omega);
    print_bmatrix(thred,pp[0],pp[0]);
    for (j = 0; j < dc_max_iter[0]; j++) {
        sparsity_pattern_changed = false;
        for (i = 0; i < p_square; i++) {
            if (thred[i] != (abs(Omega[i]) < tau[0])) {
                sparsity_pattern_changed = true;
                thred[i] = (abs(Omega[i]) < tau[0]);
            }
        }
        print_bmatrix(thred,pp[0],pp[0]);
        if (!sparsity_pattern_changed) { printf("DC iteration converges at %d. \n", j); break; }
        clime_sym_gen(Sigma_hat,lam,pp,max_sig_val_inv,num_iter,eps,thred,Rho,Omega);
    }
    if (thred != NULL) free(thred);
}



void clime_sym_path(const double* Sigma_hat, const double* lam,
                    const int* lam_length,
                    const int* pp, const double* max_sig_val_inv,
                    const int* num_iter, const double* eps,
                    double* Rho, double* Omega)
{
    int i, p = pp[0];
    // double* tmp = (double *) malloc(sizeof(double) * p * p);
    for (int i = 0; i < lam_length[0]; i++)
    {
        printf("lambda is: %f\n", lam[i]);
        if (i == 0) {
            clime_sym(Sigma_hat,lam,pp,max_sig_val_inv,num_iter,eps,Rho,Omega);
        } else {
            dmat_vcopy(p*p,Omega+(i-1)*p*p,Omega+i*p*p);
            clime_sym(Sigma_hat,lam+i,pp,max_sig_val_inv,num_iter,eps,Rho,Omega+i*p*p);
        }
    }
    
    // if (tmp != NULL) free(tmp);
}


void group_clime_sym(const double* Sigma_hat, const double* lam,
               const double* alpha, const int* num_iter,
               const int* pp, const int* KK,
               const double* max_sig_val_inv,
               const double* eps, double* Rho, double* Omega)
{
    int p = pp[0], K = KK[0];
    int max_iter = num_iter[0], p_square = p*p;
    int i,j,k,l,index;
    double rho = Rho[0], nu = 1.0-alpha[0], primal_res, dual_res,tmp;
    double* Gamma = (double *) malloc(sizeof(double) * p_square * K);
    double* Gamma_old = (double *) malloc(sizeof(double) * p_square * K);
    double* tmp1 = (double *) malloc(sizeof(double) * p_square);
    double* tmp2 = (double *) malloc(sizeof(double) * p_square);
    double* a = (double *) malloc(sizeof(double) * K);
    double* b = (double *) malloc(sizeof(double) * K);
    dmat_vset(K,alpha[0],b);
    double* c = (double *) malloc(sizeof(double) * K);
    for (k = 0; k < K; k++) { c[k] = rho/max_sig_val_inv[k]; }
    double* x = (double *) malloc(sizeof(double) * K);
    
    dmat_vset(p_square*K,.0,Gamma);
    dmat_vset(p_square*K,.0,Gamma_old);
    dmat_vset(p_square,.0,tmp2);
    for (i = 0; i < max_iter; i++) {
        
        // update Omega: store new Omega in Gamma_old //
        for (k = 0; k < K; k++) {
            dmat_waxpby(p_square,2.0,Gamma+k*p_square,-1.0,Gamma_old+k*p_square,tmp2);
            dmat_C_ATB(p,p,p,&dzero,Sigma_hat+k*p_square,tmp2,tmp1);
            dmat_C_ATB(p,p,p,&done,tmp2,Sigma_hat+k*p_square,tmp1);
            dmat_waxpby(p_square,-.5*max_sig_val_inv[k],tmp1,1.0,Omega+k*p_square,Gamma_old+k*p_square);
        }
        for (j = 1; j < p; j++){
            for (l = 0; l < j; l++) {
                index = j*p+l;
                for (k = 0; k < K; k++) { a[k] = Gamma_old[k*p_square+index]; }
                // print_dmatrix(a,1,K); print_dmatrix(b,1,K); print_dmatrix(c,1,K);
                eval_prox_gen_group(a,b,c,&nu,&K,&max_iter,x);
                // print_dmatrix(x,1,K);
                for (k = 0; k < K; k++) {
                    Gamma_old[k*p_square+index] = x[k];
                    Gamma_old[k*p_square+l*p+j] = x[k];
                }
            }
        }
        
        // calculate dual residual //
        dmat_waxpby(p_square*K,1.0,Omega,-1.0,Gamma_old,Omega);
        dual_res = rho * dmat_norm2(p_square*K,Omega) / (p * sqrt(K));
        
        // update Omega //
        dmat_vcopy(p_square*K, Gamma_old, Omega);
        
        // store Gamma to Gamma_old //
        dmat_vcopy(p_square*K,Gamma,Gamma_old);
        
        primal_res = .0;
        // update Gamma //
        for (k = 0; k < K; k++) {
            dmat_C_ATB(p,p,p,&dzero,Sigma_hat+k*p_square,Omega+k*p_square,tmp1);
            dmat_C_ATB(p,p,p,&done,Omega+k*p_square,Sigma_hat+k*p_square,tmp1);
            dmat_waxpby(p_square,.5,tmp1,.0,tmp1,tmp1);
            for (j = 0; j < p; j++) { tmp1[j*p+j] -= 1.0; }
            dmat_waxpby(p_square,1.0,tmp1,1.0,Gamma+k*p_square,Gamma+k*p_square);
            soft_threshold(Gamma+p_square*k, p_square, lam[0]);
            dmat_waxpby(p_square,1.0,Gamma_old+p_square*k,-1.0,Gamma+p_square*k,tmp1);
            tmp = dmat_norm2(p_square, tmp1); tmp = tmp * tmp;
            primal_res += tmp;
        }
        primal_res = sqrt(primal_res) / (p * sqrt(K));
        
        // printf("#iter: %d, primal residual: %f, dual residual: %f. \n", i, primal_res, dual_res);
        // stopping criterion //
        if (fmax(primal_res,dual_res) < eps[0]) {
            printf("converges in %d steps! \n", i);
            // printf("#iter: %d, primal residual: %f, dual residual: %f. \n", i, primal_res, dual_res);
            break;
        }
        
    }
    if (Gamma != NULL) free(Gamma); if (Gamma_old != NULL) free(Gamma_old);
    if (tmp1 != NULL) free(tmp1); if (tmp2 != NULL) free(tmp2);
    if (a != NULL) free(a); if (b != NULL) free(b);
    if (c != NULL) free(c); if (x != NULL) free(x);
}

void group_clime_sym_path(const double* Sigma_hat, const double* lam,
                     const int* lam_length,
                     const double* alpha, const int* alpha_length,
                     const int* num_iter,
                     const int* pp, const int* KK,
                     const double* max_sig_val_inv,
                     const double* eps, double* Rho, double* Omega)
{
    int i, j, p = pp[0], K = KK[0];
    for (i = 0; i < lam_length[0]; i++) {
        for (j = 0; j < alpha_length[0]; j++) {
            printf("lambda is: %f, %f. \n", lam[i], alpha[j]);
            if (i == 0 && j == 0){
                group_clime_sym(Sigma_hat,lam,alpha,num_iter,pp,KK,max_sig_val_inv,eps,Rho,Omega);
            }
            else if (j == 0){
                dmat_vcopy(p*p*K, Omega+p*p*K*(i-1)*alpha_length[0],Omega+p*p*K*i*alpha_length[0]);
                group_clime_sym(Sigma_hat,lam+i,alpha+j,num_iter,pp,KK,max_sig_val_inv,eps,Rho,Omega+p*p*K*(i*alpha_length[0]+j));
            } else {
                dmat_vcopy(p*p*K, Omega+p*p*K*(i*alpha_length[0]+j-1),Omega+p*p*K*(i*alpha_length[0]+j));
                group_clime_sym(Sigma_hat,lam+i,alpha+j,num_iter,pp,KK,max_sig_val_inv,eps,Rho,Omega+p*p*K*(i*alpha_length[0]+j));
            }
        }
    }
}

