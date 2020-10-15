#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>


#ifndef Rpackage
#include "blas.h"
#include "lapack.h"
#else
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#endif

#include "dmatrix.h"


int     izero       =  0;
int     ione        =  1;
double  dzero       =  0.0;
double  done        =  1.0;
double  dminusone   = -1.0;

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

void dmat_elemprod(const int n, const double *x, const double *y, double *z)
{
    if (y != z) {
        F77_CALL(dcopy)(&n, y, &ione, z, &ione);
    }
    F77_CALL(dtbmv)("U", "N", "N", &n, &izero, x, &ione, z, &ione);
}


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

double dmat_norminf(const int n, const double *x)
{
    return fabs(x[F77_CALL(idamax)(&n, x, &ione)-1]);
}

void dmat_C_ATB(int m, int n1, int n2,
                const double* A, const double* B, double* C)
{
    F77_CALL(dgemm)("N","T",&n2,&n1,&m,&done,B,&n2,A,&n1,&dzero,C,&n2);
}

// not tested //
void dmat_C_ABT(int m, int n1, int n2, double* beta,
                const double* A, const double* B, double* C)
{
    F77_CALL(dgemm)("T","N",&n2,&n1,&m,&done,B,&n2,A,&n1,beta,C,&n2);
}

void dmat_B_AAT(int m, int n, double *A, double *B)
{
    F77_CALL(dsyrk)("L","T",&m,&n,&done,A,&n,&dzero,B,&m);
}

// not tested //
void dmat_B_ADAT(int m, int n, double *A, double* D, double *B)
{
    int i;
    double* Tmp = (double*) malloc(sizeof(double)*m*n);
    dmat_vset(m*n,.0,Tmp);
    // scale cols of A by D and stores the scaled matrix in Tmp //
    for (i = 0; i < n; i++) {
        // dmat_elemprod(n, A+i*n, D, Tmp+i*n);
        dmat_waxpby(m,D[i], A+i*m,.0, Tmp,Tmp+i*m);
    }
    dmat_C_ATB(m, n, n, A, Tmp, B);
    free(Tmp);
}

// not tested //
/* A n by m, D n by n */
//void dmat_B_ATDA(int m, int n, double *A, double* D, double *B)
//{
//    int i;
//    double* Tmp = (double*) malloc(sizeof(double)*m*n);
//    dmat_vset(m*n,.0,Tmp);
//    for (i = 0; i < n; i++) {
//        // dmat_elemprod(n, A+i*n, D, Tmp+i*n);
//        dmat_waxpby(m,D[i], A+i*m,.0,Tmp+i*m,Tmp+i*m);
//    }
//    // if (isnan_mat(A,n*m)) { printf("A is nan! \n"); }
//    // if (isnan_mat(Tmp,n*m)) { printf("Tmp is nan! \n"); }
//    dmat_C_ATB(n, m, m, A, Tmp, B);
//    // if (isnan_mat(B,m*m)) { printf("B is nan! \n"); }
//    free(Tmp);
//}

void dmat_B_ATDA(int m, int n, double *A, double* D, double *B)
{
    int i;
    double* Tmp = (double*) malloc(sizeof(double)*m*n);
    dmat_vset(m*n,.0,Tmp);
    for (i = 0; i < n; i++) {
        // dmat_elemprod(n, A+i*n, D, Tmp+i*n);
        dmat_waxpby(m,D[i], A+i*m,.0,Tmp+i*m,Tmp+i*m);
    }
    // if (isnan_mat(A,n*m)) { printf("A is nan! \n"); }
    // if (isnan_mat(Tmp,n*m)) { printf("Tmp is nan! \n"); }
    dmat_C_ATB(n, m, m, A, Tmp, B);
    // if (isnan_mat(B,m*m)) { printf("B is nan! \n"); }
    free(Tmp);
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


void soft_threshold_vec(double *x, const int n, const double* lam_mat)
{
    int i;
    // #pragma omp parallel for
    for (i = 0; i < n; i++) {
        x[i] = fmax(0, x[i] - lam_mat[i]) - fmax(0, -x[i] - lam_mat[i]);
    }
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
            if (fabs(tmp1 - 1.0) < 1e-8) { converge = 1; break; }
        }
        I = .5*(min + max);
        for (i = 0; i < p; i++)
            x[i] = I * tmp[i] / (c[i]*I + lam[0]);
    }
    
    /* Deallocate memory */
    if (tmp != NULL) free(tmp);
    return converge;
}



