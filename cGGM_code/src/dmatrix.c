#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>


#include "dmatrix.h"
#include "def.h"
#include "cblas.h"


#ifndef Rpackage
#include "blas.h"
#include "lapack.h"
#else
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#endif



// for debugging //

bool isnan_mat(const double* mat, int n){
    bool tmp = false;
    for (int i = 0; i < n; i++){
        if(isnan(mat[i])) { tmp = true; break; }
    }
    return tmp;
}

void print_bmatrix(const bool* matrix,int m,int n)
{
    int i,j;
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            printf("%d, ", matrix[i*n+j]);
        }
        printf("\n");
    }
    printf("\n");
}

void print_dmatrix(const double* matrix,int m,int n)
{
    int i,j;
//    #ifndef Rpackage
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            printf("%5.4f, ", matrix[i*n+j]);
        }
        printf("\n");
    }
    printf("\n");
//    #else
//    for (i = 0; i < m; i++) {
//        for (j = 0; j < n; j++) {
//            Rprintf("%5.4f, ", matrix[i*n+j]);
//        }
//        Rprintf("\n");
//    }
//    Rprintf("\n");
//    #endif
}



// end for debugging //

/** \brief \f$ x^Ty \f$
 *
 *  Returns dot product of a vector x and y.
 *
 *  @param  n   length of a vector x.
 *  @param  x   pointer to a vector x.
 *  @param  y   pointer to a vector y.
 *  @return     result.
 */
double dmat_dot(const int n, const double *x, const double *y)
{
    return F77_CALL(ddot)(&n, x, &ione, y, &ione);
}



void dmat_yinvx(const int n, const double *x, double *y)
{
#ifdef MKL
    vdInv(n, x, y);
#else
    const double *xi;
    double *yi;
    xi = x+n-1;
    yi = y+n-1;
    do {
        *yi-- = 1/(*xi--);
    } while (xi >= x);
#endif
}


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

void dmat_bset(int n, const bool val, bool *dst)
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

/** \brief \f$ \|x\|_1 \f$
 *
 *  Returns 1-norm of a vector x.
 *
 *  @param  n   length of a vector x.
 *  @param  x   pointer to a vector x.
 *  @return     result.
 */
double dmat_norm1(const int n, const double *x)
{
    return F77_CALL(dasum)(&n, x, &ione);
}


double dmat_norminf(const int n, const double *x)
{
    return fabs(x[F77_CALL(idamax)(&n, x, &ione)-1]);
}

double dmat_norm2(const int n, const double *x)
{
    return F77_CALL(dnrm2)(&n, x, &ione);
}



/** \brief \f$ C = beta * AB \f$
 *
 *  Computes dense matrix-transpose-matrix product.
 *
 *  @param  A       pointer to m1 by n matrix.
 *  @param  B       pointer to n by m2 matrix.
 *  @param  C       pointer to m1 by m2 matrix.
 */
void dmat_C_AB(int n, int m1, int m2, double* beta, const double* A, const double* B, double* C)
{
    F77_CALL(dgemm)("N","N",&m2,&m1,&n,beta,B,&m2,A,&n,&dzero,C,&m2);
}


/** \brief \f$ C = beta * AB^T \f$
 *
 *  Computes dense matrix-transpose-matrix product.
 *
 *  @param  A       pointer to m1 by n matrix.
 *  @param  B       pointer to m2 by n matrix.
 *  @param  C       pointer to m1 by m2 matrix.
 */
void dmat_C_ABT(int n, int m1, int m2, double* beta,
                const double* A, const double* B, double* C)
{
    F77_CALL(dgemm)("T","N",&m2,&m1,&n,beta,B,&n,A,&n,&dzero,C,&m2);
}

/** \brief \f$ C = beta * A^TB \f$
 *
 *  Computes dense matrix-transpose-matrix product.
 *
 *  @param  A       pointer to m by n1 matrix.
 *  @param  B       pointer to m by n2 matrix.
 *  @param  C       pointer to n1 by n2 matrix.
 */
void dmat_C_ATB(int m, int n1, int n2, double* beta, double* A, double* B, double* C)
{
    F77_CALL(dgemm)("N","T",&n2,&n1,&m,beta,B,&n2,A,&n1,&dzero,C,&n2);
}

// different implementation using cblas. TO DO: test their differences! //
// void dmat_C_ATB_new(int m, int n1, int n2, double* beta, double* A, double* B, double* C)
// {
//     cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n1, n2, m, dzero, A, n1, B, n2, beta[0], C, n2);
// }


/** \brief \f$ B = A^TA \f$
 *
 *  Computes dense matrix-transpose-matrix product; return a
 *  row-oriented upper-diag symmetric matrix.
 *
 *  @param  A       pointer to m by n matrix.
 *  @param  B       pointer to a result matrix.
 */
void dmat_B_ATA(int m, int n, double *A, double *B)
{
    F77_CALL(dsyrk)("L","N",&n,&m,&done,A,&n,&dzero,B,&n);
}


void dmat_C_ATA(int m, int n, double* beta, double *A, double* alpha, double *C)
{
    F77_CALL(dsyrk)("L","N",&n,&m,beta,A,&n,alpha,C,&n);
}


/** \brief \f$ B = AA^T \f$
 *
 *  Computes dense matrix-transpose-matrix product; return a
 *  row-oriented upper-diag symmetric matrix.
 *
 *  @param  A       pointer to m by n matrix.
 *  @param  B       pointer to a result matrix.
 */
void dmat_B_AAT(int m, int n, double *A, double *B)
{
    F77_CALL(dsyrk)("L","T",&m,&n,&done,A,&n,&dzero,B,&m);
}


/** \brief \f$ B = A^TDA \f$
 *
 *  Computes dense matrix-transpose-matrix product.
 *
 *  @param  A       pointer to m by n matrix.
 *  @param  D       pointer to a m dimensional vector.
 *  @param  B       pointer to a result matrix.
 */
void dmat_B_ATDA(int m, int n, double *A, double* D, double *B)
{
    int i;
    double* Tmp = (double*) malloc(sizeof(double)*m*n);
    // scale rows of A by D and stores the scaled matrix in Tmp //
    for (i = 0; i < m; i++) {
        dmat_waxpby(n, sqrt(D[i]), A+i*n, 0.0, Tmp, Tmp+i*n);
    }
    dmat_B_ATA(m, n, Tmp, B);
    free(Tmp);
}


void proj_L1_gen(double *x, const double *c, const double *z, const int* nn)
{
    int n = nn[0];
    double f_0 = -z[0];
    int converge = 1;
    for (int i = 0; i < n; ++i)
    {
        f_0 = f_0 + c[i] * fabs(x[i]);
    }
    if (f_0 > .0)
    {
        converge = 0;
        double lambda_min = .0, lambda_max = .0, obj = .0, new_lam = .0;
        
        for (int i = 0; i < n; ++i)
        {
            lambda_max = fmax(lambda_max , fabs(x[i]));
        }
        
        for (int i = 0; i < 1e4; ++i)
        {
            new_lam = .5 * (lambda_max + lambda_min);
            // printf("current solution is %f at iteration %d. \n", new_lam, i);
            obj = -z[0];
            for (int j = 0; j < n; ++j)
            {
                obj += c[j] * fmax(fabs(x[j]) - c[j]*new_lam, .0);
            }
            if (obj > 0)
            {
                lambda_min = new_lam;
            } else
            {
                lambda_max = new_lam;
            }
            if (fabs(lambda_max - lambda_min)/(lambda_max+1e-10) < 1e-10)
            {
                // printf("converges in %d steps. \n", i);
                converge = 1;
                break;
            }
        }
        for (int i = 0; i < n; ++i)
        {
            x[i] = sign(x[i]) * fmax(fabs(x[i]) - c[i]*lambda_max, .0);
        }
    }
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

void proj_offdiag_InftyBall(double* x, const int n, const double* lam, const bool* ind)
{
    int i, j;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (!ind[i*n+j]) { x[i*n + j] = .0; continue; }
            if (x[i*n + j] > lam[0]) x[i*n + j] = lam[0];
            else if (x[i*n + j] < -lam[0]) x[i*n + j] = -lam[0];
        }
    }
}

void make_sym(double* x, const int n)
{
    int i, j; double tmp;
    for (i = 1; i < n; i++) {
        for (j = 0; j < i; j++) {
            tmp = .5 * (x[i*n + j] + x[j*n + i]);
            x[i*n + j] = tmp; x[j*n + i] = tmp; 
        }
    }
}

/*
 computes the Cholesky factorization of a real symmetric
 positive definite matrix A.
 */
void dmat_potrf(double *A, const int m)
{
    int info;
    F77_CALL(dpotrf)("L", &m, A, &m, &info);
}

bool dmat_ispd(const int n, double* A)
{
    int info;
    F77_CALL(dpotrf)("L", &n, A, &n, &info);
    if (info == 0) return true; else return false;
}



/*    X = eigvec' * eigval * eigvec    */
// warning: it changes the input matrix as well // 
void eigen_decomp(int n, double* X, double *eigvec, double *eigval)
{
    
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



/** \brief \f$ M = D_L M D_R \f$
 *
 *  Computes left and right diagonal scaleing of M.
 \f{eqnarray*}
 M &=& D_L M \quad\mbox{ ,if dr = NULL }\\
 M &=& M D_R \quad\mbox{ ,if dl = NULL }\\
 M &=& D_L M D_R \quad\mbox{ ,otherwise }
 \f}
 *
 *  @param  M       pointer to a matrix.
 *  @param  dl      pointer to a diagonal vector.
 *  @param  invl    inverse flag of left diagonal matrix.
 \f{eqnarray*}
 D_L &=& \mbox{diag}(dl)    \quad\mbox{ ,if invl} = 0\\
 D_L &=& \mbox{diag}(1./dl) \quad\mbox{ ,if invl} = 1
 \f}
 *  @param  dr      pointer to a diagonal vector.
 *  @param  invr    inverse flag of right diagonal matrix.
 \f{eqnarray*}
 D_R &=& \mbox{diag}(dr)    \quad\mbox{ ,if invr} = 0\\
 D_R &=& \mbox{diag}(1./dr) \quad\mbox{ ,if invr} = 1
 \f}
 */
void dmat_diagscale(double *M, int m, int n, const double *dl, const int invl,
                    const double *dr, const int invr)
{
    int i, j;
    double *val = M;
    
    if (dl != NULL && dr != NULL)
    {
        if (invl == FALSE && invr == FALSE)
        {
            for (i = 0; i < m; i++)
                for (j = 0; j < n; j++)
                    val[i*n+j] *= dl[i]*dr[j];
        }
        else if (invl == FALSE && invr == TRUE)
        {
            for (i = 0; i < m; i++)
                for (j = 0; j < n; j++)
                    val[i*n+j] *= dl[i]/dr[j];
        }
        else if (invl == TRUE && invr == FALSE)
        {
            for (i = 0; i < m; i++)
                for (j = 0; j < n; j++)
                    val[i*n+j] *= dr[j]/dl[i];
        }
        else
        {
            for (i = 0; i < m; i++)
                for (j = 0; j < n; j++)
                    val[i*n+j] /= dl[i]/dr[j];
        }
    }
    else if (dl != NULL)
    {
        if (invl == FALSE)
        {
            for (i = 0; i < m; i++)
                for (j = 0; j < n; j++)
                    val[i*n+j] *= dl[i];
        }
        else
        {
            for (i = 0; i < m; i++)
                for (j = 0; j < n; j++)
                    val[i*n+j] /= dl[i];
        }
    }
    else if (dr != NULL)
    {
        if (invr == FALSE)
        {
            for (i = 0; i < m; i++)
                for (j = 0; j < n; j++)
                    val[i*n+j] *= dr[j];
        }
        else
        {
            for (i = 0; i < m; i++)
                for (j = 0; j < n; j++)
                    val[i*n+j] /= dr[j];
        }
        
    }
}


/** \brief \f$ B = A^TDA \f$
 *
 *  Computes dense matrix-transpose-matrix product.
 *
 *  @param  A       pointer to m by n matrix.
 *  @param  D       pointer to a m dimensional vector.
 *  @param  B       pointer to a (nonsymmetric) result matrix.
 */
void dmat_uB_ATDA(int m, int n, double *A, double* D, double *B)
{
    double* Tmp = (double*) malloc(sizeof(double)*m*n);
    dmat_vcopy(m*n, A, Tmp);
    dmat_diagscale(Tmp, m, n, D, 0, NULL, 0);
    dmat_C_ATB(m, n, n, &done, A, Tmp, B);
    free(Tmp);
}


// C = alpha * C + beta * A' * A, where both upper and lower diagonals of C are stored //
void dmat_uC_ATA(int m, int n, double* beta, double *A, double* alpha, double *C)
{
    double* Tmp = (double*) malloc(sizeof(double)*n*n);
    dmat_C_ATB(m, n, n, beta, A, A, Tmp);
    dmat_waxpby(n*n, done, Tmp, alpha[0], C, C);
    free(Tmp);
    make_sym(C, n);
}












