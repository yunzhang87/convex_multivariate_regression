#ifndef DMATRIX_H
#define DMATRIX_H

#include <stdbool.h>


#ifdef __cplusplus
extern "C" {
#endif

    extern int     izero;
    extern int     ione;
    extern double  dzero;
    extern double  done;
    extern double  dminusone;
    
    
    void print_bmatrix(const bool* matrix,int m,int n);
    void print_dmatrix(const double* matrix,int m,int n);
    void dmat_elemprod(const int n, const double *x, const double *y, double *z);
    void dmat_vcopy(const int n, const double *src, double *dst);
    void dmat_vset(int n, const double val, double *dst);
    void dmat_waxpby(int n, double alpha, const double *x, double beta,
                     const double *y, double *w);
    double dmat_norm2(const int n, const double *x);
    double dmat_norminf(const int n, const double *x);
    void dmat_C_ATB(int m, int n1, int n2,
                    const double* A, const double* B, double* C);
    void dmat_C_ABT(int m, int n1, int n2, double* beta,
                    const double* A, const double* B, double* C);
    void dmat_B_AAT(int m, int n, double *A, double *B);
    void dmat_B_ADAT(int m, int n, double *A, double* D, double *B);
    void dmat_B_ATDA(int m, int n, double *A, double* D, double *B);
    void eigen_decomp(int n, double* X, double *eigvec, double *eigval);
    
    void soft_threshold_gen(double *x, bool* thred, const int n, const double lam);
    
    void soft_threshold(double *x, const int n, const double lam);
    void soft_threshold_vec(double *x, const int n, const double* lam_mat);
    
#ifdef __cplusplus
}
#endif

#endif /* DMATRIX_H */