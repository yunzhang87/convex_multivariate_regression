#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif
    
    
#define     max(x,y)                (((x)>(y))?(x):(y))
#define     min(x,y)                (((x)<(y))?(x):(y))
#define     abs(x  )                (((x) > 0)?(x):(-x))
#define     sign(x )                (((x) > 0)?(1):(-1))
#define     zero(x )                (((x) == 0)?(1):(0))
    
    
    enum penalty_types
    {
        L1,
        Truncated_L1,
        MCP
    };
    
    typedef struct
    {
        bool *BRowSparsityInd, *OmegaSparsityInd;
        double  *sigmaXX, *sigmaXY, *sigmaYY, *B, *Omega, *Sigma, *Z, *Delta, *H, *G,
        *Omega_tmp, *B_tmp, *Z_tmp, *Z_prime, *H_eig_vectors, *H_eig_values, 
        *Sigma_eig_vectors, *Sigma_eig_values, *Omega_eig_vectors, *Omega_eig_values,
        *Z_cur, *Z_prev, *Gamma_cur, *Gamma_prev, *Gamma_diff, *tmp1, *tmp2, *tmp3, *tmp_values;
        double pri_res, dual_res, obj_cur, obj_prev, duality_gap, newton_step_size;
        int num_newton_steps, dc_iter;
    } tmpvars;
    

    typedef struct
    {
        int p, q, n, lambdaB_grid, lambdaOmega_grid;
        double *Y, *X, *lambdaB, *lambdaOmega, *gamma,
        *fun_vals, *Z_ini, *sol_path, *sol_path_nc;
    } problem_data;


    typedef struct 
    {
        int newton_iter, admm_iter;
        double newton_tol, alpha, rho, eps_abs, eps_rel, mu, tau;
    } parameters;
    
    
    void print_dmatrix(const double* matrix,int m,int n);
    void print_bmatrix(const bool* matrix,int m,int n);
    void dmat_yinvx(const int n, const double *x, double *y);
    void dmat_vcopy(const int n, const double *src, double *dst);
    double dmat_dot(const int n, const double *x, const double *y);
    void dmat_vset(int n, const double val, double *dst);
    void dmat_bset(int n, const bool val, bool *dst);
    void dmat_waxpby(int n, double alpha, const double *x, double beta, const double *y, double *w);
    double dmat_norm1(const int n, const double *x);
    double dmat_norminf(const int n, const double *x);
    double dmat_norm2(const int n, const double *x);
    void dmat_C_AB(int n, int m1, int m2, double* beta, const double* A, const double* B, double* C);
    
    
    void dmat_C_ATB(int m, int n1, int n2, double* beta, double* A, double* B, double* C);
    // void dmat_C_ATB_new(int m, int n1, int n2, double* beta, double* A, double* B, double* C);
    
    
    void dmat_C_ABT(int n, int m1, int m2, double* beta, const double* A, const double* B, double* C);
    void dmat_B_ATA(int m, int n, double *A, double *B);
    void dmat_C_ATA(int m, int n, double* beta, double *A, double* alpha, double *C);
    void dmat_B_AAT(int m, int n, double *A, double *B);
    void dmat_B_ATDA(int m, int n, double *A, double* D, double *B);
    void dmat_potrf(double *A, const int m);
    bool dmat_ispd(const int n, double* A);
    void eigen_decomp(int n, double* X, double *eigvec, double *eigval);
    void proj_offdiag_InftyBall(double* x, const int n, const double* lam, const bool* ind);
    void make_sym(double* x, const int n);
    void dmat_diagscale(double *M, int m, int n, const double *dl, const int invl, const double *dr, const int invr);
    void dmat_uB_ATDA(int m, int n, double *A, double* D, double *B); 
    void dmat_uC_ATA(int m, int n, double* beta, double *A, double* alpha, double *C);
    
#ifdef __cplusplus
}
#endif
