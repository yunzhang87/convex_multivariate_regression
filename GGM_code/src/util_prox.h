
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
    
    
    void group_clime_sym(const double* Sigma_hat, const double* lam,
                         const double* alpha, const int* num_iter,
                         const int* pp, const int* KK,
                         const double* max_sig_val_inv,
                         const double* eps, double* Rho, double* Omega);
    
    
    void clime_sym(const double* Sigma_hat, const double* lam,
                   const int* pp, const double* max_sig_val_inv,
                   const int* num_iter, const double* eps,
                   double* Rho, double* Omega);
    
    
#ifdef __cplusplus
}
#endif