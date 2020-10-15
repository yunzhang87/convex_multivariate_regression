//
//  opt.h
//  cGGM
//
//  Created by Yunzhang Zhu on 4/18/16.
//  Copyright © 2016 OSU. All rights reserved.
//

#ifndef opt_h
#define opt_h


#ifdef __cplusplus
extern "C" {
#endif
    
void cGGM(double* Y, double*X, double* Z_ini, double* lambda_B,
          double* lambda_Omega, double* gamma, double* solution_path,
          double* solution_path_nc, int* p, int* q, int* n,
          int* lambdaB_grid, int* lambdaOmega_grid, double* eps_abs,
          double* eps_rel, double* newton_tol, double* mu, double* tau,
          double* rho, double* alpha, int* max_newton_iter,
          int* max_admm_iter, int* max_dc_iter, int* nonconvex);

    
#ifdef __cplusplus
}
#endif
    
#endif /* opt_h */
