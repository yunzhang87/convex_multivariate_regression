//
//  main.cpp
//  cGGM
//
//  Created by Yunzhang Zhu on 4/12/16.
//  Copyright © 2016 OSU. All rights reserved.
//

#include <iostream>
#include <random>
#include <math.h>

#include "test_liag.h"
#include "opt.h"
#include "dmatrix.h"
#include "def.h"

void test_cGGM();


int main(int argc, const char * argv[]) {
    // test_liag();
    test_cGGM();
    
    return 0;
}


void test_cGGM()
{
    int i, j;
    std::srand(1234);
    std::default_random_engine generator;
    std::normal_distribution<double> distribution(.0,1.0);
    int n = 50, p = 10, q = 10;
    double *X = new double[n*p];
    double *Y = new double[n*q];
    double *B = new double[p*q];
    double *E = new double[n*q];
    double *Omega = new double[q*q];
    double *eigvec_omega = new double[q*q];
    double *eigval_omega = new double[q];
    double *eigval_sqrt_sigma = new double[q];
    double *tmp1 = new double[q*q]; double *tmp2 = new double[n*q];
    double *tmp3 = new double[q*q];
    
    dmat_vset(q*q, .0, Omega);
    for (i = 0; i < q; i++) Omega[i*q + i] = 1.0;
    Omega[1] = .5; Omega[q+2] = .5;
    
    for (i = 0; i < n; i++) {
        for (j = 0; j < p; j++) {
            X[i*p + j] = distribution(generator);
        }
        for (j = 0; j < q; j++) {
            E[i*q + j] = distribution(generator);
        }
    }
    
    eigen_decomp(q, Omega, eigvec_omega, eigval_omega);
    for (i = 0; i < q; i++) eigval_sqrt_sigma[i] = 1.0 / sqrt(eigval_omega[i]);
    // tmp2 = E * U' * D^{-1/2} * U //
    dmat_uB_ATDA(q, q, eigvec_omega, eigval_sqrt_sigma, tmp3);
    dmat_C_AB(q, n, q, &done, E, tmp3, tmp2);
    
    dmat_vset(p*q, .0, B);
    for (i = 0; i < 2; i++) {
        for (j = 0; j < q; j++) {
            B[i*q+j] = distribution(generator);
        }
    }
    
    // print_dmatrix(B, p, q);
    // print_dmatrix(tmp2, n, q);
    
    // Y = X * B + E //
    dmat_C_AB(p, n, q, &done, X, B, Y);
    dmat_waxpby(n*q, 1.0, Y, 1.0, tmp2, Y);
    
    int nonconvex = 1;
    double lamB[2] = {1.5, 2.0};
    double lamOmega[2] = {.6, 1.0};
    double gamma = .01;
    int lamB_grid = 2;
    int lamOmega_grid = 2;
    double *sol_path = new double[lamB_grid*lamOmega_grid*(p+q)*q];
    double *sol_path_nc = new double[lamB_grid*lamOmega_grid*(p+q)*q];
    double *Z_ini = new double[(p+q)*q]; dmat_vset((p+q)*q, .0, Z_ini);
    for (i = 0; i < q; i++) Z_ini[(p+i)*q + i] = 1.0;
    
    // ADMM tolerance has to be set to be very small, otherwise, the line search would be problematic. // 
    double eps_abs = 1e-5, eps_rel=1e-6, newton_tol=1e-8, mu = 10.0, tau = 2.0, rho = 1.0, alpha = .4;
    int max_newton_iter = 50, max_admm_iter = 5e3, max_dc_iter = 10;
    
    cGGM(Y, X, Z_ini, lamB, lamOmega, &gamma, sol_path, sol_path_nc, &p, &q, &n, &lamB_grid, &lamOmega_grid, &eps_abs, &eps_rel, &newton_tol, &mu, &tau, &rho, &alpha, &max_newton_iter, &max_admm_iter, &max_dc_iter, &nonconvex);
    
    print_dmatrix(sol_path+(p+q)*q, p, q);
    printf("\n \n");
    print_dmatrix(sol_path+(p+q)*q+p*q, q, q);
    printf("\n \n \n \n ");
    print_dmatrix(sol_path_nc+(p+q)*q, p, q);
    printf("\n \n");
    print_dmatrix(sol_path_nc+(p+q)*q+p*q, q, q);
}