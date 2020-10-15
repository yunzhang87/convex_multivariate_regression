#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#include "dmatrix.h"
#include "def.h"
#include "blas.h"
#include "lapack.h"


#ifdef __cplusplus
extern "C" {
#endif

    void form_H_and_G(problem_data* data, tmpvars* tmp);
    void proximal_newton_step(problem_data* data, parameters* para,tmpvars* tmp, double lamB, double lamOmega);
    void proximal_newton_step_nonconvex(problem_data* data, parameters* para, tmpvars* tmp, double lamB, double lamOmega, double tau);
    void line_search(problem_data* data, parameters* para,tmpvars* tmp, double lamB, double lamOmega);
    double penalty(int p, int q, double* B, double* Omega, double lamB, double lamOmega, bool* BInd, bool* OmegaInd);
    double obj_cGGM(problem_data* data, parameters* para, tmpvars* tmp, double* sol, double lamB, double lamOmega);
    void free_tmpvars(tmpvars* tmp, int* nonconvex);

    void ini_problem(double* Y, double*X, double* Z_ini, double* lambdaB,
                     double* lambdaOmega, double* gamma, double* solution_path,
                     double* solution_path_nc, int* p, int* q, int* n,
                     int* lambdaB_grid, int* lambdaOmega_grid, double* eps_abs,
                     double* eps_rel, double* newton_tol, double* mu, double* tau,
                     double* rho, double* alpha, int* max_newton_iter, int* max_admm_iter,
                     int* max_dc_iter, int* nonconvex, tmpvars** tmp_,
                     problem_data** data_, parameters** para_);



    // initialization must be provided, convex solutions will not be returned. //
    void cGGM_nonconvex(double* Y, double*X, double* Z_ini, double* lambda_B,
              double* lambda_Omega, double* gamma, double* initialization_path,
              double* solution_path_nc, int* p, int* q, int* n,
              int* lambdaB_grid, int* lambdaOmega_grid, double* eps_abs,
              double* eps_rel, double* newton_tol, double* mu, double* tau,
              double* rho, double* alpha, int* max_newton_iter,
              int* max_admm_iter, int* max_dc_iter, int* nonconvex)
    {
        tmpvars* tmp;
        problem_data* data;
        parameters* para;

        //    print_dmatrix(Y, n[0], q[0]);
        //    print_dmatrix(X, n[0], p[0]);
        //    print_dmatrix(Z_ini, p[0] + q[0], q[0]);

        ini_problem(Y, X, Z_ini, lambda_B, lambda_Omega, gamma, initialization_path,
                    solution_path_nc, p, q, n, lambdaB_grid, lambdaOmega_grid,
                    eps_abs,  eps_rel,  newton_tol,  mu,  tau, rho,
                    alpha, max_newton_iter, max_admm_iter, max_dc_iter,
                    nonconvex, &tmp, &data, &para);

        int qtq = data->q * data->q, ptq = data->p * data->q;
        int ppq = data->p + data->q, ppq_t_q = ppq * data->q;
        int i, j, k, m;
        double diff;

        // printf("Initial solution is: \n");
        // print_dmatrix(tmp->Z, ppq, data->q);
        for (i = 0; i < lambdaB_grid[0]; i++)
        {
            for (j = 0; j < lambdaOmega_grid[0]; j++)
            {
                // printf("\n \n lambdaB is %f and lambdaOmega is %f. \n", lambda_B[i], lambda_Omega[j]);

                // warm start solution from provided initilization //
                dmat_vcopy(ppq_t_q, initialization_path + ppq_t_q * (i*lambdaOmega_grid[0] + j), tmp->Z);

                // printf("initial solutions is: \n");
                // print_dmatrix(tmp->Z, ppq, data->q);
                for (k = 0; k < tmp->dc_iter; k++)
                {
                    // printf("\n DC iteration #%d starts...... \n", k+1);
                    if (k == 0) {
                        dmat_bset(data->p, true, tmp->BRowSparsityInd);
                        dmat_bset(data->q*data->q, true, tmp->OmegaSparsityInd);
                        for (int s = 0; s < data->q; s++) tmp->OmegaSparsityInd[s*data->q+s] = false;
                    } else {
                        for (int s = 0; s < data->p; s++)
                        {
                            tmp->BRowSparsityInd[s] = (dmat_norm2(data->q, tmp->Z+s*data->q) < sqrt(data->q) * gamma[0]);
                        }
                        for (int s = 0; s < data->q; s++) {
                            for (int t = 0; t < data->q; t++) {
                                if (s != t) tmp->OmegaSparsityInd[s*data->q+t] = (abs(tmp->Z[(data->p+s)*data->q+t]) < gamma[0]);
                            }
                        }
                    }

                    // print_bmatrix(tmp->BRowSparsityInd, 1, data->p);
                    // print_bmatrix(tmp->OmegaSparsityInd, data->q, data->q);

                    for (m = 0; m < para->newton_iter; m++)
                    {
                        // printf("\n Newton iteration #%d starts...... \n", m+1);

                        form_H_and_G(data, tmp);

                        //                print_dmatrix(tmp->Z, ppq, data->q);
                        //                print_dmatrix(tmp->sigmaXX, data->p, data->p);
                        //                print_dmatrix(tmp->sigmaXY, data->p, data->q);
                        //                print_dmatrix(tmp->sigmaYY, data->q, data->q);
                        //                print_dmatrix(tmp->H, ppq, ppq);
                        //                print_dmatrix(tmp->G, ppq, data->q);
                        //                print_dmatrix(tmp->Sigma, data->q, data->q);

                        proximal_newton_step(data,para,tmp,lambda_B[i], lambda_Omega[j]);


                        // print newton search direction //
                        // printf("Newton search direction is: \n");
                        // print_dmatrix(tmp->Delta, ppq, q[0]);

                        // compute delta, where delta should be negative. //
                        dmat_waxpby(ppq_t_q, 1.0, tmp->Z, 1.0, tmp->Delta, tmp->Z_tmp);
                        diff = penalty(p[0], q[0], tmp->Z_tmp, tmp->Z_tmp + ptq, data->lambdaB[i], data->lambdaOmega[j], tmp->BRowSparsityInd, tmp->OmegaSparsityInd) - penalty(p[0], q[0], tmp->Z, tmp->Z + ptq, data->lambdaB[i], data->lambdaOmega[j], tmp->BRowSparsityInd, tmp->OmegaSparsityInd);
                        diff = diff + dmat_dot(qtq, tmp->Sigma, tmp->Delta + ptq) - dmat_dot(ppq_t_q, tmp->G, tmp->Delta);

                        // printf("delta is: %10.20f. \n", diff);

                        if (diff > -1e-10) {
                            // printf("warning: delta is close to positive! \n");
                            // dmat_vcopy(ppq_t_q, tmp->Z, data->sol_path + (i*lambdaOmega_grid[0] + j)*ppq_t_q);
                            // printf("Newton step converges at %d-th iterations! \n \n", m+1);
                            break;
                        }

                        // perform line search for proximal newton method //
                        line_search(data, para, tmp, data->lambdaB[i], data->lambdaOmega[j]);

                        //                printf("solution is: \n");
                        //                print_dmatrix(tmp->Z, ppq, data->q);

                        // printf("objective diff is: %10.20f \n", tmp->obj_cur - tmp->obj_prev);
                        // printf("newton termination tolerence: %10.20f \n", para->newton_tol);

                        // stopping criterion //
                        if (abs(tmp->obj_prev - tmp->obj_cur) < para->newton_tol) {
                            // dmat_vcopy(ppq_t_q, tmp->Z, data->sol_path + (i*lambdaOmega_grid[0] + j)*ppq_t_q);
                            // printf("Newton step converges at %d-th iterations! \n \n", m+1);
                            break;
                        }

                        if (m == para->newton_iter - 1) {
                            // dmat_vcopy(ppq_t_q, tmp->Z, data->sol_path + (i*lambdaOmega_grid[0] + j)*ppq_t_q);
                            // printf("warning: proximal newton not yet converges! Please increase newton iteration number. \n");
                        }
                    }
                    // store convex solution //
                    if (k==0) dmat_vcopy(ppq_t_q, tmp->Z, data->sol_path + (i*lambdaOmega_grid[0] + j)*ppq_t_q);
                    if (nonconvex[0] == 0) { break; }
                }
                // store nonconvex solution //
                if (nonconvex[0] == 1) dmat_vcopy(ppq_t_q, tmp->Z, data->sol_path_nc + (i*lambdaOmega_grid[0] + j)*ppq_t_q);
            }
        }
        free_tmpvars(tmp, nonconvex);
    }


    void cGGM(double* Y, double*X, double* Z_ini, double* lambda_B,
              double* lambda_Omega, double* gamma, double* solution_path,
              double* solution_path_nc, int* p, int* q, int* n,
              int* lambdaB_grid, int* lambdaOmega_grid, double* eps_abs,
              double* eps_rel, double* newton_tol, double* mu, double* tau,
              double* rho, double* alpha, int* max_newton_iter,
              int* max_admm_iter, int* max_dc_iter, int* nonconvex)
    {
        tmpvars* tmp;
        problem_data* data;
        parameters* para;

        //    print_dmatrix(Y, n[0], q[0]);
        //    print_dmatrix(X, n[0], p[0]);
        //    print_dmatrix(Z_ini, p[0] + q[0], q[0]);

        ini_problem(Y, X, Z_ini, lambda_B, lambda_Omega, gamma, solution_path,
                    solution_path_nc, p, q, n, lambdaB_grid, lambdaOmega_grid,
                    eps_abs,  eps_rel,  newton_tol,  mu,  tau, rho,
                    alpha, max_newton_iter, max_admm_iter, max_dc_iter,
                    nonconvex, &tmp, &data, &para);

        int qtq = data->q * data->q, ptq = data->p * data->q;
        int ppq = data->p + data->q, ppq_t_q = ppq * data->q;
        int i, j, k, m;
        double diff;

        // printf("Initial solution is: \n");
        // print_dmatrix(tmp->Z, ppq, data->q);
        for (i = 0; i < lambdaB_grid[0]; i++)
        {
            for (j = 0; j < lambdaOmega_grid[0]; j++)
            {
                // printf("\n \n lambdaB is %f and lambdaOmega is %f. \n", lambda_B[i], lambda_Omega[j]);

                // warm start solution from previous grid //
                if (i > 0 && j == 0) {
                    dmat_vcopy(ppq_t_q, solution_path + ppq_t_q * ((i-1)*lambdaOmega_grid[0]), tmp->Z);
                }
                if (j > 0) {
                    dmat_vcopy(ppq_t_q, solution_path + ppq_t_q * (i*lambdaOmega_grid[0] + j-1), tmp->Z);
                }

                //            printf("initial solutions is: \n");
                //            print_dmatrix(solution_path, ppq, data->q);
                for (k = 0; k < tmp->dc_iter; k++)
                {
                    // printf("\n DC iteration #%d starts...... \n", k+1);
                    if (k == 0) {
                        dmat_bset(data->p, true, tmp->BRowSparsityInd);
                        dmat_bset(data->q*data->q, true, tmp->OmegaSparsityInd);
                        for (int s = 0; s < data->q; s++) tmp->OmegaSparsityInd[s*data->q+s] = false;
                    } else {
                        for (int s = 0; s < data->p; s++)
                        {
                            tmp->BRowSparsityInd[s] = (dmat_norm2(data->q, tmp->Z+s*data->q) < sqrt(data->q) * gamma[0]);
                        }
                        for (int s = 0; s < data->q; s++) {
                            for (int t = 0; t < data->q; t++) {
                                if (s != t) tmp->OmegaSparsityInd[s*data->q+t] = (abs(tmp->Z[(data->p+s)*data->q+t]) < gamma[0]);
                            }
                        }
                    }

//                    print_bmatrix(tmp->BRowSparsityInd, 1, data->p);
//                    print_bmatrix(tmp->OmegaSparsityInd, data->q, data->q);

                    for (m = 0; m < para->newton_iter; m++)
                    {
                        // printf("\n Newton iteration #%d starts...... \n", m+1);

                        form_H_and_G(data, tmp);

                        //                print_dmatrix(tmp->Z, ppq, data->q);
                        //                print_dmatrix(tmp->sigmaXX, data->p, data->p);
                        //                print_dmatrix(tmp->sigmaXY, data->p, data->q);
                        //                print_dmatrix(tmp->sigmaYY, data->q, data->q);
                        //                print_dmatrix(tmp->H, ppq, ppq);
                        //                print_dmatrix(tmp->G, ppq, data->q);
                        //                print_dmatrix(tmp->Sigma, data->q, data->q);

                        proximal_newton_step(data,para,tmp,lambda_B[i], lambda_Omega[j]);


                        // print newton search direction //
                        // printf("Newton search direction is: \n");
                        // print_dmatrix(tmp->Delta, ppq, q[0]);

                        // compute delta, where delta should be negative. //
                        dmat_waxpby(ppq_t_q, 1.0, tmp->Z, 1.0, tmp->Delta, tmp->Z_tmp);
                        diff = penalty(p[0], q[0], tmp->Z_tmp, tmp->Z_tmp + ptq, data->lambdaB[i], data->lambdaOmega[j], tmp->BRowSparsityInd, tmp->OmegaSparsityInd) - penalty(p[0], q[0], tmp->Z, tmp->Z + ptq, data->lambdaB[i], data->lambdaOmega[j], tmp->BRowSparsityInd, tmp->OmegaSparsityInd);
                        diff = diff + dmat_dot(qtq, tmp->Sigma, tmp->Delta + ptq) - dmat_dot(ppq_t_q, tmp->G, tmp->Delta);

                        // printf("delta is: %10.20f. \n", diff);

                        if (diff > -1e-10) {
                            // printf("warning: delta is close to positive! \n");
                            // dmat_vcopy(ppq_t_q, tmp->Z, data->sol_path + (i*lambdaOmega_grid[0] + j)*ppq_t_q);
                            // printf("Newton step converges at %d-th iterations! \n \n", m+1);
                            break;
                        }

                        // perform line search for proximal newton method //
                        line_search(data, para, tmp, data->lambdaB[i], data->lambdaOmega[j]);

                        //                printf("solution is: \n");
                        //                print_dmatrix(tmp->Z, ppq, data->q);

                        // printf("objective diff is: %10.20f \n", tmp->obj_cur - tmp->obj_prev);
                        // printf("newton termination tolerence: %10.20f \n", para->newton_tol);

                        // stopping criterion //
                        if (abs(tmp->obj_prev - tmp->obj_cur) < para->newton_tol) {
                            // dmat_vcopy(ppq_t_q, tmp->Z, data->sol_path + (i*lambdaOmega_grid[0] + j)*ppq_t_q);
                            // printf("Newton step converges at %d-th iterations! \n \n", m+1);
                            break;
                        }

                        if (m == para->newton_iter - 1) {
                            // dmat_vcopy(ppq_t_q, tmp->Z, data->sol_path + (i*lambdaOmega_grid[0] + j)*ppq_t_q);
                            // printf("warning: proximal newton not yet converges! Please increase newton iteration number. \n");
                        }
                    }
                    // store convex solution //
                    if (k==0) dmat_vcopy(ppq_t_q, tmp->Z, data->sol_path + (i*lambdaOmega_grid[0] + j)*ppq_t_q);
                }
                // store nonconvex solution //
                if (nonconvex[0] == 1) dmat_vcopy(ppq_t_q, tmp->Z, data->sol_path_nc + (i*lambdaOmega_grid[0] + j)*ppq_t_q);
            }
        }
        free_tmpvars(tmp, nonconvex);
    }

    void ini_problem(double* Y, double*X, double* Z_ini, double* lambdaB,
                     double* lambdaOmega, double* gamma, double* solution_path,
                     double* solution_path_nc, int* p, int* q, int* n,
                     int* lambdaB_grid, int* lambdaOmega_grid, double* eps_abs,
                     double* eps_rel, double* newton_tol, double* mu, double* tau,
                     double* rho, double* alpha, int* max_newton_iter,
                     int* max_admm_iter, int* max_dc_iter,
                     int* nonconvex, tmpvars** tmp_, problem_data** data_, parameters** para_)
    {
        tmpvars* tmp = (tmpvars*) malloc(sizeof(tmpvars));
        problem_data* data = (problem_data*) malloc(sizeof(problem_data));
        parameters* para = (parameters*) malloc(sizeof(parameters));

        data->p = p[0];
        data->n = n[0];
        data->q = q[0];
        data->X = X;
        data->Y = Y;
        data->lambdaB_grid = lambdaB_grid[0];
        data->lambdaOmega_grid = lambdaOmega_grid[0];
        data->lambdaB = lambdaB;
        data->lambdaOmega = lambdaOmega;
        if (nonconvex[0]==1) data->gamma = gamma;
        data->Z_ini = Z_ini;
        data->sol_path = solution_path;
        if (nonconvex[0]==1) data->sol_path_nc = solution_path_nc;

        para->newton_iter = max_newton_iter[0];
        para->newton_tol = newton_tol[0];
        para->admm_iter = max_admm_iter[0];
        para->alpha = alpha[0];
        para->rho = rho[0];
        para->eps_abs = eps_abs[0];
        para->eps_rel = eps_rel[0];
        para->mu = mu[0];
        para->tau = tau[0];

        if (nonconvex[0]==1)
        {
            tmp->dc_iter = max_dc_iter[0];
        } else {
            tmp->dc_iter = 1;
        }
        tmp->BRowSparsityInd    = (bool*) malloc(sizeof(bool) * data->p);
        tmp->OmegaSparsityInd   = (bool*) malloc(sizeof(bool) * data->q * data->q);
        tmp->pri_res                = .0;
        tmp->dual_res               = .0;
        tmp->obj_cur                = .0;
        tmp->obj_prev           = .0;
        tmp->duality_gap        = .0;
        tmp->newton_step_size   = 1.0;
        tmp->num_newton_steps   = 0;
        tmp->sigmaXX            = (double*) malloc(sizeof(double) * data->p * data->p);
        tmp->sigmaXY            = (double*) malloc(sizeof(double) * data->p * data->q);
        tmp->sigmaYY            = (double*) malloc(sizeof(double) * data->q * data->q);
        tmp->B                  = (double*) malloc(sizeof(double) * data->p * data->q);
        tmp->Omega              = (double*) malloc(sizeof(double) * data->q * data->q);
        tmp->Sigma              = (double*) malloc(sizeof(double) * data->q * data->q);
        tmp->Z                  = (double*) malloc(sizeof(double) * (data->p + data->q) * data->q);
        tmp->Delta              = (double*) malloc(sizeof(double) * (data->p + data->q) * data->q);
        tmp->H                  = (double*) malloc(sizeof(double) * (data->p + data->q) * (data->p + data->q));
        tmp->G                  = (double*) malloc(sizeof(double) * (data->p + data->q) * data->q);
        tmp->Omega_tmp          = (double*) malloc(sizeof(double) * data->q * data->q);
        tmp->B_tmp              = (double*) malloc(sizeof(double) * data->p * data->q);
        tmp->Z_tmp              = (double*) malloc(sizeof(double) * (data->p + data->q) * data->q);
        tmp->Z_prime            = (double*) malloc(sizeof(double) * (data->p + data->q) * data->q);
        tmp->H_eig_vectors      = (double*) malloc(sizeof(double) * (data->p + data->q) * (data->p + data->q));
        tmp->H_eig_values       = (double*) malloc(sizeof(double) * (data->p + data->q));
        tmp->Sigma_eig_vectors  = (double*) malloc(sizeof(double) * data->q * data->q);
        tmp->Sigma_eig_values   = (double*) malloc(sizeof(double) * data->q);
        tmp->Omega_eig_vectors  = (double*) malloc(sizeof(double) * data->q * data->q);
        tmp->Omega_eig_values   = (double*) malloc(sizeof(double) * data->q);
        tmp->Z_cur              = (double*) malloc(sizeof(double) * (data->p + data->q) * data->q);
        tmp->Z_prev             = (double*) malloc(sizeof(double) * (data->p + data->q) * data->q);
        tmp->Gamma_cur          = (double*) malloc(sizeof(double) * (data->p + data->q) * data->q);
        tmp->Gamma_prev         = (double*) malloc(sizeof(double) * (data->p + data->q) * data->q);
        tmp->Gamma_diff         = (double*) malloc(sizeof(double) * (data->p + data->q) * data->q);
        tmp->tmp1               = (double*) malloc(sizeof(double) * (data->p + data->q) * max(data->n, data->q + data->p));
        tmp->tmp2               = (double*) malloc(sizeof(double) * data->p * data->q);
        tmp->tmp3               = (double*) malloc(sizeof(double) * (data->q + data->p) * (data->q + data->p));
        tmp->tmp_values         = (double*) malloc(sizeof(double) * data->q);

        dmat_vset(q[0]*(q[0]+p[0]), .0, tmp->Gamma_prev);
        dmat_vset(q[0]*(q[0]+p[0]), .0, tmp->Gamma_cur);

        // dmat_vcopy(q[0]*p[0], data->Z_ini, tmp->B);
        //dmat_vcopy(q[0]*q[0], data->Z_ini + q[0]*p[0], tmp->Omega);
        dmat_vcopy(q[0]*(q[0]+p[0]), data->Z_ini, tmp->Z);

        // print_dmatrix(X, data->n, data->p);
        // print_dmatrix(Y, data->n, data->q);


        // compuate XtX/n, XtY/n, and YtY/n //
        double n_inv = 1.0 / n[0];
        dmat_C_ATB(n[0], p[0], q[0], &n_inv, X, Y, tmp->sigmaXY);
        dmat_C_ATB(n[0], p[0], p[0], &n_inv, X, X, tmp->sigmaXX);
        dmat_C_ATB(n[0], q[0], q[0], &n_inv, Y, Y, tmp->sigmaYY);


        // dmat_B_ATA(n[0], p[0], X, tmp->sigmaXX);
        // dmat_B_ATA(n[0], q[0], Y, tmp->sigmaYY);
        // dmat_waxpby(p[0]*q[0], 1.0/n[0], tmp->sigmaXY, dzero, tmp->sigmaXY, tmp->sigmaXY);
        // dmat_waxpby(p[0]*p[0], 1.0/n[0], tmp->sigmaXX, dzero, tmp->sigmaXX, tmp->sigmaXX);
        // dmat_waxpby(q[0]*q[0], 1.0/n[0], tmp->sigmaYY, dzero, tmp->sigmaYY, tmp->sigmaYY);

        *data_  = data;
        *para_  = para;
        *tmp_   = tmp;
    }

    void free_tmpvars(tmpvars* tmp, int* nonconvex)
    {
        if (tmp) {
            if (tmp->BRowSparsityInd) free(tmp->BRowSparsityInd);
            if (tmp->OmegaSparsityInd) free(tmp->OmegaSparsityInd);
            if (tmp->sigmaXX) free(tmp->sigmaXX);
            if (tmp->sigmaXY) free(tmp->sigmaXY);
            if (tmp->sigmaYY) free(tmp->sigmaYY);
            if (tmp->B) free(tmp->B);
            if (tmp->Omega) free(tmp->Omega);
            if (tmp->Sigma) free(tmp->Sigma);
            if (tmp->Z) free(tmp->Z);
            if (tmp->Delta) free(tmp->Delta);
            if (tmp->H) free(tmp->H);
            if (tmp->G) free(tmp->G);
            if (tmp->Omega_tmp) free(tmp->Omega_tmp);
            if (tmp->B_tmp) free(tmp->B_tmp);
            if (tmp->Z_tmp) free(tmp->Z_tmp);
            if (tmp->Z_prime) free(tmp->Z_prime);
            if (tmp->H_eig_vectors) free(tmp->H_eig_vectors);
            if (tmp->H_eig_values) free(tmp->H_eig_values);
            if (tmp->Sigma_eig_vectors) free(tmp->Sigma_eig_vectors);
            if (tmp->Sigma_eig_values) free(tmp->Sigma_eig_values);
            if (tmp->Omega_eig_vectors) free(tmp->Omega_eig_vectors);
            if (tmp->Omega_eig_values) free(tmp->Omega_eig_values);
            if (tmp->Z_cur) free(tmp->Z_cur);
            if (tmp->Z_prev) free(tmp->Z_prev);
            if (tmp->Gamma_cur) free(tmp->Gamma_cur);
            if (tmp->Gamma_prev) free(tmp->Gamma_prev);
            if (tmp->Gamma_diff) free(tmp->Gamma_diff);
            if (tmp->tmp1) free(tmp->tmp1);
            if (tmp->tmp2) free(tmp->tmp2);
            if (tmp->tmp3) free(tmp->tmp3);
            if (tmp->tmp_values) free(tmp->tmp_values);
            free(tmp);
        }
    }

    // form H and G, where only upper diagonals of H are stored //
    void form_H_and_G(problem_data* data, tmpvars* tmp)
    {
        int p = data->p, q = data->q, qtq = q*q, ptq = p*q, ppq = p+q, n = data->n;
        int i, j;
        double two_n_inv = 2.0 / n, minus_two_n_inv = - two_n_inv;
        dmat_vcopy(ptq, tmp->Z, tmp->B);
        dmat_vcopy(qtq, tmp->Z + ptq, tmp->Omega);

        eigen_decomp(q, tmp->Omega, tmp->Omega_eig_vectors, tmp->Omega_eig_values);
        dmat_yinvx(q, tmp->Omega_eig_values, tmp->Sigma_eig_values);
        // dmat_B_ATDA(data->q, data->q, tmp->Omega_eig_vectors, tmp->Sigma_eig_values, tmp->Sigma);
        dmat_uB_ATDA(q, q, tmp->Omega_eig_vectors, tmp->Sigma_eig_values, tmp->Sigma);

        // tmp1 = X * B * Sigma //
        dmat_C_ABT(q, q, p, &done, tmp->Sigma, tmp->B, tmp->tmp2);
        dmat_C_ABT(p, n, q, &done, data->X, tmp->tmp2, tmp->tmp1);


        // tmp3 = 2 * tmp1^T * tmp1 / n + Sigma = H_22 //
        dmat_vcopy(qtq, tmp->Sigma, tmp->tmp3);
        dmat_uC_ATA(n, q, &two_n_inv, tmp->tmp1, &done, tmp->tmp3);
        // dmat_C_ATA(n, q, &two_n_inv, tmp->tmp1, &done, tmp->tmp3);

        // tmp2 = H_12 //
        dmat_C_ATB(n, p, q, &minus_two_n_inv, data->X, tmp->tmp1, tmp->tmp2);


        // form H, first store upper diagonal //
        for (i = 0; i < ppq; i++) {
            for (j = i; j < ppq; j++) {
                if (i < p) {
                    if (j < p) tmp->H[i*ppq + j] = 2.0 * tmp->sigmaXX[i*p + j];
                    else tmp->H[i*ppq + j] = tmp->tmp2[i*q + j-p];
                } else {
                    tmp->H[i*ppq + j] = tmp->tmp3[(i-p)*q + j-p];
                }

            }
        }


        // Second, store lower diagonal //
        for (i = 0; i < ppq; i++) {
            for (j = 0; j < i; j++) {
                tmp->H[i*ppq + j] = tmp->H[j*ppq + i];
            }
        }

        // form G //
        // G_11 = tmp2 + 2.0*sigmaXY //
        dmat_waxpby(ptq, 1.0, tmp->tmp2, 2.0, tmp->sigmaXY, tmp->G);

        // tmp3 = .5 * tmp3 + 1.5 * Sigma and G_21 = tmp3 - SigmaYY //
        dmat_waxpby(qtq, .5, tmp->tmp3, 1.5, tmp->Sigma, tmp->tmp3);
        dmat_waxpby(qtq, 1.0, tmp->tmp3, -1.0, tmp->sigmaYY, tmp->G + p*q);

        //    if (!dmat_ispd(p+q, tmp->H))
        //    {
        //        printf("H is not positive definite! \n");
        //        exit(1);
        //    }
    }


    // solve minimize_{Z} .5 * tr(Z * Sigma * Z' * H) - tr(Z' * G) + p(Z); and return Delta = Z_cur - Z //
    void proximal_newton_step(problem_data* data, parameters* para, tmpvars* tmp, double lamB, double lamOmega)
    {
        int p = data->p, q = data->q, ppq = p + q, s = 0;
        int p_plus_q_t_q = (p+q)*q;
        double pri_eps, dual_eps;
        //    printf("hello\n");
        //    print_dmatrix(tmp->Gamma_cur, ppq, q);
        //    print_dmatrix(tmp->Gamma_prev, ppq, q);

        dmat_vset(p_plus_q_t_q, .0, tmp->Delta);
        dmat_vcopy(p_plus_q_t_q, tmp->Z, tmp->Z_cur);
        dmat_vcopy(ppq*ppq, tmp->H, tmp->tmp1);

        // printf("hello \n");
        // print_dmatrix(tmp->tmp1, ppq, ppq);
        eigen_decomp(ppq, tmp->tmp1, tmp->H_eig_vectors, tmp->H_eig_values);


        // printf("hello \n");
        dmat_vcopy(q*q, tmp->Sigma, tmp->tmp3);
        // print_dmatrix(tmp->tmp3, q, q);
        eigen_decomp(q, tmp->tmp3, tmp->Sigma_eig_vectors, tmp->Sigma_eig_values);
        // printf("hello \n");
        tmp->duality_gap = .0;
        int i, j, k; double tmp_real;
        for (i = 0; i < para->admm_iter; i++)
        {
            // update primal variable //
            dmat_vcopy(p_plus_q_t_q, tmp->Z_cur, tmp->Z_prev);
            dmat_waxpby(p_plus_q_t_q, -2.0, tmp->Gamma_cur, para->rho, tmp->Z_prev, tmp->tmp1);
            dmat_waxpby(p_plus_q_t_q, 1.0, tmp->G, 1.0, tmp->tmp1, tmp->tmp1);
            dmat_waxpby(p_plus_q_t_q, 1.0, tmp->Gamma_prev, 1.0, tmp->tmp1, tmp->tmp1);

            dmat_C_ABT(q, q, ppq, &done, tmp->Sigma_eig_vectors, tmp->tmp1, tmp->Z_prime);
            dmat_C_ABT(ppq, ppq, q, &done, tmp->H_eig_vectors, tmp->Z_prime, tmp->tmp1);
            for (j = 0; j < ppq; j++)
                for (k = 0; k < q; k++)
                    tmp->Z_prime[j*q + k] = tmp->tmp1[j*q + k] / (para->rho + tmp->H_eig_values[j] * tmp->Sigma_eig_values[k]);

            // Z_cur = U_h' * Z_prime * U_s //
            dmat_C_AB(q, ppq, q, &done, tmp->Z_prime, tmp->Sigma_eig_vectors, tmp->tmp1);
            dmat_C_ATB(ppq, ppq, q, &done, tmp->H_eig_vectors, tmp->tmp1, tmp->Z_cur);


            // check primal update //
            //        dmat_C_AB(ppq,ppq,q, &done, tmp->H, tmp->Z_cur, tmp->tmp2);
            //        dmat_C_AB(q, ppq, q, &done, tmp->tmp2, tmp->Sigma, tmp->tmp3);


            // udpate dual variable //
            dmat_waxpby(p_plus_q_t_q, 1.0, tmp->Gamma_cur, -1.0, tmp->Gamma_prev, tmp->Gamma_diff);
            dmat_vcopy(p_plus_q_t_q, tmp->Gamma_cur, tmp->Gamma_prev);
            dmat_waxpby(p_plus_q_t_q, 1.0, tmp->Gamma_prev, para->rho, tmp->Z_cur, tmp->Gamma_cur);
            for (j = 0; j < p; j++) {
                if (tmp->BRowSparsityInd[j]) {
                    tmp_real = dmat_norm2(q, tmp->Gamma_cur + j*q);
                    if (tmp_real > lamB) {
                        tmp_real = lamB / tmp_real;
                        dmat_waxpby(q, tmp_real, tmp->Gamma_cur + j*q, .0, tmp->Gamma_cur + j*q, tmp->Gamma_cur + j*q);
                    }
                } else {
                    dmat_vset(q, .0, tmp->Gamma_cur + j*q);
                }
            }
            proj_offdiag_InftyBall(tmp->Gamma_cur + p*q, q, &lamOmega, tmp->OmegaSparsityInd);

            // stopping criterion //
            dmat_waxpby(p_plus_q_t_q, 1.0, tmp->Gamma_cur, -1.0, tmp->Gamma_prev, tmp->tmp1);
            tmp->pri_res = dmat_norm2(p_plus_q_t_q, tmp->tmp1) / para->rho;

            dmat_waxpby(p_plus_q_t_q, done, tmp->Z_cur, dminusone, tmp->Z_prev, tmp->tmp3);
            dmat_waxpby(p_plus_q_t_q, para->rho, tmp->tmp3, 1.0, tmp->Gamma_diff, tmp->tmp3);
            dmat_waxpby(p_plus_q_t_q, done, tmp->tmp3, dminusone, tmp->tmp1, tmp->tmp3);
            tmp->dual_res = dmat_norm2(p_plus_q_t_q, tmp->tmp3);

            // printf("#%d: primal residual is %f and dual residual is %f .\n", i+1, tmp->pri_res, tmp->dual_res);

            dmat_waxpby(p_plus_q_t_q, done, tmp->Z_cur, -1.0/para->rho, tmp->tmp1, tmp->tmp3);
            pri_eps = sqrt(p_plus_q_t_q) * para->eps_abs + para->eps_rel * max(dmat_norm2(p_plus_q_t_q, tmp->tmp3), dmat_norm2(p_plus_q_t_q, tmp->Z_cur));
            dual_eps =sqrt(p_plus_q_t_q) * para->eps_abs + para->eps_rel * dmat_norm2(p_plus_q_t_q, tmp->Gamma_cur);

            if ((tmp->pri_res < pri_eps) && (tmp->dual_res < dual_eps)) {
                // printf("#%d: primal residual is %f and dual residual is %f .\n", i+1, tmp->pri_res, tmp->dual_res);
                break;
            }

            if (i == (s+1)*s*50)
            {
                s+=1;
                if (tmp->pri_res > (para->mu * tmp->dual_res))
                    para->rho = para->rho * para->tau;
                if (tmp->dual_res > (para->mu * tmp->pri_res))
                    para->rho = para->rho / para->tau;
            }

            if (i == (para->admm_iter - 1)) {
                // printf("warning: ADMM reaches maximum number of iteration! \n");
            }

            // print obj. of subproblem //
            //        dmat_C_AB(q, ppq, q, &done, tmp->Z_cur, tmp->Sigma, tmp->tmp1);
            //        dmat_C_AB(ppq, ppq, q, &done, tmp->H, tmp->tmp1, tmp->tmp3);
            //        double primalVal = .5 * dmat_dot(p_plus_q_t_q, tmp->Z_cur, tmp->tmp3) - dmat_dot(p_plus_q_t_q, tmp->Z_cur, tmp->G);
            //        primalVal += penalty(p, q, tmp->Z_cur, tmp->Z_cur+p*q, lamB, lamOmega);
            //        printf("objective value is: %10.20f \n", primalVal);
        }

        // TO DO --- compute final duality gap (current version not correct) //

        //    // tmp3 = H * Z_cur * Sigma //
        //    dmat_C_AB(q, ppq, q, &done, tmp->Z_cur, tmp->Sigma, tmp->tmp1);
        //    dmat_C_AB(ppq, ppq, q, &done, tmp->H, tmp->tmp1, tmp->tmp3);
        //
        //    double primalVal = .5 * dmat_dot(p_plus_q_t_q, tmp->Z_cur, tmp->tmp3) - dmat_dot(p_plus_q_t_q, tmp->Z_cur, tmp->G);
        //    primalVal += penalty(p, q, tmp->Z_cur, tmp->Z_cur+p*q, lamB, lamOmega);
        //
        //    double dualVal = .0;
        //    // tmp1 = U_h * (G - Gamma) * U_s'
        //    dmat_waxpby(p_plus_q_t_q, done, tmp->G, dminusone, tmp->Gamma_cur, tmp->tmp1);
        //    dmat_C_ABT(q, ppq, q, &done, tmp->tmp1, tmp->Sigma_eig_vectors, tmp->tmp3);
        //    dmat_C_AB(ppq, ppq, q, &done, tmp->H_eig_vectors, tmp->tmp3, tmp->tmp1);
        //
        //    for (j = 0; j < ppq; j++)
        //        for (k = 0; k < q; k++)
        //            dualVal += tmp->tmp3[j*q + k] * tmp->tmp3[j*q + k] / (tmp->H_eig_values[j] * tmp->Sigma_eig_values[k]);
        //    dualVal = - .5 * dualVal;
        //    tmp->duality_gap = primalVal - dualVal;
        //    printf("Final duality gap is %f. \n", tmp->duality_gap);


        // Compute gradient and gap //
        //    dmat_C_AB(q, ppq, q, &done, tmp->Z_cur, tmp->Sigma, tmp->tmp1);
        //    dmat_C_AB(ppq, ppq, q, &done, tmp->H, tmp->tmp1, tmp->tmp3);
        //    dmat_waxpby(p_plus_q_t_q, 1.0, tmp->tmp3, -1.0, tmp->G, tmp->tmp3);
        //    dmat_waxpby(p_plus_q_t_q, 1.0, tmp->tmp3, 1.0, tmp->Gamma_cur, tmp->tmp3);
        //    printf("Norm of the final gradient is: %10.20f. \n", dmat_norm2(p_plus_q_t_q, tmp->tmp3)/sqrt(p_plus_q_t_q));
        //    printf("gap is: %10.20f. \n", penalty(p, q, tmp->Z_cur, tmp->Z_cur+p*q, lamB, lamOmega) - dmat_dot(p_plus_q_t_q, tmp->Z_cur, tmp->Gamma_cur));



        // adjust final solution slightly to make small entry exactly zero //

        //    dmat_waxpby(p_plus_q_t_q, done, tmp->Z_cur, dminusone, tmp->Z, tmp->tmp1);
        //    dmat_waxpby(p_plus_q_t_q, done, tmp->Gamma_cur, dminusone, tmp->Gamma_prev, tmp->tmp3);
        //    dmat_waxpby(p_plus_q_t_q, done, tmp->tmp1, 1.0/para->rho, tmp->tmp3, tmp->Delta);

        dmat_waxpby(p_plus_q_t_q, 1.0, tmp->Gamma_prev, para->rho, tmp->Z_cur, tmp->tmp1);
        for (j = 0; j < p; j++) {
            if (tmp->BRowSparsityInd[j] && (dmat_norm2(q, tmp->tmp1 + j*q) < lamB)) {
                dmat_vset(q, .0, tmp->Z_cur + j*q);
            }
        }
        for (i = 0; i < q; i++) {
            for (j = 0; j < q; j++) {
                if ((tmp->OmegaSparsityInd[i*q+j]) && (abs(tmp->tmp1[p*q + i*q + j]) < lamOmega)) tmp->Z_cur[p*q+i*q+j] = .0;
            }
        }
        dmat_waxpby(p_plus_q_t_q, done, tmp->Z_cur, dminusone, tmp->Z, tmp->Delta);
        // make Omega sysmmetric //
        make_sym(tmp->Delta + p*q, q);
    }


    // update tmp->Z to tmp->Z + t * tmp->Delta for some t > 0 //
    void line_search(problem_data* data, parameters* para,tmpvars* tmp, double lamB, double lamOmega)
    {
        int q = data->q, p = data->p;
        int p_plus_q_t_q = (p + q) * q;
        double diff = .0;
        tmp->newton_step_size = 1.0; tmp->num_newton_steps = 0;
        tmp->obj_prev = obj_cGGM(data, para, tmp, tmp->Z, lamB, lamOmega);

        while (1)
        {
            if (tmp->num_newton_steps > 40)
            {
                // printf("Newton step size is %f and termination criterion is %f. \n", diff, tmp->newton_step_size);
                // printf("warning: too many line search attempted! \n");
                // printf("obj diff is: %10.20f. \n", tmp->obj_cur - tmp->obj_prev);
                break;
            }
            // printf("t is %f \n", tmp->newton_step_size);
            // printf("Line search iteration %d \n", tmp->num_newton_steps);
            dmat_waxpby(p_plus_q_t_q, 1.0, tmp->Z, tmp->newton_step_size, tmp->Delta, tmp->Z_tmp);
            dmat_vcopy(q*q, tmp->Z_tmp + p*q, tmp->tmp3);
            if (!dmat_ispd(q, tmp->tmp3)) { tmp->newton_step_size /= 2.0; ++tmp->num_newton_steps; continue; }
            tmp->obj_cur = obj_cGGM(data, para, tmp, tmp->Z_tmp, lamB, lamOmega);

            // obj diff is //
            // printf("obj. diff is: %10.20f. \n", tmp->obj_cur - tmp->obj_prev);

            // remainder term is //
            // printf("remainder term is: %10.20f. \n", tmp->obj_cur - tmp->obj_prev - tmp->newton_step_size * (dmat_dot(q*q, tmp->Sigma, tmp->Delta + p*q) - dmat_dot(p_plus_q_t_q, tmp->G, tmp->Delta)) - penalty(p, q, tmp->Z_tmp, tmp->Z_tmp + p*q, lamB, lamOmega) + penalty(p, q, tmp->Z, tmp->Z + p*q, lamB, lamOmega));

            diff = para->alpha * tmp->newton_step_size * (dmat_dot(q*q, tmp->Sigma, tmp->Delta + p*q) - dmat_dot(p_plus_q_t_q, tmp->G, tmp->Delta));
            diff = diff + para->alpha * (penalty(p, q, tmp->Z_tmp, tmp->Z_tmp + p*q, lamB, lamOmega, tmp->BRowSparsityInd, tmp->OmegaSparsityInd) - penalty(p, q, tmp->Z, tmp->Z + p*q, lamB, lamOmega, tmp->BRowSparsityInd, tmp->OmegaSparsityInd));
            diff = tmp->obj_cur - tmp->obj_prev - diff;
            // printf("diff is: %10.20f. \n", diff);

            if (diff < 0)
            {
                dmat_vcopy(p_plus_q_t_q, tmp->Z_tmp, tmp->Z);
                // printf("line search finish at %d steps! \n", tmp->num_newton_steps);
                // printf("Newton steps: \n");
                // print_dmatrix(tmp->Delta, p+q,q);
                // printf("solution is: \n");
                // print_dmatrix(tmp->Z_tmp, p+q, data->q);
                // printf("objective values is: %f. \n", tmp->obj_cur);
                // printf("prev obj: %10.20f, curr obj: %10.20f. \n", tmp->obj_prev, tmp->obj_cur);
                // printf("obj diff is: %10.20f. \n", tmp->obj_cur - tmp->obj_prev);
                break;
            }
            else { tmp->newton_step_size /= 2.0; ++tmp->num_newton_steps; }
        }
    }

    double penalty(int p, int q, double* B, double* Omega, double lamB, double lamOmega, bool* BInd, bool* OmegaInd)
    {
        double penB = .0, penOmega = .0; int i,j;
        for (i = 0; i < p; i++) if (BInd[i]) penB += dmat_norm2(q, B+i*q);
        penB = lamB * penB;
//         penOmega = dmat_norm1(q*q, Omega);
//         for (i = 0; i < q; i++) penOmega -= Omega[i*q + i];
        for (i = 0; i < q; i++) for (j = 0; j < q; j++) if (OmegaInd[i*q+j]) penOmega += abs(Omega[i*q+j]);
        return (penB + penOmega * lamOmega);
    }

    double obj_cGGM(problem_data* data, parameters* para, tmpvars* tmp, double* sol, double lamB, double lamOmega)
    {
        double obj = .0;
        int q = data->q, p = data->p;
        dmat_vcopy(q*q, sol+p*q, tmp->tmp3);
        eigen_decomp(q, tmp->tmp3, tmp->Omega_eig_vectors, tmp->Omega_eig_values);

        // printf("eigen values are: \n ");
        // print_dmatrix(tmp->Omega_eig_values, 1, q);

        obj = penalty(p, q, sol, sol+p*q, lamB, lamOmega, tmp->BRowSparsityInd, tmp->OmegaSparsityInd);

        // tmp1 = sigmaXX * B //
        dmat_C_ATB(p, p, q, &done, tmp->sigmaXX, sol, tmp->tmp1);
        // tmp3 = B' * tmp1 = B' * sigmaXX * B //
        dmat_C_ATB(p, q, q, &done, sol, tmp->tmp1, tmp->tmp3);

        // tmp1 = Sigma //
        dmat_yinvx(q, tmp->Omega_eig_values, tmp->Sigma_eig_values);
        dmat_uB_ATDA(q, q, tmp->Omega_eig_vectors, tmp->Sigma_eig_values, tmp->tmp1);

        obj = obj + dmat_dot(q*q, tmp->sigmaYY, sol+p*q) - 2.0 * dmat_dot(p*q, sol, tmp->sigmaXY) + dmat_dot(q*q, tmp->tmp1, tmp->tmp3);

        // printf("first term is %f. \n", dmat_dot(q*q, tmp->sigmaYY, sol+p*q));
        // printf("second term is %f. \n", dmat_dot(p*q, sol, tmp->sigmaXY));
        // printf("third term is %f. \n", dmat_dot(q*q, tmp->tmp1, tmp->tmp3));
        // print_dmatrix(tmp->tmp1, q, q);
        // print_dmatrix(tmp->tmp3, q, q);

        int i;
        for (i = 0; i < q; i++) {
            obj -= log(tmp->Omega_eig_values[i]);
        }
        return obj;
    }

#ifdef __cplusplus
}
#endif
