#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include <float.h>

#include "glasso.h"
#include "dmatrix.h"


/*  single glasso training */
void create_tmp_vars(tmpvars** tmp_,int p){
    int p_square = p*p;
    tmpvars* tmp    = (tmpvars*) malloc(sizeof(tmpvars));
    tmp->Delta      = (double *) malloc(sizeof(double) * p_square);
    tmp->Delta_old  = (double *) malloc(sizeof(double) * p_square);
    tmp->tmp1       = (double *) malloc(sizeof(double) * p_square);
    tmp->tmp2       = (double *) malloc(sizeof(double) * p_square);
    tmp->Lambda     = (double *) malloc(sizeof(double) * p);
    tmp->lam_mat    = (double *) malloc(sizeof(double) * p_square);
    *tmp_ = tmp;
}

void free_tmp_vars(tmpvars* tmp){
    if (tmp->Delta != NULL) free(tmp->Delta);
    if (tmp->Delta_old != NULL) free(tmp->Delta_old);
    if (tmp->tmp1 != NULL) free(tmp->tmp1);
    if (tmp->tmp2 != NULL) free(tmp->tmp2);
    if (tmp->Lambda != NULL) free(tmp->Lambda);
    if (tmp->lam_mat  != NULL) free(tmp->lam_mat);
    if (tmp != NULL) free(tmp);
}

void glasso_gen(const double* Sigma_hat, const int* pp,
                const int* num_iter, const double* eps_abs,
                const double* eps_rel,
                double* Rho, double* Alpha, tmpvars* tmp,
                double* Gamma, double* Omega)
{
    int p = pp[0];
    int max_iter = num_iter[0], p_square = p*p;
    int i,j;
    double rho = Rho[0], alpha=Alpha[0],primal_res,dual_res,eps_primal,eps_dual;
    double *Delta = tmp->Delta, *Delta_old = tmp->Delta_old, *tmp1 = tmp->tmp1, *tmp2 = tmp->tmp2, *Lambda = tmp->Lambda, *lam_mat=tmp->lam_mat;
    dmat_waxpby(p_square,1.0,Omega,1.0,Gamma,Delta);
    dmat_waxpby(p_square,1.0/rho,lam_mat,.0,tmp1,tmp1);
    soft_threshold_vec(Delta,p_square,tmp1);
    
    for (j = 0; j < max_iter; j++) {
        // Omega udpate //
        dmat_waxpby(p_square,rho,Delta,-rho,Gamma,tmp1);
        dmat_waxpby(p_square,1.0,tmp1,-1.0,Sigma_hat,tmp1);
        
        eigen_decomp(p,tmp1,tmp2,Lambda);
        for (i = 0; i < p; i++) { Lambda[i] = (Lambda[i] + sqrt(Lambda[i]*Lambda[i]+4*rho))/(2*rho); }
        // if (isnan_mat(Lambda,p)) { printf("Lambda is nan! \n"); }
        // if (isnan_mat(tmp2,p_square)) { printf("tmp2 is nan! \n"); }
        dmat_B_ATDA(p,p,tmp2,Lambda,Omega);
        // if (isnan_mat(Omega,p_square)) { printf("Omega is nan! \n"); }
        
        // Delta update //
        dmat_vcopy(p_square,Delta,Delta_old);
        dmat_waxpby(p_square,alpha,Omega,1.0-alpha,Delta,Delta);
        dmat_waxpby(p_square,1.0,Delta,1.0,Gamma,Delta);
        dmat_waxpby(p_square,1.0/rho,lam_mat,.0,tmp1,tmp1);
        soft_threshold_vec(Delta,p_square,tmp1);
        
        // Gamma update //
        dmat_waxpby(p_square,1.0,Delta_old,-1.0,Delta,tmp1);
        dmat_waxpby(p_square,1.0,Omega,-1.0,Delta,tmp2);
        dmat_waxpby(p_square,1.0,Gamma,1.0-alpha,tmp1,Gamma);
        dmat_waxpby(p_square,1.0,Gamma,alpha,tmp2,Gamma);
        
        // check stopping rule //
        dual_res = rho * dmat_norm2(p_square,tmp1);
        primal_res = dmat_norm2(p_square,tmp2);
        
        eps_primal = eps_abs[0]*p + eps_rel[0]*fmax(dmat_norm2(p_square,Omega),dmat_norm2(p_square,Delta));
        eps_dual = eps_abs[0]*p + eps_rel[0]*dmat_norm2(p_square,Gamma)*rho;
        
        // printf("#%d, primal res: %5.4f, dual res: %5.4f. \n", j, primal_res, dual_res);
        if ((primal_res < eps_primal) && (dual_res < eps_dual)) {
            dmat_vcopy(p_square,Delta,Omega);
            // printf("#%d, primal res: %5.4f, dual res: %5.4f. \n", j, primal_res, dual_res);
            break;
        }
    }
}

void glasso_nonconvex(const double* Sigma_hat, const double* lam,
                      const double* tau, 
                      const int* pp,
                      const int* num_iter, const int* dc_max_iter,
                      const double* eps_abs,
                      const double* eps_rel,
                      double* Rho, double* Alpha, tmpvars* tmp,
                      enum penalty_types* pen_type,
                      double* Gamma_L1, double* Gamma_trL1,
                      double* Omega_L1, double* Omega_trL1)
{
    int i,j,p_square = pp[0]*pp[0];
    dmat_vset(p_square,lam[0],tmp->lam_mat);
    for (i = 0; i < pp[0]; i++) { tmp->lam_mat[i*pp[0]+i] = .0; }
    glasso_gen(Sigma_hat,pp,num_iter,eps_abs,eps_rel,Rho,Alpha,tmp,Gamma_L1,Omega_L1);
    dmat_vcopy(p_square,Omega_L1,Omega_trL1);
    dmat_vcopy(p_square,Gamma_L1,Gamma_trL1);
    for (j = 0; j < dc_max_iter[0]; j++) {
        dmat_vcopy(p_square,tmp->lam_mat,tmp->tmp1);
        if (pen_type[0] == MCP) { for (i = 0; i < p_square; i++) { tmp->lam_mat[i] = lam[0]*fmax(.0, 1.0 - fabs(Omega_trL1[i])/(2*tau[0])); } }
        if (pen_type[0] == Truncated_L1)
        {
            for (i = 0; i < p_square; i++)
            {
                if (fabs(Omega_trL1[i]) < tau[0]) tmp->lam_mat[i] = lam[0];
                else tmp->lam_mat[i] = .0;
            }
        }
        for (i = 0; i < pp[0]; i++) { tmp->lam_mat[i*pp[0]+i] = .0; }
        dmat_waxpby(p_square,1.0,tmp->lam_mat,-1.0,tmp->tmp1,tmp->tmp1);
        // print_dmatrix(tmp->lam_mat,pp[0],pp[0]);
        // print_dmatrix(Omega_trL1,pp[0],pp[0]);
        // printf("DC gap: %5.20f. \n",dmat_norminf(p_square,tmp->tmp1));
        if (dmat_norminf(p_square,tmp->tmp1) < (lam[0]*eps_abs[0])) {
            // printf("DC iteration used: %d. \n", j);
            break;
        }
        glasso_gen(Sigma_hat,pp,num_iter,eps_abs,eps_rel,Rho,Alpha,tmp,Gamma_trL1,Omega_trL1);
    }
}


void glasso_nonconvex_path_new(const double* Sigma_hat, const double* lam,
                               const int* lam_length,const double* tau,
                               const int* tau_length,const int* pp,
                               const int* num_iter, const int* dc_max_iter,
                               const double* eps_abs,
                               const double* eps_rel,
                               double* Rho, double* Alpha,
                               const int* method,
                               double* Omega_L1, double* Omega_trL1)
{
    int i,j, p = pp[0]; tmpvars* tmp;
    enum penalty_types pen_type;
    if (method[0] == 1) { pen_type = Truncated_L1; printf("TLP is used. \n");}
    else if (method[0] == 2) { pen_type = MCP; printf("MCP is used. \n");}
    else { error("nonconvex penlaty not supported! \n"); }
    double* Gamma_L1 = (double *) malloc(sizeof(double) * p*p);
    double* Gamma_trL1 = (double *) malloc(sizeof(double) * p*p);
    create_tmp_vars(&tmp,p);
    // double* tmp = (double *) malloc(sizeof(double) * p * p);
    for (i = 0; i < lam_length[0]; i++){
        // printf("lambda is: %f.\n", lam[i]);
        // get solution for maximum tau //
        if (i == 0) {
            dmat_vset(p*p,.0,Gamma_L1);
            glasso_nonconvex(Sigma_hat,lam,tau,pp,num_iter,dc_max_iter,eps_abs,eps_rel,Rho,Alpha,tmp,&pen_type,Gamma_L1,Gamma_trL1,Omega_L1,Omega_trL1);
        } else {
            dmat_vcopy(p*p,Omega_L1+(i-1)*p*p,Omega_L1+i*p*p);
            glasso_nonconvex(Sigma_hat,lam+i,tau,pp,num_iter,dc_max_iter,eps_abs,eps_rel,Rho,Alpha,tmp,&pen_type,Gamma_L1,Gamma_trL1,Omega_L1+i*p*p,Omega_trL1+i*p*p);
        }
        for (j = 1; j < tau_length[0]; j++){
            // Rho[0] = lam[i]; // adaptively varying rho //
            glasso_nonconvex(Sigma_hat,lam+i,tau+j,pp,num_iter,dc_max_iter,eps_abs,eps_rel,Rho,Alpha,tmp,&pen_type,Gamma_L1,Gamma_trL1,Omega_L1+i*p*p,Omega_trL1+i*tau_length[0]*p*p+j*p*p);
        }
    }
    if (Gamma_L1 != NULL) free(Gamma_L1);
    if (Gamma_trL1 != NULL) free(Gamma_trL1);
    free_tmp_vars(tmp);
}






// x is a n*p dim vector //
// genenral group lasso proximal operator for p vectors (of dimension n) //
void group_lasso_gen(int n,int p,double* x,const double* rho,
                     const double* lam_mat_L1, const double* lam_mat_group)
{
    soft_threshold_vec(x,n*p,lam_mat_L1);
    int i,j;
    double tmp[p], tmp1, tmp2;
    for (i = 0; i < n; i++) {
        if (fabs(lam_mat_group[i])>1e-14) {
            for (j = 0; j < p; j++) { tmp[j] = x[j*n+i]; }
            tmp1 = dmat_norm2(p,tmp);
            tmp2 = lam_mat_group[i];
            for (j = 0; j < p; j++) { x[j*n+i] = fmax(1.0 - tmp2/tmp1, .0) * tmp[j]; }
            // if (tmp1 < tmp2) { for (j = 0; j < p; j++) { x[j*n+i] = .0; } }
            // else { for (j = 0; j < p; j++) { x[j*n+i] = (1.0 - tmp2/tmp1) * tmp[j]; } }
        }
    }
}



void create_tmp_vars_group(tmpvars_group** tmp_,int p, int K){
    int p_square = p*p;
    tmpvars_group* tmp    = (tmpvars_group*) malloc(sizeof(tmpvars_group));
    tmp->Delta      = (double *) malloc(sizeof(double) * p_square*K);
    tmp->Delta_old  = (double *) malloc(sizeof(double) * p_square*K);
    tmp->tmp1       = (double *) malloc(sizeof(double) * p_square*K);
    tmp->tmp2       = (double *) malloc(sizeof(double) * p_square*K);
    tmp->Lambda     = (double *) malloc(sizeof(double) * p);
    // tmp->thred_L1   = (bool   *) malloc(sizeof(bool)   * p_square*K);
    // tmp->thred_group= (bool   *) malloc(sizeof(bool)   * p_square);
    tmp->lam_mat_L1 = (double *) malloc(sizeof(double) * p_square*K);
    tmp->lam_mat_group = (double *) malloc(sizeof(double) * p_square);
    *tmp_ = tmp;
}

void free_tmp_vars_group(tmpvars_group* tmp){
    if (tmp->Delta != NULL) free(tmp->Delta);
    if (tmp->Delta_old != NULL) free(tmp->Delta_old);
    if (tmp->tmp1 != NULL) free(tmp->tmp1);
    if (tmp->tmp2 != NULL) free(tmp->tmp2);
    if (tmp->Lambda != NULL) free(tmp->Lambda);
    // if (tmp->thred_L1  != NULL) free(tmp->thred_L1);
    // if (tmp->thred_group  != NULL) free(tmp->thred_group);
    if (tmp->lam_mat_L1 != NULL) free(tmp->lam_mat_L1);
    if (tmp->lam_mat_group != NULL) free(tmp->lam_mat_group);
    if (tmp != NULL) free(tmp);
}


void group_gen(const double* Sigma_hat,
               const int* pp, const int* KK,
               const int* num_iter, const double* eps_abs,
               const double* eps_rel,
               double* Rho, double* Alpha, tmpvars_group* tmp,
               double* Gamma, double* Omega)
{
    int p = pp[0]; int K = KK[0];
    int max_iter = num_iter[0], p_square = p*p;
    int i,j,k;
    double rho = Rho[0], alpha=Alpha[0],primal_res,dual_res,eps_primal,eps_dual;
    double *Delta = tmp->Delta, *Delta_old = tmp->Delta_old, *tmp1 = tmp->tmp1,
            *tmp2 = tmp->tmp2, *Lambda = tmp->Lambda, *lam_mat_L1 = tmp->lam_mat_L1,
            *lam_mat_group = tmp->lam_mat_group;
    dmat_waxpby(K*p_square,1.0,Omega,1.0,Gamma,Delta);
    dmat_waxpby(K*p_square,1.0/rho,lam_mat_L1,.0,tmp1,tmp1);
    dmat_waxpby(p_square,1.0/rho,lam_mat_group,.0,tmp2,tmp2);
    group_lasso_gen(p_square,K,Delta,Rho,tmp1,tmp2);
    for (j = 0; j < max_iter; j++) {
        // Omega_k udpate //
        for (k = 0; k < K; k++) {
            dmat_waxpby(p_square,rho,Delta+k*p_square,-rho,Gamma+k*p_square,tmp1);
            dmat_waxpby(p_square,1.0,tmp1,-1.0,Sigma_hat+k*p_square,tmp1);
            eigen_decomp(p,tmp1,tmp2,Lambda);
            
            for (i = 0; i < p; i++){ Lambda[i] = (Lambda[i] + sqrt(Lambda[i]*Lambda[i]+4*rho))/(2*rho); }
            // dmat_B_ADAT(p,p,tmp2,Lambda,Omega+k*p_square);
            dmat_B_ATDA(p,p,tmp2,Lambda,Omega+k*p_square);
        }
        
        // Delta update //
        dmat_vcopy(p_square*K,Delta,Delta_old);
        dmat_waxpby(p_square*K,alpha,Omega,1.0-alpha,Delta,Delta);
        dmat_waxpby(p_square*K,1.0,Delta,1.0,Gamma,Delta);
        dmat_waxpby(K*p_square,1.0/rho,lam_mat_L1,.0,tmp1,tmp1);
        dmat_waxpby(p_square,1.0/rho,lam_mat_group,.0,tmp2,tmp2);
        group_lasso_gen(p_square,K,Delta,Rho,tmp1,tmp2);
        
        // Gamma update //
        dmat_waxpby(p_square*K,1.0,Delta_old,-1.0,Delta,tmp1);
        dmat_waxpby(p_square*K,1.0,Omega,-1.0,Delta,tmp2);
        dmat_waxpby(p_square*K,1.0,Gamma,1.0-alpha,tmp1,Gamma);
        dmat_waxpby(p_square*K,1.0,Gamma,alpha,tmp2,Gamma);
        
        // check stopping rule //
        dual_res = rho * dmat_norm2(p_square*K,tmp1);
        primal_res = dmat_norm2(p_square*K,tmp2);
        
        eps_primal = eps_abs[0]*p*sqrt(K) + eps_rel[0]*fmax(dmat_norm2(p_square*K,Omega),dmat_norm2(p_square*K,Delta));
        eps_dual = eps_abs[0]*p*sqrt(K) + eps_rel[0]*dmat_norm2(p_square*K,Gamma)*rho;
        
        
        // printf("#%d, primal res: %5.4f, dual res: %5.4f. \n", j, primal_res, dual_res);
        if ((primal_res < eps_primal) && (dual_res < eps_dual)) {
            dmat_vcopy(p_square*K,Delta,Omega);
            // printf("#%d, primal res: %5.4f, dual res: %5.4f. \n", j, primal_res, dual_res);
            // printf("ADMM iterations is %d \n", j);
            break;
        }
    }
}

void group_nonconvex(const double* Sigma_hat, const double* lam,
                    const double* nu, const double* tau,
                    const int* pp, const int* KK,
                    const int* num_iter, const int* dc_max_iter,
                    const double* eps_abs,
                    const double* eps_rel,
                    double* Rho, double* Alpha, tmpvars_group* tmp,
                     enum penalty_types* pen_type,
                    double* Gamma_L1, double* Gamma_trL1,
                    double* Omega_L1, double* Omega_trL1)
{
    int i,j,k,K = KK[0], p_square = pp[0]*pp[0];
    dmat_vset(p_square*K,lam[0],tmp->lam_mat_L1);
    dmat_vset(p_square,nu[0],tmp->lam_mat_group);
    // double rho = Rho[0], alpha = Alpha[0];
    // if (fmax(lam[0],nu[0]) < (1e-5*pp[0])) { rho = .01; alpha = 1.0; }
    double* tmp3 = tmp->tmp2 + p_square;
    for (i = 0; i < pp[0]; i++) {
        tmp->lam_mat_group[i*pp[0]+i] = .0;
        for (k = 0; k < K; k++){
            tmp->lam_mat_L1[k*p_square+i*pp[0]+i] = .0;
        }
    }
    group_gen(Sigma_hat,pp,KK,num_iter,eps_abs,eps_rel,Rho,Alpha,tmp,Gamma_L1,Omega_L1);
    dmat_vcopy(p_square*K,Omega_L1,Omega_trL1);
    dmat_vcopy(p_square*K,Gamma_L1,Gamma_trL1);
    
    // print_bmatrix(thred,pp[0],pp[0]);
    for (j = 0; j < dc_max_iter[0]; j++) {
        dmat_vcopy(p_square*K,tmp->lam_mat_L1,tmp->tmp1);
        dmat_vcopy(p_square,tmp->lam_mat_group,tmp->tmp2);
        for (i = 0; i < p_square; i++) {
            tmp3[i] = .0;
            for (k = 0; k < K; k++){ tmp3[i] += (Omega_trL1[k*p_square+i] * Omega_trL1[k*p_square+i]); }
            tmp3[i] = sqrt(tmp3[i]);
        }
        if (pen_type[0] == MCP)
        {
            for (i = 0; i < K*p_square; i++) { tmp->lam_mat_L1[i] = lam[0]*fmax(.0, 1.0 - fabs(Omega_trL1[i])/(2*tau[0])); }
            for (i = 0; i < p_square; i++) { tmp->lam_mat_group[i] = lam[0]*fmax(.0,1.0 - tmp3[i]/(2*tau[0])); }
        }
        if (pen_type[0] == Truncated_L1)
        {
            for (i = 0; i < K*p_square; i++)
            {
                if (fabs(Omega_trL1[i]) < tau[0]) { tmp->lam_mat_L1[i] = lam[0]; }
                else { tmp->lam_mat_L1[i] = .0; }
            }
            for (i = 0; i < p_square; i++){
                if (tmp3[i] < tau[0]) tmp->lam_mat_group[i] = lam[0];
                else tmp->lam_mat_group[i] = .0;
            }
        }
        for (i = 0; i < pp[0]; i++) {
            tmp->lam_mat_group[i*pp[0]+i] = .0;
            for (k = 0; k < K; k++){
                tmp->lam_mat_L1[k*p_square+i*pp[0]+i] = .0;
            }
        }
        dmat_waxpby(p_square*K,1.0,tmp->lam_mat_L1,-1.0,tmp->tmp1,tmp->tmp1);
        dmat_waxpby(p_square,1.0,tmp->lam_mat_group,-1.0,tmp->tmp2,tmp->tmp2);
        // print_bmatrix(thred,pp[0],pp[0]);
        if (fmax(dmat_norminf(p_square*K,tmp->tmp1),dmat_norminf(p_square,tmp->tmp2)) < (lam[0]*1e-10)) {
            // printf("DC iteration used: %d. \n", j);
            break;
        }
        group_gen(Sigma_hat,pp,KK,num_iter,eps_abs,eps_rel,Rho,Alpha,tmp,Gamma_trL1,Omega_trL1);
    }
}

void group_nonconvex_path(const double* Sigma_hat, const double* lam,
                         const int* lam_length,const double* nu,
                         const int* nu_length, const double* tau,
                         const int* pp,const int* KK,
                         const int* num_iter, const int* dc_max_iter,
                         const double* eps_abs,
                         const double* eps_rel,
                         double* Rho, double* Alpha,
                           const int* method,
                         double* Omega_L1, double* Omega_trL1)
{
    int i, j, K = KK[0], p = pp[0]; tmpvars_group* tmp;
    enum penalty_types pen_type;
    if (method[0] == 1) { pen_type = Truncated_L1; }
    else if (method[0] == 2) { pen_type = MCP; }
    else { error("nonconvex penlaty not supported! \n"); }
    double* Gamma_L1 = (double *) malloc(sizeof(double) * p*p*K);
    double* Gamma_trL1 = (double *) malloc(sizeof(double) * p*p*K);
    create_tmp_vars_group(&tmp,p,K);
    
    // double* tmp = (double *) malloc(sizeof(double) * p * p);
    for (i = 0; i < lam_length[0]; i++)
    {
        for (j = 0; j < nu_length[0]; j++){
            // Rho[0] = lam[i]; // adaptively varying rho //
            // printf("lambda is: %f, nu is: %f\n", lam[i],nu[j]);
            if (i == 0 && j == 0) {
                dmat_vset(p*p*K,.0,Gamma_L1);
                group_nonconvex(Sigma_hat,lam,nu,tau,pp,KK,num_iter,dc_max_iter,eps_abs,eps_rel,Rho,Alpha,
                               tmp,&pen_type,Gamma_L1,Gamma_trL1,Omega_L1,Omega_trL1);
            } else if (j == 0) {
                dmat_vcopy(p*p*K,Omega_L1+(i-1)*nu_length[0]*p*p*K,Omega_L1+i*nu_length[0]*p*p*K);
                group_nonconvex(Sigma_hat,lam+i,nu,tau,pp,KK,num_iter,dc_max_iter,eps_abs,eps_rel,Rho,Alpha,
                               tmp,&pen_type,Gamma_L1,Gamma_trL1,Omega_L1+i*nu_length[0]*p*p*K,Omega_trL1+i*nu_length[0]*p*p*K);
            } else
            {
                dmat_vcopy(p*p*K,Omega_L1+i*nu_length[0]*K*p*p+(j-1)*K*p*p,Omega_L1+i*nu_length[0]*K*p*p+j*K*p*p);
                group_nonconvex(Sigma_hat,lam+i,nu+j,tau,pp,KK,num_iter,dc_max_iter,eps_abs,eps_rel,Rho,Alpha,
                               tmp,&pen_type,Gamma_L1,Gamma_trL1,Omega_L1+i*nu_length[0]*K*p*p+j*K*p*p,
                           Omega_trL1+i*nu_length[0]*K*p*p+j*K*p*p);
            }
        }
    }
    if (Gamma_L1 != NULL) free(Gamma_L1);
    if (Gamma_trL1 != NULL) free(Gamma_trL1);
    free_tmp_vars_group(tmp);
}



// lam is a K dimensional vector //
void group_nonconvex_new(const double* Sigma_hat, const double* lam,
                     const double* nu, const double* tau,
                     const int* pp, const int* KK,
                     const int* num_iter, const int* dc_max_iter,
                     const double* eps_abs,
                     const double* eps_rel,
                     double* Rho, double* Alpha, tmpvars_group* tmp,
                     enum penalty_types* pen_type,
                     double* Gamma_L1, double* Gamma_trL1,
                     double* Omega_L1, double* Omega_trL1)
{
    int i,j,k,K = KK[0], p_square = pp[0]*pp[0];
    for (k = 1; k < K; k++){ dmat_vset(p_square,lam[k],tmp->lam_mat_L1+k*p_square); }
    dmat_vset(p_square,nu[0],tmp->lam_mat_group);
    // double rho = Rho[0], alpha = Alpha[0];
    // if (fmax(lam[0],nu[0]) < (1e-5*pp[0])) { rho = .01; alpha = 1.0; }
    double* tmp3 = tmp->tmp2 + p_square;
    for (i = 0; i < pp[0]; i++) {
        tmp->lam_mat_group[i*pp[0]+i] = .0;
        for (k = 0; k < K; k++){
            tmp->lam_mat_L1[k*p_square+i*pp[0]+i] = .0;
        }
    }
    group_gen(Sigma_hat,pp,KK,num_iter,eps_abs,eps_rel,Rho,Alpha,tmp,Gamma_L1,Omega_L1);
    dmat_vcopy(p_square*K,Omega_L1,Omega_trL1);
    dmat_vcopy(p_square*K,Gamma_L1,Gamma_trL1);
    
    // print_bmatrix(thred,pp[0],pp[0]);
    for (j = 0; j < dc_max_iter[0]; j++) {
        dmat_vcopy(p_square*K,tmp->lam_mat_L1,tmp->tmp1);
        dmat_vcopy(p_square,tmp->lam_mat_group,tmp->tmp2);
        for (i = 0; i < p_square; i++) {
            tmp3[i] = .0;
            for (k = 0; k < K; k++){ tmp3[i] += (Omega_trL1[k*p_square+i] * Omega_trL1[k*p_square+i]); }
            tmp3[i] = sqrt(tmp3[i]);
        }
        if (pen_type[0] == MCP)
        {
            for (k = 0; k < K; k++) {
                for (i = 0; i < p_square; i++) { tmp->lam_mat_L1[i+k*p_square] = lam[k]*fmax(.0, 1.0 - fabs(Omega_trL1[i+k*p_square])/(2*tau[0])); }
            }
            for (i = 0; i < p_square; i++) { tmp->lam_mat_group[i] = lam[0]*fmax(.0,1.0 - tmp3[i]/(2*tau[0])); }
        }
        if (pen_type[0] == Truncated_L1)
        {
            for (k = 0; k < K; k++) {
                for (i = 0; i < p_square; i++)
                {
                    if (fabs(Omega_trL1[i+k*p_square]) < tau[0]) { tmp->lam_mat_L1[i+k*p_square] = lam[k]; }
                    else { tmp->lam_mat_L1[i+k*p_square] = .0; }
                }
            }
            for (i = 0; i < p_square; i++){
                if (tmp3[i] < tau[0]) tmp->lam_mat_group[i] = lam[0];
                else tmp->lam_mat_group[i] = .0;
            }
        }
        for (i = 0; i < pp[0]; i++) {
            tmp->lam_mat_group[i*pp[0]+i] = .0;
            for (k = 0; k < K; k++){
                tmp->lam_mat_L1[k*p_square+i*pp[0]+i] = .0;
            }
        }
        dmat_waxpby(p_square*K,1.0,tmp->lam_mat_L1,-1.0,tmp->tmp1,tmp->tmp1);
        dmat_waxpby(p_square,1.0,tmp->lam_mat_group,-1.0,tmp->tmp2,tmp->tmp2);
        // print_bmatrix(thred,pp[0],pp[0]);
        if (fmax(dmat_norminf(p_square*K,tmp->tmp1),dmat_norminf(p_square,tmp->tmp2)) < (lam[0]*1e-10)) {
            // printf("DC iteration used: %d. \n", j);
            break;
        }
        group_gen(Sigma_hat,pp,KK,num_iter,eps_abs,eps_rel,Rho,Alpha,tmp,Gamma_trL1,Omega_trL1);
    }
}













// Group lasso version: adding maximum eigenvalue constraint R = \sqrt{2.0 * tau / lam}: MCP only // 
void group_gen_constrain(const double* Sigma_hat, const double* R, 
               const int* pp, const int* KK,
               const int* num_iter, const double* eps_abs,
               const double* eps_rel,
               double* Rho, double* Alpha, tmpvars_group* tmp,
               double* Gamma, double* Omega)
{
    int p = pp[0]; int K = KK[0];
    int max_iter = num_iter[0], p_square = p*p;
    int i,j,k;
    double rho = Rho[0], alpha=Alpha[0],primal_res,dual_res,eps_primal,eps_dual;
    double *Delta = tmp->Delta, *Delta_old = tmp->Delta_old, *tmp1 = tmp->tmp1,
            *tmp2 = tmp->tmp2, *Lambda = tmp->Lambda, *lam_mat_L1 = tmp->lam_mat_L1,
            *lam_mat_group = tmp->lam_mat_group;
    dmat_waxpby(K*p_square,1.0,Omega,1.0,Gamma,Delta);
    dmat_waxpby(K*p_square,1.0/rho,lam_mat_L1,.0,tmp1,tmp1);
    dmat_waxpby(p_square,1.0/rho,lam_mat_group,.0,tmp2,tmp2);
    group_lasso_gen(p_square,K,Delta,Rho,tmp1,tmp2);
    for (j = 0; j < max_iter; j++) {
        // Omega_k udpate //
        for (k = 0; k < K; k++) {
            dmat_waxpby(p_square,rho,Delta+k*p_square,-rho,Gamma+k*p_square,tmp1);
            dmat_waxpby(p_square,1.0,tmp1,-1.0,Sigma_hat+k*p_square,tmp1);
            eigen_decomp(p,tmp1,tmp2,Lambda);
            
            for (i = 0; i < p; i++){ Lambda[i] = fmin(R[0], (Lambda[i] + sqrt(Lambda[i]*Lambda[i]+4*rho))/(2*rho)); }
            // dmat_B_ADAT(p,p,tmp2,Lambda,Omega+k*p_square);
            dmat_B_ATDA(p,p,tmp2,Lambda,Omega+k*p_square);
        }
        
        // Delta update //
        dmat_vcopy(p_square*K,Delta,Delta_old);
        dmat_waxpby(p_square*K,alpha,Omega,1.0-alpha,Delta,Delta);
        dmat_waxpby(p_square*K,1.0,Delta,1.0,Gamma,Delta);
        dmat_waxpby(K*p_square,1.0/rho,lam_mat_L1,.0,tmp1,tmp1);
        dmat_waxpby(p_square,1.0/rho,lam_mat_group,.0,tmp2,tmp2);
        group_lasso_gen(p_square,K,Delta,Rho,tmp1,tmp2);
        
        // Gamma update //
        dmat_waxpby(p_square*K,1.0,Delta_old,-1.0,Delta,tmp1);
        dmat_waxpby(p_square*K,1.0,Omega,-1.0,Delta,tmp2);
        dmat_waxpby(p_square*K,1.0,Gamma,1.0-alpha,tmp1,Gamma);
        dmat_waxpby(p_square*K,1.0,Gamma,alpha,tmp2,Gamma);
        
        // check stopping rule //
        dual_res = rho * dmat_norm2(p_square*K,tmp1);
        primal_res = dmat_norm2(p_square*K,tmp2);
        
        eps_primal = eps_abs[0]*p*sqrt(K) + eps_rel[0]*fmax(dmat_norm2(p_square*K,Omega),dmat_norm2(p_square*K,Delta));
        eps_dual = eps_abs[0]*p*sqrt(K) + eps_rel[0]*dmat_norm2(p_square*K,Gamma)*rho;
        
        
        // printf("#%d, primal res: %5.4f, dual res: %5.4f. \n", j, primal_res, dual_res);
        if ((primal_res < eps_primal) && (dual_res < eps_dual)) {
            dmat_vcopy(p_square*K,Delta,Omega);
            // printf("#%d, primal res: %5.4f, dual res: %5.4f. \n", j, primal_res, dual_res);
            // printf("ADMM iterations is %d \n", j);
            break;
        }
    }
}


void group_nonconvex_constrain(const double* Sigma_hat, const double* lam,
                    const double* nu, const double* tau, const double* R, 
                    const int* pp, const int* KK,
                    const int* num_iter, const int* dc_max_iter,
                    const double* eps_abs,
                    const double* eps_rel,
                    double* Rho, double* Alpha, tmpvars_group* tmp,
                     enum penalty_types* pen_type,
                    double* Gamma_L1, double* Gamma_trL1,
                    double* Omega_L1, double* Omega_trL1)
{
    int i,j,k,K = KK[0], p_square = pp[0]*pp[0];
    dmat_vset(p_square*K,lam[0],tmp->lam_mat_L1);
    dmat_vset(p_square,nu[0],tmp->lam_mat_group);
    // double rho = Rho[0], alpha = Alpha[0];
    // if (fmax(lam[0],nu[0]) < (1e-5*pp[0])) { rho = .01; alpha = 1.0; }
    double* tmp3 = tmp->tmp2 + p_square;
    for (i = 0; i < pp[0]; i++) {
        tmp->lam_mat_group[i*pp[0]+i] = .0;
        for (k = 0; k < K; k++){
            tmp->lam_mat_L1[k*p_square+i*pp[0]+i] = .0;
        }
    }
    group_gen_constrain(Sigma_hat,R,pp,KK,num_iter,eps_abs,eps_rel,Rho,Alpha,tmp,Gamma_L1,Omega_L1);
    dmat_vcopy(p_square*K,Omega_L1,Omega_trL1);
    dmat_vcopy(p_square*K,Gamma_L1,Gamma_trL1);
    
    // print_bmatrix(thred,pp[0],pp[0]);
    for (j = 0; j < dc_max_iter[0]; j++) {
        dmat_vcopy(p_square*K,tmp->lam_mat_L1,tmp->tmp1);
        dmat_vcopy(p_square,tmp->lam_mat_group,tmp->tmp2);
        for (i = 0; i < p_square; i++) {
            tmp3[i] = .0;
            for (k = 0; k < K; k++){ tmp3[i] += (Omega_trL1[k*p_square+i] * Omega_trL1[k*p_square+i]); }
            tmp3[i] = sqrt(tmp3[i]);
        }
        if (pen_type[0] == MCP)
        {
            for (i = 0; i < K*p_square; i++) { tmp->lam_mat_L1[i] = lam[0]*fmax(.0, 1.0 - fabs(Omega_trL1[i])/(2*tau[0])); }
            for (i = 0; i < p_square; i++) { tmp->lam_mat_group[i] = lam[0]*fmax(.0,1.0 - tmp3[i]/(2*tau[0])); }
        }
        if (pen_type[0] == Truncated_L1)
        {
            for (i = 0; i < K*p_square; i++)
            {
                if (fabs(Omega_trL1[i]) < tau[0]) { tmp->lam_mat_L1[i] = lam[0]; }
                else { tmp->lam_mat_L1[i] = .0; }
            }
            for (i = 0; i < p_square; i++){
                if (tmp3[i] < tau[0]) tmp->lam_mat_group[i] = lam[0];
                else tmp->lam_mat_group[i] = .0;
            }
        }
        for (i = 0; i < pp[0]; i++) {
            tmp->lam_mat_group[i*pp[0]+i] = .0;
            for (k = 0; k < K; k++){
                tmp->lam_mat_L1[k*p_square+i*pp[0]+i] = .0;
            }
        }
        dmat_waxpby(p_square*K,1.0,tmp->lam_mat_L1,-1.0,tmp->tmp1,tmp->tmp1);
        dmat_waxpby(p_square,1.0,tmp->lam_mat_group,-1.0,tmp->tmp2,tmp->tmp2);
        // print_bmatrix(thred,pp[0],pp[0]);
        if (fmax(dmat_norminf(p_square*K,tmp->tmp1),dmat_norminf(p_square,tmp->tmp2)) < (lam[0]*1e-10)) {
            // printf("DC iteration used: %d. \n", j);
            break;
        }
        group_gen_constrain(Sigma_hat,R,pp,KK,num_iter,eps_abs,eps_rel,Rho,Alpha,tmp,Gamma_trL1,Omega_trL1);
    }
}

void group_nonconvex_constrain_path(const double* Sigma_hat, const double* lam,
                         const int* lam_length,const double* nu,
                         const int* nu_length, const double* R_, 
                         const double* tau, const int* pp,const int* KK,
                         const int* num_iter, const int* dc_max_iter,
                         const double* eps_abs,
                         const double* eps_rel,
                         double* Rho, double* Alpha,
                           const int* method,
                         double* Omega_L1, double* Omega_trL1)
{
    int i, j, K = KK[0], p = pp[0]; tmpvars_group* tmp;
    enum penalty_types pen_type;
    if (method[0] == 1) { pen_type = Truncated_L1; }
    else if (method[0] == 2) { pen_type = MCP; }
    else { error("nonconvex penlaty not supported! \n"); }
    double* Gamma_L1 = (double *) malloc(sizeof(double) * p*p*K);
    double* Gamma_trL1 = (double *) malloc(sizeof(double) * p*p*K);
    create_tmp_vars_group(&tmp,p,K);
    double R = R_[0];
    printf("R is: %f \n", R);
    // double* tmp = (double *) malloc(sizeof(double) * p * p);
    for (i = 0; i < lam_length[0]; i++)
    {
        // R = sqrt(200.0 * tau[0] / lam[i]); 
        for (j = 0; j < nu_length[0]; j++){
            // Rho[0] = lam[i]; // adaptively varying rho //
            // printf("lambda is: %f, nu is: %f\n", lam[i],nu[j]);
            if (i == 0 && j == 0) {
                dmat_vset(p*p*K,.0,Gamma_L1);
                group_nonconvex_constrain(Sigma_hat,lam,nu,tau,&R,pp,KK,num_iter,dc_max_iter,eps_abs,eps_rel,Rho,Alpha,
                               tmp,&pen_type,Gamma_L1,Gamma_trL1,Omega_L1,Omega_trL1);
            } else if (j == 0) {
                dmat_vcopy(p*p*K,Omega_L1+(i-1)*nu_length[0]*p*p*K,Omega_L1+i*nu_length[0]*p*p*K);
                group_nonconvex_constrain(Sigma_hat,lam+i,nu,tau,&R,pp,KK,num_iter,dc_max_iter,eps_abs,eps_rel,Rho,Alpha,
                               tmp,&pen_type,Gamma_L1,Gamma_trL1,Omega_L1+i*nu_length[0]*p*p*K,Omega_trL1+i*nu_length[0]*p*p*K);
            } else
            {
                dmat_vcopy(p*p*K,Omega_L1+i*nu_length[0]*K*p*p+(j-1)*K*p*p,Omega_L1+i*nu_length[0]*K*p*p+j*K*p*p);
                group_nonconvex_constrain(Sigma_hat,lam+i,nu+j,tau,&R,pp,KK,num_iter,dc_max_iter,eps_abs,eps_rel,Rho,Alpha,
                               tmp,&pen_type,Gamma_L1,Gamma_trL1,Omega_L1+i*nu_length[0]*K*p*p+j*K*p*p,
                           Omega_trL1+i*nu_length[0]*K*p*p+j*K*p*p);
            }
        }
    }
    if (Gamma_L1 != NULL) free(Gamma_L1);
    if (Gamma_trL1 != NULL) free(Gamma_trL1);
    free_tmp_vars_group(tmp);
}







