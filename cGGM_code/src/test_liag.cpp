#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#include "dmatrix.h"
#include "def.h"

#define N 5
#define NRHS 3
#define LDA N
#define LDB NRHS


void test_C_ATB()
{
    double a[N*LDA] = {
        6.80, -6.05, -0.45,  8.32, -9.67,
        -2.11, -3.30,  2.58,  2.71, -5.14,
        5.66, 5.36, -2.70,  4.35, -7.26,
        5.97, -4.44,  0.27, -7.17, 6.08,
        8.23, 1.08,  9.04,  2.14, -6.87
    };
    double b[N*LDB] = {
        4.02, -1.56, 9.81,
        6.19,  4.00, -4.09,
        -8.22, -8.67, -4.57,
        -7.57, 1.75, -8.61,
        -3.03,  2.86, 8.99
    };
    
    double c[LDA*LDB];
    
    printf("Entry matrix A: \n");
    print_dmatrix(a, N, LDA);
    printf("Entry matrix B: \n");
    print_dmatrix(b, N, LDB);
    
    dmat_C_ATB(N, LDA, LDB, &done, a, b, c);
    
    printf("AT * B is: \n");
    print_dmatrix(c, LDA, LDB);
    
    dmat_C_ATB_new(N, LDA, LDB, &done, a, b, c);
    
    printf("AT * B is: \n");
    print_dmatrix(c, LDA, LDB);
    
}



void test_C_ABT()
{
    double a[LDA*N] = {
        6.80, -6.05, -0.45,  8.32, -9.67,
        -2.11, -3.30,  2.58,  2.71, -5.14,
        5.66, 5.36, -2.70,  4.35, -7.26,
        5.97, -4.44,  0.27, -7.17, 6.08,
        8.23, 1.08,  9.04,  2.14, -6.87
    };
    double b[LDB*N] = {
        4.02, -1.56, 9.81, 6.19,  4.00,
        -4.09, -8.22, -8.67, -4.57, -7.57,
        1.75, -8.61, -3.03,  2.86, 8.99
    };
    
    double c[LDA*LDB];
    
    printf("Entry matrix A: \n");
    print_dmatrix(a, LDA, N);
    printf("Entry matrix B: \n");
    print_dmatrix(b, LDB, N);
    
    dmat_C_ABT(N, LDA, LDB, &done, a, b, c);
    
    printf("A * BT is: \n");
    print_dmatrix(c, LDA, LDB);
}


void test_eigen()
{
    double a[NRHS*NRHS] = {
        6.80, -1.05, -0.45,
        -1.05, 5.5,  2.58,
         -0.45, 2.58, 2.0
    };
    
    printf("Entry matrix A: \n");
    print_dmatrix(a, NRHS, NRHS);

    
    double eigvec[NRHS*NRHS];
    double eigval[NRHS];
    
    eigen_decomp(NRHS, a, eigvec, eigval);
    
    printf("eigenvalues(D) are: \n");
    print_dmatrix(eigval, 1, NRHS);
    
    printf("eigenvectors(U) are: \n");
    print_dmatrix(eigvec, NRHS, NRHS);
    
    
    double tmp_a[NRHS*NRHS];
    dmat_B_ATDA(NRHS, NRHS, eigvec, eigval, tmp_a);
    
    
    printf("U' * D * U is: \n ");
    print_dmatrix(tmp_a, NRHS, NRHS);
    
    
}

void test_C_ATA()
{
    double dtwo = 2.0;
    double a[LDB*N] = {
        4.02, -1.56, 9.81, 6.19,  4.00,
        -4.09, -8.22, -8.67, -4.57, -7.57,
        1.75, -8.61, -3.03,  2.86, 8.99
    };
    
    printf("Entry matrix A: \n");
    print_dmatrix(a, LDB, N);
    
    double c[N*N];
    dmat_C_ATA(LDB, N, &dtwo, a, &dtwo, c);
    
    printf("C is : \n");
    print_dmatrix(c, N, N);
    
    printf("A after entry is: \n");
    print_dmatrix(a, LDB, N);
}

void test_C_AB()
{
    double a[LDA*N] = {
        6.80, -6.05, -0.45,  8.32, -9.67,
        -2.11, -3.30,  2.58,  2.71, -5.14,
        5.66, 5.36, -2.70,  4.35, -7.26,
        5.97, -4.44,  0.27, -7.17, 6.08,
        8.23, 1.08,  9.04,  2.14, -6.87
    };
    double b[N*LDB] = {
        4.02, -1.56, 9.81,
        6.19,  4.00, -4.09,
        -8.22, -8.67, -4.57,
        -7.57, 1.75, -8.61,
        -3.03,  2.86, 8.99
    };
    
    double c[LDA*LDB];
    
    printf("Entry matrix A: \n");
    print_dmatrix(a, LDA, N);
    printf("Entry matrix B: \n");
    print_dmatrix(b, N, LDB);
    
    dmat_C_AB(N, LDA, LDB, &done, a, b, c);
    
    printf("A * B is: \n");
    print_dmatrix(c, LDA, LDB);
}


void test_eigen_decomp()
{
    double a[LDA*N] = {
        6.80, -6.05, -0.45,  8.32, -9.67,
        -2.11, -3.30,  2.58,  2.71, -5.14,
        5.66, 5.36, -2.70,  4.35, -7.26,
        5.97, -4.44,  0.27, -7.17, 6.08,
        8.23, 1.08,  9.04,  2.14, -6.87
    };
    print_dmatrix(a, 5, 5);
    double c[25], b[25];
    eigen_decomp(5, a, b, c);
    
    print_dmatrix(a, 5, 5);
}

void test_liag()
{
    // test_C_ATB();
    // test_C_ABT();
    // test_eigen();
    // test_C_ATA();
    // test_C_AB();
    test_eigen_decomp();
}