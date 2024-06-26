// 小课题1 有限差分法 rho1 = 0
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <lapacke.h>
#include <gsl/gsl_sf_bessel.h>

const double error = 1e-10;
const int N = 99;
const double length = 1;
const double PI = 3.14159265358979323;

int main(void)
{
    FILE *eigenvalue;
    FILE *x;
    FILE *u1;
    FILE *u2;
    FILE *u3;
    FILE *u4;
    FILE *u5;
    FILE *u6;
    FILE *u7;
    FILE *u8;


    eigenvalue = fopen("data/eigenvalue.txt", "w");
    x = fopen("data/x.txt", "w");
    u1 = fopen("data/u1.txt", "w");
    u2 = fopen("data/u2.txt", "w");
    u3 = fopen("data/u3.txt", "w");
    u4 = fopen("data/u4.txt", "w");
    u5 = fopen("data/u5.txt", "w");
    u6 = fopen("data/u6.txt", "w");
    u7 = fopen("data/u7.txt", "w");
    u8 = fopen("data/u8.txt", "w");

    double dx = length/(N+1);
    double *A;
    A = (double *)calloc(N*N,sizeof(double));
    A[0] = 2.0;
    A[1] = -1.0;
    A[N*N-1] = 2.0;
    A[N*N-2] = -1.0;
    for(int i=1; i<N-1; i++){
        A[i*N+i] = 2.0;
        A[i*N+i-1] = -1.0;
        A[i*N+i+1] = -1.0;
    }
    double *lambda;
    lambda = (double*)calloc(N,sizeof(double));
/*
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            printf("%lf ", A[i*N+j]);
        }printf("\n");
    }
*/
    LAPACKE_dsyev(LAPACK_COL_MAJOR, 'V', 'U', N, A, N, lambda);
    for(int i=0; i<N; i++){
        fprintf(eigenvalue, "%lf\n", sqrt(lambda[i]/(dx*dx)));
        fprintf(x, "%lf\n", i*dx+dx);
        fprintf(u1, "%lf\n", A[0*N+i]);
        fprintf(u2, "%lf\n", A[1*N+i]);
        fprintf(u3, "%lf\n", A[2*N+i]);
        fprintf(u4, "%lf\n", A[3*N+i]);
        fprintf(u5, "%lf\n", A[4*N+i]);
        fprintf(u6, "%lf\n", A[5*N+i]);
        fprintf(u7, "%lf\n", A[6*N+i]);
        fprintf(u8, "%lf\n", A[7*N+i]);
    }

    fclose(eigenvalue);
    fclose(x);
    fclose(u1);
    fclose(u2);
    fclose(u3);
    fclose(u4);
    fclose(u5);
    fclose(u6);
    fclose(u7);
    fclose(u8);

    free(A);
    free(lambda);
    return 0;
}