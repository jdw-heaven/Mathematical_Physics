// 小课题1 谱方法 rho1 = 0 excat method
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <lapacke.h>
#include <gsl/gsl_sf_bessel.h>

const double error = 1e-10;
const int N = 1999;
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

    for(int i=0; i<N; i++){
        A[i*N+i] = ((i+1)*PI)*((i+1)*PI);
    }
    double *lambda;
    lambda = (double *)calloc(N,sizeof(double));

    LAPACKE_dsyev(LAPACK_COL_MAJOR,'V','U',N,A,N,lambda);
    for(int i=0; i<8; i++){
        fprintf(eigenvalue, "%lf\n", sqrt(lambda[i]));
    }
    double ax = 0.0;
    double *u;
    u = (double *)calloc(8,sizeof(double));
    for(int i=1; i<N+1; i++){
        ax = i*dx;
        for(int j=0; j<8; j++){
            u[j] = 0.0;
        }
        for(int j=0; j<N; j++){
            u[0] += A[0*N+j]*sin((j+1)*PI*ax);
            u[1] += A[1*N+j]*sin((j+1)*PI*ax);
            u[2] += A[2*N+j]*sin((j+1)*PI*ax);
            u[3] += A[3*N+j]*sin((j+1)*PI*ax);
            u[4] += A[4*N+j]*sin((j+1)*PI*ax);
            u[5] += A[5*N+j]*sin((j+1)*PI*ax);
            u[6] += A[6*N+j]*sin((j+1)*PI*ax);
            u[7] += A[7*N+j]*sin((j+1)*PI*ax);
        }
        fprintf(x, "%lf\n", i*dx);
        fprintf(u1, "%lf\n", u[0]);
        fprintf(u2, "%lf\n", u[1]);
        fprintf(u3, "%lf\n", u[2]);
        fprintf(u4, "%lf\n", u[3]);
        fprintf(u5, "%lf\n", u[4]);
        fprintf(u6, "%lf\n", u[5]);
        fprintf(u7, "%lf\n", u[6]);
        fprintf(u8, "%lf\n", u[7]);
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
    free(u);
    return 0;
}