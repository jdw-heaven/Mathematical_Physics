// 小课题1 有限差分法 rho1=0.3*sin(pi*x)
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <lapacke.h>
#include <gsl/gsl_sf_bessel.h>

const double error = 1e-10;
const int N = 99;
const double length = 1;
const double PI = 3.14159265358979323;

double m_f(double x);   // f(x)，非线性Helmholtz方程
void m_search(double *lambda, int *num, const int N); //查找最小的8个特征值

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
    double *Ap;
    Ap = (double *)calloc(N*N,sizeof(double));
    A[0] = 2.0*m_f(dx);
    A[1] = -1.0*m_f(dx);
    A[N*N-1] = 2.0*m_f(N*dx);
    A[N*N-2] = -1.0*m_f(N*dx);
    for(int i=1; i<N-1; i++){
        A[i*N+i] = 2.0*m_f((i+1)*dx);
        A[i*N+i-1] = -A[i*N+i]/2;
        A[i*N+i+1] = A[i*N+i-1];
    }
    double *lambda;
    lambda = (double*)calloc(N,sizeof(double));
    double *lambda_i;
    lambda_i = (double*)calloc(N,sizeof(double));
/*
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            printf("%lf ", A[i*N+j]);
        }printf("\n");
    }
*/
    LAPACKE_dgeev(LAPACK_COL_MAJOR, 'N', 'V', N, A, N, lambda, lambda_i, NULL, N, Ap, N);
    //对特征值进行排序，找到前8小的
    int *num;
    num = (int *)calloc(8, sizeof(int));
    m_search(lambda,num,N);

    for(int i=0; i<8; i++){
        fprintf(eigenvalue, "%lf\n", sqrt(lambda[num[i]]/(dx*dx)));
    }

    for(int i=0; i<N; i++){
        fprintf(x, "%lf\n", i*dx+dx);
        fprintf(u1, "%lf\n", Ap[num[0]*N+i]);
        fprintf(u2, "%lf\n", Ap[num[1]*N+i]);
        fprintf(u3, "%lf\n", Ap[num[2]*N+i]);
        fprintf(u4, "%lf\n", Ap[num[3]*N+i]);
        fprintf(u5, "%lf\n", Ap[num[4]*N+i]);
        fprintf(u6, "%lf\n", Ap[num[5]*N+i]);
        fprintf(u7, "%lf\n", Ap[num[6]*N+i]);
        fprintf(u8, "%lf\n", Ap[num[7]*N+i]);
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
    free(Ap);
    free(num);
    free(lambda);
    free(lambda_i);
    return 0;
}

double m_f(double x){
    return 1.0/(1+0.3*sin(PI*x));
}

void m_search(double *lambda, int *num, const int N){
    int m=0;
    for(int i=0; i<8; i++){
        num[i] = i;
    }
    for(int i=1; i<8; i++){
        for(int j=i-1; j>=0; j--){
            if(lambda[num[j]]>lambda[num[j+1]]){
                m = num[j];
                num[j] = num[j+1];
                num[j+1] = m;
            }
        }
    }
    for(int i=8; i<N; i++){
        if(lambda[i]<lambda[num[7]]){
            num[7] = i;
            for(int j=6; j>=0; j--){
                if(lambda[num[j]]>lambda[num[j+1]]){
                    m = num[j];
                    num[j] = num[j+1];
                    num[j+1] = m;
                }
            }
        }
    }
}