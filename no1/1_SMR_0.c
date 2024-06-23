// 小课题1 谱方法 rho1 = 0 Romberg method
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <lapacke.h>
#include <gsl/gsl_sf_bessel.h>

const double error = 1e-15;
const int N = 9;
const double length = 1;
const double PI = 3.14159265358979323;

double m_f(int m, int n, double x);
double m_f_add(int m, int n, int k, double initial_point, double t);
double m_romberg(int m, int n, double initial_point, double terminal_point, const double error);

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
        for(int j=0; j<N; j++){
            A[i*N+j] = ((j+1)*PI)*((j+1)*PI)*m_romberg(i+1,j+1,0,1,error)/m_romberg(i+1,i+1,0,1,error);
            //printf("%lf\n", m_romberg(i+1,j+1,0,1,error));
            //printf("%d,%lf\n", i+1, m_romberg(i+1,i+1,0,1,error));
        }
        //printf("%lf\n", A[i*N+i]);
    }
    //printf("%d\n", (3+1)/2);
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            if(abs(A[i*N+j])<1e-6)A[i*N+j] = 0;
        }
    }
/*
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            printf("%lf ", A[i*N+j]);
        }printf("\n");
    }
*/
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

double m_f(int m, int n, double x){
    return sin(m*PI*x)*sin(n*PI*x);
}

double m_f_add(int m, int n, int k, double initial_point, double t){
    double a = 0;
    for(int i = 1; i <= (k>2? pow(2,k-2):1); i++){
        a += m_f(m,n,initial_point+(2*i-1)*t/pow(2,k-1));
    }
    return a*t/(k>2? pow(2,k-2):1);
}

double m_romberg(int m, int n, double initial_point, double terminal_point, double error){
    double answer = 0;
    double t = terminal_point-initial_point;
    int M = 40;     //the initial M, if it's not enough, we'll give you an alert and you need to enlarge M
    double *Px;
    Px = (double *)malloc(M*M*sizeof(double));

    Px[0] = ( m_f(m,n,initial_point)+m_f(m,n,terminal_point) )*t/2.0;
    int i = 1, j = 1;
    while(i){
        Px[i*M+0] = (Px[(i-1)*M+0]+m_f_add(m,n,i+1,initial_point,t))/2;
        while(j){
            Px[i*M+j] = Px[i*M+j-1]+(Px[i*M+j-1]-Px[(i-1)*M+j-1])/(pow(4,j)-1);
            if(j == i){
                j = 1;
                break;
            }else{
                j++;
            }
        }
        double e = fabs(Px[i*M+i]-Px[(i-1)*M+i-1]);
        if(e<error&&i>10){
            answer = Px[i*M+i];
            break;
        }else if(i == M-1){
            printf("Please enlarge Px!!!\n");
        }else{
            i++;
        }
    }

    free(Px);
    return answer;
}