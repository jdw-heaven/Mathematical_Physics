// 小课题 2 谱方法 tho1 = 0.2 phi_2 有简并
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <lapacke.h>
#include <gsl/gsl_sf_bessel.h>

const double error = 1e-10;
const int N = 16;    // N可取1-16,即本征解的个数，受限于谱方法
const double length = 1;
const double PI = 3.14159265358979323;

double m_B(int E_m1, double E_x1, int E_m2, double E_x2, double r);
double m_B_add(int E_m1, double E_x1, int E_m2, double E_x2, int k, double initial_point, double t);
double m_B_romberg(int E_m1, double E_x1, int E_m2, double E_x2, double initial_point, double terminal_point, double error);

double m_T(int E_m1, int E_m2, double th);
double m_T_add(int E_m1, int E_m2, int k, double initial_point, double t);
double m_T_romberg(int E_m1, int E_m2, double initial_point, double terminal_point, double error);

double m_mmn(int E_m1, double E_x1, int E_m2, double E_x2, double r, double th);
double m_Mmn_th_add(int E_m1, double E_x1, int E_m2, double E_x2, double r, int k, double th1, double t);
double m_Mmn_th(int E_m1, double E_x1, int E_m2, double E_x2, double r, double error);
double m_Mmn_add(int E_m1, double E_x1, int E_m2, double E_x2, int k, double r1, double t, double error);
double m_Mmn(int E_m1, double E_x1, int E_m2, double E_x2, double error);

double m_f(int E_m, double th);

void m_search(double *lambda, int *num, int N);

int main(void)
{
    FILE *eigenvalue;
    FILE *x;
    FILE *y;
    const int num_files = N;
    FILE *u[num_files];
    for (int i = 0; i < num_files; i++) {
        char filename[20];  // 假设文件名长度不超过19个字符
        sprintf(filename, "data/u%d.txt", i+1); // 生成文件名
        u[i] = fopen(filename, "w"); // 打开文件
        if (u[i] == NULL) {
            perror("Error opening file");
            exit(EXIT_FAILURE);
        }
    }


    eigenvalue = fopen("data/eigenvalue.txt", "w");
    x = fopen("data/x.txt", "w");
    y = fopen("data/y.txt", "w");

    int E_m[16] = {0,1,1,2,2,0,3,3,1,1,4,4,2,2,0,5};
    double E_x[16] = {2.404826,3.831706,3.831706,5.135622,5.135622,5.520078,6.380162,6.380162,7.015587,7.015587
    ,7.588342,7.588342,8.417244,8.417244,8.653728,8.771484};
/*
    for(int i=0; i<16; i++){
        printf("%d %lf\n", E_m[i], E_x[i]);
    }
*/
    double *A;
    A = (double *)calloc(N*N,sizeof(double));
    double *Ap;
    Ap = (double *)calloc(N*N,sizeof(double));

    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            A[i*N+j] = m_Mmn(E_m[i],E_x[i],E_m[j],E_x[j],error)*E_x[i]*E_x[i]/m_B_romberg(E_m[i],E_x[i],E_m[i],E_x[i],0,1,error)/m_T_romberg(E_m[i],E_m[i],0,2*PI,error);
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
    double *lambda_i;
    lambda_i = (double *)calloc(N,sizeof(double));
    int *num;
    num = (int *)calloc(N,sizeof(int));
    LAPACKE_dgeev(LAPACK_COL_MAJOR, 'N', 'V', N, A, N, lambda, lambda_i, NULL, N, Ap, N);
    //对特征值进行排序
    m_search(lambda,num,N);
/*
    for(int i=0; i<N; i++){
        printf("%lf ", lambda[num[i]]);
    }
*/
    double r[101] = {0,};
    double th[301] = {0,};
    for(int i=0; i<101; i++){
        r[i] = i*(1/100.0);
    }
    for(int i=0; i<301; i++){
        th[i] = i*(2*PI/300.0);
    }
    for(int j=0; j<101; j++){
        for(int k=0; k<301; k++){
            fprintf(x, "%lf\n", r[j]*cos(th[k]));
            fprintf(y, "%lf\n", r[j]*sin(th[k]));
        }
    }
    double *nu;
    nu = (double *)calloc(N,sizeof(double));
    for(int i=0; i<N; i++){
        fprintf(eigenvalue, "%lf\n", sqrt(lambda[num[i]]));
        for(int j=0; j<101; j++){
            for(int k=0; k<301; k++){
                for(int l=0; l<N; l++){
                    nu[i] += Ap[num[i]*N+l]*gsl_sf_bessel_Jn(E_m[l],E_x[l]*r[j])*m_f(E_m[l],th[k]);
                }
                fprintf(u[i], "%lf\n", nu[i]);
                nu[i] = 0;
            }
        }
    }

    free(A);
    free(Ap);
    free(num);
    free(lambda);
    free(lambda_i);
    free(nu);
    return 0;
}

double m_B(int E_m1, double E_x1, int E_m2, double E_x2, double r){
    return gsl_sf_bessel_Jn(E_m1,E_x1*r)*gsl_sf_bessel_Jn(E_m2,E_x2*r)*r;
}
double m_B_add(int E_m1, double E_x1, int E_m2, double E_x2, int k, double initial_point, double t){
    double a = 0;
    for(int i = 1; i <= (k>2? pow(2,k-2):1); i++){
        a += m_B(E_m1,E_x1,E_m2,E_x2,initial_point+(2*i-1)*t/pow(2,k-1));
    }
    return a*t/(k>2? pow(2,k-2):1);
}
double m_B_romberg(int E_m1, double E_x1, int E_m2, double E_x2, double initial_point, double terminal_point, double error){
    double answer = 0;
    double t = terminal_point-initial_point;
    int M = 40;     //the initial M, if it's not enough, we'll give you an alert and you need to enlarge M
    double *Px;
    Px = (double *)malloc(M*M*sizeof(double));

    Px[0] = ( m_B(E_m1,E_x1,E_m2,E_x2,initial_point)+m_B(E_m1,E_x1,E_m2,E_x2,terminal_point) )*t/2.0;
    int i = 1, j = 1;
    while(i){
        Px[i*M+0] = (Px[(i-1)*M+0]+m_B_add(E_m1,E_x1,E_m2,E_x2,i+1,initial_point,t))/2;
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
double m_f(int E_m, double th){
    double a = cos(E_m*th);
    if(E_m==2||E_m==4||E_m==7||E_m==9||E_m==11||E_m==13)a = sin(E_m*th);
    return a;
}
double m_T(int E_m1, int E_m2, double th){
    double a = m_f(E_m1,th)*m_f(E_m2,th);
    return a;
}
double m_T_add(int E_m1, int E_m2, int k, double initial_point, double t){
    double a = 0;
    for(int i = 1; i <= (k>2? pow(2,k-2):1); i++){
        a += m_T(E_m1,E_m2,initial_point+(2*i-1)*t/pow(2,k-1));
    }
    return a*t/(k>2? pow(2,k-2):1);
}
double m_T_romberg(int E_m1, int E_m2, double initial_point, double terminal_point, double error){
    double answer = 0;
    double t = terminal_point-initial_point;
    int M = 40;     //the initial M, if it's not enough, we'll give you an alert and you need to enlarge M
    double *Px;
    Px = (double *)malloc(M*M*sizeof(double));

    Px[0] = ( m_T(E_m1,E_m2,initial_point)+m_T(E_m1,E_m2,terminal_point) )*t/2.0;
    int i = 1, j = 1;
    while(i){
        Px[i*M+0] = (Px[(i-1)*M+0]+m_T_add(E_m1,E_m2,i+1,initial_point,t))/2;
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

void m_search(double *lambda, int *num, int N){
    int m=0;
    for(int i=0; i<N; i++){
        num[i] = i;
    }
    for(int i=1; i<N; i++){
        for(int j=i-1; j>=0; j--){
            if(lambda[num[j]]>lambda[num[j+1]]){
                m = num[j];
                num[j] = num[j+1];
                num[j+1] = m;
            }
        }
    }
}

double m_mmn(int E_m1, double E_x1, int E_m2, double E_x2, double r, double th){
    return gsl_sf_bessel_Jn(E_m1,E_x1*r)*m_f(E_m1,th)*(1/(1+0.2*gsl_sf_bessel_Jn(1,3.831706*r)*cos(th)))*gsl_sf_bessel_Jn(E_m2,E_x2*r)*m_f(E_m2,th)*r;
}
//The inner layer of the double integral
double m_Mmn_th_add(int E_m1, double E_x1, int E_m2, double E_x2, double r, int k, double th1, double t){
    double a = 0;
    for(int i = 1; i <= (k>2? pow(2,k-2):1); i++){
        a += m_mmn(E_m1,E_x1,E_m2,E_x2,r,th1+(2*i-1)*t/pow(2,k-1));
    }
    return a*t/(k>2? pow(2,k-2):1);
}
double m_Mmn_th(int E_m1, double E_x1, int E_m2, double E_x2, double r, double error){
    double answer = 0;
    double th1 = 0, th2 = 2*PI;
    double t = th2-th1;
    int M = 40;     //the initial M, if it's not enough, we'll give you an alert and you need to enlarge M
    double *Px;
    Px = (double *)malloc(M*M*sizeof(double));

    Px[0] = ( m_mmn(E_m1,E_x1,E_m2,E_x2,r,th1)+m_mmn(E_m1,E_x1,E_m2,E_x2,r,th2) )*t/2;
    int i = 1, j = 1;
    while(i){
        Px[i*M+0] = (Px[(i-1)*M+0]+m_Mmn_th_add(E_m1,E_x1,E_m2,E_x2,r,i+1,th1,t))/2;
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
        if(e<error&&i>5){
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
double m_Mmn_add(int E_m1, double E_x1, int E_m2, double E_x2, int k, double r1, double t, double error){
    double a = 0;
    for(int i = 1; i <= (k>2? pow(2,k-2):1); i++){
        a += m_Mmn_th(E_m1,E_x1,E_m2,E_x2,r1+(2*i-1)*t/pow(2,k-1),error);
    }
    return a*t/(k>2? pow(2,k-2):1);
}
double m_Mmn(int E_m1, double E_x1, int E_m2, double E_x2, double error){
    double answer = 0;
    double r1 = 0.0, r2 = 1.0;
    double t = r2-r1;
    int M = 40;
    double *Hx;
    Hx = (double *)malloc(M*M*sizeof(double));

    Hx[0] = ( m_Mmn_th(E_m1,E_x1,E_m2,E_x2,r1,error)+m_Mmn_th(E_m1,E_x1,E_m2,E_x2,r2,error) )*t/2;
    int i = 1, j = 1;
    while(i){
        Hx[i*M+0] = (Hx[(i-1)*M+0]+m_Mmn_add(E_m1,E_x1,E_m2,E_x2,i+1,r1,t,error))/2;
        while(j){
            Hx[i*M+j] = Hx[i*M+j-1]+(Hx[i*M+j-1]-Hx[(i-1)*M+j-1])/(pow(4,j)-1);
            if(j == i){
                j = 1;
                break;
            }else{
                j++;
            }
        }
        double e = fabs(Hx[i*M+i]-Hx[(i-1)*M+i-1]);
        if(e<error&&i>5){
            answer = Hx[i*M+i];
            break;
        }else if(i == M-1){
            printf("Please enlarge Hx!!!\n");
        }else{
            i++;
        }
    }

    free(Hx);
    return answer;
}