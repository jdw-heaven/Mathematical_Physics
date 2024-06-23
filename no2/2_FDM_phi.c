// 小课题 2 有限差分法 \rho _{1} = 0.2\phi _{2}(\rho , \theta)
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <lapacke.h>
#include <gsl/gsl_sf_bessel.h>

const double error = 1e-10;
const int N = 31;   //为了避免奇点，取奇数
const double length = 1;
const double PI = 3.14159265358979323;

void m_search(double *lambda, int *num, int N);

int main(void)
{
    FILE *eigenvalue;
    FILE *x;
    FILE *y;
    const int num_files = 16;
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


    int D = (N+2)*(N+2);
    double dx = 2*length/(N+1);
    double *A;
    A = (double *)calloc(D*D,sizeof(double));
    //判断点是否在院内并赋值
    int row=0, c1=0, c2=0, c3=0, c4=0;
    int *nz;
    nz = (int *)calloc(D,sizeof(int));
    int nzn = 0;
    double ax = 0.0,ay = 0.0;
    for(int i=1; i<=N; i++){
        for(int j=1; j<=N; j++){
            ax = i-N/2.0;
            ay = j-N/2.0;
            //printf("%lf,%lf\n",ax,ay);
            double r = sqrt(ax*ax+ay*ay);
            double t=0.0,c=0.0;
            if(r<N/2){
                //printf("%lf\n",r);
                row = i*(N+2)+j;
                t = ax/r;
                r = sqrt(ax*ax+ay*ay)*dx;
                c = 1/(1+0.2*gsl_sf_bessel_J1(3.831706*r)*t);
                //printf("%lf\n",0.2*gsl_sf_bessel_J1(3.831706*r)*t);
                //if(abs(c)>3)printf("%lf,%lf\n",ax,ay);
                c1 = (i-1)*(N+2)+j;
                c2 = (i+1)*(N+2)+j;
                c3 = i*(N+2)+j-1;
                c4 = i*(N+2)+j+1;
                A[row*D+row] = 4*c;
                A[row*D+c1] = -1*c;
                A[row*D+c2] = -1*c;
                A[row*D+c3] = -1*c;
                A[row*D+c4] = -1*c;
                nz[nzn] = row;
                nzn += 1;
                fprintf(x, "%lf\n", (i-N/2)*dx);
                fprintf(y, "%lf\n", (j-N/2)*dx);
            }
        }
    }
/*
    for(int i=0; i<D; i++){
        for(int j=0; j<D; j++){
            printf("%lf ", A[i*D+j]);
        }printf("\n");
    }
*/

    double *Ap;
    Ap = (double *)calloc(nzn*nzn,sizeof(double));
    for(int i=0; i<nzn; i++){
        for(int j=0; j<nzn; j++){
            Ap[i*nzn+j] = A[nz[i]*D+nz[j]];
        }
    }
/*
    for(int i=0; i<nzn; i++){
        for(int j=0; j<nzn; j++){
            printf("%lf ", Ap[i*nzn+j]);
        }printf("\n");
    }
*/
    double *App;
    App = (double *)calloc(nzn*nzn,sizeof(double));
    double *lambda;
    lambda = (double *)calloc(nzn,sizeof(double));
    double *lambda_i;
    lambda_i = (double *)calloc(nzn,sizeof(double));
    LAPACKE_dgeev(LAPACK_COL_MAJOR, 'N', 'V', nzn, Ap, nzn, lambda, lambda_i, NULL, nzn, App, nzn);
    //LAPACKE_dsyev(LAPACK_COL_MAJOR,'V','U',nzn,Ap,nzn,lambda);
    int *num;
    num = (int *)calloc(16, sizeof(int));
    m_search(lambda,num,nzn);

    for(int i=0; i<16; i++){
        fprintf(eigenvalue, "%lf\n", sqrt(lambda[num[i]]/dx/dx));
        for(int j=0; j<nzn; j++){
            fprintf(u[i], "%lf\n", App[num[i]*nzn+j]);
        }
    }

    free(A);
    free(Ap);
    free(App);
    free(lambda);
    free(lambda_i);
    free(num);
    return 0;
}

void m_search(double *lambda, int *num, int N){
    int m=0;
    int z=0,zp=0;
    do{
        if(lambda[z]>0){
            num[zp] = z;
            zp += 1;
        }
        z+=1;
    }while(zp<16);
    for(int i=1; i<16; i++){
        for(int j=i-1; j>=0; j--){
            if(lambda[num[j]]>lambda[num[j+1]]){
                m = num[j];
                num[j] = num[j+1];
                num[j+1] = m;
            }
        }
    }
    for(int i=z; i<N; i++){
        if(lambda[i]<lambda[num[15]]&&lambda[i]>0){
            num[15] = i;
            for(int j=14; j>=0; j--){
                if(lambda[num[j]]>lambda[num[j+1]]){
                    m = num[j];
                    num[j] = num[j+1];
                    num[j+1] = m;
                }
            }
        }
    }
}