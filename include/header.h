#ifndef HEADER_H
#define HEADER_H

// packages
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// functions
/*Romberg_integral*/
double m_function(double x);// 函数定义在.c文件中。
double m_function_add(int k, double initial_point, double t){
    double a = 0;
    for(int i = 1; i <= (k>2? pow(2,k-2):1); i++){
        a += m_function(initial_point+(2*i-1)*t/pow(2,k-1));
    }
    return a*t/(k>2? pow(2,k-2):1);
}
double m_romberg(double initial_point, double terminal_point, double error){
    double answer = 0;
    double t = terminal_point-initial_point;
    int M = 40;     //the initial M, if it's not enough, we'll give you an alert and you need to enlarge M
    double *Px;
    Px = (double *)malloc(M*M*sizeof(double));

    Px[0] = ( m_function(initial_point)+m_function(terminal_point) )*t/2;
    int i = 1, j = 1;
    while(i){
        Px[i*M+0] = (Px[(i-1)*M+0]+m_function_add(i+1,initial_point,t))/2;
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
        if(e<error&&i!=1){
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

#endif // 宏定义的头文件名