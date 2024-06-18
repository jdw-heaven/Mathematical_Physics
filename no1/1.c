#include "../include/header.h"
#include <lapacke.h>
#include <gsl/gsl_sf_bessel.h>

const double error = 1e-10;

int main(void)
{
    double x = 3.831706; // 可以替换为任意实数
    printf("Bessel function of the first kind of order 0: %f\n", gsl_sf_bessel_J0(x));
    printf("Bessel function of the first kind of order 1: %f\n", j1(x));
}

/*被积函数
函数声明在头文件中
警告：可以改变返回值的形式，但是不能改变函数原型
参数请使用全局变量的形式加入*/
double m_function(double x){
    return cos(x)*cos(x);
}