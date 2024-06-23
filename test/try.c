double m_mmn(int E_m1, double E_x1, int E_m2, double E_x2, double r, double th){
    return gsl_sf_bessel_Jn(E_m1,E_x1*r)*cos(E_m1*th)*(1/(1+gsl_sf_bessel_Jn(1,3.831706*r)*cos(th)))*gsl_sf_bessel_Jn(E_m2,E_x2*r)*cos(E_m2*th);
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
        if(e<error&&i>10){
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