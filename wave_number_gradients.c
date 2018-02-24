#include <iostream>
#include <stdlib.h> // for exit
double wave_equation_dfdk(double k,double gamma2,double a,double b,bool is_sub_critical)
{
    double a2 = a*a;
    double a4 = a2*a2;
    double b2 = b*b;
    double b4 = b2*b2;
    double g2 = gamma2;
    double g4 = g2*g2;
    double dfdk = 0.;
    if( is_sub_critical )
        dfdk = -2.*k*(g2*a4 +(g4+1.)*a2*b2+ g2*b4);
    else
        dfdk = -2.*k*(g2*a4 -(g4+1.)*a2*b2+ g2*b4);
    return dfdk;
}
double wave_equation_dfdb(double k,double gamma2,double a,double b,bool is_sub_critical)
{
    double a2 = a*a;
    double b2 = b*b;
    double g2 = gamma2;
    double g4 = gamma2*gamma2;
    double k2 = k*k;
    double dfdb = 0.;
    if( is_sub_critical )
        dfdb = -2.*b*(2.*k2*g2*b2+g2+1. + k2*a2*(g4+1));
    else
        dfdb =  2.*b*(-2.*k2*g2*b2+g2+1.+ k2*a2*(g4+1));
    return dfdb;
}
double wave_equation_dfda(double k,double gamma2,double a,double b,bool is_sub_critical)
{
    double a2 = a*a;
    double b2 = b*b;
    double g2 = gamma2;
    double g4 = gamma2*gamma2;
    double k2 = k*k;
    double dfda = 0.;
    if( is_sub_critical )
        dfda = 2.*a*(-2.*k2*g2*a2+(g2+1.-k2*(g4+1.)*b2));
    else
        dfda = 2.*a*(-2.*k2*g2*a2+(g2+1.+k2*(g4+1.)*b2));
    return dfda;
}

double wave_equation_dbdk(double dfdk, double dfdb )
{
    if( dfdb == 0.) std::cout << "wave_equation dbdk fatal error dfdb vanishes\n";
    if( dfdb == 0.) exit(-1);
    double dbdk = -dfdk/dfdb;
    return dbdk;
}

double wave_equation_dbda(double dfda, double dfdb )
{
    if( dfdb == 0.) std::cout << "wave_equation dbda fatal error dfdb vanishes\n";
    if( dfdb == 0.) exit(-1);
    double dbda = -dfda/dfdb;
    return dbda;
}
