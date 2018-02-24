#include <iostream>
#include <iomanip>
//        std::cout << std::scientific << std::setprecision(14);
#include <limits>  // std::numeric_limits
#include <assert.h>
#include <stdlib.h>
#include <cmath>
#include "wave_number.h"
#include "secular.h"
#include "frequency_equation.h"
//timoshenko beam,  implicit relationship between the wave numbers
//p.951, s, gamma2 are derived from the parameters in equations (77) & (81)
// This is the residual for the wave equation constraint ...
// that must be satisfied for (a,b) pairs in terms of k and gamma2.
double wave_number_residual(double k,double gamma2,double a,double b, bool subcritical)
{
    double a2 = a*a;
    double b2 = b*b;
    double k2 = k*k;
    double g2 = gamma2;
    double g4 = gamma2*gamma2;
    double a4 = a2*a2;
    double b4 = b2*b2;
    double residual = 0.;
    if( subcritical )
        residual = -k2*g2*a4 - k2*(g4+1.)*a2*b2 + a2*(g2+1.) - k2*g2*b4 - (g2+1.)*b2;
    else
        residual = -k2*g2*a4 + k2*(g4+1.)*a2*b2 + a2*(g2+1.) - k2*g2*b4 + (g2+1.)*b2;
    return residual;
}

double euler_wave_number( unsigned id )
{
    double freefree[] = {4.730040744862704, 7.853204624095838,
     10.99560783800169, 14.13716549125746, 17.27875965739948};
    if( id < 5 ) return freefree[id];

    double lambda = (1.5+ static_cast<double>(id))*M_PI;

    unsigned max_id = static_cast<unsigned>( log( std::numeric_limits<double>::max() )/M_PI  - 1.5);
    if( id >= max_id ) return lambda;

    double r = 1.;
    double dr =1.;
    double h = 1.;
    unsigned max_iter = id;
    unsigned iteration = 0;

    while( (iteration < max_iter) && (fabs(r) > fabs(dr)*1.e-13) && (fabs(h) > 1.e-14) ) 
    {
        r = cos(lambda)*cosh(lambda)-1.;
        dr = cos(lambda)*sinh(lambda) - sin(lambda)*cosh(lambda);
        h = -r/dr;
        lambda += h;
        iteration = iteration + 1;
    }
    if( iteration == max_iter ) 
    {
        std::cerr << " fatal error: unable to accurately determine Euler wave number\n";
        exit(-1);
    }
    return lambda;
}
