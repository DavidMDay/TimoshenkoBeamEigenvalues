#include <math.h>
#include <limits>  // std::numeric_limits
#include <stdlib.h>
#include "frequency_equation_cc.h"
#include "frequency_equation_cf.h"
#include "frequency_equation.h"   // for boundary_condition enum
//timoshenko beam,  a<ac, implicit relationship between the wave numbers
// page 951, a_c  equation (79), g2 : equation (81),
// subcritical  a < ac
double frequency_equation_sub(double a, double b, double gamma2, double* F_a, double* F_b)
{
    double a2 = a*a;
    double b2 = b*b;
    double ab = a*b;
    double g2 = gamma2;

    double numerator = (a2 - b2)*(a2 + b2 + (g2-1.)*ab)*(a2 + b2 - (g2-1.)*ab);
    double num_a=
    2*a*(a2 + b2 + (g2-1.)*ab)*(a2 + b2 - (g2-1.)*ab)+
    (a2 - b2)*(2*a + (g2-1.)*b)*(a2 + b2 - (g2-1.)*ab)+
    (a2 - b2)*(a2 + b2 + (g2-1.)*ab)*(2.*a - (g2-1.)*b);

    double num_b = 
    -2*b*(a2 + b2 + (g2-1.)*ab)*(a2 + b2 - (g2-1.)*ab)+
    (a2 - b2)*(2*b + (g2-1.)*a)*(a2 + b2 - (g2-1.)*ab)+
    (a2 - b2)*(a2 + b2 + (g2-1.)*ab)*(2*b - (g2-1.)*a);

    double denominator = 2.*ab*(g2*a2+b2)*(a2 + g2*b2);  // 2 g2 a5 /b

    double den_a = 2.*b*(
    (b2 + g2* a2)*(a2 + g2*b2)+
    2.*a2*g2*(a2 + g2*b2)+
    2.*a2*(b2 + g2*a2));

    double den_b = 2.*a*(
     (g2*a2+b2)*(a2 + g2*b2)+
     2.*b2*(a2 + g2*b2)+
     2.*g2*b2*(b2 + g2*a2));   // 2 g2 a5

    double quotient = numerator/denominator;

    double quotient_a = (num_a*denominator - den_a*numerator)/( denominator*denominator );
    double quotient_b = (num_b*denominator - den_b*numerator)/( denominator*denominator );

    double F =      quotient*sin(a)*sinh(b) - cos(a)*cosh(b)+1.;
    *F_a =  quotient_a*sin(a)*sinh(b) + quotient*cos(a)*sinh(b) + sin(a)*cosh(b);
    *F_b =  quotient_b*sin(a)*sinh(b) + quotient*sin(a)*cosh(b) - cos(a)*sinh(b);

    return F;
}
//timoshenko beam,  a<ac, implicit relationship between the wave numbers
// page 951, a_c  equation (79), g2 : equation (81),
// super critical


double frequency_equation_super(double a, double b, double gamma2, double* F_a, double* F_b)
{
    double a2 = a*a;
    double b2 = b*b;
    double ab = a*b;
    double g2 = gamma2;

    double numerator = (a2+b2)*( (a2-b2)*(a2-b2)  + (g2-1.)*(g2-1.)*a2*b2 );

    double num_a=
    (2.*a)*( (a2-b2)*(a2-b2)  + (g2-1.)*(g2-1.)*a2*b2 )+
    (a2+b2)*( 2.*(a2-b2)*2.*a  + 2.*(g2-1.)*(g2-1.)*a*b2 );

    double num_b = 
    (2.*b)*( (a2-b2)*(a2-b2)  + (g2-1.)*(g2-1.)*a2*b2 ) +
    (a2+b2)*( -2.*(a2-b2)*2.*b  + 2.*(g2-1.)*(g2-1.)*a2*b );

    double denominator = 2.*ab*(g2*a2-b2)*(a2 - g2*b2);

    double den_a = 2.*b*(
    (g2*a2-b2)*(a2 - g2*b2)+
    (2.*g2*a2)*(a2 - g2*b2)+
    (g2*a2-b2)*(2.*a2));

    double den_b = 2.*a*(
    (g2*a2-b2)*(a2 - g2*b2)
    -2.*b2*(a2 - g2*b2)
    -2.*(g2*a2-b2)*g2*b2);

    double quotient = numerator/denominator;
    double quotient_a = (num_a*denominator - den_a*numerator)/( denominator*denominator );
    double quotient_b = (num_b*denominator - den_b*numerator)/( denominator*denominator );

    double F =      quotient*sin(a)*sin(b) - cos(a)*cos(b)+1.;
    *F_a =  quotient_a*sin(a)*sin(b) + quotient*cos(a)*sin(b) + sin(a)*cos(b);
    *F_b =  quotient_b*sin(a)*sin(b) + quotient*sin(a)*cos(b) + cos(a)*sin(b);
    return F;
}
//timoshenko beam,  implicit relationship between the wave numbers

double frequency_equation(double a,
                          double b,
                          double gamma2,
                          bool is_sub_critical,
                          boundary_condition bc,
                          double* F_a,
                          double* F_b)
{
    double F = 0;
    if( bc == freefree )
    {
        if( is_sub_critical )
        {
            F = frequency_equation_sub(a, b, gamma2, F_a, F_b);
        }
        else
        {
            F = frequency_equation_super(a, b, gamma2, F_a, F_b);
        }
    }
    else
    {
        if( bc == clampedclamped )
            F = frequency_equation_cc(a, b, gamma2, is_sub_critical,F_a, F_b);
        else // clampedfree
            F = frequency_equation_cf(a, b, gamma2, is_sub_critical,F_a, F_b);
    }
    return F;
}


double euler_wave_number( unsigned id, boundary_condition bc )
{
    double eulerfreefree[] = {4.730040744862704, 7.853204624095838,
     10.99560783800169, 14.13716549125746, 17.27875965739948};
    double eulerclampedfree[] = {1.87510406871196, 4.694091132974175,
     7.854757438237613, 10.99554073487547, 14.13716839104647};

    if( id < 5 )
    {
        if( bc != clampedfree ) return eulerfreefree[id];
        else return eulerclampedfree[id];
    }

    double offset = 0.;
    if( bc == clampedfree ) 
        offset = .5;
    else
        offset = 1.5;

    double lambda = (offset +  static_cast<double>(id))*M_PI;

    unsigned max_id = static_cast<unsigned>( log( std::numeric_limits<double>::max() )/M_PI  - 1.5);
    if( id >= max_id ) return lambda;

    double r = 1.;
    double dr =1.;
    double h = 1.;
    unsigned max_iter = id;
    unsigned iteration = 0;

    double sign = 0.;
    if( bc == clampedfree ) 
        sign = 1.;
    else
        sign = -1;

    while( (iteration < max_iter) && (fabs(r) > fabs(dr)*1.e-13) && (fabs(h) > 1.e-14) )
    {
        
        r = cos(lambda)*cosh(lambda) + sign;
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
