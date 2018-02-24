#include <iostream>
#include <math.h>
#include "frequency_equation_cc.h"
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

double frequency_equation(double a, double b, double gamma2, bool is_sub_critical, bool is_ff, double* F_a, double* F_b)
{
    double F = 0;
    if( is_ff )
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
        F = frequency_equation_cc(a, b, gamma2, is_sub_critical,F_a, F_b);
    }
    return F;
}
