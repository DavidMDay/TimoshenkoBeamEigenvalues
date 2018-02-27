#include <iostream>
#include <math.h>
//timoshenko beam,  a<ac, implicit relationship between the wave numbers
// page 951, a_c  equation (79), g2 : equation (81),
// subcritical  a < ac

double frequency_equation_cf_sub(double a, double b, double gamma2, double* F_a, double* F_b)
{
    double a2 = a*a;
    double b2 = b*b;
    double ab = a*b;
    double a4 = a2*a2;
    double b4 = b2*b2;
    double g2 = gamma2;
    double g4 = g2*g2;

    double numerator = ab*((1.+g4)*a4 + 4*g2*a2*b2 + (1.+g4)*b4);

    double b3 = b*b2;
    double b5 = b*b4;

    double num_a= 5.*(1.+g4)*a4*b + 12.*g2*a2*b3 + (1.+g4)*b5;

    double a3 = a*a2; 
    double a5 = a*a4;

    double num_b = (1.+g4)*a5 + 12.*g2*a3*b2 + 5.*(1.+g4)*a*b4;

    double denominator = (b2 + g2*a2)*(a2 + g2*b2);
    double den_a = 2.*a*((1.+g4)*b2+2.*g2*a2);
    double den_b = 2.*b*((1.+g4)*a2+2.*g2*b2);

    double quotient = numerator/denominator;
    double quotient_a = (num_a*denominator - den_a*numerator)/( denominator*denominator );
    double quotient_b = (num_b*denominator - den_b*numerator)/( denominator*denominator );

    double F = 2.*a*b+(b2-a2)*sin(a)*sinh(b) + quotient*cos(a)*cosh(b);

    *F_a = 2.*b -2.*a*sin(a)*sinh(b) + (b2-a2)*cos(a)*sinh(b) -quotient*sin(a)*cosh(b)+quotient_a*cos(a)*cosh(b);
    *F_b = 2.*a + 2.*b*sin(a)*sinh(b)  + quotient *cos(a)*sinh(b)+ (b2-a2)*sin(a)*cosh(b) + quotient_b*cos(a)*cosh(b);
    return F;
}
double frequency_equation_cf_super(double a, double b, double gamma2, double* F_a, double* F_b)
{
    double a2 = a*a;
    double b2 = b*b;
    double ab = a*b;
    double a4 = a2*a2;
    double b4 = b2*b2;
    double g2 = gamma2;
    double g4 = g2*g2;

    double numerator = ab*((1.+g4)*a4 - 4*g2*a2*b2 + (1.+g4)*b4);

    double b3 = b*b2;
    double b5 = b*b4;

    double num_a= (5.*(1.+g4)*a4*b - 12.*g2*a2*b3 + (1.+g4)*b5);

    double a3 = a*a2; 
    double a5 = a*a4;

    double num_b = ((1.+g4)*a5 - 12.*g2*a3*b2 + 5.*(1.+g4)*a*b4);

    double denominator = (g2*a2 - b2)*(a2 - g2*b2);
    double den_a = 2.*a*(-(1.+g4)*b2+2.*g2*a2);
    double den_b = 2.*b*(-(1.+g4)*a2+2.*g2*b2);


    double quotient = numerator/denominator;
    double quotient_a = (num_a*denominator - den_a*numerator)/( denominator*denominator );
    double quotient_b = (num_b*denominator - den_b*numerator)/( denominator*denominator );

    double F = 2.*a*b-(b2+a2)*sin(a)*sin(b) + quotient*cos(a)*cos(b);


    *F_a = 2.*b -2.*a*sin(a)*sin(b) -(b2+a2)*cos(a)*sin(b) -quotient*sin(a)*cos(b)+quotient_a*cos(a)*cos(b);
    *F_b = 2.*a -2.*b*sin(a)*sin(b)-quotient*cos(a)*sin(b)-  (b2+a2)*sin(a)*cos(b)+quotient_b*cos(a)*cos(b);
    return F;
}

//timoshenko beam,  a<ac, implicit relationship between the wave numbers
// page 951, a_c  equation (79), g2 : equation (81),
// super critical
//timoshenko beam,  implicit relationship between the wave numbers

double frequency_equation_cf(double a, double b, double gamma2, bool is_sub_critical, double* F_a, double* F_b)
{
    double F = 0;
    if( is_sub_critical )
    {
        F = frequency_equation_cf_sub(a, b, gamma2, F_a, F_b);
    }
    else
    {
        F = frequency_equation_cf_super(a, b, gamma2, F_a, F_b);
    }
    return F;
}
