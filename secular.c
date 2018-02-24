#include <math.h>
#include <iostream>
#include <assert.h>

double secular_equation_ff(double a, double gamma2)
{
      double sec  = .5*sin(a)*a/gamma2 - cos(a)  + 1;
      return sec;
}

double secular_equation_derivative_ff(double a, double gamma2)
{
      double dsecda = .5*cos(a)*a/gamma2 + sin(a)*(.5/gamma2 + 1.);
      return dsecda;
}

double secular_equation_cc(double a, double gamma2)
{
      double sec  = .5*sin(a)*a*gamma2 - cos(a)  + 1;
      return sec;
}

double secular_equation_derivative_cc(double a, double gamma2)
{
      double dsecda = .5*cos(a)*a*gamma2 + sin(a)*(.5*gamma2 + 1.);
      return dsecda;
}

double secular_equation(double a, double gamma2, bool is_ff)
{
    double sec  = 0;
    if( is_ff )
      sec  = secular_equation_ff(a, gamma2);
    else
      sec  = secular_equation_cc(a, gamma2);
}

double secular_equation_derivative(double a, double gamma2, bool is_ff)
{
    double dsecda  = 0;
    if( is_ff )
      dsecda  = secular_equation_derivative_ff(a, gamma2);
    else
      dsecda  = secular_equation_derivative_cc(a, gamma2);
}

double get_approximate_critical_point( unsigned mode, double gamma2, bool is_ff )
{
    double pi = 4.*atan(1.);
    if( is_ff && mode%2 == 1 )
    {
        double odd = static_cast<double>(mode);
        if( mode > 1 )
        {
            double estimate = odd*pi  + 4.*gamma2/(pi*odd);
            return estimate;
        }
        else
        {
            double left = 2.; // secular_equation_ff(pi, gamma2);
            double scale = 1.5;
            double right = secular_equation_ff(scale*pi, gamma2);
            while (right > 0.)
            {
               scale = .5*(2.+scale);
               right = secular_equation_ff(scale*pi, gamma2);
            }
            assert( left * right  < 0);
            double s = left/( left-right);
            double estimate  = pi*(1. + .5*s);
            return estimate;
        }
    }
    else
    {
        double even = static_cast<double>(mode);
        double exact = even*pi;
        return exact;
    }
}

double get_critical_point( unsigned mode, double gamma2, bool is_ff )
{
    double a = get_approximate_critical_point(mode,gamma2,is_ff);
    double p = secular_equation(a,gamma2,is_ff);
    bool converged = false;
    unsigned max_step = 10;
    unsigned step = 0; 
    while( step < max_step && !converged )
    for( unsigned step = 0; step < 10; step++)
    {
        double dpda = secular_equation_derivative(a, gamma2,is_ff);
        double d = -p/dpda ;
        a = a + d;
        p = secular_equation(a,gamma2,is_ff);
        if( fabs(d) < 1.e-14 ) converged = true;
        if( fabs(p) < 1.e-14 ) converged = true;
    }
    if(!converged) std::cout << " get critical point did not converge\n";
    return a;
}
