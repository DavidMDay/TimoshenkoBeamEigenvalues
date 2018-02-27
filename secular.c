#include <math.h>
#include <iostream>
#include <assert.h>
#include "frequency_equation.h"

double secular_equation_ff(double a, double gamma2)
{
      double sec  = .5*sin(a)*a/gamma2 - cos(a)  + 1;
      return sec;
// a = k pi + t,   .5 cos(k pi) t k pi/g2 - cos(kpi) + 1 = 0
// t = 2 g2 (1- cos(k pi))/(k pi)
//  5.2, 6.28, 10.88,12.56,
// odd k:   1.7 pi, 3.5 pi, 5.3 pi, 7.2 pi, 9.2 pi,
// even k     2 pi    4 pi    6 pi       8 pi      10 pi
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
// a = k pi + t    .5 cos(k pi) sin(t) (k pi) g2 - cos(k pi) + 1 = 0
//    cos(k pi) - 1 = .5 cos(k pi) sin(t) (k pi) g2 
//    1 - cos(k pi) = .5 t (k pi) g2 
//    t = 2 (1- cos(k pi))/(g2 k pi)
// odd k:   1.1 pi, 3.0 pi, 5.0 pi, 7.0 pi, 9.0 pi,
// even k     2 pi    4 pi    6 pi       8 pi      10 pi
}

double secular_equation_derivative_cc(double a, double gamma2)
{
      double dsecda = .5*cos(a)*a*gamma2 + sin(a)*(.5*gamma2 + 1.);
      return dsecda;
}


double secular_equation_cf(double a, double gamma2)
{
      double g4   = gamma2*gamma2;
      double sec  = 2. -a*sin(a) +  cos(a)*(1.+g4)/gamma2;
      return sec;
// a = pi/2, pi, 2 pi , 3 pi , 4 pi , 5 pi
// a-> infty,   sin(a) -> 0   a = k pi + t
//  sin( k pi +t) =  cos(k pi) sin(t)
//  2 - k pi cos(k pi) t  + cos(k pi) (1+g4)/g2
//   k pi cos(k pi) t = 2 + cos(k pi) (1+g4)/g2
//   t = (2 cos(k pi) + (1+g4)/g2)/(k pi)
}

double secular_equation_derivative_cf(double a, double gamma2)
{
      double g4 = gamma2*gamma2;
      double dsecda = -a*cos(a) - sin(a)*(1.+gamma2+g4)/gamma2;
      return dsecda;
}


double secular_equation(double a, double gamma2, boundary_condition bc)
{
    double sec  = 0;
    if( bc == freefree )
      sec  = secular_equation_ff(a, gamma2);
    else
    {
      if( bc == clampedclamped )
          sec  = secular_equation_cc(a, gamma2);
      else
      {
          sec  = secular_equation_cf(a, gamma2);
      }
    }
}

double secular_equation_derivative(double a, double gamma2, boundary_condition bc )
{
    double dsecda  = 0;
    if( bc == freefree )
      dsecda  = secular_equation_derivative_ff(a, gamma2);
    else
      if( bc == clampedclamped )
          dsecda  = secular_equation_derivative_cc(a, gamma2);
      else
      {
          dsecda  = secular_equation_derivative_cf(a, gamma2);
      }
}

double get_approximate_critical_point( unsigned mode, double gamma2, boundary_condition bc )
{
    double pi = 4.*atan(1.);
    if( bc == freefree && mode%2 == 1 )
    {
        double odd = static_cast<double>(mode);
        if( mode > 1 )
        {
            double estimate = odd*pi  + 4.*gamma2/(pi*odd);
            return estimate;
        }
        else  // mode == 1
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
        double even = static_cast<double>(mode); // cc odd or even,
        double exact = even*pi;
        return exact;
    }
}

double get_critical_point( unsigned mode, double gamma2, boundary_condition bc)
{
    double a = get_approximate_critical_point(mode,gamma2,bc);
    double p = secular_equation(a,gamma2,bc);
    bool converged = false;
    unsigned max_step = 10;
    unsigned step = 0; 
    while( step < max_step && !converged )
    for( unsigned step = 0; step < 10; step++)
    {
        double dpda = secular_equation_derivative(a, gamma2,bc);
        double d = -p/dpda ;
        a = a + d;
        p = secular_equation(a,gamma2,bc);
        if( fabs(d) < 1.e-14 ) converged = true;
        if( fabs(p) < 1.e-14 ) converged = true;
    }
    if(!converged) std::cout << " get critical point did not converge\n";
    return a;
}


double get_critical_point2( double a, double gamma2, boundary_condition bc)
{

    //std::cout << "  get_critical_point2  a = " << a << "\n";
    double p = secular_equation(a,gamma2,bc);
    bool converged = false;
    unsigned max_step = 10;
    unsigned step = 0;
    while( step < max_step && !converged )
    {
        double dpda = secular_equation_derivative(a, gamma2,bc);
        double d = -p/dpda ;
        a = a + d;
        p = secular_equation(a,gamma2,bc);
        //std::cout << "  get_critical_point2  a = " << a << "  "  << p << "\n";
        if( fabs(d) < 1.e-14 ) converged = true;
        if( fabs(p) < 1.e-14 ) converged = true;
    }
    if(!converged) std::cout << " get critical point did not converge\n";
    return a;
}


