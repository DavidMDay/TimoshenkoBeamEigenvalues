#include <cmath>
#include <iostream>
#include <limits>  // std::numeric_limits
#include <cstdlib>
#include "frequency_equation.h"
//timoshenko beam,  a<ac, implicit relationship between the wave numbers
// page 951, a_c  equation (79), g2 : equation (81),
// subcritical  a < ac

frequnency_equation_values frequency_equation_sub(frequnency_equation_parameters p) {
    double a2 = p.a*p.a;
    double b2 = p.b*p.b;
    double ab = p.a*p.b;
    double g2 = p.gamma2;

    double numerator = (a2 - b2)*(a2 + b2 + (g2-1.0)*ab)*(a2 + b2 - (g2-1.0)*ab);
    double num_a=
    2.0*p.a*(a2 + b2 + (g2-1.0)*ab)*(a2 + b2 - (g2-1.0)*ab)+
    (a2 - b2)*(2*p.a + (g2-1.0)*p.b)*(a2 + b2 - (g2-1.0)*ab)+
    (a2 - b2)*(a2 + b2 + (g2-1.0)*ab)*(2.0*p.a - (g2-1.0)*p.b);

    double num_b =
    -2.0*p.b*(a2 + b2 + (g2-1.0)*ab)*(a2 + b2 - (g2-1.0)*ab)+
    (a2 - b2)*(2*p.b + (g2-1.0)*p.a)*(a2 + b2 - (g2-1.0)*ab)+
    (a2 - b2)*(a2 + b2 + (g2-1.0)*ab)*(2.0*p.b - (g2-1.0)*p.a);

    double denominator = 2.0*ab*(g2*a2+b2)*(a2 + g2*b2);  // 2 g2 a5 /b

    double den_a = 2.0*p.b*(
    (b2 + g2* a2)*(a2 + g2*b2)+
    2.0*a2*g2*(a2 + g2*b2)+
    2.0*a2*(b2 + g2*a2));

    double den_b = 2.0*p.a*(
     (g2*a2+b2)*(a2 + g2*b2)+
     2.0*b2*(a2 + g2*b2)+
     2.0*g2*b2*(b2 + g2*a2));   // 2 g2 a5

    double quotient = numerator/denominator;

    double quotient_a = (num_a*denominator - den_a*numerator)/(denominator*denominator);
    double quotient_b = (num_b*denominator - den_b*numerator)/(denominator*denominator);

    frequnency_equation_values f;

    f.F =      quotient*sin(p.a)*sinh(p.b) - cos(p.a)*cosh(p.b)+1.0;
    f.F_a =  quotient_a*sin(p.a)*sinh(p.b) + quotient*cos(p.a)*sinh(p.b) + sin(p.a)*cosh(p.b);
    f.F_b =  quotient_b*sin(p.a)*sinh(p.b) + quotient*sin(p.a)*cosh(p.b) - cos(p.a)*sinh(p.b);

    return f;
}
//timoshenko beam,  a<ac, implicit relationship between the wave numbers
// page 951, a_c  equation (79), g2 : equation (81),
// super critical


frequnency_equation_values frequency_equation_super(frequnency_equation_parameters p) {
    double a2 = p.a*p.a;
    double b2 = p.b*p.b;
    double ab = p.a*p.b;
    double g2 = p.gamma2;

    double numerator = (a2+b2)*( (a2-b2)*(a2-b2)  + (g2-1.)*(g2-1.)*a2*b2 );

    double num_a =
        (2.0 * p.a) * ((a2 - b2) * (a2 - b2) + (g2 - 1.) * (g2 - 1.) * a2 * b2) +
        (a2 + b2) * (2.0 * (a2 - b2) * 2.0 * p.a + 2.0 * (g2 - 1.) * (g2 - 1.) * p.a * b2);

    double num_b =
        (2.0 * p.b) * ((a2 - b2) * (a2 - b2) + (g2 - 1.) * (g2 - 1.) * a2 * b2) +
        (a2 + b2) * (-2.0 * (a2 - b2) * 2.0 * p.b + 2.0 * (g2 - 1.) * (g2 - 1.) * a2 * p.b);

    double denominator = 2.0*ab*(g2*a2-b2)*(a2 - g2*b2);

    double den_a = 2.0 * p.b *
                   ((g2 * a2 - b2) * (a2 - g2 * b2) + (2.0 * g2 * a2) * (a2 - g2 * b2) +
                    (g2 * a2 - b2) * (2.0 * a2));

    double den_b = 2.0 * p.a *
                   ((g2 * a2 - b2) * (a2 - g2 * b2) - 2.0 * b2 * (a2 - g2 * b2) -
                    2.0 * (g2 * a2 - b2) * g2 * b2);

    double quotient = numerator/denominator;
    double quotient_a = (num_a*denominator - den_a*numerator)/( denominator*denominator );
    double quotient_b = (num_b*denominator - den_b*numerator)/( denominator*denominator );
    frequnency_equation_values f;

    f.F = quotient * sin(p.a) * sin(p.b) - cos(p.a) * cos(p.b) + 1.;
    f.F_a = quotient_a * sin(p.a) * sin(p.b) + quotient * cos(p.a) * sin(p.b) + sin(p.a) * cos(p.b);
    f.F_b = quotient_b * sin(p.a) * sin(p.b) + quotient * sin(p.a) * cos(p.b) + cos(p.a) * sin(p.b);
    return f;
}
//timoshenko beam,  implicit relationship between the wave numbers
frequnency_equation_values frequency_equation(frequnency_equation_parameters p) {
    if (p.bc == freefree) {
      if (p.is_sub_critical) {
        return frequency_equation_sub(p);
      } else {
        return frequency_equation_super(p);
      }
    } else {
      if (p.bc == clampedclamped) {
        return clamped_clamped(p);
      } else {
        return clamped_free(p);
      }
    }
    return {};
}

double euler_wave_number( size_t id, BoundaryCondition bc ) {
    double eulerfreefree[] = {4.730040744862704, 7.853204624095838,
     10.99560783800169, 14.13716549125746, 17.27875965739948};
    double eulerclampedfree[] = {1.87510406871196, 4.694091132974175,
     7.854757438237613, 10.99554073487547, 14.13716839104647};

    if( id < 5 ) {
      return ( bc == clampedfree ?  eulerclampedfree[id] : eulerfreefree[id] );
    }

    double offset = ( bc == clampedfree ? .5 : 1.5);

    double lambda = (offset +  static_cast<double>(id))*M_PI;

    auto max_id = static_cast<size_t>( log( std::numeric_limits<double>::max() )/M_PI  - 1.5);
    if( id >= max_id ) {
        return lambda;
    }

    double r = 1.;
    double dr =1.;
    double h = 1.;
    size_t max_iter = id;
    size_t iteration = 0;

    double sign = ( bc == clampedfree ? 1. : -1.);

    while( (iteration < max_iter) && (fabs(r) > fabs(dr)*1.e-13) && (fabs(h) > 1.e-14) ) {
        r = cos(lambda)*cosh(lambda) + sign;
        dr = cos(lambda)*sinh(lambda) - sin(lambda)*cosh(lambda);
        h = -r/dr;
        lambda += h;
        iteration = iteration + 1;
    }
    if( iteration == max_iter ) {
        std::cerr << " fatal error: unable to accurately determine Euler wave number\n";
        exit(-1);
    }
    return lambda;
}
