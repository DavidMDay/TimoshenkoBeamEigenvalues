#include <cmath>
#include "frequency_equation.h"
//timoshenko beam,  a<ac, implicit relationship between the wave numbers
// page 951, a_c  equation (79), g2 : equation (81),
// subcritical  a < ac

frequnency_equation_values clamped_free_sub(frequnency_equation_parameters p) {
  double a2 = p.a * p.a;
  double a3 = p.a * a2;
  double a4 = a2 * a2;
  double a5 = p.a * a4;
  double ab = p.a * p.b;
  double b2 = p.b * p.b;
  double b3 = p.b * b2;
  double b4 = b2 * b2;
  double b5 = p.b * b4;
  double g2 = p.gamma2;
  double g4 = g2 * g2;

  double numerator = ab * ((1.0 + g4) * a4 + 4 * g2 * a2 * b2 + (1.0 + g4) * b4);

  double num_a = 5. * (1.0 + g4) * a4 * p.b + 12.0 * g2 * a2 * b3 + (1.0 + g4) * b5;

  double num_b = (1.0 + g4) * a5 + 12.0 * g2 * a3 * b2 + 5. * (1.0 + g4) * p.a * b4;

  double denominator = (b2 + g2 * a2) * (a2 + g2 * b2);
  double den_a = 2.0 * p.a * ((1.0 + g4) * b2 + 2.0 * g2 * a2);
  double den_b = 2.0 * p.b * ((1.0 + g4) * a2 + 2.0 * g2 * b2);

  double quotient = numerator / denominator;
  double quotient_a = (num_a * denominator - den_a * numerator) / (denominator * denominator);
  double quotient_b = (num_b * denominator - den_b * numerator) / (denominator * denominator);
  frequnency_equation_values f;
  f.F = 2.0 * p.a * p.b + (b2 - a2) * sin(p.a) * sinh(p.b) + quotient * cos(p.a) * cosh(p.b);

  f.F_a = 2.0 * p.b - 2.0 * p.a * sin(p.a) * sinh(p.b) + (b2 - a2) * cos(p.a) * sinh(p.b) -
         quotient * sin(p.a) * cosh(p.b) + quotient_a * cos(p.a) * cosh(p.b);
  f.F_b = 2.0 * p.a + 2.0 * p.b * sin(p.a) * sinh(p.b) + quotient * cos(p.a) * sinh(p.b) +
         (b2 - a2) * sin(p.a) * cosh(p.b) + quotient_b * cos(p.a) * cosh(p.b);
  return f;
}

frequnency_equation_values clamped_free_super(frequnency_equation_parameters p) {
    double a2 = p.a * p.a;
    double a3 = p.a * a2;
    double a4 = a2 * a2;
    double a5 = p.a * a4;
    double ab = p.a * p.b;
    double b2 = p.b * p.b;
    double b3 = p.b * b2;
    double b4 = b2 * b2;
    double b5 = p.b * b4;
    double g2 = p.gamma2;
    double g4 = g2 * g2;

    double numerator = ab * ((1.0 + g4) * a4 - 4 * g2 * a2 * b2 + (1.0 + g4) * b4);

    double num_a = (5. * (1.0 + g4) * a4 * p.b - 12.0 * g2 * a2 * b3 + (1.0 + g4) * b5);

    double num_b = ((1.0 + g4) * a5 - 12.0 * g2 * a3 * b2 + 5. * (1.0 + g4) * p.a * b4);

    double denominator = (g2 * a2 - b2) * (a2 - g2 * b2);
    double den_a = 2.0 * p.a * (-(1.0 + g4) * b2 + 2.0 * g2 * a2);
    double den_b = 2.0 * p.b * (-(1.0 + g4) * a2 + 2.0 * g2 * b2);

    double quotient = numerator / denominator;
    double quotient_a = (num_a * denominator - den_a * numerator) / (denominator * denominator);
    double quotient_b = (num_b * denominator - den_b * numerator) / (denominator * denominator);
    frequnency_equation_values f;
    f.F = 2.0 * p.a * p.b - (b2 + a2) * sin(p.a) * sin(p.b) + quotient * cos(p.a) * cos(p.b);

    f.F_a = 2.0 * p.b - 2.0 * p.a * sin(p.a) * sin(p.b) - (b2 + a2) * cos(p.a) * sin(p.b) -
            quotient * sin(p.a) * cos(p.b) + quotient_a * cos(p.a) * cos(p.b);
    f.F_b = 2.0 * p.a - 2.0 * p.b * sin(p.a) * sin(p.b) - quotient * cos(p.a) * sin(p.b) -
            (b2 + a2) * sin(p.a) * cos(p.b) + quotient_b * cos(p.a) * cos(p.b);
    return f;
}

//timoshenko beam,  a<ac, implicit relationship between the wave numbers
// page 951, a_c  equation (79), g2 : equation (81),
// super critical
//timoshenko beam,  implicit relationship between the wave numbers
frequnency_equation_values clamped_free(frequnency_equation_parameters p) {
    if (p.is_sub_critical) {
      return clamped_free_sub(p);
    } else {
      return clamped_free_super(p);
    }
    return {};
}
