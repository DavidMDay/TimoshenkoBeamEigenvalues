#include <cmath>
#include "frequency_equation.hpp"
//timoshenko beam,  a<ac, implicit relationship between the wave numbers
// page 951, a_c  equation (79), g2 : equation (81),
// subcritical  a < ac

frequnency_equation_values clamped_clamped_sub(frequnency_equation_parameters p) {
  double a2 = p.a * p.a;
  double b2 = p.b * p.b;
  double ab = p.a * p.b;
  double g2 = p.gamma2;
  double numerator =
      (a2 - b2) * (g2 * a2 + g2 * b2 + (g2 - 1.0) * ab) * (g2 * a2 + g2 * b2 - (g2 - 1.0) * ab);
  double num_a =
      2.0 * p.a * (g2 * a2 + g2 * b2 + (g2 - 1.0) * ab) * (g2 * a2 + g2 * b2 - (g2 - 1.0) * ab) +
      (a2 - b2) * (2.0 * g2 * p.a + (g2 - 1.0) * p.b) * (g2 * a2 + g2 * b2 - (g2 - 1.0) * ab) +
      (a2 - b2) * (g2 * a2 + g2 * b2 + (g2 - 1.0) * ab) * (2.0 * g2 * p.a - (g2 - 1.0) * p.b);

  double num_b =
      -2.0 * p.b * (g2 * a2 + g2 * b2 + (g2 - 1.0) * ab) * (g2 * a2 + g2 * b2 - (g2 - 1.0) * ab) +
      (a2 - b2) * (2.0 * g2 * p.b + (g2 - 1.0) * p.a) * (g2 * a2 + g2 * b2 - (g2 - 1.0) * ab) +
      (a2 - b2) * (g2 * a2 + g2 * b2 + (g2 - 1.0) * ab) * (2.0 * g2 * p.b - (g2 - 1.0) * p.a);

  double denominator = 2.0 * ab * (b2 + g2 * a2) * (a2 + g2 * b2);

  double den_a = 2.0 * p.b *
                 ((b2 + g2 * a2) * (a2 + g2 * b2) + 2.0 * a2 * g2 * (a2 + g2 * b2) +
                  2.0 * a2 * (b2 + g2 * a2));

  double den_b = 2.0 * p.a *
                 ((b2 + g2 * a2) * (a2 + g2 * b2) + 2.0 * b2 * (a2 + g2 * b2) +
                  2.0 * g2 * b2 * (b2 + g2 * a2));

  double quotient = numerator / denominator;
  double quotient_a = (num_a * denominator - den_a * numerator) / (denominator * denominator);
  double quotient_b = (num_b * denominator - den_b * numerator) / (denominator * denominator);

  frequnency_equation_values f;
  f.F = quotient * sin(p.a) * sinh(p.b) - cos(p.a) * cosh(p.b) + 1.0;
  f.F_a = quotient_a * sin(p.a) * sinh(p.b) + quotient * cos(p.a) * sinh(p.b) + sin(p.a) * cosh(p.b);
  f.F_b = quotient_b * sin(p.a) * sinh(p.b) + quotient * sin(p.a) * cosh(p.b) - cos(p.a) * sinh(p.b);
  return f;
}
// timoshenko beam,  a<ac, implicit relationship between the wave numbers
//  page 951, a_c  equation (79), g2 : equation (81),
//  super critical

frequnency_equation_values clamped_clamped_super(frequnency_equation_parameters p) {
  double a2 = p.a * p.a;
  double b2 = p.b * p.b;
  double ab = p.a * p.b;
  double g2 = p.gamma2;
  double g4 = g2 * g2;

  double numerator = (a2 + b2) * (g4 * (a2 - b2) * (a2 - b2) + (g2 - 1.0) * (g2 - 1.0) * a2 * b2);

  double num_a =
      (2.0 * p.a) * (g4 * (a2 - b2) * (a2 - b2) + (g2 - 1.0) * (g2 - 1.0) * a2 * b2) +
      (a2 + b2) * (2.0 * g4 * (a2 - b2) * 2.0 * p.a + 2.0 * (g2 - 1.0) * (g2 - 1.0) * p.a * b2);

  double num_b =
      (2.0 * p.b) * (g4 * (a2 - b2) * (a2 - b2) + (g2 - 1.0) * (g2 - 1.0) * a2 * b2) +
      (a2 + b2) * (-2.0 * g4 * (a2 - b2) * 2.0 * p.b + 2.0 * (g2 - 1.0) * (g2 - 1.0) * a2 * p.b);

  double denominator = 2.0 * ab * (g2 * a2 - b2) * (a2 - g2 * b2);

  double den_a = 2.0 * p.b *
                 ((g2 * a2 - b2) * (a2 - g2 * b2) + (2.0 * g2 * a2) * (a2 - g2 * b2) +
                  (g2 * a2 - b2) * (2.0 * a2));

  double den_b = 2.0 * p.a *
                 ((g2 * a2 - b2) * (a2 - g2 * b2) - 2.0 * b2 * (a2 - g2 * b2) -
                  2.0 * (g2 * a2 - b2) * g2 * b2);

  double quotient = numerator / denominator;
  double quotient_a = (num_a * denominator - den_a * numerator) / (denominator * denominator);
  double quotient_b = (num_b * denominator - den_b * numerator) / (denominator * denominator);

  frequnency_equation_values f;
  f.F = quotient * sin(p.a) * sin(p.b) - cos(p.a) * cos(p.b) + 1.0;
  f.F_a = quotient_a * sin(p.a) * sin(p.b) + quotient * cos(p.a) * sin(p.b) + sin(p.a) * cos(p.b);
  f.F_b = quotient_b * sin(p.a) * sin(p.b) + quotient * sin(p.a) * cos(p.b) + cos(p.a) * sin(p.b);
  return f;
}
// timoshenko beam,  implicit relationship between the wave numbers
frequnency_equation_values clamped_clamped(frequnency_equation_parameters p) {
  if (p.is_sub_critical) {
    return clamped_clamped_sub(p);
  } else {
    return clamped_clamped_super(p);
  }
  return {};
}
