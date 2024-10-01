#ifndef TIMOSHENKO_BEAM_EIGENVALUES_FREQUENCY_EQUATION
#define TIMOSHENKO_BEAM_EIGENVALUES_FREQUENCY_EQUATION

#include <cstddef>
enum BoundaryCondition { freefree, clampedclamped, clampedfree };

struct continuation_state {
  continuation_state( double alpha, double beta, double kappa, double det) :
  a(alpha),
  b(beta),
  k(kappa),
  jacobian(det)
  {}
  double a = 0.0;
  double b = 0.0;
  double k = 0.0;
  double jacobian = 0.0;
};

struct frequnency_equation_parameters {
  double a = 0.0;
  double b = 0.0;       // end points?
  double gamma2 = 0.0;  // physical
  bool is_sub_critical = false;
  BoundaryCondition bc = freefree;
};

struct frequnency_equation_values {
  double F = 0.0;
  double F_a = 0.0;
  double F_b = 0.0;
};

frequnency_equation_values clamped_free(frequnency_equation_parameters p);
frequnency_equation_values clamped_clamped(frequnency_equation_parameters p);
frequnency_equation_values frequency_equation(frequnency_equation_parameters p);
double euler_wave_number(size_t id, enum BoundaryCondition bc);

double wave_frequency_equation(double* aptr,
                              double* bptr,
                              double gamma2,
                              bool subcritical,
                              BoundaryCondition bc,
                              double k);

continuation_state
arclengthcontinuation(continuation_state current,
                              bool subcritical,
                              double arclength,
                              double gamma2,
                              double kmax,
                              BoundaryCondition bc,
                              continuation_state previous);

#endif  // TIMOSHENKO_BEAM_EIGENVALUES_FREQUENCY_EQUATION

