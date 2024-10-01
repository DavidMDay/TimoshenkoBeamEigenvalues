#ifndef TIMOSHENKO_BEAM_EIGENVALUES_SECULAR_EQUATION
#define TIMOSHENKO_BEAM_EIGENVALUES_SECULAR_EQUATION

#include <cstddef>
#include "frequency_equation.h"
double get_critical_point(size_t mode, double gamma2, enum BoundaryCondition bc);
double get_critical_point2(double a, double gamma2, enum BoundaryCondition bc);

#endif
