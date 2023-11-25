#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include "wave_number.h"
#include "frequency_equation.h"
#include "secular.h"

void timoshenko_wave_numbers(
    double kmax, double gamma2, size_t mode, BoundaryCondition bc, double* aptr, double* bptr) {
  double arclength = .01;
  double factor = sqrt((gamma2 + 1.) / gamma2);
  // double jacobian = 0.;
  double ao = euler_wave_number(mode, bc);
  continuation_state current(ao, ao, 0.0, 0.0);
  bool subcritical = true;
  bool verbose = false;
  continuation_state previous = current;
  continuation_state before_last = previous;
  while (current.k < kmax) {
    before_last = previous;
    previous = current;
    current =
        arclengthcontinuation(previous, subcritical, arclength, gamma2, kmax, bc, before_last);

    bool next_subcritical = (factor - current.a * current.k > 0.);
    if (current.b == 0.0) {
      double a_critical = get_critical_point2(current.a, gamma2, bc);
      std::cout << "a_critical= " << a_critical << "\n";
      current.a = a_critical;
      next_subcritical = false;
    }
    if (subcritical != next_subcritical) {  // b=0 singularity
      double a_critical = get_critical_point2(current.a, gamma2, bc);
      current.a = a_critical;
      current.b = 0.;
      current.k = factor / a_critical;
      std::cout << "vanishing b singularity  k " << current.k << "  a " << current.a << "\n";
      subcritical = false;
    }
    if (current.b < 0.0) exit(-1);
    if (verbose) {
      std::cout << std::scientific << std::setprecision(6) << "k " << current.k << "  a "
                << current.a << "  b " << current.b << "\n";
    }
  }  // k
  *aptr = current.a;
  *bptr = current.b;
  wave_frequency_equation(aptr, bptr, gamma2, subcritical, bc, kmax);
}
