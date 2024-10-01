#include <iostream>
#include <iomanip>
#include <cmath>
#include "frequency_equation.h"
#include "timoshenko_wave_numbers.h"

int main() {
  double beam_thickness = 1.0;
  double beam_length = 10.0;
  double poissons_ratio = 0.3;
  double density = 0.3;
  double E = 3.0e+7;
  BoundaryCondition bc = freefree;
  size_t num_mode = 4;

  double nu = poissons_ratio;
  double kprime = 10.0* (1.0+ nu) / (12.0+ 11.0* nu);  // rectangular
  double h = beam_thickness;
  double area = h * h;  //  square cross section
  double I = h * h * h * h / 12.0;
  double kmax = sqrt(I / area) / beam_length;
  double gamma2 = 2.0* (1.0+ nu) / kprime;
  std::cout << "gamma2 = " << gamma2 << "\n";
  double factor = sqrt((gamma2 + 1.) / gamma2);

  for (size_t mode = 0; mode < num_mode; mode++) {
    double a = 0.0;
    double b = 0.0;
    timoshenko_wave_numbers(kmax, gamma2, mode, bc, &a, &b);
    bool subcritical = (factor - a * kmax > 0.0);
    double a2 = a * a;
    double b2 = b * b;
    double scaled_frequency =
        (subcritical ? sqrt((a2 - b2) / (1 + gamma2)) : sqrt((a2 + b2) / (1 + gamma2)));
    std::cout << std::scientific << std::setprecision(14);
    double frequency = scaled_frequency * sqrt(E / density) / (beam_length * 2.0* M_PI);
    std::cout << " mode " << mode << "   frequency " << frequency << "  a " << a << "  b " << b
              << "\n";
  }
  return 0;
}
