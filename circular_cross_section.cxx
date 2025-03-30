#include <iostream>
#include <iomanip>
#include <cmath>
#include "src/frequency_equation.hpp"
#include "src/timoshenko_wave_numbers.hpp"

struct cylinder {
  double outer_radius;
  double inner_radius;
  double beam_length;
};

struct beam_parameters {
  double mm;
  double mm2p12;
  double kmax;
};

beam_parameters cylindrical_beam_parameters( cylinder geo ) {
  double mm = geo.inner_radius / geo.outer_radius;
  double mm2p12 = (1 + mm * mm) * (1 + mm * mm);
  double area = M_PI * (geo.outer_radius * geo.outer_radius - geo.inner_radius * geo.inner_radius);
  std::cout << "cross section area " << area << "\n";

  double ir = geo.inner_radius;
  double rad = geo.outer_radius;
  double I = (rad * rad * rad * rad - ir * ir * ir * ir) * .25 * M_PI;
  std::cout << " I " << I << "\n";
  double kmax = sqrt(I / area) / geo.beam_length;
  std::cout << "kmax = " << kmax << "\n";
  return {mm,mm2p12,kmax};
}

struct material {
  double E;
  double density;
  double poissons_ratio;
};

struct mode {
  unsigned num_mode;
  BoundaryCondition bc;
};



int main() {
  cylinder geo = {.16, .15, 1.};
  auto mk = cylindrical_beam_parameters(geo );

  material al = {2.e+11, 7830., .29};
  double nu = al.poissons_ratio;
  double kprime = 6. * (1. + nu) * mk.mm2p12 / ((7. + 6. * nu) * mk.mm2p12 + (20. + 12. * nu) * mk.mm * mk.mm);
  double gamma2 = 2. * (1. + nu) / kprime;
  std::cout << "gamma = " << sqrt(gamma2) << "\n";
  double factor = sqrt((gamma2 + 1.) / gamma2);

  mode md = {12, clampedfree};
  std::cout << std::scientific << std::setprecision(14);
  for (unsigned mode = 0; mode < md.num_mode; mode++) {
    double a = 0;
    double b = 0;
    timoshenko_wave_numbers(mk.kmax, gamma2, mode, md.bc, &a, &b);
    bool subcritical = (factor - a * mk.kmax > 0.);

    double a2 = a * a;
    double b2 = b * b;
    double scaled_frequency =
        (subcritical ? sqrt((a2 - b2) / (1 + gamma2)) : sqrt((a2 + b2) / (1 + gamma2)));

    double frequency = scaled_frequency * sqrt(al.E / al.density) / (geo.beam_length);
    std::cout << "mode " << mode << "  f  " << frequency << "\n";
  }
  return 0;
}
