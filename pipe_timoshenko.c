#include <iostream>
#include <iomanip>
#include <math.h>  // sqrt
void timoshenko_wave_numbers( double kmax,
                              double gamma2,
                              unsigned mode,
                              bool is_ff,
                              double* aptr,
                              double* bptr);
int main()
{
  double beam_length = 1.;
  double poissons_ratio = .29;
  double density= 7830.;
  double E = 2.e+11;
  bool is_ff = true;
  unsigned num_mode = 4;
  double outer_radius = .16;
  double inner_radius = .15;

  double nu = poissons_ratio;
  double mm = inner_radius/outer_radius;
  double mm2p12 = (1+mm*mm)*(1+mm*mm);
  double kprime = 6*(1+nu)*mm2p12/(  (7+6*nu)*mm2p12 + (20+12*nu)*mm*mm  );

std::cout << " kprime = " << kprime << "\n";
  double area = M_PI*(outer_radius*outer_radius - inner_radius*inner_radius );
// area moment of inertia of the cross section about the neutral axis

  double ir = inner_radius;
  double rad = outer_radius;
  double I = (rad*rad*rad*rad - ir*ir*ir*ir)*.25*M_PI;
  double kmax = sqrt(I/area)/beam_length;
  double gamma2 = 2.*(1.+nu)/kprime;
  std::cout << "gamma2 = " << gamma2 << "\n";
  double factor = sqrt( (gamma2 + 1.)/gamma2 );

  for( unsigned mode = 0; mode<num_mode; mode++)
  {
      double a = 0;
      double b = 0;
      timoshenko_wave_numbers( kmax, gamma2, mode, is_ff, &a, &b);
      bool subcritical = ( factor-a*kmax > 0. );
      double scaled_frequency = 0.;
      if( subcritical )
               scaled_frequency = sqrt(( a*a - b*b )/(1+gamma2));
      else
               scaled_frequency = sqrt(( a*a + b*b )/(1+gamma2));

      std::cout << std::scientific << std::setprecision(14);
      //double frequency = scaled_frequency*sqrt(E/density)/(beam_length*2.*M_PI);
      std::cout << mode<< "    "<<scaled_frequency<<"   "<< a << " " << b << "\n";
  }
  return 0;
}
