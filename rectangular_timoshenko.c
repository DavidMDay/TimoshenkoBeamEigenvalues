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
  double beam_thickness = 1.;
  double beam_length = 10.;
  double poissons_ratio = .3;
  double density= .3;
  double E = 3.e+7;
  bool is_ff = true;
  unsigned num_mode = 4;

  double nu = poissons_ratio;
  double kprime = 10. * (1. + nu)/(12.+11.*nu); // rectangular
  double h = beam_thickness;
  double area = h*h; //  square cross section
  double I = h*h*h*h/12.;
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
      double frequency = scaled_frequency*sqrt(E/density)/(beam_length*2.*M_PI);
      std::cout << mode<< "    "<<frequency<<"   "<< a << " " << b << "\n";
  }
  return 0;
}
