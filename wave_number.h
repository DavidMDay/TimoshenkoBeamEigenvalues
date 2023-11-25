
#ifndef TIMOSHENKO_BEAM_EIGENVALUES_wave_EQUATION
#define TIMOSHENKO_BEAM_EIGENVALUES_wave_EQUATION

//#include <cstddef>
//#include "frequency_equation.h"

double wave_number_residual(double k,double gamma2,double a,double b,bool is_sub_critical);
double wave_equation_dfdk(double k,double gamma2,double a,double b,bool is_sub_critical);
double wave_equation_dfdb(double k,double gamma2,double a,double b,bool is_sub_critical);
double wave_equation_dfda(double k,double gamma2,double a,double b,bool is_sub_critical);
double wave_equation_dbdk(double dfdk, double dfdb);
double wave_equation_dbda(double dfda, double dfdb);

#endif
