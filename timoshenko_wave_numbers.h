
#ifndef TIMOSHENKO_BEAM_EIGENVALUES_wave_numbers
#define TIMOSHENKO_BEAM_EIGENVALUES_wave_numbers

#include <cstddef>
#include "frequency_equation.h"


void timoshenko_wave_numbers( double kmax,
                              double gamma2,
                              size_t mode,
                              BoundaryCondition bc,
                              double* aptr,
                              double* bptr);

#endif
