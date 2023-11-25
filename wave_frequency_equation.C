#include <cassert>
#include <cmath>
#include "wave_number.h"
#include "frequency_equation.h"

// output Jacobian, new (a,b) pair
// input  k and frequnency_equation_parameters 
double wave_frequency_equation(double* aptr,
                              double* bptr,
                              double gamma2,
                              bool subcritical,
                              BoundaryCondition bc,
                              double k) {
  frequnency_equation_parameters parameter;
  parameter.a = *aptr; // (a,b): F=0 & f=0
  parameter.b = *bptr;
  parameter.gamma2 = gamma2;
  parameter.is_sub_critical = subcritical;
  parameter.bc = bc;
  double jacobian = 0.0;
  assert(parameter.b > 0.0);
  unsigned num_iter = 10;
  for (unsigned iteration = 0; iteration < num_iter; iteration++) {
    auto value = frequency_equation(parameter);
    double dfda = wave_equation_dfda(k, parameter.gamma2, parameter.a, parameter.b,
                                     parameter.is_sub_critical);
    double dfdb = wave_equation_dfdb(k, parameter.gamma2, parameter.a, parameter.b,
                                     parameter.is_sub_critical);
    //J = [Fa Fb;...   rhs = -[F;f];
    //     fa fb];     J [da;db] = rhs
    jacobian = value.F_a*dfdb-value.F_b*dfda;
    double residual = wave_number_residual(k, parameter.gamma2, parameter.a, parameter.b,
                                           parameter.is_sub_critical);  // f
    double da = (value.F_b*residual - dfdb*value.F)/jacobian;
    double db = (dfda*value.F-value.F_a*residual)/jacobian;
    parameter.a += da;
    parameter.b += db;
    double error = sqrt(value.F * value.F + residual * residual);
    if (iteration == num_iter-1 ) {
      assert(error < .001 );
    }
    if (error < 1.e-8 ) break;
 }
 *aptr = parameter.a;
 *bptr = parameter.b;
 return jacobian;
}

