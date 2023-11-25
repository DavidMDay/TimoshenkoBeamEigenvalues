#include <cassert>
#include <cmath>
#include "wave_number.h"
#include "frequency_equation.h"

// output Jacobian,  new (a,b) pair
// input  k + frequnency_equation_parameters 
double wave_frequency_equation(double* aptr,
                              double* bptr,
                              double gamma2,
                              bool subcritical,
                              BoundaryCondition bc,
                              double k) {
  frequnency_equation_parameters state;
  state.a = *aptr; // (a,b): F=0 & f=0
  state.b = *bptr;
  state.gamma2 = gamma2;
  state.is_sub_critical = subcritical;
  state.bc = bc;
  double jacobian = 0.0;
  assert(state.b > 0.0);
  unsigned num_iter = 10;
  for (unsigned iteration = 0; iteration < num_iter; iteration++) {
    auto value = frequency_equation(state);
    double dfda = wave_equation_dfda(k,state.gamma2,state.a,state.b,state.is_sub_critical);
    double dfdb = wave_equation_dfdb(k,state.gamma2,state.a,state.b,state.is_sub_critical);
    //J = [Fa Fb;...   rhs = -[F;f];
    //     fa fb];     J [da;db] = rhs
    jacobian = value.F_a*dfdb-value.F_b*dfda;
    double residual = wave_number_residual(k,state.gamma2,state.a,state.b,state.is_sub_critical); // f
    double da = (value.F_b*residual - dfdb*value.F)/jacobian;
    double db = (dfda*value.F-value.F_a*residual)/jacobian;
    state.a += da;
    state.b += db;
    double error = sqrt(value.F * value.F  + residual * residual );
    if (iteration == num_iter-1 ) {
      assert(error < .001 );
    }
    if (error < 1.e-8 ) break;
 }
 *aptr = state.a;
 *bptr = state.b;
 return jacobian;
}

