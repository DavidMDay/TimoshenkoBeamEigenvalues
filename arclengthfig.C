// The original purpose was to plot the number of iterations
// taken by the nonlinear solver per step of the arc length
// continuation algorithm.  Files 'param.data' and 'iter.data' were created.  
// Currently it seems to plot a(k), b(k), k and Jacobian(k)
// as k increases from 0 to kmax.  a(0)=b(0) = 17.2
// a(kmax) = 12.5,  b(kmax) = 5.5
// At kx = 0.065,  b(kx)=0, jacobian(kx)=0,  a(kx) = 16.7
// The Jacobian changes sign often without b changing sign.
#include <iostream>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include "wave_number.h"
#include "frequency_equation.h"
#include "secular.h"

// mode 4
int main() {
  double kmax = .4;// .2273; // .4;// physics
  double gamma = 2.205;
  BoundaryCondition bc = freefree;
  double arclength = .01;
  double gamma2 = gamma*gamma;
  double factor = sqrt( (gamma2 + 1.0)/gamma2);
  double jacobian = 0.0;
  for(unsigned mode = 4; mode < 5; mode++) {
    double a = euler_wave_number(mode, bc);
    double b = a;
    double a_critical = get_critical_point(mode+1, gamma2, bc);
    double k = 0.0;
    bool subcritical = true;
    bool verbose = false;
    continuation_state current(a,b,k,jacobian);
    continuation_state previous = current;
    continuation_state before_last = previous;

    int max_iter_per_mode = 100000;
    int iter = 0;
    std::cout << "    a      b       k     jacobian" << std::endl;
    while( current.k < kmax && iter < max_iter_per_mode) {
        before_last = previous;
        previous = current;
        current = arclengthcontinuation(previous,subcritical,arclength,
                                                gamma2,kmax,bc,before_last);
        bool next_subcritical = ( factor-current.a* current.k > 0.0);
        if ( current.b == 0.0) {
            if (verbose) { 
              std::cout<<"a=a_critical= "<<a_critical<<"\n";
            }
            current.a=a_critical;
            next_subcritical = false;
        }
        if( subcritical != next_subcritical) { // b=0 singularity
            current.a = a_critical;
            current.b = 0.;
            current.k = factor/a_critical;
            subcritical = false;
            //verbose = false;
        }
        //    double a2 = current.a * current.a, b2 = current.b * current.b;
        //    double scaled_frequency = (subcritical ? sqrt((a2 - b2)/(1 + gamma2))
        //                                           : sqrt((a2 + b2)/(1 + gamma2)));
        if ( current.k < kmax) {
            std::cout <<current.a<<"  "<<current.b<<"  " <<current.k<<"  "<<current.jacobian<< std::endl;
        }
        if ( current.b < 0) exit(-1);
        ++iter;
    } // k
    a = current.a;
    b = current.b;
    std::cout <<current.a<<"  "<<current.b<<"  " <<current.k<<"  "<<current.jacobian<< std::endl;
    current.jacobian = wave_frequency_equation(&a,&b,gamma2,subcritical,bc,kmax);
    current.a = a;
    current.b = b;
    std::cout <<current.a<<"  "<<current.b<<"  " <<current.k<<"  "<<current.jacobian<< std::endl;
  } // mode
  return 0;
}
