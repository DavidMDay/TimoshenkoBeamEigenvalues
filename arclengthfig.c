#include <iostream>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "wave_number.h"
#include "frequency_equation.h"
#include "secular.h"
double wave_frequency_equation( double* aptr,
                              double* bptr,
                              double gamma2,
                              bool subcritical,
                              boundary_condition bc,
                              double k);


double arclengthcontinuation( double *a,
                              double *b, 
                              double *k, 
                              bool subcritical,
                              double arclength,
                              double gamma2,
                              double kmax,
                              boundary_condition bc,
                              double *previous_state );
int main()
{
  double kmax = .4;// .2273; // .4;// physics
  double gamma = 2.205;
  boundary_condition bc = freefree;
  double arclength = .01;
  double gamma2 = gamma*gamma;
  double factor = sqrt( (gamma2 + 1.)/gamma2 );
  double jacobian = 0.;
  for(unsigned mode = 4; mode < 5; mode++)
  {
    double a = euler_wave_number( mode, bc );
    double b = a;
    double a_critical = get_critical_point( mode+1, gamma2, bc);
    double k = 0.;
    bool subcritical = true;
    bool verbose = false;
    double state_before_last[] = {a,b,k};
    double previous_state[] = {a,b,k};
    while( k < kmax ) 
    {
        for(unsigned i=0; i<3;i++) state_before_last[i] = previous_state[i];
        previous_state[0] = a; previous_state[1] = b; previous_state[2] = k;

        double jacobian = arclengthcontinuation(&a,&b,&k,subcritical,arclength,
                                                gamma2,kmax,bc,state_before_last);
        bool next_subcritical = ( factor-a*k > 0. );
        if( b == 0. )
        {
            if( verbose ) std::cout << "a := a_critical\n";
            a=a_critical;
            next_subcritical = false;
        }
        if( subcritical != next_subcritical ) // b=0 singularity
        {
            a = a_critical;
            b = 0.;
            k = factor/a;
            subcritical = false;
            verbose = false;
        }
        double scaled_frequency = 0.;
        //if( subcritical )
        //   scaled_frequency = sqrt(( a*a - b*b )/(1+gamma2));// output state
        //else
        //   scaled_frequency = sqrt(( a*a + b*b )/(1+gamma2));
        if( k < kmax )
            std::cout <<k<<"  "<<a<<"  "<<b<<"  "<<jacobian<<"\n";
        if( b < 0 ) exit(-1);
    } // k
    jacobian = wave_frequency_equation(&a,&b,gamma2,subcritical,bc,kmax);
    std::cout <<kmax<<"  "<<a<<"  "<<b<<"  "<<jacobian<<"\n";
  } // mode
  return 0;
}
