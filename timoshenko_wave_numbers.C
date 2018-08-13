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

void timoshenko_wave_numbers( double kmax,
                              double gamma2,
                              unsigned mode,
                              boundary_condition bc,
                              double* aptr,
                              double* bptr)
{
    double arclength = .01;
    double factor = sqrt( (gamma2 + 1.)/gamma2 );
    double jacobian = 0.;
    double a = euler_wave_number( mode, bc );
    double b = a;
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
            //verbose = true;
            double a_critical = get_critical_point2( a, gamma2, bc);
            //std::cout << "a_critical= " << a_critical << "\n";
            if( verbose ) std::cout << "a := a_critical\n";
            a=a_critical;
            next_subcritical = false;
        }
        if( subcritical != next_subcritical ) // b=0 singularity
        {
            double a_critical = get_critical_point2( a, gamma2, bc);
            a = a_critical;
            b = 0.;
            k = factor/a;
            //std::cout << "singularity " << k <<  "  "  << a << "  " << b << "\n";
            subcritical = false;
            verbose = true;
        }
        double scaled_frequency = 0.;
        if( b < 0. ) exit(-1);
        //std::cout << k << " " <<  a << "  " << b << "\n";
    } // k
    wave_frequency_equation(&a,&b,gamma2,subcritical,bc,kmax);
    *aptr =  a;
    *bptr =  b;
}
