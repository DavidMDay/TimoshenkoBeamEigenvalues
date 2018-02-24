#include <iostream>
#include <assert.h>
#include <math.h>
#include "wave_number.h"
#include "frequency_equation.h"

double wave_frequency_equation( double* aptr, 
                              double* bptr,
                              double gamma2,
                              bool subcritical,
                              bool is_ff,
                              double k)
{
    double a = *aptr; // (a,b): F=0 & f=0
    double b = *bptr;
    double jacobian = 0.;
    assert( b > 0. );
    unsigned num_iter = 10;
    for( unsigned iteration = 0; iteration < num_iter; iteration++)
    {
            double Fa, Fb; // frequency equation gradients
            double F = frequency_equation(a,b,gamma2,subcritical,is_ff,&Fa,&Fb);

            double dfda = wave_equation_dfda(k,gamma2,a,b,subcritical);
            double dfdb = wave_equation_dfdb(k,gamma2,a,b,subcritical);
            //J = [Fa Fb;...   rhs = -[F;f];
            //     fa fb];     J [da;db] = rhs
            jacobian = Fa*dfdb-Fb*dfda;
            double residual = wave_number_residual(k,gamma2,a,b,subcritical); // f
            double da = (Fb*residual - dfdb*F)/jacobian;
            double db = (dfda*F-Fa*residual)/jacobian;
            a += da;
            b += db;
            double f = residual;
            if( iteration == num_iter-1 ) assert(  sqrt(F*F+f*f) < .001 );
            if( sqrt(F*F+f*f) < 1.e-8 ) break;

    }
    *aptr = a;
    *bptr = b;
    return jacobian;
}

