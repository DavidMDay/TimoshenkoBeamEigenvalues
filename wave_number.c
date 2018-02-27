//        std::cout << std::scientific << std::setprecision(14);
#include "wave_number.h"
//timoshenko beam,  implicit relationship between the wave numbers
//p.951, s, gamma2 are derived from the parameters in equations (77) & (81)
// This is the residual for the wave equation constraint ...
// that must be satisfied for (a,b) pairs in terms of k and gamma2.
double wave_number_residual(double k,double gamma2,double a,double b, bool subcritical)
{
    double a2 = a*a;
    double b2 = b*b;
    double k2 = k*k;
    double g2 = gamma2;
    double g4 = gamma2*gamma2;
    double a4 = a2*a2;
    double b4 = b2*b2;
    double residual = 0.;
    if( subcritical )
        residual = -k2*g2*a4 - k2*(g4+1.)*a2*b2 + a2*(g2+1.) - k2*g2*b4 - (g2+1.)*b2;
    else
        residual = -k2*g2*a4 + k2*(g4+1.)*a2*b2 + a2*(g2+1.) - k2*g2*b4 + (g2+1.)*b2;
    return residual;
}
