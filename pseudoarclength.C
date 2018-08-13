#include <iostream>
#include <assert.h>
#include <algorithm>    // std::swap
#include <math.h>
#include <stdlib.h>
#include "wave_number.h"
#include "frequency_equation.h"

double getdadkanddbdk(double k,double gamma2,double a,double b,
                      bool is_sub_critical, boundary_condition bc, double* dbdkptr);

void linearsolve3x3( const double input_matrix[][3], const double* input_rhs, double * lhs);

double arclengthpredictor( double a, double b, double k, bool subcritical, double arclength,
                           double gamma2, boundary_condition bc, double *previous_state, double *dads,
                           double *dbds );

double maxnorm3x3( const double matrix[][3])
{
    double maxrowsum = 0.;
    for( unsigned row = 0; row<3; row++ )
    {
        double row_sum = 0.;
        for( unsigned col = 0; col<3; col++ ) row_sum += fabs( matrix[row][col] );
        maxrowsum =  ( row_sum > maxrowsum ) ? row_sum : maxrowsum ;
    }
    return maxrowsum;
}


double arclengthcontinuation( double *a,
                              double *b, 
                              double *k, 
                              bool subcritical,
                              double arclength,
                              double gamma2,
                              double kmax,
                              boundary_condition bc,
                              double *previous_state )
{
    if( *b < 0. ) *b = 0.;
    double dads = 0.;
    double dbds = 0.;
    double a_previous = *a;
    double b_previous = *b;
    double k_previous = *k;
    assert( k_previous < kmax );
    double dkds = arclengthpredictor( *a, *b, *k, subcritical, arclength,
                           gamma2, bc, previous_state, &dads, &dbds );
    if( k_previous + dkds > kmax )
    { 
        double ratio = (kmax-k_previous )/dkds;
        dads *= ratio;
        dbds *= ratio;
        dkds *= ratio;
        arclength = arclength * ratio;
    } 
    *k += dkds;  // predictor
    *a += dads;
    *b += dbds;
    double jacobian = 0.;

    double rhs[] = { 1., 1., 1.};
    double sol[] = { 1., 1., 1.};
    double rhs_norm  = sqrt( rhs[0]*rhs[0]+ rhs[1]*rhs[1]+ rhs[2]*rhs[2]);
    double lhs_norm  = sqrt( sol[0]*sol[0]+ sol[1]*sol[1]+ sol[2]*sol[2]);
    double jacobian_matrix_norm = 1.;

    unsigned max_iter = static_cast<unsigned>(*a);
    max_iter =  (max_iter>20) ?  max_iter : 20;

    unsigned iteration = 0;
    while( lhs_norm > 1.e-13 && rhs_norm > jacobian_matrix_norm*1.e-15 && iteration < max_iter )
    {
        double Fa, Fb; // frequency equation gradients
        double F = frequency_equation(*a,*b,gamma2,subcritical,bc,&Fa,&Fb);
        double dFdk = 0; // corrector
        double f = wave_number_residual(*k,gamma2,*a,*b,subcritical);
        double dfda = wave_equation_dfda(*k,gamma2,*a,*b,subcritical);
        double dfdb = wave_equation_dfdb(*k,gamma2,*a,*b,subcritical);
        double dfdk = wave_equation_dfdk(*k,gamma2,*a,*b,subcritical);
        jacobian = Fa*dfdb-Fb*dfda;
        double jacobianmatrix[3][3];
        jacobianmatrix[0][0] = Fa;   jacobianmatrix[0][1] = Fb;   jacobianmatrix[0][2] = dFdk;
        jacobianmatrix[1][0] = dfda; jacobianmatrix[1][1] = dfdb; jacobianmatrix[1][2] = dfdk;
        jacobianmatrix[2][0] = *a - a_previous; 
        jacobianmatrix[2][1] = *b - b_previous; 
        jacobianmatrix[2][2] = *k - k_previous;
        jacobian_matrix_norm = maxnorm3x3(jacobianmatrix);
        rhs[0] = -F;
        rhs[1] = -f;
        rhs[2] = .5*(arclength*arclength -(*a - a_previous)*(*a - a_previous)
                                         -(*b - b_previous)*(*b - b_previous)
                                         -(*k - k_previous)*(*k - k_previous));
        linearsolve3x3( jacobianmatrix, rhs, sol);

        *a = *a + sol[0];            *b = *b + sol[1];             *k = *k + sol[2];
        rhs_norm  = sqrt( rhs[0]*rhs[0]+ rhs[1]*rhs[1]+ rhs[2]*rhs[2]);
        lhs_norm  = sqrt( sol[0]*sol[0]+ sol[1]*sol[1]+ sol[2]*sol[2]);

        iteration = iteration + 1;

    }

    if( iteration == max_iter ) 
    {
        std::cerr << " fatal error: nonlinear solve diverged \n";
        exit(-1);
    }
    if( *b < 0. ) *b = 0.;
    return jacobian;
}
double arclengthpredictor( double a,
                           double b, 
                           double k, 
                           bool subcritical,
                           double arclength,
                           double gamma2,
                           boundary_condition bc,
                           double *previous_state,
                           double *dads,
                           double *dbds )
{
    double dkds = 0.;
    if( b == 0. )
    {
        double dadk = 0.;
        double dbdk = 1./(arclength*arclength);
        dkds = arclength/sqrt(1.+dadk*dadk+dbdk*dbdk);
       *dads = dadk * dkds;
       *dbds = dbdk * dkds;
    }
    else
    {
        double prev[] = {a - previous_state[0],b - previous_state[1],k - previous_state[2]};
        double ptp = 0.;
        for( unsigned i=0;i<3;i++) ptp += prev[i]*prev[i];
        if( ptp == 0.)
        {
            double dbdk = 0.;
            double dadk = getdadkanddbdk(k,gamma2,a,b,subcritical,bc,&dbdk);
            dkds = arclength/sqrt(1.+dadk*dadk+dbdk*dbdk); // arclength
           *dads = dadk * dkds; // 0
           *dbds = dbdk * dkds; // 0
        }
        else
        {
             double norm = sqrt(ptp);
            *dads = prev[0]*arclength/norm;
            *dbds = prev[1]*arclength/norm;
             dkds = prev[2]*arclength/norm;
        }
    }
    return dkds;
}

void linearsolve3x3( const double input_matrix[][3], const double* input_rhs, double * lhs)
{
    const unsigned numcol = 3;
    const unsigned numrow = 3;
    double matrix[numrow][numcol];
    double rhs[numrow];

    for( unsigned row = 0; row<numrow; row++)
       for( unsigned col = 0; col<numcol; col++)
           matrix[row][col] = input_matrix[row][col];

    for( unsigned row = 0; row<numrow; row++) rhs[row] = input_rhs[row];

    unsigned max = ((fabs(matrix[0][0]) > fabs(matrix[1][0])) ? 0 : 1); // gepp 
    unsigned perm[] = { 0, 1};
    if( fabs( matrix[2][0] ) > fabs( matrix[max][0] ) ) max = 2;
    if( max > 0 )
    {
        for( unsigned i=0;i<3;i++ ) std::swap(matrix[max][i], matrix[0][i]);
        perm[0] = max;
        std::swap(rhs[max], rhs[0]);
    }
    assert (fabs(matrix[0][0]) > 0.);
    matrix[1][0]= matrix[1][0]/matrix[0][0];
    matrix[2][0]= matrix[2][0]/matrix[0][0];

    for( unsigned i=1;i<3;i++)
        for( unsigned j=1;j<3;j++)
            matrix[i][j] -= matrix[i][0]*matrix[0][j];

    if( fabs(matrix[2][1]) > fabs(matrix[1][1]) )
    {
        for( unsigned i=0;i<3;i++ ) std::swap(matrix[1][i], matrix[2][i]);
        
        perm[1] = 2;
        std::swap(rhs[1], rhs[2]);
    }
    assert( fabs( matrix[1][1] ) > 0. );
    matrix[2][1]  = matrix[2][1]/matrix[1][1];
    matrix[2][2] -= matrix[2][1]*matrix[1][2]; 

    assert( fabs(matrix[2][2]) > 0. );
    double y[] = {rhs[0],rhs[1]-matrix[1][0]*rhs[0],rhs[2]-matrix[2][0]*rhs[0]};
    y[2] = y[2] - matrix[2][1]*y[1]; //[L\U,b] Ly=b,Ux=y

    for( unsigned row = 0; row <3; row++ ) lhs[row] = y[row];

    for( int row = 2; row >= 0; row-- )
    {
       for( unsigned col = row+1; col < 3; col++)
           lhs[row] = lhs[row] - matrix[row][col]*lhs[col];

       lhs[row] = lhs[row]/matrix[row][row];
    }
}
//timoshenko beam,  han et al 1999
//implicit relationship between the wave numbers
//p.951, s, gamma2 are derived from the parameters in equations (77) & (81)
//From the wave number equation, b is a function of a and k, b(a,k).

double getdadkanddbdk(double k,
                      double gamma2,
                      double a,
                      double b,
                      bool is_sub_critical,
                      boundary_condition bc,
                      double* dbdkptr)
{
    double dadk= 0.;
    if( b > 0. )
    {
        double dfda = wave_equation_dfda(k,gamma2,a,b,is_sub_critical);
        double dfdb = wave_equation_dfdb(k,gamma2,a,b,is_sub_critical); // <-
        double dfdk = wave_equation_dfdk(k,gamma2,a,b,is_sub_critical); // <-
        double dbda = wave_equation_dbda(dfda, dfdb );

        double dbdk = wave_equation_dbdk(dfdk, dfdb );
        *dbdkptr = dbdk;
        double Fa, Fb;
        double F = frequency_equation(a,b,gamma2,is_sub_critical,bc,&Fa,&Fb);

        dadk= - Fb * dbdk/( Fa + Fb*dbda);
    }
    else
    {
        double dbdk = 0.;
        *dbdkptr = dbdk;
    }
    return dadk;
}
