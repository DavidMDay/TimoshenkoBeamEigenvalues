#include <iostream>
#include <cassert>
#include <algorithm> // std::swap
#include <cmath>
#include <cstdlib>
#include "wave_number.h"
#include "frequency_equation.h"

namespace {

struct dds {
  double dads = 0.0;
  double dbds = 0.0;
  double dkds = 0.0;
};

//timoshenko beam,  han et al 1999
//implicit relationship between the wave numbers
//p.951, s, gamma2 are derived from the parameters in equations (77) & (81)
//From the wave number equation, b is a function of a and k, b(a,k).
double getdadkanddbdk(double k,
                      double gamma2,
                      double a,
                      double b,
                      bool is_sub_critical,
                      BoundaryCondition bc,
                      double* dbdkptr) {

  frequnency_equation_parameters state;
  state.a = a;
  state.b = b;
  state.gamma2 = gamma2;
  state.is_sub_critical = is_sub_critical;
  state.bc = bc;
  double dadk= 0.0;
  if( state.b > 0. ) {
        double dfda = wave_equation_dfda(k,state.gamma2,state.a,state.b,state.is_sub_critical);
        double dfdb = wave_equation_dfdb(k,state.gamma2,state.a,state.b,state.is_sub_critical); // <-
        double dfdk = wave_equation_dfdk(k,state.gamma2,state.a,state.b,state.is_sub_critical); // <-
        double dbda = wave_equation_dbda(dfda, dfdb );

        double dbdk = wave_equation_dbdk(dfdk, dfdb );
        *dbdkptr = dbdk;
        auto value = frequency_equation(state);
        dadk= - value.F_b * dbdk/( value.F_a + value.F_b*dbda);
    } else {
        double dbdk = 0.;
        *dbdkptr = dbdk;
    }
    return dadk;
}

dds arclengthpredictor( continuation_state current,
                           bool subcritical,
                           double arclength,
                           double gamma2,
                           BoundaryCondition bc,
                           continuation_state previous) {
    dds x;
    x.dkds = 0.0;
    if (current.b == 0.0) {
        double dadk = 0.;
        double dbdk = 1. / (arclength * arclength);
        x.dkds = arclength / sqrt(1. + dadk * dadk + dbdk * dbdk);
        x.dads = dadk * x.dkds;
        x.dbds = dbdk * x.dkds;
    } else {
        double delta[] = {current.a - previous.a, current.b - previous.b, current.k - previous.k};
        double ptp = 0.0;
        for (double i : delta) {
          ptp += i * i;
        }
        if (ptp == 0.0) {
          double dbdk = 0.0;
          double dadk = getdadkanddbdk(current.k, gamma2, current.a, current.b, subcritical, bc, &dbdk);
          x.dkds = arclength / sqrt(1. + dadk * dadk + dbdk * dbdk);  // arclength
          x.dads = dadk * x.dkds;                                      // 0
          x.dbds = dbdk * x.dkds;                                      // 0
        } else {
          double norm = sqrt(ptp);
          x.dads = delta[0] * arclength / norm;
          x.dbds = delta[1] * arclength / norm;
          x.dkds = delta[2] * arclength / norm;
        }
    }
    return x;
}

void linearsolve3x3( const double input_matrix[][3], const double* input_rhs, double * lhs) {
    const size_t numcol = 3;
    const size_t numrow = 3;
    double matrix[numrow][numcol];
    double rhs[numrow];

    for( size_t row = 0; row<numrow; row++) {
       for( size_t col = 0; col<numcol; col++) {
           matrix[row][col] = input_matrix[row][col];
       }
    }

    for( size_t row = 0; row<numrow; row++) rhs[row] = input_rhs[row];

    size_t max = ((fabs(matrix[0][0]) > fabs(matrix[1][0])) ? 0 : 1); // gepp
    size_t perm[] = { 0, 1};
    if( fabs( matrix[2][0] ) > fabs( matrix[max][0] ) ) max = 2;
    if( max > 0 ) {
        for( size_t i=0;i<3;i++ ) std::swap(matrix[max][i], matrix[0][i]);
        perm[0] = max;
        std::swap(rhs[max], rhs[0]);
    }
    assert (fabs(matrix[0][0]) > 0.);
    matrix[1][0]= matrix[1][0]/matrix[0][0];
    matrix[2][0]= matrix[2][0]/matrix[0][0];

    for( size_t i=1;i<3;i++) {
        for( size_t j=1;j<3;j++) {
            matrix[i][j] -= matrix[i][0]*matrix[0][j];
        }
    }

    if( fabs(matrix[2][1]) > fabs(matrix[1][1]) ) {
        for( size_t i=0;i<3;i++ ) {
          std::swap(matrix[1][i], matrix[2][i]);
        }
        perm[1] = 2;
        std::swap(rhs[1], rhs[2]);
    }
    assert( fabs( matrix[1][1] ) > 0. );
    matrix[2][1]  = matrix[2][1]/matrix[1][1];
    matrix[2][2] -= matrix[2][1]*matrix[1][2];

    assert( fabs(matrix[2][2]) > 0. );
    double y[] = {rhs[0],rhs[1]-matrix[1][0]*rhs[0],rhs[2]-matrix[2][0]*rhs[0]};
    y[2] = y[2] - matrix[2][1]*y[1]; //[L\U,b] Ly=b,Ux=y

    for( size_t row = 0; row <3; row++ ) lhs[row] = y[row];

    for( int row = 2; row >= 0; row-- ) {
       for( size_t col = row+1; col < 3; col++)
           lhs[row] = lhs[row] - matrix[row][col]*lhs[col];

       lhs[row] = lhs[row]/matrix[row][row];
    }
}

double maxnorm3x3( const double matrix[][3]) {
    double maxrowsum = 0.;
    for( size_t row = 0; row<3; row++ ) {
        double row_sum = 0.;
        for( size_t col = 0; col<3; col++ ) row_sum += fabs( matrix[row][col] );
        maxrowsum =  ( row_sum > maxrowsum ) ? row_sum : maxrowsum ;
    }
    return maxrowsum;
}
} // anonymous namespace


continuation_state
arclengthcontinuation( continuation_state current,
                              bool subcritical,
                              double arclength,
                              double gamma2,
                              double kmax,
                              BoundaryCondition bc,
                              continuation_state previous) {
  assert( current.k < kmax );
  frequnency_equation_parameters parameter;
  parameter.gamma2 = gamma2;
  parameter.is_sub_critical = subcritical;
  parameter.bc = bc; // not used again until while loop
  auto ds = arclengthpredictor( current, parameter.is_sub_critical, arclength,
                           parameter.gamma2, parameter.bc, previous);
  continuation_state prior = current;
  if( prior.k + ds.dkds > kmax ) {
        double ratio = (kmax-prior.k)/ds.dkds;
        ds.dads *= ratio;
        ds.dbds *= ratio;
        ds.dkds *= ratio;
        arclength *= ratio;
  }
  current.a += ds.dads;
  current.b += ds.dbds;
  current.k += ds.dkds;  // predictor
  current.jacobian = 0.;

  double rhs[] = { 1., 1., 1.};
  double lhs[] = { 1., 1., 1.};
  double rhs_norm  = sqrt(3.0);
  double lhs_norm  = rhs_norm;
  double jacobian_matrix_norm = 1.;

    auto max_iter = static_cast<size_t>(current.a);
    size_t lower_bound = 20;
    max_iter = std::max(max_iter, lower_bound);
    size_t iteration = 0;

    while( lhs_norm > 1.e-13 && rhs_norm > jacobian_matrix_norm*1.e-15 && iteration < max_iter ) {

        //std::cout << iteration <<" lhs "<< lhs_norm <<" rhs "<<rhs_norm<<" J "<<jacobian_matrix_norm<< std::endl;
        parameter.a = current.a;
        parameter.b = current.b;
        frequnency_equation_values v  = frequency_equation(parameter);

        double dFdk = 0; // corrector
        double f = wave_number_residual(current.k,parameter.gamma2,parameter.a,parameter.b,parameter.is_sub_critical);

        //std::cout << "       value " << v.F << " residual "  << f << std::endl;
        double dfda = wave_equation_dfda(current.k,parameter.gamma2,parameter.a,parameter.b,parameter.is_sub_critical);
        double dfdb = wave_equation_dfdb(current.k,parameter.gamma2,parameter.a,parameter.b,parameter.is_sub_critical);
        double dfdk = wave_equation_dfdk(current.k,parameter.gamma2,parameter.a,parameter.b,parameter.is_sub_critical);
        current.jacobian = v.F_a*dfdb-v.F_b*dfda;
        double jacobianmatrix[3][3];
        jacobianmatrix[0][0] = v.F_a;   jacobianmatrix[0][1] = v.F_b;   jacobianmatrix[0][2] = dFdk;
        jacobianmatrix[1][0] = dfda; jacobianmatrix[1][1] = dfdb; jacobianmatrix[1][2] = dfdk;
        jacobianmatrix[2][0] = current.a - prior.a;
        jacobianmatrix[2][1] = current.b - prior.b;
        jacobianmatrix[2][2] = current.k - prior.k;
        jacobian_matrix_norm = maxnorm3x3(jacobianmatrix);
        rhs[0] = -v.F;
        rhs[1] = -f;
        rhs[2] = .5*(arclength * arclength -(current.a - prior.a)*(current.a - prior.a)
                                         -(current.b - prior.b)*(current.b - prior.b)
                                         -(current.k - prior.k)*(current.k - prior.k));
        linearsolve3x3(jacobianmatrix, rhs, lhs);

        current.a += lhs[0];
        current.b += lhs[1];
        current.k += lhs[2];
        rhs_norm = sqrt( rhs[0]*rhs[0]+ rhs[1]*rhs[1]+ rhs[2]*rhs[2]);
        lhs_norm = sqrt( lhs[0]*lhs[0]+ lhs[1]*lhs[1]+ lhs[2]*lhs[2]);
        iteration++;
    }

    if( iteration == max_iter ) {
        std::cerr << " fatal error: nonlinear solve diverged \n";
        exit(-1);
    }
    if( current.b < 0. ) current.b = 0.;
    return current;
}
