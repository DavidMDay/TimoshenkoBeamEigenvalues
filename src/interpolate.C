#include <cassert>
#include <iostream>
// f(xi) = interpolate(n, x(0:n-1), f(0:n-1), xi),  0 <= n <= 4
//  Low order polynomial interpolation by Newton interpolation.
//  Purpose: to interpolate the eigenvalues at specified aspect ratios(?)
//  from near by values determined by continuation.
//  What if higher order interpolation is needed?
//  For high order barycentric Lagrange interpolation is preferable; otoh
//  one must not interpolate through a singularity.
double interpolate(size_t n, const double* x, double* f, double xi) {
  assert(n < 5);          // low order
  assert(x[0] < xi);      // inter..
  assert(xi < x[n - 1]);  // ... polation
  if (n == 0) return f[0];
  double fx0x1 = (f[1] - f[0]) / (x[1] - x[0]);
  double y = f[0];
  double L = (xi - x[0]);
  y = y + fx0x1 * L;
  if (n == 1) return y;
  double fx1x2 = (f[2] - f[1]) / (x[2] - x[1]);
  double fx0x1x2 = (fx1x2 - fx0x1) / (x[2] - x[0]);
  L *= (xi - x[1]);
  y += fx0x1x2 * L;
  if (n == 2) return y;
  double fx2x3 = (f[3] - f[2]) / (x[3] - x[2]);
  double fx1x2x3 = (fx2x3 - fx1x2) / (x[3] - x[1]);
  double fx0x1x2x3 = (fx1x2x3 - fx0x1x2) / (x[3] - x[0]);
  L *= (xi - x[2]);
  y += fx0x1x2x3 * L;
  if (n == 3) return y;

  double fx3x4 = (f[4] - f[3]) / (x[4] - x[3]);
  double fx2x3x4 = (fx3x4 - fx2x3) / (x[4] - x[2]);
  double fx1x2x3x4 = (fx2x3x4 - fx1x2x3) / (x[4] - x[1]);
  double fx0x1x2x3x4 = (fx1x2x3x4 - fx0x1x2x3) / (x[4] - x[0]);
  L *= (xi - x[3]);
  y += fx0x1x2x3x4 * L;
  return y;
}
