#include <iostream>
#include <gsl/gsl_sf_bessel.h>

using namespace std;

int
main (void)
{
  double x = 5.0;
  double y = gsl_sf_bessel_J0 (x);
  cout << "J0(" << x << ") = " << y << endl;
  return 0;
}
