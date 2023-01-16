//#define XSTR(x) #x
//#define STR(x) XSTR(x)
//
//#define gsl C:/gsl/include/gsl
//
//#include STR(gsl/gsl_const.h)

#include <iostream>
#include <gsl/gsl_sf_bessel.h>

using namespace std;

int main(){
  double x = 5.0;
  double y = gsl_sf_bessel_J0 (x);
  cout << y << endl;
  return 0;
}
