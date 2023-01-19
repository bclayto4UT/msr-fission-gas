#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

#include "gaussElim.h"

//#define TEST
//#define LARGE_TEST     // to test runtime

#ifdef TEST

int main(void) {

  int n;
  cout << "Enter size of n by n Matrix A: " << flush;
  cin >> n;

  Matrix A(n,n), AA(n,n);
  Vector x(n);
  Vector b(n);
  Permutation p(n);

  cout << "Enter A by rows: " << flush;
  cin >> A;
  // cout << "A = " << endl << A << endl;

  cout << "Enter b: " << flush;
  cin >> b;
  // cout << "b = " << b << endl;

  x  = b;
  AA = A;

  ge_state s = solve(A,p,x);

  switch(s) {
  case GE_SINGULAR:
    cout << "ERROR: Singular matrix" << endl;
    return 1;
  case GE_BADDATA:
    cout << "ERROR: Bad data (wrong sizes)" << endl;
    return 1;
  case GE_SUCCESS:
    cout << endl << "The solution is:" << endl << x << endl;

    cout << "The permutation is:" << endl << p << endl;
    cout << "The parity is: " << p.parity() << endl;
    cout << "The factored matrix is:" << endl << AA << endl;

    // Compute AA*x - b (= 0?)
    Vector r(n);
    matVecMult(AA,x,r);
    r -= b;
    cout << "The norm of the error is: " << l2norm(r) << endl;
    return 0;
  }
  return 2;
}
#endif

#ifdef LARGE_TEST
#include <ctime>

int main(void) {
  int n;
  cout << "Enter n: ";
  cin >> n;

  cout << "SETUP:" << endl;

  Matrix A(n,n);
  Vector x(n);
  Permutation p(n);

  A = 0.0;
  for(int j=0; j<n; j++) {
    A(j,j) = 1;
    x(j) = j;
  }

  cout << "SOLVE:" << endl;

  clock_t t = clock();
  ge_state s = solve(A,p,x);

  switch(s) {
  case GE_SINGULAR:
    cout << "ERROR: Singular matrix" << endl;
    return 1;
  case GE_BADDATA:
    cout << "ERROR: Bad data (wrong sizes)" << endl;
    return 1;
  case GE_SUCCESS:
    t = clock() - t;
    cout << "FINISHED.  Total time = " << ((float)t)/CLOCKS_PER_SEC << endl;
    return 0;
  }
  return 2;
}
#endif
