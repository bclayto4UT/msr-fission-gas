/********************************************************************************
GAUSS ELIMINATION

This package contains four sets of functions to perform:
1. LU Factorization
    luFactorize:
    Input:
        a: a Square matrix
        p: a Permutation array (to be modified)
        inplace: whether the factorized Matrix would replace the original Matrix
    Output:
        a std::map<Matrix, Permutation>
2.  LU Solve
    luSolve only works on a factorized Matrix
    solve works on any Matrix
3.  Determinant
    detFactoredMatrix only works on a factorized Matrix
    det works on any Matrix
4.  Inverse Matrix
    invFactoredMatrix only works on a factorized Matrix
    inv works on any Matrix
********************************************************************************/

#ifndef GAUSSELIM_INCLUDED
#define GAUSSELIM_INCLUDED

#include "matrix.h"

//enum ge_state {GE_SUCCESS, GE_SINGULAR, GE_BADDATA};

Matrix luFactorize(Matrix& a, Permutation& p, bool inplace=false);

Vector luSolve(const Matrix& a, const Permutation& p, Vector& x, bool inplace=false);
Vector solve(Matrix& a, Permutation& p, Vector& x, bool inplace=false);

double detFactoredMatrix(const Matrix& a, const Permutation& p) noexcept;
double det(Matrix& a) noexcept;

Matrix invFactoredMatrix(const Matrix& a, const Permutation& p);
Matrix inv(Matrix& a);

#endif
