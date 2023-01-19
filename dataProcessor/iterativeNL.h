/**
File: iterativeNL.h
**/

#ifndef IterativeNL_Included
#define IterativeNL_Included

#include <functional>
#include "matrix.h"

Vector fixedpt(std::function<Vector(const Vector&)> g,
               Vector& x, int maxIter, double tol);

Vector newton(std::function<Vector(const Vector&)> f,
              std::function<Matrix(const Vector&)> df,
              Vector& x, int maxIter, double tol);

Vector newton(std::function<Vector(const Vector&)> f,
              Vector& x, int maxIter, double tol);

Vector broyden1(std::function<Vector(const Vector&)> f,
                Matrix& A, Vector& x, int maxIter, double tol);

Vector broyden2(std::function<Vector(const Vector&)> f,
                Matrix& B, Vector& x, int maxIter, double tol);

#endif
