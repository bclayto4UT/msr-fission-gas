#ifndef CALCULUS_H_INCLUDED
#define CALCULUS_H_INCLUDED

#include <cmath>
#include <functional>

#include "matrix.h"

// First derivatives
double deriv(std::function<double(double)>, const double);
Vector deriv(std::function<Vector(double)>, const double);
double pderiv(std::function<double(Vector&)>, const Vector&, int);
Vector pderiv(std::function<Vector(Vector&)>, const Vector&, int);

// First derivative operators
Vector grad(std::function<double(Vector&)>, const Vector&);
double div(std::function<Vector(const Vector&)>, const Vector&);
Vector curl(std::function<Vector(const Vector&)>, const Vector&);
Matrix jacobian(std::function<Vector(const Vector&)>, const Vector&);

// Second derivative
double deriv2(std::function<double(double)>, const double);
Vector deriv2(std::function<Vector(double)>, const double);
double pderiv2(std::function<double(Vector&)>, const Vector&, int, int);
Vector pderiv2(std::function<Vector(Vector&)>, const Vector&, int, int);

// Second derivative operators
double laplace(std::function<double(Vector&)>, const Vector&);
Vector laplace(std::function<Vector(const Vector&)>, const Vector&);
Matrix hessian(std::function<double(const Vector&)>, const Vector&);

// Single integrals
double intClose(std::function<double(double)>, double, double,
                double tol=1e-6, double mag = 1.0);
double intOpen (std::function<double(double)>, double, double,
                double tol=1e-6, double mag = 1.0);
double intGauss(std::function<double(double)>, double, double,
                double tol=1e-6, double mag = 1.0);
double integral(std::function<double(double)>, double, double,
                double tol=1e-6, double mag = 1.0);

#endif // CALCULUS_H_INCLUDED
