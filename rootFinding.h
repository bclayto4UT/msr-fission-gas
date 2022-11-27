#ifndef ROOTFINDING_H_INCLUDED
#define ROOTFINDING_H_INCLUDED

#include <functional>

double bisection(std::function<double(double)> f,
                 double x0, double x1, int maxIter=20, double tol=1e-6);
double fixedpt(std::function<double(double)> f,
               double x0, int maxIter=20, double tol=1e-6);
double newton(std::function<double(double)> f,
              std::function<double(double)> df,
              double x0, int maxIter=20, double tol=1e-6);
double newton(std::function<double(double)> f,
              double x0, int maxIter=20, double tol=1e-6);
double secant(std::function<double(double)> f,
              double x0, double x1, int maxIter=20, double tol=1e-6);

#endif // ROOTFINDING_H_INCLUDED
