#include <math.h>

#include "calculus.h"
#include "miscellaneous.h"
#include "rootFinding.h"

//#define MONITOR

double bisection(std::function<double(double)> f,
                 double x0, double x1, int maxIter, double tol){
    if (fabs(f(x0)) <= tol) return x0;
    if (fabs(f(x1)) <= tol) return x1;
    if (tol <= 0) throw std::invalid_argument("ERROR: Invalid tolerance.");
    if (maxIter <= 0) maxIter = 1;
    if (x0 == x1) throw std::invalid_argument("ERROR: Two guesses needed.");
    if (sgn(f(x0)) * sgn(f(x1)) == 1)
        throw std::invalid_argument("ERROR: Two guesses must give opposite signs.");

    for (int i = 0; i < maxIter; i++){
        double mid = (x0+x1)/2;
        if (fabs(f(mid)) <= tol) return mid;
        if (sgn(f(mid)) == sgn(f(x0))){
            x0 = mid;
        } else x1 = mid;
    }
    return (x0+x1)/2;
}

double fixedpt(std::function<double(double)> f,
               double x0, int maxIter, double tol){
    if (fabs(x0 - f(x0)) <= tol) return x0;
    if (tol <= 0) throw std::invalid_argument("ERROR: Invalid tolerance.");
    if (maxIter <= 0) maxIter = 1;

    for (int i = 0; i < maxIter; i++){
        x0 = f(x0);
        if (!std::isfinite(x0)) throw std::domain_error("ERROR: Unable to evaluate function.");
        if (fabs(x0 - f(x0)) <= tol) return x0;
    }
    std::string err = "ERROR: Unable to converge within " + std::to_string(maxIter) + " iteration(s).";
    throw std::runtime_error(err);
}

double newton(std::function<double(double)> f,
              std::function<double(double)> df,
              double x0, int maxIter, double tol){
    if (fabs(f(x0)) <= tol) return x0;
    if (tol <= 0) throw std::invalid_argument("ERROR: Invalid tolerance.");
    if (maxIter <= 0) maxIter = 1;

    for (int i = 0; i < maxIter; i++){
        double dfx = df(x0);
        if (dfx == 0.0) throw std::runtime_error("ERROR: Zero derivative.");
        double dx = -f(x0)/dfx;
        if (fabs(dx) <= tol) return x0;
        x0 += dx;
        if (!std::isfinite(f(x0)))
            throw std::domain_error("ERROR: Unable to evaluate function.");
    }
    std::string err = "ERROR: Unable to converge within " + std::to_string(maxIter) + " iteration(s).";
    throw std::runtime_error(err);
}

double newton(std::function<double(double)> f,
              double x0, int maxIter, double tol){
    if (tol <= 0) throw std::invalid_argument("ERROR: Invalid tolerance.");
    if (maxIter <= 0) maxIter = 1;

    for (int i = 0; i < maxIter; i++){
        double dfx = deriv(f, x0);
        if (dfx == 0.0) throw std::runtime_error("ERROR: Zero derivative.");
        double dx = -f(x0)/dfx;

#ifdef MONITOR
        std::cout << "iter = " << i << ", x = " << x0
                  << ", f(x) = " << f(x0) << ", dx = " << dx
                  << ", f'(x) = " << dfx << std::endl;
#endif // MONITOR

        if (fabs(dx) <= tol) return x0;
        x0 += dx;
        if (!std::isfinite(f(x0)))
            throw std::runtime_error("ERROR: Unable to evaluate function.");
    }
    std::string err = "ERROR: Unable to converge within " + std::to_string(maxIter) + " iteration(s).";
    throw std::runtime_error(err);
}

double secant(std::function<double(double)> f,
              double x0, double x1, int maxIter, double tol){
    if (fabs(f(x0)) <= tol) return x0;
    if (fabs(f(x1)) <= tol) return x1;
    if (tol <= 0) throw std::invalid_argument("ERROR: Invalid tolerance.");
    if (maxIter <= 0) maxIter = 1;
    if (x0 == x1) throw std::invalid_argument("ERROR: Two guesses needed.");

    for (int i = 0; i < maxIter; i++){
        double dfx = f(x1)-f(x0);
        if (dfx == 0.0) throw std::runtime_error("ERROR: Zero denominator.");
        double dx = -f(x1)*(x1-x0)/dfx;

#ifdef MONITOR
        std::cout << "iter = " << i
                  << ", x0 = " << x0 << ", x1 = " << x1
                  << ", f(x) = " << f(x0) << ", dx = " << dx
                  << ", f'(x) = " << dfx << std::endl;
#endif // MONITOR

        if (fabs(dx) <= tol) return x0;
        x0 = x1;
        x1 += dx;
        if (!std::isfinite(f(x1)))
            throw std::runtime_error("ERROR: Unable to evaluate function.");
    }
    std::string err = "ERROR: Unable to converge within " + std::to_string(maxIter) + " iteration(s).";
    throw std::runtime_error(err);
}

