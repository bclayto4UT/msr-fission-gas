#include "calculus.h"
#include "matrix.h"
#include "miscellaneous.h"
//#include "polynomial.h"
//#include "rootFinding.h"

#include <cmath>
#include <map>

double deriv(std::function<double(double)> f, const double x){
    double h = pow(3e-16, 1.0/3) * pow(10, ordMag(x));
    double fxp = f(x+h), fxm = f(x-h);
    if (!std::isfinite(fxp) || !std::isfinite(fxm))
        throw std::runtime_error("ERROR in deriv: Unable to evaluate function.");
    return (fxp-fxm)/(2*h);
}

Vector deriv(std::function<Vector(double)> f, const double x){
    double h = pow(3e-16, 1.0/3) * pow(10, ordMag(x));
    Vector fxp = f(x+h), fxm = f(x-h);
    if (!fxp.isFinite() || !fxm.isFinite())
        throw std::runtime_error("ERROR in deriv: Unable to evaluate function.");
    return (fxp-fxm)/(2*h);
}

double pderiv(std::function<double(Vector&)> f, const Vector& x, int var){
    if (var < 0 || var >= x.n())
        throw std::length_error("ERROR in pderiv: Invalid variable index.");

    double h = pow(3e-16, 1.0/3) * pow(10, ordMag(x(var)));
    Vector xp(x), xm(x);
    xp(var) += h; xm(var) -= h;
    double fxp = f(xp), fxm = f(xm);
    if (!std::isfinite(fxp) || !std::isfinite(fxm))
        throw std::runtime_error("ERROR in pderiv: Unable to evaluate function.");
    return (fxp-fxm)/(2*h);
}

Vector pderiv(std::function<Vector(Vector&)> f, const Vector& x, int var){
    if (var < 0 || var >= x.n())
        throw std::length_error("ERROR in pderiv: Invalid variable index.");

    double h = pow(3e-16, 1.0/3) * pow(10, ordMag(x(var)));
    Vector xp(x), xm(x);
    xp(var) += h; xm(var) -= h;
    Vector fxp = f(xp), fxm = f(xm);
    if (!fxp.isFinite() || !fxm.isFinite())
        throw std::runtime_error("ERROR in pderiv: Unable to evaluate function.");
    return (fxp-fxm)/(2*h);
}

Vector grad(std::function<double(Vector&)> f, const Vector& x){
    Vector g(x.n());
    for (int i = 0; i < x.n(); i++){
        try{
            g(i) = pderiv(f, x, i);
        } catch (const std::runtime_error& ex){
            throw ex;
        }
    }
    return g;
}

double div(std::function<Vector(const Vector&)> f, const Vector& x){
    Vector fx = f(x);
    if (fx.n() != x.n()) throw std::domain_error("ERROR: Incompatible sizes.");

    double d = 0.0;
    for (int i = 0; i < x.n(); i++){
        try{
            Vector df = pderiv(f, x, i);
            d += df(i);
        } catch (const std::runtime_error& ex){
            throw ex;
        }
    }
    return d;
}

Vector curl(std::function<Vector(const Vector&)> f, const Vector& x){
    Vector fx = f(x);
    if (fx.n() > 3 || x.n() > 3)
        throw std::domain_error("ERROR: Size greater than 3.");

    Vector c(3);
    Vector x0(x); x0.resize(3);

    try{
        Vector f0 = pderiv(f, x0, 0); f0.resize(3);
        Vector f1 = pderiv(f, x0, 1); f1.resize(3);
        Vector f2 = pderiv(f, x0, 2); f2.resize(3);
        c(0) = f1(2) - f2(1);
        c(1) = f2(0) - f0(2);
        c(2) = f0(1) - f1(0);
    } catch (const std::runtime_error& ex){
        throw ex;
    }
    return c;
}

Matrix jacobian(std::function<Vector(const Vector&)> f, const Vector& x){
/** Inputs:
        f              The name of the function whose Jacobian matrix is sought.
        x              The initial guess at the solution.
    Outputs:
        Matrix         The Jacobian matrix of f estimated at x
*/

    Vector fx = f(x);
    if (!fx.isFinite()) throw std::domain_error("ERROR: Unable to evaluate function.");
    int n = x.n();
    int m = fx.n();

    Matrix jac(m, n);
    for (int i = 0; i < n; i++){
        try{
            Vector df = pderiv(f, x, i);
            jac.assign(df, i);
        } catch (const std::logic_error& ex){
            throw ex;
        } catch (const std::runtime_error& ex){
            throw ex;
        }
    }
    return jac;
}

double deriv2(std::function<double(double)> f, const double x){
    double h = pow(6e-16, 1.0/3) * pow(10, ordMag(x));
    double fxp = f(x+h), fx = f(x), fxm = f(x-h);
    if (!std::isfinite(fxp) || !std::isfinite(fx) || !std::isfinite(fxm))
        throw std::runtime_error("ERROR: Unable to evaluate function.");
    return (fxp-2*fx+fxm)/(h*h);
}

Vector deriv2(std::function<Vector(double)> f, const double x){
    double h = pow(6e-16, 1.0/3) * pow(10, ordMag(x));
    Vector fxp = f(x+h), fx = f(x), fxm = f(x-h);
    if (!fxp.isFinite() || !fx.isFinite() || !fxm.isFinite())
        throw std::runtime_error("ERROR: Unable to evaluate function.");
    return (fxp-fx*2+fxm)/(h*h);
}

double pderiv2(std::function<double(Vector&)> f, const Vector& x, int var1, int var2){
    if (var1 < 0 || var2 < 0 || var1 >= x.n() || var2 >= x.n())
        throw std::length_error("ERROR: Invalid variable index(es).");

    if (var1 == var2){ // Unmixed partial derivatives
        double h = pow(6e-16, 1.0/3) * pow(10, ordMag(x(var1)));
        Vector xp(x), x0(x), xm(x);
        xp(var1) += h; xm(var1) -= h;
        double fxp = f(xp), fx = f(x0), fxm = f(xm);
        if (!std::isfinite(fxp) || !std::isfinite(fx) || !std::isfinite(fxm))
            throw std::runtime_error("ERROR: Unable to evaluate function.");
        return (fxp-fx*2+fxm)/(h*h);

    } else{ // Mixed partial derivatives
        double h1 = pow(6e-16, 1.0/3) * pow(10, ordMag(x(var1)));
        double h2 = pow(6e-16, 1.0/3) * pow(10, ordMag(x(var2)));
        Vector xpp(x), xpm(x), xmp(x), xmm(x);
        xpp(var1) += h1; xpp(var2) += h2;
        xpm(var1) += h1; xpm(var2) -= h2;
        xmp(var1) -= h1; xmp(var2) += h2;
        xmm(var1) -= h1; xmm(var2) -= h2;

        double fpp = f(xpp), fpm = f(xpm), fmp = f(xmp), fmm = f(xmm);
        if (!std::isfinite(fpp) || !std::isfinite(fpm)
            || !std::isfinite(fmp) || !std::isfinite(fmm))
            throw std::runtime_error("ERROR: Unable to evaluate function.");
        return (fpp-fpm-fmp+fmm)/(4*h1*h2);
    }
}

Vector pderiv2(std::function<Vector(Vector&)> f, const Vector& x, int var1, int var2){
    if (var1 < 0 || var2 < 0 || var1 >= x.n() || var2 >= x.n())
        throw std::length_error("ERROR: Invalid variable index(es).");

    if (var1 == var2){ // Unmixed partial derivatives
        double h = pow(6e-16, 1.0/3);
        Vector xp(x), x0(x), xm(x);
        xp(var1) += h; xm(var1) -= h;
        Vector fxp = f(xp), fx = f(x0), fxm = f(xm);
        if (!fxp.isFinite() || !fx.isFinite() || !fxm.isFinite())
            throw std::runtime_error("ERROR: Unable to evaluate function.");
        return (fxp-fx*2+fxm)/(h*h);

    } else{ // Mixed partial derivatives
        double h1 = pow(6e-16, 1.0/3) * pow(10, ordMag(x(var1)));
        double h2 = pow(6e-16, 1.0/3) * pow(10, ordMag(x(var2)));
        Vector xpp(x), xpm(x), xmp(x), xmm(x);
        xpp(var1) += h1; xpp(var2) += h2;
        xpm(var1) += h1; xpm(var2) -= h2;
        xmp(var1) -= h1; xmp(var2) += h2;
        xmm(var1) -= h1; xmm(var2) -= h2;

        Vector fpp = f(xpp), fpm = f(xpm), fmp = f(xmp), fmm = f(xmm);
        if (!fpp.isFinite() || !fpm.isFinite()
            || !fmp.isFinite() || !fmm.isFinite())
            throw std::runtime_error("ERROR: Unable to evaluate function.");
        return (fpp-fpm-fmp+fmm)/(4*h1*h2);
    }
}

double laplace(std::function<double(Vector&)> f, const Vector& x){
    double d = 0.0;
    for (int i = 0; i < x.n(); i++){
        try{
            Vector df = pderiv2(f, x, i, i);
            d += df(i);
        } catch (const std::runtime_error& ex){
            throw ex;
        }
    }
    return d;
}

Vector laplace(std::function<Vector(const Vector&)> f, const Vector& x){
    Vector fx = f(x);
    if (fx.n() != x.n()) throw std::domain_error("ERROR: Incompatible sizes.");

    Vector l(x.n());
    for (int i = 0; i < l.n(); i++){
        try{
            l += pderiv2(f, x, i, i);
        } catch (const std::runtime_error& ex){
            throw ex;
        }
    }
    return l;
}

Matrix hessian(std::function<double(const Vector&)> f, const Vector& x){
    int n = x.n();
    Matrix H(n,n);
    for (int i = 0; i < n; i++){
        try{
            //H(i,i) = pderiv2(f, x, i, i);
        } catch (const std::runtime_error& ex){
                throw ex;
        }
        for (int j = i+1; j < n; j++){
            try{
                //H(i,j) = pderiv2(f, x, i, j);
                //H(j,i) = H(i,i);
            } catch (const std::runtime_error& ex){
                throw ex;
            }
        }
    }
    return H;
}


//double intClose(std::function<double(double)> f,
//                double x0, double x1, double tol, double mag){
//
//    if (tol<0) throw std::invalid_argument("ERROR: Invalid tolerance.");
//    if (!std::isfinite(f(x0)) || !std::isfinite(f(x1)))
//        throw std::domain_error("ERROR: Cannot evaluate function.");
//    if (x0 == x1) return 0;
//
//    double h = pow(tol*180/(x1-x0)/mag, 0.25);
//    int step = ceil((x1-x0)/h);
//    while (step % 2 != 0) step++;
//    h = (x1-x0)/step;
//
//    double sum;
//    for (int i = 0; i < step/2; i++){
//        sum += h/3 * (f(x0) + 4*f(x0+h) + f(x0+2*h));
//        if (!std::isfinite(sum))
//            throw std::runtime_error("ERROR: Unable to evaluate function.");
//        x0 += 2*h;
//    }
//    return sum;
//}
//
//double intOpen(std::function<double(double)> f,
//               double x0, double x1, double tol, double mag){
//
//    if (tol<0) throw std::invalid_argument("ERROR: Invalid tolerance.");
//    if (x0 == x1) return 0;
//
//    double h = pow(tol*90/7/(x1-x0)/mag, 0.25);
//    int step = ceil((x1-x0)/h);
//    while (step % 4 != 0) step++;
//    h = (x1-x0)/step;
//
//    double sum;
//    for (int i = 0; i < step/4; i++){
//        sum += 4*h/3 * (2*f(x0+h) - f(x0+2*h) + 2*f(x0+3*h));
//        if (!std::isfinite(sum))
//            throw std::runtime_error("ERROR: Unable to evaluate function.");
//        x0 += 4*h;
//    }
//    return sum;
//}
//
//// Gaussian Quadrature and helper functions
//static std::vector<std::map<double, double>> legendre;
////static int maxPanelDeg = 5;
//
//double intGauss(std::function<double(double)> f,
//                double x0, double x1, double tol, double mag){
//    // Find the Legendre polynomial degree
//    auto gaussErr = [&](double n)
//    {
//        double error = pow(x1-x0, 2*n+1) * pow(factorial(n), 4) * mag;
//        return error / (2*n+1) / pow(factorial(2*n), 3) - tol;
//    };
//
//    int n = 2;
//    while (gaussErr(n) > 0){ n = round(n*1.3);}
//
//    /*if (n > maxPanelDeg){
//        n = n/maxPanelDeg+1;
//        double sum = 0;
//        double h = (x1-x0)/n;
//        for (int i = 0; i < n; i++){
//            sum += intGauss(f, x0, x0+h, tol, mag);
//            x0 += h;
//        }
//        return sum;
//    }*/
//
//    // Finds Legendre roots and corresponding coefficients
//    if (legendre.size() < n || legendre[n-1].empty()){
//
//        if (legendre.size() < n) legendre.resize(n);
//        polynomial p{-1, 0, 1};
//        p = p.expo(n).deriv(n) / pow(2,n) / factorial(n);
//
//        for (double r = -1.0; r <= 1.0; r += 1.0/n){
//            p.root(r);
//            if (p.getRootValidity()) break;
//        }
//
//        auto lRoots = p.getRoots();
//
//        for (double i: lRoots){
//
//            polynomial L(1);
//            for (auto j: lRoots){
//                if (i != j) L *= polynomial{-j, 1} / (i-j);
//            }
//            legendre[n-1].emplace(i, L.integ(-1, 1));
//        }
//    }
//
//    // Carries out integration
//    double sum;
//    for (auto i: legendre[n-1]){
//        double x = i.first*(x1-x0)/2.0 + (x1+x0)/2.0;
//        sum += i.second*f(x);
//    }
//    return sum*(x1-x0)/2.0;
//}
//
//double integral(std::function<double(double)> f,
//                double x0, double x1, double tol, double mag){
//    try{
//        return intGauss(f, x0, x1, tol, mag);
//    } catch (...){
//    }
//}
