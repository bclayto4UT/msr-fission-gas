#include <cmath>
#include <iostream>

#include "calculus.h"
#include "iterativeNL.h"
#include "matrix.h"
#include "gaussElim.h"

//#define MONITOR  // output history of solution and l_inf error

static bool checkTol(const Vector& x, const Vector& dx, const double tol){
    for (int i = 0; i < x.n(); i++){
        if (fabs(dx(i)/x(i)) > tol && fabs(x(i)) > 1e-100) return false;
    }
    return true;
}

///////////////////////////////////////////////////////////////////////////////
Vector fixedpt(std::function<Vector(const Vector&)> g,
               Vector& x, int maxIter, double tol) {
/**  Inputs:
        g              The name of the function for which a fixed point is sought.
        x              The initial guess at the solution.
        tol            The convergence tolerance (must be > 0).
        maxIter        The maximum number of iterations that can be taken.
    Outputs:
        Vector         The solution vector.
*/

    if(maxIter < 1) maxIter = 1;
    if(tol <= 0) throw std::invalid_argument("ERROR: Invalid tolerance.");

#ifdef MONITOR
    std::cout << "Guess: x= " << x << std::endl;
#endif

    Vector y(x);
    for(int iter = 1; iter <= maxIter; iter++) {
        x = g(y);      // x = g(y) is the new guess
        if (!x.isFinite()) throw std::runtime_error("ERROR: Unable to evaluate function.");
        y -= x;

#ifdef MONITOR
    std::cout << "Iter " << iter << ": x= " << x << ", err = " << maxNorm(y) << std::endl;
#endif
        if (checkTol(x, y, tol)) return x;
        y = x;
    }
    std::string err = "ERROR: Unable to converge within " + std::to_string(maxIter) + " iteration(s).";
    throw std::runtime_error(err);
}


///////////////////////////////////////////////////////////////////////////////
Vector newton(std::function<Vector(const Vector&)> f,
              std::function<Matrix(const Vector&)> df,
              Vector& x, int maxIter, double tol){
/**  Inputs:
        f              The name of the function for which a root is sought.
        df             The name of the Jacobian of the function f.
        x              The initial guess at the solution.
        tol            The convergence tolerance (must be > 0).
        maxIter        The maximum number of iterations that can be taken.
    Outputs:
        Vector         The solution Vector.
*/

    if(maxIter < 1) maxIter = 1;
    if(tol <= 0) throw std::invalid_argument("ERROR: Invalid tolerance.");

#ifdef MONITOR
    std::cout << "Guess: x= " << x << std::endl;
#endif

    int n = x.n();
    for(int iter = 1; iter <= maxIter; iter++) {
        Matrix jac = df(x);
        Vector y = f(x);
        if (!y.isFinite()) throw std::runtime_error("ERROR: Unable to evaluate function.");
#ifdef MONITOR
    std::cout << std::endl << "Iter " << iter << ": fx=" << y
              << ", jac= " << std::endl << jac << std::endl;
#endif

        try{
            Permutation p(n);
            solve(jac, p, y, 1);
            if (!y.isFinite()) throw std::runtime_error("ERROR: Unable to solve for new x.");
            x -= y;  // x is the new guess

#ifdef MONITOR
    std::cout << "Iter " << iter << ": step= " << y << std::endl;
    std::cout << "Iter " << iter << ": x= " << x << ", err = " << maxNorm(y) << std::endl;
#endif
            if (checkTol(x, y, tol)) return x;
        } catch (const std::logic_error& ex){
            throw ex;
        } catch (const std::runtime_error& ex){
            throw ex;
        }
    }

    std::string err = "ERROR: Unable to converge within " + std::to_string(maxIter) + " iteration(s).";
    throw std::runtime_error(err);
}


///////////////////////////////////////////////////////////////////////////////
Vector newton(std::function<Vector(const Vector&)> f,
              Vector& x, int maxIter, double tol){
/**  Inputs:
        f              The name of the function for which a root is sought.
        x              The initial guess at the solution.
        tol            The convergence tolerance (must be > 0).
        maxIter        The maximum number of iterations that can be taken.
    Outputs:
        Vector         The solution Vector.
*/

    if(maxIter < 1) maxIter = 1;
    if(tol <= 0) throw std::invalid_argument("ERROR: Invalid tolerance.");

#ifdef MONITOR
    std::cout << "Guess: x= " << x << std::endl;
#endif

    int n = x.n();
    for(int iter = 1; iter <= maxIter; iter++) {
        try{
            Matrix jac = jacobian(f,x);
            Vector y = f(x);
            if (!y.isFinite()) throw std::runtime_error("ERROR: Unable to evaluate function.");

#ifdef MONITOR
    std::cout << std::endl << "Iter " << iter << ": fx=" << y
              << ", jac= " << std::endl << jac << std::endl;
#endif
            Permutation P(n);
            solve(jac, P, y, 1);
            if (!y.isFinite()) throw std::runtime_error("ERROR: Unable to solve for new x.");

            x -= y;  // x is the new guess
#ifdef MONITOR
    std::cout << "Iter " << iter << ": step= " << y << std::endl;
    std::cout << "Iter " << iter << ": x= " << x << ", err = " << maxNorm(y) << std::endl;
#endif

            if (checkTol(x, y, tol)) return x;
        } catch (const std::logic_error& ex){
            throw ex;
        } catch (const std::runtime_error& ex){
            throw ex;
        }
    }

    std::string err = "ERROR: Unable to converge within " + std::to_string(maxIter) + " iteration(s).";
    throw std::runtime_error(err);
}


///////////////////////////////////////////////////////////////////////////////
Vector broyden1(std::function<Vector(const Vector&)> f,
                Matrix& A, Vector& x, int maxIter, double tol){
/**  Inputs:
        f              The name of the function for which a root is sought.
        A              The initial matrix, approximation of the Jacobi.
        x              The initial guess at the solution.
        tol            The convergence tolerance (must be > 0).
        maxIter        The maximum number of iterations that can be taken.
     Outputs:
        Vector         The solution Vector.
*/


    if(maxIter < 1) maxIter = 1;
    if(tol <= 0) throw std::invalid_argument("ERROR: Invalid tolerance.");
    int n = x.n();
    if (A.n(0) != n || A.n(1) != n || f(x).n() != n)
        throw std::invalid_argument("ERROR: Incompatible matrix-vector sizes.");

#ifdef MONITOR
    std::cout << "Guess: x= " << x << ", A= " << std::endl << A << std::endl;
#endif

    for(int iter = 1; iter <= maxIter; iter++) {
        Vector fx = f(x);
        if (!fx.isFinite()) throw std::runtime_error("ERROR: Unable to evaluate function.");
#ifdef MONITOR
    std::cout << std::endl << "Iter " << iter << ": fx=" << fx
              << ", A= " << std::endl << A << std::endl;
#endif

        try{
            Permutation P(n);
            Vector del = -solve(A, P, fx);
            if (!del.isFinite()) throw std::runtime_error("ERROR: Unable to solve for new x.");
            x += del;
#ifdef MONITOR
    std::cout << "Iter " << iter << ": step= " << del << std::endl;
    std::cout << "Iter " << iter << ": x= " << x << ", err = " << maxNorm(del) << std::endl;
#endif
            if (checkTol(x, del, tol)) return x;

            fx = f(x) - fx; // fx is now the change in f(x) (delta_f)
            if (!fx.isFinite()) throw std::runtime_error("ERROR: Unable to evaluate function.");
            fx -= A*del; // fx holds delta_f - A*delta_x
            A += outer(fx, del)/scDot(del, del);
        } catch (const std::logic_error& ex){
            throw ex;
        } catch (const std::runtime_error& ex){
            throw ex;
        }
    }

    std::string err = "ERROR: Unable to converge within " + std::to_string(maxIter) + " iteration(s).";
    throw std::runtime_error(err);
}


///////////////////////////////////////////////////////////////////////////////
Vector broyden2(std::function<Vector(const Vector&)> f,
                Matrix& B, Vector& x, int maxIter, double tol){
/**  Inputs:
        f              The name of the function for which a root is sought.
        B              The initial matrix, approximation of the inverse Jacobi.
        x              The initial guess at the solution.
        tol            The convergence tolerance (must be > 0).
        maxIter        The maximum number of iterations that can be taken.
    Outputs:
        Vector         The solution Vector.
*/


    if(maxIter < 1) maxIter = 1;
    if(tol <= 0) throw std::invalid_argument("ERROR: Invalid tolerance.");
    int n = x.n();
    if (B.n(0) != n || B.n(1) != n)
        throw std::invalid_argument("ERROR: Incompatible matrix-vector sizes.");

#ifdef MONITOR
    std::cout << "Guess: x= " << x << ", B= " << std::endl << B << std::endl;
#endif

    for(int iter = 1; iter <= maxIter; iter++) {
        Vector fx = f(x);
        if (!fx.isFinite()) throw std::runtime_error("ERROR: Unable to evaluate function.");
#ifdef MONITOR
    std::cout << std::endl << "Iter " << iter << ": fx=" << fx
              << ", B= " << std::endl << B << std::endl;
#endif
        Vector del = -B*fx; // del = -B*f(x_i) = x_i+1 - x_i
        x += del;

#ifdef MONITOR
    std::cout << "Iter " << iter << ": step= " << del << std::endl;
    std::cout << "Iter " << iter << ": x= " << x << ", err = " << maxNorm(del) << std::endl;
#endif
        if (checkTol(x, del, tol)) return x;

        fx = f(x) - fx; // fx = f(x_i+1) - f(x_i);
        if (!fx.isFinite()) throw std::runtime_error("ERROR: Unable to evaluate function.");

        Vector BF = B*fx;
        auto z = scDot(del, BF);
        BF = del-BF;
        B += outer(BF, del) * B / z;
    }

    std::string err = "ERROR: Unable to converge within " + std::to_string(maxIter) + " iteration(s).";
    throw std::runtime_error(err);
}
