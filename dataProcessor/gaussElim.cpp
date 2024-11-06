#include "gaussElim.h"
#include <cmath>

static void swapRows(Matrix& a, int i, int j) {
    // E_i <--> E_j
    for(int k=0; k<a.n(1); k++) {
        auto aa = a(i,k);
        a(i,k) = a(j,k);
        a(j,k) = aa;
    }
}

static void swap(Vector& v, int i, int j) {
    auto vv = v(i);
    v(i) = v(j);
    v(j) = vv;
}

static void rowReplacement(Matrix& a, int i, int j) {
  // E_j --> E_j - (a(j,i)/a(i,i))E_i
    auto f = a(j,i)/a(i,i);
    a(j,i) = f;
    for(int k=i+1; k<a.n(1); k++) {
        a(j,k) -= f*a(i,k);
    }
}

Matrix luFactorize(Matrix& a, Permutation& p, bool inplace) {
    if (!inplace){
        Matrix M(a);
        try{
            luFactorize(M, p, 1);
            return M;
        } catch (const std::logic_error& ex){
            throw ex;
        } catch (const std::runtime_error& ex){
            throw ex;
        }
    }

    int n = a.n(0);
    int i,j,k;

    if(a.n(1) != n || p.n() != n)
        throw std::invalid_argument("ERROR: Incompatible sizes.");

    // Set up permutation
    p.identity();

    // Determine scale factors for scaled partial pivoting
    Vector s(n);
    for(int i=0; i<n; i++) {
        s(i) = abs(a(i,0));
        for(int j=1; j<n; j++) {
            if( s(i) < abs(a(i,j)) ) s(i) = abs(a(i,j));
        }
    }

    // Loop on Columns
    for(j=0; j<n; j++) {

    // Get nonzero pivot (use scaled partial pivoting)
        double pivot = fabs(a(j,j))/s(j);
        i = j;
        for(k=j+1; k<n; k++) {
            double q = fabs(a(k,j))/s(k);
            if(q > pivot) {
                pivot = q;
                i = k;
            }
        }
        if(pivot == 0) throw std::runtime_error("ERROR: Zero pivot encountered.");
        if(i != j) {
            swapRows(a,i,j);
            p.swap(i,j);
            swap(s,i,j);
        }

    // Loop on rows
        for(i=j+1; i<n; i++) rowReplacement(a,j,i);
    }
    return a;
}

Vector luSolve(const Matrix& a, const Permutation& p, Vector& x, bool inplace) {
    if (!inplace){
        Matrix M(a);
        Vector b(x);
        try{
            luSolve(M, p, b, 1);
            return b;
        } catch (const std::invalid_argument& ex){
            throw ex;
        }
    }

    int n = a.n(0);
    int i,j;

    if(a.n(1) != n || p.n() != n || x.n() != n)
        throw std::invalid_argument("ERROR: Incompatible sizes.");

    // Apply permutation to x
    p.permute(x);

    // FORWARD SUBSTITUTION

    // Loop on columns and rows of a
    for(j=0; j<n; j++) {
        for(i=j+1; i<n; i++) {
            x(i) -= a(i,j)*x(j);
        }
    }

    // BACKWARD SUBSTITUTION

    for(i=n-1; i>=0; i--) {
        for(j=i+1; j<n; j++) {
            x(i) -= a(i,j)*x(j);
        }
        x(i) /= a(i,i);
    }
    return x;
}

Vector solve(Matrix& a, Permutation& p, Vector& x, bool inplace) {
    if (!inplace){
        Matrix M(a);
        Vector b(x);
        try{
            solve(M, p, b, 1);
            return b;
        } catch (const std::logic_error& ex){
            throw ex;
        } catch (const std::runtime_error& ex){
            throw ex;
        }
    }

    try{
        luFactorize(a, p, 1);
        luSolve(a, p, x, 1);
        return x;
    } catch (const std::logic_error& ex){
        throw ex;
    } catch (const std::runtime_error& ex){
        throw ex;
    }
}

double detFactoredMatrix(const Matrix& a, const Permutation& p) noexcept{
    int n = a.n(0);
    std::cout << a.n(0) << " " << a.n(1) << " " << p.n() << std::endl;
    if(a.n(1) != n || p.n() != n) return NAN;

    double detA = a(0, 0);
    for(int i = 1; i < n; i++) {
        detA *= a(i, i);
    }
    return p.parity() * detA;
}

double det(Matrix& a) noexcept{
    int n = a.n(0);
    if(a.n(1) != n) return NAN;

    Permutation p(n);
    Matrix M = luFactorize(a, p);
    return detFactoredMatrix(M, p);
}

Matrix invFactoredMatrix(const Matrix& a, const Permutation& p) {
    int n = a.n(0);
    if (a.n(1) != n || p.n() != n)
        throw std::invalid_argument("ERROR: Incompatible sizes.");

    Matrix invA(n,n);
    for(int j = 0; j < n; j++){
        // initialize j'th column of identity matrix, e_j
        Vector y(n); y(j)=1;

        // solve matrix equation where b=e_j
        try{
            luSolve(a, p, y, 1);
        } catch (const std::exception& ex){
            throw ex;
        }
        invA.assign(y, j);
    }
    return invA;
}

Matrix inv(Matrix& a){
    int n = a.n(0);
    if (a.n(1) != n) throw std::invalid_argument("ERROR: Incompatible sizes.");

    Permutation p(n);
    try{
        Matrix M = luFactorize(a, p);
        return invFactoredMatrix(M, p);
    } catch (const std::logic_error& ex){
        throw ex;
    } catch (const std::runtime_error& ex){
        throw ex;
    }
}
