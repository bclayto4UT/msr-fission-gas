/******************************************************************************
Vector, Matrix, and Permutation classes, and related functions

For example:
  Vector v(5) gives avector v of 5 doubles, indexed from 0 to 4.
  Matrix m(5,12) gives a matrix of size 5X12 doubles,
                 indexed from 0 to 4 and 0 to 11.
******************************************************************************/

/*****************************************************************************
File: matrix.h
CHANGELOG as of 08/27/2021:
 * New functions Vector::count and Matrix::count, to count the number of times
 a value appears in a Vector/Matrix.

CHANGELOG as of 03/04/2021:
 * Modification of Matrix::Matrix(const Vector&)
 * Update of the Permutation class
 * Addition of a move assignment operator in the Vector and Matrix classes
 * Addition of operator Matrix() and assign functions in Vector class

CHANGELOG as of 03/29/2021:
 * Addition of operator functions +, -, *, and / for Vector and Matrix classes,
   and redefinition of operator functions +=, -=, *=, and /=.
 * matVecMult and matMatMult functions changed to member operator functions for
   more ease coding functions using Matrix arithmetics.
 * Wrapped container type changed to std::vector<double> from double[] for more
   ease coding overloaded functions.
 * Addition of new functions in the Matrix class (implementations in matrix.cpp):
   1. Matrix(const Vector&): typecasting constructor
   2. Vector row(const int) const: returns a row in Matrix
   3. Vector col(const int) const: returns a column in Matrix
   4. Matrix T(bool inplace=false) noexcept: transpose, does not return if inplace
   5. void assign(const Vector&, int, bool colAssign=true): assign a row/column
        by a Vector
   6. void assign(const Matrix&, const int rbegin=0, const int cbegin=0): assign
        a submatrix by another matrix
   7. void resize(const int r, const int c, const int rbegin=0, const int cbegin=0):
        resize current Matrix object
*******************************************************************************/

#ifndef Matrix_Included
#define Matrix_Included

#include <array>
//#include <complex>
#include <iostream>
#include <vector>

//using cdouble = std::complex<double>;

class Vector {
    private:
        int size;
        std::vector<double> array;

    public:
        Vector(int n=0);
        Vector(const Vector&);
        Vector(std::initializer_list<double> il);

        int n(int=0) const { return size; }  // return size of the vector

        double operator() (int i) const { return array[i]; }
         // access to return element of vector
        double& operator() (int i) { return array[i]; }
         // access to modify element of vector

        Vector& operator=(double);         // fill value into vector
        Vector& operator=(double*);        // copy array to vector. NOTE: unsafe
        Vector& operator=(std::vector<double>&);  // copy std::vector into vector. NOTE: safe
        Vector& operator=(const Vector&);  // copy assignment
        Vector& operator=(Vector&&);       // move assignment

        bool operator==(const Vector&) const noexcept;

        Vector operator+(const Vector&) const; // addition
        Vector operator+(const double d) const;
        Vector& operator+=(const Vector& v){ *this = *this + v; return *this;}
        Vector& operator+=(const double d){ *this = *this + d; return *this;}

        Vector operator* (const double) const; // scalar multiplication
        Vector& operator*=(const double d){ *this = *this * d; return *this;}
        double operator*(const Vector& v) const; // scalar dot product

        Vector operator-() const{ return *this * -1;} // Unary minus operator
        Vector operator- (const Vector& v) const{ return *this + -v;}
        Vector operator- (const double d) const{ return *this + -d;}
        Vector& operator-=(const Vector& v){ *this = *this - v; return *this;}
        Vector& operator-=(const double d){ *this = *this - d; return *this;}

        Vector operator/ (const double d) const{ return *this * (1/d);}
        Vector& operator/=(const double d){ *this = *this / d; return *this;}

        //operator Matrix(){ return Matrix(*this);}

        void assign(const Vector&, const int begin=0);
        void resize(const int n, const int begin=0);
        int count(double val, double tol=0) const noexcept;

        bool isFinite() const;
};

class Matrix {
    private:
        std::array<size_t, 2> size;
        std::vector<double> array;

    public:
        Matrix(int n0, int n1);
        Matrix(const Matrix&);
        Matrix(const Vector&); // Typecasting constructor

        int n(int i) const { return size[i]; }
        bool isDiagDominant() const noexcept;
        bool isSymmetric() const noexcept;

        double operator() (int i, int j) const { return array[i + size[0]*j]; }
        // access to return element of matrix
        double& operator() (int i, int j) { return array[i + size[0]*j]; }
        // access to modify element of matrix

        Matrix& operator=(double);         // fill value into matrix
        Matrix& operator=(const Matrix&);  // copy assignment
        Matrix& operator=(Matrix&&);       // move assignment

        Matrix operator+(const Matrix&) const; // addition
        Matrix operator+(const double d) const;
        Matrix& operator+=(const Matrix& M){ *this = *this + M; return *this;}
        Matrix& operator+=(const double d){ *this = *this + d; return *this;}

        Matrix operator*(const double) const; // scalar multiplication
        Vector operator*(const Vector&) const; // matrix-vector multiplication
        Matrix operator*(const Matrix&) const; // matrix-matrix multiplication
        Matrix& operator*=(const double d){ *this = *this * d; return *this;}
        Matrix& operator*=(const Matrix& M){ *this = *this * M; return *this;}

        Matrix operator-() const{ return *this * -1;} // Unary minus operator
        Matrix operator-(const Matrix& M) const{ return *this + -M;}
        Matrix operator-(const double d) const{ return *this + -d;}
        Matrix& operator-=(const Matrix& M){ *this = *this - M; return *this;}
        Matrix& operator-=(const double d){ *this = *this - d; return *this;}

        Matrix operator/(const double d) const{ return *this * (1/d);}
        Matrix& operator/=(const double d){ *this = *this / d; return *this;}

        Vector row(const int) const; // returns a row in Matrix
        Vector col(const int) const; // returns a column in Matrix
        Matrix T(bool inplace=false) noexcept; // transpose; does not return if inplace

        void assign(const Vector&, int, bool colAssign=true);
        void assign(const Matrix&, const int rbegin=0, const int cbegin=0);
        void resize(const int r, const int c, const int rbegin=0, const int cbegin=0);
        int count(double val, double tol=0) const noexcept;
};

class Permutation {
    private:
        int size;
        std::vector<int> array;
        int my_parity;

    public:
        Permutation(int n) { size = n; identity(); };

        int n(int=0) const { return size; }
        int operator() (int i) const { return array[i]; }

        void identity();
        void swap(int i, int j);
        double parity() const { return my_parity; }

        void permute(Vector& b) const;
        void permute(Matrix& M) const;
};

std::ostream& operator<< (std::ostream&, const Vector&); // output vector
std::ostream& operator<< (std::ostream&, const Matrix&); // output matrix by rows
std::ostream& operator<< (std::ostream&, const Permutation&); // output permutation
std::istream& operator>> (std::istream&, Vector&);       // input vector
std::istream& operator>> (std::istream&, Matrix&);       // input matrix by rows

double l2norm(const Vector&);  // L2-norm of the vector
double maxNorm(const Vector&); // L-infinity (max) norm of the vector
double maxNorm(const Matrix&); // L-infinity (max) norm of the matrix

Matrix eye(int); // Identity Matrix (not Permutation)
double scDot(const Vector&, const Vector&); // Scalar dot product of vectors
Vector cross(const Vector&, const Vector&); // Vector cross product of vectors
Matrix outer(const Vector&, const Vector&); // Matrix outer product of vectors

#endif
