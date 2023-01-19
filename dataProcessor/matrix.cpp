/****************************************************************************
File: matrix.cpp
****************************************************************************/
#include "math.h"
#include "matrix.h"
// Vector /////////////////////////////////////////////////////////////////////

Vector::Vector(int n){
    size = n;
    for(int i = 0;i < n; i++) array.push_back(0);
}

Vector::Vector(const Vector& v) {
    size = v.size;
    for(int i=0; i<size; i++) array.push_back(v(i));
}

Vector::Vector(std::initializer_list<double> il){
    array = std::move(il);
    size = array.size();
}

Vector& Vector::operator=(double x) {
    for(int i=0; i<size; i++) array[i] = x;
    return *this;
}

// NOTE: this is an unsafe assignment, as the size of the array that will be copied cannot be checked.
Vector& Vector::operator=(double* a){
    for(int i=0; i<size; i++) array[i] = a[i];
    return *this;
}

// NOTE: this is safe
Vector& Vector::operator=(std::vector<double>& a){
    if(size != a.size()) throw std::invalid_argument("ERROR: Incompatible sizes.");
    array = a;
}

Vector& Vector::operator=(const Vector& v) {
    if(size != v.size) throw std::invalid_argument("ERROR: Incompatible sizes.");
    array = v.array;
    return *this;
}

Vector& Vector::operator=(Vector&& v){
    if(size != v.size) throw std::invalid_argument("ERROR: Incompatible sizes.");
    array = std::move(v.array);
    return *this;
}

bool Vector::operator==(const Vector& v) const noexcept{
    if (size != v.size) return false;
    for (int i = 0; i < size; i++){
        if (array[i] != v.array[i]) return false;
    }
    return true;
}

Vector Vector::operator+(const Vector& v) const{
    if(size != v.size) throw std::invalid_argument("ERROR: Incompatible sizes.");
    Vector res(size);
    for(int i=0; i<size; i++) res(i) = array[i] + v.array[i];
    return res;
}

Vector Vector::operator+(const double d) const{
    Vector res(size);
    for(int i=0; i<size; i++) res(i) = array[i] + d;
    return res;
}

Vector Vector::operator*(const double x) const{
    Vector res(size);
    for(int i=0; i<size; i++) res(i) = array[i] * x;
    return res;
}

double Vector::operator*(const Vector& v) const {
  int n = size;
  if(size > v.size) n = v.size;

  double dot = 0;
  for(int i=0; i<n; i++) dot += array[i]*v.array[i];
  return dot;
}

void Vector::assign(const Vector& v, const int begin){
/**
    Input:
        v: a Vector containing values to modify
        begin: the index where modification begins (DEFAULT=0).
    Output: none
            This Vector object is modified beginning at the specified index.
            If an out-of-bound index occurs, a std::exception is throw.
**/
    int n = v.n();
    if (n+begin > size) throw std::out_of_range("ERROR: Index out of range.");
    for (int i = 0; i < n; i++) array[i+begin] = v(i);
}

void Vector::resize(const int n, const int begin) {
/**
    Input: n: the number of entries in the resized Vector.
           begin: the index at which the resized Vector begins.
    Output: none
            This Vector object is resized starting at begin.
            If an index is out of bound, a 0. entry is filled in instead.
**/
    if(n<1 || begin<0)
        throw std::length_error("ERROR: Invalid indices.");
    if (n != size || begin != 0) {
        Vector b(n);
        for (int i=0; i < n; i++){
            b(i) = (i+begin < size) ? array[i+begin] : 0;
        }
        array = b.array;
        size = n;
        *this = b;
    }
}

int Vector::count(double val, double tol) const noexcept{
/**
    Input: val: the value to be counted
           tol: the acceptable deviation from val (DEFAULT=0)
    Output: the number of appearances of val += tol
**/
    int c = 0;
    if (isfinite(val)){
        for (int i = 0; i < size; i++){
            if (abs(array[i]-val) <= tol) c++;
        }
    } else{
        for (int i = 0; i < size; i++){
            if (!isfinite(array[i])) c++;
        }
    }
    return c;
}

bool Vector::isFinite() const{
    for (int i = 0; i < size; i++){
        if (!isfinite(array[i])) return false;
    }
    return true;
}

// Matrix /////////////////////////////////////////////////////////////////////

Matrix::Matrix(int n0, int n1){
    size[0] = n0;
    size[1] = n1;
    for(int i = 0; i < n0*n1; i++) array.push_back(0);
}

Matrix::Matrix(const Matrix& m) {
    size[0] = m.size[0];
    size[1] = m.size[1];
    array = m.array;
}

Matrix::Matrix(const Vector& v) {
    size[0] = v.n();
    size[1] = 1;
    for (int i = 0; i < size[0]; i++) array.push_back(v(i));
}

bool Matrix::isDiagDominant() const noexcept{
    if (size[0] != size[1]) return false;
    for (int i=0; i<size[0]; i++) {
        for (int j=0; j <size[0]; j++){
            if ((j*size[0] + j) == (i*size[0])){continue;}
            if (array[size[0]*i + i] < array[j*size[0] + i]) { return false;}
        }
    }
    return true;
}

bool Matrix::isSymmetric() const noexcept{
    if (size[0] != size[1]) return false;
    for (int i=0; i<size[0]; i++) {
        for (int j=i+1; j <size[0]; j++){
            if ((*this)(i,j) != (*this)(j, i)) return false;
        }
    }
    return true;
}

Matrix& Matrix::operator=(double x) {
    for(int i=0; i<size[0]*size[1]; i++) array[i] = x;
    return *this;
}

Matrix& Matrix::operator=(const Matrix& m) {
    if(size[0] != m.size[0] || size[1] != m.size[1]){
        throw std::invalid_argument("ERROR: Incompatible sizes.");
    }
    array = m.array;
    return *this;
}

Matrix& Matrix::operator=(Matrix&& m) {
    if(size[0] != m.size[0] || size[1] != m.size[1]){
        throw std::invalid_argument("ERROR: Incompatible sizes.");
    }
    array = std::move(m.array);
    return *this;
}

Matrix Matrix::operator+(const Matrix& M) const{
    if(size[0] != M.size[0] || size[1] != M.size[1]){
        throw std::invalid_argument("ERROR: Incompatible sizes.");
    }
    Matrix res(size[0], size[1]);
    for(int i=0; i<size[0]*size[1]; i++) res.array[i] = this->array[i] + M.array[i];
    return res;
}

Matrix Matrix::operator+(const double d) const{
    Matrix res(size[0], size[1]);
    for(int i=0; i<size[0]*size[1]; i++) res.array[i] = this->array[i] + d;
    return res;
}

Matrix Matrix::operator*(const double d) const{
    Matrix res(size[0], size[1]);
    for(int i=0; i<size[0]*size[1]; i++) res.array[i] = this->array[i] * d;
    return res;
}

Vector Matrix::operator*(const Vector& x) const{
    if(this->size[1] != x.n()) throw std::invalid_argument("ERROR: Incompatible sizes.");

    Vector y(size[0]);
    for(int i=0; i<y.n(); i++) {
        double sum = 0;
        for(int j=0; j<x.n(); j++) {
            sum += (*this)(i,j)*x(j);
        }
        y(i) = sum;
    }
    return y;
}

Matrix Matrix::operator*(const Matrix& M) const{
    if(this->size[1] != M.size[0]) throw std::invalid_argument("ERROR: Incompatible sizes.");

    Matrix C(this->size[0], M.size[1]);
    for(int j=0; j<C.n(1); j++){
        for(int i=0; i<C.n(0); i++){
            double sum = 0;
            for(int k=0; k<this->size[1]; k++) {
                sum += (*this)(i,k)*M(k,j);
            }
            C(i,j) = sum;
        }
    }
    return C;
}

Vector Matrix::row(int i) const{
/**
    Input: the index of the desired row
    Output: A Vector containing entries in that row.
**/
    if (i < 0 || i >= size[0]){
        throw std::out_of_range("ERROR: Index out of range.");
    }
    Vector r(size[1]);
    for (int c = 0; c < size[1]; c++) r(c) = (*this)(i,c);
    return r;
}

Vector Matrix::col(int i) const{
/**
    Input: the index of the desired column
    Output: A Vector containing entries in that column.
**/
    if (i < 0 || i >= size[1]){
        throw std::out_of_range("ERROR: Index out of range.");
    }

    Vector c(size[0]);
    for (int r = 0; r < size[0]; r++) c(r) = (*this)(r,i);
    return c;
}

Matrix Matrix::T(bool inplace) noexcept{
/**
    Input: a bool indicating whether this Matrix object is transposed
    Output: if inplace, no return (this Matrix object is transposed);
            else, a Matrix that is the transpose of this Matrix object.
**/
    int m = size[1], n = size[0];
    Matrix M(m, n);
    for (int i = 0; i < n; i++){
        for (int j = 0; j < m; j++){
            M(j,i) = (*this)(i,j);
        }
    }
    if (inplace){
        size[0] = m, size[1] = n;
        (*this) = M;
    } return M;
}

void Matrix::assign(const Vector& v, int pos, bool colAssign){
/**
    Input:
        v: a Vector containing values to modify
        pos: the specified position to modified
        colAssign: whether a row or column is being modified (DEFAULT=TRUE)
    Output: none
            This Matrix object is modified in the specified row/column.
            If an out-of-bound index occurs, a std::exception is throw.
**/
    if (colAssign){
        if (pos < 0 || pos >= size[1] || v.n() > size[0]){
            throw std::out_of_range("ERROR: Index out of range.");
        }
        for (int i = 0; i < v.n(); i++) (*this)(i, pos) = v(i);

    } else{
        if (pos < 0 || pos >= size[0] || v.n() > size[1]){
            throw std::out_of_range("ERROR: Index out of range.");
        }
        for (int i = 0; i < v.n(); i++) (*this)(pos, i) = v(i);
    }
}

void Matrix::assign(const Matrix& M, const int rbegin, const int cbegin){
/**
    Input:
        M: a Matrix containing values to modify
        rbegin, cbegin: the row and column indices where modification begins (DEFAULT=0).
    Output: none
            This Matrix object is modified beginning at the specified row/column.
            If an out-of-bound index occurs, a std::exception is throw.
**/
    int m = M.n(0), n = M.n(1);
    if (m+rbegin > size[0] || n+cbegin > size[1]){
        throw std::out_of_range("ERROR: Index out of range.");
    }
    for (int i = 0; i < m; i++){
        for (int j = 0; j < n; j++){
            (*this)(i+rbegin, j+cbegin) = M(i,j);
        }
    }
}

void Matrix::resize(const int r, const int c, const int rbegin, const int cbegin){
/**
    Input: r, c: the number of rows and columns in the resized Matrix.
           rbegin, cbegin: the row and column indices where the resized Matrix begins (DEFAULT=0).
    Output: none
            This Matrix object is resized starting at (rbegin, cbegin).
            If an index isvout of bound, a 0. entry is filled in instead.
**/
    if (r < 1 || c < 1 || rbegin < 0 || cbegin < 0){
        throw std::length_error("ERROR: Invalid indices.");
    }
    if (r != size[0] || c != size[1] || rbegin != 0 || cbegin != 0){
        Matrix M(r, c);
        for (int i = 0; i < r; i++){
            for (int j = 0; j < c; j++){
                M(i,j) = (i+rbegin < size[0] && j+cbegin < size[1]) ? (*this)(i+rbegin,j+cbegin) : 0.;
            }
        }
        size[0] = r; size[1] = c;
        array = M.array;
    }
}

int Matrix::count(double val, double tol) const noexcept{
/**
    Input: val: the value to be counted
           tol: the acceptable deviation from val (DEFAULT=0)
    Output: the number of appearances of val += tol
**/
    int c = 0;
    if (isfinite(val)){
        for (int i = 0; i < size[0]; i++){
            for (int j = 0; j < size[0]; j++){
                if (abs((*this)(i,j)-val) <= tol) c++;
            }
        }
    } else{
        for (int i = 0; i < size[0]; i++){
            for (int j = 0; j < size[0]; j++){
                if (!isfinite((*this)(i,j))) c++;
            }
        }
    }
    return c;
}

// Permutation ////////////////////////////////////////////////////////////////

void Permutation::identity() {
    for(int i=0; i<size; i++) array.push_back(i);
    my_parity = 1;
}

void Permutation::swap(int i, int j) {
    int k = array[i];
    array[i] = array[j];
    array[j] = k;
    my_parity *= -1;
}

void Permutation::permute(Vector& b) const {
    if(size > b.n()) return;

    Vector c(b.n());

    for(int i=0; i<size; i++) c(i) = b(array[i]);
    for(int i=0; i<size; i++) b(i) = c(i);
}

void Permutation::permute(Matrix& A) const{
    /** Added function to swap the rows in a matrix. **/
    if (size > A.n(0) || size > A.n(1)) return;
    Matrix X(A.n(0), A.n(1));
    for (int r = 0; r < size; r++){
        for (int c = 0; c < size; c++){
            X(r, c) = A(array[r], c);
        }
    }
    for (int r = 0; r < size; r++){
        for (int c = 0; c < size; c++){
            A(r, c) = X(r, c);
        }
    }
}

// Miscellaneous //////////////////////////////////////////////////////////////

std::ostream& operator<< (std::ostream& os, const Vector& v) {
    for(int i=0; i<v.n(); i++){
        if (i%10 == 0 && i>0) os << "\n";
        os << v(i) << " ";
    }
    os << "\n";
    return os;
}

std::ostream& operator<< (std::ostream& os, const Matrix& m) {
    bool flag = false;
    for(int i=0; i<m.n(0); i++) {
        for(int j=0; j<m.n(1); j++){
            if (j%10 == 0 && j>0){
                os << "\n";
                flag = true;
            }
            os << m(i,j) << " ";
        }
        os << "\n";
        if (m.n(1) > 10) os << "\n";
    }
    return os;
}

std::ostream& operator<< (std::ostream& os, const Permutation& p) {
    for(int i=0; i<p.n(); i++){
        if (i%10 == 0 && i>0) os << "\n";
        os << p(i) << "  ";
    }
    os << "\n";
    return os;
}
std::istream& operator>> (std::istream& is, Vector& v) {
  for(int i=0; i<v.n(); i++) is >> v(i);
  return is;
}

std::istream& operator>> (std::istream& is, Matrix& m) {
  for(int i=0; i<m.n(0); i++)
  for(int j=0; j<m.n(1); j++) {
    is >> m(i,j);
  }
  return is;
}

double l2norm(const Vector& v) {
  double norm = 0;
  for(int i=0; i<v.n(); i++) norm += v(i)*v(i);
  return sqrt(norm);
}

double maxNorm(const Vector& v) {
  double norm = 0;
  for(int i=0; i<v.n(); i++) {
    double a = fabs(v(i));
    if(norm < a) norm = a;
  }
  return norm;
}

double maxNorm(const Matrix& m) {
  double norm = 0;
  for(int i=0; i<m.n(0); i++) {
    double sum=0;
    for(int j=0; j<m.n(1); j++) sum += fabs(m(i,j));
    if(norm < sum) norm = sum;
  }
  return norm;
}

Matrix eye(int n){
    if (n < 1) throw std::domain_error("Invalid matrix size");
    Matrix I(n,n);
    for (int i = 0; i < n; i++){ I(i,i) = 1;}
    return I;
}

double scDot(const Vector& v1, const Vector& v2) { return v1*v2; }

Vector cross(const Vector& v1, const Vector& v2){
    if (v1.n() > 3 || v2.n() > 3) throw std::invalid_argument("ERROR: Incompatible sizes.");
    Vector w1(v1.n()), w2(v2.n()), v3(3);
    w1.resize(3); w2.resize(3);
    v3(0) = w1(2)*w2(3) - w1(3)*w2(2);
    v3(1) = w1(3)*w2(1) - w1(1)*w2(3);
    v3(2) = w1(1)*w2(2) - w1(2)*w2(1);
    return v3;
}

Matrix outer(const Vector& v1, const Vector& v2){
    int n = v1.n();
    if (v2.n() != n) throw std::invalid_argument("ERROR: Incompatible sizes.");
    Matrix m(n,n);
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            m(i,j) = v1(i)*v2(j);
        }
    }
    return m;
}
