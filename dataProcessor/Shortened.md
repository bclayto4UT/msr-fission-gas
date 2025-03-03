# Thermochimica Data Processing Code

## Core Files
- **dataProcessor.h**: Defines functions for processing thermochemical data
- **dataProcessor.cpp**: Implements the processing functions

## Main Functions

### `scaleToVector(const std::string& inFile) -> concMap`
- Parses output files to extract element concentration data over time
- Returns a map of element names to concentration vectors

### `vectToTherm(const concMap& dataMap, const std::string& outFile, std::string& strT, std::string& strP, bool includesSurr=true)`
- Converts concentration data to Thermochimica input files (.F90)
- Configures temperature and pressure conditions

### `textToExcel(const std::string& inFile, const std::string& outFile, std::string& dataType)`
- Extracts specified data from Thermochimica output
- Creates tabulated text files for analysis

### `mergeTherm(const strVect& inFiles, const std::string& outFile)`
- Combines multiple Thermochimica output files
- Sums species concentrations across salt and gas phases

### `decoupleSurr(const concMap& scaleData, const std::string& thermoRes, const std::string& thermoOut, const bool includesSS=false, const bool includesFP=false)`
- Main function that decouples surrogate elements
- Updates files with additional elements and compounds
- Optional parameters for including fluorides of structural metals and fission products

## Helper Functions
- `getIonPair(std::string str)`: Extracts cation and anion from chemical formulas
- `compVector(std::pair<std::string, double> p, std::pair<std::string, double> q)`: Sorting function for comparing pairs
- `del(const Vector& zeta)`: Calculates UF4 reduction to UF3
- `thermoFunc(const Vector& zeta)`: Calculates fluoride amounts

## Dependencies
- Standard libraries: `<algorithm>`, `<iostream>`, `<fstream>`, `<unordered_map>`
- Custom headers: "miscellaneous.h", "iterativeNL.h", "rootFinding.h", "thermoElectroChem.h"

# ThermoElectroChem Code Summary

## Core Files
- **thermoElectroChem.h**: Defines data structures and declarations for thermodynamic calculations
- **thermoElectroChem.cpp**: Implements Gibbs free energy calculations for metal fluorides

## Key Data Structures

### Element Property Maps
- `atomNumMap`: Element symbols → atomic numbers (H to Og)
- `atomWeightMap`: Element symbols → atomic weights in g/mol
- `oxiStateMap`: Element symbols → common oxidation states
- `surrogateMap`: Groups elements by surrogate categories (I, Ca, La, Pu, Th)
- `surrogateMapInv`: Individual elements → their surrogate representative

### Thermodynamic Data
- `heatData`: Stores thermodynamic data for compounds as arrays containing:
  - ΔHf: Enthalpy of formation at 298.15K (kJ/mol)
  - ΔS: Entropy at 298.15K (J/mol·K)
  - Coefficients for heat capacity calculation
  - T_ref: Reference temperature (K)

### Constants
- `R = 8.314`: Universal gas constant (J/mol·K)
- Activity coefficients for specific compounds:
  - `gamma_UF3 = 50`
  - `gamma_UF4 = 0.55`
  - `gamma_Inf_CrF2 = 0.5`
  - `gamma_Inf_FeF2 = 1.6`

## Main Functions

### Core Thermodynamic Calculations
- `calc_H(const thermoArray& data, const double T)`: Calculates enthalpy
- `calc_S(const thermoArray& data, const double T)`: Calculates entropy
- `calc_G(const thermoArray& data, const double T)`: Calculates Gibbs free energy

### Activity Coefficient Function
- `gamma_Inf_NiF2(double xLiF)`: Calculates NiF2 activity coefficient based on LiF mole fraction

### Gibbs Free Energy Functions
Functions calculating Gibbs free energy for various fluoride compounds:
- `G_HF(const double T)`: HF formation from H₂ and F₂
- `G_UF4(const double T)`: UF₄ formation from UF₃ and F₂
- `G_F(const double xUF3, const double xUF4, const double T)`: Fluorine potential based on UF₃/UF₄ ratio

### Metal Fluoride Formation Energies
Multiple functions calculating Gibbs free energies for:
- Chromium fluorides: `G_CrF2`, `G_CrF3`
- Iron fluorides: `G_FeF2`, `G_FeF3`
- Cobalt fluorides: `G_CoF2`, `G_CoF3`
- Nickel fluoride: `G_NiF2`
- Manganese fluoride: `G_MnF2`
- Niobium fluorides: `G_NbF`, `G_NbF2`, `G_NbF3`, `G_NbF4`, `G_NbF5`, `G_Nb2F10`, `G_Nb3F15`
- Molybdenum fluorides: `G_MoF4`, `G_MoF5`, `G_MoF6`

## Dependencies
- Standard C++ math library
- "miscellaneous.h" (likely contains utility functions and the `strVect` type)

# Miscellaneous Utility Code

## Files
- **miscellaneous.h**: Defines utility functions and mathematical constants
- **miscellaneous.cpp**: Implements the utility functions

## Key Components

### Mathematical Constants
- Mathematical constants like e, phi, and pi

### String Processing Functions

- `elementSymb(const std::string& str)`
  - Converts element names to proper case (first letter uppercase, rest lowercase)
  - Mimics Python's string casing functions

- `strToVect(std::string& data, char deli=' ')`
  - Splits a string into a vector of strings using the specified delimiter
  - Handles both spaces and tabs

- `strToVectDouble(const std::string& str)`
  - Converts a string to a vector of doubles
  - Supports space, colon, and comma-separated formats

### String Validation Functions

- `containsNumber(const std::string& str) noexcept`
  - Checks if a string contains any numeric digits

- `isNumeric(const std::string& str) noexcept`
  - Determines if a string can be converted to a number
  - Handles decimals, scientific notation, and signs

- `convertibleNum(std::string& str)`
  - Removes leading characters until the string can be converted to a double

### Mathematical Utility Functions

- `template<typename T> inline constexpr int sgn(T x) noexcept`
  - Returns the sign of a number as an integer (-1, 0, or 1)

- `inline int factorial(int n)`
  - Calculates the factorial of a number using recursion

- `inline int ordMag(double x)`
  - Calculates the order of magnitude of a number

## Dependencies
- C++ Standard Library: `<string>`, `<vector>`, `<cmath>`, `<stdexcept>`
- Standard character functions: `toupper()`, `tolower()`, `isdigit()`
- Conversion functions: `std::stod()`
- Math functions: `exp()`, `sqrt()`, `atan()`, `log10()`, `fabs()`, `ceil()`

# Numerical Methods Code

## Linear Algebra (GaussElim)

### Files
- **gaussElim.h**: Declares matrix operation functions
- **gaussElim.cpp**: Implements Gaussian elimination algorithms

### Key Functions
- **LU Factorization**
  - `luFactorize(Matrix& a, Permutation& p, bool inplace=false)`: Performs LU factorization

- **Linear System Solving**
  - `luSolve(const Matrix& a, const Permutation& p, Vector& x, bool inplace=false)`: Solves using factorized matrix
  - `solve(Matrix& a, Permutation& p, Vector& x, bool inplace=false)`: Factorizes and solves in one step

- **Matrix Operations**
  - `detFactoredMatrix(const Matrix& a, const Permutation& p)`: Calculates determinant from factorized matrix
  - `det(Matrix& a)`: Calculates determinant
  - `invFactoredMatrix(const Matrix& a, const Permutation& p)`: Calculates inverse from factorized matrix
  - `inv(Matrix& a)`: Calculates matrix inverse

- **Helper Functions**
  - `swapRows(Matrix& a, int i, int j)`: Swaps matrix rows
  - `swap(Vector& v, int i, int j)`: Swaps vector elements
  - `rowReplacement(Matrix& a, int i, int j)`: Performs row replacement operation

## Nonlinear Systems (IterativeNL)

### Files
- **iterativeNL.h**: Declares nonlinear system solvers
- **iterativeNL.cpp**: Implements iterative methods

### Key Functions
- `fixedpt(std::function<Vector(const Vector&)> g, Vector& x, int maxIter, double tol)`: Fixed point iteration
- `newton(std::function<Vector(const Vector&)> f, std::function<Matrix(const Vector&)> df, Vector& x, int maxIter, double tol)`: Newton with analytical Jacobian
- `newton(std::function<Vector(const Vector&)> f, Vector& x, int maxIter, double tol)`: Newton with numerical Jacobian
- `broyden1(std::function<Vector(const Vector&)> f, Matrix& A, Vector& x, int maxIter, double tol)`: Broyden's first method
- `broyden2(std::function<Vector(const Vector&)> f, Matrix& B, Vector& x, int maxIter, double tol)`: Broyden's second method
- `checkTol(const Vector& x, const Vector& dx, const double tol)`: Convergence check

## Root Finding (RootFinding)

### Files
- **rootFinding.h**: Declares scalar root-finding methods
- **rootFinding.cpp**: Implements root-finding algorithms

### Key Functions
- `bisection(std::function<double(double)> f, double x0, double x1, int maxIter=20, double tol=1e-6)`: Bisection method
- `fixedpt(std::function<double(double)> f, double x0, int maxIter=20, double tol=1e-6)`: Fixed point iteration
- `newton(std::function<double(double)> f, std::function<double(double)> df, double x0, int maxIter=20, double tol=1e-6)`: Newton with analytical derivative
- `newton(std::function<double(double)> f, double x0, int maxIter=20, double tol=1e-6)`: Newton with numerical derivative
- `secant(std::function<double(double)> f, double x0, double x1, int maxIter=20, double tol=1e-6)`: Secant method

## Main Program

### File
- **main.cpp**: Contains menu-driven interface for data processing operations

### Features
- Menu-based interface for:
  - Converting SCALE output to vector format
  - Converting vector data to Thermochimica input
  - Extracting data from Thermochimica output
  - Combining multiple Thermochimica outputs
  - Decoupling surrogate elements

### Dependencies
- **dataProcessor.h**: For all data processing functions

## Common Dependencies
- **matrix.h**: Matrix and Vector classes
- **calculus.h**: Numerical differentiation
- **miscellaneous.h**: Utility functions
- Standard C++ libraries: `<functional>`, `<cmath>`

# LinearAlgebra Library

## Overview
This library defines Vector, Matrix, and Permutation classes for mathematical operations, providing comprehensive functionality for linear algebra operations.

## Classes

### Vector
A dynamic vector of double values with standard operations.

**Constructors:**
- `Vector(int n=0)`: Creates n-sized vector (zeros)
- `Vector(const Vector&)`: Copy constructor
- `Vector(std::initializer_list<double>)`: From initializer list

**Key Methods:**
- `int n(int=0) const`: Returns vector size
- `double operator()(int i) const`/`double& operator()(int i)`: Element access
- Assignment operators for various data sources
- Arithmetic operators: `+`, `+=`, `*` (scalar mult), `*=`, `*` (dot product), `-` (negation)
- `void resize(const int n, const int begin=0)`: Resize preserving values
- `int count(double val, double tol=0) const`: Count matching elements
- `bool isFinite() const`: Check if all elements are finite

### Matrix
A dynamic 2D matrix of double values with standard operations.

**Constructors:**
- `Matrix(int n0, int n1)`: Creates n0×n1 matrix (zeros)
- `Matrix(const Matrix&)`: Copy constructor
- `Matrix(const Vector&)`: Creates column matrix from Vector

**Key Methods:**
- `int n(int i) const`: Returns row (i=0) or column (i=1) count
- `bool isDiagDominant()/isSymmetric() const`: Property checks
- Access operators for reading/writing elements
- `Vector row(const int)/col(const int) const`: Extract row/column
- `Matrix T(bool inplace=false)`: Transpose matrix
- `void assign(...)`: Assign values to row/column or submatrix
- `void resize(const int r, const int c, ...)`: Resize preserving values

### Permutation
Represents a permutation of integers with parity tracking.

**Constructor:**
- `Permutation(int n)`: Creates identity permutation of size n

**Key Methods:**
- `int n(int=0) const`: Returns permutation size
- `int operator()(int i) const`: Element access
- `void identity()`: Reset to identity permutation
- `void swap(int i, int j)`: Swap elements and update parity
- `double parity() const`: Returns parity (+1 for even, -1 for odd)
- `void permute(Vector&/Matrix&) const`: Apply permutation

## Global Functions

**I/O Operations:**
- Stream operators (`<<`, `>>`) for Vector, Matrix, and Permutation

**Mathematical Operations:**
- `double l2norm(const Vector&)`: Euclidean norm
- `double maxNorm(const Vector&/Matrix&)`: Maximum norm
- `Matrix eye(int)`: Identity matrix generator
- `double scDot(const Vector&, const Vector&)`: Scalar dot product
- `Vector cross(const Vector&, const Vector&)`: Cross product
- `Matrix outer(const Vector&, const Vector&)`: Outer product

I'll create a more concise version of the calculus library documentation by removing redundancies:

## Calculus Library (`calculus.h`)

This header defines numerical calculus functions for both scalar and vector functions.

### Function Categories:

1. **First Derivatives**
   - `deriv()`: Calculates derivative of scalar or vector function at point x
   - `pderiv()`: Calculates partial derivative of multivariate scalar or vector function

2. **Differential Operators**
   - `grad()`: Calculates gradient of scalar function
   - `div()`: Calculates divergence of vector field
   - `curl()`: Calculates curl of vector field
   - `jacobian()`: Calculates Jacobian matrix of vector function

3. **Second Derivatives**
   - `deriv2()`: Calculates second derivative of scalar or vector function
   - `pderiv2()`: Calculates second partial derivative of multivariate function
   - `laplace()`: Calculates Laplacian of scalar function or vector field
   - `hessian()`: Calculates Hessian matrix of scalar function

4. **Integrals**
   - `intClose()`: Definite integral using closed method
   - `intOpen()`: Definite integral using open method
   - `intGauss()`: Definite integral using Gaussian quadrature
   - `integral()`: Definite integral using most appropriate method

### Implementation Notes:
- Uses finite difference methods for derivatives
- Uses Vector and Matrix classes from "matrix.h"
- Integral calculation functions are declared but commented out in implementation file
- Includes error handling for invalid inputs and numerical issues
