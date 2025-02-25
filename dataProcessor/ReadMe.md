# DataProcessor
_________________
# dataProcessor.h
# Description
This header file defines functions for processing data related to Thermochimica, which appears to be a thermochemical equilibrium calculation software. The main functionality focuses on decoupling surrogate elements and updating Thermochimica result files to include additional elements and compounds.

 ## Functions Defined
```c++
concMap scaleToVector(const std::string&)
```
* **Parameters:** A string (likely a file path)
* **Returns:** A concentration map (concMap, which is an unordered map of strings to vector of doubles)
* **Summary::** Converts data from a file into a concentration map structure.

  
```c++
void vectToTherm(const concMap&, const std::string&, std::string&, std::string&, bool=true)
```

* **Parameters**: A concentration map, a string (likely an input file path), two string references (likely output paths), and an optional boolean (defaulted to true)
* **Returns:** void
* **Summary::** Converts a concentration vector to Thermochimica format, storing results in the output string paths.

```c++
void textToExcel(const std::string&, const std::string&, std::string&)
```

* **Parameters:** Two strings (likely file paths) and a string reference (likely an output path)
* **Returns:** void
* **Summary::** Converts text data to Excel format, using the provided paths.

```c++
void mergeTherm(const strVect&, const std::string&)
```

* **Parameters:** A vector of strings and a string (likely file paths)
* **Returns:** void
* **Summary::** Merges multiple Thermochimica result files.

```c++
void decoupleSurr(const concMap&, const std::string&, const std::string&, const bool=false, const bool=false)
```

* **Parameters:** A concentration map, two strings (likely file paths), and two optional boolean parameters (both defaulted to false)
* **Returns:** void
* **Summary::** Main function that decouples surrogate elements and updates Thermochimica files with additional elements and compounds.

## External Function Calls
The file includes the following external dependencies:

```miscellaneous.h```: Likely provides the strVect type used in the code
```Standard C++ library's <unordered_map>```: Provides the unordered_map data structure used to define the concMap type

# dataProcessor.cpp

## Description
This C++ code appears to be part of a thermochemical data processing tool. It handles concentration data from simulation outputs, converts them to Thermochimica input format, and extracts specific data from Thermochimica output for analysis and visualization.

 ## Functions Defined
```c++
scaleToVector(const std::string& inFile) -> concMap
```
* **Parameters:** inFile (string) - path to input file
* **Returns:** concMap - map of element names to concentration vectors
* **Summary::** Parses an output file to extract concentration data for chemical elements across multiple time intervals.

```c++
vectToTherm(const concMap& dataMap, const std::string& outFile, std::string& strT, std::string& strP, bool includesSurr)
```
* **Parameters:**

```dataMap (concMap)``` - concentration map from scaleToVector

```outFile (string)``` - output file name

```strT (string)``` - system temperature in K

```strP (string)``` - system pressure in atm

```includesSurr (bool)``` - whether to include surrogate elements


* **Returns:** None (creates output files)
* **Summary::** Converts concentration data to Thermochimica input file format (.F90) with specified temperature and pressure conditions.

```c++
textToExcel(const std::string& inFile, const std::string& outFile, std::string& dataType)
```
* **Parameters:**

```inFile (string)``` - input file with Thermochimica output

```outFile (string)``` - output file name

```dataType (string)``` - types of data to extract


* **Returns:** None (creates tabulated text file)
* **Summary::** Extracts specified data types from Thermochimica output and formats them into a tabulated text file for analysis.

```c++
mergeTherm(const strVect& inFiles, const std::string& outFile)
```
* **Parameters:**

```inFiles (strVect)``` - vector of input file names containing Thermochimica output

```outFile (string)``` - output file name


Returns: None (creates a merged output file)
Summary: Combines multiple Thermochimica output files by summing species concentrations across salt and gas phases.
```c++
decoupleSurr(const concMap& scaleData, const std::string& thermoRes, const std::string& thermoOut, const bool includesSS, const bool includesFP)

```
* **Parameters:**

```scaleData (const concMap&):``` Original concentration map obtained using the scaleToVector function

```thermoRes (const std::string&):``` Filename containing Thermochimica results with surrogates

```thermoOut (const std::string&):``` Output filename

```includesSS (const bool):``` Whether to consider fluorides of metals in reactor structure (default: false)

```includesFP (const bool):``` Whether to consider fluorides of fission products (H, Nb) (default: false)


* **Purpose:** Processes thermochemical data to decouple surrogate elements into their original counterparts

### Lambda Functions Defined
```getIonPair(std::string str)```
* **Parameters:** str (std::string)
* **Purpose:** Extracts cation and anion from a chemical formula

```compVector(std::pair<std::string, double> p, std::pair<std::string, double> q)```
* **Parameters:** str (std::string)
* **Purpose:** Sorting function that compares pairs based on their second element (value)

```del(const Vector& zeta)```
* **Parameters:** str (std::string)
* **Purpose:** Calculates the total amount of UF4 being reduced to UF3

```thermoFunc(const Vector& zeta)```
* **Parameters:** zeta (const Vector&)
* **Purpose:** Calculates the amount of fluorides from structural metals and/or fission products

## External Functions Called

```getline``` - From <iostream>, reads a line from a stream

```stod``` - From <string>, converts string to double

```strToVect``` - Likely defined in "dataProcessor.h"

```strToVectDouble``` - Likely defined in "dataProcessor.h"

```surrogateMapInv.at``` - References a map likely defined in "thermoElectroChem.h"

```atomNumMap.at``` - References a map likely defined in "thermoElectroChem.h"

```containsNumber``` - Likely defined in "dataProcessor.h"

```convertibleNum``` - Likely defined in "dataProcessor.h"

```newton``` likely defined in "rootFinding.h": Newton's method for solving nonlinear equations

Functions prefixed with G_ and gamma_Inf_ likely defined in "thermoElectroChem.h": Calculate thermodynamic properties

Various standard library functions from <algorithm>, <iostream>, and <fstream>


## Header Files Included

### Standard libraries:

```<algorithm>``` - For standard algorithms

```<iostream> ```- For input/output operations

```<fstream>``` - For file operations


### Custom headers:

```"dataProcessor.h" ```- Likely contains data processing utilities

```"iterativeNL.h" ```- Possibly contains numerical methods

```"rootFinding.h" ```- Likely contains root-finding algorithms

```"thermoElectroChem.h"``` - Probably contains thermodynamic and electrochemical definitions

# ThermoElectroChem
_________________
# thermoElectroChem.h
## Description
This header file defines thermodynamic and electrochemical data structures and function declarations used for calculations involving metals and fluorides. It contains maps for element properties (atomic numbers, weights, oxidation states) and thermodynamic data for heat capacity calculations.


## Data Structures and Maps

```c++
const std::unordered_map<std::string, int> atomNumMap
```
* **Summary:** Maps element symbols to their atomic numbers
* **Content:** All elements from H (1) to Og (118)

```c++
const std::unordered_map<std::string, double> atomWeightMap
```
* **Summary:** Maps element symbols to their atomic weights in g/mol
* **Content:** Elements from H (1.0080) to U (238.03)
* **Note:** Comment suggests this map may not be used anywhere

```c++
const std::unordered_map<std::string, int> oxiStateMap
```
* **Summary:** Maps element symbols to their common oxidation states
* **Content:** Includes alkali metals (+1), alkaline earth metals (+2), halogens (-1), and various transition/rare earth elements

```c++
const std::unordered_map<std::string, strVect> surrogateMap
```
* **Summary:** Groups elements by surrogate categories
* **Content:** 5 surrogate groups (I, Ca, La, Pu, Th) with their respective element members
* **Purpose:** Used for substituting elements with similar chemical properties

```c++
const std::unordered_map<std::string, std::string> surrogateMapInv
```
* **Summary:** Inverse mapping of surrogateMap
* **Content:** Maps individual elements to their surrogate group representative
* **Example:** "Pr" → "La", "Zr" → "Th"

```c++
const size_t CP_TERMS = 6
using thermoArray = std::array<double, CP_TERMS+3>
const std::unordered_map<std::string, thermoArray> heatData
```
* **Summary:** Stores thermodynamic data for various compounds
* **Content:** Each thermoArray contains:
  1. ΔHf: Enthalpy of formation at 298.15 K (kJ/mol)
  2. ΔS: Entropy at 298.15 K (J/mol·K)
  3-8. Six coefficients for heat capacity calculation: Cp = a + bT + cT² + dT³ + eT⁻¹ + fT⁻²
  9. T_ref: Reference temperature (K)

```c++
const double R = 8.314
```
* **Summary:** Universal gas constant in J/mol·K

```c++
const double gamma_UF3 = 50
const double gamma_UF4 = 0.55
const double gamma_Inf_CrF2 = 0.5
const double gamma_Inf_FeF2 = 1.6
```
* **Summary:** Activity coefficients for specific compounds


## Functions Defined
```c++
gamma_Inf_NiF2(double xLiF)
```
* **Parameters:** LiF mole fraction (xLiF)
* **Returns:** Activity coefficient for NiF2
* **Summary:** Calculates activity coefficient for NiF2 based on LiF mole fraction

```c++
G_HF(const double T)
```
* **Parameters:** Temperature T in Kelvin
* **Returns:** Gibbs free energy in J/mol
* **Summary:** Calculates Gibbs free energy of HF formation

```c++
G_CrF2(const double T)
```
* **Parameters:** Temperature T in Kelvin
* **Returns:** Gibbs free energy in J/mol
* **Summary:** Calculates Gibbs free energy of CrF2 formation

```c++
G_CrF3(const double T)
```
* **Parameters:** Temperature T in Kelvin
* **Returns:** Gibbs free energy in J/mol
* **Summary:** Calculates Gibbs free energy of CrF3 formation reaction

```c++
G_FeF2(const double T)
```
* **Parameters:** Temperature T in Kelvin
* **Returns:** Gibbs free energy in J/mol
* **Summary:** Calculates Gibbs free energy of FeF2 formation

```c++
G_FeF3(const double T)
```
* **Parameters:** Temperature T in Kelvin
* **Returns:** Gibbs free energy in J/mol
* **Summary:** Calculates Gibbs free energy of FeF3 formation reaction

```c++
G_CoF2(const double T)
```
* **Parameters:** Temperature T in Kelvin
* **Returns:** Gibbs free energy in J/mol
* **Summary:** Calculates Gibbs free energy of CoF2 formation

```c++
G_CoF3(const double T)
```
* **Parameters:** Temperature T in Kelvin
* **Returns:** Gibbs free energy in J/mol
* **Summary:** Calculates Gibbs free energy of CoF3 formation reaction

```c++
G_MnF2(const double T)
```
* **Parameters:** Temperature T in Kelvin
* **Returns:** Gibbs free energy in J/mol
* **Summary:** Calculates Gibbs free energy of MnF2 (inline function)

```c++
G_NiF2(const double T)
```
* **Parameters:** Temperature T in Kelvin
* **Returns:** Gibbs free energy in J/mol
* **Summary:** Calculates Gibbs free energy of NiF2 formation

```c++
G_NbF(const double T)
```
* **Parameters:** Temperature T in Kelvin
* **Returns:** Gibbs free energy in J/mol
* **Summary:** Calculates Gibbs free energy of NbF formation

```c++
G_NbF2(const double T)
```
* **Parameters:** Temperature T in Kelvin
* **Returns:** Gibbs free energy in J/mol
* **Summary:** Calculates Gibbs free energy of NbF2 formation

```c++
G_NbF3(const double T)
```
* **Parameters:** Temperature T in Kelvin
* **Returns:** Gibbs free energy in J/mol
* **Summary:** Calculates Gibbs free energy of NbF3 formation

```c++
G_NbF4(const double T)
```
* **Parameters:** Temperature T in Kelvin
* **Returns:** Gibbs free energy in J/mol
* **Summary:** Calculates Gibbs free energy of NbF4 formation

```c++
G_NbF5(const double T)
```
* **Parameters:** Temperature T in Kelvin
* **Returns:** Gibbs free energy in J/mol
* **Summary:** Calculates Gibbs free energy of NbF5 formation

```c++
G_Nb2F10(const double T)
```
* **Parameters:** Temperature T in Kelvin
* **Returns:** Gibbs free energy in J/mol
* **Summary:** Calculates Gibbs free energy of Nb2F10 formation reaction

```c++
G_Nb3F15(const double T)
```
* **Parameters:** Temperature T in Kelvin
* **Returns:** Gibbs free energy in J/mol
* **Summary:** Calculates Gibbs free energy of Nb3F15 formation reaction

```c++
G_MoF4(const double T)
```
* **Parameters:** Temperature T in Kelvin
* **Returns:** Gibbs free energy in J/mol
* **Summary:** Calculates Gibbs free energy of MoF4 formation

```c++
G_MoF5(const double T)
```
* **Parameters:** Temperature T in Kelvin
* **Returns:** Gibbs free energy in J/mol
* **Summary:** Calculates Gibbs free energy of MoF5 formation reaction

```c++
G_MoF6(const double T)
```
* **Parameters:** Temperature T in Kelvin
* **Returns:** Gibbs free energy in J/mol
* **Summary:** Calculates Gibbs free energy of MoF6 formation reaction

```c++
G_UF4(const double T)
```
* **Parameters:** Temperature T in Kelvin
* **Returns:** Gibbs free energy in J/mol
* **Summary:** Calculates Gibbs free energy of UF4 formation reaction

```c++
G_F(const double xUF3, const double xUF4, const double T)
```
* **Parameters:** UF3 mole fraction, UF4 mole fraction, Temperature in Kelvin
* **Returns:** Fluorine potential in J/mol
* **Summary:** Calculates fluorine potential based on UF3/UF4 ratio and temperature

## External Function Calls
* Includes "miscellaneous.h" which is not provided but likely contains utility functions and the `strVect` type

# thermoElectroChem.cpp
## Description
This implementation file provides Gibbs free energy calculations for various metal fluoride compounds based on thermodynamic data. It uses enthalpy and entropy calculations to determine free energies of formation and reaction at specified temperatures.

## Functions Defined
```c++
gamma_Inf_NiF2(double xLiF)
```
* **Parameters:** LiF mole fraction (xLiF)
* **Returns:** Activity coefficient for NiF2
* **Summary:** Calculates activity coefficient for NiF2 using an exponential formula based on LiF concentration

```c++
static double calc_H(const thermoArray& data, const double T)
```
* **Parameters:** Thermodynamic data array, temperature in Kelvin
* **Returns:** Enthalpy in J/mol
* **Summary:** Calculates enthalpy from heat capacity data and reference values

```c++
static double calc_S(const thermoArray& data, const double T)
```
* **Parameters:** Thermodynamic data array, temperature in Kelvin
* **Returns:** Entropy in J/mol·K
* **Summary:** Calculates entropy from heat capacity data and reference values

```c++
static double calc_G(const thermoArray& data, const double T)
```
* **Parameters:** Thermodynamic data array, temperature in Kelvin
* **Returns:** Gibbs free energy in J/mol
* **Summary:** Calculates Gibbs free energy from enthalpy and entropy (G = H - TS)

```c++
double G_HF(const double T)
```
* **Parameters:** Temperature T in Kelvin
* **Returns:** Gibbs free energy in J/mol
* **Summary:** Calculates free energy of HF formation from H2 and F2

```c++
double G_UF4(const double T)
```
* **Parameters:** Temperature T in Kelvin
* **Returns:** Gibbs free energy in J/mol
* **Summary:** Calculates free energy of UF4 formation reaction from UF3 and F2

```c++
double G_CrF2(const double T)
```
* **Parameters:** Temperature T in Kelvin
* **Returns:** Gibbs free energy in J/mol
* **Summary:** Calculates free energy of CrF2 formation with temperature-dependent Cr properties

```c++
double G_CrF3(const double T)
```
* **Parameters:** Temperature T in Kelvin
* **Returns:** Gibbs free energy in J/mol
* **Summary:** Calculates free energy of CrF3 formation reaction from CrF2 and F2

```c++
double G_FeF2(const double T)
```
* **Parameters:** Temperature T in Kelvin
* **Returns:** Gibbs free energy in J/mol
* **Summary:** Calculates free energy of FeF2 formation with temperature-dependent Fe properties

```c++
double G_FeF3(const double T)
```
* **Parameters:** Temperature T in Kelvin
* **Returns:** Gibbs free energy in J/mol
* **Summary:** Calculates free energy of FeF3 formation reaction with multiple phase transitions

```c++
double G_CoF2(const double T)
```
* **Parameters:** Temperature T in Kelvin
* **Returns:** Gibbs free energy in J/mol
* **Summary:** Calculates free energy of CoF2 formation with temperature-dependent Co properties

```c++
double G_CoF3(const double T)
```
* **Parameters:** Temperature T in Kelvin
* **Returns:** Gibbs free energy in J/mol
* **Summary:** Calculates free energy of CoF3 formation reaction from CoF2 and F2

```c++
double G_NiF2(const double T)
```
* **Parameters:** Temperature T in Kelvin
* **Returns:** Gibbs free energy in J/mol
* **Summary:** Calculates free energy of NiF2 formation with multiple Ni phase transitions

```c++
double G_NbF5(const double T)
```
* **Parameters:** Temperature T in Kelvin
* **Returns:** Gibbs free energy in J/mol
* **Summary:** Calculates free energy of NbF5 formation from Nb and F2

```c++
double G_Nb2F10(const double T)
```
* **Parameters:** Temperature T in Kelvin
* **Returns:** Gibbs free energy in J/mol
* **Summary:** Calculates free energy of Nb2F10 formation reaction from two NbF5 molecules

```c++
double G_Nb3F15(const double T)
```
* **Parameters:** Temperature T in Kelvin
* **Returns:** Gibbs free energy in J/mol
* **Summary:** Calculates free energy of Nb3F15 formation reaction from NbF5 and Nb2F10

```c++
double G_NbF4(const double T)
```
* **Parameters:** Temperature T in Kelvin
* **Returns:** Gibbs free energy in J/mol
* **Summary:** Calculates free energy of reaction from NbF4 to NbF5 with F2

```c++
double G_NbF3(const double T)
```
* **Parameters:** Temperature T in Kelvin
* **Returns:** Gibbs free energy in J/mol
* **Summary:** Calculates free energy of reaction from NbF3 to NbF4 with F2

```c++
double G_NbF2(const double T)
```
* **Parameters:** Temperature T in Kelvin
* **Returns:** Gibbs free energy in J/mol
* **Summary:** Calculates free energy of reaction from NbF2 to NbF3 with F2

```c++
double G_NbF(const double T)
```
* **Parameters:** Temperature T in Kelvin
* **Returns:** Gibbs free energy in J/mol
* **Summary:** Calculates free energy of reaction from NbF to NbF2 with F2

```c++
double G_MoF4(const double T)
```
* **Parameters:** Temperature T in Kelvin
* **Returns:** Gibbs free energy in J/mol
* **Summary:** Calculates free energy of MoF4 formation from Mo and F2

```c++
double G_MoF5(const double T)
```
* **Parameters:** Temperature T in Kelvin
* **Returns:** Gibbs free energy in J/mol
* **Summary:** Calculates free energy of MoF5 formation reaction from MoF4 and F2

```c++
double G_MoF6(const double T)
```
* **Parameters:** Temperature T in Kelvin
* **Returns:** Gibbs free energy in J/mol
* **Summary:** Calculates free energy of MoF6 formation reaction from MoF5 and F2

```c++
double G_F(const double xUF3, const double xUF4, const double T)
```
* **Parameters:** UF3 mole fraction, UF4 mole fraction, Temperature in Kelvin
* **Returns:** Fluorine potential in J/mol
* **Summary:** Calculates fluorine potential as defined by Olander (2001) based on UF3/UF4 ratio

## External Function Calls
* Uses `pow()` and mathematical functions from the C++ standard math library
* Accesses data from `heatData` map and constants defined in the header file
* Uses the gas constant `R` from the header file
  
# Miscellaneous
_________________
# miscellaneous.h
## Description
This header file defines various utility functions and constants for mathematical and string manipulation operations. It includes mathematical constants (e, phi, pi), string processing functions, and template functions for numeric operations that appear to be used across other parts of the application.

## Functions Defined
```c++
std::string elementSymb(const std::string&)
```
* **Parameters:** A string reference
* **Returns:** A string
* **Summary:** Converts an element name to proper case (first letter uppercase, rest lowercase).

```c++
strVect strToVect(std::string&, char=' ')
```
* **Parameters:** A string reference and an optional character delimiter (defaults to space)
* **Returns:** A vector of strings (strVect)
* **Summary:** Splits a string into a vector of strings using the specified delimiter.

```c++
std::vector<double> strToVectDouble(const std::string&)
```
* **Parameters:** A constant string reference
* **Returns:** A vector of doubles
* **Summary:** Converts a string to a vector of doubles, supporting various parsing options.

```c++
bool containsNumber(const std::string&) noexcept
```
* **Parameters:** A constant string reference
* **Returns:** A boolean
* **Summary:** Checks if a string contains any numeric digits.

```c++
bool isNumeric(const std::string&) noexcept
```
* **Parameters:** A constant string reference
* **Returns:** A boolean
* **Summary:** Determines if a string can be converted to a number.

```c++
void convertibleNum(std::string&)
```
* **Parameters:** A string reference
* **Returns:** void
* **Summary:** Trims leading characters until string can be converted to a double.

```c++
template<typename T> inline constexpr int sgn(T x) noexcept
```
* **Parameters:** A value of type T
* **Returns:** An integer (-1, 0, or 1)
* **Summary:** Returns the sign of a number as an integer.

```c++
inline int factorial(int n)
```
* **Parameters:** An integer
* **Returns:** An integer
* **Summary:** Calculates the factorial of a number using recursion.

```c++
inline int ordMag(double x)
```
* **Parameters:** A double
* **Returns:** An integer
* **Summary:** Calculates the order of magnitude of a number.

## External Function Calls
The file includes the following external dependencies:
* C++ Standard Library's `<string>`: For string operations
* C++ Standard Library's `<vector>`: For vector data structure
* C++ Standard Library's `<cmath>`: For mathematical functions like exp(), sqrt(), atan(), log10(), and fabs()

# miscellaneous.cpp
## Description
This implementation file provides the code for the functions declared in miscellaneous.h. It includes string manipulation utilities for parsing and converting data between different formats, with special handling for numeric string processing.

## Functions Defined
```c++
std::string elementSymb(const std::string& str)
```
* **Parameters:** A constant string reference
* **Returns:** A string
* **Summary:** Converts first character to uppercase and rest to lowercase, mimicking Python's string casing functions.

```c++
strVect strToVect(std::string& data, char deli)
```
* **Parameters:** A string reference and a character delimiter
* **Returns:** A vector of strings (strVect)
* **Summary:** Splits a string into a vector of strings by a delimiter, handling both spaces and tabs.

```c++
std::vector<double> strToVectDouble(const std::string& str)
```
* **Parameters:** A constant string reference
* **Returns:** A vector of doubles
* **Summary:** Parses a string into a vector of doubles, supporting space, colon, and comma-separated formats.

```c++
bool isNumeric(const std::string& str) noexcept
```
* **Parameters:** A constant string reference
* **Returns:** A boolean
* **Summary:** Checks if a string represents a valid number, handling decimals, scientific notation, and signs.

```c++
bool containsNumber(const std::string& str) noexcept
```
* **Parameters:** A constant string reference
* **Returns:** A boolean
* **Summary:** Determines if a string contains any numeric digit.

```c++
void convertibleNum(std::string& str)
```
* **Parameters:** A string reference
* **Returns:** void
* **Summary:** Removes leading characters until the string can be converted to a double.

## External Function Calls
The file includes the following external dependencies:
* `miscellaneous.h`: For declarations of functions implemented in this file
* C++ Standard Library's `<stdexcept>`: For exception handling classes like std::invalid_argument and std::domain_error
* C++ Standard Library character functions: toupper(), tolower(), isdigit()
* C++ Standard Library conversion functions: std::stod() for string to double conversion
* C++ Standard Library math functions: ceil()

# GaussElim
_________________
# gaussElim.h

## Description
This header file defines functions for performing Gaussian elimination and related matrix operations. It declares functions for performing LU factorization, solving linear equations, calculating determinants, and finding inverse matrices.

## Functions Defined
```cpp
Matrix luFactorize(Matrix& a, Permutation& p, bool inplace=false)
```
* **Parameters:** A matrix `a`, a permutation array `p`, and a boolean `inplace` (default: false)
* **Returns:** The factorized matrix
* **Summary:** Performs LU factorization on a square matrix, optionally modifying the original matrix

```cpp
Vector luSolve(const Matrix& a, const Permutation& p, Vector& x, bool inplace=false)
```
* **Parameters:** A factorized matrix `a`, permutation array `p`, vector `x`, and boolean `inplace` (default: false)
* **Returns:** Solution vector
* **Summary:** Solves a linear system using a pre-factorized matrix

```cpp
Vector solve(Matrix& a, Permutation& p, Vector& x, bool inplace=false)
```
* **Parameters:** A matrix `a`, permutation array `p`, vector `x`, and boolean `inplace` (default: false)
* **Returns:** Solution vector
* **Summary:** Factorizes the matrix and solves the linear system

```cpp
double detFactoredMatrix(const Matrix& a, const Permutation& p) noexcept
```
* **Parameters:** A factorized matrix `a` and permutation array `p`
* **Returns:** Determinant of the matrix
* **Summary:** Calculates determinant from a pre-factorized matrix

```cpp
double det(Matrix& a) noexcept
```
* **Parameters:** A matrix `a`
* **Returns:** Determinant of the matrix
* **Summary:** Factorizes the matrix and calculates its determinant

```cpp
Matrix invFactoredMatrix(const Matrix& a, const Permutation& p)
```
* **Parameters:** A factorized matrix `a` and permutation array `p`
* **Returns:** Inverse matrix
* **Summary:** Calculates inverse of a pre-factorized matrix

```cpp
Matrix inv(Matrix& a)
```
* **Parameters:** A matrix `a`
* **Returns:** Inverse matrix
* **Summary:** Factorizes the matrix and calculates its inverse

## External Dependencies
* `matrix.h` - Provides Matrix and Permutation classes

# gaussElim.cpp

## Description
This file implements the functions declared in gaussElim.h, providing algorithms for Gaussian elimination with scaled partial pivoting to perform LU factorization, solve linear systems, calculate determinants, and find inverse matrices.

## Functions Defined

```cpp
static void swapRows(Matrix& a, int i, int j)
```
* **Parameters:** A matrix `a`, row indices `i` and `j` 
* **Returns:** void
* **Summary:** Swaps two rows in a matrix

```cpp
static void swap(Vector& v, int i, int j)
```
* **Parameters:** A vector `v`, indices `i` and `j`
* **Returns:** void
* **Summary:** Swaps two elements in a vector

```cpp
static void rowReplacement(Matrix& a, int i, int j)
```
* **Parameters:** A matrix `a`, row indices `i` and `j`
* **Returns:** void
* **Summary:** Performs row replacement operation in Gaussian elimination

## External Dependencies
* `matrix.h` - Provides Matrix, Vector, and Permutation classes
* `cmath` - For mathematical functions like `abs` and `fabs`

# IterativeNL
_________________
# iterativeNL.h

## Description
This header file declares functions for solving nonlinear equations using various iterative methods including fixed point iteration, Newton's method, and Broyden's methods.

## Functions Defined

```cpp
Vector fixedpt(std::function<Vector(const Vector&)> g, Vector& x, int maxIter, double tol)
```
* **Parameters:** Function `g`, initial guess vector `x`, maximum iterations `maxIter`, tolerance `tol`
* **Returns:** Solution vector
* **Summary:** Implements fixed point iteration method for vector functions

```cpp
Vector newton(std::function<Vector(const Vector&)> f, std::function<Matrix(const Vector&)> df, Vector& x, int maxIter, double tol)
```
* **Parameters:** Function `f`, Jacobian `df`, initial guess `x`, maximum iterations `maxIter`, tolerance `tol`
* **Returns:** Solution vector
* **Summary:** Implements Newton's method with analytical Jacobian for vector functions

```cpp
Vector newton(std::function<Vector(const Vector&)> f, Vector& x, int maxIter, double tol)
```
* **Parameters:** Function `f`, initial guess `x`, maximum iterations `maxIter`, tolerance `tol`
* **Returns:** Solution vector
* **Summary:** Implements Newton's method with numerical Jacobian for vector functions

```cpp
Vector broyden1(std::function<Vector(const Vector&)> f, Matrix& A, Vector& x, int maxIter, double tol)
```
* **Parameters:** Function `f`, initial Jacobian approximation `A`, initial guess `x`, maximum iterations `maxIter`, tolerance `tol`
* **Returns:** Solution vector
* **Summary:** Implements Broyden's first method for vector functions

```cpp
Vector broyden2(std::function<Vector(const Vector&)> f, Matrix& B, Vector& x, int maxIter, double tol)
```
* **Parameters:** Function `f`, initial inverse Jacobian approximation `B`, initial guess `x`, maximum iterations `maxIter`, tolerance `tol`
* **Returns:** Solution vector
* **Summary:** Implements Broyden's second method for vector functions

## External Dependencies
* `functional` - For std::function
* `matrix.h` - For Matrix and Vector classes

# iterativeNL.cpp

## Description
This file implements the functions declared in iterativeNL.h, providing algorithms for various iterative methods to solve systems of nonlinear equations, including fixed point iteration, Newton's method, and Broyden's methods.

## Functions Defined

```cpp
static bool checkTol(const Vector& x, const Vector& dx, const double tol)
```
* **Parameters:** Current solution `x`, step vector `dx`, tolerance `tol`
* **Returns:** Boolean indicating if convergence is achieved
* **Summary:** Checks if relative change in solution is within tolerance

## External Dependencies
* `calculus.h` - Provides numerical differentiation functions for the Jacobian
* `gaussElim.h` - For solving linear systems in each iteration
* `matrix.h` - For Matrix and Vector classes
* `cmath` - For mathematical functions

# RootFinding
_________________
# rootFinding.h

## Description
This header file declares functions for finding roots of scalar nonlinear equations using various methods including bisection, fixed point iteration, Newton's method, and the secant method.

## Functions Defined

```cpp
double bisection(std::function<double(double)> f, double x0, double x1, int maxIter=20, double tol=1e-6)
```
* **Parameters:** Function `f`, interval endpoints `x0` and `x1`, maximum iterations `maxIter` (default: 20), tolerance `tol` (default: 1e-6)
* **Returns:** Root approximation
* **Summary:** Implements the bisection method for finding roots

```cpp
double fixedpt(std::function<double(double)> f, double x0, int maxIter=20, double tol=1e-6)
```
* **Parameters:** Function `f`, initial guess `x0`, maximum iterations `maxIter` (default: 20), tolerance `tol` (default: 1e-6)
* **Returns:** Fixed point
* **Summary:** Implements fixed point iteration for scalar functions

```cpp
double newton(std::function<double(double)> f, std::function<double(double)> df, double x0, int maxIter=20, double tol=1e-6)
```
* **Parameters:** Function `f`, derivative `df`, initial guess `x0`, maximum iterations `maxIter` (default: 20), tolerance `tol` (default: 1e-6)
* **Returns:** Root approximation
* **Summary:** Implements Newton's method with analytical derivative

```cpp
double newton(std::function<double(double)> f, double x0, int maxIter=20, double tol=1e-6)
```
* **Parameters:** Function `f`, initial guess `x0`, maximum iterations `maxIter` (default: 20), tolerance `tol` (default: 1e-6)
* **Returns:** Root approximation
* **Summary:** Implements Newton's method with numerical derivative

```cpp
double secant(std::function<double(double)> f, double x0, double x1, int maxIter=20, double tol=1e-6)
```
* **Parameters:** Function `f`, initial guesses `x0` and `x1`, maximum iterations `maxIter` (default: 20), tolerance `tol` (default: 1e-6)
* **Returns:** Root approximation
* **Summary:** Implements the secant method for finding roots

## External Dependencies
* `functional` - For std::function

# rootFinding.cpp

## Description
This file implements the functions declared in rootFinding.h, providing algorithms for root-finding methods including bisection, fixed point iteration, Newton's method, and the secant method for scalar equations.

## External Dependencies
* `calculus.h` - For numerical differentiation (deriv function)
* `miscellaneous.h` - For utility functions (sgn function)
* `math.h` - For mathematical functions

# Main
_________________
# main.cpp

## Description
This file contains the main program with a menu-driven interface for data processing operations related to the SCALE nuclear code and Thermochimica thermochemical calculation software. It allows users to convert data between formats, extract information, combine files, and decouple surrogate elements.

## Functions Defined
The file primarily contains a `main()` function with no other function definitions.

## External Dependencies
* `dataProcessor.h` - For all data processing functions used in the menu options:
  - `scaleToVector` - Converts SCALE output to vector format
  - `vectToTherm` - Converts vector data to Thermochimica input
  - `textToExcel` - Extracts data from Thermochimica output
  - `mergeTherm` - Combines multiple Thermochimica outputs
  - `decoupleSurr` - Decouples surrogate elements
