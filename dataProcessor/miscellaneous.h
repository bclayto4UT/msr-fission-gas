#ifndef MISCELLANEOUS_H_INCLUDED
#define MISCELLANEOUS_H_INCLUDED

#include <string>
#include <vector>
#include <cmath>

/* Miscellaneous functions, some may even have a Python equivalent that does the same thing. */

constexpr double e = exp(1);
constexpr double phi = (1+sqrt(5))/2;
constexpr double pi = atan(1)*4;

std::string elementSymb (const std::string&);

using strVect = std::vector<std::string>;
strVect strToVect (std::string&, char=' ');
std::vector<double> strToVectDouble (const std::string&);

bool containsNumber (const std::string&) noexcept;
bool isNumeric (const std::string&) noexcept;
void convertibleNum (std::string&);

template<typename T>
inline constexpr int sgn(T x) noexcept{ return (x>T(0))-(x<T(0));}

inline int factorial(int n){ return (n == 0 || n == 1) ? 1 : n*factorial(n-1);}
inline int ordMag(double x){ return log10(fabs(x));}

#endif // MISCELLANEOUS_H_INCLUDED
