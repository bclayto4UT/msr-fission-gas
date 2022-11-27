#ifndef MISCELLANEOUS_H_INCLUDED
#define MISCELLANEOUS_H_INCLUDED

#include <array>
#include <map>
#include <string>
#include <vector>
#include <cmath>

constexpr double e = exp(1);
constexpr double phi = (1+sqrt(5))/2;
constexpr double pi = atan(1)*4;

std::string findElement (int);
std::string elementSymb (const std::string&);

using strVect = std::vector<std::string>;
strVect strToVect (std::string&);

bool containsNumber (const std::string&) noexcept;
void convertibleNum (std::string&);

template<typename T>
inline constexpr int sgn(T x) noexcept{ return (x>T(0))-(x<T(0));}

inline int factorial(int n){ return (n == 0 || n == 1) ? 1 : n*factorial(n-1);}
inline int ordMag(double x){ return log10(fabs(x));}

//inline bool isSorted(const Vector& v, bool strict=false) noexcept{
//    if (v.n() < 3) return true;
//    int dir = 0;
//    int start = 0;
//    for (int i = 0; i < v.n()-1; i++){
//        dir = sgn(v(i+1) - v(i));
//        if (dir == 0 && strict) return false;
//        if (dir != 0) break;
//        start++;
//    }
//    for (int i = start; i < v.n()-1; i++){
//        if (sgn(v(i+1) - v(i)) == -dir) return false;
//        if (strict)
//            if (sgn(v(i+1) - v(i)) == 0) return false;
//    }
//    return true;
//}

#endif // MISCELLANEOUS_H_INCLUDED
