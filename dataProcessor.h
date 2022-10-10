#ifndef DATAPROCESSOR_H_INCLUDED
#define DATAPROCESSOR_H_INCLUDED

#include "miscellaneous.h"

std::vector<strVect> scaleToVector(const std::string&);
void vectToThermI(const std::vector<strVect>&, const std::string&);
void vectToThermO(const std::vector<strVect>&, const std::string&,
                  const strVect*);

void rmvSpecies(const std::vector<strVect>&, const std::string&,
                const std::map<int, double>&, bool);
void rmvSpecies(const std::vector<strVect>&, const std::string&,
                const double, bool);

void textToExcel(const std::string&, const std::string&, std::string&);
void mergeTherm(const strVect&, const std::string&);

void massToMole(const std::string&, const std::map<std::string, double>&,
                double mass=1.0);

#endif // DATAPROCESSOR_H_INCLUDED
