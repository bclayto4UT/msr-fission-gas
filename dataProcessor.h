#ifndef DATAPROCESSOR_H_INCLUDED
#define DATAPROCESSOR_H_INCLUDED

#include "miscellaneous.h"

std::vector<strVect> scaleToVector(const std::string&, bool=true);
void vectToThermI(const std::vector<strVect>&, const std::string&);

void textToExcel(const std::string&, const std::string&, std::string&);

void mergeTherm(const strVect&, const std::string&);

void coupleSurr(const std::vector<strVect>&, const std::string&);
void decoupleSurr(const std::vector<strVect>&,
                  const std::string&, const std::string&,
                  const bool=false, const bool=false);

/*void massToMole(const std::string&, const std::map<std::string, double>&,
                double mass=1.0);

void vectToThermO(const std::vector<strVect>&, const std::string&,
                  const strVect*);

void rmvSpecies(const std::vector<strVect>&, const std::string&,
                const std::map<int, double>&, bool);
void rmvSpecies(const std::vector<strVect>&, const std::string&,
                const double, bool);*/


#endif // DATAPROCESSOR_H_INCLUDED
