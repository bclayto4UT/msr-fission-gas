#ifndef DATAPROCESSOR_H_INCLUDED
#define DATAPROCESSOR_H_INCLUDED

#include "miscellaneous.h"

std::vector<strVect> scaleToVector(const std::string&);
void vectToTherm(const std::vector<strVect>&, const std::string&,
                 std::string&, std::string&, bool=true);

void textToExcel(const std::string&, const std::string&, std::string&);

void mergeTherm(const strVect&, const std::string&);

void decoupleSurr(const std::vector<strVect>&,
                  const std::string&, const std::string&,
                  const bool=false, const bool=false);

#endif // DATAPROCESSOR_H_INCLUDED
