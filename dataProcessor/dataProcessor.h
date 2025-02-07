#ifndef DATAPROCESSOR_H_INCLUDED
#define DATAPROCESSOR_H_INCLUDED

#include "miscellaneous.h"
#include <unordered_map>

/********************************************************************
 * These functions are directly called in main. The main function here is decoupleSurr(),
 * which updates a Thermochimica result file to include more elements and compounds.
********************************************************************/

using concMap = std::unordered_map<std::string, std::vector<double>>;

concMap scaleToVector(const std::string&);
void vectToTherm(const concMap&, const std::string&,
                 std::string&, std::string&, bool=true);

void textToExcel(const std::string&, const std::string&, std::string&);

void mergeTherm(const strVect&, const std::string&);

void decoupleSurr(const concMap&,
                  const std::string&, const std::string&,
                  const bool=false, const bool=false);

#endif // DATAPROCESSOR_H_INCLUDED
