# Code Analysis

## dataProcessor.h
### Description
This header file defines functions for processing data related to Thermochimica, which appears to be a thermochemical equilibrium calculation software. The main functionality focuses on decoupling surrogate elements and updating Thermochimica result files to include additional elements and compounds.

 ### Functions Defined
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

### External Function Calls
The file includes the following external dependencies:

```miscellaneous.h```: Likely provides the strVect type used in the code
```Standard C++ library's <unordered_map>```: Provides the unordered_map data structure used to define the concMap type

## dataProcessor.cpp
### Description
This C++ code appears to be part of a thermochemical data processing tool. It handles concentration data from simulation outputs, converts them to Thermochimica input format, and extracts specific data from Thermochimica output for analysis and visualization.

 ### Functions Defined
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

#### Lambda Functions Defined
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

### External Functions Called

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


### Header Files Included

#### Standard libraries:

```<algorithm>``` - For standard algorithms

```<iostream> ```- For input/output operations

```<fstream>``` - For file operations


#### Custom headers:

```"dataProcessor.h" ```- Likely contains data processing utilities

```"iterativeNL.h" ```- Possibly contains numerical methods

```"rootFinding.h" ```- Likely contains root-finding algorithms

```"thermoElectroChem.h"``` - Probably contains thermodynamic and electrochemical definitions
