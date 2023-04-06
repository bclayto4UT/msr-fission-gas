#include <algorithm>
#include <iostream>
#include <fstream>

#include "dataProcessor.h"
#include "iterativeNL.h"
#include "rootFinding.h"
#include "thermoElectroChem.h"

static const double FRACTION_CUTOFF = 1.0E-99;

concMap scaleToVector(const std::string& inFile)
/****************************************************************
Input(s):
    inFile:     the name of the file with .out file to be read
Output:
    a map that arranges the concentration values under the element
    names
*****************************************************************/
{
    // Opens input file
    std::ifstream input(inFile);
    std::string line;

    if (input.is_open()) getline(input,line);
    else throw std::invalid_argument("Cannot open input file.");

    std::string beginStr = "relative cutoff;";
    std::string endStr = "------";
    bool beginStrIsFound = false;

    while (line.find(beginStr) == std::string::npos) getline(input,line);
    beginStrIsFound = line.find(beginStr) >= 0;
    if (!beginStrIsFound) throw std::invalid_argument("Input file does not contain an element table.");

    getline(input,line);
    getline(input,line); // The line with the the intervals is the second one after beginStr
    concMap concen;

    // Inputs data into vector
    while(getline(input,line) && line.find(endStr) == std::string::npos){
        // Converts element symbol to atomic number and stores it
        auto pos = std::min(line.find("\t"), line.find(" "));
        while (pos == 0){
            line.erase(0, 1);
            pos = std::min(line.find("\t"), line.find(" "));
        }
        std::string ele = line.substr(0, pos);
        ele[0] = toupper(ele[0]);
        line.erase(0, pos);

        std::vector<std::string> vs = strToVect(line);
        std::vector<double> vd(vs.size());
        for (auto it = vs.cbegin(); it != vs.cend(); ++it){
            auto it2 = vd.begin() + (it-vs.cbegin());
            try{
                *it2 = stod(*it);
            } catch (std::invalid_argument& ex){
                throw std::invalid_argument("Non-numerical value encountered in concentration table.");
            }
        }
        concen[ele] = vd;
    }

    input.close();
    return concen;
}

void vectToTherm(const concMap& dataMap,
                 const std::string& outFile,
                 std::string& strT,
                 std::string& strP,
                 bool includesSurr)
/****************************************************************
Input(s):
    dataVect:     the concentration map given in scaleToVector
    outFile:      the output file name
    strT:         system temperature (in K)
    strP:         system pressure (in atm)
    includesSurr: whether surrogate elements are to be included
                  (DEFAULT: true)
Output: none
    a text file formatted as a Thermochimica input is created
    if there are multiple rows, multiple files are created
*****************************************************************/
{
    std::ofstream output;

    size_t pos = outFile.find("."); // Finds the file extension
    std::string fileStem;
    if (pos < outFile.size()) fileStem = outFile.substr(0, pos);
    else fileStem = outFile;

    pos = fileStem.rfind("/");
    std::string relFileStem;
    if (pos < fileStem.size()) relFileStem = fileStem.substr(pos+1, fileStem.size()-pos-1);
    else relFileStem = fileStem;

    std::size_t m = dataMap.cbegin()->second.size(); // m = number of time intervals

    std::vector<double> T, P;
    try{
        T = strToVectDouble(strT);
    } catch(const std::domain_error&){
        auto vs = strToVect(strT,':');
        double start = std::stod(vs[0]);
        double stop = std::stod(vs[1]);

        if (start == stop) T.push_back(start);
        else{
            double step = (stop-start)/m;
            T.reserve(m);
            for (size_t i = 0; i < m; i++) T.push_back(start+step*i);
        }

    } catch(const std::invalid_argument& ex){
        throw ex;
    } catch(const std::bad_alloc& ex){
        throw ex;
    }

    if (T.size() != 1 && T.size() < m)
        throw std::out_of_range("Temperature array out of range.");
    for (auto d: T){
        if (d <= 0) throw std::invalid_argument("Temperature must be positive.");
    }

    try{
        P = strToVectDouble(strP);
    } catch(const std::domain_error&){
        auto vs = strToVect(strP,':');
        double start = std::stod(vs[0]);
        double stop = std::stod(vs[1]);

        if (start == stop) P.push_back(start);
        else{
            double step = (stop-start)/m;
            P.reserve(m);
            for (size_t i = 0; i < m; i++) P.push_back(start+step*i);
        }
    } catch(const std::invalid_argument& ex){
        throw ex;
    } catch(const std::bad_alloc& ex){
        throw ex;
    }
    if (P.size() != 1 && P.size() < m)
    throw std::out_of_range("Pressure array out of range.");
    for (auto d: P){
        if (d <= 0) throw std::invalid_argument("Pressure must be positive.");
    }

//     Recalculates species amounts for surrogated elements
    const concMap* actualMap;
    concMap dataMap2;

    if (includesSurr){
        dataMap2 = dataMap;

        for (auto it = dataMap2.begin(); it != dataMap2.end(); ++it){
            std::string ele = (*it).first;
            std::string surr;

            // Checks to see if the element needs a surrogate
            try{ surr = surrogateMapInv.at(ele);}
            catch(const std::out_of_range& ex){ continue;}

            // Checks to see if the surrogate exists within dataMap
            try{ dataMap2.at(surr);}
            catch(const std::out_of_range& ex){ dataMap2[surr] = std::vector<double>(m);}
            for (size_t i = 0; i < m; i++) dataMap2[surr][i] += (*it).second[i];

            // Erases the element now that it is taken into account by the surrogate
            dataMap2.erase(it);
        }
        actualMap = &dataMap2;
    }
    else actualMap = &dataMap;

    for (size_t i = 0; i < m; i++){
        std::string fileName = fileStem + std::to_string(i) + ".F90";
        output.open(fileName);

        output << "program " + relFileStem + std::to_string(i) + "\n";
        output << "USE ModuleThermoIO\n";
        output << "implicit none\n";
        output << "cInputUnitTemperature = 'K'\n";
        output << "cInputUnitPressure = 'atm'\n";
        output << "cInputUnitMass = 'moles'\n";
        output << "cThermoFileName = DATA_DIRECTORY // 'MSTDB-TC_V2.0_Fluorides_8-0.dat'\n";

        if (P.size() == 1) output << "dPressure = " << P[0] << "\n";
        else output << "dPressure = " << P[i] << "\n";
        if (T.size() == 1) output << "dTemperature = " << T[0] << "\n";
        else output << "dTemperature = " << T[i] << "\n";


        for (auto it: *actualMap){
            if (it.second[i] > 0.0){
                output << "dElementMass(";
                output << atomNumMap.at(it.first);
                output << ") = ";
                output << it.second[i];
                output << "\n";
            }
        }

        output << "iPrintResultsMode = 1\n";
        output << "call ParseCSDataFile(cThermoFileName)\n";
        output << "if (INFOThermo == 0) call Thermochimica\n";
        output << "if (iPrintResultsMode > 0) call PrintResults\n";
        output << "if (INFOThermo == 0) call ResetThermoAll\n";
        output << "call ThermoDebug\n";
        output << "end program " + relFileStem + std::to_string(i) + "\n";

        output.close();
    }
}

void textToExcel(const std::string& inFile,
                 const std::string& outFile,
                 std::string& dataType)
/****************************************************************
Input(s):
    inFile:   the name of the input file that contains a
              Thermochimica output
    outFile:  the output file name
    dataType: a string that contains the types of data that the
              user wishes to extract. It can include as many
              types as necessary. The accepted types include:
     * ni:    the total amount (moles) of ions
     * nx:    the total amount (moles) of salt
     * ny:    the total amount (moles) of gas
     * x_ABC: the mole fraction of species ABC in the salt phase
     * y_ABC: the mole fraction of species ABC in the gas phase
              "ABC" can be "all", which will output all salts
              (not ions) and all gases
     * T:     the system temperature
       For example, dataType = "x_UF3 x_UF4 y ny_UF5 y_UF6"
Output(s): none
    a textfile in tabulated form, with each data type on the
    columns, is created
*****************************************************************/
{
    auto v = strToVect(dataType);
    const bool has_x_all = std::find(v.cbegin(), v.cend(), "x_all") != v.cend();
    const bool has_y_all = std::find(v.cbegin(), v.cend(), "y_all") != v.cend();

    // Initializes two std::map's that map a species to its amount
    std::unordered_map<std::string, double> vx;
    std::unordered_map<std::string, double> vy;
    double T = -1;
    for (auto s: v){
        if (s == "ni" || s == "nx") vx[s] = 0;
        else if (s == "ny") vy[s] = 0;
        else if (s.substr(0,1) == "x"){
            if (s == "x_U4+"){
                vx["U2"] = 0;
                vx["U[VI]"] = 0;
                vx["U[VII]"] = 0;
            } if (s == "x_UF4" || has_x_all){
                vx["U2F8"] = 0;
                vx["U[VI]-F4"] = 0;
                vx["U[VII]-F4"] = 0;
            } if (s == "x_UI4" || has_x_all){
                vx["U2I8"] = 0;
                vx["U[VI]-I4"] = 0;
                vx["U[VII]-I4"] = 0;
            } if (s == "x_Be"){
                vx["Be[IV]"] = 0;
                vx["Be2"] = 0;
            } if (s == "x_BeF2" || has_x_all){
                vx["BeF2"] = 0;
                vx["Be2F4"] = 0;
            } if (s == "x_BeI2" || has_x_all){
                vx["BeI2"] = 0;
                vx["Be2I4"] = 0;
            } if (!has_x_all) vx[s.erase(0,2)] = 0;
        }
        else if (s.substr(0,1) == "y" && !has_y_all) vy[s.erase(0,2)] = 0;
        else if (s == "T") T = 0;
    }

    // Checks if a key in the map is one of the data to be extracted
    // This is because more than one entry was created for U4+ and Be2+ species
    auto isLookedFor = [&](const std::string& str) -> bool
        {
            // List of species to be excluded, unless explicitly looked for
            // THIS NEEDS FIXING - Not excluding the species as intended.
            std::vector<std::string> configs{"x_U2F8", "x_U[VI]-F4", "x_U[VII]-F4",
                                             "x_U2I8", "x_U[VI]-I4", "x_U[VII]-I4",
                                             "x_Be2F4", "x_Be2I4"};

            auto it1 = std::find(v.cbegin(), v.cend(), str);
            auto it2 = std::find(configs.cbegin(), configs.cend(), str);
            if (it2 != configs.cend() && it1 == v.cend()) return true;
            return has_x_all || it1 != v.cend();
        };

    std::ofstream output(outFile);
    if (output.is_open()){
        if (!has_x_all){
            for (auto const& it: vx){
                std::string label = "";
                if (it.first != "nx" && it.first != "ni") label = "x_";
                label += it.first;
                if (isLookedFor(label)) output << label << " ";
            }
        }
        if (!has_y_all){
            for (auto const& it: vy){
                std::string prefix = "";
                if (it.first != "ny") prefix = "y_";
                output << prefix << it.first << " ";
            }
        }
        if (T >= 0) output << "T" << " ";
        output << std::endl;
    } else throw std::invalid_argument("Cannot open output file.");

    std::ifstream input(inFile);
    std::string line;
    auto* phaseRef = &vx;
    bool hasAll = false;
    bool vxIsFilled = false;
    bool vyIsFilled = false;

    if (input.is_open()){
        bool firstRowChanged = false;
        while (getline(input, line)){
            if (containsNumber(line)){
                if (!vxIsFilled && line.find("mol MSFL") < line.size()){
                    phaseRef = &vx;
                    if (vx.count("ni")){
                        line.erase(line.find("mol MSFL"), line.size());
                        convertibleNum(line);
                        vx["ni"] = std::stod(line);
                    }
                }

                else if (!vxIsFilled && line.find("Moles of pairs") < line.size()){
                    phaseRef = &vx;
                    hasAll = has_x_all;
                    if (vx.count("nx")){
                        line.erase(line.find("Moles of pairs"), line.size());
                        convertibleNum(line);
                        vx["nx"] = std::stod(line);
                    }
                }

                else if (!vyIsFilled && line.find("mol gas") < line.size()){
                    phaseRef = &vy;
                    hasAll = has_y_all;

                    if (vy.count("ny")){
                        line.erase(line.find("mol gas"), line.size());
                        convertibleNum(line);
                        vy["ny"] = std::stod(line);
                    }
                }

                else if (T >= 0 && line.find("Temperature") < line.size()){
                    line.erase(line.size()-3, line.size());
                    convertibleNum(line);
                    T = std::stod(line);
                }

                else{
                    auto vect = strToVect(line);
                    if (vect.back() == "}") vect.pop_back();
                    if (hasAll){
                    /* If the species is not current in phaseRef, it means it
                        is a new species to be added, and the label row has to
                        be re-printed. */
                        if (!phaseRef->count(vect.back())) firstRowChanged = true;
                        (*phaseRef)[vect.back()] = std::stod(vect[vect.size()-2]);
                    }

                    else{
                        for (auto const& it: *phaseRef){
                            if (vect.back() == it.first && it.second == 0){
                                (*phaseRef)[it.first] = std::stod(vect[vect.size()-2]);
                                break;
                            }
                        }
                    }
                }

                auto isFilled = [](const std::unordered_map<std::string, double>& m) -> bool
                {
                    if (m.empty()) return false;
                    for (auto const& it: m) if (it.second == 0) return false;
                    return true;
                };
                vxIsFilled = isFilled(vx);
                vyIsFilled = isFilled(vy);

            }

            else if (line.find("Quadruplet") < line.size()
                     || line.find("System properties") < line.size()){
                hasAll = false;
                // So that amount of quadruplet fractions are not extracted.
            }

            else if (line.find("DEBUG") < line.size()){
                if (firstRowChanged){
                    for (auto const& it: vx){
                        std::string label = "";
                        if (it.first != "nx" && it.first != "ni") label = "x_";
                        label += it.first;
                        if (isLookedFor(label)) output << label << " ";
                    }
                    for (auto const& it: vy){
                        std::string prefix = "";
                        if (it.first != "ny") prefix = "y_";
                        output << prefix << it.first << " ";
                    }
                    if (T >= 0) output << "T" << " ";
                    output << std::endl;
                    firstRowChanged = false;
                }

                for (auto const& it: vx){
                    std::string label = "";
                    if (it.first != "nx" && it.first != "ni") label = "x_";
                    label += it.first;

                    if (isLookedFor(label)){
                        double sum;
                        if (it.first == "U4+" && it.second == 0){
                            sum = vx["U2"]*2 + vx["U[VI]"] + vx["U[VII]"];
                        } else if (it.first == "UF4" && it.second == 0){
                            sum = vx["U2F8"]*2 + vx["U[VI]-F4"] + vx["U[VII]-F4"];
                        } else if (it.first == "UI4" && it.second == 0){
                            sum = vx["U2I8"]*2 + vx["U[VI]-I4"] + vx["U[VII]-I4"];
                        } else if (it.first == "Be"){
                            sum = vx["Be[IV]"] + vx["Be2"]*2;
                        } else if (it.first == "BeF2"){
                            sum = it.second + vx["Be2F4"]*2;
                        } else if (it.first == "BeI2"){
                            sum = it.second + vx["Be2I4"]*2;
                        } else sum = it.second;

                        if (it.second <= 1) output << sum << " ";
                        /* it.second > 1 usually means it is actually 1e-100 or less,
                        and the small exponential tends not to be left out. */
                        else output << 0.00 << " ";
                    }
                }

                for (auto const&it: vx) vx[it.first] = 0;

                for (auto const& it: vy){
                    if (it.second <= 1) output << it.second << " ";
                    else output << 0.00 << " ";
                    vy[it.first] = 0;
                }

                if (T >= 0) output << T << " ";
                output << std::endl;
                vxIsFilled = false;
                vyIsFilled = false;
            }
        }
    } else throw std::invalid_argument("Cannot open input file.");



    input.close();
    output.close();
}

void mergeTherm(const strVect& inFiles, const std::string& outFile)
/****************************************************************
Input(s):
    inFiles: the name(s) of the input file(s) that contains a
             Thermochimica output, in a vector
    outFile: the output file name
Output(s): none
    a textfile in the format of a Thermochimica output (only with
    moles of pairs and mol gas_ideal are supported) that sums up
    all species in each phase
Warning:
    this function is the most accurate if there is no phase other
    than a salt and an ideal gas phase, and the temperature and
    pressure across all files are uniform
*****************************************************************/
{
    double nx, ny, n;
    std::unordered_map<std::string, double> mx, my;
    std::unordered_map<std::string, double>* mapRef = nullptr;

    std::ifstream input;
    std::string line;

    for (auto& s: inFiles){
        input.open(s);

        if (input.is_open()){
            while (getline(input, line)){
                if (line.find("==") < line.size()){
                    mapRef = nullptr;

                } else if (containsNumber(line)){
                    auto v = strToVect(line);
                    if (v.back() == "pairs"){
                        // At the salt phrase
                        if (v.front() == "+") n = std::stod(v[1]);
                        else n = std::stod(v.front());
                        nx += n;
                        mapRef = &mx;

                    } else if (v.back() == "gas_ideal"){
                        // At the gas phrase
                        if (v.front() == "+") n = std::stod(v[1]);
                        else n = std::stod(v.front());
                        ny += n;
                        mapRef = &my;

                    } else if (mapRef){
                        // Now reads the components
                        std::string key = v.back() == "}" ? v[v.size()-2] : v.back();
                        bool check = v.front() == "+" || v.front() == "{";
                        double val = check ? std::stod(v[1]) : std::stod(v.front());
                        (*mapRef)[key] += n*val;
                        if (v.back() == "}") mapRef = nullptr;
                    }
                }
            }
            input.close();

        } else continue; // Should an exception be thrown?

    }

    // Compute new mole fractions and output results
    std::ofstream output(outFile);
    if (nx > 0){
        output << "\t" << nx << " Moles of pairs\n\t {\n ";
        for (auto const& it: mx){
            if (it.second/nx <= 1.0){
                output << "\t\t" << it.second/nx << "\t" << it.first;
                output << std::endl;
            }
        }
        output << "\t }" << std::endl;
    }
    if (ny > 0){
        output << "\t";
        if (nx > 0) output << "+ ";
        output << ny << " mol gas_ideal\n\t {\n ";
        for (auto const& it: my){
            if (it.second/ny <= 1.0){
                output << "\t\t" << it.second/ny << "\t" << it.first;
                output << std::endl;
            }
        }
        output << "\t }" << std::endl;
    }
}

void decoupleSurr(const concMap& scaleData,
                  const std::string& thermoRes,
                  const std::string& thermoOut,
                  const bool includesSS,
                  const bool includesFP)
/****************************************************************
Input(s):
    scaleData:    the original concentration map obtained using the
                  scaleToVector function)
    thermoRes:    the name of the file containing Thermochimica re-
                  results that have been ran on data with surrogates
    thermoOut:    the output file name
    includesSS:   whether fluorides of metals in the reactor struc-
                  ture are considered (DEFAULT: false).
    includesFP:   whether fluorides of fission products (H, Nb) are
                  considered (DEFAULT: false).
Output: none
    New Thermochimica results are written in the output file,
    in which each salt containing a surrogate metal or nonmetal
    has been broken further down to its original counterpart(s),
    with the ratio of the corresponding elements preserved.
Warning:
    only a prediction, works best on hypo-stoichiometric systems.
*****************************************************************/
{
    std::ifstream input(thermoRes);
    std::string line;
    if (input.is_open()) getline(input,line);
    else throw std::invalid_argument("Cannot open input file.");
    std::ofstream output(thermoOut);

    std::size_t m = scaleData.cbegin()->second.size(); // number of time intervals
    if (m == 0) return;

    const std::array<std::string, 6> gases{"H", "He", "Ne", "Ar", "Kr", "Xe"};
    const std::array<std::string, 6> metals{"Cr", "Fe", "Co", "Ni", "Nb", "Mo"};

    // lambda to extract the cation and anion from a species formula
    auto getIonPair = [](std::string str) -> std::pair<std::string, std::string>
    {
        if (str.size() < 2) throw std::domain_error("String is too short.");
        std::string cat;
        cat.push_back(str[0]);
        for (size_t i = 1; i < str.size(); i++){
            if (islower(str[i])) cat += str[i];
            else break;
        }

        std::string an = "";
        for (int i = str.size()-1; i >= 0; i--){
            if (islower(str[i])) an = str[i] + an;
            else if (isupper(str[i])){
                an = str[i] + an;
                break;
            }
        }

        return std::pair<std::string, std::string>{cat, an};
    };

    auto compVector = [](std::pair<std::string, double> p,
                         std::pair<std::string, double> q)
                        { return p.second > q.second; };

    using mapStrDbl = std::unordered_map<std::string, double>;
    using pairStrDbl = std::pair<std::string, double>;

    for (size_t i = 0; i < m; i++){
        mapStrDbl specMap; // To access and modify values
        std::vector<pairStrDbl> specSol, specPair, specGas; // vectors to be sorted
        double T, P;

        std::vector<mapStrDbl> surrElemMaps;
        surrElemMaps.reserve(6);
        std::unordered_map<std::string, size_t> mapIndex;
        // To indicate the group index in surrElemMaps

        // Populate surrElemMap with the mole fraction of each surrogated
        // element in each group.
        for (auto it: surrogateMap){
            mapStrDbl concmp;
            double sumSurr = 0.0;
            static size_t index = 0;

            for (const std::string& ele: it.second){
                if (scaleData.find(ele) != scaleData.end()){
                    // If the element is found in scaleData
                    concmp[ele] = scaleData.at(ele)[i];
                    sumSurr += concmp.at(ele);
                } else concmp[ele] = 0;
            }

            if (sumSurr > 0.0){
                for (auto& elePair: concmp) elePair.second /= sumSurr;
                surrElemMaps.push_back(concmp);
                std::string surrogate = it.first;
                mapIndex[surrogate] = index;
                index++;
                if (index == surrogateMap.size()) index = 0;
            }
        }

        // Similar, but for gases
        mapStrDbl concmp;
        double sumGas = 0.0;
        for (const std::string& ele: gases){
            if (scaleData.find(ele) != scaleData.end()){
                if (ele == "H"){
                    concmp["H2"] = scaleData.at(ele)[i]/2;
                    sumGas += concmp.at("H2");
                } else{
                    concmp[ele] = scaleData.at(ele)[i];
                    sumGas += concmp.at(ele);
                }
            } else concmp[ele] = 0;
        }
        if (sumGas > 0.0){
            specGas.reserve(gases.size());
            for (auto& elePair: concmp){
                if (elePair.second > FRACTION_CUTOFF)
                    specGas.push_back(pairStrDbl(elePair.first, elePair.second));
            }
        }

        while (line.find("Moles of pairs") >= line.size()) getline(input,line);

        auto v = strToVect(line);
        double n = std::stod(v[0]); // Total amount (moles) of salts
        double xUF4 = 0;
        double xUI4 = 0;

        while (line.find("Quadruplet") >= line.size()){
            getline(input,line);
            if (!containsNumber(line)) continue;

            auto v = strToVect(line);
            if (v.back() == "}") v.pop_back();
            double x = std::stod(v[v.size()-2]);

            if (v.back() == "U2F8") xUF4 += x*2;
            else if (v.back() == "U2I8") xUI4 += x*2;
            else if (v.back() == "U[VI]-F4" || v.back() == "U[VII]-F4")
                xUF4 += x;
            else if (v.back() == "U[VI]-I4" || v.back() == "U[VII]-I4")
                xUI4 += x;
            else if (v.back() == "Be2F4") specMap["BeF2"] += x*2;
            else if (v.back() == "Be2I4") specMap["BeI2"] += x*2;
            else specMap[v.back()] = x;
        }

        if (xUF4 > 0) specMap["UF4"] = xUF4; // Sum of all configurations
        if (xUI4 > 0) specMap["UI4"] = xUI4;

        while (line.find("Temperature") >= line.size()) getline(input,line);
        v = strToVect(line);
        T = std::stod(v[v.size()-2]); // System temperature

        while (line.find("Pressure") >= line.size()) getline(input,line);
        v = strToVect(line);
        P = std::stod(v[v.size()-2]);
        P *= 1.01325; // System pressure in bar

        double sumSol = 0.0;
        specSol.reserve(metals.size());
        for (std::string m: metals){
            try{
                double amount = scaleData.at(m)[i];
                if (amount >= FRACTION_CUTOFF){
                    sumSol += amount;
                    specSol.push_back(pairStrDbl(m, amount));
                }
            } catch (const std::out_of_range& ex){}
        }

        // Decouples data in the map and put the results in the vector
        for (auto& it: specMap){
            if (it.second <= FRACTION_CUTOFF) continue;
            it.second *= n;

            std::pair<std::string, std::string> ions;
            try{
                ions =  getIonPair(it.first);
            } catch (const std::domain_error& ex){
                continue;
            }

            if (includesSS && ions.first == "Ni"){}

            else if (ions.second == "I"){
                if (ions.first == "Ca" || ions.first == "La" ||
                    ions.first == "Pu" || ions.first == "Th"){

                    char oxiState;
                    if (ions.first == "Ca"){ oxiState = '2'; }
                    else if (ions.first == "La"){ oxiState = '3'; }
                    else if (ions.first == "Pu"){ oxiState = '3'; }
                    else{ oxiState = '4'; }


                    size_t catIndex = mapIndex[ions.first];
                    size_t anIndex = mapIndex["I"];
                    for (auto cat: surrElemMaps[catIndex]){
                        for (auto an: surrElemMaps[anIndex]){
                            std::string specName = cat.first + an.first + oxiState;
                            double amount = it.second * surrElemMaps[catIndex][cat.first]
                                                      * surrElemMaps[anIndex][an.first];
                            if (amount == 0) continue;
                            specPair.push_back(pairStrDbl{specName, amount});
                        }
                    }

                } else{
                    char oxiState = it.first.back() == 'I' ? ' ' : it.first.back();
                    size_t anIndex = mapIndex["I"];
                    for (auto an: surrElemMaps[anIndex]){
                        std::string specName = ions.first + an.first + oxiState;
                        double amount = it.second * surrElemMaps[anIndex][an.first];
                        if (amount == 0) continue;
                        specPair.push_back(pairStrDbl{specName, amount});
                    }
                }

            } else if (ions.first == "Ca" || ions.first == "La" ||
                       ions.first == "Pu" || ions.first == "Th"){

                char oxiState;
                if (ions.first == "Ca"){ oxiState = '2'; }
                else if (ions.first == "La"){ oxiState = '3'; }
                else if (ions.first == "Pu"){ oxiState = '3'; }
                else{ oxiState = '4'; }

                size_t catIndex = mapIndex[ions.first];
                for (auto cat: surrElemMaps[catIndex]){
                    std::string specName = cat.first + ions.second + oxiState; // is ions.second just 'F'?
                    double amount = it.second * surrElemMaps[catIndex][cat.first];
                    if (amount == 0) continue;
                    specPair.push_back(pairStrDbl{specName, amount});
                }
            } else if (!(includesSS && ions.first == "Ni")){
                    specPair.push_back(pairStrDbl{it.first, it.second});
            }
        }

        double sumSalt = 0.0;
        std::string message;
        for (auto it: specPair) sumSalt += it.second;

        if (includesSS || includesFP){
            auto del = [](const Vector& zeta) -> double
            {
                if (zeta.n() != 18) throw std::domain_error("Zeta's size is incorrect");
                const Vector coefficient{1,2,1,2,1,2,1,2,5,0,0,-1,-1,-1,-1,4,1,1};
                return coefficient*zeta;
            };

        // Calculates the amount of fluorides from structural metals and/or fission products
            std::function<Vector(const Vector&)> thermoFunc =
            [&](const Vector& zeta) -> Vector
            {
                Vector y(zeta);
                double nSol = sumSol-zeta(1)-zeta(3)-zeta(5)-zeta(7)-zeta(8)-zeta(15);
                double nSalt = sumSalt+zeta(1)+zeta(3)+zeta(5)+zeta(7);
                double nGas = sumGas;
                if (nGas > 0.0) nGas += zeta(0)/2+zeta(8)-zeta(9)-zeta(10)+zeta(15);

                double xUF3 = (specMap["UF3"]+del(zeta))/nSalt;
                double xUF4 = (specMap["UF4"]-del(zeta))/nSalt;
                double GF2 = G_F(xUF3, xUF4, T);

                double pMoF4, pMoF5, pMoF6;
                // pMoF5 is expected to be tiny and that leads to loss of precision?
                double aCrF2 = gamma_Inf_CrF2*(zeta(1)-zeta(2))/nSalt;
                double aCrF3 = 1.0*zeta(2)/nSalt;
                double aFeF2 = gamma_Inf_FeF2*(zeta(3)-zeta(4))/nSalt;
                double aFeF3 = 1.0*zeta(4)/nSalt;
                double aCoF2 = 1.0*(zeta(5)-zeta(6))/nSalt;
                double aCoF3 = 1.0*zeta(6)/nSalt;
                double aNiF2 = gamma_Inf_NiF2(specMap["LiF"]/sumSalt)*zeta(7)/nSalt;

                double XCr, XFe, XCo, XNi, XMo;
                for (auto it = specSol.cbegin(); it != specSol.cend(); ++it){
                    if (it->first == "Cr") XCr = (it->second-zeta(1))/nSol;
                    else if (it->first == "Fe") XFe = (it->second-zeta(3))/nSol;
                    else if (it->first == "Co") XCo = (it->second-zeta(5))/nSol;
                    else if (it->first == "Ni") XNi = (it->second-zeta(7))/nSol;
                    else if (it->first == "Mo") XMo = (it->second-zeta(15))/nSol;
                }

                if (sumGas > 0.0 && includesFP){
                    double nH2 = scaleData.at("H")[i]/2;
                    if (nH2 > 0.0){
                        nGas += zeta(0)/2;
                        double pH2 = (nH2-zeta(0)/2)/nGas*P;
                        double pHF = zeta(0)/nGas*P;
                        y(0) = pHF*pHF - pH2*exp((GF2 - 2*G_HF(T))/(R*T));
                    }

                    double nNb = scaleData.at("Nb")[i];
                    if (nNb > 0.0){
                        double XNb = (nNb-zeta(8))/nSol;
                        double pNbF5 = (zeta(8)-2*zeta(9)-zeta(10)-zeta(11))/nGas*P;
                        double pNb2F10 = (zeta(9)-zeta(10))/nGas*P;
                        double pNb3F15 = zeta(10)/nGas*P;
                        double pNbF4 = (zeta(11)-zeta(12))/nGas*P;
                        double pNbF3 = (zeta(12)-zeta(13))/nGas*P;
                        double pNbF2 = (zeta(13)-zeta(14))/nGas*P;
                        double pNbF = zeta(14)/nGas*P;

                        y(8) = pNbF5 - XNb*exp((2.5*GF2 - G_NbF5(T))/(R*T));
                        y(9) = pNb2F10 - pNbF5*pNbF5*exp(-G_Nb2F10(T)/(R*T));
                        y(10) = pNb3F15 - pNbF5*pNb2F10*exp(-G_Nb3F15(T)/(R*T));
                        y(11) = pNbF4 - pNbF5*exp((-0.5*GF2 + G_NbF4(T))/(R*T));
                        y(12) = pNbF3 - pNbF4*exp((-0.5*GF2 + G_NbF3(T))/(R*T));
                        y(13) = pNbF2 - pNbF3*exp((-0.5*GF2 + G_NbF2(T))/(R*T));
                        y(14) = pNbF - pNbF2*exp((-0.5*GF2 + G_NbF(T))/(R*T));
                    }
                }

                if (includesSS){
                    if (XCr > FRACTION_CUTOFF){
                        y(1) = aCrF2 - XCr*exp((GF2 - G_CrF2(T))/(R*T));
                        y(2) = aCrF3 - aCrF2*exp((0.5*GF2 - G_CrF3(T))/(R*T));
                    } if (XFe > FRACTION_CUTOFF){
                        y(3) = aFeF2 - XFe*exp((GF2 - G_FeF2(T))/(R*T));
                        y(4) = aFeF3 - aFeF2*exp((0.5*GF2 - G_FeF3(T))/(R*T));
                    } if (XCo > FRACTION_CUTOFF){
                        y(5) = aCoF2 - XCo*exp((GF2 - G_CoF2(T))/(R*T));
                        y(6) = aCoF3 - aCoF2*exp((0.5*GF2 - G_CoF3(T))/(R*T));
                    } if (XNi > FRACTION_CUTOFF){
                        y(7) = aNiF2 - XNi*exp((GF2 - G_NiF2(T))/(R*T));
                    } if (XMo > FRACTION_CUTOFF && sumGas > 0.0){
                        pMoF4 = (zeta(15)-zeta(16))/nGas*P;
                        pMoF5 = (zeta(16)-zeta(17))/nGas*P;
                        pMoF6 = zeta(17)/nGas*P;
                        y(15) = pMoF4 - XMo*exp((2*GF2 - G_MoF4(T))/(R*T));
                        y(16) = pMoF5 - pMoF4*exp((0.5*GF2 - G_MoF5(T))/(R*T));
                        y(17) = pMoF6 - pMoF5*exp((0.5*GF2 - G_MoF6(T))/(R*T));
                    }
                }

                return y;
            };

            Vector zeta{1e-6*std::max(sumGas,1e-10), // HF
                        1e-4*sumSalt, // CrF2
                        1e-9*sumSalt, // CrF3
                        1e-9*sumSalt, // FeF2
                        1e-10*sumSalt, // FeF3
                        1e-12*sumSalt, // CoF2
                        1e-30*sumSalt, // CoF3
                        1e-12*sumSalt, // NiF2
                        1e-20*std::max(sumGas,1e-10), //NbF5
                        1e-30*std::max(sumGas,1e-10), //Nb2F10
                        1e-40*std::max(sumGas,1e-10), //Nb3F15
                        1e-22*std::max(sumGas,1e-10), //NbF4
                        1e-25*std::max(sumGas,1e-10), //NbF3
                        1e-30*std::max(sumGas,1e-10), //NbF2
                        1e-35*std::max(sumGas,1e-10), //NbF
                        1e-30*std::max(sumGas,1e-10), // MoF4
                        1e-35*std::max(sumGas,1e-10), // MoF5
                        1e-40*std::max(sumGas,1e-10)}; // MoF6
            /* sumGas being 0 breaks the code, so if I don't want to solve for
            the gas phase it will be taken cared of by the bool includesHF. */

            try{
                zeta = newton(thermoFunc, zeta, 100, 1e-6);
                for (size_t j = 0; j < zeta.n(); j++){ // Write an iterator function for Vector?
                    if (zeta(j) <= FRACTION_CUTOFF) zeta(j) = 0.0;
                }

                sumSol -= (zeta(1)+zeta(3)+zeta(5)+zeta(7)+zeta(8)+zeta(15));
                sumSalt += zeta(1)+zeta(3)+zeta(5)+zeta(7);
                sumGas += zeta(0)/2+zeta(8)-zeta(9)-zeta(10)+zeta(15);
                specPair.reserve(specPair.size()+7);
                specGas.reserve(specGas.size()+11);

                if (includesFP){
                    if (sumGas > 0.0){
                        if (zeta(0)/sumGas >= FRACTION_CUTOFF)
                            specGas.push_back({"HF", zeta(0)});
                        for (auto it = specGas.begin(); it != specGas.end(); ++it){
                            if (it->first == "H2"){
                                double frac = (it->second-zeta(0)/2)/sumGas;
                                if (frac > FRACTION_CUTOFF) it->second -= zeta(0)/2;
                                else specGas.erase(it);
                                break;
                            }
                        }
                    }

                    if (sumSol > 0.0){
                        if ((zeta(8)-2*zeta(9)-zeta(10)-zeta(11))/sumGas >= FRACTION_CUTOFF)
                            specGas.push_back({"NbF5", zeta(8)-2*zeta(9)-zeta(10)-zeta(11)});
                        if ((zeta(9)-zeta(10))/sumGas >= FRACTION_CUTOFF)
                            specGas.push_back({"Nb2F10", zeta(9)-zeta(10)});
                        if (zeta(10)/sumGas >= FRACTION_CUTOFF)
                            specGas.push_back({"Nb3F15", zeta(10)});
                        if ((zeta(11)-zeta(12))/sumGas >= FRACTION_CUTOFF)
                            specGas.push_back({"NbF4", zeta(11)-zeta(12)});
                        if ((zeta(12)-zeta(13))/sumGas >= FRACTION_CUTOFF)
                            specGas.push_back({"NbF3", zeta(12)-zeta(13)});
                        if ((zeta(13)-zeta(14))/sumGas >= FRACTION_CUTOFF)
                            specGas.push_back({"NbF2", zeta(13)-zeta(14)});
                        if (zeta(14)/sumGas >= FRACTION_CUTOFF)
                            specGas.push_back({"NbF", zeta(14)});
                        for (auto it = specSol.begin(); it != specSol.end(); ++it){
                            if (it->first == "Nb"){
                                double frac = (it->second-zeta(8))/sumSol;
                                if (frac > FRACTION_CUTOFF) it->second -= zeta(8);
                                else specSol.erase(it);
                                break;
                            }
                        }
                    }
                }

                if(includesSS && sumSol > 0.0){
                    for (auto it = specSol.begin(); it != specSol.end(); ++it){
                        if (it->first == "Cr"){
                            double frac = (it->second-zeta(1))/sumSol;
                            if (frac >= FRACTION_CUTOFF){
                                it->second -= zeta(1);
                                if ((zeta(1)-zeta(2))/sumSol >= FRACTION_CUTOFF)
                                    specPair.push_back({"CrF2", zeta(1)-zeta(2)});
                                if (zeta(2)/sumSol >= FRACTION_CUTOFF)
                                    specPair.push_back({"CrF3", zeta(2)});
                            } else specSol.erase(it);
                        } else if (it->first == "Fe"){
                            double frac = (it->second-zeta(3))/sumSol;
                            if (frac >= FRACTION_CUTOFF){
                                it->second -= zeta(3);
                                if ((zeta(3)-zeta(4))/sumSol >= FRACTION_CUTOFF)
                                    specPair.push_back({"FeF2", zeta(3)-zeta(4)});
                                if (zeta(4)/sumSol >= FRACTION_CUTOFF)
                                    specPair.push_back({"FeF3", zeta(4)});
                            } else specSol.erase(it);
                        } else if (it->first == "Co"){
                            double frac = (it->second-zeta(5))/sumSol;
                            if (frac >= FRACTION_CUTOFF){
                                it->second -= zeta(5);
                                if ((zeta(5)-zeta(6))/sumSol >= FRACTION_CUTOFF)
                                    specPair.push_back({"CoF2", zeta(5)-zeta(6)});
                                if (zeta(6)/sumSol >= FRACTION_CUTOFF)
                                    specPair.push_back({"CoF3", zeta(6)});
                            } else specSol.erase(it);
                        } else if (it->first == "Ni"){
                            double frac = (it->second-zeta(7))/sumSol;
                            if (frac >= FRACTION_CUTOFF){
                                it->second -= zeta(7);
                                if (zeta(7)/sumSol >= FRACTION_CUTOFF)
                                    specPair.push_back({"NiF2", zeta(7)});
                            } else specSol.erase(it);
                        } else if (it->first == "Mo"){
                            double frac = (it->second-zeta(15))/sumSol;
                            if (frac >= FRACTION_CUTOFF){
                                it->second -= zeta(15);
                                if ((zeta(15)-zeta(16))/sumGas >= FRACTION_CUTOFF)
                                    specGas.push_back({"MoF4", zeta(15)-zeta(16)});
                                if ((zeta(16)-zeta(17))/sumGas >= FRACTION_CUTOFF)
                                    specGas.push_back({"MoF5", zeta(16)-zeta(17)});
                                if (zeta(17)/sumGas >= FRACTION_CUTOFF)
                                    specGas.push_back({"MoF6", zeta(17)});
                            } else specSol.erase(it);
    //                           specGas.push_back(pairStrDbl{"MoF5", zeta(9)});
    //                           specGas.push_back(pairStrDbl{"MoF6", zeta(8)});
                        }
                    }
                }

                for (auto& it: specPair){
                    if (it.first == "UF3") it.second += del(zeta);
                    if (it.first == "UF4") it.second -= del(zeta);
                }

            } catch(const std::exception& ex){
                message += std::string(ex.what()) + "\n";
                message += "WARNING: Unable to solve for zeta.\n";
            }
        }

        // Sorts the vector and outputs its data
        if (sumSol > 0.0){
            std::sort(specSol.begin(), specSol.end(), compVector);
            output << sumSol << " Moles of solid solution\n";

            for (auto it = specSol.cbegin(); it != specSol.cend(); ++it){
                std::string start = it == specSol.cbegin() ? "\t{" : "\t+";
                if (it->second >= 1.0E-4){
                    output << start << " " << std::defaultfloat << it->second/sumSol;
                    if (it->second >= 1.0E-2) output << "\t";
                } else output << start << " " << std::scientific << it->second/sumSol;
                output << "\t\t" << it->first;
                if (*it == specSol.back()) output << "\t}";
                output << "\n";
            }
        }

        if (sumSalt > 0.0){
            std::sort(specPair.begin(), specPair.end(), compVector);
            output << std::defaultfloat << sumSalt << " Moles of pairs\n";
            output << "Pair fractions:\n";

            for (auto it = specPair.cbegin(); it != specPair.cend(); ++it){
                std::string start = it == specPair.cbegin() ? "\t{" : "\t+";
                double frac = it->second/sumSalt;
                if (frac >= 1.0E-4){
                    output << start << " " << std::defaultfloat << frac;
                    if (frac >= 1.0E-2) output << "\t";
                } else output << start << " " << std::scientific << frac;
                output << "\t\t" << it->first;
                if (*it == specPair.back()) output << "\t}";
                output << "\n";
            }
        }

        if (sumGas > 0.0){
            std::sort(specGas.begin(), specGas.end(), compVector);
            output << std::defaultfloat << sumGas << " mol gas_ideal\n";

            for (auto it = specGas.cbegin(); it != specGas.cend(); ++it){
                std::string start = it == specGas.cbegin() ? "\t{" : "\t+";
                double frac = it->second/sumGas;
                if (frac >= 1.0E-4){
                    output << start << " " << std::defaultfloat << frac;
                    if (frac >= 1.0E-2) output << "\t";
                } else output << start << " " << std::scientific << frac;
                output << "\t\t" << it->first;
                if (*it == specGas.back()) output << "\t}";
                output << "\n";
            }
        }

        output << "Temperature: T = " << std::defaultfloat << T << " K.\n";
        output << "Pressure: P = " << P/1.01325 << " atm.\n";

        if (message.empty()) message = "DEBUG: Successful exit.\n\n";
        else message += "WARNING: Erroneous calculations may have occurred.\n\n";
        output << message;
    }

    input.close();
    output.close();
}
