#include <algorithm>
#include <iostream>
#include <fstream>

#include "dataProcessor.h"
#include "iterativeNL.h"
#include "rootFinding.h"
#include "thermoElectroChem.h"

std::vector<strVect> scaleToVector(const std::string& inFile,
                                   bool timeOnRows)
/****************************************************************
Input(s):
    inFile: the name of the file with the SCALE output to be read
    timeOnRows: if the orientation of the file is such that the
                elements are on the columns and the time intervals
                are on the rows (DEFAULT: true).
Output:
    a 2D vector of strings that tabulates the data such that the
    elements are on the columns and the time intervals are on the
    rows.
*****************************************************************/
{
    // Opens input file
    std::ifstream input(inFile);
    std::string line;
    strVect v;
    std::vector<strVect> vect;

    if (input.is_open()) getline(input,line);
    else throw std::invalid_argument("Cannot open input file.");

    if (timeOnRows){
        v = strToVect(line);
        for (int i = 0; i < v.size(); i++) v[i][0] = toupper(v[i][0]);
        vect.push_back(v);

         while(getline(input, line)){
            v = strToVect(line);
            v.erase(v.cbegin());
            // Consider changing into a deque for the more efficient pop_front function
            vect.push_back(v);
         }

    } else{
        // Counts the days
        std::size_t n;
        n = strToVect(line).size();
        for (int i = 0; i <= n; i++) vect.push_back(v);

        // Inputs data into vector
        while(getline(input, line)){
            // Converts element symbol to atomic number and stores it
            auto pos = std::min(line.find("\t"), line.find(" "));
            while (pos == 0){
                line.erase(0, 1);
                pos = std::min(line.find("\t"), line.find(" "));
            }
            std::string ele = line.substr(0, pos);
            ele[0] = toupper(ele[0]);
            line.erase(0, pos);

            try{
                atomNumMap.at(ele); // To check if ele is a chemical element
                vect[0].push_back(ele);
            } catch (std::out_of_range& ex){
                break;
            }

            // Stores the rest of the data
            v = strToVect(line);
            for (int i = 0; i < v.size(); i++){
                vect[i+1].push_back(v[i]);
            }
        }
    }

    input.close();
    return vect;
}

void vectToThermI(const std::vector<strVect>& dataVect,
                  const std::string& outFile)
/****************************************************************
Input(s):
    dataVect: the tabulated data in the form of a 2D vector
    outFile: the output file name
Output: none
    a text file formatted as a Thermochimica input is created
    if there are multiple rows, multiple files are created
*****************************************************************/
{
    std::ofstream output;
    auto pos = outFile.find(".");
    std::string fileEnding;
    std::string fileStem;
    std::string fileName;

    if (pos < outFile.size()){
        fileEnding = outFile.substr(pos, outFile.size());
        fileStem = outFile.substr(0, pos);
    } else{
        fileStem = outFile;
    }

    for (int i = 1; i < dataVect.size(); i++){
        fileName = fileStem + std::to_string(i) + fileEnding;
        output.open(fileName);
        for (int j = 0; j < dataVect[i].size(); j++){
            output << "dElementMass(";
            output << atomNumMap.at(dataVect[0][j]);
            output << ") = ";
            output << dataVect[i][j];
            output << "\n";
        }
        output.close();
    }

}


void textToExcel(const std::string& inFile,
                 const std::string& outFile,
                 std::string& dataType)
/****************************************************************
Input(s):
    inFile: the name of the input file that contains a
            Thermochimica output
    outFile: the output file name
    dataType: a string that contains the types of data that the
              user wishes to extract. It can include as many
              types as necessary. The accepted types include:
     * ni: the total amount (moles) of ions
     * nx: the total amount (moles) of salt
     * ny: the total amount (moles) of gas
     * x_ABC: the mole fraction of species ABC in the salt phase
     * y_ABC: the mole fraction of species ABC in the gas phase
       "ABC" can be "all", which will output all salts (not ions)
       and/or all gases
     * T: the system temperature
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
    std::map<std::string, double> vx;
    std::map<std::string, double> vy;
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
        else if (s.substr(0,1) == "y" & !has_y_all) vy[s.erase(0,2)] = 0;
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

                auto isFilled = [](const std::map<std::string, double>& m) -> bool
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
    std::map<std::string, double> mx, my;
    std::map<std::string, double>* mapRef = nullptr;

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

void decoupleSurr(const std::vector<strVect>& scaleData,
                  const std::string& thermoRes,
                  const std::string& thermoOut,
                  const bool includesSSol,
                  const bool includesHF)
/****************************************************************
Input(s):
    scaleData:    the 2-D vector (table) containing the amount of
                  each species as a function of time (can be ob-
                  tained using the scaleToVector function)
    thermoRes:    the name of the file containing Thermochimica re-
                  results that have been ran on data with surrogates
    thermoOut:    the output file name
    includesSSol: whether the solid solution of the reactor struc-
                  ture is considered (DEFAULT: false).
                  The input file must contain the mole fraction of
                  Cr, Fe, and Ni in the solid solution (not their
                  amount in moles).
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

    std::vector<std::vector<std::string>> surrElem; // List of surrogated elements
    surrElem.push_back(strVect{"Ca", "Ba", "Sr"});
    surrElem.push_back(strVect{"La", "Pr", "Pm", "Sm", "Eu", "Gd", "Tb",
                               "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Y"});
    surrElem.push_back(strVect{"Th", "Zr", "Hf"});
    surrElem.push_back(strVect{"I", "Br"});
    surrElem.push_back(strVect{"H", "He", "Ne", "Ar", "Kr", "Xe"});

    std::vector<std::map<std::string, double>> surrElemMaps(surrElem.size());

    // lambda to extract the cation and anion from a species formula
    auto getIonPair = [](std::string str) -> std::pair<std::string, std::string>
    {
        if (str.size() < 2) throw std::domain_error("String is too short.");
        std::string cat;
        cat.push_back(str[0]);
        for (int i = 1; i < str.size(); i++){
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
    for (int i = 1; i < scaleData.size(); i++){
        std::map<std::string, double> specMap; // To access and modify values
        std::vector<std::pair<std::string, double>> specPair, specGas;
        // specPair and specGas are vectors so that they can be sorted
        double T, P, XCr, XFe, XNi;
        bool validSolids = false;

        // Populate surrElemMap with the mole fraction of each surrogated
        // element in each group.
        for (int group = 0; group < surrElem.size(); group++){
            std::map<std::string, double> m;
            double sumSurr = 0;

            for (int j = 0; j < surrElem[group].size(); j++){
                std::string ele = surrElem[group][j];
                auto it = std::find(scaleData[0].cbegin(), scaleData[0].cend(), ele);
                if (it != scaleData[0].cend()){ // If the element is found
                    m[ele] = std::stod(scaleData[i][it - scaleData[0].cbegin()]);
                    if (ele == "H") m[ele] /= 2;
                    sumSurr += m.at(ele);
                } else m[ele] = 0;
            }

            if (sumSurr > 0 && group < 4){
                for (auto& it: m) it.second /= sumSurr;
                surrElemMaps[group] = m;
            } else if (group == 4){
                // Group 4 is a group of gases;
                for (auto it: m){
                    if (it.second == 0) continue;
                    std::pair<std::string, double> spec{it.first, it.second};
                    if (it.first == "H") spec.first = "H2";
                    specGas.push_back(spec);
                }
            } // Group 5 does not count
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

        if (includesHF || includesSSol){
            while (line.find("Temperature") >= line.size()) getline(input,line);
            v = strToVect(line);
            T = std::stod(v[v.size()-2]); // System temperature
        }

        if (includesHF){
            while (line.find("Pressure") >= line.size()) getline(input,line);
            v = strToVect(line);
            P = std::stod(v[v.size()-2]);
            P *= 1.01325;
            // System pressure in bar
        }

        if (includesSSol){
            auto it = std::find(scaleData[0].cbegin(), scaleData[0].cend(), "Cr");
            if (it != scaleData[0].cend()){
                XCr = std::stod(scaleData[i][it - scaleData[0].cbegin()]);
            }
            it = std::find(scaleData[0].cbegin(), scaleData[0].cend(), "Fe");
            if (it != scaleData[0].cend()){
                XFe = std::stod(scaleData[i][it - scaleData[0].cbegin()]);
            }
            it = std::find(scaleData[0].cbegin(), scaleData[0].cend(), "Ni");
            if (it != scaleData[0].cend()){
                XNi = std::stod(scaleData[i][it - scaleData[0].cbegin()]);
            }

            if (XCr >= 0 && XFe >= 0 && XNi >= 0 && XCr + XFe + XNi > 0){
                validSolids = true;
                if (XCr > 1 || XFe > 1 || XNi > 1 || XCr + XFe + XNi > 1){
                    double sum = XCr + XFe + XNi;
                    XCr /= sum;
                    XFe /= sum;
                    XNi /= sum;
                }
            }
        }

        // Decouples data in the map and put the results in the vector
        for (auto& it: specMap){
            if (it.second == 0) continue;
            it.second *= n;

            std::pair<std::string, std::string> ions;
            try{
                ions =  getIonPair(it.first);
            } catch (const std::domain_error& ex){
                continue;
            }

            if (includesSSol & ions.first == "Ni"){
            }

            if (ions.second == "I"){
                if (ions.first == "Ca" || ions.first == "La" ||
                    ions.first == "Th"){

                    char oxiState;
                    if (ions.first == "Ca") oxiState = '2';
                    else if (ions.first == "La") oxiState = '3';
                    else oxiState = '4';
                    int mapIndex = oxiState - '2';

                    for (auto cat: surrElem[mapIndex]){
                        for (auto an: surrElem[3]){
                            std::string specName = cat + an + oxiState;
                            double amount = it.second * surrElemMaps[mapIndex][cat]
                                                      * surrElemMaps[3][an];
                            if (amount == 0) continue;
                            std::pair<std::string, double> spec{specName, amount};
                            specPair.push_back(spec);
                        }
                    }

                } else{
                    char oxiState = it.first.back() == 'I' ? ' ' : it.first.back();
                    if (includesSSol && ions.first == "Ni") it.second *= XNi;
                    for (auto an: surrElem[3]){
                        std::string specName = ions.first + an + oxiState;
                        double amount = it.second * surrElemMaps[3][an];
                        if (amount == 0) continue;
                        std::pair<std::string, double> spec{specName, amount};
                        specPair.push_back(spec);
                    }
                }

            } else if (ions.first == "Ca" || ions.first == "La" ||
                       ions.first == "Th"){

                char oxiState;
                if (ions.first == "Ca") oxiState = '2';
                else if (ions.first == "La") oxiState = '3';
                else oxiState = '4';
                int mapIndex = oxiState - '2';

                for (auto cat: surrElem[mapIndex]){
                    std::string specName = cat + ions.second + oxiState;
                    double amount = it.second * surrElemMaps[mapIndex][cat];
                    if (amount == 0) continue;
                    std::pair<std::string, double> spec{specName, amount};
                    specPair.push_back(spec);
                }
            } else{
                if (!(includesSSol && ions.first == "Ni")){
                    std::pair<std::string, double> spec{it.first, it.second};
                    specPair.push_back(spec);
                }
            }
        }

        double sumSalt = 0.0;
        double sumGas = 0.0;
        std::string message;
        for (auto it: specPair) sumSalt += it.second;
        for (auto it: specGas) sumGas += it.second;

        if (includesSSol){
            // Codes to check for mole fraction in solid
            if (validSolids){
                // Outputs solid phase
                std::vector<std::pair<std::string, double>> specSol;
                specSol.push_back(std::pair<std::string, double>{"Cr", XCr});
                specSol.push_back(std::pair<std::string, double>{"Fe", XFe});
                specSol.push_back(std::pair<std::string, double>{"Ni", XNi});

                std::sort(specGas.begin(), specGas.end(), compVector);

                output << "Solid solution\n";
                for (int i = 0; i < specSol.size(); i++){
                    std::string start = i == 0 ? "\t{" : "\t+";
                    output << start << " " << specSol[i].second;
                    output << "\t\t" << specSol[i].first;
                    if (i == specSol.size()-1) output << "\t}";
                    output << "\n";
                }

                    // Calculates the fraction of corrosion products
                std::function<Vector(const Vector&)> thermoFunc =
                [&](const Vector& zeta) -> Vector
                {
                    double xUF3 = specMap["UF3"]*sumSalt;
                    xUF4 = specMap["UF4"]*sumSalt;

                    double nSalt = sumSalt+zeta(1)+zeta(3)+zeta(5);
                    double aCrF2 = 0.5*(zeta(1)-zeta(2))/nSalt;
                    double aCrF3 = 1.0*zeta(2)/nSalt;
                    double aFeF2 = 1.6*(zeta(3)-zeta(4))/nSalt;
                    double aFeF3 = 1.0*zeta(4)/nSalt;
                    double aNiF2 = 1.0*zeta(5)/nSalt;
                    xUF3 = (xUF3+zeta(0)+2*zeta(1)+zeta(2)+2*zeta(3)+zeta(4)+2*zeta(5))/nSalt;
                    xUF4 = (xUF4-zeta(0)-2*zeta(1)-zeta(2)-2*zeta(3)-zeta(4)-2*zeta(5))/nSalt;
                    double GF2 = G_F(xUF3, xUF4, T);

                    Vector y(6);
                    if (!includesHF || sumGas == 0) y(0) = zeta(0);
                    else{
                        double nH2;
                        for (auto& it: specGas){
                            if (it.first == "H2"){
                                nH2 = it.second;
                                break;
                            }
                        }

                        if (nH2 > 0.0){
                            double nGas = sumGas+zeta(0)/2;
                            double pH2 = (nH2-zeta(0)/2)/nGas*P;
                            double pHF = zeta(0)/nGas*P;
                            y(0) = pHF*pHF/pH2 - exp((GF2 - 2*G_HF(T))/(R*T));
                        } else y(0) = zeta(0);
                    }

                    y(1) = aCrF2/XCr - exp((GF2 - G_CrF2(T))/(R*T));
                    y(2) = aCrF3/aCrF2 - exp((0.5*GF2 - G_CrF3(T))/(R*T));
                    y(3) = aFeF2/XFe - exp((GF2 - G_FeF2(T))/(R*T));
                    y(4) = aFeF3/aFeF2 - exp((0.5*GF2 - G_FeF3(T))/(R*T));
                    y(5) = aNiF2/XNi - exp((GF2 - G_NiF2(T))/(R*T));
                    return y;
                };

                Vector zeta{1e-6*std::max(sumGas,1e-10), // sumGas being 0 breaks the code.
                            1e-4*sumSalt,
                            1e-9*sumSalt,
                            1e-9*sumSalt,
                            1e-12*sumSalt,
                            1e-12*sumSalt};

                try{
                    zeta = newton(thermoFunc, zeta, 100, 1e-6);
                    sumSalt += zeta(1) + zeta(3);
                    specPair.push_back(std::pair<std::string, double>{"CrF2", zeta(1)-zeta(2)});
                    specPair.push_back(std::pair<std::string, double>{"CrF3", zeta(2)});
                    specPair.push_back(std::pair<std::string, double>{"FeF2", zeta(3)-zeta(4)});
                    specPair.push_back(std::pair<std::string, double>{"FeF3", zeta(4)});
                    specPair.push_back(std::pair<std::string, double>{"NiF2", zeta(5)});
                    for (auto& it: specPair){
                        if (it.first == "UF3") it.second += (zeta(0)+2*zeta(1)+zeta(2)+2*zeta(3)+zeta(4)+2*zeta(5));
                        if (it.first == "UF4") it.second -= (zeta(0)+2*zeta(1)+zeta(2)+2*zeta(3)+zeta(4)+2*zeta(5));
                    }

                    if (includesHF && sumGas > 0){
                        sumGas += zeta(0)/2;
                        specGas.push_back(std::pair<std::string, double>{"HF", zeta(0)});
                        for (auto& it: specGas){
                            if (it.first == "H2") it.second -= zeta(0)/2;
                        }
                    }
                } catch(const std::exception& ex){
                    zeta = 0.0;
                    message += "WARNING: Unable to solve for zeta.\n";
                }

        } else{
            message += "WARNING: Invalid solid solution specification at number ";
            message += std::to_string(i);
            message += ".\n";
        }
    }

        else if (includesHF && sumGas > 0){
            double nH2;
            for (auto& it: specGas){
                if (it.first == "H2"){
                    nH2 = it.second;
                    break;
                }
            }

            if (nH2 > 0.0){
                std::function<double(const double)> thermoFunc =
                [&](const double zeta) -> double
                {
                    double xUF3 = specMap["UF3"]*sumSalt + zeta;
                    xUF4 = specMap["UF4"]*sumSalt + zeta;
                    double GF2 = G_F(xUF3, xUF4, T);

                    double nGas = sumGas+zeta/2;
                    double pH2 = (nH2-zeta/2)/nGas*P;
                    double pHF = zeta/nGas*P;

                    return pHF*pHF/pH2 - exp((GF2 - 2*G_HF(T))/(R*T));
                };

                try{
                    double zeta = 1e-6*std::max(sumGas, 1e-10);
                    zeta = newton(thermoFunc, zeta, 100, 1e-12);
                    sumGas += zeta/2;
                    specGas.push_back(std::pair<std::string, double>{"HF", zeta});
                    for (auto& it: specGas){
                        if (it.first == "H2"){
                            it.second -= zeta/2;
                            break;
                        }
                    }
                    for (auto& it: specPair){
                        if (it.first == "UF3") it.second += zeta;
                        if (it.first == "UF4") it.second -= zeta;
                    }
                } catch(const std::exception& ex){
                    message += "WARNING: Unable to solve for zeta.\n";
                }
            }
        }

        // Sorts the vector and outputs its data
        std::sort(specPair.begin(), specPair.end(), compVector);
        std::sort(specGas.begin(), specGas.end(), compVector);

        if (sumSalt > 0){
            output << sumSalt << " Moles of pairs\n";
            output << "Pair fractions:\n";
            for (int i = 0; i < specPair.size(); i++){
                std::string start = i == 0 ? "\t{" : "\t+";
                output << start << " " << specPair[i].second/sumSalt;
                output << "\t\t" << specPair[i].first;
                if (i == specPair.size()-1) output << "\t}";
                output << "\n";
            }
        }

        if (sumGas > 0){
            output << sumGas << " mol gas_ideal\n";
            for (int i = 0; i < specGas.size(); i++){
                std::string start = i == 0 ? "\t{" : "\t+";
                output << start << " " << specGas[i].second/sumGas;
                output << "\t\t" << specGas[i].first;
                if (i == specGas.size()-1) output << "\t}";
                output << "\n";
            }
        }

        if (message.empty()) message = "DEBUG: Successful exit.\n\n";
        else message += "WARNING: Erroneous calculations may have occurred.\n\n";
        output << message;
    }

    input.close();
    output.close();
}

/*
void massToMole(const std::string& outFile,
                const std::map<std::string, double>& comp,
                double mass)
{
    std::vector<strVect> v(2);
    for (auto p: comp){
        try{
            auto m = p.second * mass / elementMap.at(p.first);
            v[0].push_back(std::to_string(elementMap.at(p.first)));
            v[1].push_back(std::to_string(m));
        } catch(std::out_of_range& ex){
            throw ex;
        }
    }
    vectToThermI(v, outFile);
} */
