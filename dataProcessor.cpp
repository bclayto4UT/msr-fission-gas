#include <algorithm>
#include <iostream>
#include <fstream>

#include "dataProcessor.h"

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
                elementMap.at(ele); // To check if ele is a chemical element
                vect[0].push_back(ele);
            } catch (std::out_of_range ex){
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
            output << elementMap.at(dataVect[0][j]);
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
                  const std::string& thermoOut)
/****************************************************************
Input(s):
    scaleData: the 2-D vector (table) containing the amount of
               each species as a function of time (can be obtained
               using the scaleToVector function)
    thermoRes: the name of the file containing Thermochimica re-
               results that have been ran on data with surrogates
    thermoOut: the output file name
Output: none
    New Thermochimica results are written in the output file,
    in which each salt containing a surrogate metal or nonmetal
    has been broken further down to its original counterpart(s),
    with the ratio of the corresponding elements preserved.
Warning:
    only a prediction, works best on hypo-stoichiometric systems.
*****************************************************************/
{

}

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
}
