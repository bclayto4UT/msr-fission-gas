#include "miscellaneous.h"
#include <stdexcept>

std::string elementSymb(const std::string& str){
    std::string ele;
    ele.reserve(str.size());
    ele += toupper(str.front());
//    for (size_t i = 1; i < str.size(); i++) ele[i] += tolower(str[i]);
    for (auto it = str.begin()+1; it != str.end(); ++it) ele += tolower(*it);
    return ele;
}

strVect strToVect (std::string& data, char deli)
{
    strVect v;
    while (!data.empty()){

        while (data.find(deli) == 0 || data.find("\t") == 0) data.erase(0, 1);

        while (data.find(deli) == data.size() || data.find("\t") == data.size()){
            data.erase(data.size()-1, data.size());
        }


        if (!data.empty()){
            int pos1 = data.find(deli);
            int pos2 = data.find("\t");

            int index;
            if (pos1 > 0){
                if (pos2 > 0) index = std::min(pos1, pos2);
                else index = pos1;
            } else{
                if (pos2 > 0) index = pos2;
                else index = -1;
            }

            if (index == -1){
                v.push_back(data);
                return v;
            } else{
                auto val = data.substr(0, index);
                v.push_back(val);
                data.erase(0, index);
            }
        }
    }

    return v;
}

std::vector<double> strToVectDouble (const std::string& str){
    std::vector<double> vd;

    std::string s = str;
    strVect vs = strToVect(s); // Delimeter is ' '

    try{
        if (vs.size() == 1){ // There is no ' '
            if (isNumeric(vs[0])){ // Checks if the only item is already a double;
                vd.push_back(std::stod(vs[0]));
                return vd;
            } else throw std::invalid_argument("");
        } else{
            for (auto i: vs){
                vd.reserve(vs.size());
                vd.push_back(std::stod(i));
            }
            return vd;
        }

    } catch(const std::invalid_argument&){
        try{
            vd = std::vector<double>{};
            std::string s = str;
            strVect vs = strToVect(s, ':'); // Delimeter is ':'
            if (vs.size() == 1 || vs.size() > 3) throw std::invalid_argument("");
            else if (vs.size() == 2) throw std::domain_error("");
            // Throws a different exception so that the function can pass in another parameter

            else{
                double start = std::stod(vs[0]);
                double stop  = std::stod(vs[1]);

                if (start == stop){
                    vd.push_back(start);
                    return vd;
                }

                double step  = std::stod(vs[2]);
                if (start > stop && step > 0) step *= -1;
                int m = ceil((stop-start)/step);
                vd.reserve(m);
                for (int i = 0; i < m; i++) vd.push_back(start+step*i);
                return vd;
            }

        } catch (const std::invalid_argument&){
            vd = std::vector<double>{};
            std::string s = str;
            strVect vs = strToVect(s, ','); // Delimeter is ','
//            for (std::vector<double>::size_type i = 0; i < vs.size(); i++){
//                vd.reserve(vs.size());
//                if (isNumeric(vs[i])) vd.push_back(std::stod(vs[i]));
//                else throw std::invalid_argument("There may be invalid characters.");
//            }

            for (auto it = vs.cbegin(); it != vs.cend(); ++it){
                vd.reserve(vs.size());
                if (isNumeric(*it)) vd.push_back(std::stod(*it));
                else throw std::invalid_argument("There may be invalid characters.");
            }
            return vd;
        }
    }
}

bool isNumeric (const std::string& str) noexcept
{
    if (str.size() == 0) return false;

    bool hasDec = false;
    bool hasExp = false;
    bool hasNeg = str[0] == '-';
    bool hasPos = str[0] == '+';

    size_t start = (hasNeg || hasPos) ? 1 : 0;
    for (auto it = str.begin()+start; it != str.end(); ++it){
        char c = *it;
        if (!isdigit(c)){
            if (c == '.'){
                if (hasDec) return false;
                else hasDec = true;
            } else if (c == '-'){
                if (hasNeg) return false;
                else hasNeg = true;
            } else if (c == '+'){
                if (hasPos) return false;
                else hasPos = true;
            } else if (c == 'E' || c == 'e'){
                if (hasExp || it-str.begin() == str.size()-1) return false;
                else{
                    hasExp = true;
                    if (*(it+1) == '-') hasNeg = false;
                    else if (*(it+1) == '+') hasPos = false;
                }
            } else return false;
        }
    }
    return true;
}

bool containsNumber (const std::string& str) noexcept
{
    for (char const &c : str){
        if (std::isdigit(c) == 1) return true;
    }
    return false;
}

void convertibleNum(std::string& str){
    bool throwsError = true;
    while (throwsError){
        try{
            std::stod(str);
            return;
        } catch(const std::exception& ex){
            str.erase(0,1);
        }
    }
}

