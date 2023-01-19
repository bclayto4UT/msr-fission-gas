#include "miscellaneous.h"

//std::string findElement(int num){
//    if (num < 0 || num > 118)
//        throw std::out_of_range("Atomic number out of range");
//    for (auto const& it: atomNumMap){
//        if (it.second == num) return it.first;
//    }
//}

std::string elementSymb(const std::string& str){
    std::string ele;
    ele += toupper(str[0]);
    for (int i = 1; i < str.size(); i++) ele[i] += tolower(str[i]);
    return ele;
}

strVect strToVect (std::string& data)
{
    strVect v;
    while (!data.empty()){

        while (data.find(" ") == 0 || data.find("\t") == 0){
            data.erase(0, 1);
        }

        while (data.find(" ") == data.size() || data.find("\t") == data.size()){
            data.erase(data.size()-1, data.size());
        }

        if (!data.empty()){
            int index = std::min(data.find(" "), data.find("\t"));
            auto val = data.substr(0, index);
            v.push_back(val);
            data.erase(0, index);
        }
    }

    return v;
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

