#ifndef THERMOELECTROCHEM_H_INCLUDED
#define THERMOELECTROCHEM_H_INCLUDED

#include "miscellaneous.h"

const std::map<std::string, int> atomNumMap = {
    {"H", 1}, {"He", 2}, {"Li", 3}, {"Be", 4}, {"B", 5},
    {"C", 6}, {"N", 7}, {"O", 8}, {"F", 9}, {"Ne", 10},
    {"Na", 11}, {"Mg", 12}, {"Al", 13}, {"Si", 14}, {"P", 15},
    {"S", 16}, {"Cl", 17}, {"Ar", 18}, {"K", 19}, {"Ca", 20},
    {"Sc", 21}, {"Ti", 22}, {"V", 23}, {"Cr", 24}, {"Mn", 25},
    {"Fe", 26}, {"Co", 27}, {"Ni", 28}, {"Cu", 29}, {"Zn", 30},
    {"Ga", 31}, {"Ge", 32}, {"As", 33}, {"Se", 34}, {"Br", 35},
    {"Kr", 36}, {"Rb", 37}, {"Sr", 38}, {"Y", 39}, {"Zr", 40},
    {"Nb", 41}, {"Mo", 42}, {"Tc", 43}, {"Ru", 44}, {"Rh", 45},
    {"Pd", 46}, {"Ag", 47}, {"Cd", 48}, {"In", 49}, {"Sn", 50},
    {"Sb", 51}, {"Te", 52}, {"I", 53}, {"Xe", 54}, {"Cs", 55},
    {"Ba", 56}, {"La", 57}, {"Ce", 58}, {"Pr", 59}, {"Nd", 60},
    {"Pm", 61}, {"Sm", 62}, {"Eu", 63}, {"Gd", 64}, {"Tb", 65},
    {"Dy", 66}, {"Ho", 67}, {"Er", 68}, {"Tm", 69}, {"Yb", 70},
    {"Lu", 71}, {"Hf", 72}, {"Ta", 73}, {"W", 74}, {"Re", 75},
    {"Os", 76}, {"Ir", 77}, {"Pt", 78}, {"Au", 79}, {"Hg", 80},
    {"Tl", 81}, {"Pb", 82}, {"Bi", 83}, {"Po", 84}, {"At", 85},
    {"Rn", 86}, {"Fr", 87}, {"Ra", 88}, {"Ac", 89}, {"Th", 90},
    {"Pa", 91}, {"U", 92}, {"Np", 93}, {"Pu", 94}, {"Am", 95},
    {"Cm", 96}, {"Bk", 97}, {"Cf", 98}, {"Es", 99}, {"Fm", 100},
    {"Md", 101}, {"No", 102}, {"Lr", 103}, {"Rf", 104}, {"Db", 105},
    {"Sg", 106}, {"Bh", 107}, {"Hs", 108}, {"Mt", 109}, {"Ds", 110},
    {"Rg", 111}, {"Cn", 112}, {"Nh", 113}, {"Fl", 114}, {"Mc", 115},
    {"Lv", 116}, {"Ts", 117}, {"Og", 118}
};

const std::map<std::string, double> atomWeightMap = {
    {"H", 1.0080}, {"He", 4.0026}, {"Li", 6.9400}, {"Be", 9.0122},
    {"B", 10.810}, {"C", 12.011}, {"N", 14.007}, {"O", 15.999},
    {"F", 18.998}, {"Ne", 20.180}, {"Na", 22.990}, {"Mg", 24.305},
    {"Al", 26.982}, {"Si", 28.085}, {"P", 30.974}, {"S", 32.06},
    {"Cl", 35.45}, {"Ar", 39.95}, {"K", 39.098}, {"Ca", 40.078},
    {"Sc", 44.956}, {"Ti", 47.867}, {"V", 50.942}, {"Cr", 51.996},
    {"Mn", 54.938}, {"Fe", 55.845}, {"Co", 58.933}, {"Ni", 58.693},
    {"Cu", 63.546}, {"Zn", 65.38}, {"Ga", 69.723}, {"Ge", 72.630},
    {"As", 74.922}, {"Se", 78.971}, {"Br", 79.904}, {"Kr", 83.798},
    {"Rb", 85.468}, {"Sr", 87.62}, {"Y", 88.906}, {"Zr", 91.224},
    {"Nb", 92.906}, {"Mo", 95.95}, {"Ru", 101.07}, {"Rh", 102.91},
    {"Pd", 106.42}, {"Ag", 107.87}, {"Cd", 112.41}, {"In", 114.82},
    {"Sn", 118.71}, {"Sb", 121.76}, {"Te", 127.60}, {"I", 126.90},
    {"Xe", 131.29}, {"Cs", 132.91}, {"Ba", 137.33}, {"La", 138.91},
    {"Ce", 140.91}, {"Pr", 140.91}, {"Nd", 144.24}, {"Sm", 150.36},
    {"Eu", 151.96}, {"Gd", 157.25}, {"Tb", 158.93}, {"Dy", 162.50},
    {"Ho", 164.93}, {"Er", 167.26}, {"Tm", 168.93}, {"Yb", 173.05},
    {"Lu", 174.97}, {"Hf", 178.49}, {"Ta", 180.95}, {"W", 183.84},
    {"Re", 186.21}, {"Os", 190.23}, {"Ir", 192.22}, {"Pt", 195.08},
    {"Au", 196.97}, {"Hg", 200.59}, {"Tl", 204.38}, {"Pb", 207.2},
    {"Bi", 208.98}, {"Th", 232.04}, {"Pa", 231.04}, {"U", 238.03}
};

const std::map<std::string, int> oxiStateMap = {
    {"Li", 1}, {"Be", 2}, {"F", -1}, {"Na", 1}, {"Mg", 2},
    {"Al", 3}, {"Cl", -1}, {"K", 1}, {"Ca", 2}, {"Br", -1},
    {"Rb", 1}, {"Sr", 2}, {"Y", 3}, {"Zr", 4}, {"I", -1},
    {"Cs", 1}, {"Ba", 2}, {"La", 3}, {"Ce", 3}, {"Pr", 3},
    {"Pm", 3}, {"Sm", 3}, {"Eu", 3}, {"Gd", 3}, {"Tb", 3},
    {"Dy", 3}, {"Ho", 3}, {"Er", 3}, {"Tm", 3}, {"Hf", 4},
    {"Th", 4}, {"Pa", 4}, {"U", 4}, {"Np", 3}, {"Pu", 3},
    {"Am", 3}, {"Cm", 3}
};

const std::map<std::string, strVect> surrogateMap{
    {"I", strVect{"Br", "I"}},
    {"Ca", strVect{"Ca", "Ba", "Sr"}},
    {"La", strVect{"Y", "La", "Pr", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm"}},
    {"Pu", strVect{"Np", "Pu", "Am", "Cm"}},
    {"Th", strVect{"Zr", "Th", "Pa"}}
}; // List of surrogated elements

const std::map<std::string, std::array<double,8>> heatData = {
    /* The std::array stores numbers in this order:
    ΔHf: Enthalpy of formation at 298.15 K (kJ/mol)
    ΔS: Entropy at 298.15 K (J/mol K)
    The next five items represent the five coefficients of the specific
    heat capacity (J/mol K) as a function of temperature (K):
    Cp = a + bT + cT^2 + dT^3 + eT^(-2)
    T_ref: Reference temperature
    */
    {"F2", std::array<double, 8>{0, 202.789, 31.4451, 0.008413831, -2.7789E-06, 2.1810E-10, -211175, 298.15}},
    {"H2", std::array<double, 8>{0, 130.68, 33.066178, -0.011363417, 1.14328E-05, -2.7729E-09, -158558, 298.15}},
    {"HF", std::array<double, 8>{-272.55, 173.78, 30.11693, -0.003246612, 2.86812E-06, 4.57914E-10, -24861, 298.15}},
    {"UF3", std::array<double, 8>{-1500.897, 103.66, 130, 0, 0, 0, 0, 298.15}},
    {"UF4", std::array<double, 8>{-1914.66, 115.4, 174.74, 0, 0, 0, 0, 298.15}},
    {"Cr0", std::array<double, 8>{0, 23.618, 7.489737, 0.07150498, -9.1676E-05, 4.60445E-08, 138157, 298.15}},
    {"Cr1", std::array<double, 8>{0, 0, 18.46508, 0.005477986, 7.90433E-06, -1.1478E-09, 1265791, 600}},
    {"CrF2", std::array<double, 8>{-764.692, 86.308, 100, 0, 0, 0, 0, 298.15}},
    {"CrF3", std::array<double, 8>{-1125.281, 83.0567, 130, 0, 0, 0, 0, 298.15}},
    {"Fe0", std::array<double, 8>{0, 27.321, 18.42868, 0.02464301, -8.9137E-06, 9.66471E-09, 12643, 298.15}},
    {"Fe1", std::array<double, 8>{0, 0, -57767.65, 137.9197, -0.1227732, 3.86824E-05, 3.9931E+09, 700}},
    {"FeF2", std::array<double, 8>{-674.55, 92.54, 98.324, 0, 0, 0, 0, 298.15}},
    {"Ni0", std::array<double, 8>{0, 29.87, 13.6916, 0.08249509, -0.000174955, 1.6160E-07, -92417, 298.15}},
    {"Ni1", std::array<double, 8>{0, 0, 1248.045, -1.25751, 0, 0, -165126600, 600}},
    {"Ni2", std::array<double, 8>{0, 0, 16.49839, 0.01874913, -6.6398E-06, 1.71728E-09, 1872051, 700}},
    {"NiF2", std::array<double, 8>{-616.155, 74.5, 100, 0, 0, 0, 0, 298.15}},
};

const double R = 8.314;

const double gamma_UF3 = 50;
const double gamma_UF4 = 0.55;
const double gamma_Inf_CrF2 = 0.5;
const double gamma_Inf_FeF2 = 1.6;
const double gamma_Inf_NiF2 = 1.0;

double G_HF(const double T);
double G_UF4(const double T);
double G_CrF2(const double T);
double G_CrF3(const double T);
double G_FeF2(const double T);
double G_NiF2(const double T);
double G_F(const double xUF3, const double xUF4, const double T);

#endif // THERMOELECTROCHEM_H_INCLUDED
