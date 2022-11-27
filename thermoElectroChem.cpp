#include "thermoElectroChem.h"

#include <cmath>

static const std::array<double,8> nullArray{};

static double calc_H(const std::array<double,8>& data, const double T){
    // Enthalpy in J/mol
    double H = data[0]*1000;
    H += data[2]*(T-data[7]);
    H += data[3]*(pow(T,2)-pow(data[7],2))/2;
    H += data[4]*(pow(T,3)-pow(data[7],3))/3;
    H += data[5]*(pow(T,4)-pow(data[7],4))/4;
    H -= data[6]*(1/T-1/data[7]);
    return H;
}

static double calc_S(const std::array<double,8>& data, const double T){
    // Entropy in J/mol K
    double S = data[1];
    S += data[2]*log(T/data[7]);
    S += data[3]*(T-data[7]);
    S += data[4]*(pow(T,2)-pow(data[7],2))/2;
    S += data[5]*(pow(T,3)-pow(data[7],3))/3;
    S -= data[6]*(pow(T,-2)-pow(data[7],-2))/2;
    return S;
}

static double calc_G(const std::array<double,8>& data, const double T){
    // Gibbs free energy in J/mol
    double H = calc_H(data, T);
    double S = calc_S(data, T);
    return H-T*S;
}

double G_HF(const double T){
    // Free energy of formation of HF in J/mol
    double g_H2 = calc_G(heatData.at("H2"), T);
    double g_F2 = calc_G(heatData.at("F2"), T);
    double g_HF = calc_G(heatData.at("HF"), T);
    return g_HF - 0.5*g_H2 - 0.5*g_F2;
}

double G_UF4(const double T){
    // Free energy of reaction UF3 + 1/2 F2 <-> UF4 in J/mol
    double g_F2 = calc_G(heatData.at("F2"), T);
    double g_UF3 = calc_G(heatData.at("UF3"), T);
    double g_UF4 = calc_G(heatData.at("UF4"), T);
    return g_UF4 - 0.5*g_F2 - g_UF3;
}

double G_CrF2(const double T){
    // Free energy of formation of CrF2 in J/mol
    double g_Cr;
    if (T <= 600) g_Cr = calc_G(heatData.at("Cr0"), T);
    else{
        double h_Cr = calc_H(heatData.at("Cr0"), 600) + calc_H(heatData.at("Cr1"), T);
        double s_Cr = calc_S(heatData.at("Cr0"), 600) + calc_S(heatData.at("Cr1"), T);
        g_Cr = h_Cr-T*s_Cr;
    }
    double g_F2 = calc_G(heatData.at("F2"), T);
    double g_CrF2 = calc_G(heatData.at("CrF2"), T);
    return g_CrF2 - g_F2 - g_Cr;
}

double G_CrF3(const double T){
    // Free energy of reaction CrF2 + 1/2 F2 <-> CrF3 in J/mol
    double g_F2 = calc_G(heatData.at("F2"), T);
    double g_CrF2 = calc_G(heatData.at("CrF2"), T);
    double g_CrF3 = calc_G(heatData.at("CrF3"), T);
    return g_CrF3 - 0.5*g_F2 - g_CrF2;
}

double G_FeF2(const double T){
    // Free energy of formation of FeF2 in J/mol
    double g_Fe;
    if (T <= 700) g_Fe = calc_G(heatData.at("Fe0"), T);
    else{
        double h_Fe = calc_H(heatData.at("Fe0"), 700) + calc_H(heatData.at("Fe1"), T);
        double s_Fe = calc_S(heatData.at("Fe0"), 700) + calc_S(heatData.at("Fe1"), T);
        g_Fe = h_Fe-T*s_Fe;
    }
    double g_F2 = calc_G(heatData.at("F2"), T);
    double g_FeF2 = calc_G(heatData.at("FeF2"), T);
    return g_FeF2 - g_F2 - g_Fe;
}

double G_FeF3(const double T){
    // Free energy of reaction FeF2 + 1/2 F2 <-> FeF3 in J/mol
    double g_F2 = calc_G(heatData.at("F2"), T);
    double g_FeF2 = calc_G(heatData.at("FeF2"), T);
    double g_FeF3;

    if (T <= 367) g_FeF3 = calc_G(heatData.at("FeF30"), T);
    else if (T <= 450){
        double h_FeF3 = calc_H(heatData.at("FeF30"), 367) + calc_H(heatData.at("FeF31"), T);
        double s_FeF3 = calc_S(heatData.at("FeF30"), 367) + calc_S(heatData.at("FeF31"), T);
        g_FeF3 = h_FeF3-T*s_FeF3;
    } else{
        double h_FeF3 = calc_H(heatData.at("FeF30"), 367) + calc_H(heatData.at("FeF31"), 450)
                        + calc_H(heatData.at("FeF32"), T);
        double s_FeF3 = calc_S(heatData.at("FeF30"), 367) + calc_S(heatData.at("FeF31"), 450)
                        + calc_S(heatData.at("FeF32"), T);
        g_FeF3 = h_FeF3-T*s_FeF3;
    }

    return g_FeF3 - 0.5*g_F2 - g_FeF2;
}


double G_NiF2(const double T){
    // Free energy of formation of NiF2 in J/mol
    double g_Ni;
    if (T <= 600) g_Ni = calc_G(heatData.at("Ni0"), T);
    else if (T <= 700){
        double h_Ni = calc_H(heatData.at("Ni0"), 600) + calc_H(heatData.at("Ni1"), T);
        double s_Ni = calc_S(heatData.at("Ni0"), 600) + calc_S(heatData.at("Ni1"), T);
        g_Ni = h_Ni-T*s_Ni;
    } else{
        double h_Ni = calc_H(heatData.at("Ni0"), 600) + calc_H(heatData.at("Ni1"), 700)
                    + calc_H(heatData.at("Ni2"), T);
        double s_Ni = calc_S(heatData.at("Ni0"), 600) + calc_S(heatData.at("Ni1"), 700)
                    + calc_S(heatData.at("Ni2"), T);
        g_Ni = h_Ni-T*s_Ni;
    }
    double g_F2 = calc_G(heatData.at("F2"), T);
    double g_NiF2 = calc_G(heatData.at("NiF2"), T);
    return g_NiF2 - g_F2 - g_Ni;
}

double G_F(const double xUF3, const double xUF4, const double T){
    // Fluorine potential as defined by Olander (2001) in J/mol
    return 2*R*T*log(gamma_UF4*xUF4/gamma_UF3/xUF3) + 2*G_UF4(T);
}
