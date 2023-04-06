#include "thermoElectroChem.h"


double gamma_Inf_NiF2(double xLiF){
    return exp(-17.36486818 + 68.16877934*xLiF - 81.06301593*pow(xLiF,2) + 26.70117623*pow(xLiF,3));
}

static double calc_H(const thermoArray& data, const double T){
    // Enthalpy in J/mol
    double H = data[0]*1000;
    double T_ref = data.back();
    H += data[2]*(T-T_ref);
    H += data[3]*(pow(T,2)-pow(T_ref,2))/2;
    H += data[4]*(pow(T,3)-pow(T_ref,3))/3;
    H += data[5]*(pow(T,4)-pow(T_ref,4))/4;
    H += data[6]*log(T/T_ref);
    H -= data[7]*(1/T-1/T_ref);
    return H;
}

static double calc_S(const thermoArray& data, const double T){
    // Entropy in J/mol K
    double S = data[1];
    double T_ref = data.back();
    S += data[2]*log(T/T_ref);
    S += data[3]*(T-T_ref);
    S += data[4]*(pow(T,2)-pow(T_ref,2))/2;
    S += data[5]*(pow(T,3)-pow(T_ref,3))/3;
    S -= data[6]*(1/T-1/T_ref);
    S -= data[7]*(pow(T,-2)-pow(T_ref,-2))/2;
    return S;
}

static double calc_G(const thermoArray& data, const double T){
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

double G_CoF2(const double T){
    // Free energy of formation of CoF2 in J/mol
    double g_Co;
    if (T <= 700) g_Co = calc_G(heatData.at("Co0"), T);
    else{
        double h_Co = calc_H(heatData.at("Co0"), 700) + calc_H(heatData.at("Co1"), T);
        double s_Co = calc_S(heatData.at("Co0"), 700) + calc_S(heatData.at("Co1"), T);
        g_Co = h_Co-T*s_Co;
    }
    double g_F2 = calc_G(heatData.at("F2"), T);
    double g_CoF2 = calc_G(heatData.at("CoF2"), T);
    return g_CoF2 - g_F2 - g_Co;
}

double G_CoF3(const double T){
    // Free energy of reaction CoF2 + 1/2 F2 <-> CoF3 in J/mol
    double g_F2 = calc_G(heatData.at("F2"), T);
    double g_CoF2 = calc_G(heatData.at("CoF2"), T);
    double g_CoF3 = calc_G(heatData.at("CoF3"), T);
    return g_CoF3 - 0.5*g_F2 - g_CoF2;
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

/* Note: For niobium, the most fluorinated species (NbF5) is calculated
         directly from the metal and the lower oxidation states will
         descend from this +5 oxidation state. This is because NbF5 is
         the most chemically stable fluoride, so formulating formulae
         formula from NbF5 down will avoid loss of precision that would
         otherwise arise from subtracting near equal numbers. */

double G_NbF5(const double T){
    // Free energy of reaction Nb + 5/2 F2 <-> NbF5 in J/mol
    double g_F2 = calc_G(heatData.at("F2"), T);
    double g_Nb = calc_G(heatData.at("Nb"), T);
    double g_NbF5 = calc_G(heatData.at("NbF5"), T);
    return g_NbF5 - 2.5*g_F2 - g_Nb;
}

double G_Nb2F10(const double T){
    // Free energy of reaction 2NbF5 <-> Nb2F10 in J/mol
    double g_NbF5 = calc_G(heatData.at("NbF5"), T);
    double g_Nb2F10 = calc_G(heatData.at("Nb2F10"), T);
    return g_Nb2F10 - 2*g_NbF5;
}

double G_Nb3F15(const double T){
    // Free energy of reaction 2NbF5 <-> Nb2F10 in J/mol
    double g_NbF5 = calc_G(heatData.at("NbF5"), T);
    double g_Nb2F10 = calc_G(heatData.at("Nb2F10"), T);
    double g_Nb3F15 = calc_G(heatData.at("Nb3F15"), T);
    return g_Nb3F15 - g_Nb2F10 - g_NbF5;
}

double G_NbF4(const double T){
    // Free energy of reaction NbF4 + 1/2 F2 <-> NbF5 in J/mol
    double g_F2 = calc_G(heatData.at("F2"), T);
    double g_NbF4 = calc_G(heatData.at("NbF4"), T);
    double g_NbF5 = calc_G(heatData.at("NbF5"), T);
    return g_NbF5 - 0.5*g_F2 - g_NbF4;
}

double G_NbF3(const double T){
    // Free energy of reaction NbF3 + 1/2 F2 <-> NbF4 in J/mol
    double g_F2 = calc_G(heatData.at("F2"), T);
    double g_NbF3 = calc_G(heatData.at("NbF3"), T);
    double g_NbF4 = calc_G(heatData.at("NbF4"), T);
    return g_NbF4 - 0.5*g_F2 - g_NbF3;
}

double G_NbF2(const double T){
    // Free energy of reaction NbF2 + 1/2 F2 <-> NbF3 in J/mol
    double g_F2 = calc_G(heatData.at("F2"), T);
    double g_NbF2 = calc_G(heatData.at("NbF2"), T);
    double g_NbF3 = calc_G(heatData.at("NbF3"), T);
    return g_NbF3 - 0.5*g_F2 - g_NbF2;
}

double G_NbF(const double T){
    // Free energy of reaction NbF + 1/2 F2 <-> NbF2 in J/mol
    double g_F2 = calc_G(heatData.at("F2"), T);
    double g_NbF = calc_G(heatData.at("NbF"), T);
    double g_NbF2 = calc_G(heatData.at("NbF2"), T);
    return g_NbF2 - 0.5*g_F2 - g_NbF;
}

double G_MoF4(const double T){
    // Free energy of formation of MoF4 in J/mol
    double g_Mo = calc_G(heatData.at("Mo"), T);
    double g_F2 = calc_G(heatData.at("F2"), T);
    double g_MoF4 = calc_G(heatData.at("MoF4"), T);
    return g_MoF4 - 2*g_F2 - g_Mo;
}

double G_MoF5(const double T){
    // Free energy of reaction MoF4 + 1/2 F2 <-> MoF5 in J/mol
    double g_MoF4 = calc_G(heatData.at("MoF4"), T);
    double g_F2 = calc_G(heatData.at("F2"), T);
    double g_MoF5 = calc_G(heatData.at("MoF5"), T);
    return g_MoF5 - 0.5*g_F2 - g_MoF4;
}

double G_MoF6(const double T){
    // Free energy of reaction MoF5 + 1/2 F2 <-> MoF6 in J/mol
    double g_MoF5 = calc_G(heatData.at("MoF5"), T);
    double g_F2 = calc_G(heatData.at("F2"), T);
    double g_MoF6 = calc_G(heatData.at("MoF6"), T);
    return g_MoF6 - 0.5*g_F2 - g_MoF5;
}

double G_F(const double xUF3, const double xUF4, const double T){
    // Fluorine potential as defined by Olander (2001) in J/mol
    return 2*R*T*log(gamma_UF4*xUF4/gamma_UF3/xUF3) + 2*G_UF4(T);
}
