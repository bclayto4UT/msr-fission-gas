import numpy as np

# Constants
R = 8.314  # Gas constant (J/mol·K)

# Activity coefficients
gamma_UF3 = 50
gamma_UF4 = 0.55
gamma_Inf_CrF2 = 0.5
gamma_Inf_FeF2 = 1.6

def gamma_Inf_NiF2(x_LiF):
    """
    Calculate activity coefficient for NiF2 as function of LiF mole fraction.
    
    Args:
        x_LiF (float): Mole fraction of LiF
        
    Returns:
        float: Activity coefficient
    """
    return np.exp(-17.36486818 + 68.16877934*x_LiF - 81.06301593*x_LiF**2 + 26.70117623*x_LiF**3)

# Dictionary of thermodynamic data
# Format: [ΔHf, ΔS, a, b, c, d, e, f, T_ref]
# Cp = a + bT + cT^2 + dT^3 + eT^(-1) + fT^(-2)
heat_data = {
    "H2":     [0, 130.68, 33.066178, -0.011363417, 1.14328E-05, -2.7729E-09, 0, -158558, 298.15],
    "HF":     [-272.55, 173.78, 30.11693, -0.003246612, 2.86812E-06, 4.57914E-10, 0, -24861, 298.15],
    "F2":     [0, 202.789, 31.4451, 0.008413831, -2.7789E-06, 2.1810E-10, 0, -211175, 298.15],
    "Cr0":    [0, 23.618, 7.489737, 0.07150498, -9.1676E-05, 4.60445E-08, 0, 138157, 298.15],
    "Cr1":    [0, 0, 18.46508, 0.005477986, 7.90433E-06, -1.1478E-09, 0, 1265791, 600],
    "CrF2":   [-764.692, 86.308, 100, 0, 0, 0, 0, 0, 298.15],
    "CrF3":   [-1125.281, 83.0567, 130, 0, 0, 0, 0, 0, 298.15],
    "Fe0":    [0, 27.321, 18.42868, 0.02464301, -8.9137E-06, 9.66471E-09, 0, -12643, 298.15],
    "Fe1":    [0, 0, -57767.65, 137.9197, -0.1227732, 3.86824E-05, 0, 3.9931E+9, 700],
    "FeF2":   [-674.55, 92.54, 98.324, 0, 0, 0, 0, 0, 298.15],
    "FeF30":  [-1041.82, 98.28, -29794.50, 192.68, -0.464154, 3.96689e-4, 0, 2.83084e8, 298.15],
    "FeF31":  [0, 0, 15482.7, -77.2309, 0.144403, -9.54538e-5, 0, -2.38266e8, 367],
    "FeF32":  [0, 0, 86.66, 0.0170956, -2.794011e-6, 5.82977e-10, 0, 293204, 450],
    "Co0":    [0, 30.07, 10.99430, 0.054375, -5.55132e-5, 2.5817e-8, 0, 164533, 298.15],
    "Co1":    [0, 0, -204.5760, 0.515582, -4.2155e-4, 1.295580e-7, 0, 17926700, 700],
    "CoF2":   [-637.69, 86.24, 104.6, 0, 0, 0, 0, 0, 298.15],
    "CoF3":   [-790.36, 94.52, 118.9327, -0.04470562, 4.084630e-5, -9.721107e-9, 0, -1505467, 298.15],
    "Ni0":    [0, 29.87, 13.6916, 0.08249509, -0.000174955, 1.6160E-07, 0, -92417, 298.15],
    "Ni1":    [0, 0, 1248.045, -1.25751, 0, 0, 0, -165126600, 600],
    "Ni2":    [0, 0, 16.49839, 0.01874913, -6.6398e-6, 1.71728e-9, 0, 1872051, 700],
    "NiF2":   [-616.155, 74.5, 100, 0, 0, 0, 0, 0, 298.15],
    "Nb":     [0, 36.47, 22.01430, 9.888160e-3, -5.648530e-6, 1.759691e-9, 0, 21839, 298.15],
    "NbF":    [235.76, 239.86, 37.846, 5.410e-3, -1.845e-6, 3.045e-10, 0, -137500, 298.15],
    "NbF2":   [-338.66, 284.71, 64.463, -2.330e-3, -3.321e-7, 1.732e-10, 0, -930200, 298.15],
    "NbF3":   [-804.41, 300.25, 124.879, -0.02440, 6.799e-6, -6.737e-10, -21354.8, 1.931e6, 298.15],
    "NbF4":   [-1296.46, 331.76, 143.667, -0.01918, 4.495e-6, -3.858e-10, -21676.7, 1.905e6, 298.15],
    "NbF5":   [-1720.25, 346.86, 121.254, 0.01596, -7.436e-6, 1.142e-9, 0, -1.942e6, 298.15],
    "Nb2F10": [-3531.16, 534.92, 261.303, 0.02902, -1.352e-5, 2.078e-9, 0, -4.027e6, 298.15],
    "Nb3F15": [-5334.84, 731.65, 400.55, 0.04275, -1.971e-5, 2.993e-9, 0, -5.990e6, 298.15],
    "Mo":     [0, 28.60, 24.72736, 3.960425e-3, -1.270706e-6, 1.153065e-9, 0, -170246, 298.15],
    "MoF4":   [-947.68, 328.93, 80.13824, 0.06350308, -5.372925e-5, 1.585916e-8, 0, -847072, 298.15],
    "MoF5":   [-1241.39, 347.52, 132.6939, 1.224478e-3, -2.05483e-7, 1.4022e-11, 0, -2.767901e6, 298.15],
    "MoF6":   [-1557.66, 350.47, 150.6792, 6.156463e-3, -1.657760e-6, 1.40911e-10, 0, -2.950587e6, 298.15],
    "UF3":    [-1500.897, 103.66, 130, 0, 0, 0, 0, 0, 298.15],
    "UF4":    [-1914.66, 115.4, 174.74, 0, 0, 0, 0, 0, 298.15],
}

def calc_enthalpy(data, T):
    """
    Calculate enthalpy in J/mol at temperature T.
    
    Args:
        data (list): Thermodynamic data array
        T (float): Temperature in K
        
    Returns:
        float: Enthalpy in J/mol
    """
    H = data[0] * 1000  # Convert from kJ/mol to J/mol
    T_ref = data[8]  # Reference temperature
    H += data[2] * (T - T_ref)
    H += data[3] * (T**2 - T_ref**2) / 2
    H += data[4] * (T**3 - T_ref**3) / 3
    H += data[5] * (T**4 - T_ref**4) / 4
    H += data[6] * np.log(T / T_ref)
    H -= data[7] * (1/T - 1/T_ref)
    return H

def calc_entropy(data, T):
    """
    Calculate entropy in J/(mol·K) at temperature T.
    
    Args:
        data (list): Thermodynamic data array
        T (float): Temperature in K
        
    Returns:
        float: Entropy in J/(mol·K)
    """
    S = data[1]  # Standard entropy
    T_ref = data[8]  # Reference temperature
    S += data[2] * np.log(T / T_ref)
    S += data[3] * (T - T_ref)
    S += data[4] * (T**2 - T_ref**2) / 2
    S += data[5] * (T**3 - T_ref**3) / 3
    S -= data[6] * (1/T - 1/T_ref)
    S -= data[7] * (T**(-2) - T_ref**(-2)) / 2
    return S

def calc_gibbs_energy(data, T):
    """
    Calculate Gibbs free energy in J/mol at temperature T.
    
    Args:
        data (list): Thermodynamic data array
        T (float): Temperature in K
        
    Returns:
        float: Gibbs free energy in J/mol
    """
    H = calc_enthalpy(data, T)
    S = calc_entropy(data, T)
    return H - T * S

# Free energy functions for various fluorides

def g_hf(T):
    """
    Free energy of formation of HF in J/mol.
    
    Args:
        T (float): Temperature in K
        
    Returns:
        float: Free energy of formation in J/mol
    """
    g_H2 = calc_gibbs_energy(heat_data["H2"], T)
    g_F2 = calc_gibbs_energy(heat_data["F2"], T)
    g_HF = calc_gibbs_energy(heat_data["HF"], T)
    return g_HF - 0.5 * g_H2 - 0.5 * g_F2

def g_uf4(T):
    """
    Free energy of reaction UF3 + 1/2 F2 <-> UF4 in J/mol.
    
    Args:
        T (float): Temperature in K
        
    Returns:
        float: Free energy of reaction in J/mol
    """
    g_F2 = calc_gibbs_energy(heat_data["F2"], T)
    g_UF3 = calc_gibbs_energy(heat_data["UF3"], T)
    g_UF4 = calc_gibbs_energy(heat_data["UF4"], T)
    return g_UF4 - 0.5 * g_F2 - g_UF3

def g_crf2(T):
    """
    Free energy of formation of CrF2 in J/mol.
    
    Args:
        T (float): Temperature in K
        
    Returns:
        float: Free energy of formation in J/mol
    """
    if T <= 600:
        g_Cr = calc_gibbs_energy(heat_data["Cr0"], T)
    else:
        h_Cr = calc_enthalpy(heat_data["Cr0"], 600) + calc_enthalpy(heat_data["Cr1"], T)
        s_Cr = calc_entropy(heat_data["Cr0"], 600) + calc_entropy(heat_data["Cr1"], T)
        g_Cr = h_Cr - T * s_Cr
    
    g_F2 = calc_gibbs_energy(heat_data["F2"], T)
    g_CrF2 = calc_gibbs_energy(heat_data["CrF2"], T)
    return g_CrF2 - g_F2 - g_Cr

def g_crf3(T):
    """
    Free energy of reaction CrF2 + 1/2 F2 <-> CrF3 in J/mol.
    
    Args:
        T (float): Temperature in K
        
    Returns:
        float: Free energy of reaction in J/mol
    """
    g_F2 = calc_gibbs_energy(heat_data["F2"], T)
    g_CrF2 = calc_gibbs_energy(heat_data["CrF2"], T)
    g_CrF3 = calc_gibbs_energy(heat_data["CrF3"], T)
    return g_CrF3 - 0.5 * g_F2 - g_CrF2

def g_fef2(T):
    """
    Free energy of formation of FeF2 in J/mol.
    
    Args:
        T (float): Temperature in K
        
    Returns:
        float: Free energy of formation in J/mol
    """
    if T <= 700:
        g_Fe = calc_gibbs_energy(heat_data["Fe0"], T)
    else:
        h_Fe = calc_enthalpy(heat_data["Fe0"], 700) + calc_enthalpy(heat_data["Fe1"], T)
        s_Fe = calc_entropy(heat_data["Fe0"], 700) + calc_entropy(heat_data["Fe1"], T)
        g_Fe = h_Fe - T * s_Fe
    
    g_F2 = calc_gibbs_energy(heat_data["F2"], T)
    g_FeF2 = calc_gibbs_energy(heat_data["FeF2"], T)
    return g_FeF2 - g_F2 - g_Fe

def g_fef3(T):
    """
    Free energy of reaction FeF2 + 1/2 F2 <-> FeF3 in J/mol.
    
    Args:
        T (float): Temperature in K
        
    Returns:
        float: Free energy of reaction in J/mol
    """
    g_F2 = calc_gibbs_energy(heat_data["F2"], T)
    g_FeF2 = calc_gibbs_energy(heat_data["FeF2"], T)
    
    if T <= 367:
        g_FeF3 = calc_gibbs_energy(heat_data["FeF30"], T)
    elif T <= 450:
        h_FeF3 = calc_enthalpy(heat_data["FeF30"], 367) + calc_enthalpy(heat_data["FeF31"], T)
        s_FeF3 = calc_entropy(heat_data["FeF30"], 367) + calc_entropy(heat_data["FeF31"], T)
        g_FeF3 = h_FeF3 - T * s_FeF3
    else:
        h_FeF3 = (calc_enthalpy(heat_data["FeF30"], 367) + 
                 calc_enthalpy(heat_data["FeF31"], 450) + 
                 calc_enthalpy(heat_data["FeF32"], T))
        s_FeF3 = (calc_entropy(heat_data["FeF30"], 367) + 
                 calc_entropy(heat_data["FeF31"], 450) + 
                 calc_entropy(heat_data["FeF32"], T))
        g_FeF3 = h_FeF3 - T * s_FeF3
    
    return g_FeF3 - 0.5 * g_F2 - g_FeF2

def g_cof2(T):
    """
    Free energy of formation of CoF2 in J/mol.
    
    Args:
        T (float): Temperature in K
        
    Returns:
        float: Free energy of formation in J/mol
    """
    if T <= 700:
        g_Co = calc_gibbs_energy(heat_data["Co0"], T)
    else:
        h_Co = calc_enthalpy(heat_data["Co0"], 700) + calc_enthalpy(heat_data["Co1"], T)
        s_Co = calc_entropy(heat_data["Co0"], 700) + calc_entropy(heat_data["Co1"], T)
        g_Co = h_Co - T * s_Co
    
    g_F2 = calc_gibbs_energy(heat_data["F2"], T)
    g_CoF2 = calc_gibbs_energy(heat_data["CoF2"], T)
    return g_CoF2 - g_F2 - g_Co

def g_cof3(T):
    """
    Free energy of reaction CoF2 + 1/2 F2 <-> CoF3 in J/mol.
    
    Args:
        T (float): Temperature in K
        
    Returns:
        float: Free energy of reaction in J/mol
    """
    g_F2 = calc_gibbs_energy(heat_data["F2"], T)
    g_CoF2 = calc_gibbs_energy(heat_data["CoF2"], T)
    g_CoF3 = calc_gibbs_energy(heat_data["CoF3"], T)
    return g_CoF3 - 0.5 * g_F2 - g_CoF2

def g_mnf2(T):
    """
    Free energy of formation of MnF2 in J/mol.
    
    Args:
        T (float): Temperature in K
        
    Returns:
        float: Free energy of formation in J/mol
    """
    return (-203008 + 30.96 * T) * 4.184  # Convert from cal to J

def g_nif2(T):
    """
    Free energy of formation of NiF2 in J/mol.
    
    Args:
        T (float): Temperature in K
        
    Returns:
        float: Free energy of formation in J/mol
    """
    if T <= 600:
        g_Ni = calc_gibbs_energy(heat_data["Ni0"], T)
    elif T <= 700:
        h_Ni = calc_enthalpy(heat_data["Ni0"], 600) + calc_enthalpy(heat_data["Ni1"], T)
        s_Ni = calc_entropy(heat_data["Ni0"], 600) + calc_entropy(heat_data["Ni1"], T)
        g_Ni = h_Ni - T * s_Ni
    else:
        h_Ni = (calc_enthalpy(heat_data["Ni0"], 600) + 
               calc_enthalpy(heat_data["Ni1"], 700) + 
               calc_enthalpy(heat_data["Ni2"], T))
        s_Ni = (calc_entropy(heat_data["Ni0"], 600) + 
               calc_entropy(heat_data["Ni1"], 700) + 
               calc_entropy(heat_data["Ni2"], T))
        g_Ni = h_Ni - T * s_Ni
    
    g_F2 = calc_gibbs_energy(heat_data["F2"], T)
    g_NiF2 = calc_gibbs_energy(heat_data["NiF2"], T)
    return g_NiF2 - g_F2 - g_Ni

# Niobium fluoride functions
def g_nbf5(T):
    """
    Free energy of reaction Nb + 5/2 F2 <-> NbF5 in J/mol.
    
    Args:
        T (float): Temperature in K
        
    Returns:
        float: Free energy of reaction in J/mol
    """
    g_F2 = calc_gibbs_energy(heat_data["F2"], T)
    g_Nb = calc_gibbs_energy(heat_data["Nb"], T)
    g_NbF5 = calc_gibbs_energy(heat_data["NbF5"], T)
    return g_NbF5 - 2.5 * g_F2 - g_Nb

def g_nb2f10(T):
    """
    Free energy of reaction 2NbF5 <-> Nb2F10 in J/mol.
    
    Args:
        T (float): Temperature in K
        
    Returns:
        float: Free energy of reaction in J/mol
    """
    g_NbF5 = calc_gibbs_energy(heat_data["NbF5"], T)
    g_Nb2F10 = calc_gibbs_energy(heat_data["Nb2F10"], T)
    return g_Nb2F10 - 2 * g_NbF5

def g_nb3f15(T):
    """
    Free energy of reaction NbF5 + Nb2F10 <-> Nb3F15 in J/mol.
    
    Args:
        T (float): Temperature in K
        
    Returns:
        float: Free energy of reaction in J/mol
    """
    g_NbF5 = calc_gibbs_energy(heat_data["NbF5"], T)
    g_Nb2F10 = calc_gibbs_energy(heat_data["Nb2F10"], T)
    g_Nb3F15 = calc_gibbs_energy(heat_data["Nb3F15"], T)
    return g_Nb3F15 - g_Nb2F10 - g_NbF5

def g_nbf4(T):
    """
    Free energy of reaction NbF4 + 1/2 F2 <-> NbF5 in J/mol.
    
    Args:
        T (float): Temperature in K
        
    Returns:
        float: Free energy of reaction in J/mol
    """
    g_F2 = calc_gibbs_energy(heat_data["F2"], T)
    g_NbF4 = calc_gibbs_energy(heat_data["NbF4"], T)
    g_NbF5 = calc_gibbs_energy(heat_data["NbF5"], T)
    return g_NbF5 - 0.5 * g_F2 - g_NbF4

def g_nbf3(T):
    """
    Free energy of reaction NbF3 + 1/2 F2 <-> NbF4 in J/mol.
    
    Args:
        T (float): Temperature in K
        
    Returns:
        float: Free energy of reaction in J/mol
    """
    g_F2 = calc_gibbs_energy(heat_data["F2"], T)
    g_NbF3 = calc_gibbs_energy(heat_data["NbF3"], T)
    g_NbF4 = calc_gibbs_energy(heat_data["NbF4"], T)
    return g_NbF4 - 0.5 * g_F2 - g_NbF3

def g_nbf2(T):
    """
    Free energy of reaction NbF2 + 1/2 F2 <-> NbF3 in J/mol.
    
    Args:
        T (float): Temperature in K
        
    Returns:
        float: Free energy of reaction in J/mol
    """
    g_F2 = calc_gibbs_energy(heat_data["F2"], T)
    g_NbF2 = calc_gibbs_energy(heat_data["NbF2"], T)
    g_NbF3 = calc_gibbs_energy(heat_data["NbF3"], T)
    return g_NbF3 - 0.5 * g_F2 - g_NbF2

def g_nbf(T):
    """
    Free energy of reaction NbF + 1/2 F2 <-> NbF2 in J/mol.
    
    Args:
        T (float): Temperature in K
        
    Returns:
        float: Free energy of reaction in J/mol
    """
    g_F2 = calc_gibbs_energy(heat_data["F2"], T)
    g_NbF = calc_gibbs_energy(heat_data["NbF"], T)
    g_NbF2 = calc_gibbs_energy(heat_data["NbF2"], T)
    return g_NbF2 - 0.5 * g_F2 - g_NbF

# Molybdenum fluoride functions
def g_mof4(T):
    """
    Free energy of formation of MoF4 in J/mol.
    
    Args:
        T (float): Temperature in K
        
    Returns:
        float: Free energy of formation in J/mol
    """
    g_Mo = calc_gibbs_energy(heat_data["Mo"], T)
    g_F2 = calc_gibbs_energy(heat_data["F2"], T)
    g_MoF4 = calc_gibbs_energy(heat_data["MoF4"], T)
    return g_MoF4 - 2 * g_F2 - g_Mo

def g_mof5(T):
    """
    Free energy of reaction MoF4 + 1/2 F2 <-> MoF5 in J/mol.
    
    Args:
        T (float): Temperature in K
        
    Returns:
        float: Free energy of reaction in J/mol
    """
    g_MoF4 = calc_gibbs_energy(heat_data["MoF4"], T)
    g_F2 = calc_gibbs_energy(heat_data["F2"], T)
    g_MoF5 = calc_gibbs_energy(heat_data["MoF5"], T)
    return g_MoF5 - 0.5 * g_F2 - g_MoF4

def g_mof6(T):
    """
    Free energy of reaction MoF5 + 1/2 F2 <-> MoF6 in J/mol.
    
    Args:
        T (float): Temperature in K
        
    Returns:
        float: Free energy of reaction in J/mol
    """
    g_MoF5 = calc_gibbs_energy(heat_data["MoF5"], T)
    g_F2 = calc_gibbs_energy(heat_data["F2"], T)
    g_MoF6 = calc_gibbs_energy(heat_data["MoF6"], T)
    return g_MoF6 - 0.5 * g_F2 - g_MoF5

def g_f(x_uf3, x_uf4, T):
    """
    Fluorine potential as defined by Olander (2001) in J/mol.
    
    Args:
        x_uf3 (float): Mole fraction of UF3
        x_uf4 (float): Mole fraction of UF4
        T (float): Temperature in K
        
    Returns:
        float: Fluorine potential in J/mol
    """
    return 2 * R * T * np.log(gamma_UF4 * x_uf4 / (gamma_UF3 * x_uf3)) + 2 * g_uf4(T)
