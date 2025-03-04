"""
Core Data Structures and Utility Module for Thermochimica Data Processing.

This module provides utility functions and data structures for element properties,
surrogate mappings, and string processing capabilities required for molten salt
chemistry calculations.
"""

import pandas as pd
import numpy as np
from typing import List, Dict, Union, Optional, Any

# Constants
e = np.e
phi = (1 + np.sqrt(5)) / 2
pi = np.pi

# Element mappings from atomic number
atom_num_map = {
    "H": 1, "He": 2, "Li": 3, "Be": 4, "B": 5,
    "C": 6, "N": 7, "O": 8, "F": 9, "Ne": 10,
    "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15,
    "S": 16, "Cl": 17, "Ar": 18, "K": 19, "Ca": 20,
    "Sc": 21, "Ti": 22, "V": 23, "Cr": 24, "Mn": 25,
    "Fe": 26, "Co": 27, "Ni": 28, "Cu": 29, "Zn": 30,
    "Ga": 31, "Ge": 32, "As": 33, "Se": 34, "Br": 35,
    "Kr": 36, "Rb": 37, "Sr": 38, "Y": 39, "Zr": 40,
    "Nb": 41, "Mo": 42, "Tc": 43, "Ru": 44, "Rh": 45,
    "Pd": 46, "Ag": 47, "Cd": 48, "In": 49, "Sn": 50,
    "Sb": 51, "Te": 52, "I": 53, "Xe": 54, "Cs": 55,
    "Ba": 56, "La": 57, "Ce": 58, "Pr": 59, "Nd": 60,
    "Pm": 61, "Sm": 62, "Eu": 63, "Gd": 64, "Tb": 65,
    "Dy": 66, "Ho": 67, "Er": 68, "Tm": 69, "Yb": 70,
    "Lu": 71, "Hf": 72, "Ta": 73, "W": 74, "Re": 75,
    "Os": 76, "Ir": 77, "Pt": 78, "Au": 79, "Hg": 80,
    "Tl": 81, "Pb": 82, "Bi": 83, "Po": 84, "At": 85,
    "Rn": 86, "Fr": 87, "Ra": 88, "Ac": 89, "Th": 90,
    "Pa": 91, "U": 92, "Np": 93, "Pu": 94, "Am": 95,
    "Cm": 96, "Bk": 97, "Cf": 98, "Es": 99, "Fm": 100,
    "Md": 101, "No": 102, "Lr": 103, "Rf": 104, "Db": 105,
    "Sg": 106, "Bh": 107, "Hs": 108, "Mt": 109, "Ds": 110,
    "Rg": 111, "Cn": 112, "Nh": 113, "Fl": 114, "Mc": 115,
    "Lv": 116, "Ts": 117, "Og": 118
}

# Element atomic weights
atom_weight_map = {
    "H": 1.0080, "He": 4.0026, "Li": 6.9400, "Be": 9.0122,
    "B": 10.810, "C": 12.011, "N": 14.007, "O": 15.999,
    "F": 18.998, "Ne": 20.180, "Na": 22.990, "Mg": 24.305,
    "Al": 26.982, "Si": 28.085, "P": 30.974, "S": 32.06,
    "Cl": 35.45, "Ar": 39.95, "K": 39.098, "Ca": 40.078,
    "Sc": 44.956, "Ti": 47.867, "V": 50.942, "Cr": 51.996,
    "Mn": 54.938, "Fe": 55.845, "Co": 58.933, "Ni": 58.693,
    "Cu": 63.546, "Zn": 65.38, "Ga": 69.723, "Ge": 72.630,
    "As": 74.922, "Se": 78.971, "Br": 79.904, "Kr": 83.798,
    "Rb": 85.468, "Sr": 87.62, "Y": 88.906, "Zr": 91.224,
    "Nb": 92.906, "Mo": 95.95, "Ru": 101.07, "Rh": 102.91,
    "Pd": 106.42, "Ag": 107.87, "Cd": 112.41, "In": 114.82,
    "Sn": 118.71, "Sb": 121.76, "Te": 127.60, "I": 126.90,
    "Xe": 131.29, "Cs": 132.91, "Ba": 137.33, "La": 138.91,
    "Ce": 140.91, "Pr": 140.91, "Nd": 144.24, "Sm": 150.36,
    "Eu": 151.96, "Gd": 157.25, "Tb": 158.93, "Dy": 162.50,
    "Ho": 164.93, "Er": 167.26, "Tm": 168.93, "Yb": 173.05,
    "Lu": 174.97, "Hf": 178.49, "Ta": 180.95, "W": 183.84,
    "Re": 186.21, "Os": 190.23, "Ir": 192.22, "Pt": 195.08,
    "Au": 196.97, "Hg": 200.59, "Tl": 204.38, "Pb": 207.2,
    "Bi": 208.98, "Th": 232.04, "Pa": 231.04, "U": 238.03
}

# Element oxidation states
oxi_state_map = {
    "Li": 1, "Be": 2, "F": -1, "Na": 1, "Mg": 2,
    "Al": 3, "Cl": -1, "K": 1, "Ca": 2, "Br": -1,
    "Rb": 1, "Sr": 2, "Y": 3, "Zr": 4, "I": -1,
    "Cs": 1, "Ba": 2, "La": 3, "Ce": 3, "Pr": 3,
    "Pm": 3, "Sm": 3, "Eu": 3, "Gd": 3, "Tb": 3,
    "Dy": 3, "Ho": 3, "Er": 3, "Tm": 3, "Hf": 4,
    "Th": 4, "Pa": 4, "U": 4, "Np": 3, "Pu": 3,
    "Am": 3, "Cm": 3
}

# Surrogate element groupings
surrogate_map = {
    "I": ["Br", "I"],
    "Ca": ["Ca", "Ba", "Sr"],
    "La": ["Y", "La", "Pr", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm"],
    "Pu": ["Np", "Pu", "Am", "Cm"],
    "Th": ["Zr", "Th", "Pa"]
}

# Individual elements → surrogate representative
surrogate_map_inv = {
    "Br": "I", "Ba": "Ca", "Sr": "Ca", "Y": "La", "Pr": "La",
    "Pm": "La", "Sm": "La", "Eu": "La", "Gd": "La", "Tb": "La",
    "Dy": "La", "Ho": "La", "Er": "La", "Tm": "La", "Np": "Pu",
    "Am": "Pu", "Cm": "Pu", "Zr": "Th", "Pa": "Th"
}

# Thermodynamic data constants
CP_TERMS = 6
# Heat capacity data will be implemented in the thermochemistry module

def element_symbol(element_str: str) -> str:
    """
    Convert element name to proper case (first letter capitalized, rest lowercase).
    
    Args:
        element_str: Element name as a string
        
    Returns:
        Properly formatted element symbol
    """
    if not element_str:
        return ""
    return element_str[0].upper() + element_str[1:].lower()

def str_to_list(data: str, delimiter: str = ' ') -> List[str]:
    """
    Split a string into a list using the specified delimiter.
    Handles multiple delimiters and whitespace.
    
    Args:
        data: Input string to split
        delimiter: Character to use as delimiter (default space)
        
    Returns:
        List of strings after splitting
    """
    # Use Python's built-in split method with more robust handling
    if not data:
        return []
    
    # First handle tabs and specified delimiter
    parts = []
    current_data = data.strip()
    
    # Replace tabs with the delimiter for consistent processing
    if delimiter != '\t':
        current_data = current_data.replace('\t', delimiter)
    
    # Split and filter out empty strings
    return [part for part in current_data.split(delimiter) if part]

def is_numeric(string: str) -> bool:
    """
    Check if a string can be converted to a number.
    
    Args:
        string: String to check
        
    Returns:
        True if the string can be converted to a float, False otherwise
    """
    if not string:
        return False
    
    try:
        float(string)
        return True
    except ValueError:
        return False

def contains_number(string: str) -> bool:
    """
    Check if a string contains any numeric digits.
    
    Args:
        string: String to check
        
    Returns:
        True if the string contains any digits, False otherwise
    """
    return any(char.isdigit() for char in string)

def str_to_float_list(string: str) -> List[float]:
    """
    Convert a string to a list of floating-point numbers.
    Handles space, comma, or colon-separated values.
    For colon-separated values:
    - "start:stop" produces a linearly spaced array with steps calculated based on context
    - "start:stop:step" produces a linearly spaced array with the specified step
    
    Args:
        string: Input string containing numbers
        
    Returns:
        List of floating-point numbers
    """
    if not string:
        return []
    
    # First try space-separated
    try:
        parts = str_to_list(string)
        if len(parts) == 1 and is_numeric(parts[0]):
            return [float(parts[0])]
        return [float(part) for part in parts]
    except ValueError:
        pass
    
    # Try colon-separated for numpy-like ranges
    if ':' in string:
        parts = string.split(':')
        if len(parts) == 2 or len(parts) == 3:
            try:
                start = float(parts[0])
                stop = float(parts[1])
                
                if start == stop:
                    return [start]
                
                if len(parts) == 3:
                    step = float(parts[2])
                    if start > stop and step > 0:
                        step = -step
                    
                    # Create range from start to stop with step
                    count = int(np.ceil((stop - start) / step))
                    return [start + step * i for i in range(count)]
                else:
                    # Default to 10 steps if not specified
                    # This can be adjusted based on context if needed
                    return np.linspace(start, stop, num=10).tolist()
            except ValueError:
                pass
    
    # Try comma-separated
    try:
        parts = string.split(',')
        return [float(part.strip()) for part in parts if is_numeric(part.strip())]
    except ValueError:
        raise ValueError(f"Could not convert string '{string}' to list of floats")

def sgn(x: Union[int, float]) -> int:
    """
    Return the sign of a number: 1 for positive, -1 for negative, 0 for zero.
    
    Args:
        x: Input number
        
    Returns:
        Sign of the number as an integer
    """
    return (x > 0) - (x < 0)

def factorial(n: int) -> int:
    """
    Calculate the factorial of n.
    
    Args:
        n: Non-negative integer
        
    Returns:
        n! (n factorial)
    """
    if n < 0:
        raise ValueError("Factorial not defined for negative numbers")
    if n == 0 or n == 1:
        return 1
    return n * factorial(n - 1)

def ord_mag(x: float) -> int:
    """
    Calculate the order of magnitude of a number.
    
    Args:
        x: Input number
        
    Returns:
        Order of magnitude as an integer
    """
    if x == 0:
        return 0
    return int(np.floor(np.log10(abs(x))))

def convertible_num(string: str) -> str:
    """
    Trim leading whitespace in a string until it can be converted to a float.
    
    Args:
        string: Input string
        
    Returns:
        Modified string that can be converted to a float
    """
    if not string:
        return ""
    
    s = string
    while s:
        try:
            float(s)
            return s
        except ValueError:
            s = s[1:]
    
    return ""
