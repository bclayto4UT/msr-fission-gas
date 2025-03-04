import os
import pandas as pd
import numpy as np
from typing import Dict, List, Union

# Utility functions to be imported from other modules
from utils import str_to_list, str_to_float_list
from thermochemistry import surrogate_map, surrogate_map_inv

class DataProcessor:
    """
    Core data processing class for handling Thermochimica and SCALE data
    """
    @staticmethod
    def scale_to_vector(in_file: str) -> Dict[str, List[float]]:
        """
        Parse SCALE output file to extract element concentration data.
        
        Args:
            in_file (str): Path to the SCALE .out file
        
        Returns:
            Dict[str, List[float]]: Dictionary mapping element names to concentration vectors
        """
        try:
            # Read the file, looking for the element concentration table
            with open(in_file, 'r') as f:
                # Find the start of the concentration table
                while True:
                    line = f.readline()
                    if 'relative cutoff;' in line:
                        # Skip two lines to get to the data
                        f.readline()
                        f.readline()
                        break
                
                # Initialize dictionary to store concentrations
                concentrations = {}
                
                # Parse concentration data
                while True:
                    line = f.readline().strip()
                    
                    # Stop if end of table is reached
                    if '------' in line or not line:
                        break
                    
                    # Split line into element and concentration values
                    parts = line.split()
                    
                    # Ensure valid line
                    if len(parts) < 2:
                        continue
                    
                    # Capitalize first letter of element
                    element = parts[0].capitalize()
                    
                    # Convert concentrations to floats
                    try:
                        conc_values = [float(val) for val in parts[1:]]
                        concentrations[element] = conc_values
                    except ValueError:
                        # Skip lines with non-numeric data
                        continue
            
            return concentrations
        
        except FileNotFoundError:
            raise FileNotFoundError(f"Cannot open input file: {in_file}")
        except Exception as e:
            raise RuntimeError(f"Error processing SCALE output file: {e}")

    @staticmethod
    def vector_to_therm(data_map: Dict[str, List[float]], 
                         out_file: str, 
                         str_t: str, 
                         str_p: str, 
                         includes_surr: bool = True) -> None:
        """
        Convert concentration data to Thermochimica input files (.F90)
        
        Args:
            data_map (Dict[str, List[float]]): Dictionary with element concentration data
            out_file (str): Base name for output files
            str_t (str): Temperature specification string
            str_p (str): Pressure specification string
            includes_surr (bool): Whether to include surrogate elements
        """
        # Extract file stem and handle file path
        file_parts = out_file.rsplit('.', 1)
        file_stem = file_parts[0]
        if '/' in file_stem:
            file_stem = file_stem.rsplit('/', 1)[1]

        # Determine number of time intervals
        m = len(list(data_map.values())[0])

        # Parse temperature values
        try:
            T = DataProcessor._parse_range_values(str_t, m)
        except ValueError as e:
            raise ValueError(f"Invalid temperature specification: {e}")

        # Parse pressure values
        try:
            P = DataProcessor._parse_range_values(str_p, m)
        except ValueError as e:
            raise ValueError(f"Invalid pressure specification: {e}")

        # Handle surrogate elements if needed
        if includes_surr:
            data_map = DataProcessor._handle_surrogate_elements(data_map)

        # Generate Thermochimica input files
        for i in range(m):
            filename = f"{file_stem}{i}.F90"
            with open(filename, 'w') as output:
                # Write program header
                output.write(f"program {file_stem}{i}\n")
                output.write("USE ModuleThermoIO\n")
                output.write("implicit none\n")
                output.write("cInputUnitTemperature = 'K'\n")
                output.write("cInputUnitPressure = 'atm'\n")
                output.write("cInputUnitMass = 'moles'\n")
                output.write("cThermoFileName = DATA_DIRECTORY // 'MSTDB-TC_V2.0_Fluorides_8-0.dat'\n")

                # Write pressure and temperature
                output.write(f"dPressure = {P[0] if len(P) == 1 else P[i]}\n")
                output.write(f"dTemperature = {T[0] if len(T) == 1 else T[i]}\n")

                # Write element masses
                for element, concentrations in data_map.items():
                    if concentrations[i] > 0.0:
                        # Note: In a full implementation, you'd need a way to map element 
                        # names to atomic numbers (atomNumMap in the C++ version)
                        output.write(f"dElementMass(ATOMIC_NUMBER_{element}) = {concentrations[i]}\n")

                # Write Thermochimica execution instructions
                output.write("iPrintResultsMode = 1\n")
                output.write("call ParseCSDataFile(cThermoFileName)\n")
                output.write("if (INFOThermo == 0) call Thermochimica\n")
                output.write("if (iPrintResultsMode > 0) call PrintResults\n")
                output.write("if (INFOThermo == 0) call ResetThermoAll\n")
                output.write("call ThermoDebug\n")
                output.write(f"end program {file_stem}{i}\n")

    @staticmethod
    def _parse_range_values(range_str: str, m: int) -> List[float]:
        """
        Parse temperature or pressure range specification
        
        Args:
            range_str (str): Range specification string
            m (int): Number of time intervals
        
        Returns:
            List[float]: Parsed values
        """
        # Try direct list of values first
        try:
            return [float(x) for x in re.split(r'[ ,]', range_str)]
        except ValueError:
            pass

        # Try range specification with colon
        parts = range_str.split(':')
        if len(parts) == 2 or len(parts) == 3:
            start = float(parts[0])
            stop = float(parts[1])

            if len(parts) == 2:
                # Auto-generate steps based on time intervals
                step = (stop - start) / m
            else:
                step = float(parts[2])

            # Generate values
            values = list(np.arange(start, stop, step))
            
            # Ensure the first list of values is created
            if not values:
                values = [start]
            
            return values

        raise ValueError(f"Invalid range specification: {range_str}")

    @staticmethod
    def _handle_surrogate_elements(data_map: Dict[str, List[float]]) -> Dict[str, List[float]]:
        """
        Handle surrogate elements by redistributing their concentrations
        
        Args:
            data_map (Dict[str, List[float]]): Original concentration map
        
        Returns:
            Dict[str, List[float]]: Updated concentration map
        """
        # Create a copy of the original map to modify
        updated_map = data_map.copy()

        for element in list(updated_map.keys()):
            # Check if this element has a surrogate
            try:
                surrogate = surrogate_map_inv[element]
                
                # Ensure surrogate exists in map, if not initialize
                if surrogate not in updated_map:
                    updated_map[surrogate] = [0.0] * len(updated_map[element])
                
                # Add this element's concentration to the surrogate
                for i in range(len(updated_map[element])):
                    updated_map[surrogate][i] += updated_map[element][i]
                
                # Remove the original element
                del updated_map[element]
            
            except KeyError:
                # No surrogate for this element, skip
                continue

        return updated_map

    @staticmethod
    def text_to_excel(in_file: str, out_file: str, data_type: str) -> None:
        """
        Extract specified data from Thermochimica output file.
        
        Args:
            in_file (str): Path to Thermochimica output file
            out_file (str): Path for the output CSV file
            data_type (str): Space-separated string of data types to extract
        """
        # Parse requested data types
        requested_types = data_type.split()
        
        # Special flags
        has_x_all = 'x_all' in requested_types
        has_y_all = 'y_all' in requested_types
        
        # Initialize data collection dictionaries
        vx: Dict[str, float] = {}
        vy: Dict[str, float] = {}
        extra_data: Dict[str, float] = {
            'NUF34': -1.0,
            'T': -1.0,
            'P': -1.0
        }
        
        # Prepare data collection for specific types
        for s in requested_types:
            if s in ['ni', 'nx']:
                vx[s] = 0
            elif s == 'ny':
                vy[s] = 0
            elif s.startswith('x_'):
                if s == 'x_U4+':
                    vx['U2'] = 0
                    vx['U[VI]'] = 0
                    vx['U[VII]'] = 0
                elif s == 'x_UF4' or has_x_all:
                    vx['U2F8'] = 0
                    vx['U[VI]-F4'] = 0
                    vx['U[VII]-F4'] = 0
                elif not has_x_all:
                    vx[s[2:]] = 0  # Remove 'x_' prefix
            elif s.startswith('y_') and not has_y_all:
                vy[s[2:]] = 0  # Remove 'y_' prefix
        
        # Collect all rows of data
        data_rows = []
        current_row = {}
        
        with open(in_file, 'r') as f:
            for line in f:
                line = line.strip()
                
                # Extract specific data types
                if 'mol MSFL' in line and 'ni' in vx:
                    match = re.search(r'([\d.e-]+)\s*mol MSFL', line)
                    if match:
                        vx['ni'] = float(match.group(1))
                
                elif 'Moles of pairs' in line and 'nx' in vx:
                    match = re.search(r'([\d.e-]+)\s*Moles of pairs', line)
                    if match:
                        vx['nx'] = float(match.group(1))
                
                elif 'mol gas' in line and 'ny' in vy:
                    match = re.search(r'([\d.e-]+)\s*mol gas', line)
                    if match:
                        vy['ny'] = float(match.group(1))
                
                elif 'UF34soln' in line and extra_data['NUF34'] >= 0:
                    match = re.search(r'([\d.e-]+)\s*UF34soln', line)
                    if match:
                        extra_data['NUF34'] = float(match.group(1))
                
                elif 'Temperature' in line and extra_data['T'] >= 0:
                    match = re.search(r'([\d.e-]+)\s*Temperature', line)
                    if match:
                        extra_data['T'] = float(match.group(1))
                
                elif 'Pressure' in line and extra_data['P'] >= 0:
                    match = re.search(r'([\d.e-]+)\s*Pressure', line)
                    if match:
                        extra_data['P'] = float(match.group(1))
                
                # Collect row data
                if 'DEBUG' in line:
                    # Combine collected data
                    row_data = {}
                    for k, v in vx.items():
                        row_data[f'x_{k}' if k not in ['ni', 'nx'] else k] = v
                    for k, v in vy.items():
                        row_data[f'y_{k}' if k != 'ny' else k] = v
                    
                    # Add extra data
                    for k, v in extra_data.items():
                        if v >= 0:
                            row_data[k] = v
                    
                    # Reset data for next iteration
                    data_rows.append(row_data)
                    vx = {k: 0 for k in vx}
                    vy = {k: 0 for k in vy}
                    extra_data = {'NUF34': -1.0, 'T': -1.0, 'P': -1.0}
        
        # Convert to DataFrame and save
        df = pd.DataFrame(data_rows)
        df.to_csv(out_file, index=False)

    def merge_therm(in_files: List[str], out_file: str) -> None:
        """
        Merge multiple Thermochimica output files.
        
        Args:
            in_files (List[str]): List of input Thermochimica output files
            out_file (str): Path for the merged output file
        
        Note:
            Most accurate when files have uniform temperature and pressure,
            and contain only salt and ideal gas phases.
        """
        # Initialize accumulators for total moles
        nx, ny = 0.0, 0.0
        
        # Dictionaries to store component moles for salt and gas phases
        mx, my = {}, {}
        
        # Current processing context
        current_phase = None
        current_total_moles = 0.0
    
        for input_file in in_files:
            try:
                with open(input_file, 'r') as f:
                    for line in f:
                        line = line.strip()
                        
                        # Reset processing context when encountering phase delimiter
                        if '==' in line:
                            current_phase = None
                            continue
                        
                        # Skip empty lines
                        if not line:
                            continue
                        
                        # Detect phases and total moles
                        if 'pairs' in line:
                            current_phase = 'salt'
                            match = re.search(r'([\+\-]?\d+(?:\.\d+)?)\s*Moles of pairs', line)
                            if match:
                                n = float(match.group(1))
                                current_total_moles = abs(n)
                                nx += current_total_moles
                        
                        elif 'gas_ideal' in line:
                            current_phase = 'gas'
                            match = re.search(r'([\+\-]?\d+(?:\.\d+)?)\s*mol gas_ideal', line)
                            if match:
                                n = float(match.group(1))
                                current_total_moles = abs(n)
                                ny += current_total_moles
                        
                        # Process component moles
                        if current_phase and '{' in line:
                            continue
                        
                        if current_phase and '}' in line:
                            current_phase = None
                            current_total_moles = 0.0
                            continue
                        
                        # Parse component mole fractions
                        if current_phase:
                            # More robust parsing of component lines
                            match = re.match(r'([\+\-]?\d+(?:\.\d+)?)\s*(.+)', line)
                            if match:
                                value = float(match.group(1))
                                component = match.group(2).strip()
                                
                                if current_phase == 'salt':
                                    mx[component] = mx.get(component, 0) + value * current_total_moles
                                else:  # gas phase
                                    my[component] = my.get(component, 0) + value * current_total_moles
    
            except IOError as e:
                print(f"Error reading file {input_file}: {e}")
                continue
    
        # Write merged results
        with open(out_file, 'w') as output:
            # Salt phase
            if nx > 0:
                output.write(f"\t{nx} Moles of pairs\n\t {{\n")
                for component, total_moles in mx.items():
                    mole_fraction = total_moles / nx
                    if mole_fraction <= 1.0:
                        output.write(f"\t\t{mole_fraction:.6f}\t{component}\n")
                output.write("\t }\n")
    
            # Gas phase
            if ny > 0:
                output.write(f"\t{ny} mol gas_ideal\n\t {{\n")
                for component, total_moles in my.items():
                    mole_fraction = total_moles / ny
                    if mole_fraction <= 1.0:
                        output.write(f"\t\t{mole_fraction:.6f}\t{component}\n")
                output.write("\t }\n")

    @staticmethod
    def decouple_surr(scale_data: Dict[str, List[float]], 
                      thermo_file: str, 
                      out_file: str, 
                      calc_metals: bool = False, 
                      calc_hf: bool = False) -> None:
        """
        Decouple surrogate elements into their constituent elements.
        
        Args:
            scale_data (Dict[str, List[float]]): Element concentrations from SCALE
            thermo_file (str): Path to Thermochimica result file
            out_file (str): Path for the output file
            calc_metals (bool): Whether to calculate metal fluorides
            calc_hf (bool): Whether to calculate HF concentration
        """
        # Placeholder for implementation
        raise NotImplementedError("This method needs to be implemented")

def main():
    """
    Main entry point for data processing
    """
    print("Thermochimica Data Processing")
    print("1. Convert SCALE output to Thermochimica input files")
    print("2. Extract data from Thermochimica output")
    print("3. Merge Thermochimica output files")
    print("4. Decouple surrogate elements")
    
    option = int(input("Enter your choice (1-4): "))
    
    if option == 1:
        # Example usage for Option 1
        scale_file = input("Enter SCALE output file path: ")
        output_base = input("Enter base name for Thermochimica input files: ")
        temp_range = input("Enter temperature range/values: ")
        press_range = input("Enter pressure range/values: ")
        
        data = DataProcessor.scale_to_vector(scale_file)
        DataProcessor.vector_to_therm(data, output_base, temp_range, press_range)
    
    elif option == 2:
        # Placeholder for Option 2
        pass
    
    elif option == 3:
        # Placeholder for Option 3
        pass
    
    elif option == 4:
        # Placeholder for Option 4
        pass
    
    else:
        print("Invalid option selected.")

if __name__ == "__main__":
    main()
