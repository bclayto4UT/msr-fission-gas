# Python Implementation Plan for Thermochimica Data Processing

This plan outlines a structured approach for converting the C++ Thermochimica data processing code to Python, leveraging existing Python libraries for numerical methods and data processing. Each section includes specific information about what original code should be referenced during implementation.

## 1. Core Data Structures and Utility Module (`utils.py`)

### Implementation with Libraries:
- Use `pandas` for data processing and file handling
- Use `mendeleev` or `periodictable` library for basic element properties
- Implement custom surrogate mappings and specialized utility functions

### Required Functions:
- `element_symbol(str)`: Convert element names to proper case
- `str_to_list(data, delimiter)`: Split string to list (use Python's built-in `.split()`)
- `str_to_float_list(str)`: Convert string to list of floats with support for multiple formats

### Element Data Dictionaries:
- `surrogate_map`: Surrogate element groupings (must be implemented manually)
- `surrogate_map_inv`: Individual elements → surrogate representative (must be implemented manually)

### Code References Needed:
- `thermoElectroChem.h`: For surrogate mappings and element groupings
- `miscellaneous.cpp/h`: For string processing functions and utility methods

## 2. Thermochemistry Module (`thermochemistry.py`)

### Implementation with Libraries:
- Use NumPy for mathematical operations
- Implement custom thermodynamic functions specific to molten salt chemistry

### Required Data Structures:
- `heat_data`: Dictionary storing thermodynamic data for compounds

### Required Functions:
- Core thermodynamic calculations (may need custom implementation):
  - `calc_enthalpy(data, T)`
  - `calc_entropy(data, T)` 
  - `calc_gibbs_energy(data, T)`
- Gibbs free energy functions for fluorides:
  - `g_hf(T)`, `g_uf4(T)`, `g_f(x_uf3, x_uf4, T)`
  - Metal fluoride formation energies (Cr, Fe, Co, Ni, etc.)

### Code References Needed:
- `thermoElectroChem.cpp/h`: Complete file for thermodynamic functions and constants
- Chemical formulas and constants from the original implementation

## 3. Numerical Methods Module (`numerical.py`)

### Implementation with Libraries:
- Replace most custom numerical methods with SciPy equivalents:
  - `scipy.optimize.root_scalar()` for root finding
  - `scipy.optimize.root()` for nonlinear systems
  - `scipy.integrate` for numerical integration
  - NumPy for linear algebra operations

### Required Wrappers (thin interfaces to SciPy):
```python
# Root finding wrappers
def bisection(f, x0, x1, max_iter=20, tol=1e-6):
    return scipy.optimize.root_scalar(f, method='bisect', bracket=[x0, x1], 
                                     options={'maxiter': max_iter, 'xtol': tol}).root

def newton(f, df=None, x0=1.0, max_iter=20, tol=1e-6):
    if df is None:
        # Use numerical derivative
        return scipy.optimize.root_scalar(f, method='secant', x0=x0, x1=x0*1.0001, 
                                         options={'maxiter': max_iter, 'xtol': tol}).root
    else:
        # Use analytical derivative
        return scipy.optimize.root_scalar(f, method='newton', fprime=df, x0=x0, 
                                         options={'maxiter': max_iter, 'xtol': tol}).root

# Nonlinear system wrappers
def newton_vector(f, df=None, x0=None, max_iter=20, tol=1e-6):
    if df is None:
        # Use numerical Jacobian
        return scipy.optimize.root(f, x0, method='hybr', 
                                  options={'maxfev': max_iter, 'xtol': tol}).x
    else:
        # Use analytical Jacobian
        return scipy.optimize.root(f, x0, method='hybr', jac=df, 
                                  options={'maxfev': max_iter, 'xtol': tol}).x
```

### Specialized Functions (if needed):
- `del_function(zeta)`: Calculate UF4 reduction to UF3
- `thermo_func(zeta)`: Calculate fluoride amounts

### Code References Needed:
- `rootFinding.cpp/h`: For root-finding interface details
- `iterativeNL.cpp/h`: For nonlinear system solver interface details
- `thermoElectroChem.cpp`: For specialized thermodynamic functions

## 4. Data Processor Core Module (`data_processor.py`)

### Implementation with Libraries:
- Use `pandas` for structured data handling
- Use NumPy for array manipulations
- Use Python's file handling capabilities for I/O

### Required Functions:
- `scale_to_vector(in_file)`: Parse SCALE output files
  ```python
  def scale_to_vector(in_file):
      """
      Parse SCALE output file to extract element concentration data.
      
      Args:
          in_file (str): Path to the SCALE output file
          
      Returns:
          dict: Dictionary mapping element names to concentration vectors
      """
      # Use pandas to read and process the file
      df = pd.read_csv(in_file, delim_whitespace=True)
      # Process data into appropriate format
      # Return dictionary of element concentrations
  ```

- `vector_to_therm(data_map, out_file, str_t, str_p, includes_surr=True)`:
  ```python
  def vector_to_therm(data_map, out_file, str_t, str_p, includes_surr=True):
      """
      Convert concentration data to Thermochimica input files (.F90)
      
      Args:
          data_map (dict): Dictionary with element concentration data
          out_file (str): Base name for output files
          str_t (str): Temperature specification string
          str_p (str): Pressure specification string
          includes_surr (bool): Whether to include surrogate elements
          
      Returns:
          None: Files are written to disk
      """
      # Process temperature and pressure strings
      # Generate Fortran (.F90) files for each time step
  ```

- `text_to_excel(in_file, out_file, data_type)`:
  ```python
  def text_to_excel(in_file, out_file, data_type):
      """
      Extract specified data from Thermochimica output.
      
      Args:
          in_file (str): Path to Thermochimica output file
          out_file (str): Path for the CSV output file
          data_type (str): Space-separated string of data types to extract
          
      Returns:
          None: Extracted data is written to CSV file
      """
      # Parse data_type string to determine what to extract
      # Read and parse the Thermochimica output file
      # Extract required data and save to CSV
  ```

- `merge_therm(in_files, out_file)` and `decouple_surr(...)`:
  - Implement with similar patterns using pandas for data processing

### Code References Needed:
- `dataProcessor.cpp/h`: Complete file for full function signatures and processing logic
- Sample input/output files to understand data formats
- Any specific formats for Thermochimica input/output files

## 5. Main Application (`main.py`)

### Implementation with Libraries:
- Use `click` or `argparse` for command-line interface
- Use Python's built-in functions for file operations

### Required Functions:
```python
def main():
    """Main entry point for the application."""
    print("Thermochimica Data Processing")
    print("1. Convert SCALE output to Thermochimica input files")
    print("2. Extract data from Thermochimica output")
    print("3. Merge Thermochimica output files")
    print("4. Decouple surrogate elements")
    
    option = int(input("Enter your choice (1-4): "))
    process_option(option)

def process_option(option):
    """Process the user's menu selection."""
    if option == 1:
        run_option_1()
    elif option == 2:
        run_option_2()
    # etc.
```

### Code References Needed:
- `main.cpp`: For menu structure and program flow
- Sample commands and interactions from the original application

## Implementation Strategy and Library Dependencies

### Required Libraries:
- `numpy`: For numerical operations and array handling
- `scipy`: For optimization, root-finding, and integration
- `pandas`: For data processing and file I/O
- `mendeleev` or `periodictable`: For element properties (optional)
- `click` or `argparse`: For command-line interface (optional)

### Implementation Order:
1. Set up development environment with required libraries
2. Implement `utils.py` with element data and utility functions
3. Implement `thermochemistry.py` with molten salt thermodynamics
4. Implement `numerical.py` with SciPy wrappers
5. Implement `data_processor.py` core functionality
6. Implement `main.py` with user interface
7. Test with sample data from original application

### Testing Strategy:
- Compare outputs with original C++ implementation using identical inputs
- Start with simple cases and progress to more complex scenarios
- Validate each module independently before integration

## Sample Data and Expected Results

For effective testing, you should have:
1. Sample SCALE output files (like `msr_isotopics.out` mentioned in the documentation)
2. Sample Thermochimica input and output files
3. Expected results from the original C++ implementation for comparison

## Special Considerations

1. **Fortran File Generation**: The Python code needs to generate valid Fortran (.F90) files for Thermochimica. Pay special attention to the format of these files.

2. **File Formats**: Ensure the Python implementation can handle the same input/output formats as the original C++ code.

3. **Numeric Precision**: Be mindful of potential differences in numeric precision between C++ and Python implementations.

4. **Performance**: If performance is critical, consider using NumPy's vectorized operations where possible.
