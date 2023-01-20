# Read Me
This is the entirety of the code I used to update the molten salt's chemical composition. SCALE and Thermochimica is required. See [this page](https://github.com/ORNL-CEES/thermochimica) for installing Thermochimica.

## Using the code
All code files are under the dataProcessor folder. Most of them are numerical methods file that Dr. Pencheva and I wrote while in M 348 and M 368K (I wasn't able to install a library). 

The project can be compiled in any way you'd like. There should be a new executable file if compilation is sucessful. I liked to copy that executable file (named ``dataProcessor.exe``) to the ``thermochimica`` directory so I could run both my code and Thermochimica on Ubuntu. If the executable file cannot be opened on Ubuntu, because of a denied permission, enter
```
chmod u+x ./dataProcessor.exe # or however it is named
```
The program should look like this:
![image](https://user-images.githubusercontent.com/62024926/213612517-f9284786-c0ef-4fd8-aa60-9f4ef6022e29.png)

It has 4 options:

* Option 1: Convert a SCALE output file (.F71 table) into multiple TC input files, with each file representing the system at one time instance. A sample SCALE output file is shown in ``sample_problem/F71_output.txt``, where the elements are on the columns and the time intervals are on the rows. The answer to the prompt "Do the columns contain element symbols?" is "Y". If the reserve is true, then "N" is the answer. The system temperature (in K) and pressure (in atm) are also needed. The TC input files have the .F90 ending and should be put in the ``thermochimica/test`` directory.
* Option 2: Extract numerical data from one or multiple Thermochimica results and arrange them in table format. The accepted quantities are the total moles of ions (ni), total moles of salt (nx), total moles of gas (ny), temperature (T), mole fraction of species ABC in the salt (x_ABC), and mole fraction of species ABC in the gas (y_ABC). For example, entering "nx x_UF3 x_UF4 y_UF5" will create a table with the total moles of salt, the mole fractions of UF3 and UF4 in the salt, and the mole fraction of UF5 in the gas as a function of time.
* Option 3: Merge two TC output files into one. This function has not been updated, but is left in the code as it may be useful in future projects.
* Option 4: Decouple surrogate elements into the chemically similar elements not accounted for by MSTDB (for example, Ca into Ca, Ba, and Sr). This function requires a SCALE .F71 output (just as Option 1) and a TC result file. Optionally, it can also perform thermochemical calculations of corrosion products (Cr, Fe, and Ni) as well as HF concentration. This requires a slight modification in the .F71 file, so that the entries under Cr, Fe, and Ni represent the mole fraction of the respective metal in the alloy (as opposed to amounts in moles). The result file will be in a similar format as a TC output, so that data can be extracted through Option 2.

## Sample problem

1. Create an ORIGEN file with the appropriate starting fuel composition, flux, and irradiation time. For my project, the fuel is 5% UF4 in 2LiF-BeF2. The uranium is enriched to 20 weight-% U-235 and the lithium to 99.99 mol-% Li-7, so the specification should be:

```
    mat{
        units = MOLES
        iso = [li6 = 0.00633
               li7 = 63.32367
               be9 = 31.78905
               f19 = 146.67
               u235 = 0.98989
               u238 = 4.01011]
        }
 ```
        
2. Run the ORIGEN file, open the result .F71 file, go to Table, filter the data to display elements (preferably all elements, but an nrank of about 50-60 should be fine) in moles.
3. Copy the table to Excel and crop out the "Subtotals" and "Totals" columns. If needed, the values of the elements Cu, Fe, and Zn can be modified into their respective mole fraction in the structural metal alloy (should be the same throughout the irradiation time). Copy the table into a text file.
4. Run main.cpp (Option 1) on the text file to convert it into TC input files. There should be multiple TC inputs corresponding to however many time intervals there are.
5. Run TC on all the input files and have the results printed out in one single result file (Use a loop).
```
make
for i in $(seq 1 1 184) # the last number should be how many TC input files there are
do
./bin/TC_input$i >> TC_first_results.txt
done
```
6. Run main.cpp (Option 4) on the TC result file to decouple the surrogate elements into the actual elements they represent. Optionally, it can also calculate fission products and HF composition. This should give you a new, updated result file in roughly the same format as TC.
7. Run main.cpp (Option 2) on the result file to extract any value of interest.
