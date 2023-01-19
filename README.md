# Read Me
This is the entirety of the code I used to update the molten salt's chemical composition. SCALE and Thermochimica is required.
The main.cpp should display a few options to extract and process data from input files. Refer to dataProcessor.cpp for any details as to how each function works.
The steps are detailed below:

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
for i in $(seq 1 1 5)
do
./bin/TC_input$i >> TC_first_results.txt
done
```
6. Run main.cpp (Option 4) on the TC result file to decouple the surrogate elements into the actual elements they represent. Optionally, it can also calculate fission products and HF composition. This should give you a new, updated result file in roughly the same format as TC.
7. Run main.cpp (Option 2) on the result file to extract any value of interest.
