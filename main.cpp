#include <fstream>
#include <iostream>

#include "dataProcessor.h"

using namespace std;

int main()
{
    string inFile, outFile;
    int option = -1;
    do{
        do{
            cout << "Select option: " << endl;
            cout << "1: Convert SCALE output to Thermochimica input." << endl;
            cout << "2: Extract data from Thermochimica output." << endl;
            cout << "3: Combine multiple Thermochimica outputs." << endl;
            cout << "4: Convert mass fractions into moles." << endl;
            cout << "0: Quit." << endl;
            cin >> option;

            if (option == 1){
                cout << "Enter input file: ";
                cin >> inFile;
                cout << "Enter output file: ";
                cin >> outFile;

                try{
                    auto v2 = scaleToVector(inFile);
                    vectToThermI(v2, outFile);
                    cout << "Successful!" << endl << endl;
                } catch (const invalid_argument& ex){
                    cerr << ex.what() << endl;
                } catch (...){
                    cerr << "Error encountered!" << endl;
                }

            } else if (option == 2){
                cout << "Enter input file: ";
                cin >> inFile;
                cout << "Enter output file: ";
                cin >> outFile;
                string dataType;
                cout << "Enter the data to be extracted: ";
                getline(cin >> ws, dataType);

                try{
                    textToExcel(inFile, outFile, dataType);
                    cout << "Successful!" << endl << endl;
                } catch (const invalid_argument& ex){
                    cerr << ex.what() << endl;
                } catch (...){
                    cerr << "Error encountered!" << endl;
                }

            } else if (option == 3){
                int num;
                do{
                    cout << "Enter how many files to combine: ";
                    cin >> num;
                } while (num <= 0);

                cout << "Enter " << num << " file names, separated by Enter:" << endl;
                strVect inFiles;
                for (int i = 0; i < num; i++){
                    getline(cin >> ws, inFile);
                    inFiles.push_back(inFile);
                }

                cout << "Enter output file: ";
                cin >> outFile;

                try{
                    mergeTherm(inFiles, outFile);
                    cout << "Successful!" << endl << endl;
                } catch (const invalid_argument& ex){
                    cerr << ex.what() << endl;
                } catch (...){
                    cerr << "Error encountered!" << endl;
                }

            } else if (option == 4){
                map<string, double> m;
                string ele; double eleMass;
                double totalMass = 0;
                bool isTotal = false;

                do{
                    cout << "Is the mass absolute or relative? Enter A or R: ";
                    cin >> ele;
                } while (ele != "a" && ele != "A" && ele != "r" && ele != "R");

                if (ele == "a" || ele == "A") totalMass = 1;
                else if (ele == "r" || ele == "R") isTotal = true;

                cout << "Enter element symbols and its mass. Once done, enter QUIT 0." << endl;
                while (ele != "QUIT"){
                    cin >> ele >> eleMass;
                    if (ele == "QUIT") break;
                    m[ele] = eleMass;
                    if (isTotal) totalMass += eleMass;
                }

                cout << "Enter output file: ";
                cin >> outFile;

                try{
                    massToMole(outFile, m, totalMass);
                    cout << "Successful!" << endl << endl;
                } catch (const out_of_range& ex){
                    cerr << ex.what() << endl;
                } catch (...){
                    cerr << "Error encountered!" << endl;
                }
            }

        } while (option < 0 || option > 4);
    } while (option != 0);

    return 0;
}
