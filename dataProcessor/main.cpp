#include <fstream>
#include <iostream>
#include <iomanip>

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
            cout << "4: Decouple surrogate elements." << endl;
            cout << "0: Quit." << endl;
            cin >> option;

            if (option == 1){
                cout << "Enter input file (.out): ";
                cin >> inFile;
                cout << "Enter output file (.F90): ";
                cin >> outFile;

                try{
                    auto v2 = scaleToVector(inFile);
                    string strT, strP;
                    cout << "Enter temperature (K): ";
                    getline(cin >> ws, strT);
                    cout << "Enter pressure (atm): ";
                    getline(cin >> ws, strP);

                    vectToTherm(v2, outFile, strT, strP);
                    cout << "Successful!" << endl << endl;
                } catch (const invalid_argument& ex){
                    cerr << ex.what() << endl;
                } catch (const out_of_range& ex){
                    cerr << ex.what() << endl;
                } catch (const bad_alloc& ex){
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
                string inFile, thermoIn, thermoOut, ansYN;
                try{
                    cout << "Enter SCALE output file: ";
                    cin >> inFile;
                    auto v = scaleToVector(inFile);

                    cout << "Enter original Thermochimica output file of the above SCALE output: ";
                    cin >> thermoIn;
                    cout << "Enter new output file: ";
                    cin >> thermoOut;

                    cout << "The default setting excludes thermochemical calculations\n";
                    cout << "for any solid solution or fission products present in the system.\n";
                    do{
                        cout << "Is that how you wish to proceed? Enter Y/N: ";
                        cin >> ansYN;
                        try{
                            if (ansYN == "N" || ansYN == "n"){
                                bool includesSS, includesHF;
                                cout << "Include solid solution calculation? Type Y/N: ";
                                cin >> ansYN;
                                includesSS = (ansYN == "Y" || ansYN == "y");

                                cout << "Include fission product calculation? Type Y/N: ";
                                cin >> ansYN;
                                includesHF = (ansYN == "Y" || ansYN == "y");
                                decoupleSurr(v, thermoIn, thermoOut, includesSS, includesHF);
                                cout << "Successful!" << endl << endl;

                            } else{
                                decoupleSurr(v, thermoIn, thermoOut);
                                cout << "Successful!" << endl << endl;
                            }
                        } catch (...){
                            cerr << "Error encountered!" << endl;
                            continue;
                        }
                    } while (ansYN != "N" && ansYN != "n" && ansYN != "Y" && ansYN != "y");
                } catch(const invalid_argument& ex){
                    cerr << ex.what() << endl;
                    continue;
                } catch(...){
                    cerr << "Error encountered!" << endl;
                    continue;
                }
            }

        } while (option < 0 || option > 4);
    } while (option != 0);

    return 0;
}
