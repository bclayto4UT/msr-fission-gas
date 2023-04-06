program TC_input5
USE ModuleThermoIO
implicit none
cInputUnitTemperature = 'K'
cInputUnitPressure = 'atm'
cInputUnitMass = 'moles'
cThermoFileName = DATA_DIRECTORY // 'MSTDB-TC_V2.0_Fluorides_8-0.dat'
dPressure = 1
dTemperature = 900
dElementMass(94) = 0.000411439
dElementMass(70) = 8.0599e-013
dElementMass(58) = 0.00036921
dElementMass(57) = 0.000406438
dElementMass(20) = 0.00043092
dElementMass(6) = 1.2406e-010
dElementMass(92) = 4.9976
dElementMass(46) = 2.3644e-005
dElementMass(26) = 650
dElementMass(37) = 7.098e-005
dElementMass(44) = 0.00029961
dElementMass(23) = 0.00018856
dElementMass(8) = 1.0377e-005
dElementMass(22) = 1.8209e-007
dElementMass(31) = 2.4521e-010
dElementMass(4) = 31.789
dElementMass(5) = 3.8443e-014
dElementMass(25) = 2.5758e-005
dElementMass(34) = 1.0449e-005
dElementMass(90) = 0.000626231
dElementMass(42) = 30
dElementMass(27) = 30
dElementMass(52) = 5.292e-005
dElementMass(89) = 1.8092e-016
dElementMass(2) = 0.00011035
dElementMass(3) = 63.33
dElementMass(1) = 9.7124e-005
dElementMass(10) = 7.3326e-006
dElementMass(28) = 120
dElementMass(43) = 0.00020446
dElementMass(7) = 1.329e-007
dElementMass(33) = 2.945e-008
dElementMass(9) = 146.67
dElementMass(47) = 7.2977e-007
dElementMass(55) = 0.00034963
dElementMass(24) = 170
dElementMass(32) = 7.4987e-008
dElementMass(41) = 2.1619e-005
dElementMass(29) = 6.5393e-006
dElementMass(45) = 2.1118e-005
dElementMass(30) = 1.7625e-010
dElementMass(51) = 1.6017e-006
dElementMass(48) = 1.379e-006
dElementMass(60) = 0.00027414
dElementMass(54) = 0.00043391
dElementMass(36) = 7.5225e-005
dElementMass(49) = 2.2729e-007
dElementMass(50) = 3.1437e-006
dElementMass(21) = 1.1439e-015
dElementMass(53) = 3.39921e-005
iPrintResultsMode = 1
call ParseCSDataFile(cThermoFileName)
if (INFOThermo == 0) call Thermochimica
if (iPrintResultsMode > 0) call PrintResults
if (INFOThermo == 0) call ResetThermoAll
call ThermoDebug
end program TC_input5
