program flibe
USE ModuleThermoIO
implicit none
cInputUnitTemperature = 'K'
cInputUnitPressure = 'atm'
cInputUnitMass = 'moles'
cThermoFileName = DATA_DIRECTORY // 'MSTDB-TC_V2.0_Fluorides_8-0.dat'
dPressure = 1
dTemperature = 900
dElementMass(9) = 146.6699893
dElementMass(3) = 63.32998098
dElementMass(4) = 31.78903496
dElementMass(92) = 4.998563119
dElementMass(42) = 0.000224216
dElementMass(54) = 0.000260237
dElementMass(94) = 0.000247
dElementMass(60) = 0.000155006
dElementMass(55) = 0.000200439
dElementMass(58) = 0.000220896
dElementMass(44) = 0.000168384
dElementMass(57) = 0.000243
dElementMass(43) = 6.48E-05
dElementMass(2) = 6.39E-05
dElementMass(36) = 4.53E-05
dElementMass(37) = 4.25E-05
dElementMass(45) = 8.93E-06
dElementMass(52) = 3.52E-05
dElementMass(46) = 1.37E-05
dElementMass(1) = 1.75E-05
dElementMass(53) = 0.000026
dElementMass(8) = 6.23E-06
dElementMass(34) = 6.27E-06
dElementMass(10) = 4.40E-06
dElementMass(41) = 9.48E-06
dElementMass(50) = 1.92E-06
dElementMass(48) = 8.03E-07
dElementMass(47) = 4.56E-07
dElementMass(51) = 1.12E-06
dElementMass(49) = 1.30E-07
dElementMass(7) = 7.98E-08
dElementMass(32) = 4.62E-08
dElementMass(33) = 2.09E-08
dElementMass(90) = 0.000376
dElementMass(31) = 1.57E-10
dElementMass(6) = 7.44E-11
dElementMass(30) = 9.76E-11
dElementMass(70) = 4.30E-13
dElementMass(5) = 1.39E-14
dElementMass(29) = 7.30E-13
dElementMass(82) = 7.20E-20
dElementMass(24) = 0.17
dElementMass(28) = 0.65
dElementMass(26) = 0.12
dElementMass(20) = 0.000276
iPrintResultsMode = 1
call ParseCSDataFile(cThermoFileName)
if (INFOThermo == 0) call Thermochimica
if (iPrintResultsMode > 0) call PrintResults
if (INFOThermo == 0) call ResetThermoAll
call ThermoDebug
end program standard_flibe
