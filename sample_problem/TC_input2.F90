program flibe
USE ModuleThermoIO
implicit none
cInputUnitTemperature = 'K'
cInputUnitPressure = 'atm'
cInputUnitMass = 'moles'
cThermoFileName = DATA_DIRECTORY // 'MSTDB-TC_V2.0_Fluorides_8-0.dat'
dPressure = 1
dTemperature = 900
dElementMass(9) = 146.6699964
dElementMass(3) = 63.32999366
dElementMass(4) = 31.78904499
dElementMass(92) = 4.999520995
dElementMass(42) = 7.81E-05
dElementMass(54) = 8.59E-05
dElementMass(94) = 0.000082
dElementMass(60) = 4.66E-05
dElementMass(55) = 5.65E-05
dElementMass(58) = 7.32E-05
dElementMass(44) = 5.79E-05
dElementMass(57) = 0.000081
dElementMass(43) = 1.60E-05
dElementMass(2) = 2.13E-05
dElementMass(36) = 1.54E-05
dElementMass(37) = 1.41E-05
dElementMass(45) = 1.81E-06
dElementMass(52) = 1.65E-05
dElementMass(46) = 3.94E-06
dElementMass(1) = 5.85E-06
dElementMass(53) = 0.000015
dElementMass(8) = 2.08E-06
dElementMass(34) = 2.09E-06
dElementMass(10) = 1.47E-06
dElementMass(41) = 1.44E-06
dElementMass(50) = 6.69E-07
dElementMass(48) = 2.50E-07
dElementMass(47) = 1.67E-07
dElementMass(51) = 5.96E-07
dElementMass(49) = 3.41E-08
dElementMass(7) = 2.66E-08
dElementMass(32) = 1.73E-08
dElementMass(33) = 1.22E-08
dElementMass(90) = 0.000125
dElementMass(31) = 6.95E-11
dElementMass(6) = 2.48E-11
dElementMass(30) = 5.22E-11
dElementMass(70) = 7.57E-14
dElementMass(5) = 1.54E-15
dElementMass(29) = 6.80E-13
dElementMass(82) = 2.65E-22
dElementMass(24) = 0.17
dElementMass(28) = 0.65
dElementMass(26) = 0.12
dElementMass(20) = 0.000102
iPrintResultsMode = 1
call ParseCSDataFile(cThermoFileName)
if (INFOThermo == 0) call Thermochimica
if (iPrintResultsMode > 0) call PrintResults
if (INFOThermo == 0) call ResetThermoAll
call ThermoDebug
end program flibe
