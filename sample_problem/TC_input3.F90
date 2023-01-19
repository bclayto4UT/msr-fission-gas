program flibe
USE ModuleThermoIO
implicit none
cInputUnitTemperature = 'K'
cInputUnitPressure = 'atm'
cInputUnitMass = 'moles'
cThermoFileName = DATA_DIRECTORY // 'MSTDB-TC_V2.0_Fluorides_8-0.dat'
dPressure = 1
dTemperature = 900
dElementMass(9) = 146.6699929
dElementMass(3) = 63.32998732
dElementMass(4) = 31.78903997
dElementMass(92) = 4.999041966
dElementMass(42) = 0.000151108
dElementMass(54) = 0.000173435
dElementMass(94) = 0.000165
dElementMass(60) = 9.89E-05
dElementMass(55) = 0.000126815
dElementMass(58) = 0.000146282
dElementMass(44) = 0.000113883
dElementMass(57) = 0.000162
dElementMass(43) = 4.01E-05
dElementMass(2) = 4.26E-05
dElementMass(36) = 3.03E-05
dElementMass(37) = 2.83E-05
dElementMass(45) = 4.62E-06
dElementMass(52) = 2.63E-05
dElementMass(46) = 8.80E-06
dElementMass(1) = 1.17E-05
dElementMass(53) = 0.000022
dElementMass(8) = 4.15E-06
dElementMass(34) = 4.18E-06
dElementMass(10) = 2.93E-06
dElementMass(41) = 4.74E-06
dElementMass(50) = 1.30E-06
dElementMass(48) = 5.21E-07
dElementMass(47) = 3.16E-07
dElementMass(51) = 8.73E-07
dElementMass(49) = 8.15E-08
dElementMass(7) = 5.32E-08
dElementMass(32) = 3.17E-08
dElementMass(33) = 1.66E-08
dElementMass(90) = 0.000251
dElementMass(31) = 1.14E-10
dElementMass(6) = 4.96E-11
dElementMass(30) = 7.53E-11
dElementMass(70) = 2.45E-13
dElementMass(5) = 6.17E-15
dElementMass(29) = 7.26E-13
dElementMass(82) = 9.73E-21
dElementMass(24) = 0.17
dElementMass(28) = 0.65
dElementMass(26) = 0.12
dElementMass(20) = 0.000192
iPrintResultsMode = 1
call ParseCSDataFile(cThermoFileName)
if (INFOThermo == 0) call Thermochimica
if (iPrintResultsMode > 0) call PrintResults
if (INFOThermo == 0) call ResetThermoAll
call ThermoDebug
end program standard_flibe
