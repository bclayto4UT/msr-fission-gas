program flibe
USE ModuleThermoIO
implicit none
cInputUnitTemperature = 'K'
cInputUnitPressure = 'atm'
cInputUnitMass = 'moles'
cThermoFileName = DATA_DIRECTORY // 'MSTDB-TC_V2.0_Fluorides_8-0.dat'
dPressure = 1
dTemperature = 900
dElementMass(9) = 146.6699857
dElementMass(3) = 63.32997464
dElementMass(4) = 31.78902995
dElementMass(92) = 4.998084459
dElementMass(42) = 0.000298299
dElementMass(54) = 0.000347024
dElementMass(94) = 0.000329
dElementMass(60) = 0.00021368
dElementMass(55) = 0.000274931
dElementMass(58) = 0.000295466
dElementMass(44) = 0.000221603
dElementMass(57) = 0.000325
dElementMass(43) = 8.96E-05
dElementMass(2) = 8.52E-05
dElementMass(36) = 6.03E-05
dElementMass(37) = 5.68E-05
dElementMass(45) = 1.45E-05
dElementMass(52) = 4.41E-05
dElementMass(46) = 1.87E-05
dElementMass(1) = 2.33E-05
dElementMass(53) = 0.000030
dElementMass(8) = 8.30E-06
dElementMass(34) = 8.36E-06
dElementMass(10) = 5.87E-06
dElementMass(41) = 1.52E-05
dElementMass(50) = 2.53E-06
dElementMass(48) = 1.09E-06
dElementMass(47) = 5.93E-07
dElementMass(51) = 1.36E-06
dElementMass(49) = 1.78E-07
dElementMass(7) = 1.06E-07
dElementMass(32) = 6.06E-08
dElementMass(33) = 2.52E-08
dElementMass(90) = 0.000501
dElementMass(31) = 2.01E-10
dElementMass(6) = 9.92E-11
dElementMass(30) = 1.20E-10
dElementMass(70) = 6.17E-13
dElementMass(5) = 2.46E-14
dElementMass(29) = 7.31E-13
dElementMass(82) = 2.82E-19
dElementMass(24) = 0.17
dElementMass(28) = 0.65
dElementMass(26) = 0.12
dElementMass(20) = 0.000355
iPrintResultsMode = 1
call ParseCSDataFile(cThermoFileName)
if (INFOThermo == 0) call Thermochimica
if (iPrintResultsMode > 0) call PrintResults
if (INFOThermo == 0) call ResetThermoAll
call ThermoDebug
end program flibe
