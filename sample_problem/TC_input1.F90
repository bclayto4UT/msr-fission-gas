program TC_input1
USE ModuleThermoIO
implicit none
cInputUnitTemperature = 'K'
cInputUnitPressure = 'atm'
cInputUnitMass = 'moles'
cThermoFileName = DATA_DIRECTORY // 'MSTDB-TC_V2.0_Fluorides_8-0.dat'
dPressure = 1
dTemperature = 900
dElementMass(94) = 8.2316e-005
dElementMass(70) = 7.5658e-014
dElementMass(58) = 7.3178e-005
dElementMass(57) = 8.05339e-005
dElementMass(20) = 0.000102004
dElementMass(6) = 2.4779e-011
dElementMass(92) = 4.9995
dElementMass(46) = 3.9431e-006
dElementMass(26) = 650
dElementMass(37) = 1.4102e-005
dElementMass(44) = 6.3014e-005
dElementMass(23) = 1.0428e-005
dElementMass(8) = 2.0754e-006
dElementMass(22) = 3.6409e-008
dElementMass(31) = 6.9484e-011
dElementMass(4) = 31.789
dElementMass(5) = 1.5431e-015
dElementMass(25) = 4.3785e-006
dElementMass(34) = 2.0885e-006
dElementMass(90) = 0.00012534
dElementMass(42) = 30
dElementMass(27) = 30
dElementMass(52) = 1.6511e-005
dElementMass(89) = 5.7178e-018
dElementMass(2) = 2.2049e-005
dElementMass(3) = 63.33
dElementMass(1) = 1.9406e-005
dElementMass(10) = 1.4665e-006
dElementMass(28) = 120
dElementMass(43) = 2.8599e-005
dElementMass(7) = 2.66e-008
dElementMass(33) = 1.2169e-008
dElementMass(9) = 146.67
dElementMass(47) = 1.6723e-007
dElementMass(55) = 5.6524e-005
dElementMass(24) = 170
dElementMass(32) = 1.7329e-008
dElementMass(41) = 1.4702e-006
dElementMass(29) = 1.275e-006
dElementMass(45) = 1.8057e-006
dElementMass(30) = 5.3516e-011
dElementMass(51) = 5.9592e-007
dElementMass(48) = 2.4971e-007
dElementMass(60) = 4.6598e-005
dElementMass(54) = 8.5917e-005
dElementMass(36) = 1.5367e-005
dElementMass(49) = 3.4112e-008
dElementMass(50) = 6.6905e-007
dElementMass(21) = 5.8769e-017
dElementMass(53) = 1.50784e-005
iPrintResultsMode = 1
call ParseCSDataFile(cThermoFileName)
if (INFOThermo == 0) call Thermochimica
if (iPrintResultsMode > 0) call PrintResults
if (INFOThermo == 0) call ResetThermoAll
call ThermoDebug
end program TC_input1
