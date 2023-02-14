program test4
USE ModuleThermoIO
implicit none
cInputUnitTemperature = 'K'
cInputUnitPressure = 'atm'
cInputUnitMass = 'moles'
cThermoFileName = DATA_DIRECTORY // 'MSTDB-TC_V2.0_Fluorides_8-0.dat'
dPressure = 1
dTemperature = 900
dElementMass(20) = 0.00035486
dElementMass(94) = 0.000329316
dElementMass(70) = 6.1732e-013
dElementMass(58) = 0.00029547
dElementMass(57) = 0.000325016
dElementMass(6) = 9.9215e-011
dElementMass(92) = 4.9981
dElementMass(46) = 1.8654e-005
dElementMass(26) = 0.65
dElementMass(37) = 5.6765e-005
dElementMass(44) = 0.0002216
dElementMass(23) = 1.2954e-007
dElementMass(8) = 8.3017e-006
dElementMass(22) = 1.4567e-010
dElementMass(31) = 2.0133e-010
dElementMass(4) = 31.789
dElementMass(5) = 2.4625e-014
dElementMass(25) = 1.9831e-008
dElementMass(34) = 8.36e-006
dElementMass(90) = 0.00050105
dElementMass(42) = 0.0002983
dElementMass(27) = 3.2903e-008
dElementMass(52) = 4.4108e-005
dElementMass(32) = 6.0576e-008
dElementMass(24) = 0.17
dElementMass(89) = 1.1414e-016
dElementMass(2) = 8.5243e-005
dElementMass(3) = 63.33
dElementMass(1) = 2.3385e-005
dElementMass(10) = 5.8661e-006
dElementMass(28) = 0.12
dElementMass(43) = 8.9615e-005
dElementMass(7) = 1.0633e-007
dElementMass(33) = 2.5165e-008
dElementMass(9) = 146.67
dElementMass(47) = 5.9298e-007
dElementMass(55) = 0.00027493
dElementMass(41) = 1.52e-005
dElementMass(29) = 5.211e-009
dElementMass(45) = 1.4496e-005
dElementMass(30) = 1.1994e-010
dElementMass(51) = 1.3596e-006
dElementMass(48) = 1.0902e-006
dElementMass(60) = 0.00021368
dElementMass(54) = 0.00034702
dElementMass(36) = 6.0272e-005
dElementMass(49) = 1.7848e-007
dElementMass(50) = 2.5313e-006
dElementMass(21) = 7.8529e-019
dElementMass(53) = 3.0339e-005
iPrintResultsMode = 1
call ParseCSDataFile(cThermoFileName)
if (INFOThermo == 0) call Thermochimica
if (iPrintResultsMode > 0) call PrintResults
if (INFOThermo == 0) call ResetThermoAll
call ThermoDebug
end program test4
