program test3
USE ModuleThermoIO
implicit none
cInputUnitTemperature = 'K'
cInputUnitPressure = 'atm'
cInputUnitMass = 'moles'
cThermoFileName = DATA_DIRECTORY // 'MSTDB-TC_V2.0_Fluorides_8-0.dat'
dPressure = 1
dTemperature = 900
dElementMass(20) = 0.00027573
dElementMass(94) = 0.000247089
dElementMass(70) = 4.3008e-013
dElementMass(58) = 0.0002209
dElementMass(57) = 0.000243907
dElementMass(6) = 7.4387e-011
dElementMass(92) = 4.9986
dElementMass(46) = 1.3706e-005
dElementMass(26) = 0.65
dElementMass(37) = 4.2548e-005
dElementMass(44) = 0.00016838
dElementMass(23) = 7.8549e-008
dElementMass(8) = 6.2263e-006
dElementMass(22) = 1.0924e-010
dElementMass(31) = 1.5748e-010
dElementMass(4) = 31.789
dElementMass(5) = 1.3864e-014
dElementMass(25) = 1.4291e-008
dElementMass(34) = 6.27e-006
dElementMass(90) = 0.00037592
dElementMass(42) = 0.00022422
dElementMass(27) = 2.5397e-008
dElementMass(52) = 3.5248e-005
dElementMass(32) = 4.6163e-008
dElementMass(24) = 0.17
dElementMass(89) = 6.2653e-017
dElementMass(2) = 6.3923e-005
dElementMass(3) = 63.33
dElementMass(1) = 1.7555e-005
dElementMass(10) = 4.3996e-006
dElementMass(28) = 0.12
dElementMass(43) = 6.4837e-005
dElementMass(7) = 7.975e-008
dElementMass(33) = 2.088e-008
dElementMass(9) = 146.67
dElementMass(47) = 4.5605e-007
dElementMass(55) = 0.00020044
dElementMass(41) = 9.4832e-006
dElementMass(29) = 3.8906e-009
dElementMass(45) = 8.9291e-006
dElementMass(30) = 9.7606e-011
dElementMass(51) = 1.1177e-006
dElementMass(48) = 8.0338e-007
dElementMass(60) = 0.00015501
dElementMass(54) = 0.00026024
dElementMass(36) = 4.5311e-005
dElementMass(49) = 1.2987e-007
dElementMass(50) = 1.9167e-006
dElementMass(21) = 4.7478e-019
dElementMass(53) = 2.63827e-005
iPrintResultsMode = 1
call ParseCSDataFile(cThermoFileName)
if (INFOThermo == 0) call Thermochimica
if (iPrintResultsMode > 0) call PrintResults
if (INFOThermo == 0) call ResetThermoAll
call ThermoDebug
end program test3
