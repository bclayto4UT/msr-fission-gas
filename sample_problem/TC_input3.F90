program TC_input3
USE ModuleThermoIO
implicit none
cInputUnitTemperature = 'K'
cInputUnitPressure = 'atm'
cInputUnitMass = 'moles'
cThermoFileName = DATA_DIRECTORY // 'MSTDB-TC_V2.0_Fluorides_8-0.dat'
dPressure = 1
dTemperature = 900
dElementMass(94) = 0.000247089
dElementMass(70) = 4.3008e-013
dElementMass(58) = 0.0002209
dElementMass(57) = 0.000243908
dElementMass(20) = 0.00027573
dElementMass(6) = 7.4387e-011
dElementMass(92) = 4.9986
dElementMass(46) = 1.3706e-005
dElementMass(26) = 650
dElementMass(37) = 4.2548e-005
dElementMass(44) = 0.00018389
dElementMass(23) = 7.8549e-005
dElementMass(8) = 6.2263e-006
dElementMass(22) = 1.0924e-007
dElementMass(31) = 1.5748e-010
dElementMass(4) = 31.789
dElementMass(5) = 1.3864e-014
dElementMass(25) = 1.4292e-005
dElementMass(34) = 6.27e-006
dElementMass(90) = 0.000376
dElementMass(42) = 30
dElementMass(27) = 30
dElementMass(52) = 3.5248e-005
dElementMass(89) = 6.2653e-017
dElementMass(2) = 6.6179e-005
dElementMass(3) = 63.33
dElementMass(1) = 5.825e-005
dElementMass(10) = 4.3996e-006
dElementMass(28) = 120
dElementMass(43) = 0.00011587
dElementMass(7) = 7.975e-008
dElementMass(33) = 2.088e-008
dElementMass(9) = 146.67
dElementMass(47) = 4.5605e-007
dElementMass(55) = 0.00020044
dElementMass(24) = 170
dElementMass(32) = 4.6163e-008
dElementMass(41) = 9.5515e-006
dElementMass(29) = 3.8898e-006
dElementMass(45) = 8.9291e-006
dElementMass(30) = 1.0976e-010
dElementMass(51) = 1.1177e-006
dElementMass(48) = 8.0338e-007
dElementMass(60) = 0.00015501
dElementMass(54) = 0.00026024
dElementMass(36) = 4.5311e-005
dElementMass(49) = 1.2987e-007
dElementMass(50) = 1.9167e-006
dElementMass(21) = 4.7478e-016
dElementMass(53) = 2.63827e-005
iPrintResultsMode = 1
call ParseCSDataFile(cThermoFileName)
if (INFOThermo == 0) call Thermochimica
if (iPrintResultsMode > 0) call PrintResults
if (INFOThermo == 0) call ResetThermoAll
call ThermoDebug
end program TC_input3
