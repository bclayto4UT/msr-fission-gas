program TC_input2
USE ModuleThermoIO
implicit none
cInputUnitTemperature = 'K'
cInputUnitPressure = 'atm'
cInputUnitMass = 'moles'
cThermoFileName = DATA_DIRECTORY // 'MSTDB-TC_V2.0_Fluorides_8-0.dat'
dPressure = 1
dTemperature = 900
dElementMass(94) = 0.000164752
dElementMass(70) = 2.4542e-013
dElementMass(58) = 0.00014628
dElementMass(57) = 0.000162689
dElementMass(20) = 0.000192172
dElementMass(6) = 4.9575e-011
dElementMass(92) = 4.999
dElementMass(46) = 8.7998e-006
dElementMass(26) = 650
dElementMass(37) = 2.8326e-005
dElementMass(44) = 0.00012422
dElementMass(23) = 3.7867e-005
dElementMass(8) = 4.1509e-006
dElementMass(22) = 7.2822e-008
dElementMass(31) = 1.1364e-010
dElementMass(4) = 31.789
dElementMass(5) = 6.1671e-015
dElementMass(25) = 9.14e-006
dElementMass(34) = 4.1795e-006
dElementMass(90) = 0.00025072
dElementMass(42) = 30
dElementMass(27) = 30
dElementMass(52) = 2.6262e-005
dElementMass(89) = 2.6497e-017
dElementMass(2) = 4.4109e-005
dElementMass(3) = 63.33
dElementMass(1) = 3.8823e-005
dElementMass(10) = 2.933e-006
dElementMass(28) = 120
dElementMass(43) = 7.1654e-005
dElementMass(7) = 5.3175e-008
dElementMass(33) = 1.6592e-008
dElementMass(9) = 146.67
dElementMass(47) = 3.1635e-007
dElementMass(55) = 0.00012681
dElementMass(24) = 170
dElementMass(32) = 3.1747e-008
dElementMass(41) = 4.7928e-006
dElementMass(29) = 2.5781e-006
dElementMass(45) = 4.6203e-006
dElementMass(30) = 8.0648e-011
dElementMass(51) = 8.7284e-007
dElementMass(48) = 5.2122e-007
dElementMass(60) = 9.8906e-005
dElementMass(54) = 0.00017343
dElementMass(36) = 3.0343e-005
dElementMass(49) = 8.1544e-008
dElementMass(50) = 1.2975e-006
dElementMass(21) = 2.2625e-016
dElementMass(53) = 2.16611e-005
iPrintResultsMode = 1
call ParseCSDataFile(cThermoFileName)
if (INFOThermo == 0) call Thermochimica
if (iPrintResultsMode > 0) call PrintResults
if (INFOThermo == 0) call ResetThermoAll
call ThermoDebug
end program TC_input2
