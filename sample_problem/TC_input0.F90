program TC_input0
USE ModuleThermoIO
implicit none
cInputUnitTemperature = 'K'
cInputUnitPressure = 'atm'
cInputUnitMass = 'moles'
cThermoFileName = DATA_DIRECTORY // 'MSTDB-TC_V2.0_Fluorides_8-0.dat'
dPressure = 1
dTemperature = 900
dElementMass(92) = 5
dElementMass(26) = 650
dElementMass(4) = 31.789
dElementMass(42) = 30
dElementMass(27) = 30
dElementMass(3) = 63.33
dElementMass(28) = 120
dElementMass(9) = 146.67
dElementMass(24) = 170
iPrintResultsMode = 1
call ParseCSDataFile(cThermoFileName)
if (INFOThermo == 0) call Thermochimica
if (iPrintResultsMode > 0) call PrintResults
if (INFOThermo == 0) call ResetThermoAll
call ThermoDebug
end program TC_input0
