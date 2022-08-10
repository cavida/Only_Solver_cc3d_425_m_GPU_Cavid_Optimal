# CompuCell3D-related parameters
mediumStartCoordinate = 72 #21 #72

# Universal constants
universalGasConstant = 8.314 
bodyTemperature = 310.15 # In Kelvin
conversionPixToMeter = 1e-05 #5e-06
conversionPascalToMmHg = 133.32

# Continuum-model values
ecmDegradationByMde = 0.1#1.0#0.8 #1.2 # In molecules/MCS
initECMConcentration = 0.25 # from Hoang Ta's work --> 0.5 for volume 2 cubic micrometer. Here, 1 cubic micrometer
initCytokineConcentration = 5e-23 # Calculation in documentation
maxCytokineConcentration = 4.44e-21 # Calculation in documentation
initNumberOfCadherins = 6594 # Calculation in documentation
cadherinThreshold = 990 # Calculation in documentation
initNumberOfIntegrins = 1000
initNumberOfLigands = 1000

# Cell activity values
minimumO2PressureCancer = 5 # mmHg - lower limit for physioxia in tumor cells
minimumO2PressureEndothelial = 15 # mmHg - 15-38 mmHg is the range for physioxia in endothelial cells
maximumCO2PressureCancer = maximumCO2PressureEndothelial = 40 # mmHg - normal CO2 level
activatedCellTimeThreshold = 2880 #MCS. Equal to two days
cellCycleTime = 1200 #MCS. Equal to 20 hours
                
def defineGlobalValues():
    global mediumStartCoordinate
    global universalGasConstant
    global bodyTemperature
    global conversionPixToMeter
    global conversionPascalToMmHg
    global ecmDegradationByMde
    global initECMConcentration 
    global initCytokineConcentration
    global initNumberOfCadherins
    global cadherinThreshold
    global minimumO2PressureCancer
    global minimumO2PressureEndothelial
    global maximumCO2PressureCancer
    global maximumCO2PressureEndothelial
    global cellCycleTime
    global activatedCellTimeThreshold
    global initNumberOfIntegrins
    global initNumberOfLigands