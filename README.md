# EWB toy model
The current repository stores the code of the toy model of energy and water balance for epiphytes in tropical montane cloud forests' canopies. It uses a numerical approach with hourly timesteps. Some parameters used in the model should be adjusted for the area of interest.

# Software specifications
This code was built in MATLAB for windows, and tested from V2020a to V2022a. It will require Symbolic Math Toolbox to be installed.

# Inputs
The toy model has a built-in function to simulate a three clear-sky days in a TMCFs. However, users are encouraged to perform the simulation with their own data sets. So, the code can be run in the following way:
  - Using default data set: EPI_WBV7()
  - Using own data set: EPI_WBV7(input)<br/><br/>
The input should be a table in the workspace containing the following variables:
  - Date-Time (dt): This variable has to be in a date-time format compatible with the datetime() function.
  - Air temperature record (Ta_C): It has to be measured close to the canopy, but not in touch with the epiphytes.
  - Specific humidity (qa) or Relative Humidity (RH) record: The former has to be measured in kg/kg, and the latter in %. The relative humidity is automatically turned into specific humidity using the saturated vapor pressure at current air temperature and the mass of vapor per mass of dry air.
  - Precipitation record (RF): It may be measured with precipitation gauges in mm per hour.
  - Occult precipitation (OP): It may be measured with fog gauges in mm per hour.
  - Shortwave radiation (SW): It must be measured above canopy in Watts per square meter.<br/><br/>
The required name for the variable in the table is established in parenthesis.

# Outputs
The outputs are saved in a folder automatically created for the date and time of the running. This folder contains 2 csv and 18 PNG. The data.csv file contains all variables required to perform the model. The datags.csv contains information about stomatal conductance aperture. The following table specifies the variables in th data.csv file:

| Variable | Variable name in the csv |
| ------------- | ------------- |
| Short wave radiation reflected by TMCF soil  | SW_soil |
| Cosine of zenith angle  | cos_z |
| Net solar radiation  | PHI_net |
| Epiphyte mat volumetric water content  | THETAe |
| Water level in the epiphyte mat  | Se |
| Precipitation interception  | Irf |
| Occult precipitation interception   | Iop |
| Throughfall  | TF |
| Stemflow  | SF |
| Host tree water uptake  | WUht |
| Epiphyte evapotranspiration | ET |
| Host tree epiphyte evapotranspiration | ET |
| Dew deposition | DD |
| Host tree dew deposition | DD_ht |
| Epiphyte mat temperature in Celsius degrees | Tepi_C |
| Epiphyte mat temperature in Kelvin degrees | Tepi_K |
| Relative humidity | RH |
| Air saturated vapor pressure | es_a |
| Epiphyte saturated vapor pressure| es_epi |
| Air vapor pressure| e_a |
| Epiphyte vapor pressure| e_epi |
| Epiphyte specific humidity at saturation| qsat_epi |
| Epiphyte - Air specific humidity difference | D_a |
| Epiphyte - Air saturated vapor pressure diference| EASVD |
| Epiphyte - Air vapor pressure defficit| EAVD |
| Epiphyte water potential| PSIe |
| Soil layer 1 water potential| PSI1 |
| Soil layer 2 water potential| PSI2 |
| Epiphyte mat moisture | sm |
| Epiphyte mat latent heat | L |
| Epiphyte mat sensible heat | H |
| Epiphyte mat heat capacity | cp |
| Epiphyte mat stomatal conductance| gse |
