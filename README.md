# Epiphyte Water Balance (EWB) Model
This repository stores the code for Epiphyte Water Balance (EWB) model. Epiphytes are vascular and nonvascular plants that live on trees and have no connection to the forest floor. Epiphyte mats are a composite of living epiphytes and their detritus, including decomposed organic matter and inorganic nutrients, situated within host tree canopies (Nadkarni, 1984). The model conceptualizes the epiphyte mats as a water store inside the host tree canopy. The EWB model is uncalibrated and allows users to explore the water and energy budgets of epiphyte mats within Tropical Montane Cloud Forest (TMCF) canopies. 

![Fig1_EWSF](https://github.com/user-attachments/assets/cd95428a-04f8-4d36-82ca-59663dbecf81)
Figure 1. (a) The epiphyte water balance (EWB) model conceptualizes the epiphyte mats in TMCF canopies as a water store suspended in the canopy which includes vascular and non-vascular vegetation, mosses, and canopy soils. The host tree branch and leaves shown in black are not considered in the EWB model. (b) Photo of an epiphyte mat on a branch in a TMCF near Monteverde, Costa Rica.

# Journal article
Carchipulla-Morales, D., Corbett, H., Vaughan, D., Gotsch, S., Dawson, T., Nadkarni, N., and Lowman, L., 2025. A novel model uncovers the importance of dew deposition for canopy epiphytes in a tropical montane cloud forest.

# Methodology
The EWB model consists of a water mass balance and an energy balance to solve for water storage within the epiphyte mat and its temperature at a given point in time. In the water mass balance, the epiphyte mat is treated as a water bucket in the canopy that can be filled via rainfall and fog, and depleted via evapotranspiration, and host tree water uptake. For the energy balance, change in epiphate mat temperature is determined based on incoming net radiation, sensible and latent heat fluxes, and the epiphyte mat mass and heat capacity. Consult the accompanying journal article for further details on the equations underlying the EWB model.

![Fig2_EMB](https://github.com/user-attachments/assets/2ca5fdcc-df83-438a-bfc0-60bd876621e3)
Figure 2. Input and output terms for the water mass balance of the EWB model in the canopy that influence epiphyte mat water storage (Se). The EWB model conceptualizes all epiphytic components of the canopy into a single water store. Water enters the epiphyte mat via rainfall (R) and fog (F) interception, and dew deposition (DDe). The amount of water intercepted by R and F depends on the fraction of epiphytes in the canopy (fe) and their respective interception coefficients (cR and cF ). Water leaves the epiphyte mat through evapotranspiration (ETe) and water uptake by the host tree roots (WU).

![Fig3_EEB](https://github.com/user-attachments/assets/b8466eaa-a67e-46bf-b510-6991fe0dac21)
Figure 3. Input and output terms of the energy balance for the EWB model. The epiphyte mat warms by absorbing shortwave radiation (SW) from the above canopy, incoming solar radiation (Φ0), long wave radiation emitted from the atmosphere (LWa) and soil (LWs), latent heat (L) when dew deposition (DDe) occurs, and sensible heat (H) when the air is warmer than the epiphyte mat. The epiphyte mat cools when energy leaves through the reflection of shortwave radiation (αSW), long wave radiation emitted from the epiphyte mat (LWe), latent heat when evapotranspiration occurs (ETe) and sensible heat when the epiphyte mat temperature (Te) is warmer than the air (Ta). In addition to the latent heat term, coupling of the water and energy budgets occurs through the mass (me) and heat capacity of the epiphyte mat (cp,e), which both depend on epiphyte mat water storage (Se).

# Software specifications
This code is available as a MATLAB .m file and provided as a function. This code replicates the ideal simulation presented in the Carchipulla-Morales et al. (in review) manuscript. The function has been tested using MATLAB software versions from 2020a to 2023a. Users can modify the EWB model parameters by editing the arguments in the EWB MATLAB function version 1.0. file. Parameters are presented and described in the [parameters](Parameters.txt) file of this repository. An [user guide](https://github.com/DavidCarMor/EWB/User_guide) for the model is available in this repository.

Technical support can be provided by the first author, David Carchipulla-Morales (carcpd21@wfu.edu).

# Input Variables
The EWB model input variables are provided as a MATLAB table object with the following variables in time steps of one hour. 

| Variable | Description |
| ------------- | ------------- |
| Ta_K | Air temperature in K |
| RH | Relative humidity |
| RF_mm | Rainfall in millimeters |
| F_mm | Fog in milimeters |
| SW_Wpms | Shortwave radiation in Watts per meter squared | 

If rainfall or fog were not measured, the input data corresponds to clear-sky conditions. In this case, the user is encouraged to set all elements in RF_mm and F_mm equal to zero.

# Optional Input Variables
Two optional variables can be added to the input table for the EWB model. These variables may improve the estimation of the aerodynamic boundary layer conductane used in simulating evapotranspiration and the host tree transpiration. If these values are not included in the input variable table, the EWB model treats them as constant parameters and creates a constant vector in the input table according to the parameters values listed on [parameters](Parameters.txt).

| Variable | Description |
| ------------- | ------------- |
| AP_Pa | Atmospheric pressure in pascals|
| WS_mps | Windspeed in meters per second |

# Output Variables
The function returns three objects: an updated version of the input table, a cell array of the parameters in the model, and a table with information on the stomatal conductance of the vascular component of the epiphyte mat. Consult the [parameters](Parameters.txt) file for futher details on the parameters in the EWB model. Information on the Jarvis function used to describe the stomatal conductance of the vascular component of the epiphyte mat can be found in Lowman and Dil Godoy (2020).

| Variable | Description |
| ------------- | ------------- |
| c_pd | Specific heat capacity of the epiphyte mat in W/s/kg/K |
| clf | Estimation of cloud cover - | 
| cos_z | Cosine of zenith angle - |
| DDe | Dew deposition on the epiphyte mat in mm/h |
| EAVD | Epiphyte-atmosphere vapor pressure deficit in Pa |
| e_a | Water vapor pressure of air in Pa |
| em_a | Air emissivity in W/m2/W/m2 |
| es_a | Water vapor pressure of saturated air in Pa | 
| es_epi | Water vapor pressure of saturated epiphyte mat in Pa | 
| ETe | Evapotranspiration from epiphyte mats mm/h |
| ETht | Evapotranspiration from host tree mm/h |
| g_ba | Canopy boundary layer in m/s |
| g_e | Epiphyte mat water conductance in m/s |
| g_sve | Stomatal conductance of vascular epiphytes in m/s |
| H | Epiphyte mat sensible heat in W/m2 |
| L | Epiphyte mat latent heat in W/m2 |
| Ie | Interception of rainfall and fog in mm |
| phi_net | Net radiation received by the epiphyte mat in W/m2 |
| psi_e | Epiphyte mat water potential in m |
| Se | Epiphyte mat water content in mm |
| SF | Stemflow in mm/h |
| sm | Epiphyte mat moisture in kg/kg |
| Te_C | Epiphyte mat temperature in C |
| Te_K | Epiphyte mat temperature in K |
| TF | Throughfall in mm/h |
| theta_e | Epiphyte mat volumetric water content in m3/m3 | 
| VPD | Air Vapor Pressure Deficit in Pa |
| WUht | Host tree water uptake in mm/h |

# References
Lowman, L.E. and Godoy, L.D., 2020. Simulating stomatal response to cloud immersion for montane cloud forests in the Southern Appalachians. Agricultural and Forest Meteorology, 295, p.108165.

Nadkarni, N.M., 1984. Epiphyte biomass and nutrient capital of a neotropical elfin forest. Biotropica 16, 249–256. doi:132.174.248.95.
