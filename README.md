# Epiphyte Mat Energy and Water Mass Balance Model
The current repository stores the code for epiphyte mats' energy and water mass balance model. We consider epiphytes in the host tree as water storage inside the canopy that is filled via rainfall and fog, and depleted via evapotranspiration, and host tree water uptake. This uncalibrated model depicts the interactions of epiphyte mats in the hydrology of Tropical Montane Cloud Forest (TMCF) canopies. 
![Fig1_EWSF](https://github.com/user-attachments/assets/917ad192-96f3-42bf-a8fe-7df3b973b1ec)

# Journal article
Carchipulla-Morales, D., Corbett, D., Vaughan, D., Gotsch, D., Dawson, T., Nadkarni, N., and Lowman, L., 2024. A novel model to simulate water and energy budgets for epiphytic mats.

# Methodology
The water mass balance of the epiphyte mat is solved for a water bucket in the canopy, i.e., the epiphyte mat. In the TMCF canopies, the epiphyte mat is filled via rainfall and fog, and depleted via evapotranspiration, and host tree water uptake. Rainfall and fog interception are approximated using as a fraction of bulk water input times the fractional grid canopy area covered by epiphytes. Evapotranspiration and its reciprocal, Dew Deposition, are estimated using the Penman-Monteith equation. Finally, the host tree water uptake is a fraction of the host tree transpiration.

![Fig2_EMB](https://github.com/user-attachments/assets/495d0695-ed9e-46da-b05c-8755802f6932)

Similarly, the energy balance of the epiphyte mat is solved for the fractional grid canopy area covered by epiphytes and its interaction with the rest of the canopy. The net radiation accounts for the effective shortwave radiation that reaches the epiphyte mat, long wave radiation from the air surrounding the epiphyte mat, and the longwave radiation that the epiphyte mat emits. The sensible heat considers the temperature difference between the epiphyte mat and air. Finally, latent heat reconciles the energy required to evaporate or condense water in the epiphyte mat. 

![Fig3_EEB](https://github.com/user-attachments/assets/84cdb0be-baa5-4806-bc50-2f2973e0e823)

# Software specifications
This code was built in MATLAB and shared as a function. The function has been tested from MATLAB V2020a through MATLAB V2022a. Parameters of the model can be modified and their information can be found in [parameters](Parameters.txt). A tutorial for the model was added to this repository.

Technical support can be provided by the corresponding author Dr. Lauren E. L. Lowman at lowmanle@wfu.edu

# Input Variables
The input for the function is a MATLAB table object with the following variables in time steps of one hour. 

| Variable | Description |
| ------------- | ------------- |
| Ta_K | Air temperature in K |
| RH | Relative humidty |
| RF_mm | Rainfall in milimeters |
| F_mm | Fog in milimeters |
| SW_Wpms | Shortwave radiation in Watts per meter squared | 

If rainfall or fog were not measured, the data corresponds to clear-sky conditions. The user is encouraged to set all elements in RF_mm and F_mm equal to zero.

# Optional Input Variables
A couple of optional variables can be added to the input table for the model. They will be used to improve the estimations of the conductance of the aerodynamic boundary layer and the host tree transpiration. If these values are absent, the model creates a constant vector to add as variables in the input table according to the parameters listed on [parameters](Parameters.txt).

| Variable | Description |
| ------------- | ------------- |
| AP_Pa | Atmospheric pressure in pascals|
| WS_mps | Wwindspeed in meters per second |

# Output Variables
The function returns three objects: an updated version of the input table, a cell array of the parameters in the model, and a table with information on the stomatal conductance of the vascular component of the epiphyte mat. Go to [parameters](Parameters.txt) for details of the parameters in the model. Information on the Jarvis function used to describe the stomatal conductance of the vascular component of the epiphyte mat can be found in Lowman, Lauren EL, and Luis Dil Godoy. "Simulating stomatal response to cloud immersion for montane cloud forests in the Southern Appalachians." Agricultural and Forest Meteorology 295 (2020): 108165.

| Variable | Description |
| ------------- | ------------- |
| c_pd | Specific heat capacity of the epiphyte mat in W/s/kg/K |
| clf | Estimation of cloud cover - | 
| cos_z | Cosine of zenith angle - |
| DDe | Dew deposition on the epiphyte mat in mm/h |
| EAVD | Epiphyte-atmosphere vapor pressure deficit in Pa |
| e_a | Water vapor pressure of air in Pa |
| e_epi | Water vapor pressure of the epiphyte mat in Pa |
| em_a | Air emissivity in W/m2/W/m2 |
| es_a | Water vapor pressure of saturated air in Pa | 
| es_epi | Water vapor pressure of saturated epiphyte mat in Pa | 
| ETe | Evapotranspiration from epiphyte mats mm/h |
| ETht | Evapotranspiration from host tree mm/h |
| g_ba | Canopy boundary later in m/s |
| g_e | Epiphyte mat water conductance in m/s |
| g_sve | Stomatal conductance of vascular epiphytes in m/s |
| H | Epiphyte mat sensible heat in W/m2 |
| L | Epiphyte mat latent heat in W/m2 |
| Ie | Interception of rainfall and fog in mm |
| phi_net | Net radiation received by the epiphyte mat in W/m2 |
| psi_e | Epiphyte mat water potential in m |
| Se | Epiphyte mat water content in mm/h |
| SF | Stemflow in mm/h |
| sm | Epiphyte mat moisture in kg/kg |
| Te_C | Epiphyte mat temperature in C |
| Te_K | Epiphyte mat temperature in K |
| TF | Throughfall in mm/h |
| theta_e | Epiphyte mat volumectric water content in m3/m3 | 
| VPD | Air Vapor Pressure Deficit in Pa |
| WUht | Host tree water uptake in mm/h |
