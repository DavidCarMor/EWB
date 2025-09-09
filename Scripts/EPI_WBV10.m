function [data,opt,datags]=EPI_WBV10(data,opt)
%% Function details
% This function computes the energy and water mass balance of epiphytes in
% a TMCF canopy, where epiphytes are considered an uniform water storage
% according to Carchipulla-Morales, et al., YYYY.
% The current version of the model works with hourly weather records.
% Modifications can be applied to work at finer time resolutions or time
% resolutions smaller than 12 hours. Coarser time resolutions might
% recquire a modification of the assumptions.
% 
% Input variables:
%   -data: table with variables Date air temperature in Kelvin (Ta_K), 
%   relative humidity as fraction (RH), rainfall or vertical precipitation 
%   in mm (RF_mm), fog or horizontal precipitation in mm (F_mm), solar 
%   shortwave radiation in watts per meter squared(SW_Wpms).
%   Table data can also have atmospheric pressure in pascals (AP_Pa) and
%   windspeed in meters per second (WS_mps). If AP_Pa is not specified, the
%   code will add the variable AP_Pa to data with a fixed value according
%   to parameter AP_Pa (read list of parameters below to see more
%   information). IF WS_mps is specifiec, a variable aerodynamic conductance
%   of the canopy boundary layer is estimated. Otherwise, a fixed aerodynamic
%   conductance of the canopy boundary layer is defined according to
%   parameters ga and gb.
% 
% Input parameters:
%   - alph: Albedo of the canopy. The default value is 0.15,
%   - AP_Pa: fixed atmospheric pressure in Pascals used to create the variable AP_Pa in table data. The default value is 88000 Pa,
%   - b1: Parameter for stomatal sensitivity to incoming solar radiation.
%   The default value is 0.012 m2*s/umol with the following conversion factor
%   1 W/m2 = 2.02 umol/m2/s (dos Reis, Mariana Gonçalves, and Aristides Ribeiro, 2020),
%   - b: Texture coefficient of water potential curve. The default value is 7.75,
%   - fe: Fractional grid canopy area covered by epiphytes. The default value is 0.3 m2/m2,
%   - c_F: Fraction of fog captured in the canopy. The default value is 0.5 mm/mm,
%   - c_p: Specific heat capacity of air. The default value is 1010 W s/kg/K,
%   - c_pw: Specific heat capacity of water. The default value is 4182 W s/kg/K,
%   - c_pe: Specific heat capacity of dry epiphyte mat. The default value is 1007.75 W s/kg/K,
%   - c_RF: Fraction of rainfall captured in the canopy. The default value is 0.5 mm/mmm,
%   - day_l: Day length. The default value is 12 h,
%   - Dx: Temperature-dependent parameter for stomatal sensitivity to EAVD. The default value is 1250 Pa,
%   - epsiln_a: Night air emissivity. The default value is 0.8 W/m2/W/m2,
%   - epsiln_v: Emissivity of vegetation. The default value is 0.95 W/m2/W/m2,
%   - gamma_w: Psychrometric constant. The default value is 66.1 Pa/K,
%   - g_a: Aerodynamic conductance of canopy boundary layer. The default value is 0.02 m/s,
%   - g_b: Additional aerodynamic conductance of canopy boundary layer. The default value is 0.02 m/s,
%   - gsve_max: Max. vascular epiphyte stomatal conductance. The default value is 0.000000482654 m/s,
%   - gnve: Non-vascular epiphyte water conductance. The default value is 0.01 m/s,
%   - gsht_max: Max. host tree stomatal conductance. The default value is 0.0124 m/s,
%   - h_c: Canopy height. The default value is 30 m,
%   - h_f1: Zero displacement factor. The default value is 0.7 m/m,
%   - h_f2: Heat and vapor - momentum transfer ratio 0.2,
%   - k: von Karman's constant. The default value is 0.41,
%   - k2: Parameter for stomatal sensitivity to temperature stress. The default value is 0.0016 1/K2,
%   - k_b: Beer's law extinction coefficient. The default value is 0.5,
%   - latitude: Geographical latitude of the epiphyte mat. The default value is 10.3°,
%   - LAI: Leaf area index of the forest. The default value is 7.7 m2/m2,
%   - lambda_w: Latent heat of vaporisation of water. The default value is 2.50e6 J/kg,
%   - m_d: Dry epiphyte mat biomass in forest canopy. The default value is 49095 kg/ha,
%   - VEDB: Vascular epiphyte biomass in forest canopy. The default value is 28110 kg/ha,
%   - NVEDB: Non-vascular epiphyte biomass in forest canopy. The default value is 5331 kg/ha,
%   - porosity: Epiphyte mat porosity. Default value is 0.6 m3/m3,
%   - psi_ae: Air-Entry water potential of epiphyte mat (silty clay soil). The default value is -0.356 m,
%   - psi_fc: Water potential at epiphyte mat field capacity. The default value is -16.4113 m,
%   - psi_wp: Water potential at epiphyte mat wilting point. The default value is -1529.6164 m,
%   - rho: Air density. The default value is 1.2 kg/m3,
%   - rho_s: Soil organic matter density. The default value is 1300 kg/m3,
%   - rho_w: Density of water. The default value is 998 kg/m3,
%   - s: Pressure-temperature curve slope. The default value is 189 Pa/K,
%   - So: Initial water content. The default value is NaN mm,
%   - sigma: Stefan - Boltzmann constant. Default value is 5.67e-8 W/m2/K4,
%   - SWclr: Clear Sky Shortwave radiation. The default value is 1100 W/m2,
%   - Te_Ta: Initial guess for the difference between Te and Ta. The default value is -0.5 K,
%   - tau: Fraction of transpiration of the host tree. The default value is 0.15,
%   - Topt: Optimum air temperature for stomata conductance. The default value is 293.15 K,
%   - t0: Sunrise time. The default value is 7 am LST, and
%   - wf: Water-Biomass weight ratio. Default value is 4.
% 
% Output parameters:
%   - h: Epiphyte mat depth in m, = (opt.Se_max/1000)/opt.porosity,
%   - LAIe: Epiphyte mat leaf area index in m2/m2,
%   - mcount: Number of times we used the first day until the model stabilized,
%   - psi_0: Water potential when stomata are completely closed in m
%   - psi_1: Water potential when stomata are completely open in m
%   - rho_s: Epiphyte mat density in kg/m3,
%   - Se_max: Maximum water level of epiphyte mat in mm,
%   - SWclr: if the weather data has a shortwave radiation larger than the default SWclr, the code defines the max shortwave radiation as the new
%   Clear Sky Shortwave radiation in W/m2,
%   - theta_fc: Volumetric water content of the epiphyte mat at field capacity in m3/m3, and
%   - theta_wp: Volumetric water content of the epiphyte mat at wilting point in m3/m3
% 
% Output variables:
%   - c_pe: Specific heat capacity of the epiphyte mat in W/s/kg/K,
%   - clf: Estimation of cloud cover -, 
%   - cos_z: Cosine of zenith angle -,
%   - DDe: Dew deposition on the epiphyte mat in mm/h,
%   - EAVD: Epiphyte-atmosphere vapor pressure deficit in Pa,
%   - e_a: Water vapor pressure of air in Pa,
%   - em_a: Air emissivity in W/m2/W/m2,
%   - es_a: Water vapor pressure of saturated air in Pa, 
%   - es_epi: Water vapor pressure of saturated epiphyte mat in Pa, 
%   - ETe: Evapotranspiration from epiphyte mats mm/h,
%   - ETht: Evapotranspiration from host tree mm/h,
%   - g_ba: Canopy boundary later in m/s,
%   - g_e: Epiphyte mat water conductance in m/s,
%   - g_sve: Stomatal conductance of vascular epiphytes in m/s,
%   - H: Epiphyte mat sensible heat in W/m2,
%   - L: Epiphyte mat latent heat in W/m2,
%   - Ie: Interception of rainfall and fog in mm,
%   - phi_net: Net radiation received by the epiphyte mat in W/m2,
%   - psi_e: Epiphyte mat water potential in m,
%   - Se: Epiphyte mat water content in mm/h,
%   - SF: Stemflow in mm/h,
%   - sm: Epiphyte mat moisture in kg/kg,
%   - Te_C: Epiphyte mat temperature in °C,
%   - Te_K: Epiphyte mat temperature in K,
%   - TF: Throughfall in mm/h,
%   - theta_e: Epiphyte mat volumectric water content in m3/m3, 
%   - VPD: Air Vapor Pressure Deficit in Pa, and
%   - WUht: Host tree water uptake in mm/h.
% 
% Notes:
% 
%  If you want to run the funcion, make sure you are in the folder with all
%  the functions required to run the model. Otherwise add the path, e.g.:
% 
%           addpath(genpath(folder containing the this code));
% 
% winput() is supporting code that generates the recquired variables to run
% the model. 
% 
% An example for running the code is:
%           data=EPI_WBV9(data=winput(2),LAI=7,ce=0.25);
% 
% D. Carchipulla-Morales, H. Corbett, and L.Lowman (February, 2022)
% Last update: April, 2025

%% Inputs
arguments % Defaults
    data table % Weather data
    % opt.alph double = 0.35 % Vegetation albedo [W m-2 / W m-2]
    opt.alph double = 0.15 % Vegetation albedo [W m-2 / W m-2]
    opt.AP_Pa double = 88000 % [Pa] https://www.worldmeteo.info/es/america-central/costa-rica/monteverde/tiempo-136622/
    opt.b1 double = 0.012*2.02 % Parameter for stomatal sensitivity to incoming solar radiation [m2 W-1]
    opt.b double = 7.75 % Exponent of Soil Mositure in soil-water retention curve (Dingman, pg 338)
    opt.fe double = 0.3 % Fractional grid canopy area covered by epiphytes [m2 m-2]
    opt.c_F double = 0.5 % Fraction of fog captured in the canopy [mm/mm]
    opt.c_p double = 1010 % Specific capacity of the air [W s-1 kg-1 K-1]
    opt.c_pw double = 4182 %Specific capacity of the water [W s-1 kg-1 K-1]
    opt.c_pd double = 3433.409 % Specific capacity of the epiphyte mat [W s kg-1 K-1]
    opt.c_RF double = 0.5 % Fraction of rainfall captured in the canopy [mm/mmm]
    opt.day_l double = 12 % Day length [h]
    opt.Dx double = 1250 % Temperature dependent parameter for stomatal sensitivity to EAVD [Pa]
    opt.epsiln_a double =0.8 % Night air emissivity 
    opt.epsiln_v double = 0.95 % Emissivity of vegetation [W m-2 / W m-2]
    opt.gamma_w double = 66.1 %Psychrmetric constant [Pa K-1]
    opt.g_a double = 20*(1/1000) % Aerodynamic conductance [m s-1]
    opt.g_b double = 20*(1/1000) % Additional boundary layer conductance [m s-1]
    % opt.gsve_max double = 2.48*(1/1000) % Max. vascular epiphyte water exchange conductance [m s-1]
    opt.gsve_max double = 4.82654*(1/10000000) % Max. vascular epiphyte water exchange conductance [m s-1]
    opt.gnve double = 10*(1/1000) % Nnon-vascular epiphyte water exchange conductance [m s-1]
    opt.gsht_max double = 0.0124 % Max. host tree stomatal conductance [m s-1]
    opt.h_c double = 30 % Canopy height [m]
    opt.h_f1 double = 0.7 % Zero displacement factor [m/m]
    opt.h_f2 double = 0.2 % heat and vapour - momentum transfer ratio
    opt.k double = 0.41 % von Karman's constant
    opt.k2 double = 0.0016 % coefficient for temperature stress [K-2]
    opt.k_b double = 0.5 % Beer's law extintion coefficient
    opt.latitude double = 10.3 % Latitude of epiphyte mat [°]
    opt.LAI double = 7.7 % Leaf area index of the forest [m2 m-2]
    opt.lambda_w double = 2.50e6 % Latent heat of vaporisation of water [J kg-1]
    opt.m_d double = 49095 % Dry epiphyte mat biomass in forest [kg ha-1]
    % opt.VEDB double = 69434.57 % Vascular epiphyte dry biomass [kg ha-1]
    % opt.NVEDB double = 10585.79 % Non-vascular epiphyte dry biomass [kg ha-1]
    opt.VEDB double = 28110 % Vascular epiphyte dry biomass [kg ha-1]
    opt.NVEDB double = 5331 % Non-vascular epiphyte dry biomass [kg ha-1]
    opt.porosity double = 0.6 % Epiphyte mat porosity [m3 m-3]
    opt.psi_ae double = -0.356 % Air-Entry water potential of epiphyte mat [m]
    % opt.psi_fc double = -132.2946 % Water potential at field capacity [m] (Water is not moving anymore)
    opt.psi_fc double = -16.4113 % Water potential at field capacity [m] (Water is not moving anymore)
    opt.psi_wp double = -1529.6164 % Wilting point of plants [m]
    opt.rho double = 1.2 % Air density [kg m-3]
    opt.rho_s double = 1300 % Soil organic matter density [kg m-3]
    opt.rho_w double = 998 % Water density [kg m-3]
    opt.s double = 189 % Pressure-temperature curve slope [Pa K-1]
    opt.So double = NaN % Initial water content [mm]
    opt.sigma double = 5.67e-8 % Stefan - Boltzmann constat [W m-2 K-4]
    opt.SWclr double = 1100 % Clear Sky Shortwave radiation [W m-2]
    opt.Te0 double = 293.15 % Initial guess for Te [K]
    opt.tau double = 0.15 % Fraction of transpiration of the host tree
    opt.Topt double = 293.15 % Optimum air temperature [K]
    opt.t0 double = 7 % Sunrise time
    opt.wf double = 4 % Water-Biomass weight ratio [kg of water m-2 / kg of biomass m-2]
end

%% Vectors to store data
data.c_pe = NaN(size(data,1),1); % Specific heat capacity of the epiphyte mat in [W/s/kg/K]
data.clf = NaN(size(data,1),1); % Estimation of cloud cover 
data.DDe = NaN(size(data,1),1); % Dew deposition on the epiphyte mat [mm/h]
data.EAVD = NaN(size(data,1),1); % Epiphyte-atmosphere vapor pressure deficit [Pa]
data.e_a = NaN(size(data,1),1); % Water vapor pressure of air [Pa]
data.epsiln_a = NaN(size(data,1),1); % Air emissivity
data.es_a = NaN(size(data,1),1); % Water vapor pressure of saturated air  [Pa]
data.es_epi = NaN(size(data,1),1); % Water vapor pressure of saturated epiphyte mat [Pa]
data.ETe = NaN(size(data,1),1); % Evapotranspiration from epiphyte mats  [mm h-1]
data.ETht = NaN(size(data,1),1); % Evapotranspiration from host tree [mm h-1]
data.g_ba = zeros(size(data,1),1); % Canopy boundary later  [m s-1]
data.g_e = NaN(size(data,1),1); % Epiphyte mat water conductance  [m s-1]
data.g_sve = NaN(size(data,1),1); % Stomatal conductance of vascular epiphytes [m s-1]
data.H = NaN(size(data,1),1); % Sensible heat [W m-2]
data.L = NaN(size(data,1),1); % Latent heat [W m-2]
data.Ie = NaN(size(data,1),1); % Interception of rainfall and fog [mm]
data.phi_net = NaN(size(data,1),1); % Net radiation received by the epiphyte mat  [W m-2]
data.psi_e = NaN(size(data,1),1); % Epiphyte mat water potential [m]
data.Se = NaN(size(data,1),1); % Epiphyte mat water content [mm h-1]
data.SF = NaN(size(data,1),1); % Stemflow [mm h-1]
data.sm = NaN(size(data,1),1); % Epiphyte mat moisture [kg kg-1]
data.Te_C = NaN(size(data,1),1); % Epiphyte mat temperature [C]
data.Te_K = NaN(size(data,1),1); % Epiphyte mat temperature [K]
data.TF = NaN(size(data,1),1); % Throughfall[mm h-1]
data.theta_e = NaN(size(data,1),1); % Epiphyte mat volumectric water content [m3/m3]
data.VPD = NaN(size(data,1),1); % Air vapor pressure deficit [Pa]
data.WUht = NaN(size(data,1),1); % Host tree water uptake [mm h-1]

%% Unit change from inputs
data = table2timetable(data);
% Air temperature [K]
Ta_C = data.Ta_K-273.15;
data = addvars(data,Ta_C,'After','Ta_K');
clear Ta_C

%% Output parameters
% Epiphyte Leaf Area Index [m2 m-2]
opt.LAIe = opt.LAI*opt.fe;
% Maximum water level of epiphyte mat [mm]
opt.Se_max = opt.m_d*opt.wf/(10*opt.rho_w);
% Epiphyte mat depth [m]
opt.h = (opt.Se_max/1000)/opt.porosity;
% Water potential when stomata are completely open [m]
opt.psi_1 = opt.psi_fc;
% Water potential when stomata are completely closed [m]
opt.psi_0 = opt.psi_wp;
% Making sure SWclr is larger than max radiation
opt.SWclr = max(max(data.SW_Wpms),opt.SWclr);
% Volumetric water content of the epiphyte mat at field capacity [m3/m3]
opt.theta_fc = opt.porosity/((opt.psi_fc/opt.psi_ae)^(1/opt.b));
% Volumetric water content of the epiphyte mat at wilting point [m3/m3]
opt.theta_wp = opt.porosity/((opt.psi_wp/opt.psi_ae)^(1/opt.b));

%% Initial variables in the epiphyte mat
% Initial Water level] [mm]
if isnan(opt.So) % if an initial water level is not provided, the code use max. water level [mm]
    data.Se(1) = opt.Se_max; 
else % if an initial water level is provided,
    if opt.So <= opt.Se_max % the code checks this water level is lower than the max value [mm],
        data.Se(1) = opt.So; % then the initial water level is set as initial value in data.Se(1)
    else
        warning('Initial water level can not be greater than max. water level.')
        fprintf('Max. water level = %f', opt.Se_max)
        return
    end
end
% Initial volumetric water content in the epiphyte mat [m3/m3]
data.theta_e(1)= (data.Se(1)/1000)/opt.h; 

% Initial Epiphyte temperature
data.Te_K(1) = opt.Te0; % [K]
data.Te_C(1) = opt.Te0-273.15; % [C]

%% Stomatal conductance - Variables in Jarvis function
datags = zeros(size(data,1),5);

%% Forcing variables

% Cos(Zenith)
data = cos_zenith(data,opt.latitude,opt.t0,opt.day_l); 

% Aerodynamic conductance + Additional boundary layer conductance [m s-1]
if sum(strcmp('WS',data.Properties.VariableNames)) == 1
    data.g_ba = aerCond(opt.h_c,data.WS_mps,opt.k,opt.h_f1,opt.h_f2); % Aerodynamic conductance
    % Aerodynamic conductance for  wind speed higher than 1.5687 or nan
    data.g_ba(data.WS_mps > 1.5687 | isnan(data.WS_mps)) = (opt.g_a.*opt.g_b)./(opt.g_a+opt.g_b); 
else % Default boundary layer
    data.g_ba = data.g_ba+(opt.g_a.*opt.g_b)./(opt.g_a+opt.g_b); % Aerodynamic conductance [m s-1]
end

% Atmospheric pressure (Check if input table has pressure data)
if sum(strcmp('AP_Pa',data.Properties.VariableNames))~=1
    data.AP_Pa = zeros(size(data,1),1)+opt.AP_Pa; % Fixed Atmospheric Pressure [Pa]
end

% Water vapor pressure of saturated air [Pa] 
% Saturated vapor pressure of the air [Pa]
data.es_a = esat(data.Ta_K);

% Water vapor pressure in air [Pa]
data.e_a = data.RH.*data.es_a;

% Estimation of cloud cover 
data.clf = 1 - data.SW_Wpms./opt.SWclr;

% Air emissivity
data.epsiln_a = data.clf + (1-data.clf).*(1.24.*((data.e_a./(100*data.Ta_K)).^(1/7))); % day
data.epsiln_a(hour(data.Date) <= opt.t0 | hour(data.Date) >= opt.t0+opt.day_l) = opt.epsiln_a; % night

% Epiphyte-to-air vapor pressure deficit [Pa]
data.VPD = max(0,data.es_a - data.e_a); % Epiphyte tank - air

%% Epiphyte tank water balance

for i = 1:size(data,1)
    % Epiphyte mat moisture [kg kg-1]
    data.sm(i)=data.theta_e(i)*(opt.rho_s/opt.rho_s);

    %% Water and Energy balance
    % Atmospheric water interception [mm]
    data.Ie(i) = (opt.c_RF*data.RF_mm(i)+opt.c_F*data.F_mm(i))*opt.fe;
    if data.RF_mm(i) > 0
        % Throughfall + Stemfall [mm]
        % Partitioning literature suggest that stemflow only occurs when
        % epiphyte are saturated
        TFSF = data.RF_mm(i)*(1-opt.c_RF)*opt.fe;
    else
        TFSF = 0;
    end

    % Water vapor pressure of saturated epiphyte mat [Pa]
    data.es_epi(i) = esat(data.Te_K(i));

    % Epiphyte-atmosphere vapor pressure deficit  [Pa]
    data.EAVD(i) = max(0,data.es_epi(i) - data.e_a(i)); % Epiphyte tank - air
    
    % Epiphyte mat water potential [m]
    data.psi_e(i) = ((opt.porosity/data.theta_e(i))^(opt.b))*(opt.psi_ae);

    % Epiphyte net radiation [W m-2]
    % Air emissivity: https://www.acs.org/content/acs/en/climatescience/atmosphericwarming/singlelayermodel.html#:~:text=The%20average%20temperature%20of%20the,atmospheric%20emissivity%20of%20about%200.8.
    PHIe = data.SW_Wpms(i)*exp(-opt.k_b*opt.LAI*(1-opt.fe));
    SW_in_out = PHIe*data.cos_z(i)*(1-opt.alph);
    LW_in_out = opt.sigma*(data.epsiln_a(i)*(data.Ta_K(i)^4) - opt.epsiln_v*(data.Te_K(i)^4));
    data.phi_net(i) = (SW_in_out + LW_in_out)*opt.fe;

    % Stomatal conductance [m s-1]
    [data.g_sve(i),datags(i,:)] = sc_i(opt.gsve_max,opt.b1,SW_in_out,opt.k2,data.Ta_K(i),...
        opt.Topt,data.psi_e(i),opt.psi_0,opt.psi_1,data.EAVD(i),opt.Dx);

    % Weighted average for epiphyte mat conductance
    data.g_e(i) = (data.g_sve(i)*(opt.VEDB/(opt.VEDB+opt.NVEDB)))+...
        (opt.gnve*(opt.NVEDB/(opt.VEDB+opt.NVEDB)));

    % Dew deposition or evapotranspiration [mm h-1]
    data.L(i) = ET_DD(data.VPD(i),data.phi_net(i),data.g_ba(i),data.g_e(i),opt.lambda_w,...
        opt.gamma_w,opt.rho,opt.rho_w,opt.s,opt.c_p,opt.LAIe);

    if data.L(i) > 0 % data.es_epi(i) > data.e_a(i)
        % Evapotranspiration [mm h-1]
        data.ETe(i) = abs(data.L(i));
        data.DDe(i) = 0;
    else
        % Dew deposition [mm h-1]
        data.DDe(i) = abs(data.L(i));
        data.ETe(i) = 0;
    end

    % Latent heat [W m-2]
    data.L(i)= -1*opt.lambda_w*opt.rho_w*(data.L(i)/3600000);

    % Heat conductance [m s-1]
    gah = data.g_ba(i);

    % Sensible heat [W m-2]
    data.H(i)=(data.Ta_K(i)-data.Te_K(i))*opt.c_p*opt.rho*gah;

    % Heat capacity of the composite [W s/kg/K]
    data.c_pe(i)=cpe(data.Se(i),opt.c_pw,opt.c_pd,opt.rho_w,opt.m_d);

    % Energy Balance (Epiphyte temperature [K and °C])
    if i < size(data,1)
        data.Te_K(i+1) = data.Te_K(i) + (3600 * (data.phi_net(i)+data.L(i)+data.H(i))/ ...
            m_cp(data.Se(i),opt.c_pw,opt.c_pd,opt.rho_w,opt.m_d) ); 
        data.Te_C(i+1) = data.Te_K(i+1)-273.15; 

    end

    % Host tree evapotranspiration
    LAVD = data.EAVD(i); % It's the same value because of temperature harmony in the canopy
    data.ETht(i) = opt.rho*opt.gsht_max*0.622*LAVD/data.AP_Pa(i);

    % Host tree water uptake
    data.WUht(i) = opt.tau*data.ETht(i);

    % Water Balance (Water level [mm])
    if i < size(data,1)
        % Next water content
        data.Se(i+1) = data.Se(i)+data.Ie(i)+...
            data.DDe(i)-data.ETe(i)-data.WUht(i);

        % What if it is saturated
        if data.Se(i+1) >= opt.Se_max
            data.SF(i+1) = data.Se(i+1) - opt.Se_max; % Stemflow [mm h-1]
            data.Se(i+1) = opt.Se_max; % Adjusting water Se [mm h-1]
        end
        if TFSF > data.SF(i+1)
            data.TF(i+1) = TFSF - data.SF(i+1); % Throughfall [mm h-1]
        else
            data.TF(i+1) = TFSF; % Throughfall [mm h-1]
        end


        % Volumetric water content [mm/mm]
        data.theta_e(i+1)= (data.Se(i+1)/1000)/opt.h;
    end

    %% Stopping the model when there are data causing errors
    % Stopping model when abs(data.Te_K(i)-data.Te_K(i-1)) >= 10
    % Using a constant boundary layer prevent us to simulate low
    % maximum water contents, because sensible heat exchange is not well
    % regulated without that variable conductance
    if i < size(data,1) && abs(data.Te_K(i)-data.Te_K(i+1)) >= 10 && i > 1
        % data(~any(data.theta_e,2),:) = []; datags(~any(data.theta_e,2),:) = [];
        data(i+1, end-28:end) = {NaN};
        warning(['The simulation stopped because the epiphyte mat temperature ' ...
            'changed more than or equal to 10 K in one hour.'])
        break
    end

    % Stopping model when abs(data.Tepi_K(i)-data.Tepi_K(i-1)) >= 10
    if i < size(data,1) && data.theta_e(i+1) <= 0  && i > 1
        % data(~any(data.theta_e,2),:) = []; datags(~any(data.theta_e,2),:) = [];
        data(i+1, end-28:end) = {NaN};
        warning(['The simulation stopped because the epiphyte water tank is empty '...
            'or water conter is below wilting point.'])
        break
    end
end

%% Exporting data
% Rename variable names of datags before exporting it
datags = array2table(datags);
datags = renamevars(datags,1:size(datags,2),["g_cemax","f_sw","f_Ta","f_psi","f_eavd"]);


end

%% Supporting functions
% Zenith angle of field data
function cos_out = cos_zenith(data,latitude,t0_param,day_l)
    % data = renamevars(data,'Timestamp','Date');
    % Declination of the sun
    d_sun = asind(sind(-23.44)*cosd(((360/365.24)*(day(data.Date,'dayofyear')+10))+...
        ((360/pi)*0.0167*sind((360/365.24)*(day(data.Date,'dayofyear')-2)))));
    
    % Fractional DOY in radians
    dt = dateshift(data.Date, 'start', 'day'); uniq_dt = unique(dt);
    t0 = ones(size(uniq_dt,1),1) * t0_param; % Initial sunlight hour (hours)
    delta = ones(size(uniq_dt,1),1) * day_l; % Day length [h]
    sn = ones(size(uniq_dt,1),1) * median(t0:t0+delta); % Solar noon
    % Map daily sn values back to original data length
    sn_expanded = zeros(size(data.Date));

    for i= 1:size(uniq_dt,1)
        day_mask = dt == uniq_dt(i);
        sn(i) = find(data.SW_Wpms(day_mask) == ...
            max(data.SW_Wpms(day_mask) ), 1, 'first')-1;
        sn_expanded(day_mask) = sn(i);
    end
    
    % Now compute gamma_DOY safely
    gamma_DOY = (2*pi/365) * ((day(data.Date,'dayofyear') - 1) + ((hour(data.Date) - sn_expanded) / 24));
    
    % Equation of time
    EoT = 229.18.*(0.000075+(0.001868.*cos(gamma_DOY))-(0.032077.*sin(gamma_DOY))...
        -(0.014615.*cos(2.*gamma_DOY))-(0.040849.*sin(2.*gamma_DOY)));
    
    % Offset
    Offset = EoT +(4.*(-84.816634-(15.*-6)));
    
    % Hour angle
    h_r = 15.*((hour(data.Date)+(Offset./60))-sn_expanded);
    h_r(data.SW_Wpms==0)=0;
    % h_r = acos(-tan(d_sun).*tan(0.18)); % Simple alternative
    
    % Cosine of Zenith
    data.cos_z = (cosd(latitude).*cosd(h_r).*cosd(d_sun))+(sind(latitude).*sind(d_sun));
    data.cos_z(data.SW_Wpms==0)=0;
    data.cos_z(data.cos_z<0)=0;

    cos_out = data;
    % equations:
    % https://solarsena.com/solar-hour-angle-calculator-formula/
    % https://solarsena.com/solar-declination-angle-calculator/
end

%% Heat capacity update
function ct = cpe(Se,c_pw,c_pd,RHOw,m_d)
% Function to compute the mass of soil and water in EWT
% Inputs
%     Se = Water level in EWT [mm h-1]
% Parameters
%     c_pw % Specific capacity of the water [W s kg-1 K-1]
%     c_pd % Specific capacity of the dry epiphyte mat [W s kg-1 K-1]
%     RHOw % Water density [kg m-3]
%     m_d % Epiphyte dry biomass [kg ha-1] 

m_d = m_d/10000; % Epiphyte dry biomass [kg m-2] 
% Calculations
m_Se = (Se/1000)*RHOw; % Mass of water per surface [kg m-2]66.076
fs = m_d/(m_Se+m_d); % Influence of epiphyte mat in heat capacity
fw = m_Se/(m_Se+m_d); % Influence of water in heat capacity
ct = ((fw*c_pw)+(fs*c_pd)); % Mass times Heat capacity

end

%% Saturated vapor pressure
function es = esat(temp)
% Function to compute saturated vapor pressure in either air or leaf
%
%     input = temperature of either air or leaf [K]
%     output = saturation vapor pressure for temperature [Pa]
%
% Lowman, L. E., & Godoy, L. D. (2020). Simulating stomatal response to
% cloud immersion for montane cloud forests in the Southern Appalachians.
% Agricultural and Forest Meteorology, 295, 108165.

es = 611.71*exp((2.501/0.000461)*((1/273.15)-(1./(temp))));
end

%% Evapotranspiration Penmann-Monteith
function r = ET_DD(VPD, PHI,gba,gs,LAMBDAw,GAMMAw,RHO,RHOw,s,cp,LAIe)
% Function to compute ET or dew deposition. The result is a function of 
% the input difference between specific humidities. Use D_epi for ET and
% D_a for DD.
% 
% Inputs
%     VPD = air vapor pressure deficit [Pa].
%     PHI = Net radiation [W m-2].
%     gba = epiphyte - air boundary layer conductance [m s-1].
%     gs = epiphyte stomatal conductace [m s-1].
% Parameters
%     LAMBDAw % Latent heat of vaporisation of water [W s kg-1]
%     GAMMAw %Psychrmetric constant [Pa K-1]
%     RHO %Air density [kg m-3] 
%     RHOw %Water density [kg m-3] 
%     s % Pressure-temperature curve slope [Pa K-1]
%     cp  % Specific capacity of the air [W s kg-1 K-1]
%     LAIe % Epiphyte leaf area index [m2 m-2]
%     
% Output
%     ET_DD = Evapotranspiration or dew deposition [mm h-1].

% numerator and denominator
num1 = (s*PHI)+(RHO*cp*gba*VPD);
denom1 = RHOw*LAMBDAw*(s+(GAMMAw*(1+(gba/(gs*LAIe)))));

% Result in mm/h
r = (num1/denom1)*3600*1000;

end

%% Epiphyte leaf - atmosphere boundary layer
function ga = aerCond(hc, U, k, f1, f2)

% Function to calculate ga, aerodynamic conductance [m s-1]
%   hc = height of the canopy [m]
%   U = windspeed [m s-1]
%   k = von Karman's constant [unitless]

d = f1.*hc; % zero plane displacement height [m]
z0 = d./10; % Roughness length governing momentum transfer [m]
z0q = f2.*z0; % Roughness length governing transfer of heat and vapor [m]

ga = (U.* k.^2) ./ (log((hc-d)./z0).*log((hc-d)./z0q));

end

%% Heat capacity times mass (water + epiphyte mat)
function ct = m_cp(Se,c_pw,c_pd,RHOw,m_d)
% Function to compute the mass of soil and water in EWT
% Inputs    
%     Se = Water level in EWT [mm]
% Parameters
%     c_pw % Specific capacity of the water [W s kg-1 K-1]
%     c_pd % Specific capacity of the dry epiphyte mat [W s kg-1 K-1]  
%     RHOw % Water density [kg m-3] 
%     m_d% Epiphyte dry biomass [kg ha-1]

m_d=m_d/10000; % Epiphyte dry biomass [kg m-2]
% Calculations
m_Se = (Se/1000)*RHOw; % Mass of water per surface [kg m-2]
fs = m_d/(m_Se+m_d); % Influence of epiphyte mat in heat capacity
fw = m_Se/(m_Se+m_d); % Influence of water in heat capacity
ct = (m_Se+m_d)*((fw*c_pw)+(fs*c_pd)); % Mass per square meter times Heat capacity [W s m-2 K-1]  
    
end

%% Stomata conductance
function [sc,d2] = sc_i(gse,b1,SW_s,k2,Ta,Topt,Psi,Psi0,Psi_1,EAVD,Dx)
% Function to describe stomatal response
%     gse = maximum stomatal conductance [m s-1]
%     b1 = fitted parameter for stomatal sensitivity to SW_s [m2 W-1]
%     SW_s = incoming vertical SW radiation [W m-2]
%     k2 = coefficient for temperature stress [K-2]
%     Ta = Air temperature [K]
%     Topt = Optimum air temperature [K]
%     Psi = Water potential [m]
%     Psi0 = Minimum water potential [m]
%     Psi_1 = Maximum water potential [m]
%     EAVD = Epiphyte - air vapor pressure defficit [Pa]
%     Dx = temperature-dependent parameter that describes stomatal sensitivity to humidity [Pa]

% Stomatal sensitivity to incoming solar radiation
f_phi = min(1-exp(-b1*SW_s),1);

% Stomatal response to temperature
f_ta = min(1-(k2*((Ta-Topt)^2)),1);

% Stomatal sensitivity to water potential
if Psi < Psi0
    f_psi = 0;
elseif Psi > Psi_1
    f_psi = 1;
else
    f_psi = (Psi-Psi0)/(Psi_1-Psi0);
end

% Stomatal response to EAVD
f_eavd = min((1+(EAVD/Dx))^(-1),1);

% Stomatal response
d2=[gse,f_phi,f_ta,f_psi,f_eavd];
sc = gse*f_phi*f_ta*f_psi*f_eavd;
end

