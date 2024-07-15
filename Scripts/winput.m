function data = winput(ndays)
% Function to create a weather data set of meterological variables that
% change in diurnal cycles.
% Input:
%   - ndays: number of days to create a dataset
% Output:
%   - data: table with variables Date, time of the day in hours (t_h), air
%   temperature in Kelvin (Ta_K), relative humidity (RH), rainfall or
%   vertical precipitation (RF_mm), fog or horizontal precipitation (F_mm),
%   solar shortwave radiation (SW_Wpms).

%% Inputs that vary through time
% We are modeling data for now
n=ndays; %Number of days to model
% Date time variable in steps of 1 hour
Date = transpose(datetime(2022,04,13,0,0,0):hours(1):(datetime(2022,04,13,23,45,0)+(n-1))); % Date-time
% d = transpose(0:(60/(24*60)):n-(60/(24*60)));
%t = repmat(transpose(0:0.25:23.75),n,1); % Time of the day for n days
% t = repmat(transpose(0:1:23),n,1); % Time of the day for n days
% Time variable
t = hour(Date);
% Setting parameters of diurnal cycle
t0 = 5; % Initial sunlight hour (hours)
delta = 13; % Day length [h]
    
% Radiation modelation [W m^{-2}]
SWmax = 1100; %Maximum available energy [W m-2] 
SW_Wpms = max(0,(4*SWmax/(delta^2))*( -t.^2 + (delta+(2*t0)).*t - t0*(t0+delta) )); 

% Air temperature modelation [K]
Ta_range = 7.93; %Range in absolute temp values used to paramtrize diurnal temp (k)
Ta_K = 289.45 + max(0,(4*Ta_range/(delta^2))*( -t.^2 + (delta+(2*t0)).*t - t0*(t0+delta) ));

% Specific humidity of the air [kg kg-1]
% qrange  = 0.002476; % range in specific humidity used to parametrize diurnal humidity (kg/kg)
% q_a = 0.010844 + max(((4.*qrange)./((delta).^2)).*(-(mod(t,24).^2)+...
%     (((delta)+2.*t0).*mod(t,24))-(t0.*(t0+(delta)))),0);
% Relative humidity [-]
rh_range = 0.24;
RH = 0.94 - max(0,(4*rh_range/(delta^2))*( -t.^2 + (delta+(2*t0)).*t - t0*(t0+delta) ));
    
% Precipitation [mm]
RF_mm = zeros(length(t),1); 
    
% OP [mm]
F_mm = zeros(length(t),1); 
    
% Output as a table
data = table(Date,Ta_K,RH,RF_mm,F_mm, SW_Wpms);
    
end
