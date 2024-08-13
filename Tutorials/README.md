## Background
Details of the basis of the model and information on its peer-reviewed publication can be found on the [home page](https://github.com/DavidCarMor/EWB) of this repository. 

The code for this model was built in [MATLAB](https://www.mathworks.com/products/matlab.html), whose license usually can be accessed through your academic institution. The EWB code can be downloaded from the [Scripts](https://github.com/DavidCarMor/EWB/tree/main/Scripts) folder of this repository.

The forcing variables of the model are air temperature in Kelvin degrees, relative humidity in unit fraction, rainfall in millimeters, fog in millimeters, and shortwave radiation in Watts per meter squared. The EWB model works with a one-hour resolution, and its forcing data shall be included in the same time resolution. Additionally, atmospheric pressure in pascals and/or wind speed in meters per second can be added as forcing data for the model. Atmospheric pressure affects the host-tree water uptake from the epiphyte mat, and wind speed is used to replace the fixed boundary layer assumed in the model for a dynamic boundary layer.

## 1. Preallocating variables
1.1. Download the EPI_WBV9 and winput functions from the [Scripts](https://github.com/DavidCarMor/EWB/tree/main/Scripts) folder, open matlab and set your working directory to the same where those functions were stored.<br />
![step1](https://github.com/user-attachments/assets/f0663e70-33d0-4d1d-8a4f-f8b8fabb719d)

1.2. We will create a MATLAB script to work with the model for this tutorial.<br />
![Step2](https://github.com/user-attachments/assets/29048e14-0401-4de9-b25e-21993e9f0afa)

1.3. As mentioned in the background, the model requires air temperature, relative humidity, rainfall and fog, and shortwave radiation as forcing variables. We will create ideal forcing data using the function [winput](https://github.com/DavidCarMor/EWB/tree/main/Scripts/winput.m). Nevertheless, the user is encouraged to test the model using its own forcing data as a MATLAB table with the following variables:
    - Date: Variable storing date time information in the 'yyyy-MM-dd hh:mm:ss' format,
    - Ta_K: Air temperature in Kelvin degrees,
    - RH: Relative humidity of the air in unit fraction,
    - RF_mm: Rainfall rate in millimeters,
    - F_mm: Fog rate in millimeters, and
    - SW_Wpms: Shortwave radiation in Watts per meter squared.

  To create ideal forcing data with the winput function, use the following lines of code:
  
    n = 3; % Number of days of simulation
    data = winput(n); % Ideal forcing data
You will notice that rainfall and fog rates are constant vectors equal to zero as this ideal data re-creates a dry-down event.
 
## 2. Simulation with default parameters
2.1 The model includes default biotic and abiotic parameters that can be used to run our first simulation. There are three outputs from the model: a table with the variables used in the simulation, parameters used in the model, and variability of the stomata conductance of the vascular component of the epiphyte mat. 

The following lines of code show how to store the three outputs in our MATLAB workspace using the function [EPI_WBV9](https://github.com/DavidCarMor/EWB/tree/main/Scripts/EPI_WBV9.m). If needed, comment on the third line and uncomment one of the other two alternatives according to the desired number of outputs. 

    % [output_data] =  EPI_WBV9(data); % one output
    % [output_data, parameters] =  EPI_WBV9(data); % two outputs
    [output_data, parameters, stomatal_info] =  EPI_WBV9(data); % three outputs

2.2 As you might have noticed, the outputs are on your workspace and we can plot the new data now.
    % Plotting dew deposition and evapotranspiration
    figure
    hold on 
    plot(output_data.Date,output_data.DDe,'-','LineWidth',2)
    plot(output_data.Date,output_data.ETe,'-','LineWidth',2)
    hold off
    ylabel('DD_e and ET_e[mm/h]', 'FontSize',12)
    legend('Dew deposition', 'Evapotranspiration', 'Location', 'southoutside', ...
    'Orientation','horizontal','Fontsize',12)

    % Plotting temperature
    figure
    hold on 
    plot(output_data.Date,output_data.Ta_C,'-','LineWidth',2)
    plot(output_data.Date,output_data.Te_C,'-','LineWidth',2)
    ylabel('Temperature [°C]', 'FontSize',12)
    hold off
    legend('Air', 'Epiphyte mat', 'Location', 'southoutside', ...
        'Orientation','horizontal','Fontsize',12)


## 3. Simulation with new parameters
There are different ways to access the information regarding the default parameters in the EWB model.
- Accessing the peer-reviewed article on the model: Carchipulla-Morales, D., Corbett, H., Vaughan, D., Gotsch, S., Dawson, T., Nadkarni, N., and Lowman, L., 2024. A novel model to simulate water and energy budgets for epiphytic mats.
- Using the TEXT file listing the parameters of the mode: [Parameters.txt](https://github.com/DavidCarMor/EWB/tree/main/Parameters.txt)
- Using the assist for the EPI_WBV9 function on the MATLAB command window: help EPI_WBV9

3.1  First, we will run a couple of simulations changing the initial water level in the epiphyte mat. The following lines show how to change the initial water level using two valid coding syntaxis. 

    output_Se15 =EPI_WBV9(data,So=15); % So=15
    output_Se25 =EPI_WBV9(data,'So',25); % So=25

3.2 Now, we can plot and compare those results

    % Plotting water level
    figure
    hold on 
    plot(output_Se15.Date,output_Se15.Se,'-','LineWidth',2)
    plot(output_Se25.Date,output_Se25.Se,'-','LineWidth',2)
    ylabel('Water level [mm]', 'FontSize',12)
    hold off
    legend('S_0 = 15 mm', 'S_0 = 25 mm', 'Location', 'southoutside', ...
        'Orientation','horizontal','Fontsize',12)

3.3 It is possible to change more than one parameter at a time. The following lines of code will simulate the dry down when the initial water level equals 20 mm, and the Leaf Area Index of the forest equals 5 m2/m2 and the fractional grid canopy area covered by epiphytes equals 0.5:
    
    output_Se20_LAI_5 =EPI_WBV9(data,So=20, LAI=5); % So=20 and LAI=5
    output_Se20_fe05 =EPI_WBV9(data,So=20,fe=0.5); % So=20 and fe=0.5

3.2 Now, we can plot and compare those results

    % Plotting water level under different conditions
    figure
    hold on 
    plot(output_Se20_LAI_5.Date,output_Se20_LAI_5.Se,'-','LineWidth',2)
    plot(output_Se20_fe05.Date,output_Se20_fe05.Se,'-','LineWidth',2)
    ylabel('Water level [mm]', 'FontSize',12)
    hold off
    legend('S_0 = 20 mm \newline LAI = 5', 'S_0 = 20 mm \newline f_e = 0.5', 'Location', 'southoutside', ...
            'Orientation','horizontal','Fontsize',12)
Despite these simulations being initiated with the same water level, we can observe steeper changes in the water levels of the second curve.

## 4. Exporting data
It is possible to export the simulation results to a local machine. Nevertheless, it is suggested to turn the timetable object that stores the output data from the model into a table object before exporting the results to a local machine.

    output_table = timetable2table(output_data); % Changing the format of the object
    writetable(output_table,'./results_test1.csv') % './' saves the file in your current working directory

### Full script
    %% Ideal data
    n = 2; % Number of days of simulation
    data = winput(2); % Ideal forcing data
    
    %% Simulation with default parameters
    % [output_data] =  EPI_WBV9(data); % one output
    % [output_data, parameters] =  EPI_WBV9(data); % two outputs
    [output_data, parameters, stomatal_info] =  EPI_WBV9(data); % three outputs
    
    % Plotting dew deposition and evapotranspiration
    figure
    hold on 
    plot(output_data.Date,output_data.DDe,'-','LineWidth',2)
    plot(output_data.Date,output_data.ETe,'-','LineWidth',2)
    hold off
    ylabel('DD_e and ET_e[mm/h]', 'FontSize',12)
    legend('Dew deposition', 'Evapotranspiration', 'Location', 'southoutside', ...
        'Orientation','horizontal','Fontsize',12)
    
    % Plotting temperature
    figure
    hold on 
    plot(output_data.Date,output_data.Ta_C,'-','LineWidth',2)
    plot(output_data.Date,output_data.Te_C,'-','LineWidth',2)
    ylabel('Temperature [°C]', 'FontSize',12)
    hold off
    legend('Air', 'Epiphyte mat', 'Location', 'southoutside', ...
        'Orientation','horizontal','Fontsize',12)
    
    %% Simulation with new parameters
    output_Se15 =EPI_WBV9(data,So=15); % So=15
    output_Se25 =EPI_WBV9(data,'So',25); % So=25
    
    % Plotting water level
    figure
    hold on 
    plot(output_Se15.Date,output_Se15.Se,'-','LineWidth',2)
    plot(output_Se25.Date,output_Se25.Se,'-','LineWidth',2)
    ylabel('Water level [mm]', 'FontSize',12)
    hold off
    legend('S_0 = 15 mm', 'S_0 = 25 mm', 'Location', 'southoutside', ...
        'Orientation','horizontal','Fontsize',12)
    
    % Changing initial water level and biomass
    output_Se20_LAI_5 =EPI_WBV9(data,So=20,LAI=5); % So=20 and LAI=5
    output_Se20_fe05 =EPI_WBV9(data,So=20,fe=0.5); % So=20 and fe=0.5
    
    
    % Plotting water level under different conditions
    figure
    hold on 
    plot(output_Se20_LAI_5.Date,output_Se20_LAI_5.Se,'-','LineWidth',2)
    plot(output_Se20_fe05.Date,output_Se20_fe05.Se,'-','LineWidth',2)
    ylabel('Water level [mm]', 'FontSize',12)
    hold off
    legend('S_0 = 20 mm \newline LAI = 5', 'S_0 = 20 mm \newline f_e = 0.5', 'Location', 'southoutside', ...
            'Orientation','horizontal','Fontsize',12)

    %% Exporting data
    output_table = timetable2table(output_data); % Changin object format
    writetable(output_table,'./results_test1.csv') % './' saves the file in your current working directory
