## Background
Details of the basis of the model and information on its peer-reviewed publication can be found on the [home page](https://github.com/DavidCarMor/EWB) of this repository. 

The code for this model was built in [MATLAB](https://www.mathworks.com/products/matlab.html), whose license usually can be accessed through your academic institution. The EWB code can be downloaded from the [Scripts](https://github.com/DavidCarMor/EWB/tree/main/Scripts) folder of this repository.

The forcing variables of the model are air temperature in Kelvin degrees, relative humidity in unit fraction, rainfall in millimeters, fog in millimeters, and shortwave radiation in Watts per meter squared. Additionally, atmospheric pressure in pascals and/or wind speed in meters per second can be added as forcing data for the model. Atmospheric pressure affects the host-tree water uptake from the epiphyte mat, and wind speed is used to replace the fixed boundary layer assumed in the model for a dynamic boundary layer.

## Preallocating variables
1. Download the EPI_WBV9 and winput functions from the [Scripts](https://github.com/DavidCarMor/EWB/tree/main/Scripts) folder, open matlab and set your working directory to the same where those functions were stored.<br />
![step1](https://github.com/user-attachments/assets/f0663e70-33d0-4d1d-8a4f-f8b8fabb719d)

2. We will create a MATLAB script to work with the model for this tutorial.<br />
![Step2](https://github.com/user-attachments/assets/29048e14-0401-4de9-b25e-21993e9f0afa)

3. As mentioned in the background, the model requires air temperature, relative humidity, rainfall and fog, and shortwave radiation as forcing variables. We will create ideal forcing data by using the function [winput](https://github.com/DavidCarMor/EWB/tree/main/Scripts/winput.m). Nevertheless, the user is encouraged to test the model using its own forcing data as a MATLAB table with the following variables:
    - Date: Variable storing date time information in the 'yyyy-MM-dd hh:mm:ss' format,
    - Ta_K: Air temperature in Kelvin degrees,
    - RH: Relative humidity of the air in unit fraction,
    - RF_mm: Rainfall rate in millimeters,
    - F_mm: Fog rate in millimeters, and
    - SW_Wpms: Shortwave radiation in Watts per meter squared.

  To create ideal forcing data with the winput function, use the following lines of code:
      %% Ideal data
    n = 2; % Number of days of simulation
    data = winput(2); % Ideal forcing data
You will notice that rainfall and fog rates are defined to be equal to zero as this ideal data is re-creating a dry-down event.
 
## Simulating balance with default parameters

### Example
    python ee_Landsat_LAI_export.py -o <asset_dir> -p <path> -r <row> 
        -s <start_date> -d <end_date>
        
    Required arguments:
    -o  Earth Engine asset directory to export LAI images, can a foler or 
        an image collection
    -p  WRS Path number of the Landsat Collection 1 surface reflectance scene
    -r  WRS Row number of the Landsat Collection 1 surface reflectance scene
    -s  The start date to export image in YYYY-MM-dd
    -e  The end date (exclusive) to export image in YYYY-MM-dd

## Modifying parameters

### Example
For example, 

    python ee_Landsat_LAI_export_v0.1.1.py -o projects/ee-yanghuikang/assets/LAI_test/LAI_test_v0_1_1 -p 44 -r 33 -s 2020-06-01 -e 2020-06-15
The output will look like, 

    assetDir: projects/ee-yanghuikang/assets/LAI_test/LAI_test_v0_1_1
    WRS path: 44
    WRS row: 33
    start date: 2020-06-01
    end date: 2020-06-15
    Number of Landsat images:  2
    CFBCSQCUGNZMKSDXF3HPE7EO LAI_LC08_044033_20200606
    X6JDL76R5BVHTW7JAGHFEJOA LAI_LE07_044033_20200614

This will export two LAI images to the designated EE asset image collection. The last two lines print the Task ID and the image name.
