# Validation scripts
Validation scripts for air quality models.
Version 1.0.0, April 2022, Tereza Šedivá, tereza.sediva@gmail.com

A series of scripts for fast and easy validation of air quality (possibly other meteorological data) is provided.
These script were developed for validation of CMAQ and CAMS air quality model data, for other models there might have to be changes made to the scripts.
The validation script requires model outputs in a form of netCDF4 dataset.
The model data are validated against the observations measured at the National air quality monitoring stations in Slovakia (NMSKO).
The data from the NMSKO stations are available through the internal SHMU (Slovak hydrometeorological institute) database. 
For validation against other stations, you will need to link the validation function to your own database (it will be described below, in Set up for Validation)
All of the files for the validation are available on the master branch of this repository.
The validation is realized in a jupyter notebook called Validation_notebook.ipynb. 
Further I will explain how to use the script for your own purposes.

Disclaimer: You may want or need to adjust the scripts for your own preference. I tried to make a lot of notes to all parts of the validation script so that the user will be able to use it on their own, if the user is comfortable with working in python.
In case you have any questions, something is not working right, you need help with your data etc., let me know on tereza.sediva@gmail.com.
I'm sure we'll figure it out.

## Set up for Validation
You will need to download several files and prepare some data into a correct format for the validation to be possible.
First, make a /validation folder where you wish to save the validation data. The model outputs may be saved elsewhere.
There are several modules you need to download, which contain functions used in validation. These modules are day_of_year.py and validation_v2.py - save them into the /validation folder.

### Setting up the validation_v2.py head
The validation_v2.py file contains the validation function and you will need to specify some paths, names, metadata etc. within this file.
1. Validation folder path: Set the "vali_path" variable (the first uncommented line in the script) to your /validation folder.
2. Metadata: metadata can be provided in a form of a .csv table or from a database. In case you are using your own database, make sure to import the needed modules, libraries or files to use your queries. The metadata table needs to contain station names, station IDs, station latitude and longitude, station type (traffic, background...), station location (rural, urban....). You will need to specify names of the corresponding variables in your metadata table. Specify which station types do you wish to validate, with labels as used in metadata (in case you want to validate all station types, please include all).
3. Variables: You will need to specify names of variables in your observations, corresponding model variable names, units in which you wish to display the variables, specify constants for variables (in case you need to multiply the model values to match the units of observations).
4. Grid: Specify path to your model grid. The grid needs to be in a netcdf4 format. Specify the names of latitude and longitude fields in your grid. Preferably the latitude and longitude fields are in 2D fields. If your grid has 1D fields for latitude and longitude, you will need to make a 2D grid out of them - we will do that later in the file, in the validation function.
5. Stations blacklist: In case you want to exclude some stations from the validation (wrong meta data, missing observations etc. - some stations tend to make problems later in the validation) you can specify their codes here.
6. Time shift: The observations are usually saved as the averages of the preceeding hour - e.g. 9 AM marks the average concentrations between 8-9 AM. However, models often use the opposite approach, where 9 AM marks the beginning of the 10th hour of the day. Therefore, model 9 AM means average values for 9-10 AM. In case your model saves the data this way, you need to account for this in the validation. When you set the shift_time variable to "True" the observations are going to be taken for one hour later, but they are going to be marked as the prevoius hour, so that they match the model data at a correct time. The time dsiplayed in the figures is going to be the model time - the hours mark the beggining of an hour.

### Setting up the validation function in validation_v2.py
Scroll in the validation_v2.py file until you find the validation function. The function requires the following inputs - path to the model output, model name, grid name, start and end times of the validation, variable to be validated. For each variable, the validation function needs to be run separately.
You need to adjust a few things before running the validation function for the first time.
1. Adjusting the grid: Different models will use different grid types. In my example (in validation function), I work with models CMAQ and CAMS. The CMAQ saves LAT and LON fields in grid file with 4 dimensions (time step, layer, x and y), but time step and layer are constant, therefore it is enough to take X and Y dimensions to get 2d fields for lats and lons.
With CAMS, there were 1d fields for both latitude and longitude. I had to first combine them into 2d fields for the rest of the validation function to run smoothly. The way you do this depends on the way your data are provided. 
You want to end up with 2d field of latitudes and a second 2d field of longitudes for each grid cell.
DISCLAIMER: be sure to check the orientation of your lat and lon fields. The lat[0,0] and lon[0,0] values must both correspond to your model_variable[0,0] value possition. It is possible, especially when making a 2d field out of 1d, that you accidentaly create a transposed field (possibly only for either lats or lons or both).
There is a collumn of station's lat and lon in the model for you to reverse check whether you got values close to the real station's lat and lon.
2. Model data: You need to specify the way the model data will be loaded. This is usually different for different models. In my example, the CAMS model has a specific model output file for every day. The model CMAQ contains the whole validation period within one file. You will need to link the model time field to the date-time generated in the 'YYYY-MM-DD HH:MM:SS' format.
3. Observations: 


Ked je validovane dlhsie obdobie ako jeden rok tak nespravne mam urceny pocet dni validacie
time shift - asi by bolo lepsie to spravit opacne - teda posunut modelove data ci?
observations
date in observations







