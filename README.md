# Validation scripts
Validation scripts for air quality models.
Version 1.0.0, April 2022, Tereza Šedivá

A series of scripts for fast and easy validation of air quality (possibly other meteorological data) is provided.
These script were developed for validation of CMAQ and CAMS air quality model data, for other models there might have to be changes made to the scripts.
The validation script requires model outputs in a form of netCDF4 dataset.
The model data are validated against the observations measured at the National air quality monitoring stations in Slovakia (NMSKO).
The data from the NMSKO stations are available through the internal SHMU (Slovak hydrometeorological institute) database. 
For validation against other stations, you will need to link the validation function to your own database (it will be described below, in Set up for Validation)
All of the files for the validation are available on the master branch of this repository.
The validation is realized in a jupyter notebook called Validation_notebook.ipynb. 
Further I will explain how to use the script for your own purposes.

## Set up for Validation
You will need to download several files and prepare some data into a correct format for the validation to be possible.
First, make a /validation folder where you wish to save the validation data. The model outputs may be saved elsewhere.
There are several modules you need to download, which contain functions used in validation. These modules are day_of_year.py and validation_v2.py - save them into the /validation folder.
### Linking your own stations database into the validation function
The validation function is contained in the validation_v2.py module. To link your own stations database, do the following:


You will need a table containing data from the 








