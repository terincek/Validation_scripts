#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script for creation of stations table with x and y coordinates within the model domain

Last edit: 17.3.2021
@author: Tereza 
"""
import xarray as xr
import numpy as np
import pandas as pd

### USER INPUT ########
grid_name = ''
grid = xr.open_dataset('')   #path to grid in netCDF4 format
g_lat_name, g_lon_name = 'LAT', 'LON'  #names of LAT and LON fields in the netCDF4 file

#table with your stations, needs to contain columns with lats and lons for each station named 'lat' and 'lon'
meta_data = pd.read_csv('')  

save_path = '' #folder in which you want to save the table
table_name = 'stations_{}'.format(grid_name)

######################################################################################

# Function to find coordinates of the grid closest to the station
def getclosest_ij(lats,lons,latpt,lonpt):
    # find squared distance of every point on grid
    dist_sq = (lats-latpt)**2 + (lons-lonpt)**2  
    # 1D index of minimum dist_sq element
    minindex_flattened = dist_sq.argmin()    
    # Get 2D index for latvals and lonvals arrays from 1D index
    return np.unravel_index(minindex_flattened, lats.shape)

lats, lons = grid[g_lat_name][0,0,:,:].values, grid[g_lon_name][0,0,:,:].values 
meta_data['coord']=meta_data.apply(lambda x: getclosest_ij(lats, lons, x['lat'],x['lon']),axis=1 )
meta_data['x']=meta_data.apply(lambda x: getclosest_ij(lats, lons, x['lat'],x['lon'])[0],axis=1 )
meta_data['y']=meta_data.apply(lambda x: getclosest_ij(lats, lons, x['lat'],x['lon'])[1],axis=1 )

meta_data['lat_lon']=meta_data.apply(lambda x: (lats[x['coord']],lons[x['coord']]),axis=1 ) 

meta_data.to_csv('{}/{}'.format(save_path,table_name))






