#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Validacia vystupov z modelu s meraniami stanic NMSKO.
verzia kde vsetky pollutanty su v jednom subore

Created on Tue Jan 12 10:26:35 2021

@author: ext29843
"""
vali_path = '/data/teri/validacia'
import numpy as np
import pandas as pd
import xarray as xr
import sys
import importlib
import time
import datetime as dt
#import random
#import matplotlib.pyplot as plt
#import os

sys.path.append(vali_path)
sys.path.append('/home/KOL/p6065/python-scripts/core')
#import val_config2 as config
#importlib.reload(config)
import day_of_year as doy

from dbConnector import obs
    
# Find coordinates of the grid closest to the station
def getclosest_ij(lats,lons,latpt,lonpt):
    # find squared distance of every point on grid
    dist_sq = (lats-latpt)**2 + (lons-lonpt)**2  
    # 1D index of minimum dist_sq element
    minindex_flattened = dist_sq.argmin()    
    # Get 2D index for latvals and lonvals arrays from 1D index
    return np.unravel_index(minindex_flattened, lats.shape)

query_meta = obs.q_nmsko 
meta_data_ = obs.query(query_meta)

var = ['SO2','NO2','PM10','PM2_5','O3','CO']   # names of variables in observations
var_mod, var_const = {}, {}
var_mod['CAMS'] = {'SO2':'so2_conc', 'NO2':'no2_conc', 'PM10':'pm10_conc', 'O3':'o3_conc', 'PM2_5':'pm2p5_conc','CO':'co_conc'}
var_mod['CMAQ'] = {'SO2':'SO2', 'NO2':'NO2', 'PM10':'PM10', 'O3':'O3', 'PM2_5':'PM25', 'CO':'CO'}
var_const['CAMS'] = {'NO2':1,'O3':1, 'SO2':1, 'CO':1, 'PM10':1, 'PM2_5':1} #vystupy z CAMS su v ug/m3
var_const['CMAQ'] = {'NO2':1.92,'O3':2, 'SO2':2.66, 'CO':1.16, 'PM10':1, 'PM2_5':1}
units = 'ug/m3'
station_type = ['B','T','I']     # stations to include in the validation as named in meta data

station_blacklist = [99701,99702,99703,99501,99502,99503,99350,99223,99130,99126]
meta_data = meta_data_[~meta_data_.ii.isin(station_blacklist)].copy()

def is_leap(year):
    leap = False
    if year%4 == 0:
        leap = True
        if year%100 == 0:
            if year%400 == 0:
                leap = True
            else:
                leap = False
    
    return leap

# v nasej databaze znamena 12 hodin priemer za 11-12 kdezto v modeli 12 znamena priemer za 12-13. takze pozorovania treba zobrat o hodinu neskor a nasledne posunut casovu os aby sedela na modelovu
def obs_date(start_date, end_date):
    year = int(start_date[0:4])
    if is_leap(year):   
        month_ends = ['01-31','02-29','03-31','04-30','05-31','06-30','07-31','08-31','09-30','10-31','11-30']       
    else:
        month_ends = ['01-31','02-28','03-31','04-30','05-31','06-30','07-31','08-31','09-30','10-31','11-30']
                
    if end_date[5:10] in month_ends:
        new_month = str(int(end_date[5:7]) +1)
        if int(new_month) < 10:
            new_month = '0{}'.format(new_month)
        obs_end_date = end_date[0:5] + new_month + '-01 00:00:00'
        
    elif end_date[5:10] == '12-31':
        obs_end_date = str(int(end_date[0:4]) +1) + '-01-01 00:00:00'
        
    else:
        new_day = str(int(end_date[8:10])+1)
        if int(new_day) < 10:
            new_day = '0{}'.format(new_day) 
        obs_end_date = end_date[0:8] + new_day + ' 00:00:00' 
        
    obs_start_date = start_date[0:10]+' 01:00:00'
    
    return obs_start_date, obs_end_date
    
"""
start_date = '2021-11-01 00:00:00'
end_date = '2021-11-30 23:59:59'
#grid_name = 'cams'
#model_name = 'CAMS'
grid_name = 'ala2km'
#grid_name = 'd02'
model_name = 'CMAQ'
#path = '/data/teri/cmaq_data/COMBINE_ACONC_v531_gcc_ala_2021_202103.nc'
#path = '/data/oko/dusan/CAMS_DATA/CAMS_analysis_ensamble_med_'
#path = '/data/oko/dusan/POST/LIFEIP_Small_NO2_2015'   #d02
#path = '/data/oko/dusan/POST/WB_2020_PM_2017'  #d02
#path = '/data/oko/teri/cmaq_data/WB_2017_d02_CAMS_13_res_heat_PM_2017'  #d02
path = '/data/oko/teri/cmaq_data/COMBINE_ACONC_v533_intel_ala_2021_202111.nc'
"""
#validation('/data/teri/cmaq_data/COMBINE_ACONC_v531_gcc_ala_2021_202103.nc', 'CMAQ', 'ala2km','2021-03-05','2021-03-09 23:59:59', 'NO2')

def validation(path, model_name, grid_name, start_date, end_date, v):
    start_time = time.time()
    #print('Validating...')
        
    grid_dict = {'ala2km':'/data/teri/cmaq_data/GRIDCRO2D_2021-03-05.nc', 
            'cams':'/data/oko/dusan/CAMS_DATA/CAMS_analysis_ensamble_med_2021-03-05.nc',
            'd02':'/data/juraj/MEGAN/MEGANv2.10/Input_d02/MAP/EFMAPS.CMAQ_d02_grid.nc'}
    grid = xr.open_dataset(grid_dict[grid_name])
    
    lat_dict = {'CMAQ':'LAT','CAMS':'latitude'}
    lon_dict = {'CMAQ':'LON','CAMS':'longitude',}
    
    if model_name == 'CMAQ' and grid_name =='d02':
        lon_dict['CMAQ'] = 'LONG'
            
    if model_name == "CMAQ":
        lats, lons = grid[lat_dict[model_name]][0,0,:,:].values, grid[lon_dict[model_name]][0,0,:,:].values 

    elif model_name == "CAMS":
        #g_lat_name, g_lon_name = 'latitude', 'longitude'
        grid_lats = np.array(grid[lat_dict[model_name]])
        grid_lons = np.array(grid[lon_dict[model_name]])
        lats = np.transpose([grid_lats] * len(grid_lons))
        lons = np.tile(grid_lons, (len(grid_lats), 1))

    meta_data['coord']=meta_data.apply(lambda x: getclosest_ij(lats, lons, x['lat'],x['lon']),axis=1 )
    meta_data['lat_lon']=meta_data.apply(lambda x: (lats[x['coord']],lons[x['coord']]),axis=1 ) 
    
    
    n_start = doy.day_of_year(start_date)
    n_end = doy.day_of_year(end_date)
    n_days = n_end - n_start + 1
    
    obs_start_date, obs_end_date = obs_date(start_date, end_date)

    #for v in var:
    if 1<2:
        
        times = pd.date_range(start=start_date, periods=n_days*24, freq='H')
        len(times)

        if model_name == 'CAMS':
            model = xr.open_dataset('{0}{1}.nc'.format(path,start_date[0:10]))[var_mod[model_name][v]]
            #n_first = str(start_date[])
            for i in range(n_start+1,n_end + 1):
                date = doy.date_from_day(i, start_date[0:4])
                model_new = xr.open_dataset('{0}{1}.nc'.format(path,date))[var_mod[model_name][v]]
                model = xr.concat([model, model_new], 'time')
            
            model.coords['time'] = (('time'),times)
        
        else:                
            #model = xr.open_dataset(path)[var_mod[model_name][v]][(n_start-1)*24:n_end*24,0,:,:] #pokial je modelovy subor pre cely rok, ale validovat chceme len dake obdobie
            model = xr.open_dataset(path)[var_mod[model_name][v]][:,0,:,:]
            model.coords['TSTEP'] = (('TSTEP'),times)
        
        obs_all = pd.DataFrame()
        mod_all = pd.DataFrame()
        dictionary, dictionary2, dictionary3, dictionary4 = {}, {}, {}, {}
        dictionary['all'], dictionary2['all'], dictionary3['all'], dictionary4['all'] = 'all', 'all', 'all', 'all'
        
        #s = 99112
        for s in meta_data['ii']:
        #if 1<2:
            query = "select si.ii,si.lat,si.lon,si.elev, o.* from obs.obs_nmsko_1h o join si.si on si.id=o.si_id where ii={0} and date between '{1}' and '{2}'".format(s, obs_start_date, obs_end_date)
            obs_data = obs.query(query)
            obs_data = obs_data.set_index('date')
            #data_conc = data_conc.reindex(pd.date_range("2015-01-01 00:00:00", koniec , freq='H'))
            obs_data = obs_data.reindex(pd.date_range(obs_start_date, obs_end_date, freq = 'H'))
            #posun osi o 1 hodinu dozadu
            obs_data.index = obs_data.index + dt.timedelta(hours=-1)
            #obs_all[s] = obs_data[v][(n_start-1)*24:n_end*24]
            obs_all[s] = obs_data[v][:]
            
            coord_sta = np.array(meta_data['coord'].loc[meta_data['ii']==s].item())
            model_data = model[:,coord_sta[0],coord_sta[1]].values*var_const[model_name][v]
            mod_all[s] = obs_data[v]
            mod_all[s] = model_data
            dictionary[s] = '{0}:{1}'.format(s,meta_data['name'].loc[meta_data['ii']==s].values[0])
            dictionary2[s] = meta_data['name'].loc[meta_data['ii']==s].values[0]
            dictionary3[s] = meta_data['type'].loc[meta_data['ii']==s].values[0]
            dictionary4[s] = meta_data['loc'].loc[meta_data['ii']==s].values[0]
                
        ######### STATISTICS ##########
        stat = pd.DataFrame(index = dictionary.keys(), columns = ['name','loc','type','coverage','obs_mean','MB','MGE','RMSE','R','FAC2'])
        stat['name'] = dictionary2.values()
        stat['type'] = dictionary3.values()
        stat['loc'] = dictionary4.values()
        obs_all.fillna(value=np.nan, inplace=True)
        
        for s in stat.index:
            if s == 'all': 
                obs_all[obs_all == 0] = np.nan
                nan_len = obs_all.dropna(axis = 1, how = 'all').isna().sum().sum()
                #all_len = obs_all.dropna(axis = 1, how = 'all').count().sum()
                all_len = len(obs_all)*len(obs_all.dropna(axis = 1, how = 'all').columns)
                stat['coverage'][s] = 100*(all_len - nan_len)/all_len
                diff = (mod_all - obs_all).replace('nan',np.nan).to_numpy().flatten() 
                #obs_flat = obs_all.replace('nan',np.nan).to_numpy().flatten() 
                #obs_flat = obs_flat[~np.isnan(obs_flat)]
                diff = diff[~np.isnan(diff)]  #ponecha iba nenanove cleny                
                podiel = (mod_all/obs_all).replace('nan',np.nan).to_numpy().flatten()
                podiel = podiel[~np.isnan(podiel)]
                podiel_fac2 = podiel[podiel <= 2]
                podiel_fac2 = podiel_fac2[podiel_fac2 >= 0.5]
                
                corr_mod = mod_all.stack()
                corr_obs = obs_all.replace('nan',np.nan).stack(dropna = False)
                stat['R'][s] = corr_mod.corr(corr_obs)
                stat['obs_mean'][s] = np.nanmean(obs_all)
                #stat['FAC2'] = 100*len(podiel_fac2)/len(podiel)
        
            else:
                obs_all[s][obs_all[s] == 0] = np.nan
                nan_len = obs_all[s].isna().sum()
                all_len = len(mod_all[s])
                stat['coverage'][s] = 100*(all_len - nan_len)/all_len
                diff = np.array((mod_all[s] - obs_all[s]).dropna()) 
                #obs_flat = obs_all[s].replace('nan',np.nan).to_numpy()
                #obs_flat = obs_flat[~np.isnan(obs_flat)]
                podiel = (mod_all[s]/obs_all[s]).dropna()
                podiel_fac2 = podiel[podiel <= 2]
                podiel_fac2 = podiel_fac2[podiel_fac2 >= 0.5]
                #stat['FAC2'] = 100*len(podiel_fac2)/len(podiel)
                #R
                stat['R'][s] = mod_all[s].corr(obs_all[s])
                stat['obs_mean'][s] = np.nanmean(obs_all[s])
                
            #obs_flat[obs_flat == 0] = 0.1
                    
            if len(diff) != 0:
                #BIAS
                stat['MB'][s] = np.sum(diff)/len(diff)
                #stat['MPE'][s] = 100*np.sum((diff/obs_flat))/len(diff)
                #stat['MB%'][s] = 100*(np.sum(diff)/len(diff))/np.mean(obs_flat)
                #MGE
                stat['MGE'][s] = np.sum(np.absolute(diff))/len(diff)
                #stat['MGPE'][s] = 100*np.sum(np.absolute(diff/obs_flat))/len(diff)
                #stat['MGE%'][s] = 100*(np.sum(np.absolute(diff))/len(diff))/np.mean(obs_flat)
                #RMSE
                stat['RMSE'][s] = np.sqrt(np.sum(diff**2)/len(diff))
                #stat['RMSPE'][s] = 100*np.sqrt(np.sum((diff/obs_flat)**2)/len(diff))
                #stat['RMSE%'][s] = 100*(np.sqrt(np.sum(diff**2)/len(diff)))/np.mean(obs_flat)
                stat['FAC2'][s] = 100*len(podiel_fac2)/len(podiel)
        

        print('Validation of {0} done in {1} seconds'.format(v,time.time()-start_time))
                
    ###### SAVE DATA AS NETCDF4                                                
        val_data = xr.Dataset()            
        # data variables
        #val_data['OBS'] = (('time','station'), np.concatenate((obs_data[v].to_numpy(),np.array([s]))), axis = 1)
        val_data['OBS'.format(v)] = (('time','station'), obs_all)
        val_data['MOD'.format(v)] = (('time','station'), mod_all)
           
        # coordinates
        val_data.coords['time'] = times
        val_data.coords['station'] = obs_all.columns.values
           
        # attributes
        val_data.attrs['Description'] = 'Validation of the {0} model with the NMSKO stations. For each station, observations and model data for the selected period are stored as 1D arrays.'.format(model_name)
        val_data.attrs['pollutants'] = v
        val_data.attrs['units'] = units
        val_data.attrs['domain'] = grid_name
        val_data.attrs['Included station types'] = station_type
        val_data.attrs['Period'] = '{0}- {1}'.format(start_date[0:11], end_date[0:11])
        val_data.attrs['Stations'] = list(dictionary.values())
        val_data.attrs['Get stations as dictionary']: "stations_list = xr_data.attrs['Stations']; stations = {}; for i in range(0,len(stations_list)): stations[stations_list[i][0:5]]=stations_list[i][6:]"

    return stat, val_data

def comparison(ref_path, other_path, ref_name, other_name, ref_grid_name, other_grid_name, start_date, end_date, v):
    start_time = time.time()
    #print('Validating...')
    grid_dict = {'ala2km':'/data/teri/cmaq_data/GRIDCRO2D_2021-03-05.nc', 
            'cams':'/data/oko/dusan/CAMS_DATA/CAMS_analysis_ensamble_med_2021-03-05.nc',
            'd02':'/data/juraj/MEGAN/MEGANv2.10/Input_d02/MAP/EFMAPS.CMAQ_d02_grid.nc'}
    ref_grid = xr.open_dataset(grid_dict[ref_grid_name])
    other_grid = xr.open_dataset(grid_dict[other_grid_name])
    
    lat_dict = {'CMAQ':'LAT','CAMS':'latitude'}
    lon_dict = {'CMAQ':'LON','CAMS':'longitude',}
    
    if ref_name == 'CMAQ' and ref_grid_name =='d02':
        lon_dict['CMAQ'] = 'LONG'
        
    if other_name == 'CMAQ' and other_grid_name =='d02':
        lon_dict['CMAQ'] = 'LONG'
            
    if ref_name == "CMAQ":
        ref_lats, ref_lons = ref_grid[lat_dict[ref_name]][0,0,:,:].values, ref_grid[lon_dict[ref_name]][0,0,:,:].values 

    elif ref_name == "CAMS":
        #g_lat_name, g_lon_name = 'latitude', 'longitude'
        ref_grid_lats = np.array(ref_grid[lat_dict[ref_name]])
        ref_grid_lons = np.array(ref_grid[lon_dict[ref_name]])
        ref_lats = np.transpose([ref_grid_lats] * len(ref_grid_lons))
        ref_lons = np.tile(ref_grid_lons, (len(ref_grid_lats), 1))
        
    if other_name == "CMAQ":
        other_lats, other_lons = other_grid[lat_dict[other_name]][0,0,:,:].values, other_grid[lon_dict[other_name]][0,0,:,:].values 
        
    elif other_name == "CAMS":
        #g_lat_name, g_lon_name = 'latitude', 'longitude'
        other_grid_lats = np.array(other_grid[lat_dict[other_name]])
        other_grid_lons = np.array(other_grid[lon_dict[other_name]])
        other_lats = np.transpose([other_grid_lats] * len(other_grid_lons))
        other_lons = np.tile(other_grid_lons, (len(other_grid_lats), 1))

    meta_data['ref_coord']=meta_data.apply(lambda x: getclosest_ij(ref_lats, ref_lons, x['lat'],x['lon']),axis=1 )
    meta_data['ref_lat_lon']=meta_data.apply(lambda x: (ref_lats[x['ref_coord']],ref_lons[x['ref_coord']]),axis=1 ) 
    meta_data['other_coord']=meta_data.apply(lambda x: getclosest_ij(other_lats, other_lons, x['lat'],x['lon']),axis=1 )
    meta_data['other_lat_lon']=meta_data.apply(lambda x: (other_lats[x['other_coord']],other_lons[x['other_coord']]),axis=1 ) 
    
    n_start = doy.day_of_year(start_date)
    n_end = doy.day_of_year(end_date)
    n_days = n_end - n_start + 1

    #for v in var:
    if 1<2:
        
        times = pd.date_range(start=start_date, periods=n_days*24, freq='H')

        if ref_name == 'CAMS':
            ref_model = xr.open_dataset('{0}{1}.nc'.format(ref_path,start_date[0:10]))[var_mod[ref_name][v]]
            #n_first = str(start_date[])
            for i in range(n_start+1,n_end + 1):
                date = doy.date_from_day(i, start_date[0:4])
                ref_model_new = xr.open_dataset('{0}{1}.nc'.format(ref_path,date))[var_mod[ref_name][v]]
                ref_model = xr.concat([ref_model, ref_model_new], 'time')
            
            ref_model.coords['time'] = (('time'),times)
        
        else:                
            ref_model = xr.open_dataset(ref_path)[var_mod[ref_name][v]][(n_start-1)*24:n_end*24,0,:,:]
            ref_model.coords['TSTEP'] = (('TSTEP'),times)
            
        if other_name == 'CAMS':
            other_model = xr.open_dataset('{0}{1}.nc'.format(other_path,start_date[0:10]))[var_mod[other_name][v]]
            #n_first = str(start_date[])
            for i in range(n_start+1,n_end + 1):
                date = doy.date_from_day(i, start_date[0:4])
                other_model_new = xr.open_dataset('{0}{1}.nc'.format(other_path,date))[var_mod[other_name][v]]
                other_model = xr.concat([other_model, other_model_new], 'time')
            
            other_model.coords['time'] = (('time'),times)
        
        else:                
            other_model = xr.open_dataset(other_path)[var_mod[other_name][v]][(n_start-1)*24:n_end*24,0,:,:]
            other_model.coords['TSTEP'] = (('TSTEP'),times)
        
        ref_all = pd.DataFrame()
        other_all = pd.DataFrame()
        dictionary, dictionary2, dictionary3, dictionary4 = {}, {}, {}, {}
        dictionary['all'], dictionary2['all'], dictionary3['all'], dictionary4['all'] = 'all', 'all', 'all', 'all'
        
        #s = 11813
        for s in meta_data['ii']:
        #if 1<2:
            """
            query = "select si.ii,si.lat,si.lon,si.elev, o.* from obs.obs_nmsko_1h o join si.si on si.id=o.si_id where ii={0} and date between '{1}' and '{2}'".format(s, start_date, end_date)
            obs_data = obs.query(query)
            obs_data = obs_data.set_index('date')
            #data_conc = data_conc.reindex(pd.date_range("2015-01-01 00:00:00", koniec , freq='H'))
            obs_data = obs_data.reindex(pd.date_range(start_date, end_date, freq = 'H'))
            obs_all[s] = obs_data[v][(n_start-1)*24:n_end*24]
            """
            
            ref_coord_sta = np.array(meta_data['ref_coord'].loc[meta_data['ii']==s].item())
            ref_data = ref_model[:,ref_coord_sta[0],ref_coord_sta[1]].values*var_const[ref_name][v]            #mod_all[s] = obs_data[v]
            ref_all[s] = ref_data
            
            other_coord_sta = np.array(meta_data['other_coord'].loc[meta_data['ii']==s].item())
            other_data = other_model[:,other_coord_sta[0],other_coord_sta[1]].values*var_const[other_name][v]            #mod_all[s] = obs_data[v]
            other_all[s] = other_data
            
            dictionary[s] = '{0}:{1}'.format(s,meta_data['name'].loc[meta_data['ii']==s].values[0])
            dictionary2[s] = meta_data['name'].loc[meta_data['ii']==s].values[0]
            dictionary3[s] = meta_data['type'].loc[meta_data['ii']==s].values[0]
            dictionary4[s] = meta_data['loc'].loc[meta_data['ii']==s].values[0]
        
        date_index = pd.date_range(start_date, end_date,freq='H')
        ref_all.index = date_index
        other_all.index = date_index
        
        ######### STATISTICS ##########
        stat = pd.DataFrame(index = dictionary.keys(), columns = ['name','loc','type','ref_var_coef','other_var_coef','bias','mge','rmse','r'])
        stat['name'] = dictionary2.values()
        stat['type'] = dictionary3.values()
        stat['loc'] = dictionary4.values()
        
        for s in stat.index:
            if s == 'all': 
                #nan_len = obs_all.isna().sum().sum()
                #all_len = mod_all.count().sum()
                #stat['coverage'][s] = 100*(all_len - nan_len)/all_len
                diff = (other_all - ref_all).to_numpy().flatten() 
                #diff = diff[~np.isnan(diff)]
                corr_other = other_all.stack()
                corr_ref = ref_all.stack()
                stat['r'][s] = corr_other.corr(corr_ref)
                
                ref_std_all, other_std_all = np.std(corr_ref), np.std(corr_other)
                ref_mean_all, other_mean_all = np.mean(corr_ref), np.mean(corr_other)
                stat['ref_var_coef'][s] = ref_std_all/ref_mean_all
                stat['other_var_coef'][s] = other_std_all/other_mean_all
        
            else:
                #nan_len = obs_all[s].isna().sum()
                #all_len = len(mod_all[s])
                #stat['coverage'][s] = 100*(all_len - nan_len)/all_len
                diff = np.array((other_all[s] - ref_all[s]))    
                #R
                stat['r'][s] = other_all[s].corr(ref_all[s])
                
                ref_std, other_std = np.std(ref_all[s]), np.std(other_all[s])
                ref_mean, other_mean = np.mean(ref_all[s]), np.mean(other_all[s])
                stat['ref_var_coef'][s] = ref_std/ref_mean
                stat['other_var_coef'][s] = other_std/other_mean
        
            if len(diff) != 0:
                #BIAS
                stat['bias'][s] = np.sum(diff)/len(diff)
                #MGE
                stat['mge'][s] = np.sum(np.absolute(diff))/len(diff)
                #RMSE
                stat['rmse'][s] = np.sqrt(np.sum(diff**2)/len(diff))
                

        print('Comparison of {0} done in {1} seconds'.format(v,time.time()-start_time))
                
    ###### SAVE DATA AS NETCDF4                                                
        val_data = xr.Dataset()            
        # data variables
        #val_data['OBS'] = (('time','station'), np.concatenate((obs_data[v].to_numpy(),np.array([s]))), axis = 1)
        val_data['Reference'.format(v)] = (('time','station'), ref_all)
        val_data['Other'.format(v)] = (('time','station'), other_all)
           
        # coordinates
        val_data.coords['time'] = times
        val_data.coords['station'] = ref_all.columns.values
           
        # attributes
        val_data.attrs['Description'] = 'Comparison of the reference run and other model run at NMSKO stations.  For each station, reference and other model data for the selected period are stored as 1D arrays.'
        val_data.attrs['pollutants'] = v
        val_data.attrs['units'] = units
        val_data.attrs['ref. domain'] = ref_grid_name
        val_data.attrs['other domain'] = other_grid_name
        val_data.attrs['ref. model name'] = ref_name
        val_data.attrs['other domain name'] = other_name
        val_data.attrs['Included station types'] = station_type
        val_data.attrs['Period'] = '{0}- {1}'.format(start_date[0:11], end_date[0:11])
        val_data.attrs['Stations'] = list(dictionary.values())
        val_data.attrs['Get stations as dictionary']: "stations_list = xr_data.attrs['Stations']; stations = {}; for i in range(0,len(stations_list)): stations[stations_list[i][0:5]]=stations_list[i][6:]"

    return stat, val_data 

#def var_coef(series): 
#computes variation coefficient of a 1D array  

#AKO TAHAT DATA 
#VYBER KONKRETNU STANICU:
#val_data['MOD'].sel(station = 11813)
#DICTIONARY SO STANICAMI
#stations_list = xr_data.attrs['Stations']
#stations = {}
#for i in range(0,len(stations_list)):
#    stations[stations_list[i][0:5]]=stations_list[i][6:]



####### PLOTS
#start_time = time.time()





















