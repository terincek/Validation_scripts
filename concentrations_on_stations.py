#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script na vytiahnutie merani zo stanic vo formate netcdf4 ako pouzivame pri validacii
Pozor, casova os je tu zachovana ako na staniciach - napr. 10 hodin znamena priemer za 9-10, v modeli je 10 priemer za 10-11, 
vo validacii je to posunute na podelovy cas
Created on Tue Dec  7 09:27:59 2021

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

sys.path.append(vali_path)
sys.path.append('/home/KOL/p6065/python-scripts/core')
import day_of_year as doy

from dbConnector import obs

query_meta = obs.q_nmsko 
meta_data_ = obs.query(query_meta)


var = ['SO2','NO2','PM10','PM2_5','O3','CO'] 

units = 'ug/m3'
station_type = ['B','T','I']

station_blacklist = [99701,99702,99703,99501,99502,99503,99350,99223,99130,99126]
meta_data = meta_data_[~meta_data_.ii.isin(station_blacklist)].copy()

def pozorovania(start_date,end_date,v):
    
    obs_all = pd.DataFrame()
    dictionary, dictionary2, dictionary3, dictionary4 = {}, {}, {}, {}
    dictionary['all'], dictionary2['all'], dictionary3['all'], dictionary4['all'] = 'all', 'all', 'all', 'all'
        
    n_start = doy.day_of_year(start_date)
    n_end = doy.day_of_year(end_date)
    n_days = n_end - n_start + 1
    #times = pd.date_range(start=start_date, periods=(n_days-1)*24 + 1, freq='H')
    times = pd.date_range(start=start_date, periods=(n_days)*24, freq='H')

    
    for s in meta_data['ii']:
        query = "select si.ii,si.lat,si.lon,si.elev, o.* from obs.obs_nmsko_1h o join si.si on si.id=o.si_id where ii={0} and date between '{1}' and '{2}'".format(s, start_date, end_date)
        obs_data = obs.query(query)
        obs_data = obs_data.set_index('date')
        obs_data = obs_data.reindex(pd.date_range(start_date, end_date, freq = 'H'))
        #obs_all[s] = obs_data[v][(n_start-1)*24:n_end*24]
        obs_all[s] = obs_data[v][:]
        
        dictionary[s] = '{0}:{1}'.format(s,meta_data['name'].loc[meta_data['ii']==s].values[0])
        dictionary2[s] = meta_data['name'].loc[meta_data['ii']==s].values[0]
        dictionary3[s] = meta_data['type'].loc[meta_data['ii']==s].values[0]
        dictionary4[s] = meta_data['loc'].loc[meta_data['ii']==s].values[0]
    
    # STATS
    stat = pd.DataFrame(index = dictionary.keys(), columns = ['name','loc','type','coverage','obs_mean'])
    stat['name'] = dictionary2.values()
    stat['type'] = dictionary3.values()
    stat['loc'] = dictionary4.values()    
    obs_all.fillna(value=np.nan, inplace=True)
    
    for s in stat.index:
        print(s)
        if s == 'all': 
            obs_all[obs_all == 0] = np.nan
            nan_len = obs_all.dropna(axis = 1, how = 'all').isna().sum().sum()
            #all_len = obs_all.dropna(axis = 1, how = 'all').count().sum()
            all_len = len(obs_all)*len(obs_all.dropna(axis = 1, how = 'all').columns)
            stat['coverage'][s] = 100*(all_len - nan_len)/float(all_len)
            stat['obs_mean'][s] = np.mean(obs_all).mean()
            
        else:
            obs_all[s][obs_all[s] == 0] = np.nan
            nan_len = obs_all[s].isna().sum()
            all_len = len(obs_all[s])   # TOTO SKONTROLOVAT CI JE DOBRA HODNOTA
            stat['coverage'][s] = 100*(all_len - nan_len)/float(all_len)
            stat['obs_mean'][s] = np.mean(obs_all[s])
            
    #print(obs_all)
    #return stat, obs_all

    # PUTTING INTO DATASET FORMAT
    val_data = xr.Dataset()            
    # data variables
    #val_data['OBS'] = (('time','station'), np.concatenate((obs_data[v].to_numpy(),np.array([s]))), axis = 1)
    val_data['OBS'.format(v)] = (('time','station'), obs_all)
       
    # coordinates
    val_data.coords['time'] = times
    val_data.coords['station'] = obs_all.columns.values
       
    # attributes
    val_data.attrs['Description'] = 'Concentrations of {0} at NMSKO stations.'.format(v)
    val_data.attrs['pollutants'] = v
    val_data.attrs['units'] = units
    val_data.attrs['Included station types'] = station_type
    val_data.attrs['Period'] = '{0}- {1}'.format(start_date[0:11], end_date[0:11])
    val_data.attrs['Stations'] = list(dictionary.values())
    val_data.attrs['Get stations as dictionary']: "stations_list = xr_data.attrs['Stations']; stations = {}; for i in range(0,len(stations_list)): stations[stations_list[i][0:5]]=stations_list[i][6:]"

    return stat, val_data
    
    










