#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Vypocita podla datumu o kolkaty den v roku ide, alebo naopak podla poradia dna v roku urci datum. 
Datum je vo formate 'YYYY-MM-DD'

Created on Tue Apr  6 13:10:13 2021

@author: ext29843
"""

def day_of_year(date):
    day = int(str(date[8:10]))
    month = int(str(date[5:7]))
    year = int(date[0:4])
    leap_year = False
    
    if (year%4) == 0:
        leap_year = True
        if (year%100) == 0:
            leap_year = False
            if (year%400) == 0:
                leap_year = True

    #months = {1:31, 2:28, 3:31, 4:30, 5:31, 6:30, 7:31, 8:31, 9:30, 10:31, 11:30, 12:31} 
    
    months = {0:0, 1:31, 2:59, 3:90, 4:120, 5:151, 6:181, 7:212, 8:243, 9:273, 10:304, 11:334, 12:365}
    
    if leap_year:
        months = {0:0, 1:31, 2:60, 3:91, 4:121, 5:152, 6:182, 7:213, 8:244, 9:274, 10:305, 11:335, 12:366}
            
    n_day = months[month-1] + day
    
    return n_day

def date_from_day(n, year):
    
    leap = False
    year = int(year)
    
    if (year%4) == 0:
        leap = True
        if (year%100) == 0:
            leap = False
            if (year%400) == 0:
                leap = True
                
    if n > 365 and leap == False:
        return 'Year {} has 365 days.'.format(year)
    
    if n > 366 and leap == True:
        return 'Year {} has 366 days.'.format(year)
    
    months = {1:31, 2:59, 3:90, 4:120, 5:151, 6:181, 7:212, 8:243, 9:273, 10:304, 11:334, 12:365}
    
    if leap:
        months = {1:31, 2:60, 3:91, 4:121, 5:152, 6:182, 7:213, 8:244, 9:274, 10:305, 11:335, 12:366}
        
    if n < 32:
        day = n
        month = 1
    else:
        for m in reversed(range(1,13)):
            if (months[m]-n) < 0:
                month = m + 1
                day = -(months[m] - n)
                break
            
            elif (months[m]-n) == 0:
                month = m 
                day = (months[m] - months[m-1])
                break
        
    if day < 10:
        day = '0{}'.format(day)
        
    if month < 10:
        month = '0{}'.format(month)
        
    date = '{0}-{1}-{2}'.format(year,month,day)

    return date
        
            





