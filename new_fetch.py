# -*- coding: utf-8 -*-
"""
Created on Tue Dec  1 16:50:05 2020


variables
1)	'integral_of_surface_downwelling_shortwave_flux_in_air_wrt_time',
2)	'air_temperature_2m',
3)	'relative_humidity_2m',
4)	'x_wind_10m',
5)	'y_wind_10m',
6)	'precipitation_amount_acc'



Three cases are presented in the archive data:
1)	Basically, From 2020-02-04 on, the archive directory to find all required variables is,
for example, "https://thredds.met.no/thredds/dodsC/meps25epsarchive/2020/08/29/meps_det_2_5km_20200829T18Z.nc"
with every 3 hours for updating data (00,03,06,09,12,15,18,21).
2)	From 2019-06-01 to 2020-02-03, the directory is, for example,
“https://thredds.met.no/thredds/dodsC/meps25epsarchive/2019/07/19/meps_extracted_2_5km_20190719T03Z.nc.html"
with every 3 hours for updating data (00,03,06,09,12,15,18,21).
3)	Before 2019-05-31 the directory is, for example,
“https://thredds.met.no/thredds/dodsC/meps25epsarchive/2019/02/18/meps_subset_2_5km_20190218T18Z.nc”
with every 6 hours for updating data (00, 06,12,18).



filename = "https://thredds.met.no/thredds/dodsC/meps25epsarchive/2020/08/29/meps_det_2_5km_20200829T18Z.nc"
filename = "https://thredds.met.no/thredds/dodsC/meps25epsarchive/2019/07/19/meps_extracted_2_5km_20190719T03Z.nc"
filename= "https://thredds.met.no/thredds/dodsC/meps25epsarchive/2019/02/18/meps_subset_2_5km_20190218T18Z.nc"
filename = "https://thredds.met.no/thredds/dodsC/meps25epsarchive/YYYY/MM/DD/meps_det_2_5km_YYYYMMDDTHHZ.nc"

@author: alkragh and Li
"""
import netCDF4
import datetime
import numpy as np
import pandas as pd
import time
import multiprocessing

import os


from datetime import date, timedelta

def get_first_day(dt, d_years=0, d_months=0):
    # d_years, d_months are "deltas" to apply to dt
    y, m = dt.year + d_years, dt.month + d_months
    a, m = divmod(m-1, 12)
    d_first=date(y+a, m+1, 1)
    d_first_stmp=datetime.datetime(d_first.year,d_first.month,d_first.day,0,0,0,tzinfo=datetime.timezone.utc)
    return d_first_stmp

def get_last_day(dt):
    d_last=get_first_day(dt, 0, 1) + timedelta(-1)
    d_last_stmp=datetime.datetime(d_last.year,d_last.month,d_last.day,18,0,0,tzinfo=datetime.timezone.utc)
    return d_last_stmp

def gene_start_stop(year):
    dates_start=[datetime.datetime(year,month,1,0,0,0,tzinfo=datetime.timezone.utc) for month in range(1,13)]
    dates_stop=[get_last_day(date) for date in dates_start]

    stmps_start=[datetime.datetime.timestamp(date) for date in dates_start]
    stmps_stop=[datetime.datetime.timestamp(date) for date in dates_stop]

    return stmps_start, stmps_stop

def getclosest_ij(lats, lons, latpt, lonpt):
    # find squared distance of every point on grid
    dist_sq = (lats - latpt) ** 2 + (lons - lonpt) ** 2
    # 1D index of minimum dist_sq element
    minindex_flattened = dist_sq.argmin()
    # Get 2D index for latvals and lonvals arrays from 1D index
    return np.unravel_index(minindex_flattened, lats.shape)


def f_lan_long_mins(f):
    # for var in variables:
    lat, lon = f.variables['latitude'], f.variables['longitude']
    # extract lat/lon values (in degrees) to numpy arrays
    latvals = lat[:]
    lonvals = lon[:]
    # Nordhavn55.715833°, 12.606667  #Bornholm 55.13, 14.91
    iy_min, ix_min = getclosest_ij(latvals, lonvals, 55.716, 12.607)

    return iy_min, ix_min


def fetch_one_month(start, stop):

    filename1 = 'https://thredds.met.no/thredds/dodsC/meps25epsarchive/'
    filename2= '/meps_subset_2_5km_'
    filename3=  '/meps_extracted_2_5km_'
    filename4 = '/meps_det_2_5km_'
    variables = ['integral_of_surface_downwelling_shortwave_flux_in_air_wrt_time',
                 'air_temperature_2m', 'relative_humidity_2m',
                 'precipitation_amount_acc','x_wind_10m',
                 'y_wind_10m']


    step=21600
    secs = start
    rows = int((stop - start) / step + 1)
    cols = 49
    res = np.full((len(variables), rows, cols), -99.99)

    col_names = ['hour-ahead-' + str(k) for k in range(cols)]
    index = [(start + k * step) for k in range(rows)]

    # iy_min, ix_min = 168, 422

    i = 0
    while secs <= stop:
        print(secs)
        print(datetime.datetime.utcfromtimestamp(secs).strftime('%Y%m%dT%HZ'))
        YYYY = datetime.datetime.utcfromtimestamp(secs).strftime('%Y')
        MM = datetime.datetime.utcfromtimestamp(secs).strftime('%m')
        DD = datetime.datetime.utcfromtimestamp(secs).strftime('%d')
        HH = datetime.datetime.utcfromtimestamp(secs).strftime('%H')
        file1 = filename1 + YYYY + '/' + MM + '/' + DD + filename2 + YYYY + MM + DD + 'T' + HH + 'Z.nc'
        file2 = filename1 + YYYY + '/' + MM + '/' + DD + filename3 + YYYY + MM + DD + 'T' + HH + 'Z.nc'
        file3 = filename1 + YYYY + '/' + MM + '/' + DD + filename4 + YYYY + MM + DD + 'T' + HH + 'Z.nc'

        try:
            f = netCDF4.Dataset(file1, mode="r")
            iy_min, ix_min = f_lan_long_mins(f)
            # print(iy_min, ix_min)
            # print("I am at the case 1")
            for var_idx, var in enumerate(variables):
                temp = f.variables[var]
                loop = range(0, cols)
                for j in loop:
                    if (temp[j, 0, 0, iy_min, ix_min].mask==False):
                        res[variables.index(var), i, j] = (temp[j, 0, 0, iy_min, ix_min])

            i = i + 1
            secs = start + i * step

            # print(i, j)
        except:
            try:
                f = netCDF4.Dataset(file2, mode="r")
                iy_min, ix_min = f_lan_long_mins(f)
                # print(iy_min, ix_min)
                # print("I am at the case 2")
                for var_idx, var in enumerate(variables):
                    temp = f.variables[var]
                    loop = range(0, cols)
                    for j in loop:
                        if (temp[j, 0, 0, iy_min, ix_min].mask==False):
                            res[variables.index(var), i, j] = (temp[j, 0, 0, iy_min, ix_min])
                i = i + 1
                secs = start + i * step


            except:
                try:
                    f = netCDF4.Dataset(file3, mode="r")
                    iy_min, ix_min = f_lan_long_mins(f)
                    # print(iy_min, ix_min)
                    # print("I am at the case 3")
                    for var_idx, var in enumerate(variables):
                        temp = f.variables[var]
                        loop = range(0, cols)
                        for j in loop:
                            if (temp[j, 0, iy_min, ix_min].mask ==False):
                                res[variables.index(var), i, j] = (temp[j, 0, iy_min, ix_min])
                    i = i + 1
                    secs = start + i * step

                    # print(i, j)
                except:
                    i = i + 1
                    secs = start + i * step

    YYYY = datetime.datetime.utcfromtimestamp(start).strftime('%Y')
    MM = datetime.datetime.utcfromtimestamp(start).strftime('%m')
    for var in variables:
        df = pd.DataFrame(res[variables.index(var), :, :], index, col_names)

        filename = "./Nordhavn_weather/Nordhavn_"+var +"_"+YYYY+MM +'.csv'
        df.to_csv(filename)


if __name__ == '__main__':
    start_2018, stop_2018=gene_start_stop(year=2018)
    start_2019, stop_2019=gene_start_stop(year=2019)
    start_2020, stop_2020 = gene_start_stop(year=2020)

    # start=start_2020
    # stop=stop_2020   #
    # start=start
    # stop=list(np.array(start)+21600*2)
    # stop=stop_2020[10:11]
    if not os.path.exists('./Nordhavn_weather'):
        os.makedirs('./Nordhavn_weather')

    for k in range(1):
        if k==0:
            start = start_2018
            stop = stop_2018  #

            with multiprocessing.Pool(processes=8) as pool:
                pool.starmap(fetch_one_month, zip(start, stop))

        if k==1:
            start = start_2019
            stop = stop_2019  #

            with multiprocessing.Pool(processes=8) as pool:
                pool.starmap(fetch_one_month, zip(start, stop))

        if k==2:
            start = start_2020
            stop = stop_2020  #

            with multiprocessing.Pool(processes=8) as pool:
                pool.starmap(fetch_one_month, zip(start, stop))






# now2 = datetime.datetime.now()
# print("time cost: {}".format(now2 - now1))
#
# '''
# df_integral_of_surface_parallel_solar_flux_wrt_time = pd.DataFrame(res[],index,col_names)
# df_integral_of_surface_direct_normal_irradiance_wrt_time = pd.DataFrame(res,index,col_names)
# air_temperature_2m = pd.DataFrame(res,index,col_names)
# relative_humidity_2m = pd.DataFrame(res,index,col_names)
# air_temperature_max = pd.DataFrame(res,index,col_names)
# air_temperature_min = pd.DataFrame(res,index,col_names)
# wind_direction = pd.DataFrame(res,index,col_names)
# wind_speed = pd.DataFrame(res,index,col_names)
# precipitation_amount_acc = pd.DataFrame(res,index,col_names)
#
# '''
#

# p.map(f,variables)

# datetime.datetime(2020,2,28,0,0,0, datetime.timezone.utc)
