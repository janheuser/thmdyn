##############################################################################################################################
# Date: 2022-8-10
# Name: calctsi.py
# Author: James Anheuser
# Description: Reads in AMSR-E/AMSR2 radiances and sea ice concentration, calculates snow--ice interface temperature
# Inputs: AMSR-E/AMSR2 radiance/sea ice concentration file location
# Outputs: AMSR-E/AMSR2 radiance/sea ice concentration, snow--ice interface temperature

import numpy as np
import xarray as xr
import glob
import scipy.spatial as ss
import h5py
from pyhdf import SD
from datetime import date, timedelta
from pyproj import Transformer

transformer = Transformer.from_crs(4326, 6931)

amsrfp='/ships19/cryo/janheuser/amsr/'


def readamsr(day, chan):
    """Read in AMSR radiances
    
    Arguments:
    day: requested day in datetime date format
    chan: requested channel in string format (e.g., '06V')
    
    Returns:
    radiances in an xarray data structure
    """

    if day < date(2011, 10, 5):
        file_name = glob.glob(amsrfp + 'AMSR_E_L3_SeaIce25km*' + day.strftime("%Y%m%d") + '.hdf')

        while not file_name:
            day = day + timedelta(days=1)
            file_name = glob.glob(amsrfp + 'AMSR_E_L3_SeaIce25km*' + day.strftime("%Y%m%d") + '.hdf')
            # print('AMSR data doesn\'t exist for this day, adding day')
            if day > date(2011, 10, 4):
                print('Date is within AMSRE and AMSR2 data gap')
                return

        file = SD.SD(file_name[0], SD.SDC.READ)

        ds_name = 'SI_25km_NH_' + chan + '_DAY'
        data_obj = file.select(ds_name)  # select sds

        data = data_obj.get()  # get sds data

    elif day > date(2012, 7, 1):
        file_name = glob.glob(amsrfp + 'AMSR_U2_L3_SeaIce25km*' + day.strftime("%Y%m%d") + '.he5')

        while not file_name:
            day = day + timedelta(days=1)
            file_name = glob.glob(amsrfp + 'AMSR_U2_L3_SeaIce25km*' + day.strftime("%Y%m%d") + '.he5')
            print('AMSR data doesn\'t exist for this day, adding day')

        with h5py.File(file_name[0], 'r') as f:
            data = f['/HDFEOS/GRIDS/NpPolarGrid25km/Data Fields/SI_25km_NH_' + chan + '_DAY'][()]

    else:
        print('Date is within AMSRE and AMSR2 data gap')
        return

    lat = np.load(amsrfp + 'amsr_25km_lat.npy')
    lon = np.load(amsrfp + 'amsr_25km_lon.npy')

    return xr.Dataset({'dat': (['y', 'x'], data)}, coords={'lat': (['y', 'x'], lat), 'lon': (['y', 'x'], lon)})



def readamsr_sic(day):
    """Read in AMSR sea ice concentration
    
    Arguments:
    day: requested day in datetime date format
    chan: requested channel in string format (e.g., '06V')
    
    Returns:
    sea ice concentration in an xarray data structure
    """
        
    if day < date(2011, 10, 5):
        file_name = glob.glob(amsrfp + 'AMSR_E_L3_SeaIce25km*' + day.strftime("%Y%m%d") + '.hdf')

        while not file_name:
            day = day + timedelta(days=1)
            file_name = glob.glob(
                amsrfp + 'AMSR_E_L3_SeaIce25km*' + day.strftime("%Y%m%d") + '.hdf')
            print('AMSR data doesn\'t exist for this day, adding day')
            if day > date(2011, 10, 4):
                print('Date is within AMSRE and AMSR2 data gap')
                return

        file = SD.SD(file_name[0], SD.SDC.READ)

        ds_name = 'SI_25km_NH_ICECON_DAY'
        data_obj = file.select(ds_name)  # select sds

        data = data_obj.get()  # get sds data

    elif day > date(2012, 7, 1):
        file_name = glob.glob(amsrfp + 'AMSR_U2_L3_SeaIce25km*' + day.strftime("%Y%m%d") + '.he5')

        while not file_name:
            day = day + timedelta(days=1)
            file_name = glob.glob(amsrfp + 'AMSR_U2_L3_SeaIce25km*' + day.strftime("%Y%m%d") + '.he5')
            print('AMSR data doesn\'t exist for this day, adding day')

        with h5py.File(file_name[0], 'r') as f:
            data = f['/HDFEOS/GRIDS/NpPolarGrid25km/Data Fields/SI_25km_NH_ICECON_DAY'][()]

    else:
        print('Date is within AMSRE and AMSR2 data gap')
        return

    lat = np.load(amsrfp + 'amsr_25km_lat.npy')
    lon = np.load(amsrfp + 'amsr_25km_lon.npy')

    return xr.Dataset({'dat': (['y', 'x'], data)}, coords={'lat': (['y', 'x'], lat), 'lon': (['y', 'x'], lon)})


def calcTsi_amsr_loc(day, lati, loni):
    print('test4')
    """Calculate snow--ice interface temperature at single lat/lon point via Kilic et al. (2021)
    
    Arguments:
    day: requested day in datetime date format
    lati: latitude of requested point
    loni: longitude of requested point
    
    Returns:
    snow--ice interface temperature at requested point in float format
    """
    
    v6 = readamsr(day, '06V')/10
    v18 = readamsr(day, '18V')/10
    v36 = readamsr(day, '36V')/10
    
    Ds=1.7701+0.0175*v6-0.028*v18+0.0041*v36
    Ds=Ds.where(Ds!=1.7701)
    Ds=Ds.where(Ds>0)
    
    Tsi=1.086*v6+3.98*np.log(Ds)-10.7
    Tsi=Tsi.where(Tsi.dat>0)
        
    lon=Tsi.lon.where(Tsi.lon>0,Tsi.lon+360)
    Tsi=Tsi.assign_coords({'lon': lon})
    
    transformer = Transformer.from_crs(4326, 6931)
    new_x, new_y = transformer.transform(Tsi.lat.values.ravel(), Tsi.lon.values.ravel())
    xi, yi = transformer.transform(lati, loni)
    
    kdtree = ss.KDTree(np.stack((new_x, new_y), axis=1))
    ind = kdtree.query([xi, yi])
    
    return Tsi.dat.values.ravel()[ind[1]]

def calcTsi_amsr(day):
    """Calculate snow--ice interface temperature via Kilic et al. (2021)
    
    Arguments:
    day: requested day in datetime date format
    
    Returns:
    snow--ice interface temperature in an xarray data structure
    """
    
    v6 = readamsr(day, '06V')/10
    v18 = readamsr(day, '18V')/10
    v36 = readamsr(day, '36V')/10
    
    sic = readamsr_sic(day)
    
    v6=v6.where(sic.dat != 120)
    v6=v6.where(sic.dat > 95)
    
    v18=v18.where(sic.dat != 120)
    v18=v18.where(sic.dat > 95)
    
    v36=v36.where(sic.dat != 120)
    v36=v36.where(sic.dat > 95)
    
    Ds=1.7701+0.0175*v6-0.028*v18+0.0041*v36
    Ds=Ds.where(Ds!=1.7701)
    Ds=Ds.where(Ds>0)
    
    Tsi=1.086*v6+3.98*np.log(Ds)-10.7
    Tsi=Tsi.where(Tsi.dat>0)
    
    lon=Tsi.lon.where(Tsi.lon>0,Tsi.lon+360)
    Tsi=Tsi.assign_coords({'lon': lon})
    
    return Tsi