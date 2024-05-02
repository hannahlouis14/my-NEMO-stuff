import numpy as np
import xarray as xr
import netCDF4 as nc
import datetime as dt
import glob
import cmocean
import sys
from assign_timestamp import assign_timestamp
import matplotlib.pyplot as plt
import pandas as pd
from daterange import daterange
import scipy.signal as sp


time = daterange()
filtered_list = [date for date in time if date.year == 2002]
#print(filtered_list[2])



def lag_linregress_3D(x, y, lagx=0, lagy=0):
    """
    Input: Two xr.Datarrays of any dimensions with the first dim being time. 
    Thus the input data could be a 1D time series, or for example, have three 
    dimensions (time,lat,lon). 
    Datasets can be provided in any order, but note that the regression slope 
    and intercept will be calculated for y with respect to x.
    Output: Covariance, correlation, regression slope and intercept, p-value, 
    and standard error on regression between the two datasets along their 
    aligned time dimension.  
    Lag values can be assigned to either of the data, with lagx shifting x, and
    lagy shifting y, with the specified lag amount. 
    """ 
    #1. Ensure that the data are properly alinged to each other. 
    x,y = xr.align(x,y)

    #2. Add lag information if any, and shift the data accordingly
    if lagx!=0:

        # If x lags y by 1, x must be shifted 1 step backwards. 
        # But as the 'zero-th' value is nonexistant, xr assigns it as invalid 
        # (nan). Hence it needs to be dropped
        x   = x.shift(time = -lagx).dropna(dim='time')

        # Next important step is to re-align the two datasets so that y adjusts
        # to the changed coordinates of x
        x,y = xr.align(x,y)

    if lagy!=0:
        y   = y.shift(time = -lagy).dropna(dim='time')
        x,y = xr.align(x,y)

    #3. Compute data length, mean and standard deviation along time axis: 
    n = y.notnull().sum(dim='time')
    xmean = x.mean(axis=0)
    ymean = y.mean(axis=0)
    xstd  = x.std(axis=0)
    ystd  = y.std(axis=0)

    #4. Compute covariance along time axis
    cov   =  np.sum((x - xmean)*(y - ymean), axis=0)/(n)

    #5. Compute correlation along time axis
    cor   = cov/(xstd*ystd)

    #6. Compute regression slope and intercept:
    slope     = cov/(xstd**2)
    intercept = ymean - xmean*slope  

    #7. Compute P-value and standard error
    #Compute t-statistics
    tstats = cor*np.sqrt(n-2)/np.sqrt(1-cor**2)
    stderr = slope/tstats

    from scipy.stats import t
    pval   = t.sf(tstats, n-2)*2
    pval   = xr.DataArray(pval, dims=cor.dims, coords=cor.coords)

    return cov,cor,slope,intercept,pval,stderr

def preprocess(ds,var='ileadfra'):
    return ds[var].to_dataset()



def correlate(runid, turb_var, ice_var, two_regions, year):
    data_path = '/project/6007519/pmyers/ANHA4/ANHA4-'+runid+'-S/'
    air_data_path = '/project/6007519/hlouis/data_files/AtmForcing/'
    turb_data_path = '/project/6007519/hlouis/data_files/TurbFluxes/'+runid+'/'

    if turb_var=='latent':
        turb_files = sorted(glob.glob(turb_data_path+'EPM151_TurbLatent_no-ice_y'+str(year)+'.nc'))

    if turb_var=='sensible':
        turb_files = sorted(glob.glob(turb_data_path+'EPM151_TurbSensible_no-ice_y'+str(year)+'.nc'))

    #ice_files =  sorted(glob.glob(data_path+'*y'+str(year)+'*icemod.nc')))

    t=73
    if two_regions==True:
        shb_mask_path = '/project/6007519/hlouis/scripts/masks/shb_mask.nc'
        jb_mask_path = '/project/6007519/hlouis/scripts/masks/jb_mask.nc'
        mf1 = nc.Dataset(shb_mask_path)
        hb_mask = mf1['tmask']

        mf2 = nc.Dataset(jb_mask_path)
        jb_mask = mf2['tmask']
    
    else:
        mask_path = '/project/6007519/hlouis/scripts/masks/shbjb_mask.nc'
        mf = nc.Dataset(mask_path)
        mask = mf['tmask']
        mask = np.broadcast_to(mask,(t,)+mask.shape)
    '''
    # load the turbulent flux and ice variable
    d = xr.open_mfdataset(turb_files, preprocess=lambda ds: ds[[turb_var]], concat_dim='time_counter', combine='nested', data_vars='minimal', coords='minimal', compat='override', chunks={'time_counter': 1})
    ds = xr.open_mfdataset(ice_files, preprocess=preprocess, concat_dim='time_counter', combine='nested', data_vars='minimal', coords='minimal', compat='override')
     
    turb_flux = d[turb_var]
    ice = ds['ileadfra']
    turb_flux = turb_flux.where(mask==1) 
    ice = ice.where(mask==1) 
    
    # first, take regional nan mean of ice and turb flux
    turb_flux = turb_flux.mean(dim=['y','x'], skipna=True) 
    ice = ice.mean(dim=['y','x'], skipna=True)
    
    N = ice['time_counter'].size
    lead = 30
    
    # now, shift the turb flux data by some time window
    normal = ice #.isel({'time_counter': slice(0, N - lead)})
    #shifted = turb_flux.isel({'time_counter': slice(0 + lead, N)})
    shifted  = turb_flux.shift(time_counter=-30).dropna(dim='time_counter')
    
    #normal, shifted = xr.align(normal, shifted, join='right')
    
    ice = ice[:-30]
    print(ice)
    
    shifted['time_counter'] = ice['time_counter']
    print(shifted) 
    # now, you have to align them again along the time dimension
    #shifted['time_counter'] = normal['time_counter']
    #ice, turb_flux = xr.align(ice, turb_flux)
     
    # apply the correlation now
    #cor = xr.corr(ice, shifted, dim=('time_counter')) 
    ''' 

    coefficients = [] 
    #corrcoeff = np.corrcoef(ice[:-30].values, shifted.values)
    ''' 
    for i in range(N-lead):
        corrcoeff = xr.corr(ice[i], shifted[i])#, dim='time_counter')
        print(corrcoeff.values)
        coefficients.append(corrcoeff)
    '''
    
    for ifile, file in enumerate(turb_files):
        year = file.split("_y")[1][:4]
         
        d = xr.open_dataset(file)
        turb_flux = d[turb_var]
        turb_flux = turb_flux.where(mask==1)
        print('before', turb_flux)
        turb_flux = turb_flux[25:]#.shift(time_counter=-30).dropna('time_counter')
        
        print('after', turb_flux)
        print('after', turb_flux['time_counter'])
        
        ice_files = sorted(glob.glob(data_path+'*'+str(year)+'*icemod.nc'))
        ice = xr.open_mfdataset(ice_files, concat_dim='time_counter', combine='nested', data_vars='minimal', coords='minimal', compat='override')
        
        ice = ice[ice_var] 
        ice = ice.where(mask==1)
        ice = ice[:-25]
        print(ice)
        print(turb_flux)
        turb_flux['time_counter'] = ice['time_counter']
        print('AFTER', turb_flux['time_counter'])
        #ice = assign_timestamp(d, ice)
        #t = pd.to_datetime(ice['time_counter'].values)
        
        times = make_dates(endyear=2002, endmonth=12, endday=31, startyear=2002, startmonth=1, startday=5)
        cor = xr.corr(turb_flux, ice, dim=('y','x'))
        print(cor)
        plt.plot(times[:-25], cor)
        plt.ylabel('pearson r')
        plt.xlabel('date')
        plt.show()
        print('before nanmean', cor)
        cor = np.nanmean(cor)
        print('after nanmean', cor)
        coefficients.append(cor)

    #return coefficients
   



def make_dates(endyear, endmonth, endday, startyear, startmonth, startday):
    start_time = dt.date(startyear, startmonth, startday)

    end_time = dt.date(endyear, endmonth, endday)

    #figure out all the dates we have model files
    delta = end_time - start_time
    times = []

    i = 0
    while i < delta.days+1:
        t = start_time + dt.timedelta(days=i)
        if t.month == 2 and t.day == 29:
            t = dt.date(t.year, 3, 1)
            i = i+6
        else:
            i = i+5
        times.append(t)

    return times



def plot_corr(pearson_r):
    #years = pd.date_range(start=pd.datetime(2002, 1, 1), periods=21, freq='AS')# end='2022-12-31')
    #years = np.arange(2002,2023,1)
    years = make_dates(endyear=2022, endmonth=12, endday=31, startyear=2002, startmonth=1, startday=5)
    years = years[:-30]

    plt.plot(years, pearson_r, color='purple')
    plt.ylabel("pearson's r")
    plt.title('Correlation between Ice Concentration and Sensible Heat Flux')
    #plt.show()
    plt.savefig('IC-SenseHF.png')


#time = make_dates(endyear=2022, endmonth=12, endday=31, startyear=2002, startmonth=1, startday=5)


coeffs = correlate('EPM151', 'sensible', 'ileadfra', two_regions=False, year=2002)
#plot_corr(coeffs)
'''
sense_contra = correlate('EPM151', 'sensible', 'ileadfra', two_regions=False)
sense_thicc = correlate('EPM151', 'sensible', 'iicethic', two_regions=False)
latent_contra= correlate('EPM151', 'latent', 'ileadfra', two_regions=False)
latent_thicc = correlate('EPM151', 'latent', 'iicethic', two_regions=False)
print(np.nanmean(sense_contra))
print(np.nanmean(sense_thicc))
print(np.nanmean(latent_contra))
print(np.nanmean(latent_thicc))
years = np.arange(2002,2023,1)
plt.plot(years, sense_contra, label='Q$_{SH}$-ice conc')
plt.plot(years, sense_thicc, label='Q$_{SH}$-ice thic')
plt.plot(years, latent_contra, label='Q$_{LH}$-ice conc')
plt.plot(years, latent_thicc, label='Q$_{LH}$-ice thic')
'''
