import datetime
import xarray as xr
import netCDF4 as nc
import numpy.ma as ma
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob
import warnings



def t_profile(runid, section, endyear, endmonth, endday, startyear=2002, startmonth=1, startday=5):
    path = '/project/6007519/pmyers/ANHA4/ANHA4-'+runid+'-S/'
    output_path = '/project/6007519/hlouis/plotting/figures/ts/'
    mask_file = '/project/6007519/hlouis/scripts/HBC_mask.nc'
    
    #mdl_files = glob.glob(path+'ANHA4-'+runid+'*_gridT.nc')
    start_time = datetime.date(startyear, startmonth, startday)
    end_time = datetime.date(endyear, endmonth, endday)
 
    delta = end_time - start_time
    times = []

    i = 0
    while i < delta.days+1:
        t = start_time + datetime.timedelta(days=i)
        if t.month == 2 and t.day == 29:
            t = datetime.date(t.year, 3, 1)
            i = i+6
        else:
            i = i+5
        times.append(t)
    
    mdl_files = []
    
    for t in times:
        mdl_files.append(path+"ANHA4-"+runid+"_y"+str(t.year)+"m"+str(t.month).zfill(2)+"d"+str(t.day).zfill(2)+"_gridT.nc")
    
    #months_out = ['m01', 'm02', 'm03', 'm04', 'm05', 'm09', 'm10', 'm11', 'm12']  # for summer months
    months_out = ['m04', 'm05', 'm06', 'm07', 'm08', 'm09', 'm10', 'm11', 'm12']  # for winter months
    mdl_files = [x for x in mdl_files if not any(word in x for word in months_out)]
    
    nfiles = len(mdl_files)
    
    years = set([f.split("_y")[1][:4] for f in mdl_files])
    years = sorted(years)
    
    d = xr.open_dataset(mdl_files[0])
    ndepth = d.deptht.shape[0]

    mf = nc.Dataset(mask_file)
    rmask = mf[section][0]

    yearly_temp = {}

    for year in years:
        print(year)
        yearly_temp.update({year: [np.zeros(ndepth), 0]})

    for ifile, file in enumerate(mdl_files):
        print('reading', file, str(ifile)+'/'+str(nfiles))
        year = file.split("_y")[1][:4]
        d = xr.open_dataset(file)
        temp = d['vosaline'].values

        masked_temp = np.where(rmask==2,temp,np.nan)
        with warnings.catch_warnings():
            warnings.filterwarnings(action='ignore', message='Mean of empty slice')
            try:
                regional_temp_ts = np.nanmean(masked_temp, axis=(2,3)).flatten()
            except RuntimeWarning:
                regional_temp_ts = np.NaN
        yearly_temp[year][0] += regional_temp_ts  # running sum
        yearly_temp[year][1] +=1  # counter
    
    for year in years:
        yearly_temp[year][0] /= yearly_temp[year][1]  # compute average for each year

    fig = plt.figure(figsize=(6,5))
    nyears = len(years)
    print(nyears)
    
    for i, year in enumerate(years):
        #ax = fig.add_subplot(3,int(nyears/3), i+1, sharex=True)
        plt.plot(yearly_temp[year][0], d.deptht, 'b-')
        plt.xlabel('Salinity (psu)', fontsize=12)
        plt.ylabel('depth', fontsize=12)
        plt.ylim(270, -5)
        #ax.label_outer()
        #ax.set_title(year)
        plt.title(runid+' Winter Salinity Profiles')  #year)
    plt.savefig(output_path+'EPM151_winter_sal_profile.png')
    

t_profile(runid='EPM151', section='hbc', endyear=2022, endmonth=12, endday=31)
#t_profile(runid='EPM111', section='hbc', endyear=2009, endmonth=12, endday=31)  

