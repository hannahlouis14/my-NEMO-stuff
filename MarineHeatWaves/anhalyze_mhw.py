import numpy as np
import netCDF4 as nc
import datetime
import pandas as pd
import xarray as xr

# Plotting-related libraries
import matplotlib.pyplot as plt
import cartopy.feature as feature
import cartopy.crs as ccrs
import matplotlib.dates as mdates

# OS-specific libraries
import glob


fig_path = '/project/6007519/hlouis/plotting/figures/mhw/'
path = '/project/6007519/hlouis/scripts/MarineHeatWaves/timeseries_data/'

Hudson_bay = False           # Boolean if using Hudson Bay vs James Bay locations.
depth = 0                    # z axis location from 50 unit "depth" 

if Hudson_bay: 
    hudson_east = -75
    hudson_west = -95
    hudson_north = 65
    hudson_south = 50
    loc = 'hudson_bay'

else:
    hudson_east = -78.5
    hudson_west = -82.5
    hudson_north = 54.7
    hudson_south = 51
    loc = 'james_bay'

lat_range = (hudson_south, hudson_north)
lon_range = (hudson_west, hudson_east)



def get_row_col_range(data, lat_range=lat_range, lon_range=lon_range, grid='gridT'):
    """ Get the row and col range given lat and lon range. 
    lat_range: latitude range (tuple)
    lon_range: londitude range (tuple)
    """

    # Get all lat-lon data
    if grid == 'gridT':
        lat = data['nav_lat_grid_T'][:]
        lon = data['nav_lon_grid_T'][:]
    else:
        # Get all lat-lon data
        lat = data['nav_lat'][:]
        lon = data['nav_lon'][:]

    # Create mask given lat lon values.
    lat_mask = np.ma.filled((lat.data > lat_range[0]) & (lat.data < lat_range[1]))
    lon_mask = np.ma.filled((lon.data > lon_range[0]) & (lon.data < lon_range[1]))

    # Apply masks to data
    mask = lat
    mask[~(lat_mask & lon_mask)] = 0

    # Find the row,col range by collapsing each axis.
    row_ranges = np.where(mask.data.sum(axis=1) > 0)[0]
    col_ranges = np.where(mask.data.sum(axis=0) > 0)[0]

    # Select range
    row_range = (row_ranges[0], row_ranges[-1])
    col_range = (col_ranges[0], col_ranges[-1])

    return row_range, col_range



def get_lat_lon(data, lat_range, lon_range, cardinal=True):
    """  Getting Latitude and Longitude """

    # Given data selection range in lat-lon or row-col
    if cardinal:
        row_range, col_range = get_row_col_range(data, lat_range, lon_range)
    else:
        row_range, col_range = lat_range, lon_range

    lat = data['nav_lat_grid_T'][row_range[0]:row_range[1], col_range[0]:col_range[1]]
    lon = data['nav_lon_grid_T'][row_range[0]:row_range[1], col_range[0]:col_range[1]]

    return lat, lon



def get_date(filename, how=''):
    """  Get date from filename.
         Assuming filename format: */*/ANHA4-EPM111_y1998m04d05_gridB.nc
         Multiple output formats possible.
    """

    # Get full date from filename
    date = filename.split('_')[-2]

    # Return format specific date info
    if how == 'ymd':
        return int(date[1:5]), int(date[6:8]), int(date[9:11])
    elif how == 'y':
        return int(date[1:5])
    elif how == 'm':
        return int(date[6:8])
    else:
        return date



def get_timeseries(runid, minmax=False):
    if runid == 'EPM151':
        data_path = '/project/6007519/pmyers/ANHA4/ANHA4-EPM151-S/'
    
    if runid == 'ETW161':
        data_path = '/project/6007519/weissgib/ANHA4/ANHA4-ETW161-S/'
    
    if runid == 'ETW162':
        data_path = '/project/6007519/weissgib/ANHA4/ANHA4-ETW162-S/' 

    if runid=='EPM111':
        data_path = '/project/6007519/pmyers/ANHA4/ANHA4-EPM111-S/'

    mask_path = '/project/6007519/hlouis/scripts/masks/'
    output_path = '/project/6007519/hlouis/scripts/MarineHeatWaves/'
    depth = 0
    ANHA4 = False
    if ANHA4:
        mask = nc.Dataset(mask_path + 'ANHA4_mesh_mask.nc')
        tmask = mask['tmask'][0][depth]
    else:
        mask = nc.Dataset(mask_path + 'shbjb_mask.nc')  
        tmask = mask['tmask']  # depth is already surface and there is no time dimension in hommade masks
    
    file_list = sorted(glob.glob(data_path+'ANHA4-'+runid+'*_gridT.nc'))
    nfiles = len(file_list)
    tmask = np.broadcast_to(tmask,(1,)+tmask.shape)

    var_means = []
    var_stds = []
    var_mins = []
    var_maxs = []
    dates = []


    for filename in file_list:  
        # get the ANHA4 data
        #data = nc.Dataset(filename)
        print(filename) 
        data = xr.open_dataset(filename)
        ''' 
        cardinal = True
        if cardinal:
            row_range, col_range = get_row_col_range(data, lat_range, lon_range)
        else:
            row_range, col_range = lat_range, lon_range
       

        surf_mask = tmask[row_range[0]:row_range[1], col_range[0]:col_range[1]]
        '''
        var = 'votemper'
        
        var_data = data[var][:,0]
        var_data = var_data.rename({'y_grid_T': 'y', 'x_grid_T': 'x'})
        # Given data selection range in lat-lon or row-col for depth=0m
        #var_data = var_data[0, depth, row_range[0]:row_range[1], col_range[0]:col_range[1]]
        #tmask = np.broadcast_to(tmask,(1,)+tmask.shape) 
        # mask data
        #var_data.data[~np.ma.filled((1 == surf_mask.data))] = np.nan
        var_data = var_data.where(tmask==1)
        #plt.pcolormesh(var_data[0])
        #plt.show()
        
        var_mean = np.nanmean(var_data)
        var_std = np.nanstd(var_data)
        
        if minmax:
            var_min = np.nanmin(var_data)
            var_max = np.nanmax(var_data)
            var_mins.append(var_min)
            mar_maxs.append(var_max)

        y,m,d = get_date(filename, how='ymd')
        date = datetime.date(y, m, d)

        var_means.append(var_mean)
        var_stds.append(var_std)
        dates.append(date)
      
        data.close()
    
    timeseries_var = pd.DataFrame({'date': dates, 'var_mean': var_means, 'var_std': var_stds})
    timeseries_var.to_csv(output_path+runid+'_shbjb_timeseries_data.csv', index=False)



#get_timeseries(runid='EPM151')


def anhalyze_timeseries(timeseries,runid, mhw=True, return_series=False):
    if runid=='EPM151':
        year_standard=2006
        year_min = 2002
        year_max = 2022
        n_year = year_max - year_min + 1
    
    if runid=='EPM111':
        year_standard=2000
        year_min = 1958
        year_max = 2009
        n_year = year_max - year_min + 1

    if runid=='ETW161':
        year_standard = 2008
        year_min = 2002
        year_max = 2018
        n_year = year_max - year_min # both of tahyas experiments are missing a year

    if runid=='ETW162':
        year_standard = 2008
        year_min = 2002
        year_max = 2016
        n_year = year_max - year_min + 1

    if mhw:
        actions = ['g_quantile90', 'add_2T', 'add_3T', 'add_4T']
    else:
        actions = ['g_quantile10', 'remove_2T', 'remove_3T', 'remove_4T']

    anha_timeseries = timeseries.copy()  # why make a copy?

    # change the date formatting
    anha_timeseries['date'] = pd.to_datetime(anha_timeseries['date'], format='%Y-%m-%d')
    
    # add the year, month, and day as a column
    anha_timeseries['year'] = anha_timeseries.date.dt.year
    anha_timeseries['month'] = anha_timeseries.date.dt.month
    anha_timeseries['day'] = anha_timeseries.date.dt.day
    
     # Change date formatting
    anha_timeseries['date2'] = anha_timeseries['date'].dt.date
    anha_timeseries.drop('date', axis=1, inplace=True)  # dropping the og date column for some reason
    anha_timeseries.rename(columns={'date2': 'date'}, inplace=True)
   
    # Add wrap-day column
    anha_timeseries['wrap_day'] = anha_timeseries.apply(lambda row: datetime.date(year_standard, row.month, row.day), axis=1)
    # fold the data to get the day of year statistics
    anha_timeseries_doy_mean = anha_timeseries.groupby('wrap_day')[['var_mean','var_std']].mean().reset_index()
    anha_timeseries_doy_quantile = anha_timeseries.groupby('wrap_day')[['var_mean','var_std']].quantile(.9).reset_index()
    anha_timeseries_doy_max = anha_timeseries.groupby('wrap_day')[['var_mean','var_std']].max().reset_index()
    anha_timeseries_doy_median = anha_timeseries.groupby('wrap_day')[['var_mean','var_std']].median().reset_index()
 

    SST_timeseries = True
    if SST_timeseries:
        plt.figure(figsize=(14,7))
        plt.scatter(anha_timeseries['date'], anha_timeseries.var_mean, c='k', alpha=0.3)

        day_min = datetime.date(year_min-1,12,15)
        day_max = datetime.date(year_max+1,1,15)

        #plt.ylim([-2, 19])
        plt.xlim([day_min, day_max])
        #plt.fill_between([day_min, day_max], [-2,-2],y2=[-2,-2], alpha=0.2)
        #plt.title('James Bay Temperature Timeseries')
        #plt.xlabel('Year')
        plt.ylabel('SST ($^o$C)')
        plt.tight_layout()
        plt.show()
        plt.savefig(fig_path+runid+'_JB_SST_timeseries.png')
        #fig, ax = plt.subplots(1, 1, figsize=(9, 6))
        #ax.xaxis.set_major_formatter(mdates.DateFormatter('%b-%d'))
        #ax.xaxis.set_tick_params(rotation=30, labelsize=10)

   
    anha_timeseries['var_mean_mean'] = np.array(anha_timeseries_doy_mean['var_mean'].to_list()*(n_year))
    anha_timeseries['var_mean_quantile'] = np.array(anha_timeseries_doy_quantile['var_mean'].to_list()*(n_year))

    dT = anha_timeseries['var_mean_quantile'] - anha_timeseries['var_mean_mean']
    
    anha_timeseries['var_mean_2T'] = anha_timeseries['var_mean_mean']+2*dT
    anha_timeseries['var_mean_3T'] = anha_timeseries['var_mean_mean']+3*dT
    anha_timeseries['var_mean_4T'] = anha_timeseries['var_mean_mean']+4*dT
    

    plot_res_temp = True
    if plot_res_temp:
        plt.figure(figsize=(14,7))
        plt.xticks(fontsize=17, rotation=30)
        plt.yticks(fontsize=17)
        #plt.scatter(anha_timeseries['date'], anha_timeseries.var_mean - anha_timeseries.var_mean_mean, s=30, c='k', alpha=0.3)
        plt.plot(anha_timeseries['date'], anha_timeseries.var_mean - anha_timeseries.var_mean_mean, linewidth=2, c='k')
        plt.ylabel('SST Anomaly ($^o$C)', fontsize=17)
        #plt.xlabel('Year')
        plt.xlim([datetime.date(year_min-1,12,1), datetime.date(year_max+1,1,31)])
        #plt.tight_layout()
        plt.savefig(fig_path+runid+'_JB_SST_anomaly.png')
        #plt.show()


    year=2005
    plot_mhw(anha_timeseries, year=year, remove_mean=True, show_cat4=False)
    if return_series:
        anha_timeseries['res_SST'] = anha_timeseries['var_mean']-anha_timeseries['var_mean_mean']
        return anha_timeseries



def plot_mhw(timeseries, year, remove_mean=True, show_cat4=False, region='James Bay', mhw=True):
    timeseries_year = timeseries[timeseries['year']==year].copy()

    if remove_mean:
        labels = ['SST-T$_{c}$', '$\Delta$T', '2$\Delta$T', '3$\Delta$T', '4$\Delta$T']

        # get the freezing line
        freezing_line = timeseries_year['var_mean_mean']*-1

        # remove the yearly average/climatology, T$_{c}$, axis = 1 to apply to whole column
        timeseries_year['var_mean'] = timeseries_year.apply(lambda row: row.var_mean - row.var_mean_mean,axis=1)
        timeseries_year['var_mean_quantile'] = timeseries_year.apply(lambda row: row.var_mean_quantile - row.var_mean_mean, axis=1)  
        timeseries_year['var_mean_2T'] = timeseries_year.apply(lambda row: row.var_mean_2T - row.var_mean_mean, axis=1)
        timeseries_year['var_mean_3T'] = timeseries_year.apply(lambda row: row.var_mean_3T - row.var_mean_mean, axis=1)
        timeseries_year['var_mean_4T'] = timeseries_year.apply(lambda  row: row.var_mean_4T - row.var_mean_mean, axis=1)
        
    else:
        if mhw:
            labels = ['SST', 'T$_{90}$', 'T$_{c}$+2$\Delta$T', 'T$_{c}$+3$\Delta$T', 'T$_{c}$+4$\Delta$T']
        else:
            labels = ['SST', 'T$_{90}$', 'T$_{c}$-2$\Delta$T', 'T$_{c}$-3$\Delta$T', 'T$_{c}$-4$\Delta$T']
        
        
        freezing_line = timeseries_year['var_mean_mean'] * 0. - 2.

    # Setting up MHW/MCS related plotting variables.
    if mhw:
        colors = ['gold', 'orange', 'red', 'maroon']
        alphas = [0.2, 0.3, 0.1]
        mhw_title = 'MHW'
    else:
        colors = ['DeepSkyBlue', 'dodgerblue', 'blue', 'navy']
        alphas = [0.1, 0.5, 0.1]
        mhw_title = 'MCS'


    # Plotting data and categories
    plt.figure(figsize=(9,5))
    
    if not remove_mean:
        plt.plot(timeseries['date'], timeseries.var_mean_mean, c='gray', alpha=.5, label='T$_{c}$')
    
    
    plt.scatter(timeseries_year['date'], timeseries_year.var_mean, c='k', alpha=.3, label=labels[0])

    plt.plot(timeseries_year['date'], timeseries_year.var_mean_quantile, c=colors[0], alpha=alphas[1], label=labels[1])
    
    plt.plot(timeseries_year['date'], timeseries_year.var_mean_2T, c=colors[1], alpha=alphas[1], label=labels[2])
    
    #plt.plot(timeseries_year['date'], timeseries_year.var_mean_3T, c=colors[2], alpha=alphas[1], label=labels[3])


    if show_cat4:
        plt.plot(timeseries_year['date'], timeseries_year.var_mean_4T, c=colors[3], alpha=alphas[1], label=labels[4])
        
        plt.fill_between(timeseries_year['date'], timeseries_year.var_mean_3T, y2=timeseries_year.var_mean_4T, alpha=.1, facecolor=colors[2])


    # Fill in categories colors
    plt.fill_between(timeseries_year['date'], timeseries_year.var_mean_quantile, y2=timeseries_year.var_mean_2T, alpha=alphas[2], facecolor=colors[0])

    #plt.fill_between(timeseries_year['date'], timeseries_year.var_mean_2T,y2=timeseries_year.var_mean_3T,alpha=alphas[2], facecolor=colors[1])
    
    
    if show_cat4:
        plt.fill_between(timeseries_year['date'], timeseries_year.var_mean_3T, y2=timeseries_year.var_mean_4T, alpha=.1, facecolor=colors[2])


    # Set axis limits and labels
    if remove_mean:
        if mhw:
            plt.ylim(bottom=0)
        else:
            plt.ylim(top=2, bottom=-12)
    else:
        plt.ylim(bottom=-2)

    plt.xlim([timeseries_year['date'].iloc[0], timeseries_year['date'].iloc[-1]])
    plt.ylabel(' Residual\n Temperature (\N{DEGREE SIGN}C)', fontsize=12)
    plt.title('%s %i %s' % (region, year, mhw_title), fontsize=14)
    plt.xticks(fontsize=10, rotation=30)
    plt.yticks(fontsize=10)
    plt.legend(fontsize=10)

    plt.tight_layout()
    plt.show()
    #plt.savefig(fig_path+'EPM151_JB1_mhw_'+str(year)+'.png')



def sst_diff(ts1, ts2):
    ts1= ts1[(ts1['year'] > 2002) & (ts1['year'] < 2018)]

    plot_res_temp=True
    
    if plot_res_temp:
        plt.figure(figsize=(10,4))
        plt.scatter(ts2['date'], ts2.var_mean - ts1.var_mean, c='k', alpha=0.3)
        plt.ylabel('Residual Temperature ($^o$C)')
        plt.xlabel('Year')
        plt.tight_layout()
        #plt.savefig(fig_path+'ETW161_JB_residual_temp.png')
        plt.show()





#get_timeseries(runid='ETW162')

timeseries_151 = pd.read_csv(path+'EPM151_shbjb_timeseries_data.csv', index_col=False)
#timeseries_111 = pd.read_csv(path+'EPM111_shbjb_timeseries_data.csv', index_col=False)
#jb_timeseries_162 = jb_timeseries_151.reset_index(drop=True)
#jb_timeseries_162['date'] = pd.to_datetime(jb_timeseries_162['date'], format='%Y-%m-%d')
#jb_timeseries_162 = jb_timeseries_162[(jb_timeseries_162['date'] >= '2002-01-05') & (jb_timeseries_162['date'] < '2017-01-05')]
anhalyze_timeseries(timeseries=timeseries_151, runid='EPM151')
'''
jb_ts161 = pd.read_csv(path+'ETW161_JamesBay_all_region_timeseries_data.csv', index_col=False)
jb_ts161 = jb_ts161.reset_index(drop=True)
ts161 = anhalyze_timeseries(timeseries=jb_ts161, runid='ETW161', return_series=True)

jb_ts151 = pd.read_csv(path+'EPM151_james_bay_timeseries_data.csv', index_col=False)
jb_ts151 = jb_ts151.reset_index(drop=True)  #axis=1, inplace=True)
ts151 = anhalyze_timeseries(jb_ts151,runid='EPM151', return_series=True)
ts151 = ts151[(ts151['year'] >=2002) & (ts151['year'] < 2019)]
print(ts151)
print(ts161)
'''
#sst_diff(ts151, ts161)
#jb1_timeseries_raw = pd.read_csv(path+'EPM151_JamesBay1_timeseries_data.csv', index_col=False)
#jb1_timeseries_raw = jb1_timeseries_raw.reset_index(drop=True)

#jb2_timeseries_raw = pd.read_csv(path+'EPM151_JamesBay2_timeseries_data.csv', index_col=False)
#jb2_timeseries_raw = jb2_timeseries_raw.reset_index(drop=True)

#jb3_timeseries_raw = pd.read_csv(path+'EPM151_JamesBay3_timeseries_data.csv', index_col=False)
#jb3_timeseries_raw = jb3_timeseries_raw.reset_index(drop=True)

#jb4_timeseries_raw = pd.read_csv(path+'EPM151_JamesBay4_timeseries_data.csv', index_col=False)
#jb4_timeseries_raw = jb4_timeseries_raw.reset_index(drop=True)

#anhalyze_timeseries(timeseries=jb_timeseries_riverheat, EPM151=False)



