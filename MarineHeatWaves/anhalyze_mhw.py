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

depth = 0                    # z axis location from 50 unit "depth" 



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
    
    file_list = sorted(glob.glob(data_path+'**/*_gridT.nc'))
    nfiles = len(file_list)
    tmask = np.broadcast_to(tmask,(1,)+tmask.shape)

    var_means = []
    var_stds = []
    var_mins = []
    var_maxs = []
    dates = []


    for filename in file_list:  
        # get the ANHA4 data
        data = xr.open_dataset(filename)
        
        var = 'votemper'
        var_data = data[var][:,0]
        var_data = var_data.rename({'y_grid_T': 'y', 'x_grid_T': 'x'})
        
        # mask data
        var_data = var_data.where(tmask==1)
        
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



get_timeseries(runid='EPM111')


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







#timeseries_151 = pd.read_csv(path+'EPM151_shbjb_timeseries_data.csv', index_col=False)
#anhalyze_timeseries(timeseries=timeseries_151, runid='EPM151')

