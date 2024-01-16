import numpy as np
import netCDF4 as nc
import datetime
import pandas as pd
import xarray as xr

# Plotting-related libraries
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import cartopy.feature as feature
import cartopy.crs as ccrs
import matplotlib.dates as mdates
import matplotlib

# OS-specific libraries
import glob

###this is for running script backend###
#import matplotlib
#matplotlib.use('Agg')
###----------------------------------###

fig_path = '/project/6007519/hlouis/plotting/figures/mhw/'
path = '/project/6007519/hlouis/scripts/MarineHeatWaves/'
movie_path = '/project/rrg-pmyers-ad/hlouis/plotting/figures/mhw/movies/EPM151_SST_spatial/'
mask_path = '/project/6007519/hlouis/scripts/masks/'
output_path = '/project/6007519/hlouis/scripts/MarineHeatWaves/'



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


def surf_avg(var):
        n = len(var['depth'])
        weight = np.zeros(n)
        dz = np.zeros(n)
        dd = var['depth'][n-1]

        for i in range(n):
            if i == 0:
                weight[i] = var['depth'][i]/dd
                dz[i] = var['depth'][i]
            else:
                weight[i] = (var['depth'][i] - var['depth'][i-1])/dd
                dz[i] = var['depth'][i] - var['depth'][i-1]

        weights = xr.DataArray(weight, coords=[var['depth']], dims=['depth'])

        #and take the average
        var_weighted = var.weighted(weights)
        surface_variable = var_weighted.mean(dim='depth', skipna=True)
        
        return surface_variable



movie_path = '/project/6007519/hlouis/plotting/figures/horizontal_advection/'
export_path = '/project/6007519/hlouis/data_files/horizontal_advection/'


def advection_calc(runid):
    if runid == 'EPM151':
        data_path = '/project/6007519/pmyers/ANHA4/ANHA4-EPM151-S/'
    
    if runid == 'ETW161':
        data_path = '/project/6007519/weissgib/ANHA4/ANHA4-ETW161-S/'
    
    if runid == 'ETW162':
        data_path = '/project/6007519/weissgib/ANHA4/ANHA4-ETW162-S/' 

    # load the data for u-grid, v-grid, and t-grid
    mdl_files_u = sorted(glob.glob(data_path+'ANHA4-'+runid+'*_gridU.nc'))
    mdl_files_v = sorted(glob.glob(data_path+'ANHA4-'+runid+'*_gridV.nc'))
    mdl_files_t= sorted(glob.glob(data_path+'ANHA4-'+runid+'*_gridT.nc'))
  
    nfiles = len(mdl_files_v)
    
    dates = []
    for filename in mdl_files_t:
        y,m,d = get_date(filename,how='ymd')
        date = datetime.date(y,m,d)
        dates.append(date)
   

    mask_file = '/project/6007519/hlouis/scripts/masks/ANHA4_mesh_mask.nc'
    export_path = '/project/6007519/hlouis/data_files/'
    
    mesh = xr.open_dataset(mask_file) 
    tmask = mesh['tmask'] 
    e3t = mesh['e3t']
    e1t = mesh['e1t']
    e2t = mesh['e1t']
    mesh.close()

    #e1t = e1t.expand_dims(depth=50, axis=1)
    #e2t = e2t.expand_dims(depth=50, axis=1)
    e1t = e1t.rename({'t': 'time_counter'}) 
    e2t = e2t.rename({'t': 'time_counter'}) 


    for ifile, filename in enumerate(mdl_files_t):
        #name = filename.split('/')[-1]
        filename = filename.split('gridT')[0]
        
        du = xr.open_dataset(filename+'gridU.nc') 
        dv = xr.open_dataset(filename+'gridV.nc')
        dt = xr.open_dataset(filename+'gridT.nc')
        
        du = du.rename({'depthu': 'depth'})
        dv = dv.rename({'depthv': 'depth'})
        dt = dt.rename({'deptht': 'depth'})
        dt = dt.rename({'x_grid_T': 'x', 'y_grid_T': 'y'})
        
        lats = dt.coords['nav_lat_grid_T']
        lons = dt.coords['nav_lon_grid_T']
        
        mld = dt['somxl010']
        t = dt['votemper']
        u = du['vozocrtx']   
        v = dv['vomecrty']
        
        '''
        u_center = 0.5 * (u[:,:,:, 0:-2] + u[:,:,:, 1:-1])
        v_center = 0.5 * (v[:,:,0:-2,:] + v[:,:,1:-1,:])
       
        amountToPad = u.x.size - u_center.x.size
        u_center = u_center.pad(x=(0,amountToPad), constant_values=np.nan)
        u_center = u_center.roll(x=1, roll_coords=True)

        amountToPad = v.y.size - v_center.y.size
        v_center = v_center.pad(y=(0,amountToPad), constant_values=np.nan)
        v_center = v_center.roll(y=1,roll_coords=True)
        '''
        
        entrainment_layer = t.depth <= mld 

        
        u_mld = u.where(entrainment_layer) 
        v_mld = v.where(entrainment_layer)
        t_mld = t.where(entrainment_layer)
        
        u_mld = surf_avg(u_mld)
        v_mld = surf_avg(v_mld)
        t_mld = surf_avg(t_mld)
        
        dT_dx = np.gradient(t_mld, axis=2)/e1t 
        dT_dy = np.gradient(t_mld, axis=1)/e2t
                


        ''' 
        dTx = t_mld[:,:,:,0:-2] - t_mld[:,:,:, 1:-1]
        amountToPad = t_mld.x.size - t_mld.y.size
        dTx = dTx.pad(y=(0,amountToPad), constant_values=np.nan)
        dTx = dTx.roll(x=1,roll_coords=True)
        '''
        #print('u',u_mld)
        #print('v',v_mld)
        #print('dtdx',dT_dx)
        #print('dtdy',dT_dy)
       
        advection = -(u_mld*dT_dx + v_mld*dT_dy)
        print(advection)
        advection = advection.rename('advection')
        #advection.to_netcdf(export_path+'EPM151_'+str(dates[ifile])+'_horizontal_advection.nc')
        print('exported', filename, str(ifile)+'/'+str(nfiles))


        dv.close()
        du.close()
        dt.close()
        u.close()
        v.close()
        u_mld.close()
        v_mld.close()
        t_mld.close()
        dT_dx.close()
        dT_dy.close()
        advection.close()
    
         

#advection_calc('EPM151')



def fmt(x, pos):
    a, b = '{:.2e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)



def plot_advection(runid):
    if runid == 'EPM151':
        data_path = '/project/6007519/pmyers/ANHA4/ANHA4-EPM151-S/'

    if runid == 'ETW161':
        data_path = '/project/6007519/weissgib/ANHA4/ANHA4-ETW161-S/'
    
    if runid == 'ETW162':
        data_path = '/project/6007519/weissgib/ANHA4/ANHA4-ETW162-S/' 

    file_path = '/project/6007519/hlouis/data_files/'#horizontal_advection/'
    mask_path = '/project/6007519/hlouis/scripts/masks/'
    
    mdl_files = sorted(glob.glob(file_path+runid+'*horizontal_advection.nc'))
    mf = nc.Dataset(mask_path+'HBC_mask.nc')
    mask = mf['jb'][0,0]
    mf.close()
    file_list = sorted(glob.glob(data_path+'ANHA4-'+runid+'*_gridT.nc'))

    dates = []
    
    for filename in file_list:
        y,m,d = get_date(filename,how='ymd')
        date = datetime.date(y,m,d)
        dates.append(date)

    for ifile, file in enumerate(mdl_files):
        d = xr.open_dataset(file)
        surface_advection = d['advection']
        surface_advection = surface_advection.squeeze('t')
        print(surface_advection)
        lats = d.coords['nav_lat_grid_T']
        lons = d.coords['nav_lon_grid_T']
        '''
        n = len(advect['depth'])
        weight = np.zeros(n)
        dz = np.zeros(n)
        dd = advect['depth'][n-1]
        for i in range(n):
            if i == 0:
                weight[i] = advect['depth'][i]/dd
                dz[i] = advect['depth'][i]
            else:
                weight[i] = (advect['depth'][i] - advect['depth'][i-1])/dd
                dz[i] = advect['depth'][i] - advect['depth'][i-1]

        weights = xr.DataArray(weight, coords=[advect['depth']], dims=['depth'])

        #and take the average
        advect_weighted = advect.weighted(weights)
        surface_advection = advect_weighted.mean(dim='depth', skipna=True)
        '''

        land_50m = feature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='black', facecolor='gray', linewidth=0.5)
        projection=projection=ccrs.Mercator(central_longitude=-80)
        fig = plt.figure(figsize=(6,5))
        ax = plt.subplot(1, 1, 1, projection=projection)
        #ax.set_title(str(dates[ifile]))
        plt.title(str(dates[ifile]))
        ax.set_extent([-82.5, -78.5, 51, 55.86], crs=ccrs.PlateCarree())
        ax.add_feature(land_50m, color=[0.8, 0.8, 0.8])
        ax.coastlines(resolution='50m')
        p1 = ax.pcolormesh(lons, lats, surface_advection[0], transform=ccrs.PlateCarree(), cmap='viridis', vmin=-3*10**-6, vmax=2*10**-6)
        ax_cb = plt.axes([0.70, 0.25, 0.015, 0.5])
        fmt = '%1.2f'
        cb = plt.colorbar(p1,cax=ax_cb, orientation='vertical')
        cb.formatter.set_powerlimits((0, 0))
        cb.formatter.set_useMathText(True)
        cb.ax.yaxis.set_offset_position('left')
        cb.update_ticks()
        cb.ax.set_ylabel('Horizontal Advection in MLD ($^o$C/s)')
        #ax.gridlines()
        #plt.tight_layout()
        plt.savefig(movie_path+'EPM151-JB_horizontal_advection_'+str(dates[ifile])+'.png')
        #plt.show()
        plt.close()
         
    '''
    d = xr.open_mfdataset(mdl_files)

    lats = d.coords['nav_lat_grid_T']
    lons = d.coords['nav_lon_grid_T']

    advection = d['advection']
    
    adv_masked = np.ma.masked_where(mask==1, advection)

    land_50m = feature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='black', facecolor='gray', linewidth=0.5)
    projection=projection=ccrs.Mercator(central_longitude=-80)
    fig = plt.figure()
    ax = plt.subplot(1, 1, 1, projection=projection)
    ax.set_title(str(dates[day]))
    ax.set_extent([-82.5, -78.5, 51, 55], crs=ccrs.PlateCarree())
    ax.add_feature(land_50m, color=[0.8, 0.8, 0.8])
    ax.coastlines(resolution='50m')
    p1 = ax.pcolormesh(lons, lats, adv_masked, transform=ccrs.PlateCarree(), cmap='gist_ncar', vmin=-0.005, vmax=0.005)
    ax_cb = plt.axes([0.80, 0.25, 0.015, 0.5])
    cb = plt.colorbar(p1,cax=ax_cb, orientation='vertical')
    cb.ax.set_ylabel('Horizontal Advection in MLD ($^o$C/s)')
    ax.gridlines()
    fig.tight_layout()
    plt.savefig(movie_path+'hor_adv'+str(dates[day])+'.png')
    plt.close()
    '''

plot_advection('EPM151')




def get_timeseries(runid, minmax=False):
    if runid == 'EPM151':
        data_path = '/project/6007519/pmyers/ANHA4/ANHA4-EPM151-S/'
    
    if runid == 'ETW161':
        data_path = '/project/6007519/weissgib/ANHA4/ANHA4-ETW161-S/'
    
    if runid == 'ETW162':
        data_path = '/project/6007519/weissgib/ANHA4/ANHA4-ETW162-S/' 

    depth = 0
    ANHA4 = False
    if ANHA4:
        mask = nc.Dataset(mask_path + 'ANHA4_mesh_mask.nc')
        tmask = mask['tmask'][0][depth]
    else:
        locations = ['JamesBay1', 'JamesBay2', 'JamesBay3', 'JamesBay4', 'JamesBay_all']
        loc = locations[4]
        loc='jb'
        mask = nc.Dataset(mask_path + 'HBC_mask.nc') #'JB_all_masks.nc')  
        tmask = mask[loc][0,0]  # depth is already surface and there is no time dimension in hommade masks


    file_list = sorted(glob.glob(data_path+'ANHA4-'+runid+'*_gridT.nc'))
    nfiles = len(file_list)
    var_means = []
    var_stds = []
    var_mins = []
    var_maxs = []
    dates = []
    
    for filename in file_list:
        y,m,d = get_date(filename,how='ymd')
        date = datetime.date(y,m,d)
        dates.append(date)
    

    d = xr.open_mfdataset(file_list, concat_dim='time_counter', combine='nested', data_vars='minimal', coords='minimal', compat='override')
    
    lats = d.coords['nav_lat_grid_T']
    lons = d.coords['nav_lon_grid_T']
    
    doy_avg = d.votemper.groupby('time_counter.dayofyear').mean('time_counter')
    doy_gb = d.votemper.groupby('time_counter.dayofyear')

    d_gb = doy_gb - doy_avg #doy_gb.mean('time_counter')



    for day in range(1533): #d_gb.dims['time_counter']):
        Tc = d_gb[day,0].values  # doy_avg['votemper'][day,0].values
        Tc_masked = np.ma.masked_where(tmask==1, Tc)
        
        
        SST_Tc = Tc_masked
        
        



#get_timeseries(runid='EPM151')


