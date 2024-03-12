import numpy as np
import datetime
import xarray as xr
import netCDF4 as nc
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as feature
import glob
import numpy.ma as ma
import warnings
import matplotlib



def sea_ice_concentration(runid, endyear, endmonth,endday, startyear, startmonth, startday):
    fig_path = '/project/6007519/hlouis/plotting/figures/sea_ice/'
    obs_path = '/project/6007519/hlouis/AMSR2/'
    mdl_path = '/project/6007519/pmyers/ANHA4/ANHA4-EPM151-S/'
    grid_file = '/project/6007519/hlouis/scripts/masks/ANHA4_mesh_mask.nc'
    mask_file = '/project/6007519/hlouis/scripts/masks/shbjb_mask.nc'

    mesh = nc.Dataset(grid_file)
    lons = np.array(mesh.variables['nav_lon'])
    lats = np.array(mesh.variables['nav_lat'])
    mesh.close()
    mf = nc.Dataset(mask_file)
    mask = mf['tmask']

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
         mdl_files.append(mdl_path+"ANHA4-"+runid+"_y"+str(t.year)+"m"+str(t.month).zfill(2)+"d"+str(t.day).zfill(2)+"_icemod.nc")

    months_out = ['m04', 'm05', 'm06', 'm07', 'm08', 'm09', 'm10', 'm11', 'm12']  # for winter months
    mdl_files = [x for x in mdl_files if not any(word in x for word in months_out)]  # now, we have only JFM months for our model paths
    
    d = xr.open_mfdataset(mdl_files, concat_dim='time_counter', data_vars='minimal', coords='minimal', compat='override')
    # next, we want to take the yearly average for the winter months we have
    annual_avg = d.groupby('time_counter.year').mean('time_counter')
    yrs = annual_avg['year'].values
    ice = annual_avg['ileadfra']
    
    # now, we want to load in the AMSR2 satellite data and take the 5 day average for the whole set (2012-2022) 
    years = np.arange(2013,2023,1)
    
    for y, year in enumerate(years):
        directory = obs_path+str(year)+'/'
        files = sorted(glob.glob(directory+'*.nc'))    
        months_out = [str(year)+'04', str(year)+'05', str(year)+'06', str(year)+'07', str(year)+'08', str(year)+'09',str(year)+'10',str(year)+'11',str(year)+'12', str(year)+'0229']
        files = [x for x in files if not any(word in x for word in months_out)]
        ds = xr.open_mfdataset(files, combine='nested', concat_dim=None, data_vars='minimal', coords='minimal', compat='override')
        sic = ds['sic'].values
        masked_sic = ma.masked_where(mask==1, sic)
        masked_sic[masked_sic==0] = np.nan
        
       # now, we want to mask the model data
        sea_ice = annual_avg['ileadfra'].isel(year=y).values  # ice[y].values
        masked_ice = ma.masked_where(mask==1,sea_ice)
        masked_ice[masked_ice==0] = np.nan
        
        diff = (masked_ice - masked_sic)/masked_sic * 100  #(masked_sic - masked_ice)/masked_ice 
        vmin = -50
        vmax = 100
        #norm = mpl.colors.TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)

        land_50m = feature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='black', facecolor='gray', linewidth=0.5)
        projection=ccrs.Mercator(central_longitude=-80)
        fig = plt.figure(figsize=(6,5))
        ax = plt.subplot(1,1,1, projection=projection)
        ax.set_extent([-95, -76, 64.5, 51], crs=ccrs.PlateCarree())
        ax.add_feature(land_50m, color=[0.8, 0.8, 0.8])
        ax.coastlines(resolution='50m')
        p1 = ax.pcolormesh(lons, lats, diff, transform=ccrs.PlateCarree(), vmin=vmin, vmax=vmax, cmap='viridis')  
        #p1 = ax.pcolormesh(lons, lats, diff, transform=ccrs.PlateCarree(), norm=matplotlib.colors.SymLogNorm(vmin=vmin, vmax=vmax, linthresh=15, linscale=1), cmap='viridis')
        ax_cb = plt.axes([0.80, 0.25, 0.015, 0.5])
        cb = plt.colorbar(p1,cax=ax_cb, orientation='vertical')
        cb.ax.set_ylabel('Sea Ice Concentration Difference (%)')
        ax.set_title(runid+' Model vs Observation '+str(year))

        #fig.tight_layout()
        #plt.savefig(fig_path+'sea_ice_concentration_diff_ice-sic_normalscale_'+runid+'_'+str(year)+'.png')
        plt.show()
        plt.clf()
        
        ds.close()

if __name__ == "__main__":
    sea_ice_concentration(runid='EPM151', endyear=2022, endmonth=12, endday=31, startyear=2012, startmonth=1, startday=5) 
