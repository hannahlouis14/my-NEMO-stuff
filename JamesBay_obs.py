import datetime
import xarray as xr
import netCDF4 as nc
import numpy.ma as ma
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import gsw
import cartopy.crs as ccrs
import cartopy.feature as feature
import glob


def section_calculation(x1, x2, y1, y2):

    ii = []
    jj = []

    dx = x2-x1
    dy = y2-y1
    yi = 1

    if dy < 0:
        yi = -1
        dy = -dy

    D = (2*dy) - dx
    y = y1

    for x in range(x1,x2):
        ii.append(x)
        jj.append(y)
        if D > 0:
            y = y +yi
            D = D + (2*(dy-dx))
            ii.append(x)
            jj.append(y)
        else:
            D = D + 2*dy

    return ii, jj



def temp_depth(year, month, section):
    path = '/project/6007519/hlouis/data_files/obs/'
    fig_path = '/project/6007519/hlouis/plotting/figures/ts/'
    data_file = path+'ctd_data-'+year+'.csv'
    
    ## read in the data as a dataframe and set the variables we want
    data = pd.read_csv(data_file, keep_default_na=True)
    df = pd.DataFrame(data)
    temperature = df.iloc[:,12]
    salinity = df.iloc[:,13]
    depth = df.iloc[:,11]
    lon = df.iloc[:,3]
    lat = df.iloc[:,4]


    ## parsing out the observational data for the section we want by getting the indices at those locations
    if section == 'AB':
        lon_ind = lon.loc[(lon >= -82) & (lon <= -79.85)].index.tolist()
        lat_ind = lat.loc[(lat >= 54.2) & (lat <= 54.4)].index.tolist()
        ind = list(set(lon_ind).intersection(lat_ind))


    if section == 'CD':
        lon_ind = lon.loc[(lon >= -82) & (lon <= -79.4)].index.tolist()
        lat_ind = lat.loc[(lat >= 53.76) & (lat <= 53.83)].index.tolist()
        ind = list(set(lon_ind).intersection(lat_ind))


    if section == 'EF':
        lon_ind = lon.loc[(lon >= -80.65) & (lon <= -78.49)].index.tolist()
        lat_ind = lat.loc[(lat >= 52.2) & (lat <= 52.44)].index.tolist()
        ind = list(set(lon_ind).intersection(lat_ind))

    
    if section == 'alpha_beta':
        lon_ind = lon.loc[(lon >= -81.89) & (lon <= -81.68)].index.tolist()
        lat_ind = lat.loc[(lat >= 53.81) & (lat <= 54.3)].index.tolist()
        ind = list(set(lon_ind).intersection(lat_ind))

    
    if section == 'epsilon_zeta':
        lon_ind = lon.loc[(lon >= -80.1) & (lon <= -79.78)].index.tolist()
        lat_ind = lat.loc[(lat >= 53.838) & (lat <= 54.329)].index.tolist()
        ind = list(set(lon_ind).intersection(lat_ind))


    if section == 'eta_theta':
        lon_ind = lon.loc[(lon >= -79.6) & (lon <= -79.35)].index.tolist()
        lat_ind = lat.loc[(lat >= 52.4) & (lat <= 53.8)].index.tolist()
        ind = list(set(lon_ind).intersection(lat_ind))


    if section == 'lower_JB':
        lon_ind = lon.loc[(lon >= -80.7) & (lon <= -79.63)].index.tolist()
        lat_ind = lat.loc[(lat >= 50) & (lat <= 52.35)].index.tolist()
        ind = list(set(lon_ind).intersection(lat_ind))


    if section == 'inner_JB':
        lon_ind = lon.loc[(lon >= -80.71) & (lon <= -79.75)].index.tolist()
        lat_ind = lat.loc[(lat >= 52.38) & (lat <= 53.76)].index.tolist()
        ind = list(set(lon_ind).intersection(lat_ind))


    temperature = temperature.loc[ind].values
    salinity = salinity.loc[ind].values
    depth = depth.loc[ind].values
    lon = lon.loc[ind].values
    lat = lat.loc[ind].values

    plt.figure(figsize=(6,5))
    plt.scatter(temperature, -depth, marker='.', linewidth=1, label='observation')
    #plt.scatter(mod_sal, mod_temp, c='r', marker='s', linewidth=2, label='model')
    plt.ylabel('Depth (m)'); plt.xlabel('Temperature ($^o$C)')
    plt.title('James Bay Observation and Model Temperature Profile '+section+' August 2021')
    plt.legend(loc='upper right')
    plt.show()
    #plt.savefig(fig_path+'TS-obs_model_section-'+section+'_Aug2021.png')



def ts(year, month, section):
    path = '/project/6007519/hlouis/data_files/obs/'
    fig_path = '/project/6007519/hlouis/plotting/figures/ts/'
    data_file = path+'ctd_data-'+year+'.csv'
    
    ## read in the data as a dataframe and set the variables we want
    data = pd.read_csv(data_file, keep_default_na=True)
    df = pd.DataFrame(data)
    date = df.iloc[:,7]
    month = df.iloc[:,6]
    temperature = df.iloc[:,12]
    salinity = df.iloc[:,13]
    depth = df.iloc[:,11]
    lon = df.iloc[:,3]
    lat = df.iloc[:,4]
    

    ## parsing out the observational data for the section we want by getting the indices at those locations
    if section == 'jb_AB':
        lon_ind = lon.loc[(lon >= -82.05) & (lon <= -79.82)].index.tolist()
        lat_ind = lat.loc[(lat >= 54.2) & (lat <= 54.4)].index.tolist()
        ab_ind = list(set(lon_ind).intersection(lat_ind))
        
    if section =='jb_CD':
        lon_ind = lon.loc[(lon >= -82) & (lon <= -79.4)].index.tolist()
        lat_ind = lat.loc[(lat >= 53.76) & (lat <= 53.83)].index.tolist()
        ab_ind = list(set(lon_ind).intersection(lat_ind))
    
    if section == 'jb_EF':
        lon_ind = lon.loc[(lon >= -80.65) & (lon <= -78.49)].index.tolist()
        lat_ind = lat.loc[(lat >= 52.2) & (lat <= 52.44)].index.tolist()
        ab_ind = list(set(lon_ind).intersection(lat_ind))

    if section == 'alpha_beta':
        lon_ind = lon.loc[(lon >= -81.89) & (lon <= -81.68)].index.tolist()
        lat_ind = lat.loc[(lat >= 53.81) & (lat <= 54.3)].index.tolist()
        ab_ind = list(set(lon_ind).intersection(lat_ind))

    if section == 'epsilon_zeta':
        lon_ind = lon.loc[(lon >= -80.1) & (lon <= -79.78)].index.tolist()
        lat_ind = lat.loc[(lat >= 53.83) & (lat <= 54.329)].index.tolist()
        ab_ind = list(set(lon_ind).intersection(lat_ind))
    
    if section == 'eta_theta':
        lon_ind = lon.loc[(lon >= -79.6) & (lon <= -79.35)].index.tolist()
        lat_ind = lat.loc[(lat >= 52.4) & (lat <= 53.8)].index.tolist()
        ab_ind = list(set(lon_ind).intersection(lat_ind))

    if section == 'psi_omega':
        lon_ind = lon.loc[(lon >= -81) & (lon <= -80.17)].index.tolist()
        lat_ind = lat.loc[(lat >= 50) & (lat <= 51.6)].index.tolist()
        ab_ind = list(set(lon_ind).intersection(lat_ind))

    if section == 'lower_JB':
        lon_ind = lon.loc[(lon >= -80.7) & (lon <= -79.63)].index.tolist()
        lat_ind = lat.loc[(lat >= 50) & (lat <= 52.35)].index.tolist()
        ab_ind = list(set(lon_ind).intersection(lat_ind))


    if section == 'inner_JB':
        lon_ind = lon.loc[(lon >= -80.71) & (lon <= -79.75)].index.tolist()
        lat_ind = lat.loc[(lat >= 52.38) & (lat <= 53.76)].index.tolist()
        ab_ind = list(set(lon_ind).intersection(lat_ind))


    month = month.loc[ab_ind].values
    date = date.loc[ab_ind].values
    temperature = temperature.loc[ab_ind].values
    salinity = salinity.loc[ab_ind].values
    depth = depth.loc[ab_ind].values
    lon = lon.loc[ab_ind].values
    lat = lat.loc[ab_ind].values


    # Define the min / max values for plotting isopycnals
    t_min = temperature.min() - 5
    t_max = temperature.max() + 5
    s_min = salinity.min() - 2
    s_max = salinity.max() + 7
    
    
    # Calculate how many gridcells we need in the x and y dimensions
    #xdim = np.ceil(29 - s_min)/0.01
    #ydim = np.ceil(t_max - (-2))
    xdim = np.ceil(s_max - s_min)/0.1
    ydim = np.ceil(t_max-t_min)
    dens = np.zeros((int(ydim),int(xdim)))
    
    # Create temp and salt vectors of appropiate dimensions
    #ti = np.linspace(-4,15,int(ydim))+t_min
    #si = np.linspace(1,int(xdim),int(xdim))*0.1+s_min
    ti = np.linspace(0,int(ydim),int(ydim))+t_min
    si = np.linspace(1,int(xdim),int(xdim))*0.1+s_min
    
    
    # Loop to fill in grid with densities
    for j in range(0,int(ydim)):
        for i in range(0, int(xdim)):
            dens[j,i]=gsw.rho(si[i],ti[j],0)
    
    dens = dens - 1000
    mod_sal, mod_temp = ts_nemo('EPM151',section=section, endyear=2022, endmonth=12, endday=31)

    plt.figure(figsize=(10,6))
    CS = plt.contour(si,ti,dens, linestyles='dashed', colors='k')  # this plots the isopycnals
    plt.clabel(CS, fontsize=12, inline=1, fmt='%.2f') # Label every second level
    plt.scatter(salinity,temperature, c=depth, cmap=plt.cm.viridis, lw=2, label='observation')
    plt.scatter(mod_sal, mod_temp, c='r', marker='s', linewidth=2, label='model')
    #plt.ylim([t_min+0.75, t_max- 0.75])
    plt.ylim(bottom=-2, top=13.2)
    plt.xlim(left=20.5, right=30)
    plt.clim(0,np.max(depth)+0.5)  # bottom=-3, top=8
    plt.xlabel('Salinity (PSU)'); plt.ylabel('Temperature ($^o$C)')
    plt.title('James Bay Observation and Model T-S '+section+' August 2021')
    #plt.title('James Bay Observation T-S '+month+' '+year)
    plt.legend(loc='upper right')
    plt.colorbar(label='Depth (m)')
    #plt.show()
    #plt.savefig(fig_path+'TS-'+month+'-'+year+'.png')
    plt.savefig(fig_path+'TS-obs_model_section-'+section+'_Aug2021.png')
    

def ts_nemo(runid, section, endyear, endmonth, endday, startyear=2002, startmonth=1, startday=5):
    path = '/project/6007519/pmyers/ANHA4/ANHA4-'+runid+'-S/' 
    grid_file = '/project/6007519/weissgib/plotting/data_files/anha4_files/ANHA4_mesh_mask.nc' 
    output_path = '/project/6007519/hlouis/plotting/figures/ts/' 
    mask_file = '/project/6007519/hlouis/scripts/JB_obs2021_masks.nc'
   
    '''
    ## selecting the dates each secton coincides with
    if section == 'jb_AB':
        mdl_files = [path+'ANHA4-EPM151_y2021m08d03_gridT.nc', path+'ANHA4-EPM151_y2021m08d08_gridT.nc']
        ii, jj = section_calculation(73,83,364,362)

    if section == 'jb_CD':
        mdl_files = [path+'ANHA4-EPM151_y2021m08d03_gridT.nc', path+'ANHA4-EPM151_y2021m08d08_gridT.nc', path+'ANHA4-EPM151_y2021m08d13_gridT.nc']
        ii, jj = section_calculation(73,85,358,357)

    if section == 'jb_EF':
        mdl_files = [path+'ANHA4-EPM151_y2021m08d08_gridT.nc', path+'ANHA4-EPM151_y2021m08d13_gridT.nc']
        ii, jj = section_calculation(76,87,346,345)

    if section == 'alpha_beta':
        mdl_files = [path+'ANHA4-EPM151_y2021m08d08_gridT.nc']
        ii, jj = section_calculation(73,76,364,348)
        #ii, jj = section_calculation(74,75,359,355)

    if section == 'epsilon_zeta':
        mdl_files = [path+'ANHA4-EPM151_y2021m08d08_gridT.nc']
        ii, jj = section_calculation(82,83,362,357)

    if section == 'eta_theta':
        mdl_files = [path+'ANHA4-EPM151_y2021m08d13_gridT.nc']
        ii, jj = section_calculation(82,83,357,345)

    if section == 'lower_JB':
        mdl_files = [path+'ANHA4-EPM151_y2021m08d08_gridT.nc', path+'ANHA4-EPM151_y2021m08d13_gridT.nc']
        #ii, jj = section_calculation()

    if section == 'inner_JB':
        mdl_files = [path+'ANHA4-EPM151_y2021m08d08_gridT.nc']
        #ii, jj = section_calculation()

    if section == 'psi_omega':
        mdl_files = [path+'ANHA4-EPM151_y2021m08d08_gridT.nc', path+'ANHA4-EPM151_y2021m08d13_gridT.nc']
        ii, jj = section_calculation(78,79,337,341)
    '''
    start_time = datetime.date(startyear, startmonth, startday)

    end_time = datetime.date(endyear, endmonth, endday)

    #figure out all the dates we have model files
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


    d = xr.open_mfdataset(mdl_files, concat_dim='time_counter', combine='nested', data_vars='minimal', coords='minimal', compat='override', use_cftime=True)
    
    temp = d['votemper'].values
    sal = d['vosaline'].values
    t = d.dims['time_counter']
    print(temp)
     
    ## use the mask file I made to mask the data and get the depth?
    mf = nc.Dataset(mask_file)
    rmask = mf[section][0]
    rmask = np.broadcast_to(rmask, (t,)+rmask.shape)
    #temp = np.reshape(temp, (50,800,544))
    #sal = np.reshape(sal, (50,800,544))
    #print(np.shape(temp))

    ## masking the data
    depth = d['deptht'].values
    masked_temp = np.where(rmask==2, temp, np.nan)
    masked_sal = np.where(rmask==2, sal, np.nan)
    #masked_temp = ma.masked_where(rmask==1, temp)
    #masked_sal = ma.masked_where(rmask==1,sal)
    #masked_temp[masked_temp==0] = np.nan
    #masked_sal[masked_sal==0] = np.nan
    masked_depth = np.where(rmask==2, depth, np.nan)
    rmask2d = mf[section][0,0]
    
    ## get the lon and lat data from the ANHA4 mask
    mesh = nc.Dataset(grid_file)
    nav_lon = np.array(mesh.variables['nav_lon'])
    nav_lat = np.array(mesh.variables['nav_lat'])
    mesh.close()
    masked_lon = np.where(rmask2d==2, nav_lon, np.nan)
    masked_lat = np.where(rmask2d==2, nav_lat, np.nan)
    
    '''
    #get the grid coordinates on map coordinates
    lon = []
    lat = []
    t = []
    s = []


    for n in range(0,len(ii)-1):
        i1 = int(ii[n])
        i2 = int(ii[n+1])
        j1 = int(jj[n])
        j2 = int(jj[n+1])
        lon.append(nav_lon[j1,i1])
        lat.append(nav_lat[j1,i1]) 
        temp = d['votemper'].isel(y_grid_T=j1, x_grid_T=i1)
        sal = d['vosaline'].isel(y_grid_T=j1, x_grid_T=i1)
        sal = sal.values
        temp = temp.values

        for date in range(len(sal)):
            k = sal[date]
            k[k==0] = np.nan
            new_sal = k  # np.trim_zeros(k, 'b')
            s.append(new_sal)

        for date in range(len(temp)):
            k = temp[date]
            k[k==0] = np.nan
            new_temp = k  #np.trim_zeros(k, 'b')
            t.append(new_temp)
    '''
    '''
    print(masked_depth)
    plt.scatter(masked_sal, masked_temp)
    #plt.scatter(s,t, c='r')
    plt.show()
    '''
    '''
    land_50m = feature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='black', facecolor='gray', linewidth=0.5)
    projection=ccrs.PlateCarree()  #LambertConformal()
    fig = plt.figure(figsize=(6, 5))
    ax = plt.subplot(1, 1, 1, projection=projection)
    ax.set_extent([-83, -77.8, 50.8, 55], crs=ccrs.PlateCarree())
    ax.add_feature(land_50m, color=[0.8, 0.8, 0.8])
    ax.coastlines(resolution='50m')
    ax.gridlines()
    ax.scatter(masked_lon, masked_lat, transform=ccrs.PlateCarree())
    #ax.scatter(lon, lat, c='r', transform=ccrs.PlateCarree())
    plt.show()
    '''
    return masked_sal, masked_temp 
    


def plot_location(year, section):
    path = '/project/6007519/hlouis/data_files/obs/'
    fig_path = '/project/6007519/hlouis/plotting/figures/ts/'
    data_file = path+'ctd_data-'+year+'.csv'
    data = pd.read_csv(data_file, keep_default_na=True)
    df = pd.DataFrame(data)

    lon = df.iloc[:,3]
    lat = df.iloc[:,4]
    date = df.iloc[:,7]
    
    lon_ind = lon.loc[(lon >= -80.71) & (lon <= -79.75)].index.tolist()
    lat_ind = lat.loc[(lat >= 52.38) & (lat <= 53.76)].index.tolist()
    ind = list(set(lon_ind).intersection(lat_ind))
    
    date = date.loc[ind].values
    lon = lon.loc[ind].values
    lat = lat.loc[ind].values
    
    land_50m = feature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='black', facecolor='gray', linewidth=0.5)
    projection=ccrs.PlateCarree()  #LambertConformal()
    fig = plt.figure(figsize=(6, 5))
    ax = plt.subplot(1, 1, 1, projection=projection)
    ax.set_extent([-83, -77.8, 50.8, 55], crs=ccrs.PlateCarree())
    ax.add_feature(land_50m, color=[0.8, 0.8, 0.8])
    ax.coastlines(resolution='50m')
    ax.gridlines()
    ax.scatter(lon, lat, transform=ccrs.PlateCarree(), marker='.', lw=0)
    plt.show()
    #plt.savefig(fig_path+'JB_obs_section-'+section+'_2021.png') 


def dissO2(runid, year, endyear, endmonth, endday, startyear=2021, startmonth=8, startday=3):
    path = '/project/6007519/hlouis/data_files/obs/' 
    mod_path = '/project/6007519/pmyers/ANHA4/ANHA4-'+runid+'-S/'
    grid_file = '/project/6007519/weissgib/plotting/data_files/anha4_files/ANHA4_mesh_mask.nc'
    output_path = '/project/6007519/hlouis/plotting/figures/ts/'
    mask_file = '/project/6007519/hlouis/scripts/HBC_mask.nc'
    data_file = path+'ctd_data-'+year+'.csv'

    '''
    start_time = datetime.date(startyear, startmonth, startday)

    end_time = datetime.date(endyear, endmonth, endday)

    #figure out all the dates we have model files
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

    #and now make a list of model files to read
    mdl_files = []
    for t in times:
        mdl_files.append(mod_path+"ANHA4-"+runid+"_y"+str(t.year)+"m"+str(t.month).zfill(2)+"d"+str(t.day).zfill(2)+"_gridT.nc")

    d = xr.open_mfdataset(mdl_files, concat_dim='time_counter', combine='nested', data_vars='minimal', coords='minimal', compat='override', use_cftime=True)
    '''
    oxy = d['vooxy'].values
    mod_depth = d['deptht'].values
    t = d.dims['time_counter']
    print(mod_depth)
    mf = nc.Dataset(mask_file)
    rmask = mf['hbc'][0]
    rmask = np.broadcast_to(rmask, (t,)+rmask.shape)

    ## masking the data
    oxy = np.where(rmask==2, oxy, np.nan)
    print(np.shape(oxy))
    ## read in the data as a dataframe and set the variables we want
    data = pd.read_csv(data_file, keep_default_na=True)
    df = pd.DataFrame(data)
    date = df.iloc[:,7]
    depth = df.iloc[:,11]
    lon = df.iloc[:,3]
    lat = df.iloc[:,4]
    o2 =  df.iloc[:,14]
   
    '''
    #CS = plt.pcolormesh(lat, depth, o2, linestyles='dashed', colors='k')  # this plots the isopycnals
    #plt.clabel(CS, fontsize=12, inline=1, fmt='%.2f')
    #plt.colorbar(label='diss O$_2$')
    plt.scatter(depth,o2, lw=1)
    plt.scatter(oxy, mod_depth)
    plt.ylabel('depth (m)')
    plt.xlabel('dissolved O$_2$')
    plt.show()
    '''



#dissO2(runid='EPM151', year='2021', endyear=2021, endmonth=8, endday=28)
#plot_location('2021', section='inner_JB')
#ts_nemo(runid='EPM151', section='inner_JB', endyear=2022, endmonth=12, endday=31)
ts(year='2021',month='JULY', section='inner_JB') 
#ts(year='2022', month='AUG')
#temp_depth(year='2021', month='AUG', section='AB')
