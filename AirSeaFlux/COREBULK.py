import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import glob
import datetime
import netCDF4 as nc
import cmocean



def assign_timestamp(og, mod):
    timestamp = og.coords['time_counter']
    mod = mod.drop_vars('time_counter')
    mod = mod.assign_coords({'time_counter': timestamp})
    return mod


def cd_neutral_10m(zw10):
    rgt33 = 0.5 + 0.5 * np.sign(zw10 - 33)  # should be an ndarray for * element-wise multiplication
    cd_neutral = 10**-3 * ((1 - rgt33) * ( 2.7 / zw10 + 0.142 + zw10/13.09 - 3.14807*10**-10*zw10**6) +rgt33*2.34)
    dragcoeff = cd_neutral  #10m neutral drag coefficient
    
    return dragcoeff



def calc_psi_h(pta):
    X2 = np.sqrt(np.abs(1. - 16.*pta))
    #X2[X2 < 1] = 1
    X2 = np.maximum(X2,np.array(1))
    X = np.sqrt(X2)
    stabit = 0.5 + 0.5*np.sign(pta)
    psi_h = -5.*pta*stabit + (1. - stabit)*(2.*np.log((1. + X2)*0.5))
    return psi_h



def calc_psi_m(pta):
    rpi = np.pi 
    X2 = np.sqrt(np.abs(1. - 16.*pta))
    X2 = np.maximum(X2, np.array(1))
    X = np.sqrt(X2)
    stabit = 0.5 + 0.5*np.sign(pta)
    psi_m = -5.*pta*stabit + (1. - stabit)*(2.*np.log((1. + X)*0.5) + np.log((1. + X2)*0.5) - 2.*np.arctan(X) + rpi*0.5)
    return psi_m



def COREBULK(zt,zu,sst,T_zt,q_sat,q_zt,dU):
    # zt is height of air temperature/humiidty  in meters
    # zu is height of air speed in meters
    # sst is sea surfac temperature in Kelvin
    # T_zt is potential air temperature in kelvin
    # q_sat is sea surface specific humidity in kg/kg
    # q_zt is specific air humidity  in kg/kg
    # dU is relative wind module at wind height, m/s
    
    # change cftime to timestamps
    sst = assign_timestamp(T_zt, sst)
    q_sat = assign_timestamp(q_zt,q_sat)


    iterations = 5
    vkarmn = 0.4
    grav  = 9.80665

    if (abs(zt-zu)<0.1):
        Temp_Wind_Same_Height=1

    else:
        Temp_Wind_Same_Height=0

    U_zu=np.maximum(dU,np.array(10))
    
    # stability check
    ztmp0=T_zt*(1+0.608*q_zt)-sst*(1+0.608*q_sat)
    stab=0.5+0.5*np.sign(ztmp0)
    
    ztmp0 = cd_neutral_10m(U_zu)
    sqrt_Cd_n10=np.sqrt(ztmp0)
    Ce_n10  = 10**-3*( 34.6 * sqrt_Cd_n10 )
    Ch_n10  = 1.e-3*sqrt_Cd_n10*(18*stab + 32.7*(1 - stab))
    
    # transfer coefficients first guesses
    Cd = ztmp0  
    Ce = Ce_n10
    Ch = Ch_n10 
    sqrt_Cd = sqrt_Cd_n10 

    T_zu = T_zt 
    q_zu = q_zt
    count=0
    for i in range(1, iterations+1):
        ztmp1 = T_zu - sst
        ztmp2 = q_zu - q_sat
        
        ztmp1 = Ch / sqrt_Cd * ztmp1
        ztmp2 = Ce / sqrt_Cd * ztmp2
        ztmp0 = T_zu * (1 + 0.608*q_zu)
        ztmp0 = (vkarmn * grav / ztmp0 * (ztmp1 * (1 + 0.608 * q_zu) + 0.608*T_zu * ztmp2)) / (Cd * U_zu * U_zu)
        
        #stability parameters
        zeta_u = zu*ztmp0
        zeta_utemp=zeta_u
        zeta_u = np.sign(zeta_u) * np.minimum(np.abs(zeta_u), 10)
        #zeta_utemp = np.maximum(zeta_utemp, np.array(10)) #zeta_utemp.where(zeta_utemp > 10, 10)
        #zeta_utemp[zeta_utemp>10]=10    
        #zeta_u = zeta_utemp * np.sign(zeta_u)
        #p1 = plt.pcolormesh(zeta_u[0], cmap='cmo.thermal')
        #plt.colorbar(p1)
        #plt.show() 
        zpsi_h_u = calc_psi_h(zeta_u)
        zpsi_m_u = calc_psi_m(zeta_u)

        if Temp_Wind_Same_Height==0:
            zeta_t = zt*ztmp0
            zeta_temp = np.minimum(np.abs(zeta_t),10)
            zeta_t = np.sign(zeta_t) * zeta_temp
            stab = np.log(zu/zt) - zpsi_h_u + calc_psi_h(zeta_t)
            T_zu = T_zt + ztmp1 / vkarmn * stab
            q_zu = q_zt + ztmp2 / vkarmn * stab
            q_zu = np.maximum(np.array(0), q_zu)
        
        input1 = U_zu / (1 + sqrt_Cd_n10 / vkarmn * (np.log(zu/10) - zpsi_m_u))
        input1 = np.maximum(input1, np.array(0.5))
        ztmp0=input1
        ztmp0 = cd_neutral_10m(ztmp0)
        sqrt_Cd_n10 = np.sqrt(ztmp0)
        Ce_n10 = 1e-3 * (34.6 * sqrt_Cd_n10)
        stab = 0.5 + 0.5 * np.sign(zeta_u)
        Ch_n10 = 1.e-3 * sqrt_Cd_n10 * (18 * stab + 32.7*(1 - stab))
        ztmp1 = 1 + sqrt_Cd_n10 / vkarmn * (np.log(zu/10) - zpsi_m_u)
        Cd = ztmp0 / (ztmp1 * ztmp1)
        sqrt_Cd = np.sqrt(Cd)
        ztmp0 = (np.log(zu/10) - zpsi_h_u) / vkarmn / sqrt_Cd_n10
        ztmp2 = sqrt_Cd / sqrt_Cd_n10
        ztmp1 = 1 + Ch_n10 * ztmp0
        Ch = Ch_n10 * ztmp2 / ztmp1
        ztmp1 = 1 + Ce_n10 * ztmp0
        Ce = Ce_n10 * ztmp2 / ztmp1
        
    return Cd,Ch,Ce,T_zu,q_zu
    


def TurbulentFluxes(inputmask, AtmData, runid, year):
    if runid == 'EPM151':
        nemo_path = '/project/6007519/pmyers/ANHA4/ANHA4-EPM151-S/'
        data_path = '/project/6007519/hlouis/DATA/'
        AirHeight = 2
        WindHeight = 10
        precipmult = 1000
        snow = 0
        
        precip_file='precip_gdps_ANHA4_y'
        precip_var = 'precip'
        temp_file = 't2_gdps_ANHA4_y'
        temp_var = 'Tair'
        
        humid_file = 'q2_gdps_ANHA4_y'
        humid_var = 'Qair'
        
        u_wind_file = 'u10_gdps_ANHA4_y'
        u_wind_var = 'u_wind'
        v_wind_file = 'v10_gdps_ANHA4_y'
        v_wind_var = 'v_wind'

    # constants
    albedo = 0.066
    Lv = 2.5e6  # latent heat of vaporization [J/kg]
    lfus = 0.334e6  # latent heat of fusion [J/kg]
    cp = 1000.5  # specific heat capacity of air [J/kg/C]
    rcp = 3991.86795711963  # heat capacity
    cpic = 2067  # specific heat for ice
    Stef = 5.67e-8
    rhoa = 1.22  # air density
    zcoef_qsatw = 0.98 * 640380 / rhoa  # surface specific humidity coeffiecient? 

    t=73    
    # open the mask file
    mf = nc.Dataset(inputmask)
    rmask = mf['tmask']
    rmask = np.broadcast_to(rmask,(t,)+rmask.shape)

 
    # load on the Atmos Forcing data
    t2_files = sorted(glob.glob(data_path+temp_file+str(year)+'.nc'))
    q2_files = sorted(glob.glob(data_path+humid_file+str(year)+'.nc'))
    u10_files = sorted(glob.glob(data_path+u_wind_file+str(year)+'.nc'))
    v10_files = sorted(glob.glob(data_path+v_wind_file+str(year)+'.nc'))
    
    TempOnNemo = xr.open_mfdataset(t2_files)  # atmos forcing temp regridded to the NEMO grid
    TempOnNemo = TempOnNemo[temp_var]
   
    HumidOnNemo = xr.open_mfdataset(q2_files)
    HumidOnNemo = HumidOnNemo[humid_var]
    
    u10 = xr.open_mfdataset(u10_files)
    v10 = xr.open_mfdataset(v10_files)    
    u = u10[u_wind_var]
    v = v10[v_wind_var] 
    
    # load NEMO temp, u, and v data
    mdl_files_t = sorted(glob.glob(nemo_path+'ANHA4-'+runid+'_y'+str(year)+'*_gridT.nc')) 
    mdl_files_u = sorted(glob.glob(nemo_path+'ANHA4-'+runid+'*_y'+str(year)+'*_gridU.nc')) 
    mdl_files_v = sorted(glob.glob(nemo_path+'ANHA4-'+runid+'*_y'+str(year)+'*_gridV.nc')) 

    dt = xr.open_mfdataset(mdl_files_t)
    du = xr.open_mfdataset(mdl_files_u)
    dv = xr.open_mfdataset(mdl_files_v)
    dt = dt.rename({'x_grid_T': 'x', 'y_grid_T': 'y', 'deptht': 'depth'})  # only gridT has this weird naming convention 
    du = du.rename({'depthu': 'depth'})
    dv = dv.rename({'depthv': 'depth'})
    
    dt = dt.assign(sst=dt['votemper'] + 273.15)

    SST = dt['sst'][:,0]  #dt['votemper'][:,0]  # we just want the surface temp
    SST = assign_timestamp(TempOnNemo, SST)

    humid_mod = zcoef_qsatw * np.exp(-5107.4/SST)  # equation to calculate humidity from 


    u_mod = du['vozocrtx'][:,0]  # we just want the surface velocity
    v_mod = dv['vomecrty'][:,0]
    u_mod = assign_timestamp(u, u_mod)  # need to get the time stamp out of cftime to pandas timestamp
    v_mod = assign_timestamp(v, v_mod)
   
    # calculate the wind modulus by taking the difference between ocean and atmosphere velocity
    u_diff = u - u_mod
    v_diff = v - v_mod
    wndm = np.sqrt(u_diff**2+v_diff**2)
    
    # calculate the COREBULK formulae
    Cd,Ch,Ce,T_zu,q_zu = COREBULK(AirHeight,WindHeight,SST,TempOnNemo,humid_mod,HumidOnNemo,wndm)
    
    # now, let's apply the mask to the coefficients
    Cd = Cd.where(rmask==1)  #np.ma.masked_where(rmask==0, Cd)
    Ch = Ch.where(rmask==1)  #np.ma.masked_where(rmask==0, Ch)
    Ce = Ce.where(rmask==1)  #np.ma.masked_where(rmask==0, Ce)
    T_zu = T_zu.where(rmask==1)  #np.ma.masked_where(rmask==0, T_zu)
    q_zu = q_zu.where(rmask==1)  #np.ma.masked_where(rmask==0, q_zu)

    if AirHeight == 10:
        evap = rhoa * Ce * (humid_mod - HumidOnNemo) * wndm
        turbsens  = cp*rhoa*Ch*(SST-TempOnNemo) * wndm
    
    else:
        evap = rhoa * Ce * (humid_mod - q_zu) * wndm
        turbsens = cp * rhoa * Ch * (SST - TempOnNemo) * wndm  #np.nanmean(np.nanmean(cp * rhoa * Ch * (SST-T_zu)*wndm))

    turblat = evap*Lv  #np.nanmean(np.nanmean(evap*Lv))
    turblat = turblat.rename('latent')
    turbsens = turbsens.rename('sensible')
    
    #turblat = turblat.groupby('time_counter.month').mean('time_counter')
     
    print(turblat)
    turblat.to_netcdf('/project/6007519/hlouis/data_files/TurbFluxes/EPM151_TurbLatent_y'+str(year)+'.nc')
    turbsens.to_netcdf('/project/6007519/hlouis/data_files/TurbFluxes/EPM151_TurbSensible_y'+str(year)+'.nc')
    print('done!')
    dt.close()
    du.close()
    dv.close()
    mf.close()

region_mask = '/project/6007519/hlouis/scripts/masks/shbjb_mask.nc'

#TurbulentFluxes(region_mask, 'CGRF', runid='EPM151', year=2003)
#TurbulentFluxes(region_mask, 'CGRF', runid='EPM151', year=2004) 
#TurbulentFluxes(region_mask, 'CGRF', runid='EPM151', year=2005) 
#TurbulentFluxes(region_mask, 'CGRF', runid='EPM151', year=2006) 
#TurbulentFluxes(region_mask, 'CGRF', runid='EPM151', year=2007) 
#TurbulentFluxes(region_mask, 'CGRF', runid='EPM151', year=2008) 
#TurbulentFluxes(region_mask, 'CGRF', runid='EPM151', year=2009) 
#TurbulentFluxes(region_mask, 'CGRF', runid='EPM151', year=2010) 
#TurbulentFluxes(region_mask, 'CGRF', runid='EPM151', year=2011) 
#TurbulentFluxes(region_mask, 'CGRF', runid='EPM151', year=2012) 
#TurbulentFluxes(region_mask, 'CGRF', runid='EPM151', year=2013) 
#TurbulentFluxes(region_mask, 'CGRF', runid='EPM151', year=2014) 
#TurbulentFluxes(region_mask, 'CGRF', runid='EPM151', year=2015) 
#TurbulentFluxes(region_mask, 'CGRF', runid='EPM151', year=2016) 
#TurbulentFluxes(region_mask, 'CGRF', runid='EPM151', year=2017) 
#TurbulentFluxes(region_mask, 'CGRF', runid='EPM151', year=2018) 
#TurbulentFluxes(region_mask, 'CGRF', runid='EPM151', year=2019) 
#TurbulentFluxes(region_mask, 'CGRF', runid='EPM151', year=2020) 
#TurbulentFluxes(region_mask, 'CGRF', runid='EPM151', year=2021) 
TurbulentFluxes(region_mask, 'CGRF', runid='EPM151', year=2022) 
