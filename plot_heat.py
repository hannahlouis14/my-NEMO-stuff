import glob
from datetime import datetime
import numpy as np
import numpy.ma as ma
import xarray as xr
import netCDF4 as nc
import matplotlib.pyplot as plt
import pandas as pd

###this is for running script backend###
import matplotlib
#matplotlib.use('Agg')
###----------------------------------###

runid = ['EPM151', 'EPM152']
long_names = {'EPM151': 'HYPE, CGRF', 'EPM152': 'HYPE, ERA'}

path = '/home/hlouis/projects/rrg-pmyers-ad/hlouis/plotting/data_files/heat_contents/'
fig_path = '/home/hlouis/projects/rrg-pmyers-ad/hlouis/plotting/figures/hc_figs/'

mask_file = '/project/6007519/weissgib/plotting/data_files/anha4_files/runoff_temp_regions_mask.nc'
#mf = nc.Dataset(mask_file)


regions = {'hb_mask': 'Hudson Bay'}
def read_mask(runids, long_name):
    mf = nc.Dataset(mask_file)
    experiment = []
    region = []
    hc = []
    date = []
    
    for exp in runids:
        heat_content_file = path+exp+'test_heat_content50m.nc'
        data_file = xr.open_dataset(heat_content_file)
        heat_content_data = data_file['votemper'].values
        times = data_file['time_counter'].values  # ice concentration uses ['year']
        t = data_file.dims['time_counter']  # ice concentration uses ['year'], type = ndarray
        data_file.close()  
        #print('total size:', np.size(heat_content_data))
        #print('nan count: ', np.count_nonzero(~np.isnan(heat_content_data)))
        times = [datetime.strptime(str(x),'%Y-%m-%d %H:%M:%S') for x in list(times)]  #ice concentration uses just "%Y")


        for r in regions.keys():
            rmask = mf[r][0,0]  # what does this [0, 0] do?
            
            rmask = np.broadcast_to(rmask, (t,)+rmask.shape)  # get the mask into the correct shape
            masked_heat = ma.masked_where(rmask==1, heat_content_data)
            #print('total  of masked ice:', np.size(masked_heat))
            masked_heat[masked_heat==0] = np.nan
            plt.pcolormesh(masked_heat[0])
            plt.colorbar()
            #plt.show()
            plt.savefig('/home/hlouis/projects/rrg-pmyers-ad/hlouis/test_map.png')
            plt.clf()
            regional_hc_ts = np.nansum(masked_heat, axis=(1,2))  # masked_heat.mean(axis=(1,2))  # averaging the heat content data over the masked area
            #print('nan count of masked ice:', np.count_nonzero(~np.isnan(masked_heat[14])))
            print('regional mean', np.shape(regional_hc_ts))
'''
            for i in range(t):
                if long_name:
                    experiment.append(long_names[exp])
                else:
                    experiment.append(exp)
                region.append(r)

            hc.extend(list(regional_hc_ts))
            date.extend(times)
    
    # make a pandas database for plotting
    all_data = {'experiment': experiment, 'region':region, 'hc':hc, 'date':date}
    df = pd.DataFrame(all_data)
    
    return df
'''


def hc_region_comp(runids, long_name=False):
    df = read_mask(runids,long_name)
    
    # now, do comparision of timeseries plots
    for r in regions:
        rd = df.loc[df['region'] == r]
        rd = rd.pivot(index='date', columns='experiment', values='hc')
        rd.plot() 
        plt.title(regions[r]+' Heat Content')
        plt.ylabel('heat content (J)')  #Temperature ($^o$C)')
        plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05), ncol=2, fancybox=True, shadow=True)
        plt.legend()
        plt.tight_layout()
        plt.savefig(fig_path+r+'_heat_content_sum_timeseries.png')
        plt.clf()

def hc_region_diff(runids, long_name=False):

    df = read_mask(runids, long_name)

    #plot the difference between two runs
    for r in regions:
        rd = df.loc[df['region'] == r]
        rd = rd.pivot(index='date', columns='experiment', values='hc')
        rd['diff'] = rd[runids[0]] - rd[runids[1]]
        rd['diff'].plot()
        plt.grid(True)
        plt.title(regions[r]+' Heat Content Difference')
        plt.ylabel('heat content (J)')
        #plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05),ncol=2, fancybox=True, shadow=True)
        plt.legend()
        plt.tight_layout()
        #plt.show()
        plt.savefig(fig_path+r+'_hc_50m_timeseries_diff.png')
        plt.clf()


#if __name__ == "__main__":
hc_region_comp(runids = ['EPM151', 'EPM152'])
    #hc_region_diff(runids=['EPM151', 'EPM151'])i


