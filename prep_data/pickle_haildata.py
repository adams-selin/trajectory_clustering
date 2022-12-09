#reading in hail data, write out all trajectories greater than some threshold to a pickle file
#plot final hail size at initial location and final location

import numpy as np
import netCDF4 as nc4
import xarray as xr
import glob
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib import cm
import pickle

rstnums = np.arange(1,20)

for rstnum in rstnums:

    rstnumd = str(rstnum).zfill(3)
    print('rst'+rstnumd)


    #open in hail trajectory data, only keep defined values
    filename = '/glade/scratch/radams/hail/20120529/rst'+rstnumd+'/cm1out_haildata.nc'
    dum = xr.open_mfdataset(filename)
    ds = dum.where(dum.x>=0.)
    dum.close

    #only keep that which is bigger than 19 mm
    ds['d'] = ds.d*1.E3 #convert to mm
    big = (ds.d.max(dim='time', skipna=True) >= 19)
    ds = ds.where(big,drop=True)

    print ('done trimming by size')

    #write out the trajectory data
    hail_x = ds.x.values
    hail_y = ds.y.values
    hail_z = ds.z.values
    hail_d = ds.d.values
    hail_mass = 1E3*ds.dense*1.33333*3.1415927*(0.5*ds.d*1.E-3)**3. #g
    hail_mass = hail_mass.values
    hail_time = ds.time*1E-9
    hail_time = hail_time.values

    directory = '/'.join(filename.split('/')[:-1])
    with open(directory+'/cm1out_haildata_traj.pkl','wb') as f:
        pickle.dump(hail_x,f,pickle.HIGHEST_PROTOCOL)
        pickle.dump(hail_y,f,pickle.HIGHEST_PROTOCOL)
        pickle.dump(hail_z,f,pickle.HIGHEST_PROTOCOL)
        pickle.dump(ds.u.values,f,pickle.HIGHEST_PROTOCOL)
        pickle.dump(ds.v.values,f,pickle.HIGHEST_PROTOCOL)
        pickle.dump(ds.w.values,f,pickle.HIGHEST_PROTOCOL)
        pickle.dump(hail_d,f,pickle.HIGHEST_PROTOCOL)
        pickle.dump(ds.dense.values,f,pickle.HIGHEST_PROTOCOL)
        pickle.dump(ds.tv.values,f,pickle.HIGHEST_PROTOCOL)
        pickle.dump(ds.ts.values,f,pickle.HIGHEST_PROTOCOL)
        pickle.dump(ds.fw.values,f,pickle.HIGHEST_PROTOCOL)
        pickle.dump(ds.dice.values,f,pickle.HIGHEST_PROTOCOL)
        pickle.dump(ds.qice.values,f,pickle.HIGHEST_PROTOCOL)
        pickle.dump(ds.qliq.values,f,pickle.HIGHEST_PROTOCOL)
        pickle.dump(ds.tc.values,f,pickle.HIGHEST_PROTOCOL)
        pickle.dump(hail_mass,f,pickle.HIGHEST_PROTOCOL)
        pickle.dump(hail_time,f,pickle.HIGHEST_PROTOCOL)


    print ('done writing out trajectory pickle file')


