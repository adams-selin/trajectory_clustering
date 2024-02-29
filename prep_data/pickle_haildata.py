#reading in hail data, write out all trajectories greater than some threshold to a pickle file
#plot final hail size at initial location and final location

import numpy as np
import xarray as xr
import glob
import pickle

low_thresholds = [15, 25, 50.]
hi_thresholds = [19, 45, 150.]

rstnums = np.arange(19,37)

for rstnum in rstnums:

    rstnumd = str(rstnum).zfill(2)
    print ('')
    print('rst'+rstnumd)


    #open in hail trajectory data, only keep defined values
    filename = '/glade/derecho/scratch/radams/20120529_zshape/rst0000'+rstnumd+'/cm1out_haildata.nc'
    dum = xr.open_dataset(filename)
    ds = dum.where(dum.x>=0.)
    dum.close

    #first, make sure there was at least some moisture where the embryo started
    incloud = (ds.qliq.isel(time=[0]).values+ds.qice.isel(time=[0]).values) > 0.
    incloud = incloud[0,:] #get rid of extra dimension
    keep = xr.DataArray(incloud, dims=('xh'),coords={'xh':ds.xh})
    ds = ds.where(keep,drop=True)
    #change diameter to mm
    ds['d'] = ds.d*1.E3 #convert to mm

    print ('remaining in-cloud embryos ', ds.xh.shape[0])

    #handle some hailstones that are just getting stuck at the same spot.
    #there's probably a better way to do this, but ARGH xarray sucks sometimes
    #are the last two times the exact same?
    same = ((ds.x.isel(time=slice(-10,-1)).values == 
             ds.x.isel(time=slice(-11,-2)).values) & 
            (ds.y.isel(time=slice(-10,-1)).values == 
             ds.y.isel(time=slice(-11,-2)).values) & 
            (ds.z.isel(time=slice(-10,-1)).values ==
             ds.z.isel(time=slice(-11,-2)).values)).sum(axis=0)
    same = same==9 #whre are all these values exactly the same the whole duration (9 timesteps)? 
    keep = xr.DataArray(~same, dims=('xh'),coords={'xh':ds.xh})
    ds = ds.where(keep,drop=True)
    
    print ('removed stuck hailstones ', ds.xh.shape[0])

    #store this total array, and loop through some different thresholds
    allds = ds

    for low_threshold, hi_threshold in zip(low_thresholds,hi_thresholds):
        print (low_threshold, hi_threshold)
        #start with all the trajectories
        newds = allds
        print ('starting with: ', newds.xh.shape[0])

        #first filter by only keeping hail that is bigger than threshold at any point
        big = ((newds.d.max(dim='time', skipna=True) >= low_threshold) & 
            (newds.d.max(dim='time', skipna=True) < hi_threshold) )
        bigds = newds.where(big,drop=True)

        print ('first trim. hail trajectories remaining: ', bigds.xh.shape[0])

        #now, do an additional trimming to just those hailstones that reach the surface
        # larger than the threshold
        hail_d = bigds.d.values
        end_ind = (~np.isnan(hail_d)).cumsum(0).argmax(0)
        traj_ind = np.arange(0,hail_d.shape[1])
        endd = hail_d[end_ind,traj_ind]
        keep = xr.DataArray(((endd >= low_threshold) & (endd < hi_threshold)), 
                            dims=('xh'),coords={'xh':bigds.xh})
        endds = bigds.where(keep, drop=True)

        print ('done trimming by size. hail trajectories remaining: ', endds.xh.shape[0])

        #write out the trajectory data
        hail_x = endds.x.values
        hail_y = endds.y.values
        hail_z = endds.z.values
        hail_d = endds.d.values
        hail_mass = 1E3*endds.dense*1.33333*3.1415927*(0.5*endds.d*1.E-3)**3. #g
        hail_mass = hail_mass.values
        hail_time = endds.time*1E-9
        hail_time = np.array(hail_time.values,dtype='float')
        #round the time from immediately after the restart file down two seconds
        hail_time[0] = np.round(hail_time[0],-2)

        directory = '/'.join(filename.split('/')[:-1])
        with open(directory+'/haildata_ge'+str(int(low_threshold))+
                'lt'+str(int(hi_threshold))+'.pkl','wb') as f:
            pickle.dump(hail_x,f,pickle.HIGHEST_PROTOCOL)
            pickle.dump(hail_y,f,pickle.HIGHEST_PROTOCOL)
            pickle.dump(hail_z,f,pickle.HIGHEST_PROTOCOL)
            pickle.dump(endds.u.values,f,pickle.HIGHEST_PROTOCOL)
            pickle.dump(endds.v.values,f,pickle.HIGHEST_PROTOCOL)
            pickle.dump(endds.w.values,f,pickle.HIGHEST_PROTOCOL)
            pickle.dump(hail_d,f,pickle.HIGHEST_PROTOCOL)
            pickle.dump(endds.dense.values,f,pickle.HIGHEST_PROTOCOL)
            pickle.dump(endds.tv.values,f,pickle.HIGHEST_PROTOCOL)
            pickle.dump(endds.ts.values,f,pickle.HIGHEST_PROTOCOL)
            pickle.dump(endds.fw.values,f,pickle.HIGHEST_PROTOCOL)
            pickle.dump(endds.dice.values,f,pickle.HIGHEST_PROTOCOL)
            pickle.dump(endds.qice.values,f,pickle.HIGHEST_PROTOCOL)
            pickle.dump(endds.qliq.values,f,pickle.HIGHEST_PROTOCOL)
            pickle.dump(endds.tc.values,f,pickle.HIGHEST_PROTOCOL)
            pickle.dump(hail_mass,f,pickle.HIGHEST_PROTOCOL)
            pickle.dump(hail_time,f,pickle.HIGHEST_PROTOCOL)


        print ('done writing out trajectory pickle file for ge'+str(int(low_threshold))+
                'lt'+str(int(hi_threshold)))


