import numpy as np
import netCDF4 as nc4
import matplotlib.pyplot as plt
import glob
import pickle
import xarray as xr
from scipy import interpolate


#file to read in the left, right updraft movements, interpolate them to 30 s timesteps,
# and "convert" the hail trajectories to a left or right updraft-relative framework.

#read in the left and right updraft movements
#leftx, lefty, wtime = np.genfromtxt('left_filterwmax.txt',unpack=True)
rightx, righty, wtime = \
    np.genfromtxt('/glade/derecho/scratch/radams/20120529/clean_run/right_filterwmax.txt',
                                      unpack=True)
# These can be easily interpolated using the scipy.interpolate.interp1d command. 
# NOTE: numpy's interp command does NOT extrapolate!

#cycle through each restart run's trajectory pickle file
rstnums = np.arange(19,38)
#for each threshold grouping...
low_thresholds = [15, 25, 50.]
hi_thresholds = [19, 45, 150.]

for rstnum in rstnums:

    rstnumd = str(rstnum).zfill(3)
    print ('')
    print('rst'+rstnumd)

    filename = '/glade/derecho/scratch/radams//20120529_zshape/rst000'+rstnumd+'/cm1out_haildata.nc'
    directory = '/'.join(filename.split('/')[:-1])

    for low_threshold, hi_threshold in zip(low_thresholds,hi_thresholds):

        print (low_threshold, hi_threshold)

        with open(directory+'/haildata_ge'+str(int(low_threshold))+
                'lt'+str(int(hi_threshold))+'.pkl','rb') as f:
            hail_x = pickle.load(f)
            hail_y = pickle.load(f)
            hail_z = pickle.load(f)
            hail_u = pickle.load(f)
            hail_v = pickle.load(f)
            hail_w = pickle.load(f)
            hail_d = pickle.load(f)
            hail_dense = pickle.load(f)
            hail_tv = pickle.load(f)
            hail_ts = pickle.load(f)
            hail_fw = pickle.load(f)
            hail_dice = pickle.load(f)
            hail_qice = pickle.load(f)
            hail_qliq = pickle.load(f)
            hail_tc = pickle.load(f)
            hail_mass = pickle.load(f)
            hail_time = pickle.load(f)

        #convert this whole mess to an xarray dataset
        ds = xr.Dataset(
            data_vars=dict(
                x=(['time', 'traj'], hail_x),
                y=(['time', 'traj'], hail_y),
                z=(['time', 'traj'], hail_z),
                u=(['time', 'traj'], hail_u),
                v=(['time', 'traj'], hail_v),
                w=(['time', 'traj'], hail_w),
                d=(['time', 'traj'], hail_d),
                dense=(['time', 'traj'], hail_dense),
                tv=(['time', 'traj'], hail_tv),
                ts=(['time', 'traj'], hail_ts),
                fw=(['time', 'traj'], hail_fw),
                dice=(['time', 'traj'], hail_dice),
                qice=(['time', 'traj'], hail_qice),
                qliq=(['time', 'traj'], hail_qliq),
                tc=(['time', 'traj'], hail_tc),
                mass=(['time', 'traj'], hail_mass),
            ),
            coords=dict(
                time=(['time'], hail_time),
            ),
            attrs=dict(description='hail trajectories from restart '+rstnumd),
        )

        """
        #Following is useful if you have two supercell's worth of trajectories
        #find where the final y value is north or south of 57 km; that's the dividing line
        # b/w left and right mover.
        finald = ds.d.max(dim='time', skipna=True)
        finaly = np.nanmax(ds.y.where(ds.d == finald).values,axis=0) #keep just that value, not all the nans
        leftmover = finaly > 57000
        rightmover = finaly < 57000    
        leftds = ds.isel(traj=leftmover)
        rightds = ds.isel(traj=rightmover)


        #now...to convert the left mover x, y locations from absolute to updraft relative.
        #first, establish the x, y location at each ds timestep.
        fx = interpolate.interp1d(wtime, leftx*1.E3, fill_value='extrapolate')
        fy = interpolate.interp1d(wtime, lefty*1.E3, fill_value='extrapolate')
        all_leftx = fx(ds.time)
        all_lefty = fy(ds.time)
        # all_leftx = np.interp(ds.time,wtime,leftx) - np.interp doesn't extrapolate!
        # all_lefty = np.interp(ds.time,wtime,lefty)
        #repeat to the correct shape
        dum = np.expand_dims(all_leftx,axis=1)
        all_leftx = np.repeat(dum,leftds.x.shape[1],axis=1)
        dum = np.expand_dims(all_lefty,axis=1)
        all_lefty = np.repeat(dum,leftds.x.shape[1],axis=1)
        #convert to m from km
        all_leftx = all_leftx * 1.E3
        all_lefty = all_lefty * 1.E3
        
        #remove from the ds x and y
        leftds['x'] = leftds.x - all_leftx
        leftds['y'] = leftds.y - all_lefty

        #also want to convert the u, v motions to storm-relative.
        dtime = ds.time.values.astype('float')
        uspd = (all_leftx[1:,0]-all_leftx[0:-1,0]) / (dtime[1:] - dtime[0:-1])
        vspd = (all_lefty[1:,0]-all_lefty[0:-1,0]) / (dtime[1:] - dtime[0:-1])
        #smooth these a bit. 5-point running mean.
        smoothu = np.append(uspd[0:2],np.convolve(uspd,np.ones(5), 'valid')/5)
        smoothu = np.append(smoothu,uspd[-3:])
        smoothv = np.append(vspd[0:2],np.convolve(vspd,np.ones(5), 'valid')/5)
        smoothv = np.append(smoothv,vspd[-3:])
        #repeat to the correct shape
        dum = np.expand_dims(smoothu,axis=1)
        smoothu = np.repeat(dum,leftds.x.shape[1],axis=1)
        dum = np.expand_dims(smoothv,axis=1)
        smoothv = np.repeat(dum,leftds.x.shape[1],axis=1)

        #remove from the ds u and v
        leftds['u'] = leftds.u - smoothu
        leftds['v'] = leftds.v - smoothv
        """

        #Use this section if you just have one supercell's worth...
        #make sure there aren't any trajectories that haven't reached the surface
        #I'm using this b/c I've run the simulations for 2.5 hours apiece after embryo insertion
        if np.isnan(ds.x.isel(time=[-1]).values).sum() != ds.traj.shape[0]:
            print ('problem with trajectory not reaching surface!')
            print ((~np.isnan(ds.x.isel(time=[-1]).values)).nonzero())
            reached_surface = np.isnan(ds.x.isel(time=[-1]).values)[0,:]
            keep = xr.DataArray(reached_surface, dims=('traj'),coords={'traj':ds.traj})
            ds = ds.where(keep,drop=True)


        #make sure it isn't empty!
        if ds.traj.shape[0] > 0:

            rightds = ds
            #now...repeat for the right mover.
            #first, establish the x, y location at each ds timestep.
            fx = interpolate.interp1d(wtime, rightx*1.E3, fill_value='extrapolate')
            fy = interpolate.interp1d(wtime, righty*1.E3, fill_value='extrapolate')
            all_rightx = fx(ds.time)
            all_righty = fy(ds.time)
            # all_rightx = np.interp(ds.time,wtime,rightx) - np.interp doesn't extrapolate!
            # all_righty = np.interp(ds.time,wtime,righty)
            #repeat to the correct shape
            dum = np.expand_dims(all_rightx,axis=1)
            all_rightx = np.repeat(dum,rightds.x.shape[1],axis=1)
            dum = np.expand_dims(all_righty,axis=1)
            all_righty = np.repeat(dum,rightds.x.shape[1],axis=1)
            
            #remove from the ds x and y
            rightds['x'] = rightds.x - all_rightx
            rightds['y'] = rightds.y - all_righty

            #also want to convert the u, v motions to storm-relative.
            dtime = ds.time.values.astype('float')
            uspd = (all_rightx[1:,0]-all_rightx[0:-1,0]) / (dtime[1:] - dtime[0:-1])
            vspd = (all_righty[1:,0]-all_righty[0:-1,0]) / (dtime[1:] - dtime[0:-1])
            #smooth these a bit. 5-point running mean.
            smoothu = np.append(uspd[0:2],np.convolve(uspd,np.ones(5), 'valid')/5)
            smoothu = np.append(smoothu,uspd[-3:])
            smoothv = np.append(vspd[0:2],np.convolve(vspd,np.ones(5), 'valid')/5)
            smoothv = np.append(smoothv,vspd[-3:])
            #repeat to the correct shape
            dum = np.expand_dims(smoothu,axis=1)
            smoothu = np.repeat(dum,rightds.x.shape[1],axis=1)
            dum = np.expand_dims(smoothv,axis=1)
            smoothv = np.repeat(dum,rightds.x.shape[1],axis=1)

            #remove from the ds u and v
            rightds['u'] = rightds.u - smoothu
            rightds['v'] = rightds.v - smoothv


            #fill NaNs with -9999 to match what traj_part_4d.f90 is expecting
            #leftds = leftds.fillna(-9999.)
            rightds = rightds.fillna(-9999.)
            #move time from ns to regular s
            #leftds['time'] = leftds.time * 1E9
            rightds['time'] = rightds.time * 1E9
            #transpose to get the dimension order traj_part_4d.f90 is expecting
            #leftds = leftds.transpose()
            rightds = rightds.transpose() 

            encoding = {"x": {"dtype": "float32", '_FillValue': -9999.0}, 
                        "y": {"dtype": "float32", '_FillValue': -9999.0}, 
                        "z": {"dtype": "float32", '_FillValue': -9999.0}, 
                        "u": {"dtype": "float32", '_FillValue': -9999.0}, 
                        "v": {"dtype": "float32", '_FillValue': -9999.0}, 
                        "w": {"dtype": "float32", '_FillValue': -9999.0}, 
                        "d": {"dtype": "float32", '_FillValue': -9999.0}, 
                        "dense": {"dtype": "float32", '_FillValue': -9999.0}, 
                        "tv": {"dtype": "float32", '_FillValue': -9999.0}, 
                        "ts": {"dtype": "float32", '_FillValue': -9999.0}, 
                        "fw": {"dtype": "float32", '_FillValue': -9999.0}, 
                        "dice": {"dtype": "float32", '_FillValue': -9999.0}, 
                        "qice": {"dtype": "float32", '_FillValue': -9999.0}, 
                        "qliq": {"dtype": "float32", '_FillValue': -9999.0}, 
                        "tc": {"dtype": "float32", '_FillValue': -9999.0}, 
                        "mass": {"dtype": "float32", '_FillValue': -9999.0}, 
                        "time": {"dtype": "float32", '_FillValue': None}
                        } 

            #now, just write everything back out to separate left/right netcdf files
            #leftds.to_netcdf(directory+'/cm1out_haildata_traj_leftWrelative.nc',
            #    encoding=encoding)
            rightds.to_netcdf(directory+'/haildata_ge'+str(int(low_threshold))+
            'lt'+str(int(hi_threshold))+'_Wrel.nc',
                encoding=encoding)

        else:
            print ('array is emtpy')






































