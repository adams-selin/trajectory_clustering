import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
import glob
import sys
import subprocess, os
from matplotlib.gridspec import GridSpec

#script to read through the merged_clusters files created by 
# hyperparameter_grid_search and output # of parent trajectories per file/
# parent cluster

epss = np.arange(200,401,100)
minlnss = np.append(np.arange(3,20,2), np.arange(20,25,5))  #right mover
minlnss = np.append(np.array(1), minlnss)
dir_prefix = '../streamlined_hail_scripts/right-mover/lauren_new'
#file_prefix = 'cm1out_haildata_traj_rightW'
file_prefix = 'cntl_Wrelative_fgt45mm'

#first, pull in the total number of parent trajectories
allsubs = np.genfromtxt(dir_prefix+'/subtraj/'+file_prefix+'_subtraj.txt',
    usecols=(6),names=('parent_traj'))
print ('Total number of parent trajectories: ', np.unique(allsubs['parent_traj']).shape[0])
total_parent_trajectories = np.unique(allsubs['parent_traj']).shape[0]


#arrays to store output
entropy = np.zeros((epss.shape[0], minlnss.shape[0]))
noise = np.zeros((epss.shape[0], minlnss.shape[0]))
all_seg_clusters = np.zeros((epss.shape[0], minlnss.shape[0]))
all_seg_noise = np.zeros((epss.shape[0], minlnss.shape[0]))

for ecount, eps in enumerate(epss):
    for mcount, minlns in enumerate(minlnss):
        print ('eps: ', eps, 'minlns: ', minlns)
        long_file_prefix = file_prefix+'_'+str(eps)+'_'+str(minlns)
        #start off all parent trajs as noise, and remove as we match
        this_noise = allsubs['parent_traj']
        all_pt = np.unique(allsubs['parent_traj'])
        pt_per_cluster_list = []

        #calculate # of trajectories in each merged cluster
        cluster_files = glob.glob(dir_prefix+'/merged_clusters/'+long_file_prefix+
            '_merged_[0-9]*.txt')
        cluster_files.sort()
        pt_per_cluster = np.zeros(len(cluster_files)) #total unique parent trajs per cluster

        for i, file in enumerate(cluster_files):
            cluster = np.genfromtxt(file,skip_header=1,usecols=(6),
                names=('parent_traj'))
            pt_per_cluster[i] = np.unique(cluster['parent_traj']).shape[0]
            pt_per_cluster_list.append(np.unique(cluster['parent_traj']))

            #remove this cluster's parent trajectories from the noise array
            this_noise = np.setdiff1d(this_noise,cluster['parent_traj'])

        #what's left is noise
        num_noise_pt = np.unique(this_noise).shape[0]

        #print # of parent trajectories per parent cluster
        print (pt_per_cluster)
        print ('num parent clusters: ', len(pt_per_cluster))
        print ('noise: ', num_noise_pt)

