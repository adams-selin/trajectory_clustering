import matplotlib
matplotlib.use('agg')
import numpy as np
import matplotlib.pyplot as plt
import sys
import netCDF4 as nc4
import glob
import math
import datetime
from os.path import exists
import string

#left mover
# short_dir_prefix='/Users/rselin/Documents/NSF_hail/streamlined_hail_scripts/left-mover'
# eps='800'
# minlns='100'
# dir_prefix='/Users/rselin/Documents/NSF_hail/streamlined_hail_scripts/left-mover/4d'
# short_file_prefix = 'cm1out_haildata_traj_leftWrelative'
# useful_clusters_groupings = [['G','D'],
#                              ['E','A'],
#                              ['C','H','I','M'],
#                              ['F','K'],
#                              ['J']]
# num_sc_choices = 13

#right mover
# short_dir_prefix='/Users/rselin/Documents/NSF_hail/streamlined_hail_scripts/right-mover'
# eps='500'
# minlns='3'
# dir_prefix='/Users/rselin/Documents/NSF_hail/streamlined_hail_scripts/right-mover/4d'
# short_file_prefix = 'cm1out_haildata_traj_rightW'
# useful_clusters_groupings = [['B','A','I','H'],
#                              ['C','I','J','F'],
#                              ['D','E','G']]
# num_sc_choices = 10

#lauren
dir_prefix='/Users/rselin/Documents/NSF_hail/streamlined_hail_scripts/right-mover/lauren'
eps='300'
minlns='20'
short_file_prefix = 'hailtraj_Wrelative_fgt45mm'
label = 'Hailstones >= 45 mm'
subtitle = 'Kingfisher'
useful_clusters_groupings = [['C','F'],
                             ['G','D',],
                             ['G','E',],
                             ['A','B']]
num_sc_choices = 7


sc_colors = np.linspace(0.1,0.9,num_sc_choices)
supercluster_names = string.ascii_uppercase


file_prefix = short_file_prefix+'_'+eps+'_'+minlns

for useful_clusters in useful_clusters_groupings:

    num_traj_files = [dir_prefix+'/repr_traj/'+file_prefix+'_numtraj_bytime_'+i+'.txt' for i in useful_clusters]

    #char_names = ('x','y','z','d','dense','ts','fw','vt','ri','rw','tc','w','rate')
    #char_names = ('parent','x', 'y', 'z', 'sec', 'd', 'dense', 'ts', 'fw', 'vt', 'ri',
    #       'rw', 'tc', 'w', 'growth_rate', 'massrate', 'mass')
    #char_names =  ('parent','x','y','z','sec','d','dense','ts','fw','vt','ri','rw','tc','w','growth_rate','massrate','mass')
    #plot_variables = ('w','tc','ri','rw')

    #usecols=np.arange(1,18,1,dtype='int')


    #set up panel for hailstone vertical profiles of this trajectory
    fig = plt.figure(num=None, figsize=(8,3), dpi=100)
    fig.subplots_adjust(hspace=0.15,wspace=0.04,left=0.10,right=0.99,top=0.99,bottom=0.15)
    axis = plt.subplot(1,1,1)



    labelfont = 12
    tickfont = 10
    error_int = 20
    #colors = plt.cm.jet(np.linspace(0,1,len(num_traj_files)))
    #colors = plt.cm.Set1(np.linspace(0,1,9)[0:len(num_traj_files)])
    sc_indices = [supercluster_names.find(j) for j in useful_clusters]
    colors = plt.cm.nipy_spectral(sc_colors[sc_indices])

    #usecols=np.arange(1,18,1,dtype='int')


    for i, num_traj_file in enumerate(num_traj_files):
        print ('analyzing ', num_traj_file)
        #char_file = char_files[i]
        char_file = num_traj_file

        num_pts = np.genfromtxt(num_traj_file,max_rows=1,unpack=True).astype('int')

        num_trajs = np.genfromtxt(num_traj_file,skip_header=1,names=('num_traj'),
            max_rows=num_pts,usecols=(1),delimiter=',')
        chars_time = np.genfromtxt(char_file,skip_header=1,names=('time'),
            max_rows=num_pts,usecols=(0),delimiter=',',dtype=None,encoding=None)

        times = np.array( [datetime.datetime.strptime(i[0],'%H:%M:%S.%f') for i in chars_time])
        times = times - times[-1]
        secs = np.array([i.total_seconds() for i in times])

        #less than 3 = missing
        num_trajs['num_traj'][num_trajs['num_traj']<3] = np.nan

        #make the plot
        axis.plot(secs, num_trajs['num_traj'], label='Cluster '+useful_clusters[i],
            color=colors[i])

        axis.set_xlim([-1800,0])

    axis.set_xticks(np.arange(-1800,0,300))
    axis.set_xticklabels(np.arange(-30,0,5).astype('int'),size=tickfont)
    axis.set_xlabel('Minutes before reaching sfc', size=labelfont)
    axis.legend(loc= 'best')
    axis.set_ylabel('No. trajectories in cluster',size=labelfont)


    allcluster_string = ''.join(useful_clusters)

    plt.savefig(dir_prefix+'/plots/'+file_prefix.split('/')[-1]+
        '_trajbytime_profile_numtraj_'+allcluster_string+'.png')
    plt.close()
