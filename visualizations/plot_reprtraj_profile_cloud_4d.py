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
#useful_clusters_groupings = [['D','E','G','A']]
# num_sc_choices = 13

#right mover
short_dir_prefix='/Users/rselin/Documents/NSF_hail/streamlined_hail_scripts/right-mover'
eps='500'
minlns='3'
dir_prefix='/Users/rselin/Documents/NSF_hail/streamlined_hail_scripts/right-mover/4d'
short_file_prefix = 'cm1out_haildata_traj_rightW'
useful_clusters_groupings = [['B','A','I','H'],
                             ['C','I','J','F'],
                             ['D','E','G']]
useful_clusters_groupings = [['A','B','C']]
num_sc_choices = 10

#lauren
# dir_prefix='/Users/rselin/Documents/NSF_hail/streamlined_hail_scripts/right-mover/lauren'
# eps='300'
# minlns='20'
# short_file_prefix = 'hailtraj_Wrelative_fgt45mm'
# label = 'Hailstones >= 45 mm'
# subtitle = 'Kingfisher'
# useful_clusters_groupings = [['C','F'],
#                              ['G','D',],
#                              ['G','E',],
#                              ['A','B']]
# num_sc_choices = 7


sc_colors = np.linspace(0.1,0.9,num_sc_choices)
supercluster_names = string.ascii_uppercase


file_prefix = short_file_prefix+'_'+eps+'_'+minlns

for useful_clusters in useful_clusters_groupings:

    repr_traj_files = [dir_prefix+'/repr_traj/'+file_prefix+'_traj_bytime_'+i+'.txt' for i in useful_clusters]

    #char_names = ('x','y','z','d','dense','ts','fw','vt','ri','rw','tc','w','rate')
    #char_names = ('parent','x', 'y', 'z', 'sec', 'd', 'dense', 'ts', 'fw', 'vt', 'ri',
    #       'rw', 'tc', 'w', 'growth_rate', 'massrate', 'mass')
    char_names =  ('parent','x','y','z','sec','d','dense','ts','fw','vt','ri','rw','tc','w','growth_rate','massrate','mass')
    plot_variables = ('w','tc','ri','rw')

    usecols=np.arange(1,18,1,dtype='int')


    #set up panel for hailstone vertical profiles of this trajectory
    fig = plt.figure(num=None, figsize=(8,12), dpi=100)
    fig.subplots_adjust(hspace=0.15,wspace=0.04,left=0.10,right=0.99,top=0.99,bottom=0.07)
    ax = plt.subplot(4,1,1)
    bx = plt.subplot(4,1,2)
    cx = plt.subplot(4,1,3)
    dx = plt.subplot(4,1,4)
    axes = (ax,bx,cx,dx)


    labelfont = 12
    tickfont = 10
    error_int = 20
    #colors = plt.cm.jet(np.linspace(0,1,len(repr_traj_files)))
    #colors = plt.cm.Set1(np.linspace(0,1,9)[0:len(repr_traj_files)])

    #use these command if you want to match the cluster colors on the hailswath plot
    sc_indices = [supercluster_names.find(j) for j in useful_clusters]
    colors = plt.cm.nipy_spectral(sc_colors[sc_indices])

    #and this command if you want to just have reasonably separated values
    # sc_indices = np.arange(0,len(useful_clusters))
    # colors = np.array(['C0','C1','C2','C4','C5','C6','C7','C8','C9','C10'])
    # colors = colors[sc_indices]
    #
    usecols=np.arange(1,18,1,dtype='int')


    for i, repr_traj_file in enumerate(repr_traj_files):
        print ('analyzing ', repr_traj_file)
        #char_file = char_files[i]
        char_file = repr_traj_file

        num_pts = np.genfromtxt(repr_traj_file,max_rows=1,unpack=True).astype('int')

        chars_min = np.genfromtxt(char_file,skip_header=1,names=char_names,
            max_rows=num_pts,usecols=usecols,delimiter=',')
        chars_25p = np.genfromtxt(char_file,skip_header=1+num_pts,names=char_names,
            max_rows=num_pts,usecols=usecols,delimiter=',')
        chars_mean = np.genfromtxt(char_file,skip_header=1+num_pts*2,names=char_names,
            max_rows=num_pts,usecols=usecols,delimiter=',')
        chars_75p = np.genfromtxt(char_file,skip_header=1+num_pts*3,names=char_names,
            max_rows=num_pts,usecols=usecols,delimiter=',')
        chars_max = np.genfromtxt(char_file,skip_header=1+num_pts*4,names=char_names,
            max_rows=num_pts,usecols=usecols,delimiter=',')
        chars_time = np.genfromtxt(char_file,skip_header=1,names=('time'),
            max_rows=num_pts,usecols=(0),delimiter=',',dtype=None,encoding=None)

        #any weird missing ts values?
        chars_min['ts'][chars_min['ts']<100] = np.nan
        chars_min['dense'][chars_min['dense']<100] = np.nan

        times = np.array( [datetime.datetime.strptime(i[0],'%H:%M:%S.%f') for i in chars_time])
        times = times - times[-1]
        secs = np.array([i.total_seconds() for i in times])

        #make the plot
        for variable, axis in zip(plot_variables, axes):
            if (variable == 'rw') | (variable == 'ri'):
                chars_25p[variable] = chars_25p[variable]*1.E3
                chars_75p[variable] = chars_75p[variable]*1.E3
                chars_mean[variable] = chars_mean[variable]*1.E3

            axis.fill_between(secs, chars_25p[variable], chars_75p[variable], alpha=0.2,
                color=colors[i])
            axis.plot(secs, chars_mean[variable], label='Cluster '+useful_clusters[i],
                color=colors[i])
            #adding error bars showing the max and min points of the distribution
            # axis.errorbar(secs[::error_int], chars_mean[variable][::error_int],
            #     yerr=[chars_mean[variable][::error_int]-chars_min[variable][::error_int],
            #           chars_max[variable][::error_int]-chars_mean[variable][::error_int]],fmt='o',
            #     color=colors[i], alpha=0.5)

    for axis in axes:
        axis.set_xlabel('')
        axis.set_xticklabels('')
        axis.set_xlim([-1500,0])

    dx.set_xticks(np.arange(-1500,0,300))
    dx.set_xticklabels(np.arange(-25,0,5).astype('int'),size=tickfont)
    dx.set_xlabel('Minutes before reaching sfc', size=labelfont)
    ax.legend(loc= 'best')
    ax.set_ylabel('Updraft speed (m s$^{-1}$)',size=labelfont)
    bx.set_ylabel('Cloud temperature (K)',size=labelfont)
    cx.set_ylabel('Cloud ice (g kg$^{-1}$)',size=labelfont)
    dx.set_ylabel('Cloud water (g kg$^{-1}$)',size=labelfont)

    allcluster_string = ''.join(useful_clusters)

    plt.savefig(dir_prefix+'/plots/'+file_prefix.split('/')[-1]+
        '_trajbytime_profile_cloud_'+allcluster_string+'.png')
    plt.close()
