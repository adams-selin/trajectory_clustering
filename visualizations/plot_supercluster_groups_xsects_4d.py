#Script to plot updraft and clustered trajectory data in x-y, x-z, and y-z cross-sections
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import netCDF4 as nc4
import glob
import sys
from matplotlib import cm
import pandas as pd
import plot_xsect_functions
from os.path import exists
import string
from matplotlib.path import Path
from matplotlib.patches import PathPatch


#left mover
# short_dir_prefix='/Users/rselin/Documents/NSF_hail/streamlined_hail_scripts/left-mover'
# eps='800'
# minlns='100'
# dir_prefix='/Users/rselin/Documents/NSF_hail/streamlined_hail_scripts/left-mover/4d'
# short_file_prefix = 'cm1out_haildata_traj_leftWrelative'
# #useful_clusters_groupings = [['D','E','G','A']]
# label = 'Left mover'
# useful_clusters_groupings = [['G','D'],
#                              ['E','A'],
#                              ['C','H','I','M'],
#                              ['F','K'],
#                              ['J']]
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
label = 'Right mover'
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
#

file_prefix = short_file_prefix+'_'+eps+'_'+minlns
supercluster_names = string.ascii_uppercase
sc_colors = np.linspace(0.1,0.9,num_sc_choices)

for useful_clusters in useful_clusters_groupings:

    repr_traj_files = [dir_prefix+'/repr_traj/'+file_prefix+'_traj_bytime_'+i+'.txt' for i in useful_clusters]
    char_names =  ('parent','x','y','z','sec','d','dense','ts','fw','vt','ri','rw','tc','w','growth_rate','massrate','mass')
    plot_variables = ('d','massrate','mass','dense','fw')
    usecols=np.arange(1,18,1,dtype='int')

    # sc_indices = [supercluster_names.find(j) for j in useful_clusters]
    # colors = plt.cm.nipy_spectral(sc_colors[sc_indices])

    if dir_prefix.find('4d') >= 0:
        if dir_prefix.find('left') >= 0:
            if np.intersect1d(useful_clusters, np.array(['A','B','D','E','G','L'])).shape[0] \
                >= np.intersect1d(useful_clusters, np.array(['C','F','H','I','J','M'])).shape[0]:
                numxs = 80
                numys = 80
                xmin = -10000
                ymin = -15000
            elif np.intersect1d(useful_clusters, np.array(['A','B','D','E','G','L'])).shape[0] \
                < np.intersect1d(useful_clusters, np.array(['C','F','H','I','J','M'])).shape[0]:
                numxs = 120
                numys = 120
                xmin = -5000
                ymin = -25000
            else:
                numxs = 160
                numys = 160
                xmin = -10000
                ymin = -25000
        else:
            numxs = 100
            numys = 100
            xmin = -10000
            ymin = -10000
    else:
        numxs = 60
        numys = 60
        xmin = -7500
        ymin = -7500

    xmax=xmin+numxs*250
    ymax = ymin+numys*250
    zmax = 12000

    # 1) plot dimensions
    xh = xmin + np.arange(numxs)*250
    yh = ymin + np.arange(numys)*250
    zh = [0, 0.25, 0.5, 0.7500001, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3,
        3.25, 3.5, 3.75, 4, 4.25, 4.5, 4.75, 5, 5.25, 5.5, 5.75, 6, 6.25, 6.5,
        6.75, 7, 7.25, 7.5, 7.75, 8, 8.25, 8.5, 8.75, 9, 9.25, 9.5, 9.75, 10,
        10.25, 10.5, 10.75, 11, 11.25, 11.5, 11.75, 12, 12.25, 12.5, 12.75, 13,
        13.25, 13.5, 13.75, 14, 14.25, 14.5, 14.75, 15, 15.25, 15.5, 15.75, 16,
        16.25, 16.5, 16.75, 17, 17.25, 17.5, 17.75, 18, 18.25, 18.5, 18.75, 19,
        19.25, 19.5, 19.75, 20] * 10000
    zh = np.array(zh[0:49])

    xmeshxy,ymeshxy = np.meshgrid(xh,yh)
    xmeshxz,zmeshxz = np.meshgrid(xh,zh)
    ymeshyz,zmeshyz = np.meshgrid(yh,zh)

    if dir_prefix.find('left') >= 0:
        xticks = np.arange(xmin,xmax,5000)
        yticks =  np.arange(ymin,ymax,5000)
        zticks = np.arange(0,12000,2000)
    else:
        xticks = np.arange(xmin,xmax,5000)
        yticks =  np.arange(ymin,ymax,5000)
        zticks = np.arange(0,12000,2000)

    xticklabels = np.array(xticks/1E3,dtype='int')
    yticklabels = np.array(yticks/1E3,dtype='int')
    zticklabels =  np.array(zticks/1E3,dtype='int')


    fig = plt.figure(figsize=(11,11))
    #see https://matplotlib.org/stable/gallery/subplots_axes_and_figures/gridspec_multicolumn.html#sphx-glr-gallery-subplots-axes-and-figures-gridspec-multicolumn-py
    gs = GridSpec(3,3, figure=fig)
    gs.update(wspace=0.02,hspace=0.02,left=0.05,right=0.99,bottom=0.05,top=0.94)
    ax = fig.add_subplot(gs[1:,0:-1])  #x-y plot
    bx = fig.add_subplot(gs[0,0:-1]) #x-z plot
    cx = fig.add_subplot(gs[1:,-1])  #y-z plot


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

    ax.set_xlim(xmin,xmax)
    ax.set_ylim(ymin,ymax)
    ax.axvline(x=0,linewidth=1,color='k')
    ax.axhline(y=0,linewidth=1,color='k')
    bx.set_xlim(xmin,xmax)
    bx.set_ylim(0,zmax)
    bx.axvline(x=0,linewidth=2,color='k')
    cx.set_xlim(0,zmax)
    cx.set_ylim(ymin,ymax)
    cx.axhline(y=0,linewidth=2,color='k')



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

        #borrowed code from https://stackoverflow.com/questions/54255797/in-matplotlib-how-do-i-fill-between-two-curves-defined-by-two-different-sets-of
        #doesn't work the best, especially when curving in a circle or crossing itself
        ax.fill(np.append(chars_25p['x'],chars_mean['x'][::-1]),
                np.append(chars_25p['y'],chars_mean['y'][::-1]), alpha=0.2,
                color=colors[i])
        ax.fill(np.append(chars_mean['x'],chars_75p['x'][::-1]),
                np.append(chars_mean['y'],chars_75p['y'][::-1]), alpha=0.2,
                color=colors[i])
        ax.plot(chars_mean['x'],chars_mean['y'],label='Cluster '+useful_clusters[i],
            color=colors[i],zorder=5)

        bx.fill(np.append(chars_25p['x'],chars_mean['x'][::-1]),
                np.append(chars_25p['z'],chars_mean['z'][::-1]), alpha=0.2,
                color=colors[i])
        bx.fill(np.append(chars_mean['x'],chars_75p['x'][::-1]),
                np.append(chars_mean['z'],chars_75p['z'][::-1]), alpha=0.2,
                color=colors[i])
        bx.plot(chars_mean['x'],chars_mean['z'],color=colors[i],zorder=5)

        cx.fill(np.append(chars_25p['z'],chars_mean['z'][::-1]),
                np.append(chars_25p['y'],chars_mean['y'][::-1]), alpha=0.2,
                color=colors[i])
        cx.fill(np.append(chars_mean['z'],chars_75p['z'][::-1]),
                np.append(chars_mean['y'],chars_75p['y'][::-1]), alpha=0.2,
                color=colors[i])
        cx.plot(chars_mean['z'],chars_mean['y'],color=colors[i],zorder=5)

    ax.set_xlabel('X (km)')
    ax.set_ylabel('Y (km)')
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels,fontsize=labelfont)
    ax.set_yticks(yticks)
    ax.set_yticklabels(yticklabels,fontsize=labelfont)
    ax.set_xlim(xmin,xmax)
    ax.set_ylim(ymin,ymax)
    ax.legend(loc='best')

    bx.set_ylabel('Z (km)')
    bx.set_yticks(zticks)
    bx.set_yticklabels(zticklabels,fontsize=labelfont)
    bx.set_xticks(xticks)
    bx.set_xticklabels('')
    bx.set_xlim(xmin,xmax)
    bx.set_ylim(0,zmax)

    cx.set_xlabel('Z (km)')
    cx.set_xticks(zticks)
    cx.set_xticklabels(zticklabels,fontsize=labelfont)
    cx.set_yticks(yticks)
    cx.set_yticklabels('')
    cx.set_xlim(0,zmax)
    cx.set_ylim(ymin,ymax)

    # cbar_ax1 = fig.add_axes([0.7,0.66,0.05,0.34])
    # cbar_ax1.set_axis_off()
    # fig.colorbar(im, ax=cbar_ax1, ticks = qclvls[::2],
    #     label='Cloud water (g kg$^{-1}$)',shrink=1.25)
    # cbar_ax2 = fig.add_axes([0.85,0.66,0.05,0.34])
    # cbar_ax2.set_axis_off()
    # fig.colorbar(matplotlib.cm.ScalarMappable(cmap=cmap, norm=mass_norm),
    #     ax=cbar_ax2, ticks=mass_bins[::10], label="Hailstone mass growth rate (g s$^{-1}$)",
    #     shrink=1.25)

    allcluster_string = ''.join(useful_clusters)

    fig.suptitle(label+', Eps: '+eps+'  MinLns: '+minlns)
    #fig.text(0.75,0.7,'Updraft Speed (m s$^{-1}$)')
    plt.savefig(dir_prefix+'/plots/'+file_prefix+'_xsect_group_'+allcluster_string+'.png')
    plt.close()
