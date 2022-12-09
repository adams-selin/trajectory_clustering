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


# dir_prefix='/Users/rselin/Documents/NSF_hail/streamlined_hail_scripts/left-mover/4d'
# eps='800'
# minlns='100'
# short_file_prefix = 'cm1out_haildata_traj_leftWrelative'
# #label = 'Hailstones >= 15 mm'
# label = 'Left mover'

# dir_prefix='/Users/rselin/Documents/NSF_hail/streamlined_hail_scripts/right-mover/4d'
# eps='500'
# minlns='3'
# short_file_prefix = 'cm1out_haildata_traj_rightW'
# #label = 'Hailstones >= 19 mm'
# label = 'Right mover'

# dir_prefix='/Users/rselin/Documents/NSF_hail/streamlined_hail_scripts/right-mover/lauren'
# eps='300'
# minlns='20'
# short_file_prefix = 'hailtraj_Wrelative_fgt45mm'
# label = 'Hailstones >= 45 mm'
# subtitle = 'Kingfisher'

short_dir_prefix='/Users/rselin/Documents/NSF_hail/streamlined_hail_scripts/right-mover'
eps='500'
minlns='3'
dir_prefix='/Users/rselin/Documents/NSF_hail/streamlined_hail_scripts/right-mover/4d'
short_file_prefix = 'cm1out_haildata_traj_rightW'
useful_clusters_groupings = [['B','A','I','H'],
                             ['C','I','J','F'],
                             ['D','E','G']]
useful_clusters_groupings = [['C','D','I','G']]
useful_cluster_groupings_labels = [['A','B','C','D']]
label = 'Trajectories initialized at 95, 170 min'
num_sc_choices = 10

long_file_prefix = short_file_prefix+'_'+eps+'_'+minlns

supercluster_names = string.ascii_uppercase
if dir_prefix.find('left') >= 0:
    supercluster_choices = {'A': ('14','19'), 'B': ('15',),
                            'C': ('16','17','18','23','24','36','43','44'),
                            'D': ('1',), 'E': ('3','20'), 'F': ('22','26','27','35'),
                            'G': ('2',), 'H': ('25','29','31','32','34'),
                            'I': ('41','45','52'),
                            'J': ('49','50','55','59','61','62','63','64'),
                            'K': ('4',), 'L': ('7','8'), 'M': ('56',)}
elif dir_prefix.find('lauren') >= 0:
    supercluster_choices = {'A': ('9','19','27'),
                            'B': ('38','39','5','71','86'),
                            'C': ('80','81','82','83','79','84','87'),
                            'D': ('42','43','44','66','50','72','69','68','56','55','61','26','88'),
                            'E': ('51','52','59','67','70','74','75'),
                            'F': ('29','85',),
                            'G': ('48','49','54','64','89')}
elif label.find('170') >= 0:
    supercluster_choices = {'C': ('22','23','30','18'),
                            'D': ('58','59'),
                            'I': ('4','6','61','62','7'),
                            'G': ('39','8')}
else:
   supercluster_choices = {'A': ('11','8','67','5','10'),
                           'B': ('1','2','20','25','3','12'),
                           'C': ('22','23','30','31','32','35','36','38','42','18'),
                           'F': ('51','50'), 'E': ('63','64','65'),
                           'D': ('17','40','46','58','59','60','57'),
                           'I': ('4','6','28','29','45','26','48'),
                           'G': ('37','39','55'), 'H': ('13','24','27'),
                           'J': ('43',)}

num_sc_choices = len(supercluster_choices)
# sc_colors_inds = np.linspace(0.1,0.9,num_sc_choices)
num_clusters = num_sc_choices
# #use these command if you want to match the cluster colors on the hailswath plot
# sc_colors = plt.cm.nipy_spectral(sc_colors_inds)
#and this command if you want to just have reasonably separated values
# sc_indices = np.arange(0,num_sc_choices)
# sc_colors = np.array(['C4','C5','C2','C0','C1','C6','C3','C7','C8','C9','C0','C1','C2','C4'])
# sc_colors = sc_colors[sc_indices]

#to match the supercluster colors
sc_colors_inds = np.linspace(0.1,0.9,10)
sc_colors = plt.cm.nipy_spectral(sc_colors_inds)

char_names = ('sec1','d1','dense1','ts1','fw1','vt1','ri1','rw1','tc1','w1',
              'sec2','d2','dense2','ts2','fw2','vt2','ri2','rw2','tc2','w2')

clustersx = np.zeros(num_clusters,dtype='object')
clustersy = np.zeros(num_clusters,dtype='object')
clustersz = np.zeros(num_clusters,dtype='object')
clusterex = np.zeros(num_clusters,dtype='object')
clusterey = np.zeros(num_clusters,dtype='object')
clusterez = np.zeros(num_clusters,dtype='object')
clustermass = np.zeros(num_clusters,dtype='object')




#for count, sc_name in enumerate(supercluster_names[0:num_sc_choices]):
for count, sc_name in enumerate(['I',]):
    print ('Supercluster '+sc_name)
    clusters = supercluster_choices.get(sc_name,())

    for i, cluster in enumerate(clusters):
        print ('Plotting '+cluster)
        cluster_file = dir_prefix+'/merged_clusters/'+long_file_prefix+'_merged_'+cluster+'.txt'
        info = np.genfromtxt(cluster_file,skip_header=1,usecols=(0,1,2,3,4,5),
            names=('sx','sy','sz','ex','ey','ez'))
        cluster_chars_file = dir_prefix+'/merged_clusters/'+long_file_prefix+'_merged_chars_'+cluster+'.txt'
        chars1 = np.genfromtxt(cluster_chars_file,skip_header=1,usecols=(0,1,2,3,4,5,6,7,8,9),
         names=char_names[0:10])
        chars2 = np.genfromtxt(cluster_chars_file,skip_header=1,usecols=(10,11,12,13,14,15,16,17,18,19),
         names=char_names[10:20])
        mass1 = 1E3*chars1['dense1']*1.33333*3.1415927*(0.5*chars1['d1']*1.E-3)**3. #g
        mass2 = 1E3*chars2['dense2']*1.33333*3.1415927*(0.5*chars2['d2']*1.E-3)**3. #g
        mean_time = (chars1['sec1'] + chars2['sec2']) * 0.5
        mass_rate = (mass2 - mass1)*1.E3 / mean_time

        supersx = info['sx']
        superex = info['ex']
        supersy = info['sy']
        superey = info['ey']
        supersz = info['sz']
        superez = info['ez']
        sc_mass_rate = mass_rate

        clustersx[count] = supersx
        clustersy[count] = supersy
        clustersz[count] = supersz
        clusterex[count] = superex
        clusterey[count] = superey
        clusterez[count] = superez
        clustermass[count] = sc_mass_rate


        #create the mass color scale for the trajectories
        cmap = plt.cm.viridis_r
        max_mass = np.around(max([i.max() for i in clustermass]),1)-0.1
        mass_bins = np.linspace(0,max_mass,81)
        mass_color = cmap(np.linspace(0,1,81))
        mass_norm = matplotlib.colors.BoundaryNorm(mass_bins,cmap.N)
        labelsize = 10

        thiscount = supercluster_names.find(sc_name)


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

        ax.set_xlim(xmin,xmax)
        ax.set_ylim(ymin,ymax)
        ax.axvline(x=0,linewidth=1,color='k')
        ax.axhline(y=0,linewidth=1,color='k')
        allx = zip(clustersx[count],clusterex[count])
        ally = zip(clustersy[count],clusterey[count])
        allz = zip(clustersz[count],clusterez[count])
        lastz = 0
        for i, (x, y, z) in enumerate(zip(allx, ally, allz)):
            if (i > 16700) and (i<= 16966):
                ax.plot(x,y,color=cmap(np.minimum(clustermass[count][i]/max_mass,1)),zorder=5)
                if lastz < 50:
                    ax.plot(x,y,'o',color=sc_colors[thiscount])
                lastz = z[-1]
        ax.set_xlabel('X (km)')
        ax.set_ylabel('Y (km)')
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticklabels,fontsize=labelsize)
        ax.set_yticks(yticks)
        ax.set_yticklabels(yticklabels,fontsize=labelsize)
        ax.set_xlim(xmin,xmax)
        ax.set_ylim(ymin,ymax)

        bx.set_xlim(xmin,xmax)
        bx.set_ylim(0,zmax)
        bx.axvline(x=0,linewidth=2,color='k')
        allx = zip(clustersx[count],clusterex[count])
        allz = zip(clustersz[count],clusterez[count])
        lastz = 0
        for i, (x, y) in enumerate(zip(allx, allz)):
            if (i > 16700) and (i<= 16966):
                bx.plot(x,y,color=cmap(np.minimum(clustermass[count][i]/max_mass,1)),zorder=5)
                if lastz < 50:
                    bx.plot(x,y,'o',color=sc_colors[thiscount])
                lastz = y[-1]
        bx.set_ylabel('Z (km)')
        bx.set_yticks(zticks)
        bx.set_yticklabels(zticklabels,fontsize=labelsize)
        bx.set_xticks(xticks)
        bx.set_xticklabels('')
        bx.set_xlim(xmin,xmax)
        bx.set_ylim(0,zmax)

        cx.set_xlim(0,zmax)
        cx.set_ylim(ymin,ymax)
        cx.axhline(y=0,linewidth=2,color='k')
        allz = zip(clustersz[count],clusterez[count])
        ally = zip(clustersy[count],clusterey[count])
        lastz = 0
        for i, (x, y) in enumerate(zip(allz, ally)):
            if (i > 16700) and (i<= 16966):
                cx.plot(x,y,color=cmap(np.minimum(clustermass[count][i]/max_mass,1)),zorder=5)
                if lastz < 50:
                    cx.plot(x,y,'o',color=sc_colors[thiscount])
                lastz = x[-1]
        cx.set_xlabel('Z (km)')
        cx.set_xticks(zticks)
        cx.set_xticklabels(zticklabels,fontsize=labelsize)
        cx.set_yticks(yticks)
        cx.set_yticklabels('')
        cx.set_xlim(0,zmax)
        cx.set_ylim(ymin,ymax)

        # cbar_ax1 = fig.add_axes([0.7,0.66,0.05,0.34])
        # cbar_ax1.set_axis_off()
        # fig.colorbar(im, ax=cbar_ax1, ticks = qclvls[::2],
        #     label='Cloud water (g kg$^{-1}$)',shrink=1.25)
        cbar_ax2 = fig.add_axes([0.85,0.66,0.05,0.34])
        cbar_ax2.set_axis_off()
        fig.colorbar(matplotlib.cm.ScalarMappable(cmap=cmap, norm=mass_norm),
            ax=cbar_ax2, ticks=mass_bins[::10], label="Hailstone mass growth rate (g s$^{-1}$)",
            shrink=1.25)

        fig.suptitle(label+', Cluster '+cluster+'\n'+
            'Eps: '+eps+'  MinLns: '+minlns,multialignment='center')
        #fig.text(0.75,0.7,'Updraft Speed (m s$^{-1}$)')
        plt.savefig(dir_prefix+'/plots/'+long_file_prefix+'_'+cluster+'_xsect_19th1000.png')
        plt.close()
