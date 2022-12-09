import matplotlib
from mpl_toolkits.mplot3d import axes3d
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
from matplotlib import cm
from matplotlib import collections  as mc
import glob
import sys


dir_prefix=sys.argv[1]
file_prefix=sys.argv[2]
eps=sys.argv[3]
minlns = sys.argv[4]



short_file_prefix = file_prefix
long_file_prefix = file_prefix+'_'+eps+'_'+minlns



#set some predefined figure boundaries
if dir_prefix.find('4d') >= 0:
    if dir_prefix.find('left') >= 0:
        xmin = -10000
        ymin = -25000
        numxs = 160
        numys = 160
    else:
        xmin = -8000
        ymin = -10000
        numxs = 100
        numys = 100

elif dir_prefix.find('lauren') >= 0:
    xmin = -10000
    ymin = -10000
    numxs = 80
    numys = 80

else:  #3d steady state files
    numxs=80#160
    numys=80#160
    if dir_prefix.find('left') >= 0:
        if (time == '19'):
            xmin=50000 #32500 #30000
            xmax=70000#47500 #50000
            ymin=62000#37500 #30000
            ymax=82000#50000
        if (time == '25'):
            xmin=57000 #32500 #30000
            ymin=64000#37500 #30000
    else:
        if (time == '19'):
            xmin=32500 #30000
            ymin=32000 #30000
        elif time == '24':
            xmin=23000
            ymin=23000 #30000
        elif time == '25':
            xmin=20000
            ymin=18000 #30000
        else:
            xmin=27500 #30000
            ymin=27000 #30000

xmax=xmin+numxs*250
ymax = ymin+numys*250


print ('opening '+dir_prefix+'/subtraj/'+short_file_prefix+'_subtraj.txt')
allsubs = np.genfromtxt(dir_prefix+'/subtraj/'+short_file_prefix+'_subtraj.txt',
    usecols=(0,1,2,3,4,5,6),names=('sx','sy','sz','ex','ey','ez','parent_traj'))
#start off with all the parent trajectories as "noise"
noise = allsubs['parent_traj']
print ('Total number of parent trajectories: ', np.unique(allsubs['parent_traj']).shape[0])


#create the x,y coordinates by hand
xh = xmin + np.arange(numxs)*250
yh = ymin + np.arange(numys)*250




#Plot the merged clusters one at a time
cluster_files = glob.glob(dir_prefix+'/merged_clusters/'+long_file_prefix+'_merged_[0-9]*.txt')
cluster_files.sort()
#cluster_files = cluster_files#[start:end]
if len(cluster_files) <= 10:
    colors = ['C0','C1','C2','C3','C4','C5','C6','C7','C8','C9']
    colors = colors[0:len(cluster_files)]
else:
    #colors = plt.cm.jet(np.linspace(0,1,len(cluster_files)))
    colors = plt.cm.nipy_spectral(np.linspace(0.1,0.9,len(cluster_files)))

for i, file in enumerate(cluster_files):
    #start the figure to plot, add surface reflectivity contours
    fig = plt.figure(figsize=(11,9))
    ax = fig.add_subplot(111, projection='3d')
    fig.subplots_adjust(left=0.04,right=0.99,top=0.94,bottom=0.089)
    ax.set_xlim(xmin,xmax)
    ax.set_ylim(ymin,ymax)
    ax.set_zlim(0,12000)

    cluster_num = file.split('_')[-1].split('.')[0]

    print (cluster_num)
    cluster = np.genfromtxt(file,skip_header=1,usecols=(0,1,2,3,4,5,6),
        names=('sx','sy','sz','ex','ey','ez','parent_traj'))
    clusterx = zip(cluster['sx'],cluster['ex'])
    clustery = zip(cluster['sy'],cluster['ey'])
    clusterz = zip(cluster['sz'],cluster['ez'])

    #plot the cluster files
    for x, y, z in zip(clusterx, clustery, clusterz):
        ax.plot(x,y,z,color=colors[i])

    ax.plot([],[],[],color=colors[i],label='Cluster '+cluster_num)



    ax.set_xlim(xmin,xmax)
    ax.set_ylim(ymin,ymax)
    ax.set_zlim(0,12000)
    ax.legend(loc='best')
    ax.set_xlabel('X (m)')
    ax.set_ylabel('Y (m)')
    ax.set_zlabel('Z (m)')
    ax.set_title('Cluster '+cluster_num+'with Eps: '+eps+'  MinLns: '+minlns)

    fig.tight_layout()
    ax.view_init(elev=10,azim=-70)


    plt.savefig(dir_prefix+'/plots/'+long_file_prefix+'_merged_clusters_'+cluster_num+'_bytime.png')
    plt.close()

    #remove this cluster's parent trajectories from the noise array
    noise = np.setdiff1d(noise,cluster['parent_traj'])



#plot the remaining line segments as noise
noise_match = np.isin(allsubs['parent_traj'],noise)
#zip to get in format for plotting line segments
allx = zip(allsubs['sx'][noise_match],allsubs['ex'][noise_match])
ally = zip(allsubs['sy'][noise_match],allsubs['ey'][noise_match])
allz = zip(allsubs['sz'][noise_match],allsubs['ez'][noise_match])

fig = plt.figure(figsize=(11,9))
ax = fig.add_subplot(111, projection='3d')
fig.subplots_adjust(left=0.04,right=0.99,top=0.94,bottom=0.089)
ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin,ymax)
ax.set_zlim(0,12000)

print ('Noise')
print ('Total number of unclustered parent trajectories: ',
    np.unique(allsubs['parent_traj'][noise_match]).shape[0])

for x, y, z in zip(allx, ally, allz):
    ax.plot(x,y,z,color='C7')

ax.plot([],[],[],color='C7',label='Noise')


ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin,ymax)
ax.set_zlim(0,12000)
ax.legend(loc='best')
ax.set_xlabel('X (m)')
ax.set_ylabel('Y (m)')
ax.set_zlabel('Z (m)')
ax.set_title('Noise with Eps: '+eps+'  MinLns: '+minlns)

fig.tight_layout()
ax.view_init(elev=10,azim=-70)

plt.savefig(dir_prefix+'/plots/'+long_file_prefix+'_merged_clusters_noise_bytime.png')
plt.close()
