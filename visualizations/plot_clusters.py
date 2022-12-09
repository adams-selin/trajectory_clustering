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
file_prefix= sys.argv[2]
eps=sys.argv[3]
minlns = sys.argv[4]


#set some predefined figure boundaries
xmin=27500 #30000
xmax=47500 #50000
ymin=27000 #30000
ymax=47500

numxs=100
numys=100
xmax=xmin+numxs*250
ymax = ymin+numys*250


short_file_prefix = file_prefix
long_file_prefix = file_prefix+'_'+eps+'_'+minlns
print ('opening '+dir_prefix+'/subtraj/'+short_file_prefix+'_subtraj.txt')
allsubs = np.genfromtxt(dir_prefix+'/subtraj/'+short_file_prefix+'_subtraj.txt',
    usecols=(0,1,2,3,4,5),names=('sx','sy','sz','ex','ey','ez'))
#zip to get in format for plotting line segments
int = 1
allx = zip(allsubs['sx'][::int],allsubs['ex'][::int])
ally = zip(allsubs['sy'][::int],allsubs['ey'][::int])
allz = zip(allsubs['sz'][::int],allsubs['ez'][::int])

#start the figure to plot
fig = plt.figure(figsize=(11,9))
ax = fig.add_subplot(111, projection='3d')
fig.subplots_adjust(left=0.04,right=0.99,top=0.99,bottom=0.089)
ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin,ymax)
ax.set_zlim(0,12000)

#plot all the segments
for x, y, z in zip(allx, ally, allz):
    ax.plot(x,y,z,color='grey')


#Read in the cluster files
print ('opening '+dir_prefix+'/cluster*.txt files')
cluster_files = glob.glob(dir_prefix+'/cluster/'+long_file_prefix+'_cluster_[0-9]*.txt')
cluster_files.sort()
cluster_files = cluster_files#[start:end]
if len(cluster_files) <= 10:
    colors = ['C0','C1','C2','C3','C4','C5','C6','C8','C9']
    colors = colors[0:len(cluster_files)]
else:
    colors = plt.cm.jet(np.linspace(0,1,len(cluster_files)))

for i, file in enumerate(cluster_files):
    print (file)
    cluster_num = file.split('_')[-1].split('.')[0]
    cluster = np.genfromtxt(file,skip_header=1,usecols=(0,1,2,3,4,5),
        names=('sx','sy','sz','ex','ey','ez'))
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
ax.set_xlabel('X (m)')
ax.set_ylabel('Y (m)')
ax.set_zlabel('Z (m)')
ax.set_title('Eps: '+eps+'  MinLns: '+minlns,multialignment='center')

fig.tight_layout()
ax.view_init(elev=10,azim=-70)

plt.savefig(dir_prefix+'/plots/'+long_file_prefix+'_clusters.png')
plt.close()
