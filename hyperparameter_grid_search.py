import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
import glob
import sys
import subprocess, os
from matplotlib.gridspec import GridSpec

#script to search through a grid of hyperparameter values (eps and minlns) and
# calculate the resulting "entropy" and "noise factor"

# epss = np.arange(300,1501,100)
# minlnss = np.append(np.arange(3,20,2), np.arange(20,70,5))  #right mover
# dir_prefix = 'right-mover/4d'
# file_prefix = 'cm1out_haildata_traj_rightW'

# epss = np.arange(200,2401,200)
# minlnss = np.arange(3,20,2) 
# dir_prefix = '/glade/derecho/scratch/radams/20120529'
# file_prefix = '20120529_ge25lt45_Wrel'

#epss = np.arange(200,4401,400)
epss = np.arange(200,2401,200)
minlnss = np.arange(3,21,2) 
dir_prefix = '/glade/derecho/scratch/radams/20120529_zshape'
file_prefix = '20120529_zshape_ge25lt45_Wrel'


#Some theory...
# #p_x is the number of segments per cluster / sum off segments over all clusters
# # At edge cases, we are assuming these values are uniform.
# # For a small epsilon, number of segments per cluster is 1:
# p_x_small_e = 1/total_num_segments.
# # For a large epsilon, all sgements are in 1 cluster:
# p_x_large_e = total_num_segments / total_num_segments.
# So, resulting entropy equation can be...
# entropy = (p_x * np.log2(p_x))*total_num_segments


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

        #do the clustering
        print ('calling ','./cluster.exe',dir_prefix,file_prefix,
            str(eps), str(minlns))
        out_string = subprocess.check_output(['./cluster.exe',dir_prefix,
           file_prefix,str(eps), str(minlns)])
        seg_clusters = int(out_string.splitlines()[0].decode("utf-8").strip())
        seg_noise = int(out_string.splitlines()[1].decode("utf-8").strip())
        all_seg_clusters[ecount, mcount] = seg_clusters
        all_seg_noise[ecount, mcount] = seg_noise

        #do the merging
        print ('calling ','merge_clusters.py',dir_prefix,
            file_prefix, str(eps), str(minlns))
        subprocess.call(['python','merge_clusters.py',dir_prefix,
            file_prefix, str(eps), str(minlns)])

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

        #calculate our two factors, entropy and noise
        p_x = pt_per_cluster / (np.sum(pt_per_cluster))
        entropy[ecount, mcount] = np.sum(p_x * np.log2(p_x))
        noise[ecount, mcount] = num_noise_pt / total_parent_trajectories

        #once we've stored this, delete all the generated text files so it doesn't
        # fill up the computer
        # files = glob.glob(dir_prefix+'/cluster/'+long_file_prefix+'_cluster*.txt')
        # for file in files:
        #     os.remove(file)
        # files = glob.glob(dir_prefix+'/merged_clusters/'+long_file_prefix+'_merged*.txt')
        # for file in files:
        #     os.remove(file)

labelsize=14
ticklabelsize=12

#Fit these all on one plot
fig = plt.figure(figsize=(18,7))
gs = GridSpec(1,3, figure=fig)
gs.update(wspace=0.02,hspace=0.02,left=0.045,right=0.99,bottom=-0.07,top=0.94)
ax = fig.add_subplot(gs[0,0])
bx = fig.add_subplot(gs[0,1])
cx = fig.add_subplot(gs[0,2])

#Now we want to plot 2d plots of entropy, noise, and our combo factor
mark = np.where(entropy==entropy[np.isfinite(entropy)].min())
PB = ax.contourf(minlnss, epss, entropy, np.arange(-6.5,0.01,0.75),cmap=plt.cm.hot_r,extend='max')#, np.arange(0,1.01,0.1),extend='max')
#ax.plot(minlnss[mark[1][0]],epss[mark[0][0]],'bs',markersize=12)
ax.text(61.5,315,'(a)',color='w',fontsize=labelsize)
ax.set_xlabel('MinLns',size=labelsize)
ax.set_ylabel('eps',size=labelsize)
ax.set_title('Cluster Entropy, $H_c(X)$',size=labelsize)
ax.tick_params(axis='both', which='major', labelsize=ticklabelsize)
cbar = fig.colorbar(PB, ax=ax, orientation='horizontal',shrink = 0.9, pad=0.08)
cbar.ax.tick_params(labelsize=ticklabelsize)


mark = np.where(np.cos(noise*np.pi/2)==np.cos(noise*np.pi/2).max())
PB=bx.contourf(minlnss, epss, np.cos(noise*np.pi/2),np.arange(0.,1.01,0.1), cmap=plt.cm.hot)
#bx.plot(minlnss[mark[1][0]],epss[mark[0][0]],'bs',markersize=12)
bx.text(61.5,315,'(b)',color='w',fontsize=labelsize)
bx.set_xlabel('MinLns',size=labelsize)
#bx.set_ylabel('eps')
bx.set_yticklabels('')
bx.set_title('Noise factor, '+r'$\cos\left(\dfrac{Q\pi}{2}\right)$',size=labelsize)
bx.tick_params(axis='both', which='major', labelsize=ticklabelsize)
cbar = fig.colorbar(PB, ax=bx, orientation='horizontal',shrink = 0.9, pad=0.08)
cbar.ax.tick_params(labelsize=ticklabelsize)
#plt.savefig(dir_prefix+'/plots/hyperparameter_noise.png')
#plt.close()


combo = entropy * np.cos(noise*np.pi/2)
mark = np.where(combo == combo[np.isfinite(combo)].min())
PB=cx.contourf(minlnss, epss, combo, np.arange(-6.5,0.01,0.75), cmap=plt.cm.hot_r,extend='max')#, np.arange(0,1.01,0.1),extend='max')
#cx.plot(minlnss[mark[1][0]],epss[mark[0][0]],'bs',markersize=12)
cx.text(61.5,315,'(c)',color='w',fontsize=labelsize)
cx.set_xlabel('MinLns',size=labelsize)
#cx.set_ylabel('eps')
cx.set_yticklabels('')
cx.set_title('Parent trajectory entropy, $H_{pt}(X)$',size=labelsize)
cx.tick_params(axis='both', which='major', labelsize=ticklabelsize)
cbar = fig.colorbar(PB, ax=cx, orientation='horizontal',shrink = 0.9, pad=0.08)
cbar.ax.tick_params(labelsize=ticklabelsize)

plt.savefig(dir_prefix+'/plots/'+file_prefix+'_hyperparameter_all.png')
plt.close()


np.savetxt(dir_prefix+'/'+file_prefix+'_entropy.txt',entropy,delimiter=',')
np.savetxt(dir_prefix+'/'+file_prefix+'_noise.txt',noise,delimiter=',')
np.savetxt(dir_prefix+'/'+file_prefix+'_combo.txt',combo,delimiter=',')
