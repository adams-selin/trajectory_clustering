import numpy as np
import glob
import sys
from scipy.stats import mode

dir_prefix=sys.argv[1]
file_prefix=sys.argv[2]
# time=sys.argv[2]
# size=sys.argv[3]
# density=sys.argv[4]
eps=sys.argv[3]
minlns = sys.argv[4]
# dir_prefix =  '/Users/rselin/Documents/NSF_hail/streamlined_hail_scripts/left-mover'
# time ='19'
# size='5mm'
# density='500kg0p5dperp'
# eps='516'
# minlns='4'


# short_file_prefix='cm1out_0000'+time+'_'+size+'_'+density
# long_file_prefix = 'cm1out_0000'+time+'_'+size+'_'+density+'_'+eps+'_'+minlns
short_file_prefix = file_prefix
long_file_prefix = file_prefix+'_'+eps+'_'+minlns
#dir_prefix = dir_prefix+'/cm1out_0000'+time

# print (long_file_prefix)

cluster_files = glob.glob(dir_prefix+'/cluster/'+long_file_prefix+'_cluster_[0-9]*.txt')
cluster_files.sort()
cluster_num = np.array([int(file.split('.txt')[0].split('_')[-1]) for file in cluster_files])
num_clusters = len(cluster_files)

merged_parent_trajs = []
merged_parent_trajs_chars = []

#Loop through all cluster files, read in contents to one big array after appending cluster_num
# print ('read cluster files in to one big array')
for i,file in enumerate(cluster_files):
    cluster = np.genfromtxt(file,skip_header=1,usecols=6)
    num_subtrajs = cluster.shape[0]

    cluster_withnum = np.append(np.reshape(cluster,(num_subtrajs,1)),
        np.ones((num_subtrajs,1))*cluster_num[i],axis=1)
    if i==0:
        all_clusters = cluster_withnum
    else:
        all_clusters = np.append(all_clusters,cluster_withnum,axis=0)

#Loop through each parent trajectory within all clusters, find the clusters it contains
# print ('find clusters for each parent trajectory')
unique_parent_trajs = np.unique(all_clusters[:,0])
all_pt_clusters = []
for parent_traj in unique_parent_trajs:
    #find the clusters it contains (pt = parent trajectory)
    this_traj_inds = all_clusters[:,0] == parent_traj
    pt_clusters = np.unique(all_clusters[this_traj_inds,1])
    all_pt_clusters.append(pt_clusters)

all_pt_clusters = np.asarray(all_pt_clusters,dtype='object')


#Read in all subtrajectories (for output)
# print ('read in all subtrajectories')
subtraj_file = dir_prefix+'/subtraj/'+short_file_prefix+'_subtraj.txt'
subtraj_chars_file = dir_prefix+'/subtraj/'+short_file_prefix+'_subtraj_chars.txt'
subtrajs = np.genfromtxt(subtraj_file)
subtraj_chars = np.genfromtxt(subtraj_chars_file)


#Finally, loop through parent trajectories again, see where combinations have more
#than two clusters in common, and output those parent trajs to a merged cluster.
if num_clusters < 2:
    common_req = num_clusters
else:
    common_req = 2

#print ('Output parent trajectories with at least two clusters in common')
i=0
count=0
while unique_parent_trajs.size > 0:
    this_pt_clusters = all_pt_clusters[i]
    num_matching = np.array([np.intersect1d(this_pt_clusters,j).size for j in all_pt_clusters[i+1:]])
    pt_to_merge = num_matching >= common_req
    #print (this_pt_clusters, num_matching)

    parent_traj = unique_parent_trajs[i]
    #print ('Parent trajectory: ', parent_traj)
    if pt_to_merge.sum() > 0: #we found some that match
        parent_trajs_to_output = unique_parent_trajs[i+1:][pt_to_merge]
        #include current parent traj too
        parent_trajs_to_output = np.append(parent_traj,parent_trajs_to_output)

        #require at least three different parent trajectories to output as a
        #merged cluster
        if len(parent_trajs_to_output) > 2:
            #print ('Parent trajectories to output: ', parent_trajs_to_output)
            inds_to_output = np.isin(subtrajs[:,-1], parent_trajs_to_output)
            subtrajs_to_output = subtrajs[inds_to_output,:]
            subtraj_chars_to_output = subtraj_chars[inds_to_output,:]

            #output the data
            new_file = dir_prefix+'/merged_clusters/'+long_file_prefix+'_merged_'+\
                str(count+1)+'.txt'
            len_cluster = subtrajs_to_output.shape[0]
            np.savetxt(new_file,subtrajs_to_output,fmt=(7*'%10.1f'),header=str(len_cluster),
                comments='        ')

            new_chars_file = dir_prefix+'/merged_clusters/'+long_file_prefix+'_merged_chars_'+\
                str(count+1)+'.txt'
            np.savetxt(new_chars_file,subtraj_chars_to_output,fmt=(20*'%12.4f'),
                header=str(len_cluster),comments='        ')

            #increment the merged cluster count
            count = count +1

        #remove the parent trajectories already matched from the list
        inds_to_remove =  np.isin(unique_parent_trajs,parent_trajs_to_output)
        unique_parent_trajs = np.delete(unique_parent_trajs,inds_to_remove)
        all_pt_clusters = np.delete(all_pt_clusters,inds_to_remove)

        #increment the merged cluster count - now moved earlier, inside if statement above
        #count = count +1

    else:
        #remove the trajectory we just tested - it will be noise
        inds_to_remove = unique_parent_trajs==parent_traj
        unique_parent_trajs = np.delete(unique_parent_trajs,inds_to_remove)
        all_pt_clusters = np.delete(all_pt_clusters,inds_to_remove)
