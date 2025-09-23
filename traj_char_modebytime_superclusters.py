import numpy as np
import pandas as pd
import glob
import sys
from scipy.stats import mode
import datetime
import string

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# ! DESCRIPTION:  Given a list of files each containing a cluster of subtrajectory
# ! line segments, find the cluster characteristics (including mean repr traj) in
# ! time. Repeat for
# ! each cluster. Write the cluster characteristics out to its own file.
# !   1) Calculate "time before reaching surface" for each subtractory.
# !   2) Cycle through all times, starting from most negative,  and identify the
# !      parent trajectory closest to the "median" cluster location at each timestep.
# !      Keep a running tally of which parent trajectory is most often closest to the median.
# !   3) Return the most often identified parent trajectory as the representative trajectory.
# !
# ! INPUT: File, directory prefix for merged cluster files,  total number of cluster files,
# !   epsilon, and MinLns
# ! OUTPUT: Multiple files, each containing a representative trajectory for that cluster.
# !
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


# dir_prefix='left-mover/4d'
# eps='800'
# minlns='100'
# short_file_prefix = 'cm1out_haildata_traj_leftWrelative'

# dir_prefix='right-mover/4d'
# eps='500'
# minlns='3'
# short_file_prefix = 'cm1out_haildata_traj_rightW'

# dir_prefix='right-mover/lauren'
# eps='300'
# minlns='20'
# short_file_prefix = 'hailtraj_Wrelative_fgt45mm'

# dir_prefix='/Users/rselin/Documents/NSF_hail/streamlined_hail_scripts/20120529_sounding3'
# eps='1200'
# minlns='3'
# short_file_prefix = '20120529_ge15lt19_Wrel'
# label = 'Hailstones >= 15 mm, < 19 mm'

# dir_prefix='/Users/rselin/Documents/NSF_hail/streamlined_hail_scripts/20120529_sounding3'
# eps='800'
# minlns='11'
# short_file_prefix = '20120529_ge25lt45_Wrel'
# label = 'Hailstones >= 25 mm, < 45 mm'

dir_prefix='/Users/rselin/Documents/NSF_hail/streamlined_hail_scripts/20120529_sounding3'
eps='1000'
minlns='3'
short_file_prefix = '20120529_ge50lt150_Wrel'
label = 'Hailstones >= 50 mm'

long_file_prefix = short_file_prefix+'_'+eps+'_'+minlns
endtime = datetime.datetime(2012,5,30,0,0,0) + datetime.timedelta(seconds=7200) #dummy end time

#endtime = datetime.datetime(2012,5,29,int(time)%24,0,0) + datetime.timedelta(seconds=3600)

names = ('x1','y1','z1','x2','y2','z2','parent')
char_names = ('sec1','d1','dense1','ts1','fw1','vt1','ri1','rw1','tc1','w1','u1','v1',
              'sec2','d2','dense2','ts2','fw2','vt2','ri2','rw2','tc2','w2','u2','v2')
first = ('x1','y1','z1','sec1','d1','dense1','ts1','fw1','vt1','ri1','rw1','tc1','w1','u1','v1')
second = ('x2','y2','z2','sec2','d2','dense2','ts2','fw2','vt2','ri2','rw2','tc2','w2','u2','v2')
third =  ('x','y','z','sec','d','dense','ts','fw','vt','ri','rw','tc','w','u','v')

supercluster_names = string.ascii_uppercase
if dir_prefix.find('4d') >= 0:
    if dir_prefix.find('left') >= 0:
        supercluster_choices = {'A': ('14','19'), 'B': ('15',),
                                'C': ('16','17','18','23','24','36','43','44'),
                                'D': ('1',), 'E': ('3','20'), 'F': ('22','26','27','35'),
                                'G': ('2',), 'H': ('25','29','31','32','34'),
                                'I': ('41','45','52'),
                                'J': ('49','50','55','59','61','62','63','64'),
                                'K': ('4',), 'L': ('7','8'), 'M': ('56',)}
    else:
       supercluster_choices = {'A': ('11','8','67','5','10'),
                               'B': ('1','2','20','25','3','12'),
                               'C': ('22','23','30','31','32','35','36','38','42','18'),
                               'F': ('51','50'), 'E': ('63','64','65'),
                               'D': ('17','40','46','58','59','60','57'),
                               'I': ('4','6','28','29','45','26','48'),
                               'G': ('37','39','55'), 'H': ('24','27'),  #took out 13 - didn't really match
                               'J': ('43',), 'K': ('13',)}
elif dir_prefix.find('lauren') >= 0:
    supercluster_choices = {'A': ('9','19','27'),
                            'B': ('38','39','5','71','86'),
                            'C': ('80','81','82','83','79','84','87'),
                            'D': ('42','43','44','66','50','72','69','68','56','55','61','26','88'),
                            'E': ('51','52','59','67','70','74','75'),
                            'F': ('29','85',),
                            'G': ('48','49','54','64','89')}
elif dir_prefix.find('20120529') >= 0:
    if dir_prefix.find('zshape') >= 0:
        if short_file_prefix.find('ge15lt19') >= 0:
            supercluster_choices = {'A': ('1','2')}
                                    
        if short_file_prefix.find('ge25lt45') >= 0:
            supercluster_choices = {'A': ('21','26','40','43'),
                                    'B': ('7','8','10','13'),
                                    'C': ('27','32','42','49','50'),
                                    'D': ('22','29','56','57','60','61'),
                                    'E': ('20','24'),
                                    'F': ('1','2','3','5','6','12'),
                                    'G': ('4',),
                                    'H': ('47','51','53','54','55','59')}

        if short_file_prefix.find('ge50lt150') >= 0:
            supercluster_choices = {'A': ('1',), 
                                    'B': ('2',),
                                    'C': ('3',),
                                    'D': ('4','5')}

    else: #sounding 3
        if short_file_prefix.find('ge15lt19') >= 0:
            supercluster_choices = {'A': ('1','6','7','10','11','12','13','14'),
                                    'B': ('15','16','17','18','19'),
                                    'C': ('2','4','5','8','9'),
                                    'D': ('3',),
                                    'E': ('13',)}

        if short_file_prefix.find('ge25lt45') >= 0:
            # supercluster_choices = {'A': ('41','92',),
            #                         'B': ('1','55','72','74','82','83','84','90','95'),
            #                         'C': ('6','60','79','87','88','102'),
            #                         'D': ('7','8','35','77','89'),
            #                         'E': ('2','18'),
            #                         'F': ('9','13','28'),
            #                         'G': ('10','21','31','60'),
            #                         'H': ('97','98'),
            #                         'I': ('5','25','32','37'), 'M': ('22','26','76','81'),
            #                         'J': ('3','15'),
            #                         'K': ('107',),
            #                         'L': ('71',)}
            #merged these again after looking at xsect plots
            # supercluster_choices = {'A': ('41','92','6','60','79','87','88','102'),
            #                         'B': ('1','55','72','74','82','83','84','90','95'),
            #                         'D': ('7','8','35','77','89','22','26','76','81'),
            #                         'E': ('2','18'),
            #                         'F': ('9','13','28'),
            #                         'G': ('10','21','31','60','71'),
            #                         'H': ('97','98','107'),
            #                         'I': ('5','25','32','37'), 
            #                         'C': ('3','15')}
            #decided that was too much merging, back to just extra merging of H and K
            supercluster_choices = {'A': ('41','92'),}#,'6','60','79','87','88','102'),
                                    # 'B': ('1','55','72','74','82','83','84','90','95'),
                                    # #'C': ('6','60','79','87','88','102'),
                                    # 'C': ('7','8','35','77','89'),
                                    # 'D': ('2','18'),
                                    # 'E': ('9','13','28'),
                                    # 'F': ('10','21','31','60'),
                                    # 'G': ('97','98','107',),
                                    # 'H': ('5','25','32','37'), 'J': ('22','26','76','81'),
                                    # 'I': ('3','15'),
                                    # 'K': ('71',)}
        if short_file_prefix.find('ge50lt150') >= 0:
            supercluster_choices = {'A': ('1','5','8'),
                                    'B': ('3','6','7')}#, 'C': ('3',)}

else:
    print ('YOU HAVE THE WRONG SCRIPT. GO GET SOME REST.')

num_sc_choices = len(supercluster_choices)
for count, sc_name in enumerate(supercluster_names[0:num_sc_choices]):
    print ('Supercluster '+sc_name)
    clusters = supercluster_choices.get(sc_name,())

    for i, cluster in enumerate(clusters):
        file = dir_prefix+'/merged_clusters/'+long_file_prefix+'_merged_'+cluster+'.txt'
        char_file = dir_prefix+'/merged_clusters/'+long_file_prefix+'_merged_chars_'+cluster+'.txt'

        thistrajs = pd.read_csv(file,header=1,names=names,sep=' ',skipinitialspace=True)
        thischars = pd.read_csv(char_file,header=1,names=char_names,sep=' ',skipinitialspace=True)
        if i==0:
            trajs = thistrajs
            chars = thischars
            print (trajs.shape)
        else:
            trajs = pd.concat([trajs,thistrajs])
            chars = pd.concat([chars,thischars])
            print (trajs.shape)

    if len(clusters) > 0:
        combine = pd.concat([trajs,chars],axis=1)
        for a, b,c in zip(first, second, third):
            new = (combine[a] + combine[b])/2.
            combine[c] = new
            combine = combine.drop(columns=[a,b])

        combine['growth_rate'] = 100. * (chars['d2']-chars['d1'])/(chars['sec2']-chars['sec1']) #1E-2 mm/s
        mass1 = 1E3*chars['dense1']*1.33333*3.1415927*(0.5*chars['d1']*1.E-3)**3. #g
        mass2 = 1E3*chars['dense2']*1.33333*3.1415927*(0.5*chars['d2']*1.E-3)**3. #g
        combine['massrate'] = (mass2 - mass1) / (chars['sec2']-chars['sec1'])  #g/s
        combine['mass'] = (mass2 + mass1) / 2.

        parents = np.array(combine['parent'].values,dtype='int')
        alltrajs = list(combine.groupby('parent'))
        first_time_through = 1
        for dum in alltrajs:
            thistraj = dum[1]
            rsec = thistraj.loc[:,'sec'] - thistraj['sec'].values.max()
            td = pd.to_timedelta([i for i in rsec],unit='s') + endtime
            thistraj = thistraj.assign(rsec=td)
            #thistraj = thistraj.assign(rsec=rsec)
            thistraj = thistraj.set_index('rsec')

            #want to resample AND interpolate to every 30 s
            # empty frame with desired index. Borrowed code from
            # https://stackoverflow.com/questions/25234941/python-regularise-irregular-time-series-with-linear-interpolation
            rs = pd.DataFrame(index=thistraj.resample('30s').asfreq().iloc[1:].index)

            # array of indexes corresponding with closest timestamp after resample
            idx_after = np.searchsorted(thistraj.index.values, rs.index.values)

            # values and timestamp before/after resample
            for name in thistraj.columns:
                rs['after'] = thistraj.loc[thistraj.index[idx_after], name].values
                rs['before'] = thistraj.loc[thistraj.index[idx_after - 1], name].values
                rs['after_time'] = thistraj.index[idx_after]
                rs['before_time'] = thistraj.index[idx_after - 1]

                #calculate new weighted value
                rs['span'] = (rs['after_time'] - rs['before_time'])
                rs['after_weight'] = (rs['after_time'] - rs.index) / rs['span']
                # I got errors here unless I turn the index to a series
                rs['before_weight'] = (pd.Series(data=rs.index, index=rs.index) - rs['before_time']) / rs['span']

                rs[name] = rs.eval('after * before_weight + before * after_weight')


            #get rid of the extra columns
            rs = rs.drop(columns=['after','before','after_time','before_time','span','after_weight','before_weight'])
            thistraj = rs #resampled to every 30 sec

            if first_time_through == 1:
                alltrajs_rsec = thistraj
                first_time_through = 0
            else:
                alltrajs_rsec = pd.concat([alltrajs_rsec,thistraj])

        alltrajs_rsec = alltrajs_rsec.sort_index()

        times = alltrajs_rsec.index.unique()
        num_columns = alltrajs_rsec.shape[1]
        num_trajs_per_sc = []
        closest_to_median_trajs = []
        closest_to_min_trajs = []
        closest_to_max_trajs = []

        for time in times:
            try:
                dum = alltrajs_rsec.loc[time].shape[1] #check if only one record

                if alltrajs_rsec.loc[time].shape[0] >= 5:
                    diff_loc = alltrajs_rsec.loc[time,('x','y','z')] - \
                               alltrajs_rsec.loc[time,('x','y','z')].median()
                    diff_loc['parent'] = alltrajs_rsec.loc[time,'parent']
                    diff_loc['rmse_median'] = ( (diff_loc['x']**2. + diff_loc['y']**2. + diff_loc['z']**2.)/3.)**0.5
                    min_ptraj = diff_loc.iloc[diff_loc['rmse_median'].argmin(),3]  #parent is column 3
                    closest_to_median_trajs.append(int(min_ptraj))

                    diff_loc = alltrajs_rsec.loc[time,('x','y','z')] - \
                               alltrajs_rsec.loc[time,('x','y','z')].min()
                    diff_loc['parent'] = alltrajs_rsec.loc[time,'parent']
                    diff_loc['rmse_min'] = ( (diff_loc['x']**2. + diff_loc['y']**2. + diff_loc['z']**2.)/3.)**0.5
                    min_ptraj = diff_loc.iloc[diff_loc['rmse_min'].argmin(),3]  #parent is column 3
                    closest_to_min_trajs.append(int(min_ptraj))

                    diff_loc = alltrajs_rsec.loc[time,('x','y','z')] - \
                               alltrajs_rsec.loc[time,('x','y','z')].max()
                    diff_loc['parent'] = alltrajs_rsec.loc[time,'parent']
                    diff_loc['rmse_max'] = ( (diff_loc['x']**2. + diff_loc['y']**2. + diff_loc['z']**2.)/3.)**0.5
                    min_ptraj = diff_loc.iloc[diff_loc['rmse_max'].argmin(),3]  #parent is column 3
                    closest_to_max_trajs.append(int(min_ptraj)) 
 
                num_trajs_per_sc.append(alltrajs_rsec.loc[time].shape[0])

            except IndexError as error: # just one record, so don't include a median traj
                num_trajs_per_sc.append(1)

        #find the most frequent "median", "min", "max" trajectory.
        #Note these could all be the same....will leave it as is for now.
        unique, counts = np.unique(closest_to_median_trajs, return_counts=True)
        median_ptraj = unique[counts.argmax()]
        unique, counts = np.unique(closest_to_min_trajs, return_counts=True)
        min_ptraj = unique[counts.argmax()]
        unique, counts = np.unique(closest_to_max_trajs, return_counts=True)
        max_ptraj = unique[counts.argmax()]
        
        allmedian = alltrajs_rsec[alltrajs_rsec['parent']==median_ptraj]
        allmin = alltrajs_rsec[alltrajs_rsec['parent']==min_ptraj]
        allmax = alltrajs_rsec[alltrajs_rsec['parent']==max_ptraj]
        all_num_trajs_per_sc = pd.DataFrame(num_trajs_per_sc,index=times)
        len_allmedian = len(allmedian)
        len_allmin = len(allmin)
        len_allmax = len(allmax)

        cluster_num = file.split('_')[-1].split('.')[0]
        #print (cluster_num)

        #write results out to the repr_traj file
        repr_traj_file = dir_prefix+'/repr_traj/'+long_file_prefix+'_traj_modebytime_'+\
            sc_name+'.txt'
        f = open(repr_traj_file,'w')
        f.write(str(len_allmedian)+', '+str(len_allmin)+', '+str(len_allmax)+'\n')
        f.close()

        allmin.to_csv(repr_traj_file, date_format='%H:%M:%S.%f',mode='a',header=False)
        allmedian.to_csv(repr_traj_file, date_format='%H:%M:%S.%f',mode='a',header=False)
        allmax.to_csv(repr_traj_file, date_format='%H:%M:%S.%f',mode='a',header=False)

        #and to number of trajectories file
        numtraj_traj_file = dir_prefix+'/repr_traj/'+long_file_prefix+'_numtraj_bytime_'+\
            sc_name+'.txt'
        num_times = times.shape[0]
        f = open(numtraj_traj_file,'w')
        f.write(str(num_times)+'\n')
        f.close()
        all_num_trajs_per_sc.to_csv(numtraj_traj_file, date_format='%H:%M:%S.%f',
            mode='a',header=False)
