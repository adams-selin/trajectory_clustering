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
# !   2) Cycle through all times, starting from most negative,  and calculate the
# !       distribution of all cluster trajectory characteristics during that same
# !       time period (w/in time_int s either side.)
# !
# ! INPUT: File, directory prefix for merged cluster files,  total number of cluster files,
# !   epsilon, and MinLns
# ! OUTPUT: Multiple files, each containing a representative trajectory for that cluster.
# !
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


# dir_prefix='/Users/rselin/Documents/NSF_hail/streamlined_hail_scripts/left-mover/4d'
# eps='800'
# minlns='100'
# short_file_prefix = 'cm1out_haildata_traj_leftWrelative'

# dir_prefix='/Users/rselin/Documents/NSF_hail/streamlined_hail_scripts/right-mover/4d'
# eps='500'
# minlns='3'
# short_file_prefix = 'cm1out_haildata_traj_rightW'

dir_prefix='/Users/rselin/Documents/NSF_hail/streamlined_hail_scripts/right-mover/lauren'
eps='300'
minlns='20'
short_file_prefix = 'hailtraj_Wrelative_fgt45mm'
label = 'Hailstones >= 45 mm'
subtitle = 'Kingfisher'


long_file_prefix = short_file_prefix+'_'+eps+'_'+minlns
endtime = datetime.datetime(2012,5,30,0,0,0) + datetime.timedelta(seconds=3600) #dummy end time

#endtime = datetime.datetime(2012,5,29,int(time)%24,0,0) + datetime.timedelta(seconds=3600)

names = ('x1','y1','z1','x2','y2','z2','parent')
char_names = ('sec1','d1','dense1','ts1','fw1','vt1','ri1','rw1','tc1','w1',
              'sec2','d2','dense2','ts2','fw2','vt2','ri2','rw2','tc2','w2')
first = ('x1','y1','z1','sec1','d1','dense1','ts1','fw1','vt1','ri1','rw1','tc1','w1')
second = ('x2','y2','z2','sec2','d2','dense2','ts2','fw2','vt2','ri2','rw2','tc2','w2')
third =  ('x','y','z','sec','d','dense','ts','fw','vt','ri','rw','tc','w')

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
                               'G': ('37','39','55'), 'H': ('13','24','27'),
                               'J': ('43',)}
elif dir_prefix.find('lauren') >= 0:
    supercluster_choices = {'A': ('9','19','27'),
                            'B': ('38','39','5','71','86'),
                            'C': ('80','81','82','83','79','84','87'),
                            'D': ('42','43','44','66','50','72','69','68','56','55','61','26','88'),
                            'E': ('51','52','59','67','70','74','75'),
                            'F': ('29','85',),
                            'G': ('48','49','54','64','89')}
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
        allmedian = []
        allmin = []
        allmax = []
        all25p = []
        all75p = []
        for time in times:
            try:
                dum = alltrajs_rsec.loc[time].shape[1] #check if only one record

                if alltrajs_rsec.loc[time].shape[0] >= 3:
                    allmedian.append(alltrajs_rsec.loc[time].median())
                    allmin.append(alltrajs_rsec.loc[time].min())
                    allmax.append(alltrajs_rsec.loc[time].max())
                    all25p.append(alltrajs_rsec.loc[time].quantile(0.25,interpolation='linear'))
                    all75p.append(alltrajs_rsec.loc[time].quantile(0.75,interpolation='linear'))
                else:
                    allmedian.append(np.ones((num_columns))*np.nan)
                    allmin.append(np.ones((num_columns))*np.nan)
                    allmax.append(np.ones((num_columns))*np.nan)
                    all25p.append(np.ones((num_columns))*np.nan)
                    all75p.append(np.ones((num_columns))*np.nan)

                num_trajs_per_sc.append(alltrajs_rsec.loc[time].shape[0])

            except IndexError as error: # just one record, so set time as missing
                allmedian.append(np.ones((num_columns))*np.nan)
                allmin.append(np.ones((num_columns))*np.nan)
                allmax.append(np.ones((num_columns))*np.nan)
                all25p.append(np.ones((num_columns))*np.nan)
                all75p.append(np.ones((num_columns))*np.nan)
                num_trajs_per_sc.append(1)

        allmedian = pd.DataFrame(allmedian,index=times)
        allmin = pd.DataFrame(allmin,index=times)
        allmax = pd.DataFrame(allmax,index=times)
        all25p = pd.DataFrame(all25p,index=times)
        all75p = pd.DataFrame(all75p,index=times)
        all_num_trajs_per_sc = pd.DataFrame(num_trajs_per_sc,index=times)

        cluster_num = file.split('_')[-1].split('.')[0]
        #print (cluster_num)

        #write results out to the repr_traj file
        repr_traj_file = dir_prefix+'/repr_traj/'+long_file_prefix+'_traj_bytime_'+\
            sc_name+'.txt'
        num_times = times.shape[0]
        f = open(repr_traj_file,'w')
        f.write(str(num_times)+'\n')
        f.close()

        allmin.to_csv(repr_traj_file, date_format='%H:%M:%S.%f',mode='a',header=False)
        all25p.to_csv(repr_traj_file, date_format='%H:%M:%S.%f',mode='a',header=False)
        allmedian.to_csv(repr_traj_file, date_format='%H:%M:%S.%f',mode='a',header=False)
        all75p.to_csv(repr_traj_file, date_format='%H:%M:%S.%f',mode='a',header=False)
        allmax.to_csv(repr_traj_file, date_format='%H:%M:%S.%f',mode='a',header=False)

        #and to number of trajectories file
        numtraj_traj_file = dir_prefix+'/repr_traj/'+long_file_prefix+'_numtraj_bytime_'+\
            sc_name+'.txt'
        f = open(numtraj_traj_file,'w')
        f.write(str(num_times)+'\n')
        f.close()
        all_num_trajs_per_sc.to_csv(numtraj_traj_file, date_format='%H:%M:%S.%f',
            mode='a',header=False)
