import numpy as np
import glob as glob

directory = '/glade/scratch/radams/hail/20120529/'
cases = ['right']

for thiscase in cases:

    case_files = glob.glob(directory+'rst*/subtraj/*'+thiscase+'Wrelative_subtraj.txt')
    case_files.sort()
    case_char_files = glob.glob(directory+'rst*/subtraj/*'+thiscase+'Wrelative_subtraj_chars.txt')
    case_char_files.sort()

    for case_file, case_char_file in zip(case_files, case_char_files):
        print (case_file)
        rstnum = int(case_file.split('rst')[1][0:3])

        casedata = np.genfromtxt(case_file)
        casechardata = np.genfromtxt(case_char_file)
        rstnum_array = np.ones((casedata.shape[0])) * rstnum * 1E3
        #casedata = np.concatenate((casedata,rstnum_array),axis=1)
        casedata[:,6] = casedata[:,6] + rstnum_array

        if rstnum == 1:
            allcasedata = casedata
            allcasechardata = casechardata
        else:
            allcasedata = np.concatenate((allcasedata,casedata),axis=0)
            allcasechardata = np.concatenate((allcasechardata,casechardata),axis=0)


    #write combined data out to a new total subtraj file
    outfile = directory+'subtraj/'+case_file.split('/')[-1]
    np.savetxt(outfile,allcasedata,fmt='%10.1f %10.1f %10.1f %10.1f %10.1f %10.1f  %7.0f')

    #and new total subtraj chars file
    chars_outfile = directory+'subtraj/'+case_char_file.split('/')[-1]
    chars_fmt = '%12.4f %12.4f %12.4f %12.4f %12.4f '+\
                '%12.4f %12.4f %12.4f %12.4f %12.4f '+\
                '%12.4f %12.4f %12.4f %12.4f %12.4f '+\
                '%12.4f %12.4f %12.4f %12.4f %12.4f '
    np.savetxt(chars_outfile,allcasechardata,fmt=chars_fmt)