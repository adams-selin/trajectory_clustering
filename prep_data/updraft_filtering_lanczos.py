import numpy as np
import netCDF4 as nc4
import matplotlib.pyplot as plt
import glob
from scipy.signal import convolve2d


def low_pass_lanczos1d(nwts, cutoff):
    """Calculate weights for a low pass Lanczos filter.
    from https://scitools-iris.readthedocs.io/en/stable/generated/gallery/general/plot_SOI_filtering.html?highlight=lanczos%20filter
    Args:
    nwts: int
        Number of weights. must be odd, if not 1 is added.

    cutoff: float
        The cutoff frequency in inverse time steps.
    """
    #order = ((window - 1) // 2) + 1
    #nwts = 2 * order + 1

    if (nwts % 2) == 0:
        nwts = nwts + 1

    w = np.zeros([nwts])
    n = nwts // 2
    w[n] = 2 * cutoff
    k = np.arange(1.0, n)
    sigma = np.sin(np.pi * k / n) * n / (np.pi * k)
    firstfactor = np.sin(2.0 * np.pi * cutoff * k) / (np.pi * k)
    w[n - 1 : 0 : -1] = firstfactor * sigma
    w[n + 1 : -1] = firstfactor * sigma
    return w[1:-1]


def low_pass_lanczos2d(nwts, cutoff):
    """Calculate weights for a symmetrical low pass Lanczos filter in 2d. 
    Args:
    nwts: int
        Number of weights. must be odd, if not 1 is added. 

    cutoff: float
        The cutoff frequency in inverse time steps.
    """
    w1d = low_pass_lanczos1d(nwts,cutoff)
    w1d = np.append(0.,np.append(w1d,0.))
    w2d = np.zeros((w1d.shape[0],w1d.shape[0]))
    n = nwts // 2
    w2d[n,:] = w1d
    w2d[:,n] = w1d

    #i'm sure there is a better way to do this....but it works!!!!!!
    for j in range(n+1,nwts-1):
        for i in range(n+1,nwts-1):
            w2d[j,i] = w1d[j]*w1d[i]

    for j in range(n,0,-1):
        for i in range(n+1,nwts-1):
            w2d[j,i] = w1d[j]*w1d[i]

    for j in range(n+1,nwts-1):
        for i in range(n,0,-1):
            w2d[j,i] = w1d[j]*w1d[i]

    for j in range(n,0,-1):
        for i in range(n,0,-1):
            w2d[j,i] = w1d[j]*w1d[i]

    return w2d[1:-1,1:-1]



directory = '/glade/derecho/scratch/radams/20120529/clean_run/'
files = glob.glob(directory+'cm1out_00*.nc')
files.sort()
more_files = glob.glob('/glade/derecho/scratch/radams/20120529/rst000037/cm1out_00*.nc')
more_files.sort()
all_files = np.append(np.array(files[18:]),np.array(more_files))


#define our averaging kernel based on the number of points we want
# and the cutoff frequency

#let's say cutoff frequency of 1 cycle per 24 gridpoints
cutoff = 1/24.
#Based on Duchon 1979 Fig. 5, number of weights = 15 should be OK.
num_weights = 15
lanc2d = low_pass_lanczos2d(num_weights,cutoff)

#lists to store our filtered maximum updraft locations
# leftx = []
# lefty = []
rightx = []
righty = []
timestep = []


for filename in all_files:
    print (filename)
    time = filename.split('_')[-1][0:-3]
    direct = '/'.join(filename.split('/')[0:-1])

    nc = nc4.Dataset(filename,'r')
    xh = np.array(nc.variables['xh'])
    yh = np.array(nc.variables['yh'])
    zf = np.array(nc.variables['zf'])
    w = np.array(nc.variables['w'][0,:,:,:])
    seconds = np.array(nc.variables['time'][0])
    nc.close()
    mm, ss = divmod(seconds, 60)
    hh, mm = divmod(mm, 60)
    label = "%02d:%02d" % (hh,mm)

    wmax2d = w.max(axis=0)

    #and filter...
    wmax_filter = convolve2d(wmax2d,lanc2d,mode='valid') / np.sum(lanc2d)

    #mesh the the x, y for contouring
    halfn = (num_weights-2) //2
    xmesh,ymesh = np.meshgrid(xh,yh)
    xmeshf,ymeshf = np.meshgrid(xh[halfn:-halfn], yh[halfn:-halfn])

    # directory = '/glade/derecho/scratch/radams/20120529/clean_run'
    # with open(directory+'/cm1out_wmax_filter_'+time+'_incsum.npy','wb') as f:
    #     np.save(f, wmax_filter)
    #     np.save(f, xmeshf)
    #     np.save(f, ymeshf)

    """
    #This section is useful if you have splitting supercells (two updrafts)
    #identify the points of maximum unfiltered and filtered w for both left and right mover
    #use y=60km as a good dividing line (239)    
    left_peakw_filter = (wmax_filter[238:,:] == wmax_filter[238:,:].max()).nonzero()
    left_peakw_filterx = xmeshf[238:,:][left_peakw_filter][0]
    left_peakw_filtery = ymeshf[238:,:][left_peakw_filter][0]
    right_peakw_filter = (wmax_filter[:238,:] == wmax_filter[:238,:].max()).nonzero()
    right_peakw_filterx = xmeshf[:238,:][right_peakw_filter][0]
    right_peakw_filtery = ymeshf[:238,:][right_peakw_filter][0]
    #raw data
    left_peakw = (wmax2d[240:,:] == wmax2d[240:,:].max()).nonzero()
    left_peakwx = xmeshf[240:,:][left_peakw][0]
    left_peakwy = ymeshf[240:,:][left_peakw][0]
    right_peakw = (wmax2d[:240,:] == wmax2d[:240,:].max()).nonzero()
    right_peakwx = xmeshf[:240,:][right_peakw][0]
    right_peakwy = ymeshf[:240,:][right_peakw][0]

    leftx.append(left_peakw_filterx)
    lefty.append(left_peakw_filtery)
    rightx.append(right_peakw_filterx)
    righty.append(right_peakw_filtery)
    timestep.append(seconds[0])
    """

    #This section if you want to focus on just have one updraft
    #use the northern end of the embryo insertion box (y2) as a dividing line
    if (seconds >=9000) & (seconds <= 10800): #ignore extra split
        y1= ( (seconds-9000)/60 ) * 0.25 + 50 #65 #40#30  #0.37
        y2 = y1+32
    elif (seconds < 9000):
        y1= -( (seconds-5400)/60 ) * 0.0 + 65 #40#30   #0.31
        y2 = y1+32#27
    else:
        y2 = 50+32


    yval = yh[yh<y2].argmax()
    #raw data
    right_peakw = (wmax2d[:yval,:] == wmax2d[:yval,:].max()).nonzero()
    right_peakwx = xmeshf[:yval,:][right_peakw][0]
    right_peakwy = ymeshf[:yval,:][right_peakw][0]

    #filtered data
    yval = ymeshf[ymeshf[:,0]<y2,0].argmax()
    right_peakw_filter = (wmax_filter[:yval,:] == wmax_filter[:yval,:].max()).nonzero()
    right_peakw_filterx = xmeshf[:yval,:][right_peakw_filter][0]
    right_peakw_filtery = ymeshf[:yval,:][right_peakw_filter][0]


    rightx.append(right_peakw_filterx)
    righty.append(right_peakw_filtery)
    timestep.append(seconds)



#write updraft locations out to a text file
# with open('left_filterwmax.txt','w') as f:
#     for x, y, t in zip(leftx, lefty, timestep):
#         f.write(format(x,'>7.3f')+format(y,'>9.3f')+format(t,'>7.0f')+'\n')

with open(direct+'/right_filterwmax.txt','w') as f:
    for x, y, t in zip(rightx, righty, timestep):
        f.write(format(x,'>7.3f')+format(y,'>9.3f')+format(t,'>7.0f')+'\n')


