from ReadStats import Statistics, Pdfs 
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.style
from matplotlib import rc
from scipy import interpolate


rc('text', usetex=True)
rc('text.latex', preamble=r"\usepackage{fourier}")
rc('font', family='serif')
rc('font', size=20)
rc('axes', linewidth=1.5)
rc('lines', linewidth=2)

opath = '/scratch/local1/m300551/ForKatherine/plots/3D/Re042/'

path_1 = '/scratch/local1/m300551/ForKatherine/qCBL_3D/Re042/2560x576x2560/'
path_2 = '/scratch/local1/m300551/ForKatherine/qCBL_3D/Re042/2560x704x2560/'
path_3 = '/scratch/local1/m300551/ForKatherine/qCBL_3D/Re042/2560x896x2560/'

path_S20 = '/scratch/local1/m300551/ForKatherine/qCBL_3D/Re042/3072x960x4608-S20/'
##########################################################################
# Constants

nu = 1./25000.
B0 = 0.005
N = np.sqrt(3)
L_0 = (B0/N**3)**0.5
ceps = 0.1
cb=0.1

#######################################################################
# Calculate a running mean of a time series with a specified window size.
# Window size is number of entries on either side of the entry being averaged.
# Leaves out entries at the beginning and end of the time series such that the
# window size is always the same, but the resulting time series is shorter than
# the original.

def runningmean(timeseries,window):
    nt = len(timeseries)
    outseries = np.zeros(nt-(window*2))
    for n in range(window,nt-window):
        outseries[n-window] = np.mean(timeseries[n-window:n+window+1])
    return outseries

##########################################################################
# Stats

NS42_1 = Statistics(path_1+'stats/pdftimes/avg20500-53000.nc')
NS42_2 = Statistics(path_2+'stats/pdftimes/avg60000-74500.nc')
NS42_3 = Statistics(path_3+'stats/pdftimes/avg83000-127500.nc')

S20 = Statistics(path_S20+'stats/pdftimes/avg42000-148000.nc')

z_enc_NS42 = np.concatenate((NS42_1.z_enc,NS42_2.z_enc,NS42_3.z_enc))


############################################################################

# Pdf

vortlist_1 = [path_1+'stats/pdfs/pdf20500.LnEnstrophyW_iW_i',path_1+'stats/pdfs/pdf24000.LnEnstrophyW_iW_i',path_1+'stats/pdfs/pdf27500.LnEnstrophyW_iW_i',path_1+'stats/pdfs/pdf32000.LnEnstrophyW_iW_i',path_1+'stats/pdfs/pdf36500.LnEnstrophyW_iW_i',path_1+'stats/pdfs/pdf41500.LnEnstrophyW_iW_i',path_1+'stats/pdfs/pdf47500.LnEnstrophyW_iW_i',path_1+'stats/pdfs/pdf53000.LnEnstrophyW_iW_i']
vortlist_2 = [path_2+'stats/pdfs/pdf60000.LnEnstrophyW_iW_i',path_2+'stats/pdfs/pdf67000.LnEnstrophyW_iW_i',path_2+'stats/pdfs/pdf74500.LnEnstrophyW_iW_i']
vortlist_3 = [path_3+'stats/pdfs/pdf83000.LnEnstrophyW_iW_i',path_3+'stats/pdfs/pdf93000.LnEnstrophyW_iW_i',path_3+'stats/pdfs/pdf102500.LnEnstrophyW_iW_i',path_3+'stats/pdfs/pdf113000.LnEnstrophyW_iW_i',path_3+'stats/pdfs/pdf127500.LnEnstrophyW_iW_i']

pvlist_1 = [path_1+'stats/pdfs/pdf20500.LnPotentialEnstrophy',path_1+'stats/pdfs/pdf24000.LnPotentialEnstrophy',path_1+'stats/pdfs/pdf27500.LnPotentialEnstrophy',path_1+'stats/pdfs/pdf32000.LnPotentialEnstrophy',path_1+'stats/pdfs/pdf36500.LnPotentialEnstrophy',path_1+'stats/pdfs/pdf41500.LnPotentialEnstrophy',path_1+'stats/pdfs/pdf47500.LnPotentialEnstrophy',path_1+'stats/pdfs/pdf53000.LnPotentialEnstrophy']
pvlist_2 = [path_2+'stats/pdfs/pdf60000.LnPotentialEnstrophy',path_2+'stats/pdfs/pdf67000.LnPotentialEnstrophy',path_2+'stats/pdfs/pdf74500.LnPotentialEnstrophy']
pvlist_3 = [path_3+'stats/pdfs/pdf83000.LnPotentialEnstrophy',path_3+'stats/pdfs/pdf93000.LnPotentialEnstrophy',path_3+'stats/pdfs/pdf102500.LnPotentialEnstrophy',path_3+'stats/pdfs/pdf113000.LnPotentialEnstrophy',path_3+'stats/pdfs/pdf127500.LnPotentialEnstrophy']

vortlist_S20 = [path_S20+'stats/pdfs/pdf42000.LnEnstrophyW_iW_i',path_S20+'stats/pdfs/pdf45000.LnEnstrophyW_iW_i',path_S20+'stats/pdfs/pdf51000.LnEnstrophyW_iW_i',path_S20+'stats/pdfs/pdf58000.LnEnstrophyW_iW_i',path_S20+'stats/pdfs/pdf66000.LnEnstrophyW_iW_i',path_S20+'stats/pdfs/pdf75000.LnEnstrophyW_iW_i',path_S20+'stats/pdfs/pdf84000.LnEnstrophyW_iW_i',path_S20+'stats/pdfs/pdf94000.LnEnstrophyW_iW_i',path_S20+'stats/pdfs/pdf105000.LnEnstrophyW_iW_i',path_S20+'stats/pdfs/pdf117000.LnEnstrophyW_iW_i',path_S20+'stats/pdfs/pdf130000.LnEnstrophyW_iW_i',path_S20+'stats/pdfs/pdf148000.LnEnstrophyW_iW_i']

pvlist_S20 = [path_S20+'stats/pdfs/pdf42000.LnPotentialEnstrophy',path_S20+'stats/pdfs/pdf45000.LnPotentialEnstrophy',path_S20+'stats/pdfs/pdf51000.LnPotentialEnstrophy',path_S20+'stats/pdfs/pdf58000.LnPotentialEnstrophy',path_S20+'stats/pdfs/pdf66000.LnPotentialEnstrophy',path_S20+'stats/pdfs/pdf75000.LnPotentialEnstrophy',path_S20+'stats/pdfs/pdf84000.LnPotentialEnstrophy',path_S20+'stats/pdfs/pdf94000.LnPotentialEnstrophy',path_S20+'stats/pdfs/pdf105000.LnPotentialEnstrophy',path_S20+'stats/pdfs/pdf117000.LnPotentialEnstrophy',path_S20+'stats/pdfs/pdf130000.LnPotentialEnstrophy',path_S20+'stats/pdfs/pdf148000.LnPotentialEnstrophy']

NS42_vortpdf_1 = Pdfs(vortlist_1,path_1+'y.dat')
NS42_vortpdf_2 = Pdfs(vortlist_2,path_2+'y.dat')
NS42_vortpdf_3 = Pdfs(vortlist_3,path_3+'y.dat')

NS42_pvpdf_1 = Pdfs(pvlist_1,path_1+'y.dat')
NS42_pvpdf_2 = Pdfs(pvlist_2,path_2+'y.dat')
NS42_pvpdf_3 = Pdfs(pvlist_3,path_3+'y.dat')

vortpdf_S20 = Pdfs(vortlist_S20,path_S20+'y.dat')
pvpdf_S20 = Pdfs(pvlist_S20,path_S20+'y.dat')

# Create grid on which to interpolate pdfs

NS42_vortpdf_interp_data = Pdfs([path_3+'stats/pdfs/pdf102500.LnEnstrophyW_iW_i'],path_3+'y.dat')
NS42_pvpdf_interp_data = Pdfs([path_3+'stats/pdfs/pdf102500.LnPotentialEnstrophy'],path_3+'y.dat')

S20_vortpdf_interp_data = Pdfs([path_S20+'stats/pdfs/pdf75000.LnEnstrophyW_iW_i'],path_S20+'y.dat')
S20_pvpdf_interp_data = Pdfs([path_S20+'stats/pdfs/pdf75000.LnPotentialEnstrophy'],path_S20+'y.dat')

# Interpolate pdfs in y-direction

NS42_vortpdf_interp_1_y = np.zeros((len(vortlist_1),NS42_3.y_len,NS42_vortpdf_1.nb+2))
for n in range(len(vortlist_1)):
    for i in range(NS42_vortpdf_1.nb+2):
        NS42_vortpdf_interp_1_y[n,:,i] = np.interp(NS42_3.y,NS42_1.y,NS42_vortpdf_1.pdf[n,:-1,i])

NS42_vortpdf_interp_2_y = np.zeros((len(vortlist_2),NS42_3.y_len,NS42_vortpdf_2.nb+2))
for n in range(len(vortlist_2)):
    for i in range(NS42_vortpdf_2.nb+2):
        NS42_vortpdf_interp_2_y[n,:,i] = np.interp(NS42_3.y,NS42_2.y,NS42_vortpdf_2.pdf[n,:-1,i])

NS42_pvpdf_interp_1_y = np.zeros((len(pvlist_1),NS42_3.y_len,NS42_pvpdf_1.nb+2))
for n in range(len(pvlist_1)):
    for i in range(NS42_pvpdf_1.nb+2):
        NS42_pvpdf_interp_1_y[n,:,i] = np.interp(NS42_3.y,NS42_1.y,NS42_pvpdf_1.pdf[n,:-1,i])

NS42_pvpdf_interp_2_y = np.zeros((len(pvlist_2),NS42_3.y_len,NS42_pvpdf_2.nb+2))
for n in range(len(pvlist_2)):
    for i in range(NS42_pvpdf_2.nb+2):
        NS42_pvpdf_interp_2_y[n,:,i] = np.interp(NS42_3.y,NS42_2.y,NS42_pvpdf_2.pdf[n,:-1,i])
        
# Interpolate pdfs in x-direction

NS42_vortpdf_interp_1 = np.zeros((len(vortlist_1),NS42_3.y_len,NS42_vortpdf_1.nb))
for n in range(len(vortlist_1)):
    for j in range(NS42_3.y_len):
        NS42_vortpdf_interp_1[n,j,:] = np.interp(NS42_vortpdf_interp_data.xy[0,0,j,:],np.linspace(NS42_vortpdf_interp_1_y[n,j,NS42_vortpdf_1.nb],NS42_vortpdf_interp_1_y[n,j,NS42_vortpdf_1.nb + 1],num=NS42_vortpdf_1.nb),NS42_vortpdf_interp_1_y[n,j,:-2])

NS42_vortpdf_interp_2 = np.zeros((len(vortlist_2),NS42_3.y_len,NS42_vortpdf_2.nb))
for n in range(len(vortlist_2)):
    for j in range(NS42_3.y_len):
        NS42_vortpdf_interp_2[n,j,:] = np.interp(NS42_vortpdf_interp_data.xy[0,0,j,:],np.linspace(NS42_vortpdf_interp_2_y[n,j,NS42_vortpdf_2.nb],NS42_vortpdf_interp_2_y[n,j,NS42_vortpdf_2.nb + 1],num=NS42_vortpdf_2.nb),NS42_vortpdf_interp_2_y[n,j,:-2])

NS42_vortpdf_interp_3 = np.zeros((len(vortlist_3),NS42_3.y_len,NS42_vortpdf_3.nb))
for n in range(len(vortlist_3)):
    for j in range(NS42_3.y_len):
        NS42_vortpdf_interp_3[n,j,:] = np.interp(NS42_vortpdf_interp_data.xy[0,0,j,:],NS42_vortpdf_3.xy[0,n,j,:],NS42_vortpdf_3.pdf[n,j,:-2])

NS42_vortpdf_interp = np.concatenate((NS42_vortpdf_interp_1,NS42_vortpdf_interp_2,NS42_vortpdf_interp_3),axis=0)

NS42_pvpdf_interp_1 = np.zeros((len(pvlist_1),NS42_3.y_len,NS42_pvpdf_1.nb))
for n in range(len(pvlist_1)):
    for j in range(NS42_3.y_len):
        NS42_pvpdf_interp_1[n,j,:] = np.interp(NS42_pvpdf_interp_data.xy[0,0,j,:],np.linspace(NS42_pvpdf_interp_1_y[n,j,NS42_pvpdf_1.nb],NS42_pvpdf_interp_1_y[n,j,NS42_pvpdf_1.nb + 1],num=NS42_pvpdf_1.nb),NS42_pvpdf_interp_1_y[n,j,:-2])

NS42_pvpdf_interp_2 = np.zeros((len(pvlist_2),NS42_3.y_len,NS42_pvpdf_2.nb))
for n in range(len(pvlist_2)):
    for j in range(NS42_3.y_len):
        NS42_pvpdf_interp_2[n,j,:] = np.interp(NS42_pvpdf_interp_data.xy[0,0,j,:],np.linspace(NS42_pvpdf_interp_2_y[n,j,NS42_pvpdf_2.nb],NS42_pvpdf_interp_2_y[n,j,NS42_pvpdf_2.nb + 1],num=NS42_pvpdf_2.nb),NS42_pvpdf_interp_2_y[n,j,:-2])

NS42_pvpdf_interp_3 = np.zeros((len(pvlist_3),NS42_3.y_len,NS42_pvpdf_3.nb))
for n in range(len(pvlist_3)):
    for j in range(NS42_3.y_len):
        NS42_pvpdf_interp_3[n,j,:] = np.interp(NS42_pvpdf_interp_data.xy[0,0,j,:],NS42_pvpdf_3.xy[0,n,j,:],NS42_pvpdf_3.pdf[n,j,:-2])

NS42_pvpdf_interp = np.concatenate((NS42_pvpdf_interp_1,NS42_pvpdf_interp_2,NS42_pvpdf_interp_3),axis=0)

S20_vortpdf_interp = np.zeros((len(vortlist_S20),S20.y_len,vortpdf_S20.nb))
for n in range(len(vortlist_S20)):
    for j in range(S20.y_len):
        S20_vortpdf_interp[n,j,:] = np.interp(S20_vortpdf_interp_data.xy[0,0,j,:],vortpdf_S20.xy[0,n,j,:],vortpdf_S20.pdf[n,j,:-2])

S20_pvpdf_interp = np.zeros((len(pvlist_S20),S20.y_len,pvpdf_S20.nb))
for n in range(len(pvlist_S20)):
    for j in range(S20.y_len):
        S20_pvpdf_interp[n,j,:] = np.interp(S20_pvpdf_interp_data.xy[0,0,j,:],pvpdf_S20.xy[0,n,j,:],pvpdf_S20.pdf[n,j,:-2])

# Running mean of pdfs

NS42_vortpdf_interp_runmean = np.zeros((np.ma.size(NS42_vortpdf_interp,0)-2,NS42_3.y_len,NS42_vortpdf_1.nb))
for n in range(1,np.ma.size(NS42_vortpdf_interp,0)-1):
    NS42_vortpdf_interp_runmean[n-1,:,:] = np.mean(NS42_vortpdf_interp[n-1:n+2,:,:],axis=0)

S20_vortpdf_interp_runmean = np.zeros((np.ma.size(S20_vortpdf_interp,0)-2,S20.y_len,vortpdf_S20.nb))
for n in range(1,np.ma.size(S20_vortpdf_interp,0)-1):
    S20_vortpdf_interp_runmean[n-1,:,:] = np.mean(S20_vortpdf_interp[n-1:n+2,:,:],axis=0)

NS42_pvpdf_interp_runmean = np.zeros((np.ma.size(NS42_pvpdf_interp,0)-2,NS42_3.y_len,NS42_pvpdf_1.nb))
for n in range(1,np.ma.size(NS42_pvpdf_interp,0)-1):
    NS42_pvpdf_interp_runmean[n-1,:,:] = np.mean(NS42_pvpdf_interp[n-1:n+2,:,:],axis=0)

S20_pvpdf_interp_runmean = np.zeros((np.ma.size(S20_pvpdf_interp,0)-2,S20.y_len,pvpdf_S20.nb))
for n in range(1,np.ma.size(S20_pvpdf_interp,0)-1):
    S20_pvpdf_interp_runmean[n-1,:,:] = np.mean(S20_pvpdf_interp[n-1:n+2,:,:],axis=0)

# Find where pdf has a maximum at each height

maxvort_NS42 = np.zeros((np.ma.size(NS42_vortpdf_interp_runmean,0),NS42_3.y_len))
maxprob_vort_NS42 = np.zeros((np.ma.size(NS42_vortpdf_interp_runmean,0),NS42_3.y_len))
for t in range(0,np.ma.size(NS42_vortpdf_interp_runmean,0)):
    for j in range(0,NS42_3.y_len):
        maxvort_NS42[t,j] = NS42_vortpdf_interp_data.xy[0,0,j,np.argmax(NS42_vortpdf_interp_runmean[t,j,:])]
        maxprob_vort_NS42[t,j] = np.max(NS42_vortpdf_interp_runmean[t,j,:])
maxvort_NS42 = np.log10(np.exp(maxvort_NS42)/(ceps*B0/nu))

maxvort_S20 = np.zeros((np.ma.size(S20_vortpdf_interp_runmean,0),S20.y_len))
maxprob_vort_S20 = np.zeros((np.ma.size(S20_vortpdf_interp_runmean,0),S20.y_len))
for t in range(0,np.ma.size(S20_vortpdf_interp_runmean,0)):
    for j in range(0,S20.y_len):
        maxvort_S20[t,j] = S20_vortpdf_interp_data.xy[0,0,j,np.argmax(S20_vortpdf_interp_runmean[t,j,:])]
        maxprob_vort_S20[t,j] = np.max(S20_vortpdf_interp_runmean[t,j,:])
maxvort_S20 = np.log10(np.exp(maxvort_S20)/(ceps*B0/nu))

maxpv_NS42 = np.zeros((np.ma.size(NS42_pvpdf_interp_runmean,0),NS42_3.y_len))
maxprob_pv_NS42 = np.zeros((np.ma.size(NS42_pvpdf_interp_runmean,0),NS42_3.y_len))
for t in range(0,np.ma.size(NS42_pvpdf_interp_runmean,0)):
    for j in range(0,NS42_3.y_len):
        maxpv_NS42[t,j] = NS42_pvpdf_interp_data.xy[0,0,j,np.argmax(NS42_pvpdf_interp_runmean[t,j,:])]
        maxprob_pv_NS42[t,j] = np.max(NS42_pvpdf_interp_runmean[t,j,:])
    maxpv_NS42[t,:] = np.log10(np.exp(maxpv_NS42[t,:])/(cb*ceps*(z_enc_NS42[t+1]/L_0)**(-4./3.)*(B0/nu)**3/42))

maxpv_S20 = np.zeros((np.ma.size(S20_pvpdf_interp_runmean,0),S20.y_len))
maxprob_pv_S20 = np.zeros((np.ma.size(S20_pvpdf_interp_runmean,0),S20.y_len))
for t in range(0,np.ma.size(S20_pvpdf_interp_runmean,0)):
    for j in range(0,S20.y_len):
        maxpv_S20[t,j] = S20_pvpdf_interp_data.xy[0,0,j,np.argmax(S20_pvpdf_interp_runmean[t,j,:])]
        maxprob_pv_S20[t,j] = np.max(S20_pvpdf_interp_runmean[t,j,:])
    maxpv_S20[t,:] = np.log10(np.exp(maxpv_S20[t,:])/(cb*ceps*(S20.z_enc[t+1]/L_0)**(-4./3.)*(B0/nu)**3/42))

# Find jump in maxvort/maxpv

maxit_vort_NS42 = np.zeros(np.ma.size(NS42_vortpdf_interp_runmean,0))
y_maxit_vort_NS42 = np.zeros(np.ma.size(NS42_vortpdf_interp_runmean,0))
maxvort_it_NS42 = np.zeros(np.ma.size(NS42_vortpdf_interp_runmean,0))
for t in range(0,np.ma.size(NS42_vortpdf_interp_runmean,0)):
    for j in range(0,NS42_3.y_len):
        if np.abs(maxvort_NS42[t,j+1])-np.abs(maxvort_NS42[t,j]) > 0.2:
            maxit_vort_NS42[t] = j+1
            break
    y_maxit_vort_NS42[t] = NS42_3.y[int(maxit_vort_NS42[t])]/z_enc_NS42[t+1]
    maxvort_it_NS42[t] = (maxvort_NS42[t,int(maxit_vort_NS42[t]-1)]+maxvort_NS42[t,int(maxit_vort_NS42[t])])/2

maxit_pv_NS42 = np.zeros(np.ma.size(NS42_pvpdf_interp_runmean,0))
y_maxit_pv_NS42 = np.zeros(np.ma.size(NS42_pvpdf_interp_runmean,0))
maxpv_it_NS42 = np.zeros(np.ma.size(NS42_pvpdf_interp_runmean,0))
for t in range(0,np.ma.size(NS42_pvpdf_interp_runmean,0)):
    for j in range(NS42_1.z_enc_arg[1],NS42_3.y_len):
        if np.abs(maxpv_NS42[t,j+1])-np.abs(maxpv_NS42[t,j]) > 0.2:
            maxit_pv_NS42[t] = j+1
            break
    y_maxit_pv_NS42[t] = NS42_3.y[int(maxit_pv_NS42[t])]/z_enc_NS42[t+1]
    maxpv_it_NS42[t] = (maxpv_NS42[t,int(maxit_pv_NS42[t]-1)]+maxpv_NS42[t,int(maxit_pv_NS42[t])])/2

maxit_vort_S20 = np.zeros(np.ma.size(S20_vortpdf_interp_runmean,0))
y_maxit_vort_S20 = np.zeros(np.ma.size(S20_vortpdf_interp_runmean,0))
maxvort_it_S20 = np.zeros(np.ma.size(S20_vortpdf_interp_runmean,0))
for t in range(0,np.ma.size(S20_vortpdf_interp_runmean,0)):
    for j in range(0,S20.y_len):
        if np.abs(maxvort_S20[t,j+1])-np.abs(maxvort_S20[t,j]) > 0.5:
            maxit_vort_S20[t] = j+1
            break
    y_maxit_vort_S20[t] = S20.y[int(maxit_vort_S20[t])]/S20.z_enc[t+1]
    maxvort_it_S20[t] = (maxvort_S20[t,int(maxit_vort_S20[t]-1)]+maxvort_S20[t,int(maxit_vort_S20[t])])/2

maxit_pv_S20 = np.zeros(np.ma.size(S20_pvpdf_interp_runmean,0))
y_maxit_pv_S20 = np.zeros(np.ma.size(S20_pvpdf_interp_runmean,0))
maxpv_it_S20 = np.zeros(np.ma.size(S20_pvpdf_interp_runmean,0))
for t in range(0,np.ma.size(S20_pvpdf_interp_runmean,0)):
    for j in range(0,S20.y_len):
        if np.abs(maxpv_S20[t,j+1])-np.abs(maxpv_S20[t,j]) > 0.7:
            maxit_pv_S20[t] = j+1
            break
    y_maxit_pv_S20[t] = S20.y[int(maxit_pv_S20[t])]/S20.z_enc[t+1]
    maxpv_it_S20[t] = (maxpv_S20[t,int(maxit_pv_S20[t]-1)]+maxpv_S20[t,int(maxit_pv_S20[t])])/2

# Find saddle as point where maxprob has a minimum

y_vort_NS42_saddle = np.zeros(np.ma.size(NS42_vortpdf_interp_runmean,0))
maxvort_NS42_saddle = np.zeros(np.ma.size(NS42_vortpdf_interp_runmean,0))
for t in range(np.ma.size(NS42_vortpdf_interp_runmean,0)):
    y_vort_NS42_saddle[t] = NS42_3.y[np.argmin(maxprob_vort_NS42[t,:-150])]
    maxvort_NS42_saddle[t] = maxvort_NS42[t,np.argmin(np.abs(y_vort_NS42_saddle[t]-NS42_3.y))]
    y_vort_NS42_saddle[t] = y_vort_NS42_saddle[t]/z_enc_NS42[t+1]

y_vort_S20_saddle = np.zeros(np.ma.size(S20_vortpdf_interp_runmean,0))
maxvort_S20_saddle = np.zeros(np.ma.size(S20_vortpdf_interp_runmean,0))
for t in range(np.ma.size(S20_vortpdf_interp_runmean,0)):
    y_vort_S20_saddle[t] = S20.y[np.argmin(maxprob_vort_S20[t,:-100])]
    maxvort_S20_saddle[t] = maxvort_S20[t,np.argmin(np.abs(y_vort_S20_saddle[t]-S20.y))]
    y_vort_S20_saddle[t] = y_vort_S20_saddle[t]/S20.z_enc[t+1]

y_pv_NS42_saddle = np.zeros(np.ma.size(NS42_pvpdf_interp_runmean,0))
maxpv_NS42_saddle = np.zeros(np.ma.size(NS42_pvpdf_interp_runmean,0))
for t in range(np.ma.size(NS42_pvpdf_interp_runmean,0)):
    y_pv_NS42_saddle[t] = NS42_3.y[np.argmin(maxprob_pv_NS42[t,np.argmin(np.abs(z_enc_NS42[t+1]-NS42_3.y)):])+np.argmin(np.abs(z_enc_NS42[t+1]-NS42_3.y))]
    maxpv_NS42_saddle[t] = maxpv_NS42[t,np.argmin(np.abs(y_pv_NS42_saddle[t]-NS42_1.y))]
    y_pv_NS42_saddle[t] = y_pv_NS42_saddle[t]/z_enc_NS42[t+1]


y_pv_S20_saddle = np.zeros(np.ma.size(S20_pvpdf_interp_runmean,0))
maxpv_S20_saddle = np.zeros(np.ma.size(S20_pvpdf_interp_runmean,0))
for t in range(np.ma.size(S20_pvpdf_interp_runmean,0)):
    y_pv_S20_saddle[t] = S20.y[np.argmin(maxprob_pv_S20[t,S20.z_enc_arg[t+1]:])+S20.z_enc_arg[t+1]]
    maxpv_S20_saddle[t] = maxpv_S20[t,np.argmin(np.abs(y_pv_S20_saddle[t]-S20.y))]
    y_pv_S20_saddle[t] = y_pv_S20_saddle[t]/S20.z_enc[t+1]

# concatenate over time

time = np.concatenate((NS42_1.z_enc/L_0,NS42_2.z_enc/L_0,NS42_3.z_enc/L_0))
z_is = np.concatenate((NS42_1.z_is/NS42_1.z_enc,NS42_2.z_is/NS42_2.z_enc,NS42_3.z_is/NS42_3.z_enc))
z_if = np.concatenate((NS42_1.z_if/NS42_1.z_enc,NS42_2.z_if/NS42_2.z_enc,NS42_3.z_if/NS42_3.z_enc))


#####################################################################
# Plot

blues = matplotlib.cm.get_cmap('Blues')
oranges = matplotlib.cm.get_cmap('Oranges')
greys = matplotlib.cm.get_cmap('Greys')

def multicolor_ylabel(ax,list_of_strings,list_of_colors,axis='x',anchorpad=0,**kw):
    """this function creates axes labels with multiple colors
    ax specifies the axes object where the labels should be drawn
    list_of_strings is a list of all of the text items
    list_if_colors is a corresponding list of colors for the strings
    axis='x', 'y', or 'both' and specifies which label(s) should be drawn"""
    from matplotlib.offsetbox import AnchoredOffsetbox, TextArea, HPacker, VPacker

    # x-axis label
    if axis=='x' or axis=='both':
        boxes = [TextArea(text, textprops=dict(color=color, ha='left',va='bottom',**kw)) 
                    for text,color in zip(list_of_strings,list_of_colors) ]
        xbox = HPacker(children=boxes,align="center",pad=0, sep=5)
        anchored_xbox = AnchoredOffsetbox(loc=3, child=xbox, pad=anchorpad,frameon=False,bbox_to_anchor=(0.2, -0.09),
                                          bbox_transform=ax.transAxes, borderpad=0.)
        ax.add_artist(anchored_xbox)

    # y-axis label
    if axis=='y' or axis=='both':
        boxes = [TextArea(text, textprops=dict(color=color, ha='left',va='bottom',rotation=90,**kw)) 
                     for text,color in zip(list_of_strings[::-1],list_of_colors) ]
        ybox = VPacker(children=boxes,align="center", pad=0, sep=5)
        anchored_ybox = AnchoredOffsetbox(loc=3, child=ybox, pad=anchorpad, frameon=False, bbox_to_anchor=(-0.22, 0.05), 
                                          bbox_transform=ax.transAxes, borderpad=0.)
        ax.add_artist(anchored_ybox)


f, (ax1,ax2) = plt.subplots(1,2,figsize=(10,5))
ax1.grid(True)
ax2.grid(True)
ax1.set_ylim(1,1.55)
ax1.set_xlim(15,30)
ax2.set_ylim(-4,0)
ax2.set_xlim(15,30)
ax1.plot(time[1:-1],y_maxit_vort_NS42,c=blues(0.5))
ax1.plot(time[1:-1],y_maxit_pv_NS42,c=oranges(0.5))
ax1.plot(S20.z_enc[1:-1]/L_0,y_maxit_vort_S20,c=blues(0.9))
ax1.plot(S20.z_enc[1:-1]/L_0,y_maxit_pv_S20,c=oranges(0.9))
ax1.plot(S20.z_enc[1:-1]/L_0,runningmean(S20.z_ig/S20.z_enc,1),c=greys(0.9),ls='--',label=r'$z_{i,g}/z_\mathrm{enc}$')
ax1.plot(time[1:-1],runningmean(z_is,1),c=greys(.5),ls='-.',label=r'$z_{i,s}/z_\mathrm{enc}$')
ax1.plot(time[1:-1],runningmean(z_if,1),c=greys(.5),ls=':',label=r'$z_{i,f}/z_\mathrm{enc}$')
ax2.plot(time[1:-1],maxvort_it_NS42,c=blues(0.5))
ax2.plot(time[1:-1],maxpv_it_NS42,c=oranges(0.5))
ax2.plot(S20.z_enc[1:-1]/L_0,maxvort_it_S20,c=blues(0.9))
ax2.plot(S20.z_enc[1:-1]/L_0,maxpv_it_S20,c=oranges(0.9))
ax2.text(16,-3.4,r'$Fr_0=0$',color=greys(0.5),fontsize=20)
ax2.text(16,-3.9,r'$Fr_0=20$',color=greys(0.9),fontsize=20)
ax1.set_xlabel(r'$z_\mathrm{enc}/L_0$')
ax2.set_xlabel(r'$z_\mathrm{enc}/L_0$')
ax1.set_ylabel(r'$z_\mathrm{saddle}/z_\mathrm{enc}$')
multicolor_ylabel(ax2,('$\mathrm{log}_{10}(\omega^2/\omega_0^2),$','$\mathrm{log}_{10}(\Pi^2/\Pi_0^2)$'),(oranges(0.7),blues(0.7)),axis='y',fontsize=20)
ax1.set_title('(a)',fontsize=20,loc='left')
ax2.set_title('(b)',fontsize=20,loc='left')
ax1.legend(loc='best',fontsize=20,borderaxespad=0.1,handlelength=1.2,ncol=2,columnspacing=1.2)
plt.tight_layout()
plt.savefig(opath+'pdfs_saddle_height_time_S20_S0_interpxy.pdf')
plt.show()

f, (ax1,ax2) = plt.subplots(1,2,figsize=(10,5))
ax1.grid(True,linewidth=1.5)
ax2.grid(True,linewidth=1.5)
ax1.tick_params(bottom=False,left=False)
ax2.tick_params(bottom=False,left=False)
ax1.set_ylim(1,1.55)
ax1.set_xlim(15,30)
ax2.set_ylim(-2,0)
ax2.set_xlim(15,30)
ax1.plot(time[1:-1],y_maxit_vort_NS42,c=blues(0.5))
ax1.plot(S20.z_enc[1:-1]/L_0,y_maxit_vort_S20,c=blues(0.9))
ax1.plot(S20.z_enc[1:-1]/L_0,runningmean(S20.z_ig/S20.z_enc,1),c=greys(0.9),ls='--',label=r'$z_{i,g}/z_\mathrm{enc}$')
ax1.plot(time[1:-1],runningmean(z_if,1),c=greys(.5),ls='--',label=r'$z_{i,f}/z_\mathrm{enc}$')
ax2.plot(time[1:-1],maxvort_it_NS42,c=blues(0.5),label=r'$Fr_0=0$')
ax2.plot(S20.z_enc[1:-1]/L_0,maxvort_it_S20,c=blues(0.9),label=r'$Fr_0=20$')
# ax2.text(16,-1.7,r'$Fr_0=0$',color=blues(0.5),fontsize=20)
# ax2.text(16,-1.9,r'$Fr_0=20$',color=blues(0.9),fontsize=20)
ax1.set_xlabel(r'$z_\mathrm{enc}/L_0$')
ax2.set_xlabel(r'$z_\mathrm{enc}/L_0$')
ax1.set_ylabel(r'$z_\mathrm{saddle}/z_\mathrm{enc}$')
ax2.set_ylabel(r'$\mathrm{log}_{10}(\omega^2/\omega_0^2)$')
ax1.set_title('(a)',fontsize=20,loc='left')
ax2.set_title('(b)',fontsize=20,loc='left')
ax1.legend(loc='best',fontsize=20)
ax2.legend(loc='best',fontsize=20)
plt.tight_layout()
plt.savefig(opath+'pdfs_saddle_vort_height_time_S20_S0_interpxy.pdf',bbox_inches='tight')
plt.show()

#Presentations

# f, (ax1,ax2) = plt.subplots(1,2,sharex='all',figsize=(10,5))
# ax1.grid(True)
# ax2.grid(True)
# ax1.set_ylim(1,1.4)
# ax2.set_ylim(-4.5,0.5)
# ax1.plot(time,y_vort_saddle)
# ax1.plot(time,y_pv_saddle)
# ax1.plot(time,z_is,'k--',label=r'$z_{i,s}/z_\mathrm{enc}$')
# ax1.plot(time,z_if,'k-.',label=r'$z_{i,f}/z_\mathrm{enc}$')
# ax2.plot(time,maxvort_saddle,label=r'$\log_{10}(\omega^2/\omega_0^2)$')
# ax2.plot(time,maxpv_saddle,label=r'$\log_{10}(\Pi^2/\Pi_0^2)$')
# ax1.text(20,1.08,r'$\omega^2$',color='C0',fontsize=20)
# ax1.text(20,1.22,r'$\Pi^2$',color='C1',fontsize=20)
# ax1.set_xlabel(r'$z_\mathrm{enc}/L_0$')
# ax2.set_xlabel(r'$z_\mathrm{enc}/L_0$')
# ax1.set_ylabel(r'$z_\mathrm{saddle}/z_\mathrm{enc}$')
# ax1.legend(loc='best',fontsize=20)
# ax2.legend(loc='best',fontsize=20)
# plt.tight_layout()
# plt.savefig(opath+'Presentations/pdfs_saddle_height_time.pdf')
# plt.show()
