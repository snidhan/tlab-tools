#!/usr/bin/python3

from ReadStats import Statistics, Pdfs 
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.style
from matplotlib import rc
from matplotlib.colors import LinearSegmentedColormap
from scipy import interpolate


rc('text', usetex=True)
rc('text.latex', preamble=r"\usepackage{fourier}")
rc('font', family='serif')
rc('font', size=20)

opath = '/scratch/local1/m300551/ForKatherine/plots/3D/Re117/5120x1024x5120/'
opath_42 = '/scratch/local1/m300551/ForKatherine/plots/3D/Re042/'
colourmap_path = '/home/mpim/m300551/local/ScientificColourMaps5/'

path_117 = '/scratch/local1/m300551/ForKatherine/qCBL_3D/Re117/'
path_42 = '/scratch/local1/m300551/ForKatherine/qCBL_3D/Re042/'
path_25 = '/scratch/local1/m300551/ForKatherine/qCBL_3D/Re025/'
##########################################################################
# Constants

nu_117 = 1./70000.
nu_42 = 1./25000.
nu_25 = 1./15000.
B0 = 0.005
N = np.sqrt(3)
L0 = (B0/N**3)**0.5
ceps = 0.1
cb=0.1

##########################################################################
# Stats

NS117 = Statistics(path_117+'5120x1024x5120/stats/pdftimes/avg111100-149000.nc')
NS42 = Statistics(path_42+'2560x576x2560/stats/pdftimes/avg36500-47500.nc')
NS42_3 = Statistics(path_42+'2560x896x2560/stats/pdftimes/avg83000-127500.nc')
NS25 = Statistics(path_25+'2560x512x2560/stats/pdftimes/avg17000-21000.nc')

S20_42 = Statistics(path_42+'3072x960x4608-S20/stats/pdftimes/avg66000-84000.nc')

############################################################################
# Pdf

NS117_vortpdf = Pdfs([path_117+'5120x1024x5120/stats/pdfs/pdf111100.LnEnstrophyW_iW_i',path_117+'5120x1024x5120/stats/pdfs/pdf129500.LnEnstrophyW_iW_i',path_117+'5120x1024x5120/stats/pdfs/pdf149000.LnEnstrophyW_iW_i'],path_117+'5120x1024x5120/y.dat')
NS117_pvpdf   = Pdfs([path_117+'5120x1024x5120/stats/pdfs/pdf111100.LnPotentialEnstrophy',path_117+'5120x1024x5120/stats/pdfs/pdf129500.LnPotentialEnstrophy',path_117+'5120x1024x5120/stats/pdfs/pdf149000.LnPotentialEnstrophy'],path_117+'5120x1024x5120/y.dat')
NS42_vortpdf  = Pdfs([path_42+'2560x576x2560/stats/pdfs/pdf36500.LnEnstrophyW_iW_i',path_42+'2560x576x2560/stats/pdfs/pdf41500.LnEnstrophyW_iW_i',path_42+'2560x576x2560/stats/pdfs/pdf47500.LnEnstrophyW_iW_i'],path_42+'2560x576x2560/y.dat')
NS42_pvpdf    = Pdfs([path_42+'2560x576x2560/stats/pdfs/pdf36500.LnPotentialEnstrophy',path_42+'2560x576x2560/stats/pdfs/pdf41500.LnPotentialEnstrophy',path_42+'2560x576x2560/stats/pdfs/pdf47500.LnPotentialEnstrophy'],path_42+'2560x576x2560/y.dat')
NS42_s1gradpdf = Pdfs([path_42+'2560x576x2560/stats/pdfs/pdf36500.LnGradientG_iG_i',path_42+'2560x576x2560/stats/pdfs/pdf41500.LnGradientG_iG_i',path_42+'2560x576x2560/stats/pdfs/pdf47500.LnGradientG_iG_i'],path_42+'2560x576x2560/y.dat')
NS25_vortpdf = Pdfs([path_25+'2560x512x2560/stats/pdfs/pdf17000.LnEnstrophyW_iW_i',path_25+'2560x512x2560/stats/pdfs/pdf19000.LnEnstrophyW_iW_i',path_25+'2560x512x2560/stats/pdfs/pdf21000.LnEnstrophyW_iW_i'],path_25+'2560x512x2560/y.dat')
NS25_pvpdf = Pdfs([path_25+'2560x512x2560/stats/pdfs/pdf17000.LnPotentialEnstrophy',path_25+'2560x512x2560/stats/pdfs/pdf19000.LnPotentialEnstrophy',path_25+'2560x512x2560/stats/pdfs/pdf21000.LnPotentialEnstrophy'],path_25+'2560x512x2560/y.dat')

S20_42_vortpdf = Pdfs([path_42+'3072x960x4608-S20/stats/pdfs/pdf66000.LnEnstrophyW_iW_i',path_42+'3072x960x4608-S20/stats/pdfs/pdf75000.LnEnstrophyW_iW_i',path_42+'3072x960x4608-S20/stats/pdfs/pdf84000.LnEnstrophyW_iW_i'],path_42+'3072x960x4608-S20/y.dat')
S20_42_pvpdf = Pdfs([path_42+'3072x960x4608-S20/stats/pdfs/pdf66000.LnPotentialEnstrophy',path_42+'3072x960x4608-S20/stats/pdfs/pdf75000.LnPotentialEnstrophy',path_42+'3072x960x4608-S20/stats/pdfs/pdf84000.LnPotentialEnstrophy'],path_42+'3072x960x4608-S20/y.dat')
S20_42_s1gradpdf = Pdfs([path_42+'3072x960x4608-S20/stats/pdfs/pdf66000.LnGradientG_iG_i',path_42+'3072x960x4608-S20/stats/pdfs/pdf75000.LnGradientG_iG_i',path_42+'3072x960x4608-S20/stats/pdfs/pdf84000.LnGradientG_iG_i'],path_42+'3072x960x4608-S20/y.dat')

# Create grid on which to interpolate pdfs

NS42_vortpdf_interp_data = Pdfs([path_42+'2560x896x2560/stats/pdfs/pdf102500.LnEnstrophyW_iW_i'],path_42+'2560x896x2560/y.dat')
NS42_pvpdf_interp_data = Pdfs([path_42+'2560x896x2560/stats/pdfs/pdf102500.LnPotentialEnstrophy'],path_42+'2560x896x2560/y.dat')
S20_vortpdf_interp_data = Pdfs([path_42+'3072x960x4608-S20/stats/pdfs/pdf75000.LnEnstrophyW_iW_i'],path_42+'3072x960x4608-S20/y.dat')
S20_pvpdf_interp_data = Pdfs([path_42+'3072x960x4608-S20/stats/pdfs/pdf75000.LnPotentialEnstrophy'],path_42+'3072x960x4608-S20/y.dat')

NS117_vortpdf_interp_data = Pdfs([path_117+'5120x1024x5120/stats/pdfs/pdf149000.LnEnstrophyW_iW_i'],path_117+'5120x1024x5120/y.dat')
NS117_pvpdf_interp_data = Pdfs([path_117+'5120x1024x5120/stats/pdfs/pdf149000.LnPotentialEnstrophy'],path_117+'5120x1024x5120/y.dat')
NS25_vortpdf_interp_data = Pdfs([path_25+'2560x512x2560/stats/pdfs/pdf21000.LnEnstrophyW_iW_i'],path_25+'2560x512x2560/y.dat')
NS25_pvpdf_interp_data = Pdfs([path_25+'2560x512x2560/stats/pdfs/pdf21000.LnPotentialEnstrophy'],path_25+'2560x512x2560/y.dat')

# Interpolate pdfs in y-direction

NS42_vortpdf_interp_y = np.zeros((3,NS42_3.y_len,NS42_vortpdf.nb+2))
for n in range(3):
    for i in range(NS42_vortpdf.nb+2):
        NS42_vortpdf_interp_y[n,:,i] = np.interp(NS42_3.y,NS42.y,NS42_vortpdf.pdf[n,:-1,i])

NS42_pvpdf_interp_y = np.zeros((3,NS42_3.y_len,NS42_pvpdf.nb+2))
for n in range(3):
    for i in range(NS42_pvpdf.nb+2):
        NS42_pvpdf_interp_y[n,:,i] = np.interp(NS42_3.y,NS42.y,NS42_pvpdf.pdf[n,:-1,i])

# Interpolate pdfs in x-direction

NS42_vortpdf_interp = np.zeros((3,NS42_3.y_len,NS42_vortpdf.nb))
for n in range(3):
    for j in range(NS42_3.y_len):
        NS42_vortpdf_interp[n,j,:] = np.interp(NS42_vortpdf_interp_data.xy[0,0,j,:],np.linspace(NS42_vortpdf_interp_y[n,j,NS42_vortpdf.nb],NS42_vortpdf_interp_y[n,j,NS42_vortpdf.nb+1],num=NS42_vortpdf.nb),NS42_vortpdf_interp_y[n,j,:-2])

NS42_pvpdf_interp = np.zeros((3,NS42_3.y_len,NS42_pvpdf.nb))
for n in range(3):
    for j in range(NS42_3.y_len):
        NS42_pvpdf_interp[n,j,:] = np.interp(NS42_pvpdf_interp_data.xy[0,0,j,:],np.linspace(NS42_pvpdf_interp_y[n,j,NS42_pvpdf.nb],NS42_pvpdf_interp_y[n,j,NS42_pvpdf.nb+1],num=NS42_pvpdf.nb),NS42_pvpdf_interp_y[n,j,:-2])

S20_vortpdf_interp = np.zeros((3,S20_42.y_len,S20_42_vortpdf.nb))
for n in range(3):
    for j in range(S20_42.y_len):
        S20_vortpdf_interp[n,j,:] = np.interp(S20_vortpdf_interp_data.xy[0,0,j,:],S20_42_vortpdf.xy[0,n,j,:],S20_42_vortpdf.pdf[n,j,:-2])

S20_pvpdf_interp = np.zeros((3,S20_42.y_len,S20_42_vortpdf.nb))
for n in range(3):
    for j in range(S20_42.y_len):
        S20_pvpdf_interp[n,j,:] = np.interp(S20_pvpdf_interp_data.xy[0,0,j,:],S20_42_pvpdf.xy[0,n,j,:],S20_42_pvpdf.pdf[n,j,:-2])

NS117_vortpdf_interp = np.zeros((3,NS117.y_len,NS117_vortpdf.nb))
for n in range(3):
    for j in range(NS117.y_len):
        NS117_vortpdf_interp[n,j,:] = np.interp(NS117_vortpdf_interp_data.xy[0,0,j,:],NS117_vortpdf.xy[0,n,j,:],NS117_vortpdf.pdf[n,j,:-2])

NS117_pvpdf_interp = np.zeros((3,NS117.y_len,NS117_pvpdf.nb))
for n in range(3):
    for j in range(NS117.y_len):
        NS117_pvpdf_interp[n,j,:] = np.interp(NS117_pvpdf_interp_data.xy[0,0,j,:],NS117_pvpdf.xy[0,n,j,:],NS117_pvpdf.pdf[n,j,:-2])

NS25_vortpdf_interp = np.zeros((3,NS25.y_len,NS25_vortpdf.nb))
for n in range(3):
    for j in range(NS25.y_len):
        NS25_vortpdf_interp[n,j,:] = np.interp(NS25_vortpdf_interp_data.xy[0,0,j,:],NS25_vortpdf.xy[0,n,j,:],NS25_vortpdf.pdf[n,j,:-2])

NS25_pvpdf_interp = np.zeros((3,NS25.y_len,NS25_pvpdf.nb))
for n in range(3):
    for j in range(NS25.y_len):
        NS25_pvpdf_interp[n,j,:] = np.interp(NS25_pvpdf_interp_data.xy[0,0,j,:],NS25_pvpdf.xy[0,n,j,:],NS25_pvpdf.pdf[n,j,:-2])

# Mean of pdfs

NS42_vortpdf_interp_runmean = np.mean(NS42_vortpdf_interp,axis=0)
S20_vortpdf_interp_runmean = np.mean(S20_vortpdf_interp,axis=0)

NS42_pvpdf_interp_runmean = np.mean(NS42_pvpdf_interp,axis=0)
S20_pvpdf_interp_runmean = np.mean(S20_pvpdf_interp,axis=0)

NS117_vortpdf_interp_runmean = np.mean(NS117_vortpdf_interp,axis=0)
NS25_vortpdf_interp_runmean = np.mean(NS25_vortpdf_interp,axis=0)

NS117_pvpdf_interp_runmean = np.mean(NS117_pvpdf_interp,axis=0)
NS25_pvpdf_interp_runmean = np.mean(NS25_pvpdf_interp,axis=0)

# Find where pdf has a maximum at each height

maxvort_NS117 = np.zeros(NS117.y_len)
maxprob_vort_NS117 = np.zeros(NS117.y_len)
for i in range(0,NS117.y_len):
    maxvort_NS117[i] = NS117_vortpdf_interp_data.xy[0,0,i,np.argmax(NS117_vortpdf_interp_runmean[i,:])]
    maxprob_vort_NS117[i] = np.max(NS117_vortpdf_interp_runmean[i,:])
maxvort_NS117 = np.log10(np.exp(maxvort_NS117)/(ceps*B0/nu_117)) 

maxpv_NS117 = np.zeros(NS117.y_len)
maxprob_pv_NS117 = np.zeros(NS117.y_len)
for i in range(0,NS117.y_len):
    maxpv_NS117[i] = NS117_pvpdf_interp_data.xy[0,0,i,np.argmax(NS117_pvpdf_interp_runmean[i,:])]
    maxprob_pv_NS117[i] = np.max(NS117_pvpdf_interp_runmean[i,:])
maxpv_NS117 = np.log10(np.exp(maxpv_NS117)/(cb*ceps*N**6*117**2*(np.mean(NS117.z_enc)/L0)**(-4./3.)))

maxvort_NS42 = np.zeros(NS42_3.y_len)
maxprob_vort_NS42 = np.zeros(NS42_3.y_len)
for i in range(0,NS42_3.y_len):
    maxvort_NS42[i] = NS42_vortpdf_interp_data.xy[0,0,i,np.argmax(NS42_vortpdf_interp_runmean[i,:])]
    maxprob_vort_NS42[i] = np.max(NS42_vortpdf_interp_runmean[i,:])
maxvort_NS42 = np.log10(np.exp(maxvort_NS42)/(ceps*B0/nu_42))

maxpv_NS42 = np.zeros(NS42_3.y_len)
maxprob_pv_NS42 = np.zeros(NS42_3.y_len)
for i in range(0,NS42_3.y_len):
    maxpv_NS42[i] = NS42_pvpdf_interp_data.xy[0,0,i,np.argmax(NS42_pvpdf_interp_runmean[i,:])]
    maxprob_pv_NS42[i] = np.max(NS42_pvpdf_interp_runmean[i,:])
maxpv_NS42 = np.log10(np.exp(maxpv_NS42)/(cb*ceps*N**6*42**2*(np.mean(NS42.z_enc)/L0)**(-4./3.)))
    
maxvort_NS25 = np.zeros(NS25.y_len)
maxprob_vort_NS25 = np.zeros(NS25.y_len)
for i in range(0,NS25.y_len):
    maxvort_NS25[i] = NS25_vortpdf_interp_data.xy[0,0,i,np.argmax(NS25_vortpdf_interp_runmean[i,:])]
    maxprob_vort_NS25[i] = np.max(NS25_vortpdf_interp_runmean[i,:])
maxvort_NS25 = np.log10(np.exp(maxvort_NS25)/(ceps*B0/nu_25))
    
maxpv_NS25 = np.zeros(NS25.y_len)
maxprob_pv_NS25 = np.zeros(NS25.y_len)
for i in range(0,NS25.y_len):
    maxpv_NS25[i] = NS25_pvpdf_interp_data.xy[0,0,i,np.argmax(NS25_pvpdf_interp_runmean[i,:])]
    maxprob_pv_NS25[i] = np.max(NS25_pvpdf_interp_runmean[i,:])
maxpv_NS25 = np.log10(np.exp(maxpv_NS25)/(cb*ceps*N**6*25**2*(np.mean(NS25.z_enc)/L0)**(-4./3.)))

maxvort_S20 = np.zeros(S20_42.y_len)
maxprob_vort_S20 = np.zeros(S20_42.y_len)
for i in range(0,S20_42.y_len):
    maxvort_S20[i] = S20_vortpdf_interp_data.xy[0,0,i,np.argmax(S20_vortpdf_interp_runmean[i,:])]
    maxprob_vort_S20[i] = np.max(S20_vortpdf_interp_runmean[i,:])
maxvort_S20 = np.log10(np.exp(maxvort_S20)/(ceps*B0/nu_42))

maxpv_S20 = np.zeros(S20_42.y_len)
maxprob_pv_S20 = np.zeros(S20_42.y_len)
for i in range(0,S20_42.y_len):
    maxpv_S20[i] = S20_pvpdf_interp_data.xy[0,0,i,np.argmax(S20_pvpdf_interp_runmean[i,:])]
    maxprob_pv_S20[i] = np.max(S20_pvpdf_interp_runmean[i,:])
maxpv_S20 = np.log10(np.exp(maxpv_S20)/(cb*ceps*N**6*42**2*(np.mean(S20_42.z_enc)/L0)**(-4./3.)))

# Find jump in maxvort/maxpv

for i in range(NS42_3.y_len):
    if np.abs(maxvort_NS42[i+1])-np.abs(maxvort_NS42[i]) > 0.2:
        maxit_vort_NS42 = i+1
        break

for i in range(NS42.z_enc_arg[0],NS42_3.y_len):
    if np.abs(maxpv_NS42[i+1])-np.abs(maxpv_NS42[i]) > 0.2:
        maxit_pv_NS42 = i+1
        break
  
for i in range(S20_42.y_len):
    if np.abs(maxvort_S20[i+1])-np.abs(maxvort_S20[i]) > 0.5:
        maxit_vort_S20 = i+1
        break

for i in range(S20_42.y_len):
    if np.abs(maxpv_S20[i+1])-np.abs(maxpv_S20[i]) > 0.7:
        maxit_pv_S20 = i+1
        break

for i in range(NS117.y_len):
    if np.abs(maxvort_NS117[i+1])-np.abs(maxvort_NS117[i]) > 0.7:
        maxit_vort_NS117 = i+1
        break

for i in range(NS117.y_len):
    if np.abs(maxpv_NS117[i+1])-np.abs(maxpv_NS117[i]) > 0.7:
        maxit_pv_NS117 = i+1
        break

for i in range(NS25.y_len):
    if np.abs(maxvort_NS25[i+1])-np.abs(maxvort_NS25[i]) > 0.2:
        maxit_vort_NS25 = i+1
        break

for i in range(NS25.z_enc_arg[0],NS25.y_len):
    if np.abs(maxpv_NS25[i+1])-np.abs(maxpv_NS25[i]) > 0.2:
        maxit_pv_NS25 = i+1
        break
    
# Find saddle as point where maxprob has a minimum

y_vort_117_saddle = NS117.y[np.argmin(maxprob_vort_NS117)]
maxvort_117_saddle = maxvort_NS117[np.argmin(np.abs(y_vort_117_saddle-NS117.y))]
y_vort_117_saddle = y_vort_117_saddle/np.mean(NS117.z_enc)

y_pv_117_saddle = NS117.y[np.argmin(maxprob_pv_NS117[NS117.z_enc_arg[0]:])+NS117.z_enc_arg[0]]
maxpv_117_saddle = maxpv_NS117[np.argmin(np.abs(y_pv_117_saddle-NS117.y))]
y_pv_117_saddle = y_pv_117_saddle/np.mean(NS117.z_enc)


y_vort_NS42_saddle = NS42_3.y[np.argmin(maxprob_vort_NS42[:-100])]
maxvort_NS42_saddle = maxvort_NS42[np.argmin(np.abs(y_vort_NS42_saddle-NS42_3.y))]
y_vort_NS42_saddle = y_vort_NS42_saddle/np.mean(NS42.z_enc)

y_pv_NS42_saddle = NS42_3.y[np.argmin(maxprob_pv_NS42[NS42.z_enc_arg[0]:])+NS42.z_enc_arg[0]]
maxpv_NS42_saddle = maxpv_NS42[np.argmin(np.abs(y_pv_NS42_saddle-NS42_3.y))]
y_pv_NS42_saddle = y_pv_NS42_saddle/np.mean(NS42.z_enc)

y_vort_25_saddle = NS25.y[np.argmin(maxprob_vort_NS25)]
maxvort_25_saddle = maxvort_NS25[np.argmin(np.abs(y_vort_25_saddle-NS25.y))]
y_vort_25_saddle = y_vort_25_saddle/np.mean(NS25.z_enc)

y_pv_25_saddle = NS25.y[np.argmin(maxprob_pv_NS25[NS25.z_enc_arg[0]:])+NS25.z_enc_arg[0]]
maxpv_25_saddle = maxpv_NS25[np.argmin(np.abs(y_pv_25_saddle-NS25.y))]
y_pv_25_saddle = y_pv_25_saddle/np.mean(NS25.z_enc)

y_vort_S20_saddle = S20_42.y[np.argmin(maxprob_vort_S20)]
maxvort_S20_saddle = maxvort_S20[np.argmin(np.abs(y_vort_S20_saddle-S20_42.y))]
y_vort_S20_saddle = y_vort_S20_saddle/np.mean(S20_42.z_enc)

y_pv_S20_saddle = S20_42.y[np.argmin(maxprob_pv_S20[S20_42.z_enc_arg[0]:])+S20_42.z_enc_arg[0]]
maxpv_S20_saddle = maxpv_S20[np.argmin(np.abs(y_pv_S20_saddle-S20_42.y))]
y_pv_S20_saddle = y_pv_S20_saddle/np.mean(S20_42.z_enc)

# Normalisation of y axis

NS117_vortpdf_y_mean = NS117_vortpdf_interp_data.xy[1,0,:,:]/np.mean(NS117.z_enc)
NS117_pvpdf_y_mean = NS117_pvpdf_interp_data.xy[1,0,:,:]/np.mean(NS117.z_enc)

NS42_vortpdf_y_mean = NS42_vortpdf_interp_data.xy[1,0,:,:]/np.mean(NS42.z_enc)
NS42_pvpdf_y_mean = NS42_pvpdf_interp_data.xy[1,0,:,:]/np.mean(NS42.z_enc)
NS42_s1gradpdf_y_mean = NS42_s1gradpdf.xy[1,0,:,:]/np.mean(NS42.z_enc)

NS25_vortpdf_y_mean = NS25_vortpdf_interp_data.xy[1,0,:,:]/np.mean(NS25.z_enc)
NS25_pvpdf_y_mean = NS25_pvpdf_interp_data.xy[1,0,:,:]/np.mean(NS25.z_enc)

S20_42_vortpdf_y_mean = S20_vortpdf_interp_data.xy[1,0,:,:]/np.mean(S20_42.z_enc)
S20_42_pvpdf_y_mean = S20_pvpdf_interp_data.xy[1,0,:,:]/np.mean(S20_42.z_enc)
S20_42_s1gradpdf_y_mean = S20_42_s1gradpdf.xy[1,0,:,:]/np.mean(S20_42.z_enc)

# Normalisation of x axis

NS117_vortpdf_x_mean = np.log10(np.exp(NS117_vortpdf_interp_data.xy[0,0,:,:])/(ceps*B0/nu_117))
NS117_pvpdf_x_mean = np.log10(np.exp(NS117_pvpdf_interp_data.xy[0,0,:,:])/(cb*ceps*N**6*117**2*(np.mean(NS117.z_enc)/L0)**(-4./3.)))

NS42_vortpdf_x_mean = np.log10(np.exp(NS42_vortpdf_interp_data.xy[0,0,:,:])/(ceps*B0/nu_42))
NS42_pvpdf_x_mean = np.log10(np.exp(NS42_pvpdf_interp_data.xy[0,0,:,:])/(cb*ceps*N**6*42**2*(np.mean(NS42.z_enc)/L0)**(-4./3.)))
NS42_s1gradpdf_x_mean = np.log10(np.exp(NS42_s1gradpdf.xy[0,0,:,:])/(cb*N**4*42*(np.mean(NS42.z_enc)/L0)**(-4./3.)))

NS25_vortpdf_x_mean = np.log10(np.exp(NS25_vortpdf_interp_data.xy[0,0,:,:])/(ceps*B0/nu_25))
NS25_pvpdf_x_mean = np.log10(np.exp(NS25_pvpdf_interp_data.xy[0,0,:,:])/(cb*ceps*N**6*25**2*(np.mean(NS25.z_enc)/L0)**(-4./3.)))

S20_42_vortpdf_x_mean = np.log10(np.exp(S20_vortpdf_interp_data.xy[0,0,:,:])/(ceps*B0/nu_42))
S20_42_pvpdf_x_mean = np.log10(np.exp(S20_pvpdf_interp_data.xy[0,0,:,:])/(cb*ceps*N**6*42**2*(np.mean(S20_42.z_enc)/L0)**(-4./3.)))
S20_42_s1gradpdf_x_mean = np.log10(np.exp(S20_42_s1gradpdf.xy[0,0,:,:])/(cb*N**4*42*(np.mean(S20_42.z_enc)/L0)**(-4./3.)))

#####################################################################
# Colourmaps

imola_data = np.loadtxt(colourmap_path+'imola/imola.txt')
imola_map = LinearSegmentedColormap.from_list('imola',imola_data)

davos_data = np.loadtxt(colourmap_path+'davos/davos.txt')
davos_map = LinearSegmentedColormap.from_list('davos',davos_data)
#####################################################################
# Plot

f, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,sharex='col',sharey='all',figsize=(10,10))
ax1.tick_params(bottom=True,top=True,left=True,right=True)
ax2.tick_params(bottom=True,top=True,left=True,right=True)
ax3.tick_params(bottom=True,top=True,left=True,right=True)
ax4.tick_params(bottom=True,top=True,left=True,right=True)
ax1.set_xlim(-4,2) 
ax2.set_xlim(-7,2)
ax1.set_ylim(0,1.6)
ax1.set_yticks([0,0.25,0.5,0.75,1,1.25,1.5])
ax2.set_xticks([-6,-4,-2,0,2])
cs1 = ax1.contourf(NS117_vortpdf_x_mean,NS117_vortpdf_y_mean,NS117_vortpdf_interp_runmean,cmap=imola_map,levels=np.linspace(0,0.4,9))
ax1.plot(maxvort_NS117[:maxit_vort_NS117],NS117.y[:maxit_vort_NS117]/np.mean(NS117.z_enc),'k',lw=1)
ax1.plot(maxvort_NS117[maxit_vort_NS117:],NS117.y[maxit_vort_NS117:]/np.mean(NS117.z_enc),'k',lw=1)
ax1.scatter((maxvort_NS117[maxit_vort_NS117-1]+maxvort_NS117[maxit_vort_NS117])/2,NS117.y[maxit_vort_NS117]/np.mean(NS117.z_enc),100,color='k',marker='*')
ax1.axhline(np.mean(NS117.z_ig/NS117.z_enc),0,0.05,color='C1',linewidth=2)
ax1.axhline(np.mean(NS117.z_is/NS117.z_enc),0,0.05,color='C1',linewidth=2)
ax1.axhline(np.mean(NS117.z_if/NS117.z_enc),0,0.05,color='C1',linewidth=2)
cs2 = ax2.contourf(NS117_pvpdf_x_mean,NS117_pvpdf_y_mean,NS117_pvpdf_interp_runmean,cmap=imola_map,levels=np.linspace(0,0.24,9))
ax2.plot(maxpv_NS117[:maxit_pv_NS117],NS117.y[:maxit_pv_NS117]/np.mean(NS117.z_enc),'k',lw=1)
ax2.plot(maxpv_NS117[maxit_pv_NS117:],NS117.y[maxit_pv_NS117:]/np.mean(NS117.z_enc),'k',lw=1)
ax2.scatter((maxpv_NS117[maxit_pv_NS117-1]+maxpv_NS117[maxit_pv_NS117])/2,NS117.y[maxit_pv_NS117]/np.mean(NS117.z_enc),100,color='k',marker='*')
ax2.axhline(np.mean(NS117.z_ig/NS117.z_enc),0,0.05,color='C1',linewidth=2)
ax2.axhline(np.mean(NS117.z_is/NS117.z_enc),0,0.05,color='C1',linewidth=2)
ax2.axhline(np.mean(NS117.z_if/NS117.z_enc),0,0.05,color='C1',linewidth=2)
cs3 = ax3.contourf(NS25_vortpdf_x_mean,NS25_vortpdf_y_mean,NS25_vortpdf_interp_runmean,cmap=imola_map,levels=np.linspace(0,0.4,9))
ax3.plot(maxvort_NS25[:maxit_vort_NS25],NS25.y[:maxit_vort_NS25]/np.mean(NS25.z_enc),'k',lw=1)
ax3.plot(maxvort_NS25[maxit_vort_NS25:],NS25.y[maxit_vort_NS25:]/np.mean(NS25.z_enc),'k',lw=1)
ax3.scatter((maxvort_NS25[maxit_vort_NS25-1]+maxvort_NS25[maxit_vort_NS25])/2,NS25.y[maxit_vort_NS25]/np.mean(NS25.z_enc),100,color='k',marker='*')
ax3.axhline(np.mean(NS25.z_ig/NS25.z_enc),0,0.05,color='C1',linewidth=2)
ax3.axhline(np.mean(NS25.z_is/NS25.z_enc),0,0.05,color='C1',linewidth=2)
ax3.axhline(np.mean(NS25.z_if/NS25.z_enc),0,0.05,color='C1',linewidth=2)
cs4 = ax4.contourf(NS25_pvpdf_x_mean,NS25_pvpdf_y_mean,NS25_pvpdf_interp_runmean,cmap=imola_map,levels=np.linspace(0,0.24,9))
ax4.plot(maxpv_NS25[:maxit_pv_NS25],NS25.y[:maxit_pv_NS25]/np.mean(NS25.z_enc),'k',lw=1)
ax4.plot(maxpv_NS25[maxit_pv_NS25:],NS25.y[maxit_pv_NS25:]/np.mean(NS25.z_enc),'k',lw=1)
ax4.scatter((maxpv_NS25[maxit_pv_NS25-1]+maxpv_NS25[maxit_pv_NS25])/2,NS25.y[maxit_pv_NS25]/np.mean(NS25.z_enc),100,color='k',marker='*')
ax4.axhline(np.mean(NS25.z_ig/NS25.z_enc),0,0.05,color='C1',linewidth=2)
ax4.axhline(np.mean(NS25.z_is/NS25.z_enc),0,0.05,color='C1',linewidth=2)
ax4.axhline(np.mean(NS25.z_if/NS25.z_enc),0,0.05,color='C1',linewidth=2)
ax3.set_xlabel(r'$\log_{10}(\omega^2/\omega_0^2)$')
ax4.set_xlabel(r'$\log_{10}(\Pi^2/\Pi_0^2)$')
ax1.set_ylabel(r'$z/z_\mathrm{enc}$')
ax3.set_ylabel(r'$z/z_\mathrm{enc}$')
ax1.set_title(r'(a) $Re_0=117$',fontsize=20,loc='left')
ax2.set_title(r'(b) $Re_0=117$',fontsize=20,loc='left')
ax3.set_title(r'(c) $Re_0=25$',fontsize=20,loc='left')
ax4.set_title(r'(d) $Re_0=25$',fontsize=20,loc='left')
plt.colorbar(cs1,ax=ax1)
plt.colorbar(cs2,ax=ax2)
plt.colorbar(cs3,ax=ax3)
plt.colorbar(cs4,ax=ax4)
plt.tight_layout()
plt.savefig(opath+'pdfs_vort_pv_subplots_Re_25_117_timeavg_interp.pdf')
plt.show()


f, (ax1,ax2) = plt.subplots(1,2,sharex='all',sharey='all',figsize=(10,5))
ax1.tick_params(bottom=True,top=True,left=True,right=True)
ax2.tick_params(bottom=True,top=True,left=True,right=True)
ax1.set_xlim(-4,2) 
ax1.set_ylim(0,1.6)
ax1.set_yticks([0,0.25,0.5,0.75,1,1.25,1.5])
ax1.set_xticks([-4,-3,-2,-1,0,1,2])
cs1 = ax1.contourf(NS117_vortpdf_x_mean,NS117_vortpdf_y_mean,NS117_vortpdf_interp_runmean,cmap=imola_map,levels=np.linspace(0,0.4,9))
ax1.plot(maxvort_NS117[:maxit_vort_NS117],NS117.y[:maxit_vort_NS117]/np.mean(NS117.z_enc),'k',lw=1)
ax1.plot(maxvort_NS117[maxit_vort_NS117:],NS117.y[maxit_vort_NS117:]/np.mean(NS117.z_enc),'k',lw=1)
ax1.scatter((maxvort_NS117[maxit_vort_NS117-1]+maxvort_NS117[maxit_vort_NS117])/2,NS117.y[maxit_vort_NS117]/np.mean(NS117.z_enc),100,color='k',marker='*')
ax1.axhline(np.mean(NS117.z_ig/NS117.z_enc),0,0.05,color='C1',linewidth=2)
ax1.axhline(np.mean(NS117.z_is/NS117.z_enc),0,0.05,color='C1',linewidth=2)
ax1.axhline(np.mean(NS117.z_if/NS117.z_enc),0,0.05,color='C1',linewidth=2)
cs2 = ax2.contourf(NS25_vortpdf_x_mean,NS25_vortpdf_y_mean,NS25_vortpdf_interp_runmean,cmap=imola_map,levels=np.linspace(0,0.4,9))
ax2.plot(maxvort_NS25[:maxit_vort_NS25],NS25.y[:maxit_vort_NS25]/np.mean(NS25.z_enc),'k',lw=1)
ax2.plot(maxvort_NS25[maxit_vort_NS25:],NS25.y[maxit_vort_NS25:]/np.mean(NS25.z_enc),'k',lw=1)
ax2.scatter((maxvort_NS25[maxit_vort_NS25-1]+maxvort_NS25[maxit_vort_NS25])/2,NS25.y[maxit_vort_NS25]/np.mean(NS25.z_enc),100,color='k',marker='*')
ax2.axhline(np.mean(NS25.z_ig/NS25.z_enc),0,0.05,color='C1',linewidth=2)
ax2.axhline(np.mean(NS25.z_is/NS25.z_enc),0,0.05,color='C1',linewidth=2)
ax2.axhline(np.mean(NS25.z_if/NS25.z_enc),0,0.05,color='C1',linewidth=2)
ax1.set_xlabel(r'$\log_{10}(\omega^2/\omega_0^2)$')
ax2.set_xlabel(r'$\log_{10}(\omega^2/\omega_0^2)$')
ax1.set_ylabel(r'$z/z_\mathrm{enc}$')
ax1.set_title(r'(a) $Re_0=117$',fontsize=20,loc='left')
ax2.set_title(r'(b) $Re_0=25$',fontsize=20,loc='left')
cbar_ax = f.add_axes([0.3,0.1,0.5,0.03])
cbar = f.colorbar(cs1,cax=cbar_ax,orientation='horizontal')
plt.tight_layout(rect=[0,0.1,1,1],h_pad=2)
plt.savefig(opath+'pdfs_vort_subplots_Re_25_117_timeavg_interp.pdf',bbox_inches='tight')
plt.show()


f, (ax1,ax2) = plt.subplots(1,2,sharey='all',figsize=(10,5))
ax1.tick_params(bottom=True,top=True,left=True,right=True)
ax2.tick_params(bottom=True,top=True,left=True,right=True)
ax1.set_ylim(0,1.6)
ax1.set_xlim(-3,3)
ax2.set_xlim(-3,3)
cs1 = ax1.contourf(NS42_s1gradpdf.xy[0,1,:,:],NS42_s1gradpdf.xy[1,1,:,:],NS42_s1gradpdf.pdf[1,:NS42_s1gradpdf.ny,:NS42_s1gradpdf.nb],cmap=imola_map,levels=np.linspace(0,0.4,11),extend='max')
cs2 = ax2.contourf(S20_42_s1gradpdf.xy[0,1,:,:],S20_42_s1gradpdf.xy[1,1,:,:],S20_42_s1gradpdf.pdf[1,:S20_42_s1gradpdf.ny,:S20_42_s1gradpdf.nb],cmap=imola_map,levels=np.linspace(0,0.4,11),extend='max')
ax1.set_xlabel(r'$\log_{10}(|\nabla b|^2/|\nabla b|_0^2)$')
ax2.set_xlabel(r'$\log_{10}(|\nabla b|^2/|\nabla b|_0^2)$')
ax1.set_ylabel(r'$z/z_\mathrm{enc}$')
ax1.set_title(r'(a) $Fr_0=0$',loc='left',fontsize=20)
ax2.set_title(r'(b) $Fr_0=20$',loc='left',fontsize=20)
plt.colorbar(cs1,ax=ax1)
plt.colorbar(cs2,ax=ax2)
plt.tight_layout()
plt.savefig(opath_42+'pdfs_s1grad_S20_S0_20.pdf')
plt.show()



f, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,sharex='col',sharey='all',figsize=(10,10))
ax1.tick_params(bottom=True,top=True,left=True,right=True)
ax2.tick_params(bottom=True,top=True,left=True,right=True)
ax3.tick_params(bottom=True,top=True,left=True,right=True)
ax4.tick_params(bottom=True,top=True,left=True,right=True)
ax1.set_xlim(-4,2) 
ax2.set_xlim(-7,2)
ax1.set_ylim(0,1.6)
ax1.set_yticks([0,0.25,0.5,0.75,1,1.25,1.5])
ax2.set_xticks([-6,-4,-2,0,2])
cs1 = ax1.contourf(NS42_vortpdf_x_mean,NS42_vortpdf_y_mean,NS42_vortpdf_interp_runmean,cmap=imola_map,levels=np.linspace(0,0.4,9))
ax1.plot(maxvort_NS42[:maxit_vort_NS42],NS42_3.y[:maxit_vort_NS42]/np.mean(NS42.z_enc),'k',lw=1)
ax1.plot(maxvort_NS42[maxit_vort_NS42:],NS42_3.y[maxit_vort_NS42:]/np.mean(NS42.z_enc),'k',lw=1)
ax1.scatter((maxvort_NS42[maxit_vort_NS42-1]+maxvort_NS42[maxit_vort_NS42])/2,NS42_3.y[maxit_vort_NS42]/np.mean(NS42.z_enc),100,color='k',marker='*')
ax1.axhline(np.mean(NS42.z_ig/NS42.z_enc),0,0.05,color='C1',linewidth=2)
ax1.axhline(np.mean(NS42.z_is/NS42.z_enc),0,0.05,color='C1',linewidth=2)
ax1.axhline(np.mean(NS42.z_if/NS42.z_enc),0,0.05,color='C1',linewidth=2)
cs2 = ax2.contourf(NS42_pvpdf_x_mean,NS42_pvpdf_y_mean,NS42_pvpdf_interp_runmean,cmap=imola_map,levels=np.linspace(0,0.24,9))
ax2.plot(maxpv_NS42[:maxit_pv_NS42],NS42_3.y[:maxit_pv_NS42]/np.mean(NS42.z_enc),'k',lw=1)
ax2.plot(maxpv_NS42[maxit_pv_NS42:],NS42_3.y[maxit_pv_NS42:]/np.mean(NS42.z_enc),'k',lw=1)
ax2.scatter((maxpv_NS42[maxit_pv_NS42-1]+maxpv_NS42[maxit_pv_NS42])/2,NS42_3.y[maxit_pv_NS42]/np.mean(NS42.z_enc),100,color='k',marker='*')
ax2.axhline(np.mean(NS42.z_ig/NS42.z_enc),0,0.05,color='C1',linewidth=2)
ax2.axhline(np.mean(NS42.z_is/NS42.z_enc),0,0.05,color='C1',linewidth=2)
ax2.axhline(np.mean(NS42.z_if/NS42.z_enc),0,0.05,color='C1',linewidth=2)
cs3 = ax3.contourf(S20_42_vortpdf_x_mean,S20_42_vortpdf_y_mean,S20_vortpdf_interp_runmean,cmap=imola_map,levels=np.linspace(0,0.4,9))
ax3.plot(maxvort_S20[:maxit_vort_S20],S20_42.y[:maxit_vort_S20]/np.mean(S20_42.z_enc),'k',lw=1)
ax3.plot(maxvort_S20[maxit_vort_S20:],S20_42.y[maxit_vort_S20:]/np.mean(S20_42.z_enc),'k',lw=1)
ax3.scatter((maxvort_S20[maxit_vort_S20-1]+maxvort_S20[maxit_vort_S20])/2,S20_42.y[maxit_vort_S20]/np.mean(S20_42.z_enc),100,color='k',marker='*')
ax3.axhline(np.mean(S20_42.z_ig/S20_42.z_enc),0,0.05,color='C1',linewidth=2)
ax3.axhline(np.mean(S20_42.z_is/S20_42.z_enc),0,0.05,color='C1',linewidth=2)
ax3.axhline(np.mean(S20_42.z_if/S20_42.z_enc),0,0.05,color='C1',linewidth=2)
cs4 = ax4.contourf(S20_42_pvpdf_x_mean,S20_42_pvpdf_y_mean,S20_pvpdf_interp_runmean,cmap=imola_map,levels=np.linspace(0,0.24,9))
ax4.plot(maxpv_S20[:maxit_pv_S20],S20_42.y[:maxit_pv_S20]/np.mean(S20_42.z_enc),'k',lw=1)
ax4.plot(maxpv_S20[maxit_pv_S20:],S20_42.y[maxit_pv_S20:]/np.mean(S20_42.z_enc),'k',lw=1)
ax4.scatter((maxpv_S20[maxit_pv_S20-1]+maxpv_S20[maxit_pv_S20])/2,S20_42.y[maxit_pv_S20]/np.mean(S20_42.z_enc),100,color='k',marker='*')
ax4.axhline(np.mean(S20_42.z_ig/S20_42.z_enc),0,0.05,color='C1',linewidth=2)
ax4.axhline(np.mean(S20_42.z_is/S20_42.z_enc),0,0.05,color='C1',linewidth=2)
ax4.axhline(np.mean(S20_42.z_if/S20_42.z_enc),0,0.05,color='C1',linewidth=2)
ax1.set_ylabel(r'$z/z_\mathrm{enc}$')
ax3.set_ylabel(r'$z/z_\mathrm{enc}$')
ax3.set_xlabel(r'$\log_{10}(\omega^2/\omega_0^2)$')
ax4.set_xlabel(r'$\log_{10}(\Pi^2/\Pi_0^2)$')
ax1.set_title(r'(a)$Fr_0=0$',fontsize=20,loc='left')
ax2.set_title(r'(b)$Fr_0=0$',fontsize=20,loc='left')
ax3.set_title(r'(c)$Fr_0=20$',fontsize=20,loc='left')
ax4.set_title(r'(d)$Fr_0=20$',fontsize=20,loc='left')
plt.colorbar(cs1,ax=ax1)
plt.colorbar(cs2,ax=ax2)
plt.colorbar(cs3,ax=ax3)
plt.colorbar(cs4,ax=ax4)
plt.tight_layout()
plt.savefig(opath_42+'pdfs_vort_pv_subplots_S20_S0_timeavg_interpxy.pdf')
plt.show()

f, (ax1,ax2) = plt.subplots(1,2,sharex='all',sharey='all',figsize=(10,5))
ax1.tick_params(bottom=True,top=True,left=True,right=True)
ax2.tick_params(bottom=True,top=True,left=True,right=True)
ax1.set_xlim(-4,2) 
ax1.set_ylim(0,1.6)
ax1.set_xticks([-4,-3,-2,-1,0,1,2])
ax1.set_yticks([0,0.25,0.5,0.75,1,1.25,1.5])
cs1 = ax1.contourf(NS42_vortpdf_x_mean,NS42_vortpdf_y_mean,NS42_vortpdf_interp_runmean,cmap=imola_map,levels=np.linspace(0,0.4,9))
ax1.plot(maxvort_NS42[:maxit_vort_NS42],NS42_3.y[:maxit_vort_NS42]/np.mean(NS42.z_enc),'k',lw=1)
ax1.plot(maxvort_NS42[maxit_vort_NS42:],NS42_3.y[maxit_vort_NS42:]/np.mean(NS42.z_enc),'k',lw=1)
ax1.scatter((maxvort_NS42[maxit_vort_NS42-1]+maxvort_NS42[maxit_vort_NS42])/2,NS42_3.y[maxit_vort_NS42]/np.mean(NS42.z_enc),100,color='k',marker='*')
ax1.axhline(np.mean(NS42.z_ig/NS42.z_enc),0,0.05,color='C1',linewidth=2)
ax1.axhline(np.mean(NS42.z_is/NS42.z_enc),0,0.05,color='C1',linewidth=2)
ax1.axhline(np.mean(NS42.z_if/NS42.z_enc),0,0.05,color='C1',linewidth=2)
cs2 = ax2.contourf(S20_42_vortpdf_x_mean,S20_42_vortpdf_y_mean,S20_vortpdf_interp_runmean,cmap=imola_map,levels=np.linspace(0,0.4,9))
ax2.plot(maxvort_S20[:maxit_vort_S20],S20_42.y[:maxit_vort_S20]/np.mean(S20_42.z_enc),'k',lw=1)
ax2.plot(maxvort_S20[maxit_vort_S20:],S20_42.y[maxit_vort_S20:]/np.mean(S20_42.z_enc),'k',lw=1)
ax2.scatter((maxvort_S20[maxit_vort_S20-1]+maxvort_S20[maxit_vort_S20])/2,S20_42.y[maxit_vort_S20]/np.mean(S20_42.z_enc),100,color='k',marker='*')
ax2.axhline(np.mean(S20_42.z_ig/S20_42.z_enc),0,0.05,color='C1',linewidth=2)
ax2.axhline(np.mean(S20_42.z_is/S20_42.z_enc),0,0.05,color='C1',linewidth=2)
ax2.axhline(np.mean(S20_42.z_if/S20_42.z_enc),0,0.05,color='C1',linewidth=2)
ax1.set_ylabel(r'$z/z_\mathrm{enc}$')
ax1.set_xlabel(r'$\log_{10}(\omega^2/\omega_0^2)$')
ax2.set_xlabel(r'$\log_{10}(\omega^2/\omega_0^2)$')
ax1.set_title(r'(a)$Fr_0=0$',fontsize=20,loc='left')
ax2.set_title(r'(b)$Fr_0=20$',fontsize=20,loc='left')
cbar_ax = f.add_axes([0.3,0.1,0.5,0.03])
cbar = f.colorbar(cs1,cax=cbar_ax,orientation='horizontal')
plt.tight_layout(rect=[0,0.1,1,1],h_pad=2)
plt.savefig(opath_42+'pdfs_vort_subplots_S20_S0_timeavg_interpxy.pdf',bbox_inches='tight')
plt.show()



#Presentations

# f, (ax1,ax2) = plt.subplots(1,2,sharey='all',figsize=(10,5))
# ax1.set_xlim(-4,2) 
# ax2.set_xlim(-7,2)
# ax1.set_ylim(0,1.6)
# ax2.set_xticks([-6,-4,-2,0,2])
# cs1 = ax1.contourf(NS42_vortpdf.xy[0,0,:,:],NS42_vortpdf.xy[1,0,:,:],NS42_vortpdf.pdf_timeavg[:NS42_vortpdf.ny,:NS42_vortpdf.nb],extend='max')
# ax1.plot(maxvort_42,NS42.y/NS42.z_enc,'r')
# ax1.scatter(maxvort_42_saddle,y_vort_42_saddle,100,color='k',marker='*')
# cs2 = ax2.contourf(NS42_pvpdf.xy[0,0,:,:],NS42_pvpdf.xy[1,0,:,:],NS42_pvpdf.pdf_timeavg[:NS42_pvpdf.ny,:NS42_pvpdf.nb],levels=np.linspace(0,0.24,9),extend='max')
# ax2.plot(maxpv_42,NS42.y/NS42.z_enc,'r')
# ax2.scatter(maxpv_42_saddle,y_pv_42_saddle,100,color='k',marker='*')
# ax1.set_xlabel(r'$\log_{10}(\omega^2/\omega_0^2)$')
# ax2.set_xlabel(r'$\log_{10}(\Pi^2/\Pi_0^2)$')
# ax1.set_ylabel(r'$z/z_\mathrm{enc}$')
# plt.colorbar(cs1,ax=ax1)
# plt.colorbar(cs2,ax=ax2)
# plt.tight_layout()
# plt.savefig('/scratch/local1/m300551/ForKatherine/plots/3D/Re042/Presentations/pdfs_vort_pv_saddle.pdf')
# plt.show()

