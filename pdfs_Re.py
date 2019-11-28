#!/usr/bin/python3

from ReadStats import Statistics, Pdfs 
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.style
from matplotlib import rc
from matplotlib.colors import LinearSegmentedColormap


rc('text', usetex=True)
rc('text.latex', preamble=r"\usepackage{fourier}")
rc('font', family='serif')
rc('font', size=20)

opath = '/home/mpim/m300551/Figures/JAS2020/'
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
NS42_vortpdf  = Pdfs([path_42+'2560x576x2560/stats/pdfs/pdf36500.LnEnstrophyW_iW_i',path_42+'2560x576x2560/stats/pdfs/pdf41500.LnEnstrophyW_iW_i',path_42+'2560x576x2560/stats/pdfs/pdf47500.LnEnstrophyW_iW_i'],path_42+'2560x576x2560/y.dat')
NS25_vortpdf = Pdfs([path_25+'2560x512x2560/stats/pdfs/pdf17000.LnEnstrophyW_iW_i',path_25+'2560x512x2560/stats/pdfs/pdf19000.LnEnstrophyW_iW_i',path_25+'2560x512x2560/stats/pdfs/pdf21000.LnEnstrophyW_iW_i'],path_25+'2560x512x2560/y.dat')
S20_42_vortpdf = Pdfs([path_42+'3072x960x4608-S20/stats/pdfs/pdf66000.LnEnstrophyW_iW_i',path_42+'3072x960x4608-S20/stats/pdfs/pdf75000.LnEnstrophyW_iW_i',path_42+'3072x960x4608-S20/stats/pdfs/pdf84000.LnEnstrophyW_iW_i'],path_42+'3072x960x4608-S20/y.dat')

# Create grid on which to interpolate pdfs

NS42_vortpdf_interp_data = Pdfs([path_42+'2560x896x2560/stats/pdfs/pdf102500.LnEnstrophyW_iW_i'],path_42+'2560x896x2560/y.dat')
S20_vortpdf_interp_data = Pdfs([path_42+'3072x960x4608-S20/stats/pdfs/pdf75000.LnEnstrophyW_iW_i'],path_42+'3072x960x4608-S20/y.dat')

NS117_vortpdf_interp_data = Pdfs([path_117+'5120x1024x5120/stats/pdfs/pdf149000.LnEnstrophyW_iW_i'],path_117+'5120x1024x5120/y.dat')
NS25_vortpdf_interp_data = Pdfs([path_25+'2560x512x2560/stats/pdfs/pdf21000.LnEnstrophyW_iW_i'],path_25+'2560x512x2560/y.dat')

# Interpolate pdfs in y-direction

NS42_vortpdf_interp_y = np.zeros((3,NS42_3.y_len,NS42_vortpdf.nb+2))
for n in range(3):
    for i in range(NS42_vortpdf.nb+2):
        NS42_vortpdf_interp_y[n,:,i] = np.interp(NS42_3.y,NS42.y,NS42_vortpdf.pdf[n,:-1,i])

# Interpolate pdfs in x-direction

NS42_vortpdf_interp = np.zeros((3,NS42_3.y_len,NS42_vortpdf.nb))
for n in range(3):
    for j in range(NS42_3.y_len):
        NS42_vortpdf_interp[n,j,:] = np.interp(NS42_vortpdf_interp_data.xy[0,0,j,:],np.linspace(NS42_vortpdf_interp_y[n,j,NS42_vortpdf.nb],NS42_vortpdf_interp_y[n,j,NS42_vortpdf.nb+1],num=NS42_vortpdf.nb),NS42_vortpdf_interp_y[n,j,:-2])

S20_vortpdf_interp = np.zeros((3,S20_42.y_len,S20_42_vortpdf.nb))
for n in range(3):
    for j in range(S20_42.y_len):
        S20_vortpdf_interp[n,j,:] = np.interp(S20_vortpdf_interp_data.xy[0,0,j,:],S20_42_vortpdf.xy[0,n,j,:],S20_42_vortpdf.pdf[n,j,:-2])

NS117_vortpdf_interp = np.zeros((3,NS117.y_len,NS117_vortpdf.nb))
for n in range(3):
    for j in range(NS117.y_len):
        NS117_vortpdf_interp[n,j,:] = np.interp(NS117_vortpdf_interp_data.xy[0,0,j,:],NS117_vortpdf.xy[0,n,j,:],NS117_vortpdf.pdf[n,j,:-2])

NS25_vortpdf_interp = np.zeros((3,NS25.y_len,NS25_vortpdf.nb))
for n in range(3):
    for j in range(NS25.y_len):
        NS25_vortpdf_interp[n,j,:] = np.interp(NS25_vortpdf_interp_data.xy[0,0,j,:],NS25_vortpdf.xy[0,n,j,:],NS25_vortpdf.pdf[n,j,:-2])

# Mean of pdfs

NS42_vortpdf_interp_runmean = np.mean(NS42_vortpdf_interp,axis=0)
S20_vortpdf_interp_runmean = np.mean(S20_vortpdf_interp,axis=0)

NS117_vortpdf_interp_runmean = np.mean(NS117_vortpdf_interp,axis=0)
NS25_vortpdf_interp_runmean = np.mean(NS25_vortpdf_interp,axis=0)

# Find where pdf has a maximum at each height

maxvort_NS117 = np.zeros(NS117.y_len)
maxprob_vort_NS117 = np.zeros(NS117.y_len)
for i in range(0,NS117.y_len):
    maxvort_NS117[i] = NS117_vortpdf_interp_data.xy[0,0,i,np.argmax(NS117_vortpdf_interp_runmean[i,:])]
    maxprob_vort_NS117[i] = np.max(NS117_vortpdf_interp_runmean[i,:])
maxvort_NS117 = np.log10(np.exp(maxvort_NS117)/(ceps*B0/nu_117)) 

maxvort_NS42 = np.zeros(NS42_3.y_len)
maxprob_vort_NS42 = np.zeros(NS42_3.y_len)
for i in range(0,NS42_3.y_len):
    maxvort_NS42[i] = NS42_vortpdf_interp_data.xy[0,0,i,np.argmax(NS42_vortpdf_interp_runmean[i,:])]
    maxprob_vort_NS42[i] = np.max(NS42_vortpdf_interp_runmean[i,:])
maxvort_NS42 = np.log10(np.exp(maxvort_NS42)/(ceps*B0/nu_42))
    
maxvort_NS25 = np.zeros(NS25.y_len)
maxprob_vort_NS25 = np.zeros(NS25.y_len)
for i in range(0,NS25.y_len):
    maxvort_NS25[i] = NS25_vortpdf_interp_data.xy[0,0,i,np.argmax(NS25_vortpdf_interp_runmean[i,:])]
    maxprob_vort_NS25[i] = np.max(NS25_vortpdf_interp_runmean[i,:])
maxvort_NS25 = np.log10(np.exp(maxvort_NS25)/(ceps*B0/nu_25))

maxvort_S20 = np.zeros(S20_42.y_len)
maxprob_vort_S20 = np.zeros(S20_42.y_len)
for i in range(0,S20_42.y_len):
    maxvort_S20[i] = S20_vortpdf_interp_data.xy[0,0,i,np.argmax(S20_vortpdf_interp_runmean[i,:])]
    maxprob_vort_S20[i] = np.max(S20_vortpdf_interp_runmean[i,:])
maxvort_S20 = np.log10(np.exp(maxvort_S20)/(ceps*B0/nu_42))

# Find jump in maxvort

for i in range(NS42_3.y_len):
    if np.abs(maxvort_NS42[i+1])-np.abs(maxvort_NS42[i]) > 0.2:
        maxit_vort_NS42 = i+1
        break
  
for i in range(S20_42.y_len):
    if np.abs(maxvort_S20[i+1])-np.abs(maxvort_S20[i]) > 0.5:
        maxit_vort_S20 = i+1
        break

for i in range(NS117.y_len):
    if np.abs(maxvort_NS117[i+1])-np.abs(maxvort_NS117[i]) > 0.7:
        maxit_vort_NS117 = i+1
        break
    
for i in range(NS25.y_len):
    if np.abs(maxvort_NS25[i+1])-np.abs(maxvort_NS25[i]) > 0.2:
        maxit_vort_NS25 = i+1
        break

# Find saddle as point where maxprob has a minimum

y_vort_117_saddle = NS117.y[np.argmin(maxprob_vort_NS117)]
maxvort_117_saddle = maxvort_NS117[np.argmin(np.abs(y_vort_117_saddle-NS117.y))]
y_vort_117_saddle = y_vort_117_saddle/np.mean(NS117.z_enc)

y_vort_NS42_saddle = NS42_3.y[np.argmin(maxprob_vort_NS42[:-100])]
maxvort_NS42_saddle = maxvort_NS42[np.argmin(np.abs(y_vort_NS42_saddle-NS42_3.y))]
y_vort_NS42_saddle = y_vort_NS42_saddle/np.mean(NS42.z_enc)

y_vort_25_saddle = NS25.y[np.argmin(maxprob_vort_NS25)]
maxvort_25_saddle = maxvort_NS25[np.argmin(np.abs(y_vort_25_saddle-NS25.y))]
y_vort_25_saddle = y_vort_25_saddle/np.mean(NS25.z_enc)

y_vort_S20_saddle = S20_42.y[np.argmin(maxprob_vort_S20)]
maxvort_S20_saddle = maxvort_S20[np.argmin(np.abs(y_vort_S20_saddle-S20_42.y))]
y_vort_S20_saddle = y_vort_S20_saddle/np.mean(S20_42.z_enc)

# Normalisation of y axis

NS117_vortpdf_y_mean = NS117_vortpdf_interp_data.xy[1,0,:,:]/np.mean(NS117.z_enc)
NS42_vortpdf_y_mean = NS42_vortpdf_interp_data.xy[1,0,:,:]/np.mean(NS42.z_enc)
NS25_vortpdf_y_mean = NS25_vortpdf_interp_data.xy[1,0,:,:]/np.mean(NS25.z_enc)
S20_42_vortpdf_y_mean = S20_vortpdf_interp_data.xy[1,0,:,:]/np.mean(S20_42.z_enc)

# Normalisation of x axis

NS117_vortpdf_x_mean = np.log10(np.exp(NS117_vortpdf_interp_data.xy[0,0,:,:])/(ceps*B0/nu_117))
NS42_vortpdf_x_mean = np.log10(np.exp(NS42_vortpdf_interp_data.xy[0,0,:,:])/(ceps*B0/nu_42))
NS25_vortpdf_x_mean = np.log10(np.exp(NS25_vortpdf_interp_data.xy[0,0,:,:])/(ceps*B0/nu_25))
S20_42_vortpdf_x_mean = np.log10(np.exp(S20_vortpdf_interp_data.xy[0,0,:,:])/(ceps*B0/nu_42))

#####################################################################
# Colourmaps

imola_data = np.loadtxt(colourmap_path+'imola/imola.txt')
imola_map = LinearSegmentedColormap.from_list('imola',imola_data)

#####################################################################
# Plot

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
ax1.axhline(np.mean(NS117.z_if/NS117.z_enc),0,0.05,color='C1',linewidth=2)
cs2 = ax2.contourf(NS25_vortpdf_x_mean,NS25_vortpdf_y_mean,NS25_vortpdf_interp_runmean,cmap=imola_map,levels=np.linspace(0,0.4,9))
ax2.plot(maxvort_NS25[:maxit_vort_NS25],NS25.y[:maxit_vort_NS25]/np.mean(NS25.z_enc),'k',lw=1)
ax2.plot(maxvort_NS25[maxit_vort_NS25:],NS25.y[maxit_vort_NS25:]/np.mean(NS25.z_enc),'k',lw=1)
ax2.scatter((maxvort_NS25[maxit_vort_NS25-1]+maxvort_NS25[maxit_vort_NS25])/2,NS25.y[maxit_vort_NS25]/np.mean(NS25.z_enc),100,color='k',marker='*')
ax2.axhline(np.mean(NS25.z_ig/NS25.z_enc),0,0.05,color='C1',linewidth=2)
ax2.axhline(np.mean(NS25.z_if/NS25.z_enc),0,0.05,color='C1',linewidth=2)
ax1.set_xlabel(r'$\log_{10}(\omega^2/\omega_0^2)$')
ax2.set_xlabel(r'$\log_{10}(\omega^2/\omega_0^2)$')
ax1.set_ylabel(r'$z/z_\mathrm{enc}$')
ax1.set_title(r'(a) $Re_0=117$',fontsize=20,loc='left')
ax2.set_title(r'(b) $Re_0=25$',fontsize=20,loc='left')
cbar_ax = f.add_axes([0.3,0.1,0.5,0.03])
cbar = f.colorbar(cs1,cax=cbar_ax,orientation='horizontal')
plt.tight_layout(rect=[0,0.1,1,1],h_pad=2)
plt.savefig(opath+'Fig11.pdf',bbox_inches='tight')
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
ax1.axhline(np.mean(NS42.z_if/NS42.z_enc),0,0.05,color='C1',linewidth=2)
cs2 = ax2.contourf(S20_42_vortpdf_x_mean,S20_42_vortpdf_y_mean,S20_vortpdf_interp_runmean,cmap=imola_map,levels=np.linspace(0,0.4,9))
ax2.plot(maxvort_S20[:maxit_vort_S20],S20_42.y[:maxit_vort_S20]/np.mean(S20_42.z_enc),'k',lw=1)
ax2.plot(maxvort_S20[maxit_vort_S20:],S20_42.y[maxit_vort_S20:]/np.mean(S20_42.z_enc),'k',lw=1)
ax2.scatter((maxvort_S20[maxit_vort_S20-1]+maxvort_S20[maxit_vort_S20])/2,S20_42.y[maxit_vort_S20]/np.mean(S20_42.z_enc),100,color='k',marker='*')
ax2.axhline(np.mean(S20_42.z_ig/S20_42.z_enc),0,0.05,color='C1',linewidth=2)
ax2.axhline(np.mean(S20_42.z_if/S20_42.z_enc),0,0.05,color='C1',linewidth=2)
ax1.set_ylabel(r'$z/z_\mathrm{enc}$')
ax1.set_xlabel(r'$\log_{10}(\omega^2/\omega_0^2)$')
ax2.set_xlabel(r'$\log_{10}(\omega^2/\omega_0^2)$')
ax1.set_title(r'(a)$Fr_0=0$',fontsize=20,loc='left')
ax2.set_title(r'(b)$Fr_0=20$',fontsize=20,loc='left')
cbar_ax = f.add_axes([0.3,0.1,0.5,0.03])
cbar = f.colorbar(cs1,cax=cbar_ax,orientation='horizontal')
plt.tight_layout(rect=[0,0.1,1,1],h_pad=2)
plt.savefig(opath+'Fig1.pdf',bbox_inches='tight')
plt.show()

