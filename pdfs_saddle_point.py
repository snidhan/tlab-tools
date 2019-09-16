from ReadStats import Statistics, Pdfs 
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.style
from matplotlib import rc


rc('text', usetex=True)
rc('text.latex', preamble=r"\usepackage{fourier}")
rc('font', family='serif')
rc('font', size=20)

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

S20 = Statistics(path_S20+'stats/pdftimes/avg42000-84000.nc')

############################################################################

# Pdf 

vortlist_1 = [path_1+'stats/pdfs/pdf20500.LnEnstrophyW_iW_i',path_1+'stats/pdfs/pdf24000.LnEnstrophyW_iW_i',path_1+'stats/pdfs/pdf27500.LnEnstrophyW_iW_i',path_1+'stats/pdfs/pdf32000.LnEnstrophyW_iW_i',path_1+'stats/pdfs/pdf36500.LnEnstrophyW_iW_i',path_1+'stats/pdfs/pdf41500.LnEnstrophyW_iW_i',path_1+'stats/pdfs/pdf47500.LnEnstrophyW_iW_i',path_1+'stats/pdfs/pdf53000.LnEnstrophyW_iW_i']
vortlist_2 = [path_2+'stats/pdfs/pdf60000.LnEnstrophyW_iW_i',path_2+'stats/pdfs/pdf67000.LnEnstrophyW_iW_i',path_2+'stats/pdfs/pdf74500.LnEnstrophyW_iW_i']
vortlist_3 = [path_3+'stats/pdfs/pdf83000.LnEnstrophyW_iW_i',path_3+'stats/pdfs/pdf93000.LnEnstrophyW_iW_i',path_3+'stats/pdfs/pdf102500.LnEnstrophyW_iW_i',path_3+'stats/pdfs/pdf113000.LnEnstrophyW_iW_i',path_3+'stats/pdfs/pdf127500.LnEnstrophyW_iW_i']

pvlist_1 = [path_1+'stats/pdfs/pdf20500.LnPotentialEnstrophy',path_1+'stats/pdfs/pdf24000.LnPotentialEnstrophy',path_1+'stats/pdfs/pdf27500.LnPotentialEnstrophy',path_1+'stats/pdfs/pdf32000.LnPotentialEnstrophy',path_1+'stats/pdfs/pdf36500.LnPotentialEnstrophy',path_1+'stats/pdfs/pdf41500.LnPotentialEnstrophy',path_1+'stats/pdfs/pdf47500.LnPotentialEnstrophy',path_1+'stats/pdfs/pdf53000.LnPotentialEnstrophy']
pvlist_2 = [path_2+'stats/pdfs/pdf60000.LnPotentialEnstrophy',path_2+'stats/pdfs/pdf67000.LnPotentialEnstrophy',path_2+'stats/pdfs/pdf74500.LnPotentialEnstrophy']
pvlist_3 = [path_3+'stats/pdfs/pdf83000.LnPotentialEnstrophy',path_3+'stats/pdfs/pdf93000.LnPotentialEnstrophy',path_3+'stats/pdfs/pdf102500.LnPotentialEnstrophy',path_3+'stats/pdfs/pdf113000.LnPotentialEnstrophy',path_3+'stats/pdfs/pdf127500.LnPotentialEnstrophy']

vortlist_S20 = [path_S20+'stats/pdfs/pdf42000.LnEnstrophyW_iW_i',path_S20+'stats/pdfs/pdf45000.LnEnstrophyW_iW_i',path_S20+'stats/pdfs/pdf51000.LnEnstrophyW_iW_i',path_S20+'stats/pdfs/pdf58000.LnEnstrophyW_iW_i',path_S20+'stats/pdfs/pdf66000.LnEnstrophyW_iW_i',path_S20+'stats/pdfs/pdf75000.LnEnstrophyW_iW_i',path_S20+'stats/pdfs/pdf84000.LnEnstrophyW_iW_i']

pvlist_S20 = [path_S20+'stats/pdfs/pdf42000.LnPotentialEnstrophy',path_S20+'stats/pdfs/pdf45000.LnPotentialEnstrophy',path_S20+'stats/pdfs/pdf51000.LnPotentialEnstrophy',path_S20+'stats/pdfs/pdf58000.LnPotentialEnstrophy',path_S20+'stats/pdfs/pdf66000.LnPotentialEnstrophy',path_S20+'stats/pdfs/pdf75000.LnPotentialEnstrophy',path_S20+'stats/pdfs/pdf84000.LnPotentialEnstrophy']

NS42_vortpdf_1 = Pdfs(vortlist_1,path_1+'y.dat')
NS42_vortpdf_2 = Pdfs(vortlist_2,path_2+'y.dat')
NS42_vortpdf_3 = Pdfs(vortlist_3,path_3+'y.dat')

NS42_pvpdf_1 = Pdfs(pvlist_1,path_1+'y.dat')
NS42_pvpdf_2 = Pdfs(pvlist_2,path_2+'y.dat')
NS42_pvpdf_3 = Pdfs(pvlist_3,path_3+'y.dat')

vortpdf_S20 = Pdfs(vortlist_S20,path_S20+'y.dat')
pvpdf_S20 = Pdfs(pvlist_S20,path_S20+'y.dat')

# Running average of pdfs

NS42_vortpdf_timmean_1 = np.zeros((len(vortlist_1)-2,NS42_vortpdf_1.ny+1,NS42_vortpdf_1.nb+2))
for t in range(len(vortlist_1)-2):
    NS42_vortpdf_timmean_1[t,:,:] = np.mean(NS42_vortpdf_1.pdf_interp[t:t+3,:,:],axis=0)

# Find where pdf has a maximum at each height

maxvort_1 = np.zeros((len(vortlist_1)-2,NS42_vortpdf_1.ny))
maxprob_vort_1 = np.zeros((len(vortlist_1)-2,NS42_vortpdf_1.ny))
for t in range(len(vortlist_1)-2):
    for i in range(0,NS42_vortpdf_1.ny):
        maxvort_1[t,i] = NS42_vortpdf_1.xy[0,0,i,np.argmax(NS42_vortpdf_timmean_1[t,i,:NS42_vortpdf_1.nb])]
        maxprob_vort_1[t,i] = np.max(NS42_vortpdf_timmean_1[t,i,:NS42_vortpdf_1.nb])
maxvort_1 = np.log10(np.exp(maxvort_1)/(ceps*B0/nu))

maxvort_2 = np.zeros((len(vortlist_2),NS42_vortpdf_2.ny))
maxprob_vort_2 = np.zeros((len(vortlist_2),NS42_vortpdf_2.ny))
for t in range(len(vortlist_2)):
    for i in range(0,NS42_vortpdf_2.ny):
        maxvort_2[t,i] = NS42_vortpdf_2.xy[0,0,i,np.argmax(NS42_vortpdf_2.pdf_interp[t,i,:NS42_vortpdf_2.nb])]
        maxprob_vort_2[t,i] = np.max(NS42_vortpdf_2.pdf_interp[t,i,:NS42_vortpdf_2.nb])
maxvort_2 = np.log10(np.exp(maxvort_2)/(ceps*B0/nu))

maxvort_3 = np.zeros((len(vortlist_3),NS42_vortpdf_3.ny))
maxprob_vort_3 = np.zeros((len(vortlist_3),NS42_vortpdf_3.ny))
for t in range(len(vortlist_3)):
    for i in range(0,NS42_vortpdf_3.ny):
        maxvort_3[t,i] = NS42_vortpdf_3.xy[0,0,i,np.argmax(NS42_vortpdf_3.pdf_interp[t,i,:NS42_vortpdf_3.nb])]
        maxprob_vort_3[t,i] = np.max(NS42_vortpdf_3.pdf_interp[t,i,:NS42_vortpdf_3.nb])
maxvort_3 = np.log10(np.exp(maxvort_3)/(ceps*B0/nu))

maxvort_S20 = np.zeros((len(vortlist_S20),vortpdf_S20.ny))
maxprob_vort_S20 = np.zeros((len(vortlist_S20),vortpdf_S20.ny))
for t in range(len(vortlist_S20)):
    for i in range(0,vortpdf_S20.ny):
        maxvort_S20[t,i] = vortpdf_S20.xy[0,0,i,np.argmax(vortpdf_S20.pdf_interp[t,i,:vortpdf_S20.nb])]
        maxprob_vort_S20[t,i] = np.max(vortpdf_S20.pdf_interp[t,i,:vortpdf_S20.nb])
maxvort_S20 = np.log10(np.exp(maxvort_S20)/(ceps*B0/nu))
         
maxpv_1 = np.zeros((len(pvlist_1),NS42_pvpdf_1.ny))
maxprob_pv_1 = np.zeros((len(pvlist_1),NS42_pvpdf_1.ny))
for t in range(len(pvlist_1)):
    for i in range(0,NS42_pvpdf_1.ny):
        maxpv_1[t,i] = NS42_pvpdf_1.xy[0,0,i,np.argmax(NS42_pvpdf_1.pdf_interp[t,i,:NS42_pvpdf_1.nb])]
        maxprob_pv_1[t,i] = np.max(NS42_pvpdf_1.pdf_interp[t,i,:NS42_pvpdf_1.nb])
    maxpv_1[t,:] = np.log10(np.exp(maxpv_1[t,:])/(cb*ceps*(NS42_1.z_enc[t]/L_0)**(-4./3.)*(B0/nu)**3/42))
    
maxpv_2 = np.zeros((len(pvlist_2),NS42_pvpdf_2.ny))
maxprob_pv_2 = np.zeros((len(pvlist_2),NS42_pvpdf_2.ny))
for t in range(len(pvlist_2)):
    for i in range(0,NS42_pvpdf_2.ny):
        maxpv_2[t,i] = NS42_pvpdf_2.xy[0,0,i,np.argmax(NS42_pvpdf_2.pdf_interp[t,i,:NS42_pvpdf_2.nb])]
        maxprob_pv_2[t,i] = np.max(NS42_pvpdf_2.pdf_interp[t,i,:NS42_pvpdf_2.nb])
    maxpv_2[t,:] = np.log10(np.exp(maxpv_2[t,:])/(cb*ceps*(NS42_2.z_enc[t]/L_0)**(-4./3.)*(B0/nu)**3/42))
    
maxpv_3 = np.zeros((len(pvlist_3),NS42_pvpdf_3.ny))
maxprob_pv_3 = np.zeros((len(pvlist_3),NS42_pvpdf_3.ny))
for t in range(len(pvlist_3)):
    for i in range(0,NS42_pvpdf_3.ny):
        maxpv_3[t,i] = NS42_pvpdf_3.xy[0,0,i,np.argmax(NS42_pvpdf_3.pdf_interp[t,i,:NS42_pvpdf_3.nb])]
        maxprob_pv_3[t,i] = np.max(NS42_pvpdf_3.pdf_interp[t,i,:NS42_pvpdf_3.nb])
    maxpv_3[t,:] = np.log10(np.exp(maxpv_3[t,:])/(cb*ceps*(NS42_3.z_enc[t]/L_0)**(-4./3.)*(B0/nu)**3/42))

maxpv_S20 = np.zeros((len(pvlist_S20),pvpdf_S20.ny))
maxprob_pv_S20 = np.zeros((len(pvlist_S20),pvpdf_S20.ny))
for t in range(len(pvlist_S20)):
    for i in range(0,pvpdf_S20.ny):
        maxpv_S20[t,i] = pvpdf_S20.xy[0,0,i,np.argmax(pvpdf_S20.pdf_interp[t,i,:pvpdf_S20.nb])]
        maxprob_pv_S20[t,i] = np.max(pvpdf_S20.pdf_interp[t,i,:pvpdf_S20.nb])
    maxpv_S20[t,:] = np.log10(np.exp(maxpv_S20[t,:])/(cb*ceps*(S20.z_enc[t]/L_0)**(-4./3.)*(B0/nu)**3/42))

# Find saddle as point where maxprob has a minimum

y_vort_1_saddle = np.zeros(len(vortlist_1)-2)
maxvort_1_saddle = np.zeros(len(vortlist_1)-2)
for t in range(len(vortlist_1)-2):
    y_vort_1_saddle[t] = NS42_1.y[np.argmin(maxprob_vort_1[t,:])] 
    maxvort_1_saddle[t] = maxvort_1[t,np.argmin(np.abs(y_vort_1_saddle[t]-NS42_1.y))]
    y_vort_1_saddle[t] = y_vort_1_saddle[t]/NS42_1.z_enc[t]
    
y_vort_2_saddle = np.zeros(len(vortlist_2))
maxvort_2_saddle = np.zeros(len(vortlist_2))
for t in range(len(vortlist_2)):
    y_vort_2_saddle[t] = NS42_2.y[np.argmin(maxprob_vort_2[t,NS42_2.z_enc_arg[t]:NS42_2.z_ig_arg[t]])+NS42_2.z_enc_arg[t]]
    maxvort_2_saddle[t] = maxvort_2[t,np.argmin(np.abs(y_vort_2_saddle[t]-NS42_2.y))]
    y_vort_2_saddle[t] = y_vort_2_saddle[t]/NS42_2.z_enc[t]
     
y_vort_3_saddle = np.zeros(len(vortlist_3))
maxvort_3_saddle = np.zeros(len(vortlist_3))
for t in range(len(vortlist_3)):
    y_vort_3_saddle[t] = NS42_3.y[np.argmin(maxprob_vort_3[t,NS42_3.z_enc_arg[t]:NS42_3.z_ig_arg[t]])+NS42_3.z_enc_arg[t]]
    maxvort_3_saddle[t] = maxvort_3[t,np.argmin(np.abs(y_vort_3_saddle[t]-NS42_3.y))]
    y_vort_3_saddle[t] = y_vort_3_saddle[t]/NS42_3.z_enc[t]

y_vort_S20_saddle = np.zeros(len(vortlist_S20))
maxvort_S20_saddle = np.zeros(len(vortlist_S20))
for t in range(len(vortlist_S20)):
    y_vort_S20_saddle[t] = S20.y[np.argmin(maxprob_vort_S20[t,S20.z_enc_arg[t]:S20.z_ig_arg[t]])+S20.z_enc_arg[t]] 
    maxvort_S20_saddle[t] = maxvort_S20[t,np.argmin(np.abs(y_vort_S20_saddle[t]-S20.y))]
    y_vort_S20_saddle[t] = y_vort_S20_saddle[t]/S20.z_enc[t]

y_pv_1_saddle = np.zeros(len(pvlist_1))
maxpv_1_saddle = np.zeros(len(pvlist_1))
for t in range(len(pvlist_1)):
    y_pv_1_saddle[t] = NS42_1.y[np.argmin(maxprob_pv_1[t,NS42_1.z_enc_arg[t]:-30])+NS42_1.z_enc_arg[t]]
    maxpv_1_saddle[t] = maxpv_1[t,np.argmin(np.abs(y_pv_1_saddle[t]-NS42_1.y))]
    y_pv_1_saddle[t] = y_pv_1_saddle[t]/NS42_1.z_enc[t]
    
y_pv_2_saddle = np.zeros(len(pvlist_2))
maxpv_2_saddle = np.zeros(len(pvlist_2))
for t in range(len(pvlist_2)):
    y_pv_2_saddle[t] = NS42_2.y[np.argmin(maxprob_pv_2[t,NS42_2.z_enc_arg[t]:-30])+NS42_2.z_enc_arg[t]]
    maxpv_2_saddle[t] = maxpv_2[t,np.argmin(np.abs(y_pv_2_saddle[t]-NS42_2.y))]
    y_pv_2_saddle[t] = y_pv_2_saddle[t]/NS42_2.z_enc[t]
      
y_pv_3_saddle = np.zeros(len(pvlist_3))
maxpv_3_saddle = np.zeros(len(pvlist_3))
for t in range(len(pvlist_3)):
    y_pv_3_saddle[t] = NS42_3.y[np.argmin(maxprob_pv_3[t,NS42_3.z_enc_arg[t]:-30])+NS42_3.z_enc_arg[t]]
    maxpv_3_saddle[t] = maxpv_3[t,np.argmin(np.abs(y_pv_3_saddle[t]-NS42_3.y))]
    y_pv_3_saddle[t] = y_pv_3_saddle[t]/NS42_3.z_enc[t]

y_pv_S20_saddle = np.zeros(len(pvlist_S20))
maxpv_S20_saddle = np.zeros(len(pvlist_S20))
for t in range(len(pvlist_S20)):
    y_pv_S20_saddle[t] = S20.y[np.argmin(maxprob_pv_S20[t,S20.z_enc_arg[t]:-30])+S20.z_enc_arg[t]]
    maxpv_S20_saddle[t] = maxpv_S20[t,np.argmin(np.abs(y_pv_S20_saddle[t]-S20.y))]
    y_pv_S20_saddle[t] = y_pv_S20_saddle[t]/S20.z_enc[t]

# concatenate over time

time = np.concatenate((NS42_1.z_enc/L_0,NS42_2.z_enc/L_0,NS42_3.z_enc/L_0))
z_is = np.concatenate((NS42_1.z_is/NS42_1.z_enc,NS42_2.z_is/NS42_2.z_enc,NS42_3.z_is/NS42_3.z_enc))
z_if = np.concatenate((NS42_1.z_if/NS42_1.z_enc,NS42_2.z_if/NS42_2.z_enc,NS42_3.z_if/NS42_3.z_enc))
y_vort_saddle = np.concatenate((y_vort_1_saddle,y_vort_2_saddle,y_vort_3_saddle))
y_pv_saddle = np.concatenate((y_pv_1_saddle,y_pv_2_saddle,y_pv_3_saddle))
maxvort_saddle = np.concatenate((maxvort_1_saddle,maxvort_2_saddle,maxvort_3_saddle))
maxpv_saddle = np.concatenate((maxpv_1_saddle,maxpv_2_saddle,maxpv_3_saddle))
y_vort_saddle_mean = runningmean(y_vort_saddle,1)
y_pv_saddle_mean = runningmean(y_pv_saddle,1)
maxvort_saddle_mean = runningmean(maxvort_saddle,1)
maxpv_saddle_mean = runningmean(maxpv_saddle,1)

#####################################################################
# Plot



f, (ax1,ax2) = plt.subplots(1,2,sharex='all',figsize=(10,5))
ax1.grid(True)
ax2.grid(True)
ax1.set_ylim(1,1.4)
ax1.set_xlim(15,30)
ax2.set_ylim(-4.5,0.5)
ax2.set_xlim(15,30)
ax1.plot(time[1:-1],y_vort_saddle_mean)
ax1.plot(time[1:-1],y_pv_saddle_mean)
ax1.plot(S20.z_enc[1:-1]/L_0,runningmean(y_vort_S20_saddle,1))
ax1.plot(S20.z_enc[1:-1]/L_0,runningmean(y_pv_S20_saddle,1))
ax1.plot(S20.z_enc[1:-1]/L_0,runningmean(S20.z_ig/S20.z_enc,1),'k--',label=r'$z_{i,g}/z_\mathrm{enc}$')
ax1.plot(time[1:-1],runningmean(z_is,1),'k-.',label=r'$z_{i,s}/z_\mathrm{enc}$')
ax1.plot(time[1:-1],runningmean(z_if,1),'k:',label=r'$z_{i,f}/z_\mathrm{enc}$')
ax2.plot(time[1:-1],maxvort_saddle_mean,label=r'$Fr_0=0$')
ax2.plot(time[1:-1],maxpv_saddle_mean,label=r'$Fr_0=0$')
ax2.plot(S20.z_enc[1:-1]/L_0,runningmean(maxvort_S20_saddle,1),label=r'$Fr_0=20$')
ax2.plot(S20.z_enc[1:-1]/L_0,runningmean(maxpv_S20_saddle,1),label=r'$Fr_0=20$')
ax1.text(16,1.05,r'$Fr_0=0$',color='C0',fontsize=20)
ax1.text(16,1.02,r'$Fr_0=0$',color='C1',fontsize=20)
ax1.text(23,1.05,r'$Fr_0=20$',color='C2',fontsize=20)
ax1.text(23,1.02,r'$Fr_0=20$',color='C3',fontsize=20)
ax2.text(16,-3.5,r'$\log_{10}(\omega^2/\omega_0^2)$',color='C0',fontsize=20)
ax2.text(16,-4.1,r'$\log_{10}(\Pi^2/\Pi_0^2)$',color='C1',fontsize=20)
ax2.text(24,-3.5,r'$\log_{10}(\omega^2/\omega_0^2)$',color='C2',fontsize=20)
ax2.text(24,-4.1,r'$\log_{10}(\Pi^2/\Pi_0^2)$',color='C3',fontsize=20)
ax1.set_title('(a)',fontsize=20,loc='left')
ax2.set_title('(b)',fontsize=20,loc='left')
ax1.set_xlabel(r'$z_\mathrm{enc}/L_0$')
ax2.set_xlabel(r'$z_\mathrm{enc}/L_0$')
ax1.set_ylabel(r'$z_\mathrm{saddle}/z_\mathrm{enc}$')
ax1.legend(loc='best',fontsize=20)
plt.tight_layout()
#plt.savefig(opath+'pdfs_saddle_height_time_S20_S0.pdf')
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
