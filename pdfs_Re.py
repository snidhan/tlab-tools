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
it=1

##########################################################################
# Stats

NS117 = Statistics(path_117+'5120x1024x5120/stats/avg150000-150000.nc')
NS42 = Statistics(path_42+'2560x576x2560/stats/pdftimes/avg36500-47500.nc')
NS25 = Statistics(path_25+'2560x512x2560/stats/avg21000-21000.nc')

S20_42 = Statistics(path_42+'3072x960x4608-S20/stats/pdftimes/avg66000-84000.nc')

############################################################################
# Pdf

NS117_vortpdf = Pdfs([path_117+'5120x1024x5120/stats/pdfs/pdf150000.LnEnstrophyW_iW_i'],path_117+'5120x1024x5120/y.dat')
NS117_pvpdf   = Pdfs([path_117+'5120x1024x5120/stats/pdfs/pdf150000.LnPotentialEnstrophy'],path_117+'5120x1024x5120/y.dat')
NS117_s1gradpdf = Pdfs([path_117+'5120x1024x5120/stats/pdfs/pdf150000.LnGradientG_iG_i'],path_117+'5120x1024x5120/y.dat')
NS42_vortpdf  = Pdfs([path_42+'2560x576x2560/stats/pdfs/pdf36500.LnEnstrophyW_iW_i',path_42+'2560x576x2560/stats/pdfs/pdf41500.LnEnstrophyW_iW_i',path_42+'2560x576x2560/stats/pdfs/pdf47500.LnEnstrophyW_iW_i'],path_42+'2560x576x2560/y.dat')
NS42_pvpdf    = Pdfs([path_42+'2560x576x2560/stats/pdfs/pdf36500.LnPotentialEnstrophy',path_42+'2560x576x2560/stats/pdfs/pdf41500.LnPotentialEnstrophy',path_42+'2560x576x2560/stats/pdfs/pdf47500.LnPotentialEnstrophy'],path_42+'2560x576x2560/y.dat')
NS42_s1gradpdf = Pdfs([path_42+'2560x576x2560/stats/pdfs/pdf36500.LnGradientG_iG_i',path_42+'2560x576x2560/stats/pdfs/pdf41500.LnGradientG_iG_i',path_42+'2560x576x2560/stats/pdfs/pdf47500.LnGradientG_iG_i'],path_42+'2560x576x2560/y.dat')
NS25_vortpdf = Pdfs([path_25+'2560x512x2560/stats/pdfs/pdf21000.LnEnstrophyW_iW_i'],path_25+'2560x512x2560/y.dat')
NS25_pvpdf = Pdfs([path_25+'2560x512x2560/stats/pdfs/pdf21000.LnPotentialEnstrophy'],path_25+'2560x512x2560/y.dat')

S20_42_vortpdf = Pdfs([path_42+'3072x960x4608-S20/stats/pdfs/pdf66000.LnEnstrophyW_iW_i',path_42+'3072x960x4608-S20/stats/pdfs/pdf75000.LnEnstrophyW_iW_i',path_42+'3072x960x4608-S20/stats/pdfs/pdf84000.LnEnstrophyW_iW_i'],path_42+'3072x960x4608-S20/y.dat')
S20_42_pvpdf = Pdfs([path_42+'3072x960x4608-S20/stats/pdfs/pdf66000.LnPotentialEnstrophy',path_42+'3072x960x4608-S20/stats/pdfs/pdf75000.LnPotentialEnstrophy',path_42+'3072x960x4608-S20/stats/pdfs/pdf84000.LnPotentialEnstrophy'],path_42+'3072x960x4608-S20/y.dat')
S20_42_s1gradpdf = Pdfs([path_42+'3072x960x4608-S20/stats/pdfs/pdf66000.LnGradientG_iG_i',path_42+'3072x960x4608-S20/stats/pdfs/pdf75000.LnGradientG_iG_i',path_42+'3072x960x4608-S20/stats/pdfs/pdf84000.LnGradientG_iG_i'],path_42+'3072x960x4608-S20/y.dat')

# Find where pdf has a maximum at each height
maxvort_117 = np.zeros(NS117_vortpdf.ny)
maxprob_vort_117 = np.zeros(NS117_vortpdf.ny)
for i in range(0,NS117_vortpdf.ny):
    maxvort_117[i] = NS117_vortpdf.xy[0,i,np.argmax(NS117_vortpdf.pdf_timeavg[i,:NS117_vortpdf.nb])]
    maxprob_vort_117[i] = np.max(NS117_vortpdf.pdf_timeavg[i,:NS117_vortpdf.nb])
maxvort_117 = np.log10(np.exp(maxvort_117)/(ceps*B0/nu_117)) 

maxpv_117 = np.zeros(NS117_pvpdf.ny)
maxprob_pv_117 = np.zeros(NS117_pvpdf.ny)
for i in range(0,NS117_pvpdf.ny):
    maxpv_117[i] = NS117_pvpdf.xy[0,i,np.argmax(NS117_pvpdf.pdf_timeavg[i,:NS117_pvpdf.nb])]
    maxprob_pv_117[i] = np.max(NS117_pvpdf.pdf_timeavg[i,:NS117_pvpdf.nb])
maxpv_117 = np.log10(np.exp(maxpv_117)/(cb*ceps*(NS117.z_enc/L0)**(-4./3.)*(B0/nu_117)**3/117))
    
maxvort_42_avg = np.zeros(NS42_vortpdf.ny)
maxprob_vort_42_avg = np.zeros(NS42_vortpdf.ny)
for i in range(0,NS42_vortpdf.ny):
    maxvort_42_avg[i] = NS42_vortpdf.xy[0,i,np.argmax(NS42_vortpdf.pdf_timeavg[i,:NS42_vortpdf.nb])]
    maxprob_vort_42_avg[i] = np.max(NS42_vortpdf.pdf_timeavg[i,:NS42_vortpdf.nb])
maxvort_42_avg = np.log10(np.exp(maxvort_42_avg)/(ceps*B0/nu_42))
    
maxpv_42_avg = np.zeros(NS42_pvpdf.ny)
maxprob_pv_42_avg = np.zeros(NS42_pvpdf.ny)
for i in range(0,NS42_pvpdf.ny):
    maxpv_42_avg[i] = NS42_pvpdf.xy[0,i,np.argmax(NS42_pvpdf.pdf_timeavg[i,:NS42_pvpdf.nb])]
    maxprob_pv_42_avg[i] = np.max(NS42_pvpdf.pdf_timeavg[i,:NS42_pvpdf.nb])
maxpv_42_avg = np.log10(np.exp(maxpv_42_avg)/(cb*ceps*(np.mean(NS42.z_enc)/L0)**(-4./3.)*(B0/nu_42)**3/42))

maxvort_42 = np.zeros(NS42_vortpdf.ny)
maxprob_vort_42 = np.zeros(NS42_vortpdf.ny)
for i in range(0,NS42_vortpdf.ny):
    maxvort_42[i] = NS42_vortpdf.xy[0,i,np.argmax(NS42_vortpdf.pdf[it,i,:NS42_vortpdf.nb])]
    maxprob_vort_42[i] = np.max(NS42_vortpdf.pdf[it,i,:NS42_vortpdf.nb])
maxvort_42 = np.log10(np.exp(maxvort_42)/(ceps*B0/nu_42))
    
maxpv_42 = np.zeros(NS42_pvpdf.ny)
maxprob_pv_42 = np.zeros(NS42_pvpdf.ny)
for i in range(0,NS42_pvpdf.ny):
    maxpv_42[i] = NS42_pvpdf.xy[0,i,np.argmax(NS42_pvpdf.pdf[it,i,:NS42_pvpdf.nb])]
    maxprob_pv_42[i] = np.max(NS42_pvpdf.pdf[it,i,:NS42_pvpdf.nb])
maxpv_42 = np.log10(np.exp(maxpv_42)/(cb*ceps*(NS42.z_enc[it]/L0)**(-4./3.)*(B0/nu_42)**3/42))
    
maxvort_25 = np.zeros(NS25_vortpdf.ny)
maxprob_vort_25 = np.zeros(NS25_vortpdf.ny)
for i in range(0,NS25_vortpdf.ny):
    maxvort_25[i] = NS25_vortpdf.xy[0,i,np.argmax(NS25_vortpdf.pdf_timeavg[i,:NS25_vortpdf.nb])]
    maxprob_vort_25[i] = np.max(NS25_vortpdf.pdf_timeavg[i,:NS25_vortpdf.nb])
maxvort_25 = np.log10(np.exp(maxvort_25)/(ceps*B0/nu_25))
    
maxpv_25 = np.zeros(NS25_pvpdf.ny)
maxprob_pv_25 = np.zeros(NS25_pvpdf.ny)
for i in range(0,NS25_pvpdf.ny):
    maxpv_25[i] = NS25_pvpdf.xy[0,i,np.argmax(NS25_pvpdf.pdf_timeavg[i,:NS25_pvpdf.nb])]
    maxprob_pv_25[i] = np.max(NS25_pvpdf.pdf_timeavg[i,:NS25_pvpdf.nb])
maxpv_25 = np.log10(np.exp(maxpv_25)/(cb*ceps*(NS25.z_enc/L0)**(-4./3.)*(B0/nu_25)**3/25))

maxvort_S20_42_avg = np.zeros(S20_42_vortpdf.ny)
maxprob_vort_S20_42_avg = np.zeros(S20_42_vortpdf.ny)
for i in range(0,S20_42_vortpdf.ny):
    maxvort_S20_42_avg[i] = S20_42_vortpdf.xy[0,i,np.argmax(S20_42_vortpdf.pdf_timeavg[i,:S20_42_vortpdf.nb])]
    maxprob_vort_S20_42_avg[i] = np.max(S20_42_vortpdf.pdf_timeavg[i,:S20_42_vortpdf.nb])
maxvort_S20_42_avg = np.log10(np.exp(maxvort_S20_42_avg)/(ceps*B0/nu_42))

maxpv_S20_42_avg = np.zeros(S20_42_pvpdf.ny)
maxprob_pv_S20_42_avg = np.zeros(S20_42_pvpdf.ny)
for i in range(0,S20_42_pvpdf.ny):
    maxpv_S20_42_avg[i] = S20_42_pvpdf.xy[0,i,np.argmax(S20_42_pvpdf.pdf_timeavg[i,:S20_42_pvpdf.nb])]
    maxprob_pv_S20_42_avg[i] = np.max(S20_42_pvpdf.pdf_timeavg[i,:S20_42_pvpdf.nb])
maxpv_S20_42_avg = np.log10(np.exp(maxpv_S20_42_avg)/(cb*ceps*(np.mean(S20_42.z_enc)/L0)**(-4./3.)*(B0/nu_42)**3/42))

maxvort_S20_42 = np.zeros(S20_42_vortpdf.ny)
maxprob_vort_S20_42 = np.zeros(S20_42_vortpdf.ny)
for i in range(0,S20_42_vortpdf.ny):
    maxvort_S20_42[i] = S20_42_vortpdf.xy[0,i,np.argmax(S20_42_vortpdf.pdf[it,i,:S20_42_vortpdf.nb])]
    maxprob_vort_S20_42[i] = np.max(S20_42_vortpdf.pdf[it,i,:S20_42_vortpdf.nb])
maxvort_S20_42 = np.log10(np.exp(maxvort_S20_42)/(ceps*B0/nu_42))

maxpv_S20_42 = np.zeros(S20_42_pvpdf.ny)
maxprob_pv_S20_42 = np.zeros(S20_42_pvpdf.ny)
for i in range(0,S20_42_pvpdf.ny):
    maxpv_S20_42[i] = S20_42_pvpdf.xy[0,i,np.argmax(S20_42_pvpdf.pdf[it,i,:S20_42_pvpdf.nb])]
    maxprob_pv_S20_42[i] = np.max(S20_42_pvpdf.pdf[it,i,:S20_42_pvpdf.nb])
maxpv_S20_42 = np.log10(np.exp(maxpv_S20_42)/(cb*ceps*(S20_42.z_enc[it]/L0)**(-4./3.)*(B0/nu_42)**3/42))
    
# Find saddle as point where maxprob has a minimum

y_vort_117_saddle = NS117.y[np.argmin(maxprob_vort_117)]
maxvort_117_saddle = maxvort_117[np.argmin(np.abs(y_vort_117_saddle-NS117.y))]
y_vort_117_saddle = y_vort_117_saddle/NS117.z_enc

y_pv_117_saddle = NS117.y[np.argmin(maxprob_pv_117[NS117.z_enc_arg[0]:])+NS117.z_enc_arg[0]]
maxpv_117_saddle = maxpv_117[np.argmin(np.abs(y_pv_117_saddle-NS117.y))]
y_pv_117_saddle = y_pv_117_saddle/NS117.z_enc

y_vort_42_saddle_avg = NS42.y[np.argmin(maxprob_vort_42_avg)]
maxvort_42_saddle_avg = maxvort_42_avg[np.argmin(np.abs(y_vort_42_saddle_avg-NS42.y))]
y_vort_42_saddle_avg = y_vort_42_saddle_avg/np.mean(NS42.z_enc)

y_pv_42_saddle_avg = NS42.y[np.argmin(maxprob_pv_42_avg[NS42.z_enc_arg[0]:])+NS42.z_enc_arg[0]]
maxpv_42_saddle_avg = maxpv_42_avg[np.argmin(np.abs(y_pv_42_saddle_avg-NS42.y))]
y_pv_42_saddle_avg = y_pv_42_saddle_avg/np.mean(NS42.z_enc)

y_vort_42_saddle = NS42.y[np.argmin(maxprob_vort_42)]
maxvort_42_saddle = maxvort_42[np.argmin(np.abs(y_vort_42_saddle-NS42.y))]
y_vort_42_saddle = y_vort_42_saddle/NS42.z_enc[it]

y_pv_42_saddle = NS42.y[np.argmin(maxprob_pv_42[NS42.z_enc_arg[0]:])+NS42.z_enc_arg[0]]
maxpv_42_saddle = maxpv_42[np.argmin(np.abs(y_pv_42_saddle-NS42.y))]
y_pv_42_saddle = y_pv_42_saddle/NS42.z_enc[it]

y_vort_25_saddle = NS25.y[np.argmin(maxprob_vort_25)]
maxvort_25_saddle = maxvort_25[np.argmin(np.abs(y_vort_25_saddle-NS25.y))]
y_vort_25_saddle = y_vort_25_saddle/NS25.z_enc

y_pv_25_saddle = NS25.y[np.argmin(maxprob_pv_25[NS25.z_enc_arg[0]:])+NS25.z_enc_arg[0]]
maxpv_25_saddle = maxpv_25[np.argmin(np.abs(y_pv_25_saddle-NS25.y))]
y_pv_25_saddle = y_pv_25_saddle/NS25.z_enc

y_vort_S20_42_saddle_avg = S20_42.y[np.argmin(maxprob_vort_S20_42_avg)]
maxvort_S20_42_saddle_avg = maxvort_S20_42_avg[np.argmin(np.abs(y_vort_S20_42_saddle_avg-S20_42.y))]
y_vort_S20_42_saddle_avg = y_vort_S20_42_saddle_avg/np.mean(S20_42.z_enc)

y_pv_S20_42_saddle_avg = S20_42.y[np.argmin(maxprob_pv_S20_42_avg[S20_42.z_enc_arg[0]:])+S20_42.z_enc_arg[0]]
maxpv_S20_42_saddle_avg = maxpv_S20_42_avg[np.argmin(np.abs(y_pv_S20_42_saddle_avg-S20_42.y))]
y_pv_S20_42_saddle_avg = y_pv_S20_42_saddle_avg/np.mean(S20_42.z_enc)

y_vort_S20_42_saddle = S20_42.y[np.argmin(maxprob_vort_S20_42)]
maxvort_S20_42_saddle = maxvort_S20_42[np.argmin(np.abs(y_vort_S20_42_saddle-S20_42.y))]
y_vort_S20_42_saddle = y_vort_S20_42_saddle/S20_42.z_enc[it]

y_pv_S20_42_saddle = S20_42.y[np.argmin(maxprob_pv_S20_42[S20_42.z_enc_arg[0]:])+S20_42.z_enc_arg[0]]
maxpv_S20_42_saddle = maxpv_S20_42[np.argmin(np.abs(y_pv_S20_42_saddle-S20_42.y))]
y_pv_S20_42_saddle = y_pv_S20_42_saddle/S20_42.z_enc[it]

# Normalisation of y axis
NS117_vortpdf.xy[1,:,:] = NS117_vortpdf.xy[1,:,:]/NS117.z_enc
NS117_pvpdf.xy[1,:,:] = NS117_pvpdf.xy[1,:,:]/NS117.z_enc
NS117_s1gradpdf.xy[1,:,:] = NS117_s1gradpdf.xy[1,:,:]/NS117.z_enc

NS42_vortpdf_y_mean = NS42_vortpdf.xy[1,:,:]/np.mean(NS42.z_enc)
NS42_pvpdf_y_mean = NS42_pvpdf.xy[1,:,:]/np.mean(NS42.z_enc)
NS42_s1gradpdf_y_mean = NS42_s1gradpdf.xy[1,:,:]/np.mean(NS42.z_enc)

NS42_vortpdf.xy[1,:,:] = NS42_vortpdf.xy[1,:,:]/NS42.z_enc[it]
NS42_pvpdf.xy[1,:,:] = NS42_pvpdf.xy[1,:,:]/NS42.z_enc[it]
NS42_s1gradpdf.xy[1,:,:] = NS42_s1gradpdf.xy[1,:,:]/NS42.z_enc[it]

NS25_vortpdf.xy[1,:,:] = NS25_vortpdf.xy[1,:,:]/NS25.z_enc
NS25_pvpdf.xy[1,:,:] = NS25_pvpdf.xy[1,:,:]/NS25.z_enc

S20_42_vortpdf_y_mean = S20_42_vortpdf.xy[1,:,:]/np.mean(S20_42.z_enc)
S20_42_pvpdf_y_mean = S20_42_pvpdf.xy[1,:,:]/np.mean(S20_42.z_enc)
S20_42_s1gradpdf_y_mean = S20_42_s1gradpdf.xy[1,:,:]/np.mean(S20_42.z_enc)

S20_42_vortpdf.xy[1,:,:] = S20_42_vortpdf.xy[1,:,:]/S20_42.z_enc[it]
S20_42_pvpdf.xy[1,:,:] = S20_42_pvpdf.xy[1,:,:]/S20_42.z_enc[it]
S20_42_s1gradpdf.xy[1,:,:] = S20_42_s1gradpdf.xy[1,:,:]/S20_42.z_enc[it]

# Normalisation of x axis
NS117_vortpdf.xy[0,:,:] = np.log10(np.exp(NS117_vortpdf.xy[0,:,:])/(ceps*B0/nu_117)) 
NS117_pvpdf.xy[0,:,:] = np.log10(np.exp(NS117_pvpdf.xy[0,:,:])/(cb*ceps*(NS117.z_enc/L0)**(-4./3.)*(B0/nu_117)**3/117))
NS117_s1gradpdf.xy[0,:,:] = np.log10(np.exp(NS117_s1gradpdf.xy[0,:,:])/(cb*N**4*117*(NS117.z_enc/L0)**(-4./3.)))

NS42_vortpdf_x_mean = np.log10(np.exp(NS42_vortpdf.xy[0,:,:])/(ceps*B0/nu_42))
NS42_pvpdf_x_mean = np.log10(np.exp(NS42_pvpdf.xy[0,:,:])/(cb*ceps*(np.mean(NS42.z_enc)/L0)**(-4./3.)*(B0/nu_42)**3/42))
NS42_s1gradpdf_x_mean = np.log10(np.exp(NS42_s1gradpdf.xy[0,:,:])/(cb*N**4*42*(np.mean(NS42.z_enc)/L0)**(-4./3.)))

NS42_vortpdf.xy[0,:,:] = np.log10(np.exp(NS42_vortpdf.xy[0,:,:])/(ceps*B0/nu_42))
NS42_pvpdf.xy[0,:,:] = np.log10(np.exp(NS42_pvpdf.xy[0,:,:])/(cb*ceps*(NS42.z_enc[it]/L0)**(-4./3.)*(B0/nu_42)**3/42))
NS42_s1gradpdf.xy[0,:,:] = np.log10(np.exp(NS42_s1gradpdf.xy[0,:,:])/(cb*N**4*42*(NS42.z_enc[it]/L0)**(-4./3.)))

NS25_vortpdf.xy[0,:,:] = np.log10(np.exp(NS25_vortpdf.xy[0,:,:])/(ceps*B0/nu_25))
NS25_pvpdf.xy[0,:,:] = np.log10(np.exp(NS25_pvpdf.xy[0,:,:])/(cb*ceps*(NS25.z_enc/L0)**(-4./3.)*(B0/nu_25)**3/25))

S20_42_vortpdf_x_mean = np.log10(np.exp(S20_42_vortpdf.xy[0,:,:])/(ceps*B0/nu_42))
S20_42_pvpdf_x_mean = np.log10(np.exp(S20_42_pvpdf.xy[0,:,:])/(cb*ceps*(np.mean(S20_42.z_enc)/L0)**(-4./3.)*(B0/nu_42)**3/42))
S20_42_s1gradpdf_x_mean = np.log10(np.exp(S20_42_s1gradpdf.xy[0,:,:])/(cb*N**4*42*(np.mean(S20_42.z_enc)/L0)**(-4./3.)))

S20_42_vortpdf.xy[0,:,:] = np.log10(np.exp(S20_42_vortpdf.xy[0,:,:])/(ceps*B0/nu_42))
S20_42_pvpdf.xy[0,:,:] = np.log10(np.exp(S20_42_pvpdf.xy[0,:,:])/(cb*ceps*(S20_42.z_enc[it]/L0)**(-4./3.)*(B0/nu_42)**3/42))
S20_42_s1gradpdf.xy[0,:,:] = np.log10(np.exp(S20_42_s1gradpdf.xy[0,:,:])/(cb*N**4*42*(S20_42.z_enc[it]/L0)**(-4./3.)))

#####################################################################
# Colourmaps

imola_data = np.loadtxt(colourmap_path+'imola/imola.txt')
imola_map = LinearSegmentedColormap.from_list('imola',imola_data)

davos_data = np.loadtxt(colourmap_path+'davos/davos.txt')
davos_map = LinearSegmentedColormap.from_list('davos',davos_data)
#####################################################################
# Plot

f, ((ax1,ax2),(ax3,ax4),(ax5,ax6)) = plt.subplots(3,2,sharex='col',sharey='all',figsize=(10,16))
ax1.tick_params(bottom=True,top=True,left=True,right=True)
ax2.tick_params(bottom=True,top=True,left=True,right=True)
ax3.tick_params(bottom=True,top=True,left=True,right=True)
ax4.tick_params(bottom=True,top=True,left=True,right=True)
ax5.tick_params(bottom=True,top=True,left=True,right=True)
ax6.tick_params(bottom=True,top=True,left=True,right=True)
ax1.set_xlim(-4,2) 
ax2.set_xlim(-7,2) 
ax1.set_ylim(0,1.6)
ax2.set_xticks([-6,-4,-2,0,2])
cs1 = ax1.contourf(NS117_vortpdf.xy[0,0,:,:],NS117_vortpdf.xy[1,0,:,:],NS117_vortpdf.pdf_timeavg[:NS117_vortpdf.ny,:NS117_vortpdf.nb],cmap=imola_map,extend='max')
ax1.plot(maxvort_117,NS117.y/NS117.z_enc,'k',lw=1)
ax1.scatter(maxvort_117_saddle,y_vort_117_saddle,100,color='k',marker='*')
ax1.axhline(NS117.z_ig/NS117.z_enc,0,0.05,color='k',linewidth=2)
ax1.axhline(NS117.z_is/NS117.z_enc,0,0.05,color='k',linewidth=2)
ax1.axhline(NS117.z_if/NS117.z_enc,0,0.05,color='k',linewidth=2)
cs2 = ax2.contourf(NS117_pvpdf.xy[0,0,:,:],NS117_pvpdf.xy[1,0,:,:],NS117_pvpdf.pdf_timeavg[:NS117_pvpdf.ny,:NS117_pvpdf.nb],cmap=imola_map,extend='max')
ax2.plot(maxpv_117,NS117.y/NS117.z_enc,'k',lw=1)
ax2.scatter(maxpv_117_saddle,y_pv_117_saddle,100,color='k',marker='*')
ax2.axhline(NS117.z_ig/NS117.z_enc,0,0.05,color='k',linewidth=2)
ax2.axhline(NS117.z_is/NS117.z_enc,0,0.05,color='k',linewidth=2)
ax2.axhline(NS117.z_if/NS117.z_enc,0,0.05,color='k',linewidth=2)
cs3 = ax3.contourf(NS42_vortpdf.xy[0,0,:,:],NS42_vortpdf.xy[1,0,:,:],NS42_vortpdf.pdf_timeavg[:NS42_vortpdf.ny,:NS42_vortpdf.nb],cmap=imola_map,extend='max')
ax3.plot(maxvort_42,NS42.y/NS42.z_enc,'k',lw=1)
ax3.scatter(maxvort_42_saddle,y_vort_42_saddle,100,color='k',marker='*')
ax3.axhline(NS42.z_ig/NS42.z_enc,0,0.05,color='k',linewidth=2)
ax3.axhline(NS42.z_is/NS42.z_enc,0,0.05,color='k',linewidth=2)
ax3.axhline(NS42.z_if/NS42.z_enc,0,0.05,color='k',linewidth=2)
cs4 = ax4.contourf(NS42_pvpdf.xy[0,0,:,:],NS42_pvpdf.xy[1,0,:,:],NS42_pvpdf.pdf_timeavg[:NS42_pvpdf.ny,:NS42_pvpdf.nb],cmap=imola_map,levels=np.linspace(0,0.24,9),extend='max')
ax4.plot(maxpv_42,NS42.y/NS42.z_enc,'k',lw=1)
ax4.scatter(maxpv_42_saddle,y_pv_42_saddle,100,color='k',marker='*')
ax4.axhline(NS42.z_ig/NS42.z_enc,0,0.05,color='k',linewidth=2)
ax4.axhline(NS42.z_is/NS42.z_enc,0,0.05,color='k',linewidth=2)
ax4.axhline(NS42.z_if/NS42.z_enc,0,0.05,color='k',linewidth=2)
cs5 = ax5.contourf(NS25_vortpdf.xy[0,0,:,:],NS25_vortpdf.xy[1,0,:,:],NS25_vortpdf.pdf_timeavg[:NS25_vortpdf.ny,:NS25_vortpdf.nb],cmap=imola_map,extend='max')
ax5.plot(maxvort_25,NS25.y/NS25.z_enc,'k',lw=1)
ax5.scatter(maxvort_25_saddle,y_vort_25_saddle,100,color='k',marker='*')
ax5.axhline(NS25.z_ig/NS25.z_enc,0,0.05,color='k',linewidth=2)
ax5.axhline(NS25.z_is/NS25.z_enc,0,0.05,color='k',linewidth=2)
ax5.axhline(NS25.z_if/NS25.z_enc,0,0.05,color='k',linewidth=2)
cs6 = ax6.contourf(NS25_pvpdf.xy[0,0,:,:],NS25_pvpdf.xy[1,0,:,:],NS25_pvpdf.pdf_timeavg[:NS25_pvpdf.ny,:NS25_pvpdf.nb],cmap=imola_map,levels=np.linspace(0,0.24,9),extend='max')
ax6.plot(maxpv_25,NS25.y/NS25.z_enc,'k',lw=1)
ax6.scatter(maxpv_25_saddle,y_pv_25_saddle,100,color='k',marker='*')
ax6.axhline(NS25.z_ig/NS25.z_enc,0,0.05,color='k',linewidth=2)
ax6.axhline(NS25.z_is/NS25.z_enc,0,0.05,color='k',linewidth=2)
ax6.axhline(NS25.z_if/NS25.z_enc,0,0.05,color='k',linewidth=2)
ax5.set_xlabel(r'$\log_{10}(\omega^2/\omega_0^2)$')
ax6.set_xlabel(r'$\log_{10}(\Pi^2/\Pi_0^2)$')
ax1.set_ylabel(r'$z/z_\mathrm{enc}$')
ax3.set_ylabel(r'$z/z_\mathrm{enc}$')
ax5.set_ylabel(r'$z/z_\mathrm{enc}$')
ax1.set_title(r'(a) $Re_0=117$',fontsize=20,loc='left')
ax2.set_title(r'(b) $Re_0=117$',fontsize=20,loc='left')
ax3.set_title(r'(c) $Re_0=42$',fontsize=20,loc='left')
ax4.set_title(r'(d) $Re_0=42$',fontsize=20,loc='left')
ax5.set_title(r'(e) $Re_0=25$',fontsize=20,loc='left')
ax6.set_title(r'(f) $Re_0=25$',fontsize=20,loc='left')
plt.colorbar(cs1,ax=ax1)
plt.colorbar(cs2,ax=ax2)
plt.colorbar(cs3,ax=ax3)
plt.colorbar(cs4,ax=ax4)
plt.colorbar(cs5,ax=ax5)
plt.colorbar(cs6,ax=ax6)
plt.tight_layout()
plt.savefig(opath+'pdfs_vort_pv_subplots_Re.pdf')
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
cs1 = ax1.contourf(NS117_vortpdf.xy[0,0,:,:],NS117_vortpdf.xy[1,0,:,:],NS117_vortpdf.pdf_timeavg[:NS117_vortpdf.ny,:NS117_vortpdf.nb],cmap=imola_map,extend='max')
ax1.plot(maxvort_117,NS117.y/NS117.z_enc,'k',lw=1)
ax1.scatter(maxvort_117_saddle,y_vort_117_saddle,100,color='k',marker='*')
ax1.axhline(NS117.z_ig/NS117.z_enc,0,0.05,color='k',linewidth=2)
ax1.axhline(NS117.z_is/NS117.z_enc,0,0.05,color='k',linewidth=2)
ax1.axhline(NS117.z_if/NS117.z_enc,0,0.05,color='k',linewidth=2)
cs2 = ax2.contourf(NS117_pvpdf.xy[0,0,:,:],NS117_pvpdf.xy[1,0,:,:],NS117_pvpdf.pdf_timeavg[:NS117_pvpdf.ny,:NS117_pvpdf.nb],cmap=imola_map,extend='max')
ax2.plot(maxpv_117,NS117.y/NS117.z_enc,'k',lw=1)
ax2.scatter(maxpv_117_saddle,y_pv_117_saddle,100,color='k',marker='*')
ax2.axhline(NS117.z_ig/NS117.z_enc,0,0.05,color='k',linewidth=2)
ax2.axhline(NS117.z_is/NS117.z_enc,0,0.05,color='k',linewidth=2)
ax2.axhline(NS117.z_if/NS117.z_enc,0,0.05,color='k',linewidth=2)
cs3 = ax3.contourf(NS25_vortpdf.xy[0,0,:,:],NS25_vortpdf.xy[1,0,:,:],NS25_vortpdf.pdf_timeavg[:NS25_vortpdf.ny,:NS25_vortpdf.nb],cmap=imola_map,extend='max')
ax3.plot(maxvort_25,NS25.y/NS25.z_enc,'k',lw=1)
ax3.scatter(maxvort_25_saddle,y_vort_25_saddle,100,color='k',marker='*')
ax3.axhline(NS25.z_ig/NS25.z_enc,0,0.05,color='k',linewidth=2)
ax3.axhline(NS25.z_is/NS25.z_enc,0,0.05,color='k',linewidth=2)
ax3.axhline(NS25.z_if/NS25.z_enc,0,0.05,color='k',linewidth=2)
cs4 = ax4.contourf(NS25_pvpdf.xy[0,0,:,:],NS25_pvpdf.xy[1,0,:,:],NS25_pvpdf.pdf_timeavg[:NS25_pvpdf.ny,:NS25_pvpdf.nb],cmap=imola_map,levels=np.linspace(0,0.24,9),extend='max')
ax4.plot(maxpv_25,NS25.y/NS25.z_enc,'k',lw=1)
ax4.scatter(maxpv_25_saddle,y_pv_25_saddle,100,color='k',marker='*')
ax4.axhline(NS25.z_ig/NS25.z_enc,0,0.05,color='k',linewidth=2)
ax4.axhline(NS25.z_is/NS25.z_enc,0,0.05,color='k',linewidth=2)
ax4.axhline(NS25.z_if/NS25.z_enc,0,0.05,color='k',linewidth=2)
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
plt.savefig(opath+'pdfs_vort_pv_subplots_Re_25_117.pdf')
plt.show()

f, (ax1,ax2) = plt.subplots(1,2,sharey='all',figsize=(10,5))
ax1.tick_params(bottom=True,top=True,left=True,right=True)
ax2.tick_params(bottom=True,top=True,left=True,right=True)
ax1.set_ylim(0,1.6)
ax1.set_xlim(-3,3)
ax2.set_xlim(-3,3)
cs1 = ax1.contourf(NS42_s1gradpdf.xy[0,0,:,:],NS42_s1gradpdf.xy[1,0,:,:],NS42_s1gradpdf.pdf_timeavg[:NS42_s1gradpdf.ny,:NS42_s1gradpdf.nb],cmap=imola_map,levels=np.linspace(0,0.4,11),extend='max')
cs2 = ax2.contourf(NS117_s1gradpdf.xy[0,0,:,:],NS117_s1gradpdf.xy[1,0,:,:],NS117_s1gradpdf.pdf_timeavg[:NS117_s1gradpdf.ny,:NS117_s1gradpdf.nb],cmap=imola_map,levels=np.linspace(0,0.4,11),extend='max')
ax1.set_xlabel(r'$\log_{10}(|\nabla b|^2/|\nabla b|_0^2)$')
ax2.set_xlabel(r'$\log_{10}(|\nabla b|^2/|\nabla b|_0^2)$')
ax1.set_ylabel(r'$z/z_\mathrm{enc}$')
ax1.set_title(r'(a) $Re_0=42$',loc='left',fontsize=20)
ax2.set_title(r'(b) $Re_0=117$',loc='left',fontsize=20)
plt.colorbar(cs1,ax=ax1)
plt.colorbar(cs2,ax=ax2)
plt.tight_layout()
plt.savefig(opath+'pdfs_s1grad_Re.pdf')
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
ax3.tick_params(bottom=True,top=True,left=True,right=True)
ax1.set_xlim(-4,2) 
ax2.set_xlim(-7,2)
ax1.set_ylim(0,1.6)
ax1.set_yticks([0,0.25,0.5,0.75,1,1.25,1.5])
ax2.set_xticks([-6,-4,-2,0,2])
cs1 = ax1.contourf(NS42_vortpdf_x_mean,NS42_vortpdf_y_mean,NS42_vortpdf.pdf_timeavg[:NS42_vortpdf.ny,:NS42_vortpdf.nb],cmap=imola_map)
ax1.plot(maxvort_42_avg,NS42.y/np.mean(NS42.z_enc),'k',lw=1)
ax1.scatter(maxvort_42_saddle_avg,y_vort_42_saddle_avg,100,color='k',marker='*')
ax1.axhline(np.mean(NS42.z_ig/NS42.z_enc),0,0.05,color='k',linewidth=2)
ax1.axhline(np.mean(NS42.z_is/NS42.z_enc),0,0.05,color='k',linewidth=2)
ax1.axhline(np.mean(NS42.z_if/NS42.z_enc),0,0.05,color='k',linewidth=2)
cs2 = ax2.contourf(NS42_pvpdf_x_mean,NS42_pvpdf_y_mean,NS42_pvpdf.pdf_timeavg[:NS42_pvpdf.ny,:NS42_pvpdf.nb],cmap=imola_map,levels=np.linspace(0,0.24,9))
ax2.plot(maxpv_42_avg,NS42.y/np.mean(NS42.z_enc),'k',lw=1)
ax2.scatter(maxpv_42_saddle_avg,y_pv_42_saddle_avg,100,color='k',marker='*')
ax2.axhline(np.mean(NS42.z_ig/NS42.z_enc),0,0.05,color='k',linewidth=2)
ax2.axhline(np.mean(NS42.z_is/NS42.z_enc),0,0.05,color='k',linewidth=2)
ax2.axhline(np.mean(NS42.z_if/NS42.z_enc),0,0.05,color='k',linewidth=2)
cs3 = ax3.contourf(S20_42_vortpdf_x_mean,S20_42_vortpdf_y_mean,S20_42_vortpdf.pdf_timeavg[:S20_42_vortpdf.ny,:S20_42_vortpdf.nb],cmap=imola_map,levels=np.linspace(0,0.4,9))
ax3.plot(maxvort_S20_42_avg,S20_42.y/np.mean(S20_42.z_enc),'k',lw=1)
ax3.scatter(maxvort_S20_42_saddle_avg,y_vort_S20_42_saddle_avg,100,color='k',marker='*')
ax3.axhline(np.mean(S20_42.z_ig/S20_42.z_enc),0,0.05,color='k',linewidth=2)
ax3.axhline(np.mean(S20_42.z_is/S20_42.z_enc),0,0.05,color='k',linewidth=2)
ax3.axhline(np.mean(S20_42.z_if/S20_42.z_enc),0,0.05,color='k',linewidth=2)
cs4 = ax4.contourf(S20_42_pvpdf_x_mean,S20_42_pvpdf_y_mean,S20_42_pvpdf.pdf_timeavg[:S20_42_pvpdf.ny,:S20_42_pvpdf.nb],cmap=imola_map,levels=np.linspace(0,0.24,9))
ax4.plot(maxpv_S20_42_avg,S20_42.y/np.mean(S20_42.z_enc),'k',lw=1)
ax4.scatter(maxpv_S20_42_saddle_avg,y_pv_S20_42_saddle_avg,100,color='k',marker='*')
ax4.axhline(np.mean(S20_42.z_ig/S20_42.z_enc),0,0.05,color='k',linewidth=2)
ax4.axhline(np.mean(S20_42.z_is/S20_42.z_enc),0,0.05,color='k',linewidth=2)
ax4.axhline(np.mean(S20_42.z_if/S20_42.z_enc),0,0.05,color='k',linewidth=2)
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
#plt.savefig(opath_42+'pdfs_vort_pv_subplots_S20_S0_timeavg.pdf')
plt.show()




f, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,sharex='col',sharey='all',figsize=(10,10))
ax1.tick_params(bottom=True,top=True,left=True,right=True)
ax2.tick_params(bottom=True,top=True,left=True,right=True)
ax3.tick_params(bottom=True,top=True,left=True,right=True)
ax3.tick_params(bottom=True,top=True,left=True,right=True)
ax1.set_xlim(-4,2) 
ax2.set_xlim(-7,2)
ax1.set_ylim(0,1.6)
ax1.set_yticks([0,0.25,0.5,0.75,1,1.25,1.5])
ax2.set_xticks([-6,-4,-2,0,2])
cs1 = ax1.contourf(NS42_vortpdf.xy[0,:,:],NS42_vortpdf.xy[1,:,:],NS42_vortpdf.pdf[2,:NS42_vortpdf.ny,:NS42_vortpdf.nb],cmap=imola_map)
ax1.plot(maxvort_42,NS42.y/NS42.z_enc[2],'k',lw=1)
ax1.scatter(maxvort_42_saddle,y_vort_42_saddle,100,color='k',marker='*')
ax1.axhline(NS42.z_ig[2]/NS42.z_enc[2],0,0.05,color='k',linewidth=2)
ax1.axhline(NS42.z_is[2]/NS42.z_enc[2],0,0.05,color='k',linewidth=2)
ax1.axhline(NS42.z_if[2]/NS42.z_enc[2],0,0.05,color='k',linewidth=2)
cs2 = ax2.contourf(NS42_pvpdf.xy[0,:,:],NS42_pvpdf.xy[1,:,:],NS42_pvpdf.pdf[2,:NS42_pvpdf.ny,:NS42_pvpdf.nb],cmap=imola_map,levels=np.linspace(0,0.24,9))
ax2.plot(maxpv_42,NS42.y/NS42.z_enc[2],'k',lw=1)
ax2.scatter(maxpv_42_saddle,y_pv_42_saddle,100,color='k',marker='*')
ax2.axhline(NS42.z_ig[2]/NS42.z_enc[2],0,0.05,color='k',linewidth=2)
ax2.axhline(NS42.z_is[2]/NS42.z_enc[2],0,0.05,color='k',linewidth=2)
ax2.axhline(NS42.z_if[2]/NS42.z_enc[2],0,0.05,color='k',linewidth=2)
cs3 = ax3.contourf(S20_42_vortpdf.xy[0,:,:],S20_42_vortpdf.xy[1,:,:],S20_42_vortpdf.pdf[2,:S20_42_vortpdf.ny,:S20_42_vortpdf.nb],cmap=imola_map,levels=np.linspace(0,0.4,9))
ax3.plot(maxvort_S20_42,S20_42.y/S20_42.z_enc[2],'k',lw=1)
ax3.scatter(maxvort_S20_42_saddle,y_vort_S20_42_saddle,100,color='k',marker='*')
ax3.axhline(S20_42.z_ig[2]/S20_42.z_enc[2],0,0.05,color='k',linewidth=2)
ax3.axhline(S20_42.z_is[2]/S20_42.z_enc[2],0,0.05,color='k',linewidth=2)
ax3.axhline(S20_42.z_if[2]/S20_42.z_enc[2],0,0.05,color='k',linewidth=2)
cs4 = ax4.contourf(S20_42_pvpdf.xy[0,:,:],S20_42_pvpdf.xy[1,:,:],S20_42_pvpdf.pdf[2,:S20_42_pvpdf.ny,:S20_42_pvpdf.nb],cmap=imola_map,levels=np.linspace(0,0.24,9))
ax4.plot(maxpv_S20_42,S20_42.y/S20_42.z_enc[2],'k',lw=1)
ax4.scatter(maxpv_S20_42_saddle,y_pv_S20_42_saddle,100,color='k',marker='*')
ax4.axhline(S20_42.z_ig[2]/S20_42.z_enc[2],0,0.05,color='k',linewidth=2)
ax4.axhline(S20_42.z_is[2]/S20_42.z_enc[2],0,0.05,color='k',linewidth=2)
ax4.axhline(S20_42.z_if[2]/S20_42.z_enc[2],0,0.05,color='k',linewidth=2)
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
plt.savefig(opath_42+'pdfs_vort_pv_subplots_S20_S0_21.pdf')
plt.show()

f, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,sharex='col',sharey='all',figsize=(10,10))
ax1.tick_params(bottom=True,top=True,left=True,right=True)
ax2.tick_params(bottom=True,top=True,left=True,right=True)
ax3.tick_params(bottom=True,top=True,left=True,right=True)
ax3.tick_params(bottom=True,top=True,left=True,right=True)
ax1.set_xlim(-4,2) 
ax2.set_xlim(-7,2)
ax1.set_ylim(0,1.6)
ax1.set_yticks([0,0.25,0.5,0.75,1,1.25,1.5])
ax2.set_xticks([-6,-4,-2,0,2])
cs1 = ax1.contourf(NS42_vortpdf.xy[0,:,:],NS42_vortpdf.xy[1,:,:],NS42_vortpdf.pdf[1,:NS42_vortpdf.ny,:NS42_vortpdf.nb],cmap=imola_map)
ax1.plot(maxvort_42,NS42.y/NS42.z_enc[1],'k',lw=1)
ax1.scatter(maxvort_42_saddle,y_vort_42_saddle,100,color='k',marker='*')
ax1.axhline(NS42.z_ig[1]/NS42.z_enc[1],0,0.05,color='k',linewidth=2)
ax1.axhline(NS42.z_is[1]/NS42.z_enc[1],0,0.05,color='k',linewidth=2)
ax1.axhline(NS42.z_if[1]/NS42.z_enc[1],0,0.05,color='k',linewidth=2)
cs2 = ax2.contourf(NS42_pvpdf.xy[0,:,:],NS42_pvpdf.xy[1,:,:],NS42_pvpdf.pdf[1,:NS42_pvpdf.ny,:NS42_pvpdf.nb],cmap=imola_map,levels=np.linspace(0,0.24,9))
ax2.plot(maxpv_42,NS42.y/NS42.z_enc[1],'k',lw=1)
ax2.scatter(maxpv_42_saddle,y_pv_42_saddle,100,color='k',marker='*')
ax2.axhline(NS42.z_ig[1]/NS42.z_enc[1],0,0.05,color='k',linewidth=2)
ax2.axhline(NS42.z_is[1]/NS42.z_enc[1],0,0.05,color='k',linewidth=2)
ax2.axhline(NS42.z_if[1]/NS42.z_enc[1],0,0.05,color='k',linewidth=2)
cs3 = ax3.contourf(S20_42_vortpdf.xy[0,:,:],S20_42_vortpdf.xy[1,:,:],S20_42_vortpdf.pdf[1,:S20_42_vortpdf.ny,:S20_42_vortpdf.nb],cmap=imola_map,levels=np.linspace(0,0.4,9))
ax3.plot(maxvort_S20_42,S20_42.y/S20_42.z_enc[1],'k',lw=1)
ax3.scatter(maxvort_S20_42_saddle,y_vort_S20_42_saddle,100,color='k',marker='*')
ax3.axhline(S20_42.z_ig[1]/S20_42.z_enc[1],0,0.05,color='k',linewidth=2)
ax3.axhline(S20_42.z_is[1]/S20_42.z_enc[1],0,0.05,color='k',linewidth=2)
ax3.axhline(S20_42.z_if[1]/S20_42.z_enc[1],0,0.05,color='k',linewidth=2)
cs4 = ax4.contourf(S20_42_pvpdf.xy[0,:,:],S20_42_pvpdf.xy[1,:,:],S20_42_pvpdf.pdf[1,:S20_42_pvpdf.ny,:S20_42_pvpdf.nb],cmap=imola_map,levels=np.linspace(0,0.24,9))
ax4.plot(maxpv_S20_42,S20_42.y/S20_42.z_enc[1],'k',lw=1)
ax4.scatter(maxpv_S20_42_saddle,y_pv_S20_42_saddle,100,color='k',marker='*')
ax4.axhline(S20_42.z_ig[1]/S20_42.z_enc[1],0,0.05,color='k',linewidth=2)
ax4.axhline(S20_42.z_is[1]/S20_42.z_enc[1],0,0.05,color='k',linewidth=2)
ax4.axhline(S20_42.z_if[1]/S20_42.z_enc[1],0,0.05,color='k',linewidth=2)
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
#plt.savefig(opath_42+'pdfs_vort_pv_subplots_S20_S0_20.pdf')
plt.show()

f, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,sharex='col',sharey='all',figsize=(10,10))
ax1.tick_params(bottom=True,top=True,left=True,right=True)
ax2.tick_params(bottom=True,top=True,left=True,right=True)
ax3.tick_params(bottom=True,top=True,left=True,right=True)
ax3.tick_params(bottom=True,top=True,left=True,right=True)
ax1.set_xlim(-4,2) 
ax2.set_xlim(-7,2)
ax1.set_ylim(0,1.6)
ax1.set_yticks([0,0.25,0.5,0.75,1,1.25,1.5])
ax2.set_xticks([-6,-4,-2,0,2])
cs1 = ax1.contourf(NS42_vortpdf.xy[0,:,:],NS42_vortpdf.xy[1,:,:],NS42_vortpdf.pdf[0,:NS42_vortpdf.ny,:NS42_vortpdf.nb],cmap=imola_map)
ax1.plot(maxvort_42,NS42.y/NS42.z_enc[0],'k',lw=1)
ax1.scatter(maxvort_42_saddle,y_vort_42_saddle,100,color='k',marker='*')
ax1.axhline(NS42.z_ig[0]/NS42.z_enc[0],0,0.05,color='k',linewidth=2)
ax1.axhline(NS42.z_is[0]/NS42.z_enc[0],0,0.05,color='k',linewidth=2)
ax1.axhline(NS42.z_if[0]/NS42.z_enc[0],0,0.05,color='k',linewidth=2)
cs2 = ax2.contourf(NS42_pvpdf.xy[0,:,:],NS42_pvpdf.xy[1,:,:],NS42_pvpdf.pdf[0,:NS42_pvpdf.ny,:NS42_pvpdf.nb],cmap=imola_map,levels=np.linspace(0,0.24,9))
ax2.plot(maxpv_42,NS42.y/NS42.z_enc[0],'k',lw=1)
ax2.scatter(maxpv_42_saddle,y_pv_42_saddle,100,color='k',marker='*')
ax2.axhline(NS42.z_ig[0]/NS42.z_enc[0],0,0.05,color='k',linewidth=2)
ax2.axhline(NS42.z_is[0]/NS42.z_enc[0],0,0.05,color='k',linewidth=2)
ax2.axhline(NS42.z_if[0]/NS42.z_enc[0],0,0.05,color='k',linewidth=2)
cs3 = ax3.contourf(S20_42_vortpdf.xy[0,:,:],S20_42_vortpdf.xy[1,:,:],S20_42_vortpdf.pdf[0,:S20_42_vortpdf.ny,:S20_42_vortpdf.nb],cmap=imola_map,levels=np.linspace(0,0.4,9))
ax3.plot(maxvort_S20_42,S20_42.y/S20_42.z_enc[0],'k',lw=1)
ax3.scatter(maxvort_S20_42_saddle,y_vort_S20_42_saddle,100,color='k',marker='*')
ax3.axhline(S20_42.z_ig[0]/S20_42.z_enc[0],0,0.05,color='k',linewidth=2)
ax3.axhline(S20_42.z_is[0]/S20_42.z_enc[0],0,0.05,color='k',linewidth=2)
ax3.axhline(S20_42.z_if[0]/S20_42.z_enc[0],0,0.05,color='k',linewidth=2)
cs4 = ax4.contourf(S20_42_pvpdf.xy[0,:,:],S20_42_pvpdf.xy[1,:,:],S20_42_pvpdf.pdf[0,:S20_42_pvpdf.ny,:S20_42_pvpdf.nb],cmap=imola_map,levels=np.linspace(0,0.24,9))
ax4.plot(maxpv_S20_42,S20_42.y/S20_42.z_enc[0],'k',lw=1)
ax4.scatter(maxpv_S20_42_saddle,y_pv_S20_42_saddle,100,color='k',marker='*')
ax4.axhline(S20_42.z_ig[0]/S20_42.z_enc[0],0,0.05,color='k',linewidth=2)
ax4.axhline(S20_42.z_is[0]/S20_42.z_enc[0],0,0.05,color='k',linewidth=2)
ax4.axhline(S20_42.z_if[0]/S20_42.z_enc[0],0,0.05,color='k',linewidth=2)
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
plt.savefig(opath_42+'pdfs_vort_pv_subplots_S20_S0_19.pdf')
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

