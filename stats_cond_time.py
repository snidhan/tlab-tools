#####################################################################
# Modules

from ReadStats import Statistics, Conditional_Stats
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.style
import matplotlib.lines as mlines
from matplotlib import rc

rc('text', usetex=True)
rc('text.latex', preamble=r"\usepackage{fourier}")
rc('font', family='serif')
rc('font', size=24)
rc('axes', linewidth=1.5)
rc('axes', labelsize=24)
rc('lines', linewidth=2)

opath = '/home/mpim/m300551/Figures/JAS2020/'

#######################################################################
# Constants

nu = 1./25000.
B0 = 0.005
N = np.sqrt(3)
L0 = (B0/N**3)**0.5

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

#######################################################################
# Stats

path_1 = '/scratch/local1/m300551/ForKatherine/qCBL_3D/Re042/2560x576x2560/'
path_2 = '/scratch/local1/m300551/ForKatherine/qCBL_3D/Re042/2560x704x2560/'
path_3 = '/scratch/local1/m300551/ForKatherine/qCBL_3D/Re042/2560x896x2560/'
path_S20 = '/scratch/local1/m300551/ForKatherine/qCBL_3D/Re042/3072x960x4608-S20/'

path_vort = 'stats/gate-vorticity/gate-2-08/'
path_vort_S20 = 'stats/gate-vorticity/gate-1-24/'

# Conventional

NS42_1 = Statistics(path_1+'stats/pdftimes/avg20500-53000.nc')
NS42_2 = Statistics(path_2+'stats/pdftimes/avg60000-74500.nc')
NS42_3 = Statistics(path_3+'stats/pdftimes/avg83000-127500.nc')

S20 = Statistics(path_S20+'stats/pdftimes/avg42000-148000.nc')

NS42_s1_var_zif_1 = [NS42_1.r2S[n,NS42_1.z_if_arg[n]] for n in range(0,NS42_1.t_len)]
NS42_s1_var_zif_2 = [NS42_2.r2S[n,NS42_2.z_if_arg[n]] for n in range(0,NS42_2.t_len)]
NS42_s1_var_zif_3 = [NS42_3.r2S[n,NS42_3.z_if_arg[n]] for n in range(0,NS42_3.t_len)]

s1_var_zif = NS42_s1_var_zif_1 + NS42_s1_var_zif_2 + NS42_s1_var_zif_3

NS42_w_var_zif_1 = [NS42_1.Ryy[n,NS42_1.z_if_arg[n]] for n in range(0,NS42_1.t_len)]
NS42_w_var_zif_2 = [NS42_2.Ryy[n,NS42_2.z_if_arg[n]] for n in range(0,NS42_2.t_len)]
NS42_w_var_zif_3 = [NS42_3.Ryy[n,NS42_3.z_if_arg[n]] for n in range(0,NS42_3.t_len)]

w_var_zif = NS42_w_var_zif_1 + NS42_w_var_zif_2 + NS42_w_var_zif_3

NS42_s1_flux_zif_1 = [NS42_1.Rsv[n,NS42_1.z_if_arg[n]] for n in range(0,NS42_1.t_len)]
NS42_s1_flux_zif_2 = [NS42_2.Rsv[n,NS42_2.z_if_arg[n]] for n in range(0,NS42_2.t_len)]
NS42_s1_flux_zif_3 = [NS42_3.Rsv[n,NS42_3.z_if_arg[n]] for n in range(0,NS42_3.t_len)]

s1_flux_zif = NS42_s1_flux_zif_1 + NS42_s1_flux_zif_2 + NS42_s1_flux_zif_3

S20_s1_var_zif = [S20.r2S[n,S20.z_if_arg[n]] for n in range(0,S20.t_len)]
S20_w_var_zif = [S20.Ryy[n,S20.z_if_arg[n]] for n in range(0,S20.t_len)]
S20_s1_flux_zif = [S20.Rsv[n,S20.z_if_arg[n]] for n in range(0,S20.t_len)]

## Conditional ##

NS42_vort_int_1 = Conditional_Stats(path_1+path_vort+'int20500-53000.nc',path_1+path_vort+'Partition1/cavg20500-53000.nc',path_1+path_vort+'Partition2/cavg20500-53000.nc')
NS42_vort_int_2 = Conditional_Stats(path_2+path_vort+'int60000-74500.nc',path_2+path_vort+'Partition1/cavg60000-74500.nc',path_2+path_vort+'Partition2/cavg60000-74500.nc')
NS42_vort_int_3 = Conditional_Stats(path_3+path_vort+'int83000-127500.nc',path_3+path_vort+'Partition1/cavg83000-127500.nc',path_3+path_vort+'Partition2/cavg83000-127500.nc')

S20_vort_int = Conditional_Stats(path_S20+path_vort_S20+'int42000-148000.nc',path_S20+path_vort_S20+'Partition1/cavg42000-148000.nc',path_S20+path_vort_S20+'Partition2/cavg42000-148000.nc')

# Vorticity #

NS42_vort_turbareafrac_zif_1 = [NS42_vort_int_1.int2[n,NS42_1.z_if_arg[n]] for n in range(0,NS42_1.t_len)]
NS42_vort_turbareafrac_zif_2 = [NS42_vort_int_2.int2[n,NS42_2.z_if_arg[n]] for n in range(0,NS42_2.t_len)]
NS42_vort_turbareafrac_zif_3 = [NS42_vort_int_3.int2[n,NS42_3.z_if_arg[n]] for n in range(0,NS42_3.t_len)]

NS42_vort_turbareafrac_zif = NS42_vort_turbareafrac_zif_1 + NS42_vort_turbareafrac_zif_2 +  NS42_vort_turbareafrac_zif_3

S20_vort_turbareafrac_zif = [S20_vort_int.int2[n,S20.z_if_arg[n]] for n in range(0,S20.t_len)]

# Non-turbulent

NS42_vort_p1_s1_mean_zif_1 = [NS42_vort_int_1.P1S1Mom1[n,NS42_1.z_if_arg[n]] for n in range(0,NS42_1.t_len)]
NS42_vort_p1_s1_mean_zif_2 = [NS42_vort_int_2.P1S1Mom1[n,NS42_2.z_if_arg[n]] for n in range(0,NS42_2.t_len)]
NS42_vort_p1_s1_mean_zif_3 = [NS42_vort_int_3.P1S1Mom1[n,NS42_3.z_if_arg[n]] for n in range(0,NS42_3.t_len)]

NS42_vort_p1_s1_mean_zif = NS42_vort_p1_s1_mean_zif_1 + NS42_vort_p1_s1_mean_zif_2 + NS42_vort_p1_s1_mean_zif_3

NS42_vort_p1_w_mean_zif_1 = [NS42_vort_int_1.P1VMom1[n,NS42_1.z_if_arg[n]] for n in range(0,NS42_1.t_len)]
NS42_vort_p1_w_mean_zif_2 = [NS42_vort_int_2.P1VMom1[n,NS42_2.z_if_arg[n]] for n in range(0,NS42_2.t_len)]
NS42_vort_p1_w_mean_zif_3 = [NS42_vort_int_3.P1VMom1[n,NS42_3.z_if_arg[n]] for n in range(0,NS42_3.t_len)]

NS42_vort_p1_w_mean_zif = NS42_vort_p1_w_mean_zif_1 + NS42_vort_p1_w_mean_zif_2 + NS42_vort_p1_w_mean_zif_3

NS42_vort_p1_s1_var_zif_1 = [NS42_vort_int_1.P1S1Mom2[n,NS42_1.z_if_arg[n]] for n in range(0,NS42_1.t_len)]
NS42_vort_p1_s1_var_zif_2 = [NS42_vort_int_2.P1S1Mom2[n,NS42_2.z_if_arg[n]] for n in range(0,NS42_2.t_len)]
NS42_vort_p1_s1_var_zif_3 = [NS42_vort_int_3.P1S1Mom2[n,NS42_3.z_if_arg[n]] for n in range(0,NS42_3.t_len)]

NS42_vort_p1_s1_var_zif = NS42_vort_p1_s1_var_zif_1 + NS42_vort_p1_s1_var_zif_2 + NS42_vort_p1_s1_var_zif_3

NS42_vort_p1_w_var_zif_1 = [NS42_vort_int_1.P1VMom2[n,NS42_1.z_if_arg[n]] for n in range(0,NS42_1.t_len)]
NS42_vort_p1_w_var_zif_2 = [NS42_vort_int_2.P1VMom2[n,NS42_2.z_if_arg[n]] for n in range(0,NS42_2.t_len)]
NS42_vort_p1_w_var_zif_3 = [NS42_vort_int_3.P1VMom2[n,NS42_3.z_if_arg[n]] for n in range(0,NS42_3.t_len)]

NS42_vort_p1_w_var_zif = NS42_vort_p1_w_var_zif_1 + NS42_vort_p1_w_var_zif_2 + NS42_vort_p1_w_var_zif_3

NS42_vort_p1_v1_zif_1 = [NS42_vort_int_1.P1v1[n,NS42_1.z_if_arg[n]] for n in range(0,NS42_1.t_len)]
NS42_vort_p1_v1_zif_2 = [NS42_vort_int_2.P1v1[n,NS42_2.z_if_arg[n]] for n in range(0,NS42_2.t_len)]
NS42_vort_p1_v1_zif_3 = [NS42_vort_int_3.P1v1[n,NS42_3.z_if_arg[n]] for n in range(0,NS42_3.t_len)]

NS42_vort_p1_v1_zif = NS42_vort_p1_v1_zif_1 + NS42_vort_p1_v1_zif_2 + NS42_vort_p1_v1_zif_3

S20_vort_p1_s1_mean_zif = [S20_vort_int.P1S1Mom1[n,S20.z_if_arg[n]] for n in range(0,S20.t_len)]
S20_vort_p1_w_mean_zif = [S20_vort_int.P1VMom1[n,S20.z_if_arg[n]] for n in range(0,S20.t_len)]
S20_vort_p1_s1_var_zif = [S20_vort_int.P1S1Mom2[n,S20.z_if_arg[n]] for n in range(0,S20.t_len)]
S20_vort_p1_w_var_zif = [S20_vort_int.P1VMom2[n,S20.z_if_arg[n]] for n in range(0,S20.t_len)]
S20_vort_p1_v1_zif = [S20_vort_int.P1v1[n,S20.z_if_arg[n]] for n in range(0,S20.t_len)]

# Turbulent

NS42_vort_p2_s1_mean_zif_1 = [NS42_vort_int_1.P2S1Mom1[n,NS42_1.z_if_arg[n]] for n in range(0,NS42_1.t_len)]
NS42_vort_p2_s1_mean_zif_2 = [NS42_vort_int_2.P2S1Mom1[n,NS42_2.z_if_arg[n]] for n in range(0,NS42_2.t_len)]
NS42_vort_p2_s1_mean_zif_3 = [NS42_vort_int_3.P2S1Mom1[n,NS42_3.z_if_arg[n]] for n in range(0,NS42_3.t_len)]

NS42_vort_p2_s1_mean_zif = NS42_vort_p2_s1_mean_zif_1 + NS42_vort_p2_s1_mean_zif_2 + NS42_vort_p2_s1_mean_zif_3

NS42_vort_p2_w_mean_zif_1 = [NS42_vort_int_1.P2VMom1[n,NS42_1.z_if_arg[n]] for n in range(0,NS42_1.t_len)]
NS42_vort_p2_w_mean_zif_2 = [NS42_vort_int_2.P2VMom1[n,NS42_2.z_if_arg[n]] for n in range(0,NS42_2.t_len)]
NS42_vort_p2_w_mean_zif_3 = [NS42_vort_int_3.P2VMom1[n,NS42_3.z_if_arg[n]] for n in range(0,NS42_3.t_len)]

NS42_vort_p2_w_mean_zif = NS42_vort_p2_w_mean_zif_1 + NS42_vort_p2_w_mean_zif_2 + NS42_vort_p2_w_mean_zif_3

NS42_vort_p2_s1_var_zif_1 = [NS42_vort_int_1.P2S1Mom2[n,NS42_1.z_if_arg[n]] for n in range(0,NS42_1.t_len)]
NS42_vort_p2_s1_var_zif_2 = [NS42_vort_int_2.P2S1Mom2[n,NS42_2.z_if_arg[n]] for n in range(0,NS42_2.t_len)]
NS42_vort_p2_s1_var_zif_3 = [NS42_vort_int_3.P2S1Mom2[n,NS42_3.z_if_arg[n]] for n in range(0,NS42_3.t_len)]

NS42_vort_p2_s1_var_zif = NS42_vort_p2_s1_var_zif_1 + NS42_vort_p2_s1_var_zif_2 + NS42_vort_p2_s1_var_zif_3

NS42_vort_p2_w_var_zif_1 = [NS42_vort_int_1.P2VMom2[n,NS42_1.z_if_arg[n]] for n in range(0,NS42_1.t_len)]
NS42_vort_p2_w_var_zif_2 = [NS42_vort_int_2.P2VMom2[n,NS42_2.z_if_arg[n]] for n in range(0,NS42_2.t_len)]
NS42_vort_p2_w_var_zif_3 = [NS42_vort_int_3.P2VMom2[n,NS42_3.z_if_arg[n]] for n in range(0,NS42_3.t_len)]

NS42_vort_p2_w_var_zif = NS42_vort_p2_w_var_zif_1 + NS42_vort_p2_w_var_zif_2 + NS42_vort_p2_w_var_zif_3

NS42_vort_p2_v1_zif_1 = [NS42_vort_int_1.P2v1[n,NS42_1.z_if_arg[n]] for n in range(0,NS42_1.t_len)]
NS42_vort_p2_v1_zif_2 = [NS42_vort_int_2.P2v1[n,NS42_2.z_if_arg[n]] for n in range(0,NS42_2.t_len)]
NS42_vort_p2_v1_zif_3 = [NS42_vort_int_3.P2v1[n,NS42_3.z_if_arg[n]] for n in range(0,NS42_3.t_len)]

NS42_vort_p2_v1_zif = NS42_vort_p2_v1_zif_1 + NS42_vort_p2_v1_zif_2 + NS42_vort_p2_v1_zif_3

S20_vort_p2_s1_mean_zif = [S20_vort_int.P2S1Mom1[n,S20.z_if_arg[n]] for n in range(0,S20.t_len)]
S20_vort_p2_w_mean_zif = [S20_vort_int.P2VMom1[n,S20.z_if_arg[n]] for n in range(0,S20.t_len)]
S20_vort_p2_s1_var_zif = [S20_vort_int.P2S1Mom2[n,S20.z_if_arg[n]] for n in range(0,S20.t_len)]
S20_vort_p2_w_var_zif = [S20_vort_int.P2VMom2[n,S20.z_if_arg[n]] for n in range(0,S20.t_len)]
S20_vort_p2_v1_zif = [S20_vort_int.P2v1[n,S20.z_if_arg[n]] for n in range(0,S20.t_len)]

# time

time_1 = [NS42_1.z_enc[n]/L0 for n in range(0,NS42_1.t_len)]
time_2 = [NS42_2.z_enc[n]/L0 for n in range(0,NS42_2.t_len)]
time_3 = [NS42_3.z_enc[n]/L0 for n in range(0,NS42_3.t_len)]

time = time_1 + time_2 + time_3

z_if = np.concatenate([NS42_1.z_if,NS42_2.z_if,NS42_3.z_if])
z_ig = np.concatenate([NS42_1.z_ig,NS42_2.z_ig,NS42_3.z_ig])
z_enc = np.concatenate([NS42_1.z_enc,NS42_2.z_enc,NS42_3.z_enc])
w_enc = (B0*z_enc)**(1./3.)


rho_bw = NS42_1.Rsv/(NS42_1.r2S*NS42_1.Ryy)**0.5
rho_bw_turb = (NS42_vort_int_1.P2v1[-4:-1,:]-NS42_vort_int_1.P2S1Mom1[-4:-1,:]*NS42_vort_int_1.P2VMom1[-4:-1,:])/(NS42_vort_int_1.P2S1Mom2[-4:-1,:]*NS42_vort_int_1.P2VMom2[-4:-1,:])**0.5
rho_bw_S20 = S20.Rsv/(S20.r2S*S20.Ryy)**0.5
rho_bw_S20_turb = (S20_vort_int.P2v1[4:7,:]-S20_vort_int.P2S1Mom1[4:7,:]*S20_vort_int.P2VMom1[4:7,:])/(S20_vort_int.P2S1Mom2[4:7,:]*S20_vort_int.P2VMom2[4:7,:])**0.5
rho_bw_zif = np.array(s1_flux_zif)/(np.array(s1_var_zif)*np.array(w_var_zif))**0.5
rho_bw_zif_turb = (np.array(NS42_vort_p2_v1_zif)-np.array(NS42_vort_p2_s1_mean_zif)*np.array(NS42_vort_p2_w_mean_zif))/(np.array(NS42_vort_p2_s1_var_zif)*np.array(NS42_vort_p2_w_var_zif))**0.5
rho_bw_zif_S20 = np.array(S20_s1_flux_zif)/(np.array(S20_s1_var_zif)*np.array(S20_w_var_zif))**0.5
rho_bw_zif_S20_turb =  (np.array(S20_vort_p2_v1_zif)-np.array(S20_vort_p2_s1_mean_zif)*np.array(S20_vort_p2_w_mean_zif))/(np.array(S20_vort_p2_s1_var_zif)*np.array(S20_vort_p2_w_var_zif))**0.5
#######################################################################
# Plot

blues = matplotlib.cm.get_cmap('Blues')

f, (ax1,ax2) = plt.subplots(1,2,figsize=(10,5))
ax1.grid(True,linewidth=1.5)
ax2.grid(True,linewidth=1.5)
ax1.tick_params(bottom=False,left=False)
ax2.tick_params(bottom=False,left=False)
ax1.set_ylim(0,1.6)
ax2.set_ylim(0.4,1)
ax2.set_xlim(15,30)
ax1.plot(np.mean(NS42_vort_int_1.int2[-4:-1,:],axis=0),NS42_1.y/np.mean(NS42_1.z_enc[-4:-1]),c=blues(0.5),label=r'$Fr_0=0$')
ax1.plot(np.mean(S20_vort_int.int2[4:7,:],axis=0),S20.y/np.mean(S20.z_enc[4:7]),c=blues(0.9),label=r'$Fr_0=20$')
ax2.plot(time[1:-1],runningmean(NS42_vort_turbareafrac_zif,1),c=blues(0.5))
ax2.plot(S20.z_enc[1:-1]/L0,runningmean(S20_vort_turbareafrac_zif,1),c=blues(0.9))
ax1.axhline(np.mean(NS42_1.z_if[-4:-1]/NS42_1.z_enc[-4:-1]),0,0.05,c=blues(0.5))
ax1.axhline(np.mean(S20.z_if[4:7]/S20.z_enc[4:7]),0,0.05,c=blues(0.9))
ax1.set_xlabel(r'$a_\mathrm{T}$')
ax1.set_ylabel(r'$z/z_\mathrm{enc}$')
ax1.set_title(r'(a)',fontsize=20,loc='left')
ax2.set_xlabel(r'$z_\mathrm{enc}/L_0$')
ax2.set_ylabel(r'$(a_\mathrm{T})_{z_{i,f}}$')
ax2.set_title(r'(b)',fontsize=24,loc='left')
ax1.legend(loc='lower left',fontsize=24,handlelength=1)
plt.tight_layout()
plt.savefig(opath+'Fig5.pdf',bbox_inches='tight')
plt.show()

f, (ax1,ax2) = plt.subplots(1,2,figsize=(10,5))
ax1.grid(True,linewidth=1.5)
ax2.grid(True,linewidth=1.5)
ax1.tick_params(bottom=False,left=False)
ax2.tick_params(bottom=False,left=False)
ax1.set_xlim(1,1.4)
ax2.set_xlim(-0.15,0.15)
ax1.set_ylim(1,1.4)
ax2.set_ylim(1,1.4)
ax1.plot(np.mean(NS42_vort_int_1.P2S1Mom1[-4:-1,:],axis=0)/np.mean(N**2*NS42_1.z_enc[-4:-1]),NS42_1.y/np.mean(NS42_1.z_enc[-4:-1]),c=blues(0.5),label=r'$Fr_0=0$')
ax1.plot(np.mean(S20_vort_int.P2S1Mom1[4:7,:],axis=0)/np.mean(N**2*S20.z_enc[4:7]),S20.y/np.mean(S20.z_enc[4:7]),c=blues(0.9),label=r'$Fr_0=20$')
ax1.plot(np.mean(NS42_vort_int_1.P1S1Mom1[-4:-1,:],axis=0)/np.mean(N**2*NS42_1.z_enc[-4:-1]),NS42_1.y/np.mean(NS42_1.z_enc[-4:-1]),c=blues(0.5),ls='--')
ax1.plot(np.mean(S20_vort_int.P1S1Mom1[4:7,:],axis=0)/np.mean(N**2*S20.z_enc[4:7]),S20.y/np.mean(S20.z_enc[4:7]),c=blues(0.9),ls='--')
ax1.plot(S20.y/np.mean(S20.z_enc[4:7]),S20.y/np.mean(S20.z_enc[4:7]),'k:')
ax2.plot(np.mean(NS42_vort_int_1.P2VMom1[-4:-1,:],axis=0)/(np.mean(NS42_1.z_enc[-4:-1])*B0)**(1./3.),NS42_1.y/np.mean(NS42_1.z_enc[-4:-1]),c=blues(0.5))
ax2.plot(np.mean(S20_vort_int.P2VMom1[4:7,:],axis=0)/(np.mean(S20.z_enc[4:7])*B0)**(1./3.),S20.y/np.mean(S20.z_enc[4:7]),c=blues(0.9))
ax2.plot(np.mean(NS42_vort_int_1.P1VMom1[-4:-1,:],axis=0)/(np.mean(NS42_1.z_enc[-4:-1])*B0)**(1./3.),NS42_1.y/np.mean(NS42_1.z_enc[-4:-1]),c=blues(0.5),ls='--')
ax2.plot(np.mean(S20_vort_int.P1VMom1[4:7,:],axis=0)/(np.mean(S20.z_enc[4:7])*B0)**(1./3.),S20.y/np.mean(S20.z_enc[4:7]),c=blues(0.9),ls='--')
black_solid = mlines.Line2D([],[],c='k',label='Turbulent')
black_dashed = mlines.Line2D([],[],c='k',ls='--',label='Non- \n turbulent')
black_dot = mlines.Line2D([],[],c='k',ls=':',label=r'$N_0^2z$')
ax1.axhline(np.mean(NS42_1.z_if[-4:-1]/NS42_1.z_enc[-4:-1]),0,0.05,c=blues(0.5))
ax1.axhline(np.mean(S20.z_if[4:7]/S20.z_enc[4:7]),0,0.05,c=blues(0.9))
ax2.axhline(np.mean(NS42_1.z_if[-4:-1]/NS42_1.z_enc[-4:-1]),0,0.05,c=blues(0.5))
ax2.axhline(np.mean(S20.z_if[4:7]/S20.z_enc[4:7]),0,0.05,c=blues(0.9))
ax1.set_xlabel(r'$\langle b \rangle / b_\mathrm{enc}$')
ax1.set_ylabel(r'$z/z_\mathrm{enc}$')
ax1.set_title(r'(a)',fontsize=24,loc='left')
ax2.set_xlabel(r'$\langle w \rangle / w_*$')
ax2.set_ylabel(r'$z/z_\mathrm{enc}$')
ax2.set_title(r'(b)',fontsize=24,loc='left')
leg1 = ax1.legend(loc='upper left',fontsize=18,borderaxespad=0.1,handlelength=1)
leg2 = ax1.legend(handles=[black_solid,black_dashed,black_dot],loc='lower right',fontsize=18,borderaxespad=0.1,handlelength=1,labelspacing=0.3)
ax1.add_artist(leg1)
plt.tight_layout()
plt.savefig(opath+'Fig7.pdf',bbox_inches='tight')
plt.show()


f, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(10,10))
ax1.grid(True,linewidth=1.5)
ax2.grid(True,linewidth=1.5)
ax3.grid(True,linewidth=1.5)
ax4.grid(True,linewidth=1.5)
ax1.tick_params(bottom=False,left=False)
ax2.tick_params(bottom=False,left=False)
ax3.tick_params(bottom=False,left=False)
ax4.tick_params(bottom=False,left=False)
ax1.set_xlim(0,2.5)
ax1.set_ylim(1,1.4)
ax2.set_xlim(15,30)
ax2.set_ylim(0,2)
ax3.set_ylim(1,1.4)
ax3.set_xlim(0,0.6)
ax4.set_ylim(0.2,0.6)
ax4.set_xlim(15,30)
ax1.plot(np.mean(np.sqrt(NS42_vort_int_1.P2S1Mom2[-4:-1,:]),axis=0)/(N**2*L0),NS42_1.y/np.mean(NS42_1.z_enc[-4:-1]),c=blues(0.5),label=r'$Fr_0=0$')
ax1.plot(np.mean(np.sqrt(S20_vort_int.P2S1Mom2[4:7,:]),axis=0)/(N**2*L0),S20.y/np.mean(S20.z_enc[4:7],axis=0),c=blues(0.9),label=r'$Fr_0=20$')
ax2.plot(time[1:-1],runningmean(np.sqrt(NS42_vort_p2_s1_var_zif),1)/(N**2*L0),c=blues(0.5))
ax2.plot(S20.z_enc[1:-1]/L0,runningmean(np.sqrt(S20_vort_p2_s1_var_zif),1)/(N**2*L0),c=blues(0.9))
ax3.plot(np.mean(np.sqrt(NS42_vort_int_1.P2VMom2[-4:-1,:]),axis=0)/np.mean(w_enc[4:7]),NS42_1.y/np.mean(NS42_1.z_enc[-4:-1]),c=blues(0.5))
ax3.plot(np.mean(np.sqrt(S20_vort_int.P2VMom2[4:7,:]),axis=0)/np.mean((B0*S20.z_enc[4:7])**(1./3.)),S20.y/np.mean(S20.z_enc[4:7]),c=blues(0.9))
ax4.plot(time[1:-1],runningmean(np.sqrt(NS42_vort_p2_w_var_zif)/w_enc,1),c=blues(0.5),label=r'$Fr_0=0$')
ax4.plot(S20.z_enc[1:-1]/L0,runningmean(np.sqrt(S20_vort_p2_w_var_zif)/(B0*S20.z_enc)**(1./3.),1),c=blues(0.9),label=r'$Fr_0=20$')
ax1.axhline(np.mean(NS42_1.z_if[-4:-1]/NS42_1.z_enc[-4:-1]),0,0.05,c=blues(0.5))
ax1.axhline(np.mean(S20.z_if[4:7]/S20.z_enc[4:7]),0,0.05,c=blues(0.9))
ax3.axhline(np.mean(NS42_1.z_if[-4:-1]/NS42_1.z_enc[-4:-1]),0,0.05,c=blues(0.5))
ax3.axhline(np.mean(S20.z_if[4:7]/S20.z_enc[4:7]),0,0.05,c=blues(0.9))
ax1.set_xlabel(r'$(b_\mathrm{rms})_\mathrm{T} /(N^2L_0)$')
ax1.set_ylabel(r'$z/z_\mathrm{enc}$')
ax1.set_title('(a)',fontsize=24,loc='left')
ax2.set_xlabel(r'$z_\mathrm{enc}/L_0$')
ax2.set_ylabel(r'$((b_\mathrm{rms})_\mathrm{T})_{z_{i,f}}/(N^2L_0)$')
ax2.set_title('(b)',fontsize=24,loc='left')
ax3.set_xlabel(r'$(w_\mathrm{rms})_\mathrm{T} / w_*$')
ax3.set_ylabel(r'$z/z_\mathrm{enc}$')
ax3.set_title('(c)',fontsize=24,loc='left')
ax4.set_xlabel(r'$z_\mathrm{enc}/L_0$')
ax4.set_ylabel(r'$((w_\mathrm{rms})_\mathrm{T})_{z_{i,f}}/w_*$')
ax4.set_title('(d)',fontsize=24,loc='left')
ax1.legend(loc='best',fontsize=24,handlelength=1,borderaxespad=0.1)
plt.tight_layout()
plt.savefig(opath+'Fig8.pdf',bbox_inches='tight')
plt.show()


f, (ax1,ax2) = plt.subplots(1,2,figsize=(10,5))
ax1.grid(True,linewidth=1.5)
ax2.grid(True,linewidth=1.5)
ax1.tick_params(bottom=False,left=False)
ax2.tick_params(bottom=False,left=False)
ax1.set_ylim(1,1.4)
ax1.set_xlim(-0.2,0)
ax2.set_ylim(-0.02,0.2)
ax2.set_xlim(15,30)
ax1.plot(np.mean(rho_bw_turb,axis=0),NS42_1.y/np.mean(NS42_1.z_enc[-4:-1]),c=blues(0.5))
ax1.plot(np.mean(rho_bw_S20_turb,axis=0),S20.y/np.mean(S20.z_enc[4:7]),c=blues(0.9))
ax2.plot(time[1:-1],-runningmean(rho_bw_zif_turb,1),c=blues(0.5),label=r'$Fr_0=0$')
ax2.plot(S20.z_enc[1:-1]/L0,-runningmean(rho_bw_zif_S20_turb,1),c=blues(0.9),label=r'$Fr_0=20$')
ax1.axhline(np.mean(NS42_1.z_if[-4:-1]/NS42_1.z_enc[-4:-1]),0,0.05,c=blues(0.5))
ax1.axhline(np.mean(S20.z_if[4:7]/S20.z_enc[4:7]),0,0.05,c=blues(0.9))
ax1.set_xlabel(r'$(\rho_{bw})_\mathrm{T}$')
ax1.set_ylabel(r'$z/z_\mathrm{enc}$')
ax1.set_title('(a)',fontsize=20,loc='left')
ax2.set_xlabel(r'$z_\mathrm{enc}/L_0$')
ax2.set_ylabel(r'$-((\rho_{bw})_\mathrm{T})_{z_{i,f}}$')
ax2.set_title('(b)',fontsize=24,loc='left')
ax2.legend(loc='best',fontsize=24,handlelength=1,borderaxespad=0.1)
plt.tight_layout()
plt.savefig(opath+'Fig9.pdf',bbox_inches='tight')
plt.show()

f, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(10,10))
ax1.grid(True,linewidth=1.5)
ax2.axis('off')
ax3.grid(True,linewidth=1.5)
ax4.grid(True,linewidth=1.5)
ax1.tick_params(bottom=False,left=False)
ax3.tick_params(bottom=False,left=False)
ax4.tick_params(bottom=False,left=False)
ax1.set_xlim(-0.2,0.05)
ax1.set_ylim(1,1.4)
ax3.set_ylim(1,1.4)
ax3.set_xlim(-0.4,0.2)
ax4.set_ylim(0,0.4)
ax4.set_xlim(15,30)
ax1.plot(np.mean(NS42_vort_int_1.int2[-4:-1,:]*(NS42_vort_int_1.P2v1[-4:-1,:]-NS42_vort_int_1.P2S1Mom1[-4:-1,:]*NS42_vort_int_1.P2VMom1[-4:-1,:])/B0,axis=0),NS42_1.y/np.mean(NS42_1.z_enc[-4:-1]),c=blues(0.5))
ax1.plot(np.mean(S20_vort_int.int2[4:7,:]*(S20_vort_int.P2v1[4:7,:]-S20_vort_int.P2S1Mom1[4:7,:]*S20_vort_int.P2VMom1[4:7,:])/B0,axis=0),S20.y/np.mean(S20.z_enc[4:7]),c=blues(0.9))
ax1.plot(np.mean(NS42_vort_int_1.int1[-4:-1,:]*(NS42_vort_int_1.P1v1[-4:-1,:]-NS42_vort_int_1.P1S1Mom1[-4:-1,:]*NS42_vort_int_1.P1VMom1[-4:-1,:])/B0,axis=0),NS42_1.y/np.mean(NS42_1.z_enc[-4:-1]),c=blues(0.5),ls='--')
ax1.plot(np.mean(S20_vort_int.int1[4:7,:]*(S20_vort_int.P1v1[4:7,:]-S20_vort_int.P1S1Mom1[4:7,:]*S20_vort_int.P1VMom1[4:7,:])/B0,axis=0),S20.y/np.mean(S20.z_enc[4:7]),c=blues(0.9),ls='--')
ax1.plot(np.mean(NS42_vort_int_1.int2[-4:-1,:]*NS42_vort_int_1.int1[-4:-1,:]*(NS42_vort_int_1.P2S1Mom1[-4:-1,:]-NS42_vort_int_1.P1S1Mom1[-4:-1,:])*(NS42_vort_int_1.P2VMom1[-4:-1,:]-NS42_vort_int_1.P1VMom1[-4:-1,:])/B0,axis=0),NS42_1.y/np.mean(NS42_1.z_enc[-4:-1]),c=blues(0.5),ls='-.')
ax1.plot(np.mean(S20_vort_int.int2[4:7,:]*S20_vort_int.int1[4:7,:]*(S20_vort_int.P2S1Mom1[4:7,:]-S20_vort_int.P1S1Mom1[4:7,:])*(S20_vort_int.P2VMom1[4:7,:]-S20_vort_int.P1VMom1[4:7,:])/B0,axis=0),S20.y/np.mean(S20.z_enc[4:7]),c=blues(0.9),ls='-.')
ax1.plot(np.mean(NS42_1.Rsv[-4:-1,:]/B0,axis=0),NS42_1.y/np.mean(NS42_1.z_enc[-4:-1]),c=blues(0.5),ls=':')
ax1.plot(np.mean(S20.Rsv[4:7,:]/B0,axis=0),S20.y/np.mean(S20.z_enc[4:7]),c=blues(0.9),ls=':')
ax3.plot(np.mean(NS42_vort_int_1.P2v1[-4:-1,:]-NS42_vort_int_1.P2S1Mom1[-4:-1,:]*NS42_vort_int_1.P2VMom1[-4:-1,:],axis=0)/B0,NS42_1.y/np.mean(NS42_1.z_enc[-4:-1]),c=blues(0.5))
ax3.plot(np.mean(S20_vort_int.P2v1[4:7,:]-S20_vort_int.P2S1Mom1[4:7,:]*S20_vort_int.P2VMom1[4:7,:],axis=0)/B0,S20.y/np.mean(S20.z_enc[4:7]),c=blues(0.9))
ax4.plot(time[1:-1],-runningmean(np.array(NS42_vort_p2_v1_zif)-np.array(NS42_vort_p2_s1_mean_zif)*np.array(NS42_vort_p2_w_mean_zif),1)/B0,c=blues(0.5),label=r'$Fr_0=0$')
ax4.plot(S20.z_enc[1:-1]/L0,-runningmean(np.array(S20_vort_p2_v1_zif)-np.array(S20_vort_p2_s1_mean_zif)*np.array(S20_vort_p2_w_mean_zif),1)/B0,c=blues(0.9),label=r'$Fr_0=20$')
ax1.axhline(np.mean(NS42_1.z_if[-4:-1]/NS42_1.z_enc[-4:-1]),0,0.05,c=blues(0.5))
ax1.axhline(np.mean(S20.z_if[4:7]/S20.z_enc[4:7]),0,0.05,c=blues(0.9))
ax3.axhline(np.mean(NS42_1.z_if[-4:-1]/NS42_1.z_enc[-4:-1]),0,0.05,c=blues(0.5))
ax3.axhline(np.mean(S20.z_if[4:7]/S20.z_enc[4:7]),0,0.05,c=blues(0.9))
ax1.set_xlabel(r'$f /B_0$')
ax1.set_ylabel(r'$z/z_\mathrm{enc}$')
ax1.set_title('(a)',fontsize=24,loc='left')
ax3.set_xlabel(r'$\langle b^\prime w^\prime\rangle_\mathrm{T}/B_0$')
ax3.set_ylabel(r'$z/z_\mathrm{enc}$')
ax3.set_title(r'(b)',fontsize=24,loc='left')
ax4.set_xlabel(r'$z_\mathrm{enc}/L_0$')
ax4.set_ylabel(r'$-(\langle b^\prime w^\prime \rangle_\mathrm{T})_{z_{i,f}}/B_0$')
ax4.set_title(r'(c)',fontsize=24,loc='left')
ax4.legend(loc='best',fontsize=24,handlelength=1)
black_dot = mlines.Line2D([],[],c='k',ls=':',label=r'$\langle b^\prime w^\prime \rangle$')
black_solid = mlines.Line2D([],[],c='k',label=r'$a_\mathrm{T}\langle b^\prime w^\prime \rangle_\mathrm{T}$')
black_dashed = mlines.Line2D([],[],c='k',ls='--',label=r'$a_\mathrm{NT}\langle b^\prime w^\prime \rangle_\mathrm{NT}$')
black_dashdot = mlines.Line2D([],[],c='k',ls='-.',label=r'$a_\mathrm{T}a_\mathrm{NT}(\langle b \rangle_\mathrm{T} - \langle b \rangle_\mathrm{NT})(\langle w \rangle_\mathrm{T} - \langle w \rangle_\mathrm{NT})$')
ax2.legend(handles=[black_dot,black_solid,black_dashed,black_dashdot],loc='upper right',fontsize=20,frameon=False,handlelength=1)
plt.tight_layout()
plt.savefig(opath+'Fig6.pdf',bbox_inches='tight')
plt.show()






