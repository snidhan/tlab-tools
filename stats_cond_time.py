#####################################################################
# Modules

from ReadStats import Statistics, Conditional_Stats
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
path_pv = 'stats/gate-pv/gate-2-41/'

path_vort_S20 = 'stats/gate-vorticity/gate-1-24/'
path_pv_S20 = 'stats/gate-pv/gate-3-92/'

# Conventional

NS42_1 = Statistics(path_1+'stats/pdftimes/avg20500-53000.nc')
NS42_2 = Statistics(path_2+'stats/pdftimes/avg60000-74500.nc')
NS42_3 = Statistics(path_3+'stats/pdftimes/avg83000-127500.nc')

S20 = Statistics(path_S20+'stats/pdftimes/avg42000-148000.nc')

NS42_s1_mean_zig_1 = [NS42_1.rS[n,NS42_1.z_ig_arg[n]] for n in range(0,NS42_1.t_len)]
NS42_s1_mean_zig_2 = [NS42_2.rS[n,NS42_2.z_ig_arg[n]] for n in range(0,NS42_2.t_len)]
NS42_s1_mean_zig_3 = [NS42_3.rS[n,NS42_3.z_ig_arg[n]] for n in range(0,NS42_3.t_len)]

s1_mean_zig = NS42_s1_mean_zig_1 + NS42_s1_mean_zig_2 + NS42_s1_mean_zig_3

NS42_s1_grad_zig_1 = [NS42_1.rS_y[n,NS42_1.z_ig_arg[n]] for n in range(0,NS42_1.t_len)]
NS42_s1_grad_zig_2 = [NS42_2.rS_y[n,NS42_2.z_ig_arg[n]] for n in range(0,NS42_2.t_len)]
NS42_s1_grad_zig_3 = [NS42_3.rS_y[n,NS42_3.z_ig_arg[n]] for n in range(0,NS42_3.t_len)]

s1_grad_zig = NS42_s1_grad_zig_1 + NS42_s1_grad_zig_2 + NS42_s1_grad_zig_3

NS42_s1_var_zig_1 = [NS42_1.r2S[n,NS42_1.z_ig_arg[n]] for n in range(0,NS42_1.t_len)]
NS42_s1_var_zig_2 = [NS42_2.r2S[n,NS42_2.z_ig_arg[n]] for n in range(0,NS42_2.t_len)]
NS42_s1_var_zig_3 = [NS42_3.r2S[n,NS42_3.z_ig_arg[n]] for n in range(0,NS42_3.t_len)]

s1_var_zig = NS42_s1_var_zig_1 + NS42_s1_var_zig_2 + NS42_s1_var_zig_3

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

S20_s1_mean_zig = [S20.rS[n,S20.z_ig_arg[n]] for n in range(0,S20.t_len)]
S20_s1_grad_zig = [S20.rS_y[n,S20.z_ig_arg[n]] for n in range(0,S20.t_len)]
S20_s1_var_zig = [S20.r2S[n,S20.z_ig_arg[n]] for n in range(0,S20.t_len)]
S20_s1_var_zif = [S20.r2S[n,S20.z_if_arg[n]] for n in range(0,S20.t_len)]
S20_w_var_zif = [S20.Ryy[n,S20.z_if_arg[n]] for n in range(0,S20.t_len)]
S20_s1_flux_zif = [S20.Rsv[n,S20.z_if_arg[n]] for n in range(0,S20.t_len)]

## Conditional ##

NS42_vort_int_1 = Conditional_Stats(path_1+path_vort+'int20500-53000.nc',path_1+path_vort+'Partition1/cavg20500-53000.nc',path_1+path_vort+'Partition2/cavg20500-53000.nc')
NS42_vort_int_2 = Conditional_Stats(path_2+path_vort+'int60000-74500.nc',path_2+path_vort+'Partition1/cavg60000-74500.nc',path_2+path_vort+'Partition2/cavg60000-74500.nc')
NS42_vort_int_3 = Conditional_Stats(path_3+path_vort+'int83000-127500.nc',path_3+path_vort+'Partition1/cavg83000-127500.nc',path_3+path_vort+'Partition2/cavg83000-127500.nc')

NS42_pv_int_1 = Conditional_Stats(path_1+path_pv+'int20500-53000.nc',path_1+path_pv+'Partition1/cavg20500-53000.nc',path_1+path_pv+'Partition2/cavg20500-53000.nc')
NS42_pv_int_2 = Conditional_Stats(path_2+path_pv+'int60000-74500.nc',path_2+path_pv+'Partition1/cavg60000-74500.nc',path_2+path_pv+'Partition2/cavg60000-74500.nc')
NS42_pv_int_3 = Conditional_Stats(path_3+path_pv+'int83000-127500.nc',path_3+path_pv+'Partition1/cavg83000-127500.nc',path_3+path_pv+'Partition2/cavg83000-127500.nc')

S20_vort_int = Conditional_Stats(path_S20+path_vort_S20+'int42000-148000.nc',path_S20+path_vort_S20+'Partition1/cavg42000-148000.nc',path_S20+path_vort_S20+'Partition2/cavg42000-148000.nc')
S20_pv_int = Conditional_Stats(path_S20+path_pv_S20+'int42000-148000.nc',path_S20+path_pv_S20+'Partition1/cavg42000-148000.nc',path_S20+path_pv_S20+'Partition2/cavg42000-148000.nc')

# Vorticity #

NS42_vort_turbareafrac_zig_1 = [NS42_vort_int_1.int2[n,NS42_1.z_ig_arg[n]] for n in range(0,NS42_1.t_len)]
NS42_vort_turbareafrac_zig_2 = [NS42_vort_int_2.int2[n,NS42_2.z_ig_arg[n]] for n in range(0,NS42_2.t_len)]
NS42_vort_turbareafrac_zig_3 = [NS42_vort_int_3.int2[n,NS42_3.z_ig_arg[n]] for n in range(0,NS42_3.t_len)]

NS42_vort_turbareafrac_zig = NS42_vort_turbareafrac_zig_1 + NS42_vort_turbareafrac_zig_2 +  NS42_vort_turbareafrac_zig_3

NS42_vort_turbareafrac_zif_1 = [NS42_vort_int_1.int2[n,NS42_1.z_if_arg[n]] for n in range(0,NS42_1.t_len)]
NS42_vort_turbareafrac_zif_2 = [NS42_vort_int_2.int2[n,NS42_2.z_if_arg[n]] for n in range(0,NS42_2.t_len)]
NS42_vort_turbareafrac_zif_3 = [NS42_vort_int_3.int2[n,NS42_3.z_if_arg[n]] for n in range(0,NS42_3.t_len)]

NS42_vort_turbareafrac_zif = NS42_vort_turbareafrac_zif_1 + NS42_vort_turbareafrac_zif_2 +  NS42_vort_turbareafrac_zif_3

S20_vort_turbareafrac_zig = [S20_vort_int.int2[n,S20.z_ig_arg[n]] for n in range(0,S20.t_len)]
S20_vort_turbareafrac_zif = [S20_vort_int.int2[n,S20.z_if_arg[n]] for n in range(0,S20.t_len)]

# Non-turbulent

NS42_vort_p1_s1_mean_zig_1 = [NS42_vort_int_1.P1S1Mom1[n,NS42_1.z_ig_arg[n]] for n in range(0,NS42_1.t_len)]
NS42_vort_p1_s1_mean_zig_2 = [NS42_vort_int_2.P1S1Mom1[n,NS42_2.z_ig_arg[n]] for n in range(0,NS42_2.t_len)]
NS42_vort_p1_s1_mean_zig_3 = [NS42_vort_int_3.P1S1Mom1[n,NS42_3.z_ig_arg[n]] for n in range(0,NS42_3.t_len)]

NS42_vort_p1_s1_mean_zig = NS42_vort_p1_s1_mean_zig_1 + NS42_vort_p1_s1_mean_zig_2 + NS42_vort_p1_s1_mean_zig_3

NS42_vort_p1_s1_mean_zif_1 = [NS42_vort_int_1.P1S1Mom1[n,NS42_1.z_if_arg[n]] for n in range(0,NS42_1.t_len)]
NS42_vort_p1_s1_mean_zif_2 = [NS42_vort_int_2.P1S1Mom1[n,NS42_2.z_if_arg[n]] for n in range(0,NS42_2.t_len)]
NS42_vort_p1_s1_mean_zif_3 = [NS42_vort_int_3.P1S1Mom1[n,NS42_3.z_if_arg[n]] for n in range(0,NS42_3.t_len)]

NS42_vort_p1_s1_mean_zif = NS42_vort_p1_s1_mean_zif_1 + NS42_vort_p1_s1_mean_zif_2 + NS42_vort_p1_s1_mean_zif_3

NS42_vort_p1_w_mean_zif_1 = [NS42_vort_int_1.P1VMom1[n,NS42_1.z_if_arg[n]] for n in range(0,NS42_1.t_len)]
NS42_vort_p1_w_mean_zif_2 = [NS42_vort_int_2.P1VMom1[n,NS42_2.z_if_arg[n]] for n in range(0,NS42_2.t_len)]
NS42_vort_p1_w_mean_zif_3 = [NS42_vort_int_3.P1VMom1[n,NS42_3.z_if_arg[n]] for n in range(0,NS42_3.t_len)]

NS42_vort_p1_w_mean_zif = NS42_vort_p1_w_mean_zif_1 + NS42_vort_p1_w_mean_zif_2 + NS42_vort_p1_w_mean_zif_3

NS42_vort_p1_s1_var_zig_1 = [NS42_vort_int_1.P1S1Mom2[n,NS42_1.z_ig_arg[n]] for n in range(0,NS42_1.t_len)]
NS42_vort_p1_s1_var_zig_2 = [NS42_vort_int_2.P1S1Mom2[n,NS42_2.z_ig_arg[n]] for n in range(0,NS42_2.t_len)]
NS42_vort_p1_s1_var_zig_3 = [NS42_vort_int_3.P1S1Mom2[n,NS42_3.z_ig_arg[n]] for n in range(0,NS42_3.t_len)]

NS42_vort_p1_s1_var_zig = NS42_vort_p1_s1_var_zig_1 + NS42_vort_p1_s1_var_zig_2 + NS42_vort_p1_s1_var_zig_3

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

S20_vort_p1_s1_mean_zig = [S20_vort_int.P1S1Mom1[n,S20.z_ig_arg[n]] for n in range(0,S20.t_len)]
S20_vort_p1_s1_mean_zif = [S20_vort_int.P1S1Mom1[n,S20.z_if_arg[n]] for n in range(0,S20.t_len)]
S20_vort_p1_w_mean_zif = [S20_vort_int.P1VMom1[n,S20.z_if_arg[n]] for n in range(0,S20.t_len)]
S20_vort_p1_s1_var_zig = [S20_vort_int.P1S1Mom2[n,S20.z_ig_arg[n]] for n in range(0,S20.t_len)]
S20_vort_p1_s1_var_zif = [S20_vort_int.P1S1Mom2[n,S20.z_if_arg[n]] for n in range(0,S20.t_len)]
S20_vort_p1_w_var_zif = [S20_vort_int.P1VMom2[n,S20.z_if_arg[n]] for n in range(0,S20.t_len)]
S20_vort_p1_v1_zif = [S20_vort_int.P1v1[n,S20.z_if_arg[n]] for n in range(0,S20.t_len)]

# Turbulent

NS42_vort_p2_s1_mean_zig_1 = [NS42_vort_int_1.P2S1Mom1[n,NS42_1.z_ig_arg[n]] for n in range(0,NS42_1.t_len)]
NS42_vort_p2_s1_mean_zig_2 = [NS42_vort_int_2.P2S1Mom1[n,NS42_2.z_ig_arg[n]] for n in range(0,NS42_2.t_len)]
NS42_vort_p2_s1_mean_zig_3 = [NS42_vort_int_3.P2S1Mom1[n,NS42_3.z_ig_arg[n]] for n in range(0,NS42_3.t_len)]

NS42_vort_p2_s1_mean_zig = NS42_vort_p2_s1_mean_zig_1 + NS42_vort_p2_s1_mean_zig_2 + NS42_vort_p2_s1_mean_zig_3

NS42_vort_p2_s1_mean_zif_1 = [NS42_vort_int_1.P2S1Mom1[n,NS42_1.z_if_arg[n]] for n in range(0,NS42_1.t_len)]
NS42_vort_p2_s1_mean_zif_2 = [NS42_vort_int_2.P2S1Mom1[n,NS42_2.z_if_arg[n]] for n in range(0,NS42_2.t_len)]
NS42_vort_p2_s1_mean_zif_3 = [NS42_vort_int_3.P2S1Mom1[n,NS42_3.z_if_arg[n]] for n in range(0,NS42_3.t_len)]

NS42_vort_p2_s1_mean_zif = NS42_vort_p2_s1_mean_zif_1 + NS42_vort_p2_s1_mean_zif_2 + NS42_vort_p2_s1_mean_zif_3

NS42_vort_p2_w_mean_zif_1 = [NS42_vort_int_1.P2VMom1[n,NS42_1.z_if_arg[n]] for n in range(0,NS42_1.t_len)]
NS42_vort_p2_w_mean_zif_2 = [NS42_vort_int_2.P2VMom1[n,NS42_2.z_if_arg[n]] for n in range(0,NS42_2.t_len)]
NS42_vort_p2_w_mean_zif_3 = [NS42_vort_int_3.P2VMom1[n,NS42_3.z_if_arg[n]] for n in range(0,NS42_3.t_len)]

NS42_vort_p2_w_mean_zif = NS42_vort_p2_w_mean_zif_1 + NS42_vort_p2_w_mean_zif_2 + NS42_vort_p2_w_mean_zif_3

NS42_vort_p2_s1_var_zig_1 = [NS42_vort_int_1.P2S1Mom2[n,NS42_1.z_ig_arg[n]] for n in range(0,NS42_1.t_len)]
NS42_vort_p2_s1_var_zig_2 = [NS42_vort_int_2.P2S1Mom2[n,NS42_2.z_ig_arg[n]] for n in range(0,NS42_2.t_len)]
NS42_vort_p2_s1_var_zig_3 = [NS42_vort_int_3.P2S1Mom2[n,NS42_3.z_ig_arg[n]] for n in range(0,NS42_3.t_len)]

NS42_vort_p2_s1_var_zig = NS42_vort_p2_s1_var_zig_1 + NS42_vort_p2_s1_var_zig_2 + NS42_vort_p2_s1_var_zig_3

NS42_vort_p2_s1_var_zif_1 = [NS42_vort_int_1.P2S1Mom2[n,NS42_1.z_if_arg[n]] for n in range(0,NS42_1.t_len)]
NS42_vort_p2_s1_var_zif_2 = [NS42_vort_int_2.P2S1Mom2[n,NS42_2.z_if_arg[n]] for n in range(0,NS42_2.t_len)]
NS42_vort_p2_s1_var_zif_3 = [NS42_vort_int_3.P2S1Mom2[n,NS42_3.z_if_arg[n]] for n in range(0,NS42_3.t_len)]

NS42_vort_p2_s1_var_zif = NS42_vort_p2_s1_var_zif_1 + NS42_vort_p2_s1_var_zif_2 + NS42_vort_p2_s1_var_zif_3

NS42_vort_p2_w_var_zif_1 = [NS42_vort_int_1.P2VMom2[n,NS42_1.z_if_arg[n]] for n in range(0,NS42_1.t_len)]
NS42_vort_p2_w_var_zif_2 = [NS42_vort_int_2.P2VMom2[n,NS42_2.z_if_arg[n]] for n in range(0,NS42_2.t_len)]
NS42_vort_p2_w_var_zif_3 = [NS42_vort_int_3.P2VMom2[n,NS42_3.z_if_arg[n]] for n in range(0,NS42_3.t_len)]

NS42_vort_p2_w_var_zif = NS42_vort_p2_w_var_zif_1 + NS42_vort_p2_w_var_zif_2 + NS42_vort_p2_w_var_zif_3

NS42_vort_p2_w_var_zig_1 = [NS42_vort_int_1.P2VMom2[n,NS42_1.z_ig_arg[n]] for n in range(0,NS42_1.t_len)]
NS42_vort_p2_w_var_zig_2 = [NS42_vort_int_2.P2VMom2[n,NS42_2.z_ig_arg[n]] for n in range(0,NS42_2.t_len)]
NS42_vort_p2_w_var_zig_3 = [NS42_vort_int_3.P2VMom2[n,NS42_3.z_ig_arg[n]] for n in range(0,NS42_3.t_len)]

NS42_vort_p2_w_var_zig = NS42_vort_p2_w_var_zig_1 + NS42_vort_p2_w_var_zig_2 + NS42_vort_p2_w_var_zig_3

NS42_vort_p2_v1_zif_1 = [NS42_vort_int_1.P2v1[n,NS42_1.z_if_arg[n]] for n in range(0,NS42_1.t_len)]
NS42_vort_p2_v1_zif_2 = [NS42_vort_int_2.P2v1[n,NS42_2.z_if_arg[n]] for n in range(0,NS42_2.t_len)]
NS42_vort_p2_v1_zif_3 = [NS42_vort_int_3.P2v1[n,NS42_3.z_if_arg[n]] for n in range(0,NS42_3.t_len)]

NS42_vort_p2_v1_zif = NS42_vort_p2_v1_zif_1 + NS42_vort_p2_v1_zif_2 + NS42_vort_p2_v1_zif_3

S20_vort_p2_s1_mean_zig = [S20_vort_int.P2S1Mom1[n,S20.z_ig_arg[n]] for n in range(0,S20.t_len)]
S20_vort_p2_s1_mean_zif = [S20_vort_int.P2S1Mom1[n,S20.z_if_arg[n]] for n in range(0,S20.t_len)]
S20_vort_p2_w_mean_zif = [S20_vort_int.P2VMom1[n,S20.z_if_arg[n]] for n in range(0,S20.t_len)]
S20_vort_p2_s1_var_zig = [S20_vort_int.P2S1Mom2[n,S20.z_ig_arg[n]] for n in range(0,S20.t_len)]
S20_vort_p2_s1_var_zif = [S20_vort_int.P2S1Mom2[n,S20.z_if_arg[n]] for n in range(0,S20.t_len)]
S20_vort_p2_w_var_zif = [S20_vort_int.P2VMom2[n,S20.z_if_arg[n]] for n in range(0,S20.t_len)]
S20_vort_p2_w_var_zig = [S20_vort_int.P2VMom2[n,S20.z_ig_arg[n]] for n in range(0,S20.t_len)]
S20_vort_p2_v1_zif = [S20_vort_int.P2v1[n,S20.z_if_arg[n]] for n in range(0,S20.t_len)]

# PV #

NS42_pv_turbareafrac_zig_1 = [NS42_pv_int_1.int2[n,NS42_1.z_ig_arg[n]] for n in range(0,NS42_1.t_len)]
NS42_pv_turbareafrac_zig_2 = [NS42_pv_int_2.int2[n,NS42_2.z_ig_arg[n]] for n in range(0,NS42_2.t_len)]
NS42_pv_turbareafrac_zig_3 = [NS42_pv_int_3.int2[n,NS42_3.z_ig_arg[n]] for n in range(0,NS42_3.t_len)]

NS42_pv_turbareafrac_zig = NS42_pv_turbareafrac_zig_1 + NS42_pv_turbareafrac_zig_2 +  NS42_pv_turbareafrac_zig_3

NS42_pv_turbareafrac_zif_1 = [NS42_pv_int_1.int2[n,NS42_1.z_if_arg[n]] for n in range(0,NS42_1.t_len)]
NS42_pv_turbareafrac_zif_2 = [NS42_pv_int_2.int2[n,NS42_2.z_if_arg[n]] for n in range(0,NS42_2.t_len)]
NS42_pv_turbareafrac_zif_3 = [NS42_pv_int_3.int2[n,NS42_3.z_if_arg[n]] for n in range(0,NS42_3.t_len)]

NS42_pv_turbareafrac_zif = NS42_pv_turbareafrac_zif_1 + NS42_pv_turbareafrac_zif_2 +  NS42_pv_turbareafrac_zif_3

S20_pv_turbareafrac_zig = [S20_pv_int.int2[n,S20.z_ig_arg[n]] for n in range(0,S20.t_len)]
S20_pv_turbareafrac_zif = [S20_pv_int.int2[n,S20.z_if_arg[n]] for n in range(0,S20.t_len)]

# Non-turbulent

NS42_pv_p1_s1_mean_zig_1 = [NS42_pv_int_1.P1S1Mom1[n,NS42_1.z_ig_arg[n]] for n in range(0,NS42_1.t_len)]
NS42_pv_p1_s1_mean_zig_2 = [NS42_pv_int_2.P1S1Mom1[n,NS42_2.z_ig_arg[n]] for n in range(0,NS42_2.t_len)]
NS42_pv_p1_s1_mean_zig_3 = [NS42_pv_int_3.P1S1Mom1[n,NS42_3.z_ig_arg[n]] for n in range(0,NS42_3.t_len)]

NS42_pv_p1_s1_mean_zig = NS42_pv_p1_s1_mean_zig_1 + NS42_pv_p1_s1_mean_zig_2 + NS42_pv_p1_s1_mean_zig_3

NS42_pv_p1_s1_mean_zif_1 = [NS42_pv_int_1.P1S1Mom1[n,NS42_1.z_if_arg[n]] for n in range(0,NS42_1.t_len)]
NS42_pv_p1_s1_mean_zif_2 = [NS42_pv_int_2.P1S1Mom1[n,NS42_2.z_if_arg[n]] for n in range(0,NS42_2.t_len)]
NS42_pv_p1_s1_mean_zif_3 = [NS42_pv_int_3.P1S1Mom1[n,NS42_3.z_if_arg[n]] for n in range(0,NS42_3.t_len)]

NS42_pv_p1_s1_mean_zif = NS42_pv_p1_s1_mean_zif_1 + NS42_pv_p1_s1_mean_zif_2 + NS42_pv_p1_s1_mean_zif_3

NS42_pv_p1_w_mean_zif_1 = [NS42_pv_int_1.P1VMom1[n,NS42_1.z_if_arg[n]] for n in range(0,NS42_1.t_len)]
NS42_pv_p1_w_mean_zif_2 = [NS42_pv_int_2.P1VMom1[n,NS42_2.z_if_arg[n]] for n in range(0,NS42_2.t_len)]
NS42_pv_p1_w_mean_zif_3 = [NS42_pv_int_3.P1VMom1[n,NS42_3.z_if_arg[n]] for n in range(0,NS42_3.t_len)]

NS42_pv_p1_w_mean_zif = NS42_pv_p1_w_mean_zif_1 + NS42_pv_p1_w_mean_zif_2 + NS42_pv_p1_w_mean_zif_3

NS42_pv_p1_s1_var_zig_1 = [NS42_pv_int_1.P1S1Mom2[n,NS42_1.z_ig_arg[n]] for n in range(0,NS42_1.t_len)]
NS42_pv_p1_s1_var_zig_2 = [NS42_pv_int_2.P1S1Mom2[n,NS42_2.z_ig_arg[n]] for n in range(0,NS42_2.t_len)]
NS42_pv_p1_s1_var_zig_3 = [NS42_pv_int_3.P1S1Mom2[n,NS42_3.z_ig_arg[n]] for n in range(0,NS42_3.t_len)]

NS42_pv_p1_s1_var_zig = NS42_pv_p1_s1_var_zig_1 + NS42_pv_p1_s1_var_zig_2 + NS42_pv_p1_s1_var_zig_3

NS42_pv_p1_s1_var_zif_1 = [NS42_pv_int_1.P1S1Mom2[n,NS42_1.z_if_arg[n]] for n in range(0,NS42_1.t_len)]
NS42_pv_p1_s1_var_zif_2 = [NS42_pv_int_2.P1S1Mom2[n,NS42_2.z_if_arg[n]] for n in range(0,NS42_2.t_len)]
NS42_pv_p1_s1_var_zif_3 = [NS42_pv_int_3.P1S1Mom2[n,NS42_3.z_if_arg[n]] for n in range(0,NS42_3.t_len)]

NS42_pv_p1_s1_var_zif = NS42_pv_p1_s1_var_zif_1 + NS42_pv_p1_s1_var_zif_2 + NS42_pv_p1_s1_var_zif_3

NS42_pv_p1_v_var_zif_1 = [NS42_pv_int_1.P1VMom2[n,NS42_1.z_if_arg[n]] for n in range(0,NS42_1.t_len)]
NS42_pv_p1_v_var_zif_2 = [NS42_pv_int_2.P1VMom2[n,NS42_2.z_if_arg[n]] for n in range(0,NS42_2.t_len)]
NS42_pv_p1_v_var_zif_3 = [NS42_pv_int_3.P1VMom2[n,NS42_3.z_if_arg[n]] for n in range(0,NS42_3.t_len)]

NS42_pv_p1_v_var_zif = NS42_pv_p1_v_var_zif_1 + NS42_pv_p1_v_var_zif_2 + NS42_pv_p1_v_var_zif_3

NS42_pv_p1_v1_zif_1 = [NS42_pv_int_1.P1v1[n,NS42_1.z_if_arg[n]] for n in range(0,NS42_1.t_len)]
NS42_pv_p1_v1_zif_2 = [NS42_pv_int_2.P1v1[n,NS42_2.z_if_arg[n]] for n in range(0,NS42_2.t_len)]
NS42_pv_p1_v1_zif_3 = [NS42_pv_int_3.P1v1[n,NS42_3.z_if_arg[n]] for n in range(0,NS42_3.t_len)]

NS42_pv_p1_v1_zif = NS42_pv_p1_v1_zif_1 + NS42_pv_p1_v1_zif_2 + NS42_pv_p1_v1_zif_3

S20_pv_p1_s1_mean_zig = [S20_pv_int.P1S1Mom1[n,S20.z_ig_arg[n]] for n in range(0,S20.t_len)]
S20_pv_p1_s1_mean_zif = [S20_pv_int.P1S1Mom1[n,S20.z_if_arg[n]] for n in range(0,S20.t_len)]
S20_pv_p1_w_mean_zif = [S20_pv_int.P1VMom1[n,S20.z_if_arg[n]] for n in range(0,S20.t_len)]
S20_pv_p1_s1_var_zig = [S20_pv_int.P1S1Mom2[n,S20.z_ig_arg[n]] for n in range(0,S20.t_len)]
S20_pv_p1_s1_var_zif = [S20_pv_int.P1S1Mom2[n,S20.z_if_arg[n]] for n in range(0,S20.t_len)]
S20_pv_p1_v_var_zif = [S20_pv_int.P1VMom2[n,S20.z_if_arg[n]] for n in range(0,S20.t_len)]
S20_pv_p1_v1_zif = [S20_pv_int.P1v1[n,S20.z_if_arg[n]] for n in range(0,S20.t_len)]

# Turbulent

NS42_pv_p2_s1_mean_zig_1 = [NS42_pv_int_1.P2S1Mom1[n,NS42_1.z_ig_arg[n]] for n in range(0,NS42_1.t_len)]
NS42_pv_p2_s1_mean_zig_2 = [NS42_pv_int_2.P2S1Mom1[n,NS42_2.z_ig_arg[n]] for n in range(0,NS42_2.t_len)]
NS42_pv_p2_s1_mean_zig_3 = [NS42_pv_int_3.P2S1Mom1[n,NS42_3.z_ig_arg[n]] for n in range(0,NS42_3.t_len)]

NS42_pv_p2_s1_mean_zig = NS42_pv_p2_s1_mean_zig_1 + NS42_pv_p2_s1_mean_zig_2 + NS42_pv_p2_s1_mean_zig_3

NS42_pv_p2_s1_mean_zif_1 = [NS42_pv_int_1.P2S1Mom1[n,NS42_1.z_if_arg[n]] for n in range(0,NS42_1.t_len)]
NS42_pv_p2_s1_mean_zif_2 = [NS42_pv_int_2.P2S1Mom1[n,NS42_2.z_if_arg[n]] for n in range(0,NS42_2.t_len)]
NS42_pv_p2_s1_mean_zif_3 = [NS42_pv_int_3.P2S1Mom1[n,NS42_3.z_if_arg[n]] for n in range(0,NS42_3.t_len)]

NS42_pv_p2_s1_mean_zif = NS42_pv_p2_s1_mean_zif_1 + NS42_pv_p2_s1_mean_zif_2 + NS42_pv_p2_s1_mean_zif_3

NS42_pv_p2_w_mean_zif_1 = [NS42_pv_int_1.P2VMom1[n,NS42_1.z_if_arg[n]] for n in range(0,NS42_1.t_len)]
NS42_pv_p2_w_mean_zif_2 = [NS42_pv_int_2.P2VMom1[n,NS42_2.z_if_arg[n]] for n in range(0,NS42_2.t_len)]
NS42_pv_p2_w_mean_zif_3 = [NS42_pv_int_3.P2VMom1[n,NS42_3.z_if_arg[n]] for n in range(0,NS42_3.t_len)]

NS42_pv_p2_w_mean_zif = NS42_pv_p2_w_mean_zif_1 + NS42_pv_p2_w_mean_zif_2 + NS42_pv_p2_w_mean_zif_3

NS42_pv_p2_s1_var_zig_1 = [NS42_pv_int_1.P2S1Mom2[n,NS42_1.z_ig_arg[n]] for n in range(0,NS42_1.t_len)]
NS42_pv_p2_s1_var_zig_2 = [NS42_pv_int_2.P2S1Mom2[n,NS42_2.z_ig_arg[n]] for n in range(0,NS42_2.t_len)]
NS42_pv_p2_s1_var_zig_3 = [NS42_pv_int_3.P2S1Mom2[n,NS42_3.z_ig_arg[n]] for n in range(0,NS42_3.t_len)]

NS42_pv_p2_s1_var_zig = NS42_pv_p2_s1_var_zig_1 + NS42_pv_p2_s1_var_zig_2 + NS42_pv_p2_s1_var_zig_3

NS42_pv_p2_s1_var_zif_1 = [NS42_pv_int_1.P2S1Mom2[n,NS42_1.z_if_arg[n]] for n in range(0,NS42_1.t_len)]
NS42_pv_p2_s1_var_zif_2 = [NS42_pv_int_2.P2S1Mom2[n,NS42_2.z_if_arg[n]] for n in range(0,NS42_2.t_len)]
NS42_pv_p2_s1_var_zif_3 = [NS42_pv_int_3.P2S1Mom2[n,NS42_3.z_if_arg[n]] for n in range(0,NS42_3.t_len)]

NS42_pv_p2_s1_var_zif = NS42_pv_p2_s1_var_zif_1 + NS42_pv_p2_s1_var_zif_2 + NS42_pv_p2_s1_var_zif_3

NS42_pv_p2_v_var_zif_1 = [NS42_pv_int_1.P2VMom2[n,NS42_1.z_if_arg[n]] for n in range(0,NS42_1.t_len)]
NS42_pv_p2_v_var_zif_2 = [NS42_pv_int_2.P2VMom2[n,NS42_2.z_if_arg[n]] for n in range(0,NS42_2.t_len)]
NS42_pv_p2_v_var_zif_3 = [NS42_pv_int_3.P2VMom2[n,NS42_3.z_if_arg[n]] for n in range(0,NS42_3.t_len)]

NS42_pv_p2_v_var_zif = NS42_pv_p2_v_var_zif_1 + NS42_pv_p2_v_var_zif_2 + NS42_pv_p2_v_var_zif_3

NS42_pv_p2_v1_zif_1 = [NS42_pv_int_1.P2v1[n,NS42_1.z_if_arg[n]] for n in range(0,NS42_1.t_len)]
NS42_pv_p2_v1_zif_2 = [NS42_pv_int_2.P2v1[n,NS42_2.z_if_arg[n]] for n in range(0,NS42_2.t_len)]
NS42_pv_p2_v1_zif_3 = [NS42_pv_int_3.P2v1[n,NS42_3.z_if_arg[n]] for n in range(0,NS42_3.t_len)]

NS42_pv_p2_v1_zif = NS42_pv_p2_v1_zif_1 + NS42_pv_p2_v1_zif_2 + NS42_pv_p2_v1_zif_3


S20_pv_p2_s1_mean_zig = [S20_pv_int.P2S1Mom1[n,S20.z_ig_arg[n]] for n in range(0,S20.t_len)]
S20_pv_p2_s1_mean_zif = [S20_pv_int.P2S1Mom1[n,S20.z_if_arg[n]] for n in range(0,S20.t_len)]
S20_pv_p2_w_mean_zif = [S20_pv_int.P2VMom1[n,S20.z_if_arg[n]] for n in range(0,S20.t_len)]
S20_pv_p2_s1_var_zig = [S20_pv_int.P2S1Mom2[n,S20.z_ig_arg[n]] for n in range(0,S20.t_len)]
S20_pv_p2_s1_var_zif = [S20_pv_int.P2S1Mom2[n,S20.z_if_arg[n]] for n in range(0,S20.t_len)]
S20_pv_p2_v_var_zif = [S20_pv_int.P2VMom2[n,S20.z_if_arg[n]] for n in range(0,S20.t_len)]
S20_pv_p2_v1_zif = [S20_pv_int.P2v1[n,S20.z_if_arg[n]] for n in range(0,S20.t_len)]

# time

time_1 = [NS42_1.z_enc[n]/L0 for n in range(0,NS42_1.t_len)]
time_2 = [NS42_2.z_enc[n]/L0 for n in range(0,NS42_2.t_len)]
time_3 = [NS42_3.z_enc[n]/L0 for n in range(0,NS42_3.t_len)]

time = time_1 + time_2 + time_3

z_if = np.concatenate([NS42_1.z_if,NS42_2.z_if,NS42_3.z_if])
z_ig = np.concatenate([NS42_1.z_ig,NS42_2.z_ig,NS42_3.z_ig])
z_enc = np.concatenate([NS42_1.z_enc,NS42_2.z_enc,NS42_3.z_enc])
z_is = np.concatenate([NS42_1.z_is,NS42_2.z_is,NS42_3.z_is])
L_Oz_zif = np.concatenate([NS42_1.L_Oz_zif,NS42_2.L_Oz_zif,NS42_3.L_Oz_zif])
w_enc = (B0*z_enc)**(1./3.)

# L_Oz_param = L0*(0.23-0.85*(L0/z_enc))**0.5
# s1_mean_zig_turb_param = N**2*(1.184*z_enc+1.78*L_Oz_param)

delta_b = (N**2*z_ig - s1_mean_zig)/(s1_grad_zig - N**2)
b_delta = delta_b*s1_grad_zig
delta_b_S20 = (N**2*S20.z_ig - S20_s1_mean_zig)/(S20_s1_grad_zig - N**2)
b_delta_S20 = delta_b_S20*S20_s1_grad_zig

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
oranges = matplotlib.cm.get_cmap('Oranges')

f, (ax1,ax2) = plt.subplots(1,2,figsize=(10,5))
ax1.grid(True)
ax2.grid(True)
ax1.set_ylim(0,1.6)
ax2.set_ylim(0.4,1)
ax2.set_xlim(15,30)
ax1.plot(np.mean(NS42_vort_int_1.int2[-4:-1,:],axis=0),NS42_1.y/np.mean(NS42_1.z_enc[-4:-1]),c=blues(0.5),label=r'$Fr_0=0$')
ax1.plot(np.mean(S20_vort_int.int2[4:7,:],axis=0),S20.y/np.mean(S20.z_enc[4:7]),c=blues(0.9),label=r'$Fr_0=20$')
# ax2.plot(time[1:-1],runningmean(NS42_vort_turbareafrac_zig,1),ls='--')
# ax2.plot(S20.z_enc[1:-1]/L0,runningmean(S20_vort_turbareafrac_zig,1),ls='--')
ax2.plot(time[1:-1],runningmean(NS42_vort_turbareafrac_zif,1),c=blues(0.5))
ax2.plot(S20.z_enc[1:-1]/L0,runningmean(S20_vort_turbareafrac_zif,1),c=blues(0.9))
ax1.axhline(np.mean(NS42_1.z_if[-4:-1]/NS42_1.z_enc[-4:-1]),0,0.05,c=blues(0.5))
ax1.axhline(np.mean(S20.z_if[4:7]/S20.z_enc[4:7]),0,0.05,c=blues(0.9))
ax1.set_xlabel(r'$a_\mathrm{T}$')
ax1.set_ylabel(r'$z/z_\mathrm{enc}$')
ax1.set_title(r'(a)',fontsize=20,loc='left')
ax2.set_xlabel(r'$z_\mathrm{enc}/L_0$')
ax2.set_ylabel(r'$(a_\mathrm{T})_{z_{i,f}}$')
ax2.set_title(r'(b)',fontsize=20,loc='left')
ax1.legend(loc='lower left',fontsize=20)
plt.tight_layout()
plt.savefig(opath+'area_frac_height_time_S20_S0_zif_vort.pdf')
plt.show()

f, (ax1,ax2) = plt.subplots(1,2,figsize=(10,5))
ax1.grid(True)
ax2.grid(True)
ax1.set_ylim(0,1.6)
ax2.set_ylim(0.4,1)
ax1.plot(np.mean(NS42_pv_int_1.int2[-4:-1,:],axis=0),NS42_1.y/np.mean(NS42_1.z_enc[-4:-1]),c=oranges(0.5),label=r'$Fr_0=0$')
ax1.plot(np.mean(S20_pv_int.int2[4:7,:],axis=0),S20.y/np.mean(S20.z_enc[4:7]),c=oranges(0.9),label=r'$Fr_0=20$')
# ax2.plot(time[1:-1],runningmean(NS42_pv_turbareafrac_zig,1),'C2')
# ax2.plot(S20.z_enc[1:-1]/L0,runningmean(S20_pv_turbareafrac_zig,1),'C4')
ax2.plot(time[1:-1],runningmean(NS42_pv_turbareafrac_zif,1),c=oranges(0.5))
ax2.plot(S20.z_enc[1:-1]/L0,runningmean(S20_pv_turbareafrac_zif,1),c=oranges(0.9))
ax1.axhline(np.mean(NS42_1.z_if[-4:-1]/NS42_1.z_enc[-4:-1]),0,0.05,c=oranges(0.5))
ax1.axhline(np.mean(S20.z_if[4:7]/S20.z_enc[4:7]),0,0.05,c=oranges(0.9))
ax1.set_xlabel(r'$a_\mathrm{T}$')
ax1.set_ylabel(r'$z/z_\mathrm{enc}$')
ax1.set_title(r'(a)',fontsize=20,loc='left')
ax2.set_xlabel(r'$z_\mathrm{enc}/L_0$')
ax2.set_ylabel(r'$(a_\mathrm{T})_{z_{i,f}}$')
ax2.set_title(r'(b)',fontsize=20,loc='left')
ax1.legend(loc='lower left',fontsize=20)
plt.tight_layout()
plt.savefig(opath+'area_frac_height_time_S20_S0_zif_PV.pdf')
plt.show()


f, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(10,10))
ax1.grid(True)
ax2.grid(True)
ax3.grid(True)
ax4.grid(True)
ax1.set_ylim(1,1.4)
ax1.set_xlim(0.4,1.4)
ax2.set_ylim(1,1.2)
ax2.set_xlim(15,30)
ax3.set_ylim(1,1.4)
ax4.set_xlim(15,30)
ax4.set_ylim(-0.05,0.05)
ax1.plot(np.mean(NS42_vort_int_1.P2S1Mom1[-4:-1,:],axis=0)/np.mean(N**2*NS42_1.z_enc[-4:-1]),NS42_1.y/np.mean(NS42_1.z_enc[-4:-1]),label=r'$Fr_0=0$')
ax1.plot(np.mean(S20_vort_int.P2S1Mom1[-3:,:],axis=0)/np.mean(N**2*S20.z_enc[-3:]),S20.y/np.mean(S20.z_enc[-3:]),'C2',label=r'$Fr_0=20$')
ax1.plot(np.mean(NS42_vort_int_1.P1S1Mom1[-4:-1,:],axis=0)/np.mean(N**2*NS42_1.z_enc[-4:-1]),NS42_1.y/np.mean(NS42_1.z_enc[-4:-1]),'C0',ls='--')
ax1.plot(np.mean(S20_vort_int.P1S1Mom1[-3:,:],axis=0)/np.mean(N**2*S20.z_enc[-3:]),S20.y/np.mean(S20.z_enc[-3:]),'C2',ls='--')
ax2.plot(time[1:-1],runningmean(NS42_vort_p2_s1_mean_zif/(N**2*z_enc),1))
ax2.plot(S20.z_enc[1:-1]/L0,runningmean(S20_vort_p2_s1_mean_zif/(N**2*S20.z_enc),1),'C2')
ax2.plot(time[1:-1],runningmean(NS42_vort_p1_s1_mean_zif/(N**2*z_enc),1),'C0',ls='--')
ax2.plot(S20.z_enc[1:-1]/L0,runningmean(S20_vort_p1_s1_mean_zif/(N**2*S20.z_enc),1),'C2',ls='--')
ax3.plot(np.mean(NS42_vort_int_1.P2VMom1[-4:-1,:],axis=0)/(np.mean(NS42_1.z_enc[-4:-1])*B0)**(1./3.),NS42_1.y/np.mean(NS42_1.z_enc[-4:-1]),label=r'$Fr_0=0$')
ax3.plot(np.mean(S20_vort_int.P2VMom1[-3:,:],axis=0)/(np.mean(S20.z_enc[-3:])*B0)**(1./3.),S20.y/np.mean(S20.z_enc[-3:]),'C2',label=r'$Fr_0=20$')
ax3.plot(np.mean(NS42_vort_int_1.P1VMom1[-4:-1,:],axis=0)/(np.mean(NS42_1.z_enc[-4:-1])*B0)**(1./3.),NS42_1.y/np.mean(NS42_1.z_enc[-4:-1]),'C0',ls='--')
ax3.plot(np.mean(S20_vort_int.P1VMom1[-3:,:],axis=0)/(np.mean(S20.z_enc[-3:])*B0)**(1./3.),S20.y/np.mean(S20.z_enc[-3:]),'C2',ls='--')
ax4.plot(time[1:-1],runningmean(NS42_vort_p2_w_mean_zif/(z_enc*B0)**(1./3.),1))
ax4.plot(S20.z_enc[1:-1]/L0,runningmean(S20_vort_p2_w_mean_zif/(S20.z_enc*B0)**(1./3.),1),'C2')
ax4.plot(time[1:-1],runningmean(NS42_vort_p1_w_mean_zif/(z_enc*B0)**(1./3.),1),'C0',ls='--')
ax4.plot(S20.z_enc[1:-1]/L0,runningmean(S20_vort_p1_w_mean_zif/(S20.z_enc*B0)**(1./3.),1),'C2',ls='--')
ax1.axhline(np.mean(NS42_1.z_if[-4:-1]/NS42_1.z_enc[-4:-1]),0,0.05,c='C0')
ax1.axhline(np.mean(S20.z_if[-3:]/S20.z_enc[-3:]),0,0.05,c='C2')
ax3.axhline(np.mean(NS42_1.z_if[-4:-1]/NS42_1.z_enc[-4:-1]),0,0.05,c='C0')
ax3.axhline(np.mean(S20.z_if[-3:]/S20.z_enc[-3:]),0,0.05,c='C2')
ax1.set_xlabel(r'$\langle b \rangle / b_\mathrm{enc}$')
ax1.set_ylabel(r'$z/z_\mathrm{enc}$')
ax1.set_title(r'(a)',fontsize=20,loc='left')
ax2.set_xlabel(r'$z_\mathrm{enc}/L_0$')
ax2.set_ylabel(r'$(\langle b \rangle)_{z_{i,f}}/b_\mathrm{enc}$')
ax2.set_title(r'(b)',fontsize=20,loc='left')
ax3.set_xlabel(r'$\langle w \rangle / w_*$')
ax3.set_ylabel(r'$z/z_\mathrm{enc}$')
ax3.set_title(r'(c)',fontsize=20,loc='left')
ax4.set_xlabel(r'$z_\mathrm{enc}/L_0$')
ax4.set_ylabel(r'$(\langle w \rangle)_{z_{i,f}}/w_*$')
ax4.set_title(r'(d)',fontsize=20,loc='left')
ax1.legend(loc='best',fontsize=20)
plt.tight_layout()
plt.savefig(opath+'s1_w_mean_height_time_S20_S0_vort.pdf')
plt.show()

f, (ax1,ax2) = plt.subplots(1,2,figsize=(10,5))
ax1.grid(True)
ax2.grid(True)
ax1.set_xlim(0.8,1.4)
ax2.set_xlim(-0.15,0.15)
ax1.set_ylim(1,1.4)
ax2.set_ylim(1,1.4)
ax1.plot(np.mean(NS42_vort_int_1.P2S1Mom1[-4:-1,:],axis=0)/np.mean(N**2*NS42_1.z_enc[-4:-1]),NS42_1.y/np.mean(NS42_1.z_enc[-4:-1]),c=blues(0.5),label=r'$Fr_0=0$')
ax1.plot(np.mean(S20_vort_int.P2S1Mom1[4:7,:],axis=0)/np.mean(N**2*S20.z_enc[4:7]),S20.y/np.mean(S20.z_enc[4:7]),c=blues(0.9),label=r'$Fr_0=20$')
ax1.plot(np.mean(NS42_vort_int_1.P1S1Mom1[-4:-1,:],axis=0)/np.mean(N**2*NS42_1.z_enc[-4:-1]),NS42_1.y/np.mean(NS42_1.z_enc[-4:-1]),c=blues(0.5),ls='--')
ax1.plot(np.mean(S20_vort_int.P1S1Mom1[4:7,:],axis=0)/np.mean(N**2*S20.z_enc[4:7]),S20.y/np.mean(S20.z_enc[4:7]),c=blues(0.9),ls='--')
ax1.plot(S20.y/np.mean(S20.z_enc[4:7]),S20.y/np.mean(S20.z_enc[4:7]),'k--',label=r'$N^2z$')
ax2.plot(np.mean(NS42_vort_int_1.P2VMom1[-4:-1,:],axis=0)/(np.mean(NS42_1.z_enc[-4:-1])*B0)**(1./3.),NS42_1.y/np.mean(NS42_1.z_enc[-4:-1]),c=blues(0.5),label=r'$Fr_0=0$')
ax2.plot(np.mean(S20_vort_int.P2VMom1[4:7,:],axis=0)/(np.mean(S20.z_enc[4:7])*B0)**(1./3.),S20.y/np.mean(S20.z_enc[4:7]),c=blues(0.9),label=r'$Fr_0=20$')
ax2.plot(np.mean(NS42_vort_int_1.P1VMom1[-4:-1,:],axis=0)/(np.mean(NS42_1.z_enc[-4:-1])*B0)**(1./3.),NS42_1.y/np.mean(NS42_1.z_enc[-4:-1]),c=blues(0.5),ls='--')
ax2.plot(np.mean(S20_vort_int.P1VMom1[4:7,:],axis=0)/(np.mean(S20.z_enc[4:7])*B0)**(1./3.),S20.y/np.mean(S20.z_enc[4:7]),c=blues(0.9),ls='--')
#ax2.plot(np.mean(S20_vort_int.int2[4:7,]*S20_vort_int.P2VMom1[4:7,:],axis=0)+np.mean(S20_vort_int.int1[4:7,]*S20_vort_int.P1VMom1[4:7,:],axis=0),S20.y/np.mean(S20.z_enc[4:7]),'k')
ax1.axhline(np.mean(NS42_1.z_if[-4:-1]/NS42_1.z_enc[-4:-1]),0,0.05,c=blues(0.5))
ax1.axhline(np.mean(S20.z_if[4:7]/S20.z_enc[4:7]),0,0.05,c=blues(0.7))
ax2.axhline(np.mean(NS42_1.z_if[-4:-1]/NS42_1.z_enc[-4:-1]),0,0.05,c=blues(0.5))
ax2.axhline(np.mean(S20.z_if[4:7]/S20.z_enc[4:7]),0,0.05,c=blues(0.7))
ax1.set_xlabel(r'$\langle b \rangle / b_\mathrm{enc}$')
ax1.set_ylabel(r'$z/z_\mathrm{enc}$')
ax1.set_title(r'(a)',fontsize=20,loc='left')
ax2.set_xlabel(r'$\langle w \rangle / w_*$')
ax2.set_ylabel(r'$z/z_\mathrm{enc}$')
ax2.set_title(r'(b)',fontsize=20,loc='left')
ax1.legend(loc='best',fontsize=20)
plt.tight_layout()
plt.savefig(opath+'s1_w_mean_height_S20_S0_vort.pdf')
plt.show()


f, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(10,10))
ax1.grid(True)
ax2.grid(True)
ax3.grid(True)
ax4.grid(True)
ax1.set_ylim(1,1.4)
ax2.set_xlim(15,30)
ax2.set_ylim(0,2)
ax3.set_ylim(1,1.4)
ax3.set_xlim(0,0.6)
ax4.set_ylim(0,0.6)
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
ax1.set_title('(a)',fontsize=20,loc='left')
ax2.set_xlabel(r'$z_\mathrm{enc}/L_0$')
ax2.set_ylabel(r'$((b_\mathrm{rms})_\mathrm{T})_{z_{i,f}}/(N^2L_0)$')
ax2.set_title('(b)',fontsize=20,loc='left')
ax3.set_xlabel(r'$(w_\mathrm{rms})_\mathrm{T} / w_*$')
ax3.set_ylabel(r'$z/z_\mathrm{enc}$')
ax3.set_title('(c)',fontsize=20,loc='left')
ax4.set_xlabel(r'$z_\mathrm{enc}/L_0$')
ax4.set_ylabel(r'$((w_\mathrm{rms})_\mathrm{T})_{z_{i,f}}/w_*$')
ax4.set_title('(d)',fontsize=20,loc='left')
ax1.legend(loc='best',fontsize=20)
plt.tight_layout()
plt.savefig(opath+'s1_w_rms_height_time_S20_S0_vort_zif.pdf')
plt.show()

plt.figure(figsize=(5,5))
plt.grid(True)
plt.xlim(0,1)
plt.ylim(1,1.4)
plt.plot(np.mean((NS42_vort_int_1.int2[-4:-1,:]*NS42_vort_int_1.P2S1Mom2[-4:-1,:])/NS42_1.r2S[-4:-1,:],axis=0),NS42_1.y/np.mean(NS42_1.z_enc[-4:-1]),'C0')
plt.plot(np.mean((NS42_vort_int_1.int1[-4:-1,:]*NS42_vort_int_1.P1S1Mom2[-4:-1,:])/NS42_1.r2S[-4:-1,:],axis=0),NS42_1.y/np.mean(NS42_1.z_enc[-4:-1]),'C0',ls='--')
plt.plot(np.mean((NS42_vort_int_1.int2[-4:-1,:]*NS42_vort_int_1.int1[-4:-1,:]*(NS42_vort_int_1.P2S1Mom1[-4:-1,:]-NS42_vort_int_1.P1S1Mom1[-4:-1,:])**2)/NS42_1.r2S[-4:-1,:],axis=0),NS42_1.y/np.mean(NS42_1.z_enc[-4:-1]),'C0',ls='-.')
#plt.plot((NS42_vort_int_1.int2[-2,:]*NS42_vort_int_1.P2S1Mom2[-2,:]+NS42_vort_int_1.int1[-2,:]*NS42_vort_int_1.P1S1Mom2[-2,:]+NS42_vort_int_1.int2[-2,:]*NS42_vort_int_1.int1[-2,:]*(NS42_vort_int_1.P2S1Mom1[-2,:]-NS42_vort_int_1.P1S1Mom1[-2,:])**2)/NS42_1.r2S[-2,:],NS42_1.y/NS42_1.z_enc[-2])
# plt.plot((NS42_pv_int_1.int2[-2,:]*NS42_pv_int_1.P2S1Mom2[-2,:])/NS42_1.r2S[-2,:],NS42_1.y/NS42_1.z_enc[-2],'C1')
# plt.plot((NS42_pv_int_1.int1[-2,:]*NS42_pv_int_1.P1S1Mom2[-2,:])/NS42_1.r2S[-2,:],NS42_1.y/NS42_1.z_enc[-2],'C1',ls='--')
# plt.plot((NS42_pv_int_1.int2[-2,:]*NS42_pv_int_1.int1[-2,:]*(NS42_pv_int_1.P2S1Mom1[-2,:]-NS42_pv_int_1.P1S1Mom1[-2,:])**2)/NS42_1.r2S[-2,:],NS42_1.y/NS42_1.z_enc[-2],'C1',ls='-.')
plt.plot(np.mean((S20_vort_int.int2[-3:,:]*S20_vort_int.P2S1Mom2[-3:,:])/S20.r2S[-3:,:],axis=0),S20.y/np.mean(S20.z_enc[-3:]),'C2')
plt.plot(np.mean((S20_vort_int.int1[-3:,:]*S20_vort_int.P1S1Mom2[-3:,:])/S20.r2S[-3:,:],axis=0),S20.y/np.mean(S20.z_enc[-3:]),'C2',ls='--')
plt.plot(np.mean((S20_vort_int.int2[-3:,:]*S20_vort_int.int1[-3:,:]*(S20_vort_int.P2S1Mom1[-3:,:]-S20_vort_int.P1S1Mom1[-3:,:])**2)/S20.r2S[-3:,:],axis=0),S20.y/np.mean(S20.z_enc[-3:]),'C2',ls='-.')
# plt.plot((S20_pv_int.int2[-1,:]*S20_pv_int.P2S1Mom2[-1,:])/S20.r2S[-1,:],S20.y/S20.z_enc[-1],'C3')
# plt.plot((S20_pv_int.int1[-1,:]*S20_pv_int.P1S1Mom2[-1,:])/S20.r2S[-1,:],S20.y/S20.z_enc[-1],'C3',ls='--')
# plt.plot((S20_pv_int.int2[-1,:]*S20_pv_int.int1[-1,:]*(S20_pv_int.P2S1Mom1[-1,:]-S20_pv_int.P1S1Mom1[-1,:])**2)/S20.r2S[-1,:],S20.y/S20.z_enc[-1],'C3',ls='-.')
plt.xlabel(r'$f/\langle b^{\prime ^2} \rangle$')
plt.ylabel(r'$z/z_\mathrm{enc}$')
plt.title(r'$z_\mathrm{enc}/L_0=20$',fontsize=20,loc='left')
plt.tight_layout()
plt.savefig(opath+'s1_var_contributions_S20_S0_vort.pdf')
plt.show()


f, (ax1,ax2) = plt.subplots(1,2,figsize=(10,5))
ax1.grid(True)
ax2.grid(True)
ax1.set_ylim(1,1.4)
ax1.set_xlim(0,0.6)
ax2.set_ylim(0,0.6)
ax1.plot(np.mean(np.sqrt(NS42_vort_int_1.P2VMom2[-4:-1,:]),axis=0)/np.mean(w_enc[4:7]),NS42_1.y/np.mean(NS42_1.z_enc[-4:-1]),c=blues(0.5))
ax1.plot(np.mean(np.sqrt(S20_vort_int.P2VMom2[4:7,:]),axis=0)/np.mean((B0*S20.z_enc[4:7])**(1./3.)),S20.y/np.mean(S20.z_enc[4:7]),c=blues(0.9))
ax2.plot(time[1:-1],runningmean(np.sqrt(NS42_vort_p2_w_var_zif)/w_enc,1),c=blues(0.5),label=r'$Fr_0=0$')
ax2.plot(S20.z_enc[1:-1]/L0,runningmean(np.sqrt(S20_vort_p2_w_var_zif)/(B0*S20.z_enc)**(1./3.),1),c=blues(0.9),label=r'$Fr_0=20$')
ax1.axhline(np.mean(NS42_1.z_if[-4:-1]/NS42_1.z_enc[-4:-1]),0,0.05,c=blues(0.5))
ax1.axhline(np.mean(S20.z_if[4:7]/S20.z_enc[4:7]),0,0.05,c=blues(0.9))
ax1.set_xlabel(r'$(w_\mathrm{rms})_\mathrm{T} / w_*$')
ax1.set_ylabel(r'$z/z_\mathrm{enc}$')
ax2.set_xlabel(r'$z_\mathrm{enc}/L_0$')
ax2.set_ylabel(r'$((w_\mathrm{rms})_\mathrm{T})_{z_{i,f}}/w_*$')
ax2.legend(loc='lower left',fontsize=20)
plt.tight_layout()
plt.savefig(opath+'w_rms_height_time_S20_S0_zif.pdf')
plt.show()


f, (ax1,ax2) = plt.subplots(1,2,figsize=(10,5))
ax1.grid(True)
ax2.grid(True)
ax1.set_ylim(1,1.4)
ax1.set_xlim(-0.4,0.2)
ax2.set_ylim(0,0.4)
ax2.set_xlim(15,30)
ax1.plot(np.mean(NS42_vort_int_1.P2v1[-4:-1,:]-NS42_vort_int_1.P2S1Mom1[-4:-1,:]*NS42_vort_int_1.P2VMom1[-4:-1,:],axis=0)/B0,NS42_1.y/np.mean(NS42_1.z_enc[-4:-1]),c=blues(0.5))
ax1.plot(np.mean(S20_vort_int.P2v1[4:7,:]-S20_vort_int.P2S1Mom1[4:7,:]*S20_vort_int.P2VMom1[4:7,:],axis=0)/B0,S20.y/np.mean(S20.z_enc[4:7]),c=blues(0.9))
ax2.plot(time[1:-1],-runningmean(np.array(NS42_vort_p2_v1_zif)-np.array(NS42_vort_p2_s1_mean_zif)*np.array(NS42_vort_p2_w_mean_zif),1)/B0,c=blues(0.5),label=r'$Fr_0=0$')
ax2.plot(S20.z_enc[1:-1]/L0,-runningmean(np.array(S20_vort_p2_v1_zif)-np.array(S20_vort_p2_s1_mean_zif)*np.array(S20_vort_p2_w_mean_zif),1)/B0,c=blues(0.9),label=r'$Fr_0=20$')
ax1.axhline(np.mean(NS42_1.z_if[-4:-1]/NS42_1.z_enc[-4:-1]),0,0.05,c=blues(0.5))
ax1.axhline(np.mean(S20.z_if[4:7]/S20.z_enc[4:7]),0,0.05,c=blues(0.9))
ax1.set_xlabel(r'$\langle b^\prime w^\prime\rangle_\mathrm{T}/B_0$')
ax1.set_ylabel(r'$z/z_\mathrm{enc}$')
ax1.set_title(r'(a)',fontsize=20,loc='left')
ax2.set_xlabel(r'$z_\mathrm{enc}/L_0$')
ax2.set_ylabel(r'$-(\langle b^\prime w^\prime \rangle_\mathrm{T})_{z_{i,f}}/B_0$')
ax2.set_title(r'(b)',fontsize=20,loc='left')
ax2.legend(loc='best',fontsize=20)
plt.tight_layout()
plt.savefig(opath+'s1_vflux_turb_height_time_S20_S0_vort.pdf')
plt.show()

f, (ax1,ax2) = plt.subplots(1,2,figsize=(10,5))
ax1.grid(True)
ax2.grid(True)
ax1.set_ylim(1,1.4)
ax1.set_xlim(-0.4,0.2)
ax2.set_ylim(0,0.4)
ax2.set_xlim(15,30)
ax1.plot(np.mean(NS42_pv_int_1.P2v1[-4:-1,:]-NS42_pv_int_1.P2S1Mom1[-4:-1,:]*NS42_pv_int_1.P2VMom1[-4:-1,:],axis=0)/B0,NS42_1.y/np.mean(NS42_1.z_enc[-4:-1]),c=oranges(0.5))
ax1.plot(np.mean(S20_pv_int.P2v1[4:7,:]-S20_pv_int.P2S1Mom1[4:7,:]*S20_pv_int.P2VMom1[4:7,:],axis=0)/B0,S20.y/np.mean(S20.z_enc[4:7]),c=oranges(0.9))
ax2.plot(time[1:-1],-runningmean(np.array(NS42_pv_p2_v1_zif)-np.array(NS42_pv_p2_s1_mean_zif)*np.array(NS42_pv_p2_w_mean_zif),1)/B0,c=oranges(0.5),label=r'$Fr_0=0$')
ax2.plot(S20.z_enc[1:-1]/L0,-runningmean(np.array(S20_pv_p2_v1_zif)-np.array(S20_pv_p2_s1_mean_zif)*np.array(S20_pv_p2_w_mean_zif),1)/B0,c=oranges(0.9),label=r'$Fr_0=20$')
ax1.axhline(np.mean(NS42_1.z_if[-4:-1]/NS42_1.z_enc[-4:-1]),0,0.05,c=oranges(0.5))
ax1.axhline(np.mean(S20.z_if[4:7]/S20.z_enc[4:7]),0,0.05,c=oranges(0.9))
ax1.set_xlabel(r'$\langle b^\prime w^\prime\rangle_\mathrm{T}/B_0$')
ax1.set_ylabel(r'$z/z_\mathrm{enc}$')
ax1.set_title(r'(a)',fontsize=20,loc='left')
ax2.set_xlabel(r'$z_\mathrm{enc}/L_0$')
ax2.set_ylabel(r'$-(\langle b^\prime w^\prime \rangle_\mathrm{T})_{z_{i,f}}/B_0$')
ax2.set_title(r'(b)',fontsize=20,loc='left')
ax2.legend(loc='best',fontsize=20)
plt.tight_layout()
plt.savefig(opath+'s1_vflux_turb_height_time_S20_S0_pv.pdf')
plt.show()

f, (ax1,ax2) = plt.subplots(1,2,figsize=(10,5))
ax1.grid(True)
ax2.grid(True)
ax1.set_ylim(0,1.4)
#ax1.set_xlim(-0.5,0.5)
ax2.set_ylim(0,0.2)
ax2.set_xlim(15,30)
ax1.plot(np.mean(NS42_1.Rsv[-4:-1,:]/B0,axis=0),NS42_1.y/np.mean(NS42_1.z_enc[-4:-1]))
ax1.plot(np.mean(S20.Rsv[-3:,:]/B0,axis=0),S20.y/np.mean(S20.z_enc[-3:]),'C2')
ax2.plot(time[1:-1],-runningmean(np.array(s1_flux_zif)/B0,1),'C0',label=r'$Fr_0=0$')
ax2.plot(S20.z_enc[1:-1]/L0,-runningmean(np.array(S20_s1_flux_zif)/B0,1),'C2',label=r'$Fr_0=20$')
ax1.set_xlabel(r'$\langle b^\prime w^\prime\rangle/B_0$')
ax1.set_ylabel(r'$z/z_\mathrm{enc}$')
ax1.set_title('(a)',fontsize=20,loc='left')
ax2.set_xlabel(r'$z_\mathrm{enc}/L_0$')
ax2.set_ylabel(r'$-(\langle b^\prime w^\prime \rangle)_{z_{i,f}}/B_0$')
ax2.set_title('(b)',fontsize=20,loc='left')
ax2.legend(loc='lower left',fontsize=20)
plt.tight_layout()
plt.savefig(opath+'s1_vflux_total_height_time_S20_S0.pdf')
plt.show()


f, (ax1,ax2) = plt.subplots(1,2,figsize=(10,5))
ax1.grid(True)
ax2.grid(True)
ax1.set_ylim(1,1.4)
ax1.set_xlim(-0.2,0)
ax2.set_ylim(0,0.2)
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
ax2.set_title('(b)',fontsize=20,loc='left')
ax2.legend(loc='best',fontsize=20)
plt.tight_layout()
plt.savefig(opath+'rho_bw_turb_height_time_S20_S0.pdf')
plt.show()

f, (ax1,ax2) = plt.subplots(1,2,figsize=(10,5))
ax1.grid(True)
ax2.grid(True)
ax1.set_ylim(-0.2,0.2)
ax2.set_ylim(-0.2,0.2)
ax1.plot(time[1:-1],runningmean((np.array(NS42_vort_p2_s1_mean_zif)-np.array(NS42_vort_p1_s1_mean_zif))/(z_enc*N**2),1),label=r'$Fr_0=0$')
ax1.plot(S20.z_enc[1:-1]/L0,runningmean((np.array(S20_vort_p2_s1_mean_zif)-np.array(S20_vort_p1_s1_mean_zif))/(S20.z_enc*N**2),1),'C2',label=r'$Fr_0=20$')
ax2.plot(time[1:-1],runningmean((np.array(NS42_vort_p2_w_mean_zif)-np.array(NS42_vort_p1_w_mean_zif))/(z_enc*B0)**(1./3.),1))
ax2.plot(S20.z_enc[1:-1]/L0,runningmean((np.array(S20_vort_p2_w_mean_zif)-np.array(S20_vort_p1_w_mean_zif))/(S20.z_enc*B0)**(1./3.),1),'C2')
ax1.set_xlabel(r'$z_\mathrm{enc}/L_0$')
ax1.set_ylabel(r'$(\langle b \rangle_\mathrm{T}-\langle b \rangle_\mathrm{NT})/b_\mathrm{enc}$')
ax1.set_title(r'(a)',fontsize=20,loc='left')
ax2.set_xlabel(r'$z_\mathrm{enc}/L_0$')
ax2.set_ylabel(r'$(\langle w \rangle_\mathrm{T}-\langle w \rangle_\mathrm{NT})/w_* $')
ax2.set_title(r'(b)',fontsize=20,loc='left')
ax1.legend(loc='best',fontsize=20)
plt.tight_layout()
plt.savefig(opath+'s1_w_meandiff_zif_time_S20_S0.pdf')
plt.show()

plt.figure(figsize=(5,5))
plt.grid(True)
plt.xlim(-0.2,0.05)
plt.ylim(1,1.4)
plt.plot(np.mean(NS42_vort_int_1.int2[-4:-1,:]*(NS42_vort_int_1.P2v1[-4:-1,:]-NS42_vort_int_1.P2S1Mom1[-4:-1,:]*NS42_vort_int_1.P2VMom1[-4:-1,:])/B0,axis=0),NS42_1.y/np.mean(NS42_1.z_enc[-4:-1]),c=blues(0.5),label=r'$Fr_0=0$')
plt.plot(np.mean(S20_vort_int.int2[4:7,:]*(S20_vort_int.P2v1[4:7,:]-S20_vort_int.P2S1Mom1[4:7,:]*S20_vort_int.P2VMom1[4:7,:])/B0,axis=0),S20.y/np.mean(S20.z_enc[4:7]),c=blues(0.9),label=r'$Fr_0=20$')
plt.plot(np.mean(NS42_vort_int_1.int1[-4:-1,:]*(NS42_vort_int_1.P1v1[-4:-1,:]-NS42_vort_int_1.P1S1Mom1[-4:-1,:]*NS42_vort_int_1.P1VMom1[-4:-1,:])/B0,axis=0),NS42_1.y/np.mean(NS42_1.z_enc[-4:-1]),c=blues(0.5),ls='--')
plt.plot(np.mean(S20_vort_int.int1[4:7,:]*(S20_vort_int.P1v1[4:7,:]-S20_vort_int.P1S1Mom1[4:7,:]*S20_vort_int.P1VMom1[4:7,:])/B0,axis=0),S20.y/np.mean(S20.z_enc[4:7]),c=blues(0.9),ls='--')
plt.plot(np.mean(NS42_vort_int_1.int2[-4:-1,:]*NS42_vort_int_1.int1[-4:-1,:]*(NS42_vort_int_1.P2S1Mom1[-4:-1,:]-NS42_vort_int_1.P1S1Mom1[-4:-1,:])*(NS42_vort_int_1.P2VMom1[-4:-1,:]-NS42_vort_int_1.P1VMom1[-4:-1,:])/B0,axis=0),NS42_1.y/np.mean(NS42_1.z_enc[-4:-1]),c=blues(0.5),ls='-.')
plt.plot(np.mean(S20_vort_int.int2[4:7,:]*S20_vort_int.int1[4:7,:]*(S20_vort_int.P2S1Mom1[4:7,:]-S20_vort_int.P1S1Mom1[4:7,:])*(S20_vort_int.P2VMom1[4:7,:]-S20_vort_int.P1VMom1[4:7,:])/B0,axis=0),S20.y/np.mean(S20.z_enc[4:7]),c=blues(0.9),ls='-.')
# plt.plot(np.mean((NS42_vort_int_1.int2[-4:-1,:]*(NS42_vort_int_1.P2v1[-4:-1,:]-NS42_vort_int_1.P2S1Mom1[-4:-1,:]*NS42_vort_int_1.P2VMom1[-4:-1,:])+NS42_vort_int_1.int1[-4:-1,:]*(NS42_vort_int_1.P1v1[-4:-1,:]-NS42_vort_int_1.P1S1Mom1[-4:-1,:]*NS42_vort_int_1.P1VMom1[-4:-1,:])+NS42_vort_int_1.int2[-4:-1,:]*NS42_vort_int_1.int1[-4:-1,:]*(NS42_vort_int_1.P2S1Mom1[-4:-1,:]-NS42_vort_int_1.P1S1Mom1[-4:-1,:])*(NS42_vort_int_1.P2VMom1[-4:-1,:]-NS42_vort_int_1.P1VMom1[-4:-1,:]))/B0,axis=0),NS42_1.y/np.mean(NS42_1.z_enc[-4:-1]),'C1')
# plt.plot(np.mean(NS42_1.Rsv[-4:-1,:]/B0,axis=0),NS42_1.y/np.mean(NS42_1.z_enc[-4:-1]),'k--')
# plt.plot(np.mean((S20_vort_int.int2[4:7,:]*(S20_vort_int.P2v1[4:7,:]-S20_vort_int.P2S1Mom1[4:7,:]*S20_vort_int.P2VMom1[4:7,:])+S20_vort_int.int1[4:7,:]*(S20_vort_int.P1v1[4:7,:]-S20_vort_int.P1S1Mom1[4:7,:]*S20_vort_int.P1VMom1[4:7,:])+S20_vort_int.int2[4:7,:]*S20_vort_int.int1[4:7,:]*(S20_vort_int.P2S1Mom1[4:7,:]-S20_vort_int.P1S1Mom1[4:7,:])*(S20_vort_int.P2VMom1[4:7,:]-S20_vort_int.P1VMom1[4:7,:]))/B0,axis=0),S20.y/np.mean(S20.z_enc[4:7]),'C3')
# plt.plot(np.mean(S20.Rsv[4:7,:]/B0,axis=0),S20.y/np.mean(S20.z_enc[4:7]),'k--')
plt.xlabel(r'$f /B_0$')
plt.ylabel(r'$z/z_\mathrm{enc}$')
plt.legend(loc='best',fontsize=20)
plt.tight_layout()
plt.savefig(opath+'s1_vflux_contributions_S20_S0.pdf')
plt.show()

plt.figure(figsize=(5,5))
plt.grid(True)
plt.xlim(-0.2,0.05)
plt.ylim(1,1.4)
plt.plot(np.mean(NS42_pv_int_1.int2[-4:-1,:]*(NS42_pv_int_1.P2v1[-4:-1,:]-NS42_pv_int_1.P2S1Mom1[-4:-1,:]*NS42_pv_int_1.P2VMom1[-4:-1,:])/B0,axis=0),NS42_1.y/np.mean(NS42_1.z_enc[-4:-1]),c=oranges(0.5),label=r'$Fr_0=0$')
plt.plot(np.mean(S20_pv_int.int2[4:7,:]*(S20_pv_int.P2v1[4:7,:]-S20_pv_int.P2S1Mom1[4:7,:]*S20_pv_int.P2VMom1[4:7,:])/B0,axis=0),S20.y/np.mean(S20.z_enc[4:7]),c=oranges(0.9),label=r'$Fr_0=20$')
plt.plot(np.mean(NS42_pv_int_1.int1[-4:-1,:]*(NS42_pv_int_1.P1v1[-4:-1,:]-NS42_pv_int_1.P1S1Mom1[-4:-1,:]*NS42_pv_int_1.P1VMom1[-4:-1,:])/B0,axis=0),NS42_1.y/np.mean(NS42_1.z_enc[-4:-1]),c=oranges(0.5),ls='--')
plt.plot(np.mean(S20_pv_int.int1[4:7,:]*(S20_pv_int.P1v1[4:7,:]-S20_pv_int.P1S1Mom1[4:7,:]*S20_pv_int.P1VMom1[4:7,:])/B0,axis=0),S20.y/np.mean(S20.z_enc[4:7]),c=oranges(0.9),ls='--')
plt.plot(np.mean(NS42_pv_int_1.int2[-4:-1,:]*NS42_pv_int_1.int1[-4:-1,:]*(NS42_pv_int_1.P2S1Mom1[-4:-1,:]-NS42_pv_int_1.P1S1Mom1[-4:-1,:])*(NS42_pv_int_1.P2VMom1[-4:-1,:]-NS42_pv_int_1.P1VMom1[-4:-1,:])/B0,axis=0),NS42_1.y/np.mean(NS42_1.z_enc[-4:-1]),c=oranges(0.5),ls='-.')
plt.plot(np.mean(S20_pv_int.int2[4:7,:]*S20_pv_int.int1[4:7,:]*(S20_pv_int.P2S1Mom1[4:7,:]-S20_pv_int.P1S1Mom1[4:7,:])*(S20_pv_int.P2VMom1[4:7,:]-S20_pv_int.P1VMom1[4:7,:])/B0,axis=0),S20.y/np.mean(S20.z_enc[4:7]),c=oranges(0.9),ls='-.')
# plt.plot(np.mean((NS42_pv_int_1.int2[-4:-1,:]*(NS42_pv_int_1.P2v1[-4:-1,:]-NS42_pv_int_1.P2S1Mom1[-4:-1,:]*NS42_pv_int_1.P2VMom1[-4:-1,:])+NS42_pv_int_1.int1[-4:-1,:]*(NS42_pv_int_1.P1v1[-4:-1,:]-NS42_pv_int_1.P1S1Mom1[-4:-1,:]*NS42_pv_int_1.P1VMom1[-4:-1,:])+NS42_pv_int_1.int2[-4:-1,:]*NS42_pv_int_1.int1[-4:-1,:]*(NS42_pv_int_1.P2S1Mom1[-4:-1,:]-NS42_pv_int_1.P1S1Mom1[-4:-1,:])*(NS42_pv_int_1.P2VMom1[-4:-1,:]-NS42_pv_int_1.P1VMom1[-4:-1,:]))/B0,axis=0),NS42_1.y/np.mean(NS42_1.z_enc[-4:-1]),c=oranges(0.5))
# plt.plot(np.mean(NS42_1.Rsv[-4:-1,:]/B0,axis=0),NS42_1.y/np.mean(NS42_1.z_enc[-4:-1]),'k--')
# plt.plot(np.mean((S20_pv_int.int2[4:7,:]*(S20_pv_int.P2v1[4:7,:]-S20_pv_int.P2S1Mom1[4:7,:]*S20_pv_int.P2VMom1[4:7,:])+S20_pv_int.int1[4:7,:]*(S20_pv_int.P1v1[4:7,:]-S20_pv_int.P1S1Mom1[4:7,:]*S20_pv_int.P1VMom1[4:7,:])+S20_pv_int.int2[4:7,:]*S20_pv_int.int1[4:7,:]*(S20_pv_int.P2S1Mom1[4:7,:]-S20_pv_int.P1S1Mom1[4:7,:])*(S20_pv_int.P2VMom1[4:7,:]-S20_pv_int.P1VMom1[4:7,:]))/B0,axis=0),S20.y/np.mean(S20.z_enc[4:7]),c=oranges(0.9))
# plt.plot(np.mean(S20.Rsv[4:7,:]/B0,axis=0),S20.y/np.mean(S20.z_enc[4:7]),'k--')
plt.xlabel(r'$f /B_0$')
plt.ylabel(r'$z/z_\mathrm{enc}$')
plt.legend(loc='best',fontsize=20)
plt.tight_layout()
plt.savefig(opath+'s1_vflux_contributions_S20_S0_pv.pdf')
plt.show()


plt.figure(figsize=(5,5))
plt.grid(True)
#plt.ylim(0,0.5)
plt.plot(time,L_Oz_zif/delta_b,label=r'$(L_{Oz})_{z_{i,f}}/\delta_b$')
plt.plot(time,L_Oz_zif/L0,label=r'$(L_{Oz})_{z_{i,f}}/L_0$')
plt.plot(time,delta_b/L0,label=r'$\delta_b/L_0$')
plt.xlabel(r'$z_\mathrm{enc}/L_0$')
#plt.ylabel(r'$(L_{Oz})_{z_{i,f}}/\delta_b$')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig(opath+'L_Oz_delta_b.pdf')
plt.show()

f, (ax1,ax2) = plt.subplots(1,2,sharey='row',figsize=(10,5))
ax1.grid(True)
ax2.grid(True)
ax1.set_ylim(0.8,1.1)
ax1.set_xlim(15,30)
ax2.set_xlim(15,30)
ax1.plot(time[1:-1],runningmean(NS42_vort_p2_s1_mean_zig/(N**2*z_if),1),label=r'$Fr_0=0$')
#ax1.plot(time,NS42_pv_p2_s1_mean_zig/(N**2*z_if),label=r'$\Pi^2$, $Fr_0=0$')
ax1.plot(S20.z_enc[1:-1]/L0,runningmean(S20_vort_p2_s1_mean_zig/(N**2*S20.z_if),1),'C2',label=r'$Fr_0=20$')
#ax1.plot(S20.z_enc/L0,S20_pv_p2_s1_mean_zig/(N**2*S20.z_if),label=r'$\Pi^2$, $Fr_0=20$')
ax2.plot(time[1:-1],runningmean(NS42_vort_p1_s1_mean_zig/(N**2*z_ig),1))
#ax2.plot(time,NS42_pv_p1_s1_mean_zig/(N**2*z_ig))
ax2.plot(S20.z_enc[1:-1]/L0,runningmean(S20_vort_p1_s1_mean_zig/(N**2*S20.z_ig),1),'C2')
#ax2.plot(S20.z_enc/L0,S20_pv_p1_s1_mean_zig/(N**2*S20.z_ig))
ax1.set_xlabel(r'$z_\mathrm{enc}/L_0$')
ax2.set_xlabel(r'$z_\mathrm{enc}/L_0$')
ax1.set_ylabel(r'$(\langle b\rangle_\mathrm{T})_{z_{i,g}}/N^2z_{i,f}$')
ax2.set_ylabel(r'$(\langle b\rangle_\mathrm{NT})_{z_{i,g}}/N^2z_{i,g}$')
ax1.set_title('(a)',loc='left',fontsize=20)
ax2.set_title('(b)',loc='left',fontsize=20)
ax1.legend(loc='best',fontsize=20)
plt.tight_layout()
plt.savefig(opath+'s1_mean_zig_scaled_vort.pdf')
plt.show()

plt.figure(figsize=(5,5))
plt.grid(True)
plt.ylim(0.8,1)
plt.plot(time,NS42_vort_p2_s1_mean_zig/s1_mean_zig_turb_param,label=r'$\omega^2$')
plt.plot(time,NS42_pv_p2_s1_mean_zig/s1_mean_zig_turb_param,label=r'$\Pi^2$')
plt.xlabel(r'$z_\mathrm{enc}/L_0$')
plt.ylabel(r'$(\langle b\rangle_\mathrm{T})_{z_{i,g}}/N^2(z_{i,g})_\mathrm{model}$')
plt.title(r'$Re_0=42$',fontsize=20)
plt.legend(loc='best',fontsize=20)
plt.tight_layout()
plt.savefig(opath+'s1_mean_zig_turb_param.pdf')
plt.show()

plt.figure(figsize=(5,5))
plt.grid(True)
plt.ylim(0,1)
plt.plot(time,np.sqrt(s1_var_zig)/(N**2*(delta_b+NS42_vort_turbareafrac_zig*(z_ig-z_if))),'C4',label=r'$Fr_0=0$')
#plt.plot(time,np.sqrt(s1_var_zig)/(N**2*(delta_b+NS42_pv_turbareafrac_zig*(z_ig-z_if))),label=r'$Fr_0=0$')
plt.plot(S20.z_enc/L0,np.sqrt(S20_s1_var_zig)/(N**2*(delta_b_S20+S20_vort_turbareafrac_zig*(S20.z_ig-S20.z_if))),'C5',label=r'$Fr_0=20$')
#plt.plot(S20.z_enc/L0,np.sqrt(S20_s1_var_zig)/(N**2*(delta_b_S20+S20_pv_turbareafrac_zig*(S20.z_ig-S20.z_if))),label=r'$Fr_0=20$')
plt.xlabel(r'$z_\mathrm{enc}/L_0$')
plt.ylabel(r'$(b_\mathrm{rms})_{z_{i,g}}/N^2[\delta_b + a_\mathrm{T}(z_{i,g}-z_{i,f})]$')
plt.legend(loc='best',fontsize=20)
plt.tight_layout()
plt.show()

# Presentations

plt.figure(figsize=(5,5))
plt.grid(True)
plt.ylim(1,1.4)
plt.xlim(-0.2,0)
plt.plot(np.mean(NS42_1.Rsv[-4:-1,:]/B0,axis=0),NS42_1.y/np.mean(NS42_1.z_enc[-4:-1]),label='$Fr_0=0$')
plt.plot(np.mean(S20.Rsv[-3:,:]/B0,axis=0),S20.y/np.mean(S20.z_enc[-3:]),label='$Fr_0=20$')
plt.xlabel(r'$\langle b^\prime w^\prime \rangle /B_0$')
plt.ylabel(r'$z/z_\mathrm{enc}$')
plt.legend(loc='best',fontsize=20)
plt.tight_layout()
plt.savefig(opath+'Presentations/s1_vflux_total_S20_S0.pdf')
plt.show()

plt.figure(figsize=(5,5))
plt.grid(True)
plt.ylim(1,1.4)
plt.xlim(0,1)
plt.plot(np.mean(NS42_1.Tke[-4:-1,:],axis=0)/(B0*np.mean(NS42_1.z_enc[-4:-1]))**(2./3.),NS42_1.y/np.mean(NS42_1.z_enc[-4:-1]),label='$Fr_0=0$')
plt.plot(np.mean(S20.Tke[-3:,:],axis=0)/(B0*np.mean(S20.z_enc[-3:]))**(2./3.),S20.y/np.mean(S20.z_enc[-3:]),label='$Fr_0=20$')
plt.xlabel(r'$(\langle u_i^{\prime 2}\rangle /2)/w_*^2$')
plt.ylabel(r'$z/z_\mathrm{enc}$')
plt.legend(loc='best',fontsize=20)
plt.tight_layout()
plt.savefig(opath+'Presentations/TKE_S20_S0.pdf')
plt.show()

plt.figure(figsize=(5,5))
plt.grid(True)
plt.xlim(-0.2,0.05)
plt.ylim(1,1.4)
plt.plot(np.mean(NS42_vort_int_1.int2[-4:-1,:]*(NS42_vort_int_1.P2v1[-4:-1,:]-NS42_vort_int_1.P2S1Mom1[-4:-1,:]*NS42_vort_int_1.P2VMom1[-4:-1,:])/B0,axis=0),NS42_1.y/np.mean(NS42_1.z_enc[-4:-1]),label=r'$Fr_0=0$')
plt.plot(np.mean(S20_vort_int.int2[-3:,:]*(S20_vort_int.P2v1[-3:,:]-S20_vort_int.P2S1Mom1[-3:,:]*S20_vort_int.P2VMom1[-3:,:])/B0,axis=0),S20.y/np.mean(S20.z_enc[-3:]),label=r'$Fr_0=20$')
plt.plot(np.mean(NS42_vort_int_1.int1[-4:-1,:]*(NS42_vort_int_1.P1v1[-4:-1,:]-NS42_vort_int_1.P1S1Mom1[-4:-1,:]*NS42_vort_int_1.P1VMom1[-4:-1,:])/B0,axis=0),NS42_1.y/np.mean(NS42_1.z_enc[-4:-1]),'C0',ls='--')
plt.plot(np.mean(S20_vort_int.int1[-3:,:]*(S20_vort_int.P1v1[-3:,:]-S20_vort_int.P1S1Mom1[-3:,:]*S20_vort_int.P1VMom1[-3:,:])/B0,axis=0),S20.y/np.mean(S20.z_enc[-3:]),'C1',ls='--')
plt.plot(np.mean(NS42_vort_int_1.int2[-4:-1,:]*NS42_vort_int_1.int1[-4:-1,:]*(NS42_vort_int_1.P2S1Mom1[-4:-1,:]-NS42_vort_int_1.P1S1Mom1[-4:-1,:])*(NS42_vort_int_1.P2VMom1[-4:-1,:]-NS42_vort_int_1.P1VMom1[-4:-1,:])/B0,axis=0),NS42_1.y/np.mean(NS42_1.z_enc[-4:-1]),'C0',ls='-.')
plt.plot(np.mean(S20_vort_int.int2[-3:,:]*S20_vort_int.int1[-3:,:]*(S20_vort_int.P2S1Mom1[-3:,:]-S20_vort_int.P1S1Mom1[-3:,:])*(S20_vort_int.P2VMom1[-3:,:]-S20_vort_int.P1VMom1[-3:,:])/B0,axis=0),S20.y/np.mean(S20.z_enc[-3:]),'C1',ls='-.')
# plt.plot(np.mean((NS42_vort_int_1.int2[-4:-1,:]*(NS42_vort_int_1.P2v1[-4:-1,:]-NS42_vort_int_1.P2S1Mom1[-4:-1,:]*NS42_vort_int_1.P2VMom1[-4:-1,:])+NS42_vort_int_1.int1[-4:-1,:]*(NS42_vort_int_1.P1v1[-4:-1,:]-NS42_vort_int_1.P1S1Mom1[-4:-1,:]*NS42_vort_int_1.P1VMom1[-4:-1,:])+NS42_vort_int_1.int2[-4:-1,:]*NS42_vort_int_1.int1[-4:-1,:]*(NS42_vort_int_1.P2S1Mom1[-4:-1,:]-NS42_vort_int_1.P1S1Mom1[-4:-1,:])*(NS42_vort_int_1.P2VMom1[-4:-1,:]-NS42_vort_int_1.P1VMom1[-4:-1,:]))/B0,axis=0),NS42_1.y/np.mean(NS42_1.z_enc[-4:-1]),'C1')
# plt.plot(np.mean(NS42_1.Rsv[-4:-1,:]/B0,axis=0),NS42_1.y/np.mean(NS42_1.z_enc[-4:-1]),'k--')
# plt.plot(np.mean((S20_vort_int.int2[-3:,:]*(S20_vort_int.P2v1[-3:,:]-S20_vort_int.P2S1Mom1[-3:,:]*S20_vort_int.P2VMom1[-3:,:])+S20_vort_int.int1[-3:,:]*(S20_vort_int.P1v1[-3:,:]-S20_vort_int.P1S1Mom1[-3:,:]*S20_vort_int.P1VMom1[-3:,:])+S20_vort_int.int2[-3:,:]*S20_vort_int.int1[-3:,:]*(S20_vort_int.P2S1Mom1[-3:,:]-S20_vort_int.P1S1Mom1[-3:,:])*(S20_vort_int.P2VMom1[-3:,:]-S20_vort_int.P1VMom1[-3:,:]))/B0,axis=0),S20.y/np.mean(S20.z_enc[-3:]),'C3')
# plt.plot(np.mean(S20.Rsv[-3:,:]/B0,axis=0),S20.y/np.mean(S20.z_enc[-3:]),'k--')
plt.xlabel(r'$f /B_0$')
plt.ylabel(r'$z/z_\mathrm{enc}$')
plt.legend(loc='best',fontsize=20)
plt.tight_layout()
plt.savefig(opath+'Presentations/s1_vflux_contributions_S20_S0.pdf')
plt.show()

f, (ax1,ax2) = plt.subplots(1,2,figsize=(10,5))
ax1.grid(True)
ax2.grid(True)
ax1.set_ylim(1,1.4)
ax1.set_xlim(-0.5,0.5)
ax2.set_ylim(0,0.5)
ax2.set_xlim(15,30)
ax1.plot(np.mean(NS42_vort_int_1.P2v1[-4:-1,:]-NS42_vort_int_1.P2S1Mom1[-4:-1,:]*NS42_vort_int_1.P2VMom1[-4:-1,:],axis=0)/B0,NS42_1.y/np.mean(NS42_1.z_enc[-4:-1]))
ax1.plot(np.mean(S20_vort_int.P2v1[-3:,:]-S20_vort_int.P2S1Mom1[-3:,:]*S20_vort_int.P2VMom1[-3:,:],axis=0)/B0,S20.y/np.mean(S20.z_enc[-3:]))
ax2.plot(time[1:-1],-runningmean(np.array(NS42_vort_p2_v1_zif)-np.array(NS42_vort_p2_s1_mean_zif)*np.array(NS42_vort_p2_w_mean_zif),1)/B0,label=r'$Fr_0=0$')
ax2.plot(S20.z_enc[1:-1]/L0,-runningmean(np.array(S20_vort_p2_v1_zif)-np.array(S20_vort_p2_s1_mean_zif)*np.array(S20_vort_p2_w_mean_zif),1)/B0,'C1',label=r'$Fr_0=20$')
ax1.axhline(np.mean(NS42_1.z_if[-4:-1]/NS42_1.z_enc[-4:-1]),0,0.05,c='C0')
ax1.axhline(np.mean(S20.z_if[-3:]/S20.z_enc[-3:]),0,0.05,c='C1')
ax1.set_xlabel(r'$\langle b^\prime w^\prime\rangle_\mathrm{T}/B_0$')
ax1.set_ylabel(r'$z/z_\mathrm{enc}$')
ax1.set_title(r'(a)',fontsize=20,loc='left')
ax2.set_xlabel(r'$z_\mathrm{enc}/L_0$')
ax2.set_ylabel(r'$-(\langle b^\prime w^\prime \rangle_\mathrm{T})_{z_{i,f}}/B_0$')
ax2.set_title(r'(b)',fontsize=20,loc='left')
ax2.legend(loc='best',fontsize=20)
plt.tight_layout()
plt.savefig(opath+'Presentations/s1_vflux_turb_height_time_S20_S0_vort.pdf')
plt.show()

f, (ax1,ax2) = plt.subplots(1,2,figsize=(10,5))
ax1.grid(True)
ax2.grid(True)
ax1.set_ylim(1,1.4)
ax2.set_xlim(15,30)
ax2.set_ylim(0,2.5)
ax1.plot(np.mean(np.sqrt(NS42_vort_int_1.P2S1Mom2[-4:-1,:]),axis=0)/(N**2*L0),NS42_1.y/np.mean(NS42_1.z_enc[-4:-1]),label=r'$Fr_0=0$')
ax1.plot(np.mean(np.sqrt(S20_vort_int.P2S1Mom2[-3:,:]),axis=0)/(N**2*L0),S20.y/np.mean(S20.z_enc[-3:],axis=0),label=r'$Fr_0=20$')
ax2.plot(time[1:-1],runningmean(np.sqrt(NS42_vort_p2_s1_var_zif),1)/(N**2*L0),'C0')
ax2.plot(S20.z_enc[1:-1]/L0,runningmean(np.sqrt(S20_vort_p2_s1_var_zif),1)/(N**2*L0),'C1')
ax1.axhline(np.mean(NS42_1.z_if[-4:-1]/NS42_1.z_enc[-4:-1]),0,0.05,c='C0')
ax1.axhline(np.mean(S20.z_if[-3:]/S20.z_enc[-3:]),0,0.05,c='C1')
ax1.set_xlabel(r'$(b_\mathrm{rms})_\mathrm{T} /(N^2L_0)$')
ax1.set_ylabel(r'$z/z_\mathrm{enc}$')
ax1.set_title('(a)',fontsize=20,loc='left')
ax2.set_xlabel(r'$z_\mathrm{enc}/L_0$')
ax2.set_ylabel(r'$((b_\mathrm{rms})_\mathrm{T})_{z_{i,f}}/(N^2L_0)$')
ax2.set_title('(b)',fontsize=20,loc='left')
ax1.legend(loc='best',fontsize=20)
plt.tight_layout(w_pad=2)
plt.savefig(opath+'Presentations/s1_rms_height_time_S20_S0_vort_zif.pdf')
plt.show()

f, (ax1,ax2) = plt.subplots(1,2,figsize=(10,5))
ax1.grid(True)
ax2.grid(True)
ax1.set_ylim(1,1.4)
ax1.set_xlim(-0.2,0)
ax2.set_ylim(0,0.2)
ax2.set_xlim(15,30)
ax1.plot(np.mean(rho_bw_turb,axis=0),NS42_1.y/np.mean(NS42_1.z_enc[-4:-1]))
ax1.plot(np.mean(rho_bw_S20_turb,axis=0),S20.y/np.mean(S20.z_enc[-3:]))
ax2.plot(time[1:-1],-runningmean(rho_bw_zif_turb,1),label=r'$Fr_0=0$')
ax2.plot(S20.z_enc[1:-1]/L0,-runningmean(rho_bw_zif_S20_turb,1),'C1',label=r'$Fr_0=20$')
ax1.axhline(np.mean(NS42_1.z_if[-4:-1]/NS42_1.z_enc[-4:-1]),0,0.05,c='C0')
ax1.axhline(np.mean(S20.z_if[-3:]/S20.z_enc[-3:]),0,0.05,c='C1')
ax1.set_xlabel(r'$(\rho_{bw})_\mathrm{T}$')
ax1.set_ylabel(r'$z/z_\mathrm{enc}$')
ax1.set_title('(a)',fontsize=20,loc='left')
ax2.set_xlabel(r'$z_\mathrm{enc}/L_0$')
ax2.set_ylabel(r'$-((\rho_{bw})_\mathrm{T})_{z_{i,f}}$')
ax2.set_title('(b)',fontsize=20,loc='left')
ax2.legend(loc='best',fontsize=20)
plt.tight_layout()
plt.savefig(opath+'Presentations/rho_bw_turb_height_time_S20_S0.pdf')
plt.show()

f, (ax1,ax2) = plt.subplots(1,2,figsize=(10,5))
ax1.grid(True)
ax2.grid(True)
ax1.set_ylim(1,1.4)
ax2.set_ylim(0,1)
ax2.set_xlim(15,30)
ax1.plot(np.mean(NS42_vort_int_1.int2[-4:-1,:],axis=0),NS42_1.y/np.mean(NS42_1.z_enc[-4:-1]))
ax1.plot(np.mean(S20_vort_int.int2[-3:,:],axis=0),S20.y/np.mean(S20.z_enc[-3:]))
ax2.plot(time[1:-1],runningmean(NS42_vort_turbareafrac_zif,1),'C0',label=r'$Fr_0=0$')
ax2.plot(S20.z_enc[1:-1]/L0,runningmean(S20_vort_turbareafrac_zif,1),'C1',label=r'$Fr_0=20$')
ax1.axhline(np.mean(NS42_1.z_if[-4:-1]/NS42_1.z_enc[-4:-1]),0,0.05,c='C0')
ax1.axhline(np.mean(S20.z_if[-3:]/S20.z_enc[-3:]),0,0.05,c='C1')
ax1.set_xlabel(r'$a_\mathrm{T}$')
ax1.set_ylabel(r'$z/z_\mathrm{enc}$')
ax1.set_title(r'(a)',fontsize=20,loc='left')
ax2.set_xlabel(r'$z_\mathrm{enc}/L_0$')
ax2.set_ylabel(r'$(a_\mathrm{T})_{z_{i,f}}$')
ax2.set_title(r'(b)',fontsize=20,loc='left')
ax2.legend(loc='lower left',fontsize=20)
plt.tight_layout(w_pad=2)
plt.savefig(opath+'Presentations/area_frac_height_time_S20_S0_zif_vort_EZ.pdf')
plt.show()


f, (ax1,ax2) = plt.subplots(1,2,figsize=(10,5))
ax1.grid(True)
ax2.grid(True)
ax1.set_ylim(0,1.4)
ax2.set_ylim(0,1)
ax1.plot(NS42_vort_int_1.int2[-2,:],NS42_1.y/NS42_1.z_enc[-2],label=r'$\omega^2$')
ax1.plot(NS42_pv_int_1.int2[-2,:],NS42_1.y/NS42_1.z_enc[-2],label=r'$\Pi^2$')
ax2.plot(time,NS42_vort_turbareafrac_zig,label=r'$\omega^2$')
ax2.plot(time,NS42_pv_turbareafrac_zig,label=r'$\Pi^2$')
ax2.plot(time,NS42_vort_turbareafrac_zif,'C0',ls='--')
ax2.plot(time,NS42_pv_turbareafrac_zif,'C1',ls='--')
ax1.set_xlabel(r'$a_\mathrm{T}$')
ax1.set_ylabel(r'$z/z_\mathrm{enc}$')
#ax1.set_title(r'(a) $Re_0 = 42$, $z_\mathrm{enc}/L_0 = 21$',fontsize=20,loc='left')
ax2.set_xlabel(r'$z_\mathrm{enc}/L_0$')
#ax2.set_ylabel(r'$(a_\mathrm{T})_{z_{i,g}}$')
#ax2.set_title(r'(b) $Re_0 = 42$',fontsize=20,loc='left')
ax1.legend(loc='lower left',fontsize=20)
#ax2.legend(loc='lower left',fontsize=20)
ax2.text(20,0.3,r'$(a_\mathrm{T})_{z_{i,g}}$')
ax2.text(20,0.9,r'$(a_\mathrm{T})_{z_{i,f}}$')
plt.tight_layout()
plt.savefig(opath+'Presentations/area_frac_height_time_zig_zif.pdf')
plt.show()

plt.figure(figsize=(5,5))
plt.grid(True)
plt.ylim(0,0.5)
plt.plot(time,NS42_vort_turbareafrac_zig,label=r'$\omega^2$')
plt.plot(time,NS42_pv_turbareafrac_zig,label=r'$\Pi^2$')
plt.xlabel(r'$z_\mathrm{enc}/L_0$')
plt.ylabel(r'$(a_\mathrm{T})_{z_{i,g}}$')
plt.legend(loc='lower left',fontsize=20)
plt.tight_layout()
plt.savefig(opath+'Presentations/area_frac_time.pdf')
plt.show()

plt.figure(figsize=(5,5))
plt.grid(True)
plt.ylim(1,1.4)
plt.xlim(0.8,1.4)
plt.plot(NS42_vort_int_1.P2S1Mom1[1,:]/(N**2*NS42_1.z_enc[1]),NS42_1.y/NS42_1.z_enc[1],label=r'$\omega^2$')
plt.plot(NS42_pv_int_1.P2S1Mom1[1,:]/(N**2*NS42_1.z_enc[1]),NS42_1.y/NS42_1.z_enc[1],label=r'$\Pi^2$')
plt.plot(NS42_vort_int_1.P1S1Mom1[1,:]/(N**2*NS42_1.z_enc[1]),NS42_1.y/NS42_1.z_enc[1],'C0',ls='--')
plt.plot(NS42_pv_int_1.P1S1Mom1[1,:]/(N**2*NS42_1.z_enc[1]),NS42_1.y/NS42_1.z_enc[1],'C1',ls='--')
plt.plot(N**2*NS42_1.y/(N**2*NS42_1.z_enc[1]),NS42_1.y/NS42_1.z_enc[1],'k--')
plt.text(1.15,1.1,r'$b_\mathrm{bg}=N^2z$',fontsize=20)
plt.xlabel(r'$\langle b \rangle / b_\mathrm{enc}$')
plt.ylabel(r'$z/z_\mathrm{enc}$')
plt.legend(loc='best',fontsize=20)
plt.tight_layout()
plt.savefig(opath+'Presentations/s1_mean_height.pdf')
plt.show()

plt.figure(figsize=(5,5))
plt.grid(True)
plt.ylim(1,1.4)
plt.plot(time,NS42_vort_p2_s1_mean_zig/(N**2*z_enc),label=r'$\omega^2$')
plt.plot(time,NS42_pv_p2_s1_mean_zig/(N**2*z_enc),label=r'$\Pi^2$')
plt.plot(time,NS42_vort_p1_s1_mean_zig/(N**2*z_enc),'C0',ls='--')
plt.plot(time,NS42_pv_p1_s1_mean_zig/(N**2*z_enc),'C1',ls='--')
plt.plot(time,z_ig/z_enc,'k--')
plt.text(20,1.24,r'$N^2z_{i,g}$',fontsize=20)
plt.xlabel(r'$z_\mathrm{enc}/L_0$')
plt.ylabel(r'$(\langle b \rangle)_{z_{i,g}}/b_\mathrm{enc}$')
plt.legend(loc='best',fontsize=20)
plt.tight_layout()
plt.savefig(opath+'Presentations/s1_mean_time.pdf')
plt.show()

plt.figure(figsize=(5,5))
plt.grid(True)
plt.ylim(0.9,1.05)
plt.plot(time,NS42_vort_p2_s1_mean_zig/(N**2*z_ig),label=r'$\omega^2$')
plt.plot(time,NS42_pv_p2_s1_mean_zig/(N**2*z_ig),label=r'$\Pi^2$')
plt.plot(time,NS42_vort_p1_s1_mean_zig/(N**2*z_ig),'C0',ls='--')
plt.plot(time,NS42_pv_p1_s1_mean_zig/(N**2*z_ig),'C1',ls='--')
plt.plot(time,z_ig/z_enc,'k--')
plt.text(20,1.24,r'$N^2z_{i,g}$',fontsize=20)
plt.xlabel(r'$z_\mathrm{enc}/L_0$')
plt.ylabel(r'$(\langle b \rangle)_{z_{i,g}}/N^2z_{i,g}$')
plt.legend(loc='best',fontsize=20)
plt.tight_layout()
plt.savefig(opath+'Presentations/s1_mean_time_zig.pdf')
plt.show()


