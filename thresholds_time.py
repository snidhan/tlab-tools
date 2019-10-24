#####################################################################
# Modules

from ReadStats import Statistics
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

colourmap_path = '/home/mpim/m300551/local/ScientificColourMaps5/'
opath = '/scratch/local1/m300551/ForKatherine/plots/3D/Re042/'
blues = matplotlib.cm.get_cmap('Blues')
######################################################################
# Grid 

path_NS42 = '/scratch/local1/m300551/ForKatherine/qCBL_3D/Re042/2560x576x2560/'
path_S20 = '/scratch/local1/m300551/ForKatherine/qCBL_3D/Re042/3072x960x4608-S20/'

#######################################################################
# Constants

nu_42 = 1./25000.
B0 = 0.005
N = np.sqrt(3)
L0 = (B0/N**3)**0.5
ceps = 0.1
cb=0.1

#######################################################################
# Stats

NS42 = Statistics(path_NS42+'stats/timmean/avg36500-47500.nc')
S20 = Statistics(path_S20+'stats/avg66000-84000.nc')

delta_b_NS42 = np.zeros(NS42.t_len)
b_delta_NS42 = np.zeros(NS42.t_len)
delta_b_S20 = np.zeros(S20.t_len)
b_delta_S20 = np.zeros(S20.t_len)
for t in range(0,S20.t_len):
    delta_b_NS42[t] = (N**2*NS42.z_ig[t] - NS42.rS[t,NS42.z_ig_arg[t]])/(NS42.rS_y[t,NS42.z_ig_arg[t]] - N**2)
    b_delta_NS42[t] = delta_b_NS42[t]*NS42.rS_y[t,NS42.z_ig_arg[t]]
    delta_b_S20[t] = (N**2*S20.z_ig[t] - S20.rS[t,S20.z_ig_arg[t]])/(S20.rS_y[t,S20.z_ig_arg[t]] - N**2)
    b_delta_S20[t] = delta_b_S20[t]*S20.rS_y[t,S20.z_ig_arg[t]]

########################################################################
# Conditional Stats

vort_thresholds = np.array([0.16,0.32,0.48,0.64,0.8,0.96,1.12,1.28,1.44,1.6,3.2,8,16])
filenames = ['0-16','0-32','0-48','0-64','0-8','0-96','1-12','1-28','1-44','1-6','3-2','8','16']
vort_turb_area_fracs_NS42 = np.zeros((NS42.t_len,NS42.y_len,len(vort_thresholds)))
vort_turb_area_fracs_S20 = np.zeros((S20.t_len,S20.y_len,len(vort_thresholds)))
vort_turb_s1_mean_NS42 = np.zeros((NS42.t_len,NS42.y_len,len(vort_thresholds)))
vort_nonturb_s1_mean_NS42 = np.zeros((NS42.t_len,NS42.y_len,len(vort_thresholds)))
vort_turb_s1_mean_S20 = np.zeros((S20.t_len,S20.y_len,len(vort_thresholds)))
vort_nonturb_s1_mean_S20 = np.zeros((S20.t_len,S20.y_len,len(vort_thresholds)))
vort_turb_s1_var_NS42 = np.zeros((NS42.t_len,NS42.y_len,len(vort_thresholds)))
vort_nonturb_s1_var_NS42 = np.zeros((NS42.t_len,NS42.y_len,len(vort_thresholds)))
vort_turb_s1_var_S20 = np.zeros((S20.t_len,S20.y_len,len(vort_thresholds)))
vort_nonturb_s1_var_S20 = np.zeros((S20.t_len,S20.y_len,len(vort_thresholds)))
vort_turb_w_mean_NS42 = np.zeros((NS42.t_len,NS42.y_len,len(vort_thresholds)))
vort_nonturb_w_mean_NS42 = np.zeros((NS42.t_len,NS42.y_len,len(vort_thresholds)))
vort_turb_w_mean_S20 = np.zeros((S20.t_len,S20.y_len,len(vort_thresholds)))
vort_nonturb_w_mean_S20 = np.zeros((S20.t_len,S20.y_len,len(vort_thresholds)))
vort_turb_w_var_NS42 = np.zeros((NS42.t_len,NS42.y_len,len(vort_thresholds)))
vort_nonturb_w_var_NS42 = np.zeros((NS42.t_len,NS42.y_len,len(vort_thresholds)))
vort_turb_w_var_S20 = np.zeros((S20.t_len,S20.y_len,len(vort_thresholds)))
vort_nonturb_w_var_S20 = np.zeros((S20.t_len,S20.y_len,len(vort_thresholds)))
vort_turb_s1_flux_NS42 = np.zeros((NS42.t_len,NS42.y_len,len(vort_thresholds)))
vort_nonturb_s1_flux_NS42 = np.zeros((NS42.t_len,NS42.y_len,len(vort_thresholds)))
vort_turb_s1_flux_S20 = np.zeros((S20.t_len,S20.y_len,len(vort_thresholds)))
vort_nonturb_s1_flux_S20 = np.zeros((S20.t_len,S20.y_len,len(vort_thresholds)))
vort_turb_s1_flux_zif_NS42 = np.zeros((NS42.t_len,len(vort_thresholds)))
vort_turb_s1_flux_zif_S20 = np.zeros((S20.t_len,len(vort_thresholds)))
vort_turb_s1_mean_zif_NS42 = np.zeros((NS42.t_len,len(vort_thresholds)))
vort_turb_s1_mean_zif_S20 = np.zeros((S20.t_len,len(vort_thresholds)))
vort_turb_w_mean_zif_NS42 = np.zeros((NS42.t_len,len(vort_thresholds)))
vort_turb_w_mean_zif_S20 = np.zeros((S20.t_len,len(vort_thresholds)))
for t in range(0,NS42.t_len):
    for i in range(0,len(vort_thresholds)):
        vort_turb_file = path_NS42+'stats/gate-vorticity/gate-'+filenames[i]+'/Partition2/cavg36500-47500.nc'
        vort_turb_data = nc.Dataset(vort_turb_file,'r')
        vort_turb_s1mean = vort_turb_data.variables['Scalar1Mom1'][t,:]
        vort_turb_s1_mean_zif_NS42[t,i] = vort_turb_s1mean[NS42.z_if_arg[t]]
        vort_turb_wmean = vort_turb_data.variables['VMom1'][t,:]
        vort_turb_w_mean_zif_NS42[t,i] = vort_turb_wmean[NS42.z_if_arg[t]]
        vort_turb_flux_file = path_NS42+'stats/gate-vorticity/gate-'+filenames[i]+'/Partition2/avgMom36500-47500.nc'
        vort_turb_flux_data = nc.Dataset(vort_turb_flux_file,'r')
        vort_turb_s1flux = vort_turb_flux_data.variables['v1Mom1'][t,:]
        vort_turb_s1_flux_zif_NS42[t,i] = vort_turb_s1flux[NS42.z_if_arg[t]]
vort_turb_s1_flux_zif_NS42_timmean = np.mean(vort_turb_s1_flux_zif_NS42-(vort_turb_s1_mean_zif_NS42*vort_turb_w_mean_zif_NS42),axis=0)
for t in range(0,S20.t_len):
    for i in range(0,len(vort_thresholds)):
        vort_turb_file = path_S20+'stats/gate-vorticity/gate-'+filenames[i]+'/Partition2/cavg66000-84000.nc'
        vort_turb_data = nc.Dataset(vort_turb_file,'r')
        vort_turb_s1mean = vort_turb_data.variables['Scalar1Mom1'][t,:]
        vort_turb_s1_mean_zif_S20[t,i] = vort_turb_s1mean[S20.z_if_arg[t]]
        vort_turb_wmean = vort_turb_data.variables['VMom1'][t,:]
        vort_turb_w_mean_zif_S20[t,i] = vort_turb_wmean[S20.z_if_arg[t]]
        vort_turb_flux_file = path_S20+'stats/gate-vorticity/gate-'+filenames[i]+'/Partition2/avgMom66000-84000.nc'
        vort_turb_flux_data = nc.Dataset(vort_turb_flux_file,'r')
        vort_turb_s1flux = vort_turb_flux_data.variables['v1Mom1'][t,:]
        vort_turb_s1_flux_zif_S20[t,i] = vort_turb_s1flux[S20.z_if_arg[t]]
vort_turb_s1_flux_zif_S20_timmean = np.mean(vort_turb_s1_flux_zif_S20-(vort_turb_s1_mean_zif_S20*vort_turb_w_mean_zif_S20),axis=0)
for i in range(0,len(vort_thresholds)):
    vort_gate_file = path_NS42+'stats/gate-vorticity/gate-'+filenames[i]+'/int36500-47500.nc'
    vort_gate_data = nc.Dataset(vort_gate_file,'r')
    vort_turb_area = vort_gate_data.variables['Partition2'][0,:]
    vort_turb_area_fracs_NS42[:,:,i] = vort_turb_area
    vort_turb_area_fracs_NS42_timmean = np.mean(vort_turb_area_fracs_NS42,axis=0)
    vort_turb_file = path_NS42+'stats/gate-vorticity/gate-'+filenames[i]+'/Partition2/cavg36500-47500.nc'
    vort_turb_data = nc.Dataset(vort_turb_file,'r')
    vort_turb_s1mean = vort_turb_data.variables['Scalar1Mom1'][:,:]
    vort_turb_s1_mean_NS42[:,:,i] = vort_turb_s1mean
    vort_turb_s1_mean_NS42_timmean = np.mean(vort_turb_s1_mean_NS42,axis=0)
    vort_turb_s1var = vort_turb_data.variables['Scalar1Mom2'][:,:]
    vort_turb_s1_var_NS42[:,:,i] = vort_turb_s1var
    vort_turb_s1_var_NS42_timmean = np.mean(vort_turb_s1_var_NS42,axis=0)
    vort_turb_wmean = vort_turb_data.variables['VMom1'][:,:]
    vort_turb_w_mean_NS42[:,:,i] = vort_turb_wmean
    vort_turb_w_mean_NS42_timmean = np.mean(vort_turb_w_mean_NS42,axis=0)
    vort_turb_wvar = vort_turb_data.variables['VMom2'][:,:]
    vort_turb_w_var_NS42[:,:,i] = vort_turb_wvar
    vort_turb_w_var_NS42_timmean = np.mean(vort_turb_w_var_NS42,axis=0)
    vort_turb_flux_file = path_NS42+'stats/gate-vorticity/gate-'+filenames[i]+'/Partition2/avgMom36500-47500.nc'
    vort_turb_flux_data = nc.Dataset(vort_turb_flux_file,'r')
    vort_turb_s1flux = vort_turb_flux_data.variables['v1Mom1'][:,:]
    vort_turb_s1_flux_NS42[:,:,i] = vort_turb_s1flux
    vort_turb_s1_flux_NS42_timmean = np.mean(vort_turb_s1_flux_NS42,axis=0)
    vort_nonturb_file = path_NS42+'stats/gate-vorticity/gate-'+filenames[i]+'/Partition1/cavg36500-47500.nc'
    vort_nonturb_data = nc.Dataset(vort_nonturb_file,'r')
    vort_nonturb_s1mean = vort_nonturb_data.variables['Scalar1Mom1'][:,:]
    vort_nonturb_s1_mean_NS42[:,:,i] = vort_nonturb_s1mean
    vort_nonturb_s1_mean_NS42_timmean = np.mean(vort_nonturb_s1_mean_NS42,axis=0)  
    vort_nonturb_s1var = vort_nonturb_data.variables['Scalar1Mom2'][:,:]
    vort_nonturb_s1_var_NS42[:,:,i] = vort_nonturb_s1var
    vort_nonturb_s1_var_NS42_timmean = np.mean(vort_nonturb_s1_var_NS42,axis=0)
    vort_nonturb_wmean = vort_nonturb_data.variables['VMom1'][:,:]
    vort_nonturb_w_mean_NS42[:,:,i] = vort_nonturb_wmean
    vort_nonturb_w_mean_NS42_timmean = np.mean(vort_nonturb_w_mean_NS42,axis=0)
    vort_nonturb_wvar = vort_nonturb_data.variables['VMom2'][:,:]
    vort_nonturb_w_var_NS42[:,:,i] = vort_nonturb_wvar
    vort_nonturb_w_var_NS42_timmean = np.mean(vort_nonturb_w_var_NS42,axis=0)
    vort_nonturb_flux_file = path_NS42+'stats/gate-vorticity/gate-'+filenames[i]+'/Partition1/avgMom36500-47500.nc'
    vort_nonturb_flux_data = nc.Dataset(vort_nonturb_flux_file,'r')
    vort_nonturb_s1flux = vort_nonturb_flux_data.variables['v1Mom1'][:,:]
    vort_nonturb_s1_flux_NS42[:,:,i] = vort_nonturb_s1flux
    vort_nonturb_s1_flux_NS42_timmean = np.mean(vort_nonturb_s1_flux_NS42,axis=0)
    vort_gate_file = path_S20+'stats/gate-vorticity/gate-'+filenames[i]+'/int66000-84000.nc'
    vort_gate_data = nc.Dataset(vort_gate_file,'r')
    vort_turb_area = vort_gate_data.variables['Partition2'][0,:]
    vort_turb_area_fracs_S20[:,:,i] = vort_turb_area
    vort_turb_area_fracs_S20_timmean = np.mean(vort_turb_area_fracs_S20,axis=0)
    vort_turb_file = path_S20+'stats/gate-vorticity/gate-'+filenames[i]+'/Partition2/cavg66000-84000.nc'
    vort_turb_data = nc.Dataset(vort_turb_file,'r')
    vort_turb_s1mean = vort_turb_data.variables['Scalar1Mom1'][:,:]
    vort_turb_s1_mean_S20[:,:,i] = vort_turb_s1mean
    vort_turb_s1_mean_S20_timmean = np.mean(vort_turb_s1_mean_S20,axis=0)
    vort_turb_s1var = vort_turb_data.variables['Scalar1Mom2'][:,:]
    vort_turb_s1_var_S20[:,:,i] = vort_turb_s1var
    vort_turb_s1_var_S20_timmean = np.mean(vort_turb_s1_var_S20,axis=0)
    vort_turb_wmean = vort_turb_data.variables['VMom1'][:,:]
    vort_turb_w_mean_S20[:,:,i] = vort_turb_wmean
    vort_turb_w_mean_S20_timmean = np.mean(vort_turb_w_mean_S20,axis=0)
    vort_turb_wvar = vort_turb_data.variables['VMom2'][:,:]
    vort_turb_w_var_S20[:,:,i] = vort_turb_wvar
    vort_turb_w_var_S20_timmean = np.mean(vort_turb_w_var_S20,axis=0)
    vort_turb_flux_file = path_S20+'stats/gate-vorticity/gate-'+filenames[i]+'/Partition2/avgMom66000-84000.nc'
    vort_turb_flux_data = nc.Dataset(vort_turb_flux_file,'r')
    vort_turb_s1flux = vort_turb_flux_data.variables['v1Mom1'][:,:]
    vort_turb_s1_flux_S20[:,:,i] = vort_turb_s1flux
    vort_turb_s1_flux_S20_timmean = np.mean(vort_turb_s1_flux_S20,axis=0)
    vort_nonturb_file = path_S20+'stats/gate-vorticity/gate-'+filenames[i]+'/Partition1/cavg66000-84000.nc'
    vort_nonturb_data = nc.Dataset(vort_nonturb_file,'r')
    vort_nonturb_s1mean = vort_nonturb_data.variables['Scalar1Mom1'][:,:]
    vort_nonturb_s1_mean_S20[:,:,i] = vort_nonturb_s1mean
    vort_nonturb_s1_mean_S20_timmean = np.mean(vort_nonturb_s1_mean_S20,axis=0)
    vort_nonturb_s1var = vort_nonturb_data.variables['Scalar1Mom2'][:,:]
    vort_nonturb_s1_var_S20[:,:,i] = vort_nonturb_s1var
    vort_nonturb_s1_var_S20_timmean = np.mean(vort_nonturb_s1_var_S20,axis=0)
    vort_nonturb_wmean = vort_nonturb_data.variables['VMom1'][:,:]
    vort_nonturb_w_mean_S20[:,:,i] = vort_nonturb_wmean
    vort_nonturb_w_mean_S20_timmean = np.mean(vort_nonturb_w_mean_S20,axis=0)
    vort_nonturb_wvar = vort_nonturb_data.variables['VMom2'][:,:]
    vort_nonturb_w_var_S20[:,:,i] = vort_nonturb_wvar
    vort_nonturb_w_var_S20_timmean = np.mean(vort_nonturb_w_var_S20,axis=0)
    vort_nonturb_flux_file = path_S20+'stats/gate-vorticity/gate-'+filenames[i]+'/Partition1/avgMom66000-84000.nc'
    vort_nonturb_flux_data = nc.Dataset(vort_nonturb_flux_file,'r')
    vort_nonturb_s1flux = vort_nonturb_flux_data.variables['v1Mom1'][:,:]
    vort_nonturb_s1_flux_S20[:,:,i] = vort_nonturb_s1flux
    vort_nonturb_s1_flux_S20_timmean = np.mean(vort_nonturb_s1_flux_S20,axis=0)


pvthresholds = range(-9,5)
pv_turb_area_fracs_NS42 = np.zeros((NS42.t_len,NS42.y_len,len(pvthresholds)))
pv_turb_area_fracs_S20 = np.zeros((S20.t_len,S20.y_len,len(pvthresholds)))
pv_turb_s1_mean_NS42 = np.zeros((NS42.t_len,NS42.y_len,len(pvthresholds)))
pv_nonturb_s1_mean_NS42 = np.zeros((NS42.t_len,NS42.y_len,len(pvthresholds)))
pv_turb_s1_mean_S20 = np.zeros((S20.t_len,S20.y_len,len(pvthresholds)))
pv_nonturb_s1_mean_S20 = np.zeros((S20.t_len,S20.y_len,len(pvthresholds)))
pv_turb_s1_var_NS42 = np.zeros((NS42.t_len,NS42.y_len,len(pvthresholds)))
pv_nonturb_s1_var_NS42 = np.zeros((NS42.t_len,NS42.y_len,len(pvthresholds)))
pv_turb_s1_var_S20 = np.zeros((S20.t_len,S20.y_len,len(pvthresholds)))
pv_nonturb_s1_var_S20 = np.zeros((S20.t_len,S20.y_len,len(pvthresholds)))
pv_turb_s1_flux_NS42 =  np.zeros((NS42.t_len,NS42.y_len,len(pvthresholds)))
pv_nonturb_s1_flux_NS42 =  np.zeros((NS42.t_len,NS42.y_len,len(pvthresholds)))
pv_turb_s1_flux_S20 =  np.zeros((S20.t_len,S20.y_len,len(pvthresholds)))
pv_nonturb_s1_flux_S20 =  np.zeros((S20.t_len,S20.y_len,len(pvthresholds)))
pv_turb_w_mean_NS42 = np.zeros((NS42.t_len,NS42.y_len,len(pvthresholds)))
pv_nonturb_w_mean_NS42 = np.zeros((NS42.t_len,NS42.y_len,len(pvthresholds)))
pv_turb_w_mean_S20 = np.zeros((S20.t_len,S20.y_len,len(pvthresholds)))
pv_nonturb_w_mean_S20 = np.zeros((S20.t_len,S20.y_len,len(pvthresholds)))
for i in range(0,len(pvthresholds)):
    pv_gate_file = path_NS42+'stats/gate-pv/gate'+str(pvthresholds[i])+'/int36500-47500.nc'
    pv_gate_data = nc.Dataset(pv_gate_file,'r')
    pv_turb_area = pv_gate_data.variables['Partition2'][0,:]
    pv_turb_area_fracs_NS42[:,:,i] = pv_turb_area
    pv_turb_area_fracs_NS42_timmean = np.mean(pv_turb_area_fracs_NS42,axis=0)
    pv_turb_file = path_NS42+'stats/gate-pv/gate'+str(pvthresholds[i])+'/Partition2/cavg36500-47500.nc'
    pv_turb_data = nc.Dataset(pv_turb_file,'r')
    pv_turb_s1mean = pv_turb_data.variables['Scalar1Mom1'][:,:]
    pv_turb_s1_mean_NS42[:,:,i] = pv_turb_s1mean 
    pv_turb_s1_mean_NS42_timmean = np.mean(pv_turb_s1_mean_NS42,axis=0)
    pv_turb_s1var = pv_turb_data.variables['Scalar1Mom2'][:,:]
    pv_turb_s1_var_NS42[:,:,i] = pv_turb_s1var
    pv_turb_s1_var_NS42_timmean = np.mean(pv_turb_s1_var_NS42,axis=0)
    pv_turb_wmean =  pv_turb_data.variables['VMom1'][:,:]
    pv_turb_w_mean_NS42[:,:,i] = pv_turb_wmean
    pv_turb_w_mean_NS42_timmean = np.mean(pv_turb_w_mean_NS42,axis=0)
    pv_turb_flux_file = path_NS42+'stats/gate-pv/gate'+str(pvthresholds[i])+'/Partition2/avgMom36500-47500.nc'
    pv_turb_flux_data = nc.Dataset(pv_turb_flux_file,'r')
    pv_turb_s1flux = pv_turb_flux_data.variables['v1Mom1'][:,:]
    pv_turb_s1_flux_NS42[:,:,i] = pv_turb_s1flux
    pv_turb_s1_flux_NS42_timmean = np.mean(pv_turb_s1_flux_NS42,axis=0)
    pv_nonturb_file = path_NS42+'stats/gate-pv/gate'+str(pvthresholds[i])+'/Partition1/cavg36500-47500.nc'
    pv_nonturb_data = nc.Dataset(pv_nonturb_file,'r')
    pv_nonturb_s1mean = pv_nonturb_data.variables['Scalar1Mom1'][:,:]
    pv_nonturb_s1_mean_NS42[:,:,i] = pv_nonturb_s1mean 
    pv_nonturb_s1_mean_NS42_timmean = np.mean(pv_nonturb_s1_mean_NS42,axis=0)  
    pv_nonturb_s1var = pv_nonturb_data.variables['Scalar1Mom2'][:,:]
    pv_nonturb_s1_var_NS42[:,:,i] = pv_nonturb_s1var
    pv_nonturb_s1_var_NS42_timmean = np.mean(pv_nonturb_s1_var_NS42,axis=0)
    pv_nonturb_wmean =  pv_nonturb_data.variables['VMom1'][:,:]
    pv_nonturb_w_mean_NS42[:,:,i] = pv_nonturb_wmean
    pv_nonturb_w_mean_NS42_timmean = np.mean(pv_nonturb_w_mean_NS42,axis=0)
    pv_nonturb_flux_file = path_NS42+'stats/gate-pv/gate'+str(pvthresholds[i])+'/Partition1/avgMom36500-47500.nc'
    pv_nonturb_flux_data = nc.Dataset(pv_nonturb_flux_file,'r')
    pv_nonturb_s1flux = pv_nonturb_flux_data.variables['v1Mom1'][:,:]
    pv_nonturb_s1_flux_NS42[:,:,i] = pv_nonturb_s1flux
    pv_nonturb_s1_flux_NS42_timmean = np.mean(pv_nonturb_s1_flux_NS42,axis=0)
    pv_gate_file = path_S20+'stats/gate-pv/gate'+str(pvthresholds[i])+'/int66000-84000.nc'
    pv_gate_data = nc.Dataset(pv_gate_file,'r')
    pv_turb_area = pv_gate_data.variables['Partition2'][0,:]
    pv_turb_area_fracs_S20[:,:,i] = pv_turb_area
    pv_turb_area_fracs_S20_timmean = np.mean(pv_turb_area_fracs_S20,axis=0)
    pv_turb_file = path_S20+'stats/gate-pv/gate'+str(pvthresholds[i])+'/Partition2/cavg66000-84000.nc'
    pv_turb_data = nc.Dataset(pv_turb_file,'r')
    pv_turb_s1mean = pv_turb_data.variables['Scalar1Mom1'][:,:]
    pv_turb_s1_mean_S20[:,:,i] = pv_turb_s1mean 
    pv_turb_s1_mean_S20_timmean = np.mean(pv_turb_s1_mean_S20,axis=0)
    pv_turb_s1var = pv_turb_data.variables['Scalar1Mom2'][:,:]
    pv_turb_s1_var_S20[:,:,i] = pv_turb_s1var
    pv_turb_s1_var_S20_timmean = np.mean(pv_turb_s1_var_S20,axis=0)
    pv_turb_wmean =  pv_turb_data.variables['VMom1'][:,:]
    pv_turb_w_mean_S20[:,:,i] = pv_turb_wmean
    pv_turb_w_mean_S20_timmean = np.mean(pv_turb_w_mean_S20,axis=0)
    pv_turb_flux_file = path_S20+'stats/gate-pv/gate'+str(pvthresholds[i])+'/Partition2/avgMom66000-84000.nc'
    pv_turb_flux_data = nc.Dataset(pv_turb_flux_file,'r')
    pv_turb_s1flux = pv_turb_flux_data.variables['v1Mom1'][:,:]
    pv_turb_s1_flux_S20[:,:,i] = pv_turb_s1flux
    pv_turb_s1_flux_S20_timmean = np.mean(pv_turb_s1_flux_S20,axis=0)
    pv_nonturb_file = path_S20+'stats/gate-pv/gate'+str(pvthresholds[i])+'/Partition1/cavg66000-84000.nc'
    pv_nonturb_data = nc.Dataset(pv_nonturb_file,'r')
    pv_nonturb_s1mean = pv_nonturb_data.variables['Scalar1Mom1'][:,:]
    pv_nonturb_s1_mean_S20[:,:,i] = pv_nonturb_s1mean 
    pv_nonturb_s1_mean_S20_timmean = np.mean(pv_nonturb_s1_mean_S20,axis=0)
    pv_nonturb_s1var = pv_nonturb_data.variables['Scalar1Mom2'][:,:]
    pv_nonturb_s1_var_S20[:,:,i] = pv_nonturb_s1var
    pv_nonturb_s1_var_S20_timmean = np.mean(pv_nonturb_s1_var_S20,axis=0)
    pv_nonturb_wmean =  pv_nonturb_data.variables['VMom1'][:,:]
    pv_nonturb_w_mean_S20[:,:,i] = pv_nonturb_wmean
    pv_nonturb_w_mean_S20_timmean = np.mean(pv_nonturb_w_mean_S20,axis=0)
    pv_nonturb_flux_file = path_S20+'stats/gate-pv/gate'+str(pvthresholds[i])+'/Partition1/avgMom66000-84000.nc'
    pv_nonturb_flux_data = nc.Dataset(pv_nonturb_flux_file,'r')
    pv_nonturb_s1flux = pv_nonturb_flux_data.variables['v1Mom1'][:,:]
    pv_nonturb_s1_flux_S20[:,:,i] = pv_nonturb_s1flux
    pv_nonturb_s1_flux_S20_timmean = np.mean(pv_nonturb_s1_flux_S20,axis=0)

pvthresholds = np.exp(pvthresholds)
########################################################################
# Colourmaps

imola_data = np.loadtxt(colourmap_path+'imola/imola.txt')
imola_map = LinearSegmentedColormap.from_list('imola',imola_data)
vik_data = np.loadtxt(colourmap_path+'vik/vik.txt')
vik_map = LinearSegmentedColormap.from_list('vik',vik_data)
########################################################################
# Plot

f, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,sharex='col',sharey='row',figsize=(10,10))
ax1.tick_params(bottom=True,top=True,left=True,right=True)
ax2.tick_params(bottom=True,top=True,left=True,right=True)
ax3.tick_params(bottom=True,top=True,left=True,right=True)
ax4.tick_params(bottom=True,top=True,left=True,right=True)
ax1.set_ylim(1,1.5)
ax3.set_ylim(1,1.5)
ax2.set_xlim(np.log10(pvthresholds[0]/(cb*ceps*N**6*42**2*(np.mean(NS42.z_enc)/L0)**(-4./3.))),np.log10(pvthresholds[-1]/(cb*ceps*N**6*42**2*(np.mean(S20.z_enc)/L0)**(-4./3.))))
ax4.set_xlim(np.log10(pvthresholds[0]/(cb*ceps*N**6*42**2*(np.mean(NS42.z_enc)/L0)**(-4./3.))),np.log10(pvthresholds[-1]/(cb*ceps*N**6*42**2*(np.mean(S20.z_enc)/L0)**(-4./3.))))
ax4.set_xticks([-4,-3,-2,-1,0])
cs1 = ax1.contourf(np.log10(vort_thresholds**2/(ceps*B0/nu_42)),NS42.y/np.mean(NS42.z_enc),vort_turb_area_fracs_NS42_timmean,cmap=imola_map,levels=np.arange(0,1.1,0.1))
cs2 = ax2.contourf(np.log10(pvthresholds/(cb*ceps*N**6*42**2*(np.mean(NS42.z_enc)/L0)**(-4./3.))),NS42.y/np.mean(NS42.z_enc),pv_turb_area_fracs_NS42_timmean,cmap=imola_map,levels=np.arange(0,1.1,0.1))
cs3 = ax3.contourf(np.log10(vort_thresholds**2/(ceps*B0/nu_42)),S20.y/np.mean(S20.z_enc),vort_turb_area_fracs_S20_timmean,cmap=imola_map,levels=np.arange(0,1.1,0.1))
cs4 = ax4.contourf(np.log10(pvthresholds/(cb*ceps*N**6*42**2*(np.mean(S20.z_enc)/L0)**(-4./3.))),S20.y/np.mean(S20.z_enc),pv_turb_area_fracs_S20_timmean,cmap=imola_map,levels=np.arange(0,1.1,0.1))
cbar_ax = f.add_axes([0.15,0.08,0.8,0.02])
cbar = f.colorbar(cs1,cax=cbar_ax,orientation='horizontal')
ax1.axhline(np.mean(NS42.z_ig)/np.mean(NS42.z_enc),0,0.05,color='k',linewidth=2)
ax1.axhline(np.mean(NS42.z_if)/np.mean(NS42.z_enc),0,0.05,color='k',linewidth=2)
ax2.axhline(np.mean(NS42.z_ig)/np.mean(NS42.z_enc),0,0.05,color='k',linewidth=2)
ax2.axhline(np.mean(NS42.z_if)/np.mean(NS42.z_enc),0,0.05,color='k',linewidth=2)
ax3.axhline(np.mean(S20.z_ig)/np.mean(S20.z_enc),0,0.05,color='k',linewidth=2)
ax3.axhline(np.mean(S20.z_if)/np.mean(S20.z_enc),0,0.05,color='k',linewidth=2)
ax4.axhline(np.mean(S20.z_ig)/np.mean(S20.z_enc),0,0.05,color='k',linewidth=2)
ax4.axhline(np.mean(S20.z_if)/np.mean(S20.z_enc),0,0.05,color='k',linewidth=2)
cbar.ax.set_xlabel(r'$a_\mathrm{T}$')
ax1.set_ylabel(r'$z/z_\mathrm{enc}$')
ax3.set_ylabel(r'$z/z_\mathrm{enc}$')
ax3.set_xlabel(r'$\log_{10}(\omega_\mathrm{th}^2/\omega_0^2)$')
ax4.set_xlabel(r'$\log_{10}(\Pi^2_\mathrm{th}/\Pi_0^2)$')
ax1.set_title('(a) $Fr_0=0$',loc='left',fontsize=20)
ax2.set_title('(b) $Fr_0=0$',loc='left',fontsize=20)
ax3.set_title('(c) $Fr_0=20$',loc='left',fontsize=20)
ax4.set_title('(d) $Fr_0=20$',loc='left',fontsize=20)
ax1.axvline(-0.46,0,NS42.y[-1]/np.mean(NS42.z_enc),color='k',linestyle='--')
ax2.axvline(-1.99,0,NS42.y[-1]/np.mean(NS42.z_enc),color='k',linestyle='--')
ax3.axvline(-0.91,0,S20.y[-1]/np.mean(S20.z_enc),color='k',linestyle='--')
ax4.axvline(-2.65,0,S20.y[-1]/np.mean(S20.z_enc),color='k',linestyle='--')
plt.tight_layout(rect=[0,0.1,1,1])
plt.savefig(opath+'area_frac_threshold_S20_S0_timmean.pdf')
plt.show()


f, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,sharex='col',sharey='row',figsize=(10,10))
ax1.tick_params(bottom=True,top=True,left=True,right=True)
ax2.tick_params(bottom=True,top=True,left=True,right=True)
ax3.tick_params(bottom=True,top=True,left=True,right=True)
ax4.tick_params(bottom=True,top=True,left=True,right=True)
ax1.set_ylim(1,1.5)
ax3.set_ylim(1,1.5)
cs1 = ax1.contourf(np.log10(vort_thresholds**2/(ceps*B0/nu_42)),NS42.y/np.mean(NS42.z_enc),(vort_turb_s1_mean_NS42_timmean-vort_nonturb_s1_mean_NS42_timmean)/np.mean(NS42.b_enc),cmap=vik_map,levels=np.arange(-0.5,0.55,0.05),extend='min')
cs2 = ax2.contourf(np.log10(vort_thresholds**2/(ceps*B0/nu_42)),NS42.y/np.mean(NS42.z_enc),(vort_turb_w_mean_NS42_timmean-vort_nonturb_w_mean_NS42_timmean)/(np.mean(NS42.z_enc)*B0)**(1./3.),cmap=vik_map,levels=np.arange(-0.5,0.55,0.05),extend='min')
cs3 = ax3.contourf(np.log10(vort_thresholds**2/(ceps*B0/nu_42)),S20.y/np.mean(S20.z_enc),(vort_turb_s1_mean_S20_timmean-vort_nonturb_s1_mean_S20_timmean)/np.mean(S20.b_enc),cmap=vik_map,levels=np.arange(-0.5,0.55,0.05),extend='min')
cs4 = ax4.contourf(np.log10(vort_thresholds**2/(ceps*B0/nu_42)),S20.y/np.mean(S20.z_enc),(vort_turb_w_mean_S20_timmean-vort_nonturb_w_mean_S20_timmean)/(np.mean(S20.z_enc)*B0)**(1./3.),cmap=vik_map,levels=np.arange(-0.5,0.55,0.05),extend='min')
cbar_ax = f.add_axes([0.15,0.05,0.8,0.02])
cbar = f.colorbar(cs1,cax=cbar_ax,orientation='horizontal')
ax1.set_ylabel(r'$z/z_\mathrm{enc}$')
ax3.set_ylabel(r'$z/z_\mathrm{enc}$')
ax3.set_xlabel(r'$\log_{10}(\omega_\mathrm{th}^2/\omega_0^2)$')
ax4.set_xlabel(r'$\log_{10}(\omega_\mathrm{th}^2/\omega_0^2)$')
ax1.set_title(r'(a) $Fr_0=0$, $(\langle b \rangle_\mathrm{T}-\langle b \rangle_\mathrm{NT})/b_\mathrm{enc}$',loc='left',fontsize=20)
ax2.set_title(r'(b) $Fr_0=0$, $(\langle w \rangle_\mathrm{T}-\langle w \rangle_\mathrm{NT})/w_*$',loc='left',fontsize=20)
ax3.set_title(r'(c) $Fr_0=20$, $(\langle b \rangle_\mathrm{T}-\langle b \rangle_\mathrm{NT})/b_\mathrm{enc}$',loc='left',fontsize=20)
ax4.set_title(r'(d) $Fr_0=20$, $(\langle w \rangle_\mathrm{T}-\langle w \rangle_\mathrm{NT})/w_*$',loc='left',fontsize=20)
ax1.axhline(np.mean(NS42.z_ig)/np.mean(NS42.z_enc),0,0.05,color='k',linewidth=2)
ax1.axhline(np.mean(NS42.z_if)/np.mean(NS42.z_enc),0,0.05,color='k',linewidth=2)
ax2.axhline(np.mean(NS42.z_ig)/np.mean(NS42.z_enc),0,0.05,color='k',linewidth=2)
ax2.axhline(np.mean(NS42.z_if)/np.mean(NS42.z_enc),0,0.05,color='k',linewidth=2)
ax3.axhline(np.mean(S20.z_ig)/np.mean(S20.z_enc),0,0.05,color='k',linewidth=2)
ax3.axhline(np.mean(S20.z_if)/np.mean(S20.z_enc),0,0.05,color='k',linewidth=2)
ax4.axhline(np.mean(S20.z_ig)/np.mean(S20.z_enc),0,0.05,color='k',linewidth=2)
ax4.axhline(np.mean(S20.z_if)/np.mean(S20.z_enc),0,0.05,color='k',linewidth=2)
ax1.axvline(-0.46,0,NS42.y[-1]/np.mean(NS42.z_enc),color='k',linestyle='--')
ax2.axvline(-0.46,0,NS42.y[-1]/np.mean(NS42.z_enc),color='k',linestyle='--')
ax3.axvline(-0.91,0,S20.y[-1]/np.mean(S20.z_enc),color='k',linestyle='--')
ax4.axvline(-0.91,0,S20.y[-1]/np.mean(S20.z_enc),color='k',linestyle='--')
plt.tight_layout(rect=[0,0.08,1,1])
plt.savefig(opath+'s1wmeandiff_cond_threshold_S20_S0_vort_timmean.pdf')
plt.show()


f, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,sharey='row',sharex='col',figsize=(10,10))
ax1.tick_params(bottom=True,top=True,left=True,right=True)
ax2.tick_params(bottom=True,top=True,left=True,right=True)
ax3.tick_params(bottom=True,top=True,left=True,right=True)
ax4.tick_params(bottom=True,top=True,left=True,right=True)
ax1.set_ylim(1,1.5)
ax3.set_ylim(1,1.5)
ax1.set_yticks([1,1.1,1.2,1.3,1.4,1.5])
ax3.set_yticks([1,1.1,1.2,1.3,1.4,1.5])
cs1 = ax1.contourf(np.log10(vort_thresholds**2/(ceps*B0/nu_42)),NS42.y/np.mean(NS42.z_enc),np.sqrt(vort_nonturb_s1_var_NS42_timmean)/(N**2*L0),cmap=imola_map,levels=np.arange(0,2.1,0.1))
cs2 = ax2.contourf(np.log10(vort_thresholds**2/(ceps*B0/nu_42)),NS42.y/np.mean(NS42.z_enc),np.sqrt(vort_turb_s1_var_NS42_timmean)/(N**2*L0),cmap=imola_map,levels=np.arange(0,2.1,0.1))
cs3 = ax3.contourf(np.log10(vort_thresholds**2/(ceps*B0/nu_42)),S20.y/np.mean(S20.z_enc),np.sqrt(vort_nonturb_s1_var_S20_timmean)/(N**2*L0),cmap=imola_map,levels=np.arange(0,2.1,0.1))
cs4 = ax4.contourf(np.log10(vort_thresholds**2/(ceps*B0/nu_42)),S20.y/np.mean(S20.z_enc),np.sqrt(vort_turb_s1_var_S20_timmean)/(N**2*L0),cmap=imola_map,levels=np.arange(0,2.1,0.1))
cbar_ax = f.add_axes([0.15,0.1,0.8,0.03])
cbar = f.colorbar(cs1,cax=cbar_ax,orientation='horizontal')
ax1.set_ylabel(r'$z/z_\mathrm{enc}$')
ax3.set_ylabel(r'$z/z_\mathrm{enc}$')
ax3.set_xlabel(r'$\log_{10}(\omega_\mathrm{th}^2/\omega_0^2)$')
ax4.set_xlabel(r'$\log_{10}(\omega_\mathrm{th}^2/\omega_0^2)$')
ax1.set_title(r'(a) $Fr_0=0$, $(b_\mathrm{rms})_\mathrm{NT}/(N^2L_0)$',loc='left',fontsize=20)
ax2.set_title(r'(b) $Fr_0=0$, $(b_\mathrm{rms})_\mathrm{T}/(N^2L_0)$',loc='left',fontsize=20)
ax3.set_title(r'(c) $Fr_0=20$, $(b_\mathrm{rms})_\mathrm{NT}/(N^2L_0)$',loc='left',fontsize=20)
ax4.set_title(r'(d) $Fr_0=20$, $(b_\mathrm{rms})_\mathrm{T}/(N^2L_0)$',loc='left',fontsize=20)
ax1.axhline(np.mean(NS42.z_ig)/np.mean(NS42.z_enc),0,0.05,color='k',linewidth=2)
ax1.axhline(np.mean(NS42.z_if)/np.mean(NS42.z_enc),0,0.05,color='k',linewidth=2)
ax2.axhline(np.mean(NS42.z_ig)/np.mean(NS42.z_enc),0,0.05,color='k',linewidth=2)
ax2.axhline(np.mean(NS42.z_if)/np.mean(NS42.z_enc),0,0.05,color='k',linewidth=2)
ax3.axhline(np.mean(S20.z_ig)/np.mean(S20.z_enc),0,0.05,color='k',linewidth=2)
ax3.axhline(np.mean(S20.z_if)/np.mean(S20.z_enc),0,0.05,color='k',linewidth=2)
ax4.axhline(np.mean(S20.z_ig)/np.mean(S20.z_enc),0,0.05,color='k',linewidth=2)
ax4.axhline(np.mean(S20.z_if)/np.mean(S20.z_enc),0,0.05,color='k',linewidth=2)
ax1.axvline(-0.46,0,NS42.y[-1]/np.mean(NS42.z_enc),color='k',linestyle='--')
ax2.axvline(-0.46,0,NS42.y[-1]/np.mean(NS42.z_enc),color='k',linestyle='--')
ax3.axvline(-0.91,0,S20.y[-1]/np.mean(S20.z_enc),color='k',linestyle='--')
ax4.axvline(-0.91,0,S20.y[-1]/np.mean(S20.z_enc),color='k',linestyle='--')
plt.tight_layout(rect=[0,0.15,1,1],h_pad=2)
plt.savefig(opath+'s1rms_cond_threshold_S20_S0_vort_timmean.pdf')
plt.show()


f, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,sharey='row',sharex='col',figsize=(10,10))
ax1.tick_params(bottom=True,top=True,left=True,right=True)
ax2.tick_params(bottom=True,top=True,left=True,right=True)
ax3.tick_params(bottom=True,top=True,left=True,right=True)
ax4.tick_params(bottom=True,top=True,left=True,right=True)
ax1.set_ylim(1,1.5)
ax3.set_ylim(1,1.5)
ax1.set_yticks([1,1.1,1.2,1.3,1.4,1.5])
ax3.set_yticks([1,1.1,1.2,1.3,1.4,1.5])
cs1 = ax1.contourf(np.log10(vort_thresholds**2/(ceps*B0/nu_42)),NS42.y/np.mean(NS42.z_enc),np.sqrt(vort_nonturb_w_var_NS42_timmean)/(np.mean(NS42.z_enc)*B0)**(1./3.),cmap=imola_map,levels=np.arange(0,0.65,0.05),extend='max')
cs2 = ax2.contourf(np.log10(vort_thresholds**2/(ceps*B0/nu_42)),NS42.y/np.mean(NS42.z_enc),np.sqrt(vort_turb_w_var_NS42_timmean)/(np.mean(NS42.z_enc)*B0)**(1./3.),cmap=imola_map,levels=np.arange(0,0.65,0.05),extend='max')
cs3 = ax3.contourf(np.log10(vort_thresholds**2/(ceps*B0/nu_42)),S20.y/np.mean(S20.z_enc),np.sqrt(vort_nonturb_w_var_S20_timmean)/(np.mean(S20.z_enc)*B0)**(1./3.),cmap=imola_map,levels=np.arange(0,0.65,0.05),extend='max')
cs4 = ax4.contourf(np.log10(vort_thresholds**2/(ceps*B0/nu_42)),S20.y/np.mean(S20.z_enc),np.sqrt(vort_turb_w_var_S20_timmean)/(np.mean(S20.z_enc)*B0)**(1./3.),cmap=imola_map,levels=np.arange(0,0.65,0.05),extend='max')
cbar_ax = f.add_axes([0.15,0.1,0.8,0.03])
cbar = f.colorbar(cs1,cax=cbar_ax,orientation='horizontal')
ax1.set_ylabel(r'$z/z_\mathrm{enc}$')
ax3.set_ylabel(r'$z/z_\mathrm{enc}$')
ax3.set_xlabel(r'$\log_{10}(\omega_\mathrm{th}^2/\omega_0^2)$')
ax4.set_xlabel(r'$\log_{10}(\omega_\mathrm{th}^2/\omega_0^2)$')
ax1.set_title(r'(a) $Fr_0=0$, $(w_\mathrm{rms})_\mathrm{NT}/w_*$',loc='left',fontsize=20)
ax2.set_title(r'(b) $Fr_0=0$, $(w_\mathrm{rms})_\mathrm{T}/w_*$',loc='left',fontsize=20)
ax3.set_title(r'(c) $Fr_0=20$, $(w_\mathrm{rms})_\mathrm{NT}/w_*$',loc='left',fontsize=20)
ax4.set_title(r'(d) $Fr_0=20$, $(w_\mathrm{rms})_\mathrm{T}/w_*$',loc='left',fontsize=20)
ax1.axhline(np.mean(NS42.z_ig)/np.mean(NS42.z_enc),0,0.05,color='k',linewidth=2)
ax1.axhline(np.mean(NS42.z_if)/np.mean(NS42.z_enc),0,0.05,color='k',linewidth=2)
ax2.axhline(np.mean(NS42.z_ig)/np.mean(NS42.z_enc),0,0.05,color='k',linewidth=2)
ax2.axhline(np.mean(NS42.z_if)/np.mean(NS42.z_enc),0,0.05,color='k',linewidth=2)
ax3.axhline(np.mean(S20.z_ig)/np.mean(S20.z_enc),0,0.05,color='k',linewidth=2)
ax3.axhline(np.mean(S20.z_if)/np.mean(S20.z_enc),0,0.05,color='k',linewidth=2)
ax4.axhline(np.mean(S20.z_ig)/np.mean(S20.z_enc),0,0.05,color='k',linewidth=2)
ax4.axhline(np.mean(S20.z_if)/np.mean(S20.z_enc),0,0.05,color='k',linewidth=2)
ax1.axvline(-0.46,0,NS42.y[-1]/np.mean(NS42.z_enc),color='k',linestyle='--')
ax2.axvline(-0.46,0,NS42.y[-1]/np.mean(NS42.z_enc),color='k',linestyle='--')
ax3.axvline(-0.91,0,S20.y[-1]/np.mean(S20.z_enc),color='k',linestyle='--')
ax4.axvline(-0.91,0,S20.y[-1]/np.mean(S20.z_enc),color='k',linestyle='--')
plt.tight_layout(rect=[0,0.15,1,1],h_pad=2)
plt.savefig(opath+'wrms_cond_threshold_S20_S0_vort_timmean.pdf')
plt.show()


f, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,sharey='row',sharex='col',figsize=(10,10))
ax1.tick_params(bottom=True,top=True,left=True,right=True)
ax2.tick_params(bottom=True,top=True,left=True,right=True)
ax3.tick_params(bottom=True,top=True,left=True,right=True)
ax4.tick_params(bottom=True,top=True,left=True,right=True)
ax1.set_ylim(1,1.5)
ax3.set_ylim(1,1.5)
ax1.set_yticks([1,1.1,1.2,1.3,1.4,1.5])
ax3.set_yticks([1,1.1,1.2,1.3,1.4,1.5])
cs1 = ax1.contourf(np.log10(vort_thresholds**2/(ceps*B0/nu_42)),NS42.y/np.mean(NS42.z_enc),np.mean(vort_nonturb_s1_flux_NS42-(vort_nonturb_s1_mean_NS42*vort_nonturb_w_mean_NS42),axis=0)/B0,cmap=vik_map,levels=np.arange(-0.3,0.32,0.02),extend='min')
cs2 = ax2.contourf(np.log10(vort_thresholds**2/(ceps*B0/nu_42)),NS42.y/np.mean(NS42.z_enc),np.mean(vort_turb_s1_flux_NS42-(vort_turb_s1_mean_NS42*vort_turb_w_mean_NS42),axis=0)/B0,cmap=vik_map,levels=np.arange(-0.3,0.32,0.02),extend='min')
cs3 = ax3.contourf(np.log10(vort_thresholds**2/(ceps*B0/nu_42)),S20.y/np.mean(S20.z_enc),np.mean(vort_nonturb_s1_flux_S20-(vort_nonturb_s1_mean_S20*vort_nonturb_w_mean_S20),axis=0)/B0,cmap=vik_map,levels=np.arange(-0.3,0.32,0.02),extend='min')
cs4 = ax4.contourf(np.log10(vort_thresholds**2/(ceps*B0/nu_42)),S20.y/np.mean(S20.z_enc),np.mean(vort_turb_s1_flux_S20-(vort_turb_s1_mean_S20*vort_turb_w_mean_S20),axis=0)/B0,cmap=vik_map,levels=np.arange(-0.3,0.32,0.02),extend='min')
cbar_ax = f.add_axes([0.15,0.1,0.8,0.03])
cbar = f.colorbar(cs1,cax=cbar_ax,orientation='horizontal')
ax1.set_ylabel(r'$z/z_\mathrm{enc}$')
ax3.set_ylabel(r'$z/z_\mathrm{enc}$')
ax3.set_xlabel(r'$\log_{10}(\omega_\mathrm{th}^2/\omega_0^2)$')
ax4.set_xlabel(r'$\log_{10}(\omega_\mathrm{th}^2/\omega_0^2)$')
ax1.set_title(r'(a) $Fr_0=0$, $\langle b^\prime w^\prime \rangle_\mathrm{NT}/B_0$',loc='left',fontsize=20)
ax2.set_title(r'(b) $Fr_0=0$, $\langle b^\prime w^\prime \rangle_\mathrm{T}/B_0$',loc='left',fontsize=20)
ax3.set_title(r'(c) $Fr_0=20$, $\langle b^\prime w^\prime \rangle_\mathrm{NT}/B_0$',loc='left',fontsize=20)
ax4.set_title(r'(d) $Fr_0=20$, $\langle b^\prime w^\prime \rangle_\mathrm{T}/B_0$',loc='left',fontsize=20)
ax1.axhline(np.mean(NS42.z_ig)/np.mean(NS42.z_enc),0,0.05,color='k',linewidth=2)
ax1.axhline(np.mean(NS42.z_if)/np.mean(NS42.z_enc),0,0.05,color='k',linewidth=2)
ax2.axhline(np.mean(NS42.z_ig)/np.mean(NS42.z_enc),0,0.05,color='k',linewidth=2)
ax2.axhline(np.mean(NS42.z_if)/np.mean(NS42.z_enc),0,0.05,color='k',linewidth=2)
ax3.axhline(np.mean(S20.z_ig)/np.mean(S20.z_enc),0,0.05,color='k',linewidth=2)
ax3.axhline(np.mean(S20.z_if)/np.mean(S20.z_enc),0,0.05,color='k',linewidth=2)
ax4.axhline(np.mean(S20.z_ig)/np.mean(S20.z_enc),0,0.05,color='k',linewidth=2)
ax4.axhline(np.mean(S20.z_if)/np.mean(S20.z_enc),0,0.05,color='k',linewidth=2)
ax1.axvline(-0.46,0,NS42.y[-1]/np.mean(NS42.z_enc),color='k',linestyle='--')
ax2.axvline(-0.46,0,NS42.y[-1]/np.mean(NS42.z_enc),color='k',linestyle='--')
ax3.axvline(-0.91,0,S20.y[-1]/np.mean(S20.z_enc),color='k',linestyle='--')
ax4.axvline(-0.91,0,S20.y[-1]/np.mean(S20.z_enc),color='k',linestyle='--')
plt.tight_layout(rect=[0,0.15,1,1],h_pad=2)
plt.savefig(opath+'s1_vflux_cond_threshold_S20_S0_vort_timmean.pdf')
plt.show()

f, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,sharey='row',sharex='col',figsize=(10,10))
ax1.tick_params(bottom=True,top=True,left=True,right=True)
ax2.tick_params(bottom=True,top=True,left=True,right=True)
ax3.tick_params(bottom=True,top=True,left=True,right=True)
ax4.tick_params(bottom=True,top=True,left=True,right=True)
ax1.set_ylim(1,1.5)
ax3.set_ylim(1,1.5)
ax1.set_yticks([1,1.1,1.2,1.3,1.4,1.5])
ax3.set_yticks([1,1.1,1.2,1.3,1.4,1.5])
ax1.set_xlim(np.log10(pvthresholds[0]/(cb*ceps*N**6*42**2*(np.mean(NS42.z_enc)/L0)**(-4./3.))),np.log10(pvthresholds[-1]/(cb*ceps*N**6*42**2*(np.mean(S20.z_enc)/L0)**(-4./3.))))
ax2.set_xlim(np.log10(pvthresholds[0]/(cb*ceps*N**6*42**2*(np.mean(NS42.z_enc)/L0)**(-4./3.))),np.log10(pvthresholds[-1]/(cb*ceps*N**6*42**2*(np.mean(S20.z_enc)/L0)**(-4./3.))))
ax3.set_xlim(np.log10(pvthresholds[0]/(cb*ceps*N**6*42**2*(np.mean(NS42.z_enc)/L0)**(-4./3.))),np.log10(pvthresholds[-1]/(cb*ceps*N**6*42**2*(np.mean(S20.z_enc)/L0)**(-4./3.))))
ax4.set_xlim(np.log10(pvthresholds[0]/(cb*ceps*N**6*42**2*(np.mean(NS42.z_enc)/L0)**(-4./3.))),np.log10(pvthresholds[-1]/(cb*ceps*N**6*42**2*(np.mean(S20.z_enc)/L0)**(-4./3.))))
ax3.set_xticks([-4,-3,-2,-1,0])
ax4.set_xticks([-4,-3,-2,-1,0])
cs1 = ax1.contourf(np.log10(pvthresholds**2/(cb*ceps*N**6*42**2*(np.mean(NS42.z_enc)/L0)**(-4./3.))),NS42.y/np.mean(NS42.z_enc),np.mean(pv_nonturb_s1_flux_NS42-(pv_nonturb_s1_mean_NS42*pv_nonturb_w_mean_NS42),axis=0)/B0,cmap=vik_map,levels=np.arange(-0.3,0.32,0.02),extend='min')
cs2 = ax2.contourf(np.log10(pvthresholds**2/(cb*ceps*N**6*42**2*(np.mean(NS42.z_enc)/L0)**(-4./3.))),NS42.y/np.mean(NS42.z_enc),np.mean(pv_turb_s1_flux_NS42-(pv_turb_s1_mean_NS42*pv_turb_w_mean_NS42),axis=0)/B0,cmap=vik_map,levels=np.arange(-0.3,0.32,0.02),extend='min')
cs3 = ax3.contourf(np.log10(pvthresholds**2/(cb*ceps*N**6*42**2*(np.mean(NS42.z_enc)/L0)**(-4./3.))),S20.y/np.mean(S20.z_enc),np.mean(pv_nonturb_s1_flux_S20-(pv_nonturb_s1_mean_S20*pv_nonturb_w_mean_S20),axis=0)/B0,cmap=vik_map,levels=np.arange(-0.3,0.32,0.02),extend='min')
cs4 = ax4.contourf(np.log10(pvthresholds**2/(cb*ceps*N**6*42**2*(np.mean(NS42.z_enc)/L0)**(-4./3.))),S20.y/np.mean(S20.z_enc),np.mean(pv_turb_s1_flux_S20-(pv_turb_s1_mean_S20*pv_turb_w_mean_S20),axis=0)/B0,cmap=vik_map,levels=np.arange(-0.3,0.32,0.02),extend='min')
cbar_ax = f.add_axes([0.15,0.1,0.8,0.03])
cbar = f.colorbar(cs1,cax=cbar_ax,orientation='horizontal')
ax1.set_ylabel(r'$z/z_\mathrm{enc}$')
ax3.set_ylabel(r'$z/z_\mathrm{enc}$')
ax3.set_xlabel(r'$\log_{10}(\Pi_\mathrm{th}^2/\Pi_0^2)$')
ax4.set_xlabel(r'$\log_{10}(\Pi_\mathrm{th}^2/\Pi_0^2)$')
ax1.set_title(r'(a) $Fr_0=0$, $\langle b^\prime w^\prime \rangle_\mathrm{NT}/B_0$',loc='left',fontsize=20)
ax2.set_title(r'(b) $Fr_0=0$, $\langle b^\prime w^\prime \rangle_\mathrm{T}/B_0$',loc='left',fontsize=20)
ax3.set_title(r'(c) $Fr_0=20$, $\langle b^\prime w^\prime \rangle_\mathrm{NT}/B_0$',loc='left',fontsize=20)
ax4.set_title(r'(d) $Fr_0=20$, $\langle b^\prime w^\prime \rangle_\mathrm{T}/B_0$',loc='left',fontsize=20)
ax1.axhline(np.mean(NS42.z_ig)/np.mean(NS42.z_enc),0,0.05,color='k',linewidth=2)
ax1.axhline(np.mean(NS42.z_if)/np.mean(NS42.z_enc),0,0.05,color='k',linewidth=2)
ax2.axhline(np.mean(NS42.z_ig)/np.mean(NS42.z_enc),0,0.05,color='k',linewidth=2)
ax2.axhline(np.mean(NS42.z_if)/np.mean(NS42.z_enc),0,0.05,color='k',linewidth=2)
ax3.axhline(np.mean(S20.z_ig)/np.mean(S20.z_enc),0,0.05,color='k',linewidth=2)
ax3.axhline(np.mean(S20.z_if)/np.mean(S20.z_enc),0,0.05,color='k',linewidth=2)
ax4.axhline(np.mean(S20.z_ig)/np.mean(S20.z_enc),0,0.05,color='k',linewidth=2)
ax4.axhline(np.mean(S20.z_if)/np.mean(S20.z_enc),0,0.05,color='k',linewidth=2)
ax1.axvline(-1.99,0,NS42.y[-1]/np.mean(NS42.z_enc),color='k',linestyle='--')
ax2.axvline(-1.99,0,NS42.y[-1]/np.mean(NS42.z_enc),color='k',linestyle='--')
ax3.axvline(-2.647,0,S20.y[-1]/np.mean(S20.z_enc),color='k',linestyle='--')
ax4.axvline(-2.647,0,S20.y[-1]/np.mean(S20.z_enc),color='k',linestyle='--')
plt.tight_layout(rect=[0,0.15,1,1],h_pad=2)
plt.savefig(opath+'s1_vflux_cond_threshold_S20_S0_pv_timmean.pdf')
plt.show()
    

####################################
# Tests
####################################

plt.figure(figsize=(5,5))
plt.grid(True)
plt.xticks([-2,-1,0,1])
plt.ylim(0,0.3)
plt.plot(np.log10(vort_thresholds**2/(ceps*B0/nu_42)),-vort_turb_s1_flux_zif_NS42_timmean/B0,c=blues(0.5),label=r'$Fr_0=0$')
plt.plot(np.log10(vort_thresholds**2/(ceps*B0/nu_42)),-vort_turb_s1_flux_zif_S20_timmean/B0,c=blues(0.7),label=r'$Fr_0=20$')
plt.xlabel(r'$\log_{10}(\omega_\mathrm{th}^2/\omega_0^2)$')
plt.ylabel(r'$-(\langle b^\prime w^\prime \rangle_\mathrm{T})_{z_{i,f}}/B_0$')
plt.legend(loc='best',fontsize=20)
plt.tight_layout()
plt.savefig(opath+'s1_vflux_zif_cond_threshold_S20_S0_vort_timmean.pdf')
plt.show()

f, ((ax1,ax2),(ax3,ax4),(ax5,ax6)) = plt.subplots(3,2,sharey='row',sharex='col',figsize=(10,15))
ax1.tick_params(bottom=True,top=True,left=True,right=True)
ax2.tick_params(bottom=True,top=True,left=True,right=True)
ax1.set_ylim(1,1.5)
ax3.set_ylim(1,1.5)
ax5.set_ylim(1,1.5)
ax1.set_xlim(np.log10(vort_thresholds[0]**2/(ceps*B0/nu_42)),np.log10(vort_thresholds[-1]**2/(ceps*B0/nu_42)))
ax2.set_xlim(np.log10(vort_thresholds[0]**2/(ceps*B0/nu_42)),np.log10(vort_thresholds[-1]**2/(ceps*B0/nu_42)))
cs1 = ax1.contourf(np.log10(vort_thresholds**2/(ceps*B0/nu_42)),S20.y/S20.z_enc[0],np.sqrt(vort_nonturb_s1_var_S20[0])/(N**2*L0),cmap=imola_map,levels=np.arange(0,2.1,0.1),extend='max')
cs2 = ax2.contourf(np.log10(vort_thresholds**2/(ceps*B0/nu_42)),S20.y/S20.z_enc[0],np.sqrt(vort_turb_s1_var_S20[0])/(N**2*L0),cmap=imola_map,levels=np.arange(0,2.1,0.1),extend='max')
cs3 = ax3.contourf(np.log10(vort_thresholds**2/(ceps*B0/nu_42)),S20.y/S20.z_enc[3],np.sqrt(vort_nonturb_s1_var_S20[3])/(N**2*L0),cmap=imola_map,levels=np.arange(0,2.1,0.1),extend='max')
cs4 = ax4.contourf(np.log10(vort_thresholds**2/(ceps*B0/nu_42)),S20.y/S20.z_enc[3],np.sqrt(vort_turb_s1_var_S20[3])/(N**2*L0),cmap=imola_map,levels=np.arange(0,2.1,0.1),extend='max')
cs5 = ax5.contourf(np.log10(vort_thresholds**2/(ceps*B0/nu_42)),S20.y/S20.z_enc[-1],np.sqrt(vort_nonturb_s1_var_S20[-1])/(N**2*L0),cmap=imola_map,levels=np.arange(0,2.1,0.1),extend='max')
cs6 = ax6.contourf(np.log10(vort_thresholds**2/(ceps*B0/nu_42)),S20.y/S20.z_enc[-1],np.sqrt(vort_turb_s1_var_S20[-1])/(N**2*L0),cmap=imola_map,levels=np.arange(0,2.1,0.1),extend='max')
cbar_ax = f.add_axes([0.15,0.1,0.8,0.03])
cbar = f.colorbar(cs1,cax=cbar_ax,orientation='horizontal')
ax1.set_ylabel(r'$z/z_\mathrm{enc}$')
ax5.set_xlabel(r'$\log_{10}(\omega_\mathrm{th}^2/\omega_0^2)$')
ax6.set_xlabel(r'$\log_{10}(\omega_\mathrm{th}^2/\omega_0^2)$')
ax1.set_title(r'(a) $z_\mathrm{enc}/L_0=15$',loc='left',fontsize=20)
ax2.set_title(r'(b) $z_\mathrm{enc}/L_0=15$',loc='left',fontsize=20)
ax3.set_title(r'(c) $z_\mathrm{enc}/L_0=18$',loc='left',fontsize=20)
ax4.set_title(r'(d) $z_\mathrm{enc}/L_0=18$',loc='left',fontsize=20)
ax5.set_title(r'(e) $z_\mathrm{enc}/L_0=21$',loc='left',fontsize=20)
ax6.set_title(r'(f) $z_\mathrm{enc}/L_0=21$',loc='left',fontsize=20)
ax1.axhline(S20.z_ig[0]/S20.z_enc[0],0.95,1,color='k',linewidth=2)
ax1.axhline(S20.z_if[0]/S20.z_enc[0],0.95,1,color='k',linewidth=2)
ax2.axhline(S20.z_ig[0]/S20.z_enc[0],0.95,1,color='k',linewidth=2)
ax2.axhline(S20.z_if[0]/S20.z_enc[0],0.95,1,color='k',linewidth=2)
ax3.axhline(S20.z_ig[3]/S20.z_enc[3],0.95,1,color='k',linewidth=2)
ax3.axhline(S20.z_if[3]/S20.z_enc[3],0.95,1,color='k',linewidth=2)
ax4.axhline(S20.z_ig[3]/S20.z_enc[3],0.95,1,color='k',linewidth=2)
ax4.axhline(S20.z_if[3]/S20.z_enc[3],0.95,1,color='k',linewidth=2)
ax5.axhline(S20.z_ig[-1]/S20.z_enc[-1],0.95,1,color='k',linewidth=2)
ax5.axhline(S20.z_if[-1]/S20.z_enc[-1],0.95,1,color='k',linewidth=2)
ax6.axhline(S20.z_ig[-1]/S20.z_enc[-1],0.95,1,color='k',linewidth=2)
ax6.axhline(S20.z_if[-1]/S20.z_enc[-1],0.95,1,color='k',linewidth=2)
plt.tight_layout(rect=[0,0.15,1,1],h_pad=2)
plt.savefig(opath+'s1rms_cond_threshold_S20_S0_vort_N2L0_times.pdf')
plt.show()

f, (ax1,ax2) = plt.subplots(1,2,sharey='row',figsize=(10,5))
ax1.tick_params(bottom=True,top=True,left=True,right=True)
ax2.tick_params(bottom=True,top=True,left=True,right=True)
ax1.set_ylim(1,1.5)
ax1.set_xlim(np.log10(vort_thresholds[0]**2/(ceps*B0/nu_42)),np.log10(vort_thresholds[-1]**2/(ceps*B0/nu_42)))
ax2.set_xlim(np.log10(vort_thresholds[0]**2/(ceps*B0/nu_42)),np.log10(vort_thresholds[-1]**2/(ceps*B0/nu_42)))
cs1 = ax1.contourf(np.log10(vort_thresholds**2/(ceps*B0/nu_42)),S20.y/np.mean(S20.z_enc),np.sqrt(vort_nonturb_s1_var_S20_timmean)/np.mean(S20.b_enc),cmap=imola_map,levels=np.arange(0,0.095,0.005))
cs2 = ax2.contourf(np.log10(vort_thresholds**2/(ceps*B0/nu_42)),S20.y/np.mean(S20.z_enc),np.sqrt(vort_turb_s1_var_S20_timmean)/np.mean(S20.b_enc),cmap=imola_map,levels=np.arange(0,0.095,0.005))
cbar_ax = f.add_axes([0.15,0.1,0.8,0.03])
cbar = f.colorbar(cs1,cax=cbar_ax,orientation='horizontal')
ax1.set_ylabel(r'$z/z_\mathrm{enc}$')
ax1.set_xlabel(r'$\log_{10}(\omega_\mathrm{th}^2/\omega_0^2)$')
ax2.set_xlabel(r'$\log_{10}(\omega_\mathrm{th}^2/\omega_0^2)$')
ax1.set_title(r'(a) $Fr_0=20$, $(b_\mathrm{rms})_\mathrm{NT}/b_\mathrm{enc}$',loc='left',fontsize=20)
ax2.set_title(r'(b) $Fr_0=20$, $(b_\mathrm{rms})_\mathrm{T}/b_\mathrm{enc}$',loc='left',fontsize=20)
ax1.axhline(np.mean(S20.z_ig)/np.mean(S20.z_enc),0.95,1,color='k',linewidth=2)
ax1.axhline(np.mean(S20.z_if)/np.mean(S20.z_enc),0.95,1,color='k',linewidth=2)
ax2.axhline(np.mean(S20.z_ig)/np.mean(S20.z_enc),0.95,1,color='k',linewidth=2)
ax2.axhline(np.mean(S20.z_if)/np.mean(S20.z_enc),0.95,1,color='k',linewidth=2)
ax1.axvline(-0.8,0,S20.y[-1]/np.mean(S20.z_enc),color='k',linestyle='--')
ax2.axvline(-0.8,0,S20.y[-1]/np.mean(S20.z_enc),color='k',linestyle='--')
plt.tight_layout(rect=[0,0.15,1,1],h_pad=2)
plt.savefig(opath+'s1rms_cond_threshold_S20_S0_vort_timmean_benc.pdf')
plt.show()

f, ((ax1,ax2),(ax3,ax4),(ax5,ax6)) = plt.subplots(3,2,sharey='row',sharex='col',figsize=(10,15))
ax1.tick_params(bottom=True,top=True,left=True,right=True)
ax2.tick_params(bottom=True,top=True,left=True,right=True)
ax1.set_ylim(1,1.5)
ax3.set_ylim(1,1.5)
ax5.set_ylim(1,1.5)
ax1.set_xlim(np.log10(vort_thresholds[0]**2/(ceps*B0/nu_42)),np.log10(vort_thresholds[-1]**2/(ceps*B0/nu_42)))
ax2.set_xlim(np.log10(vort_thresholds[0]**2/(ceps*B0/nu_42)),np.log10(vort_thresholds[-1]**2/(ceps*B0/nu_42)))
cs1 = ax1.contourf(np.log10(vort_thresholds**2/(ceps*B0/nu_42)),S20.y/S20.z_enc[0],np.sqrt(vort_nonturb_s1_var_S20[0])/S20.b_enc[0],cmap=imola_map,levels=np.arange(0,0.105,0.005),extend='max')
cs2 = ax2.contourf(np.log10(vort_thresholds**2/(ceps*B0/nu_42)),S20.y/S20.z_enc[0],np.sqrt(vort_turb_s1_var_S20[0])/S20.b_enc[0],cmap=imola_map,levels=np.arange(0,0.105,0.005),extend='max')
cs3 = ax3.contourf(np.log10(vort_thresholds**2/(ceps*B0/nu_42)),S20.y/S20.z_enc[3],np.sqrt(vort_nonturb_s1_var_S20[3])/S20.b_enc[3],cmap=imola_map,levels=np.arange(0,0.105,0.005),extend='max')
cs4 = ax4.contourf(np.log10(vort_thresholds**2/(ceps*B0/nu_42)),S20.y/S20.z_enc[3],np.sqrt(vort_turb_s1_var_S20[3])/S20.b_enc[3],cmap=imola_map,levels=np.arange(0,0.105,0.005),extend='max')
cs5 = ax5.contourf(np.log10(vort_thresholds**2/(ceps*B0/nu_42)),S20.y/S20.z_enc[-1],np.sqrt(vort_nonturb_s1_var_S20[-1])/S20.b_enc[-1],cmap=imola_map,levels=np.arange(0,0.105,0.005),extend='max')
cs6 = ax6.contourf(np.log10(vort_thresholds**2/(ceps*B0/nu_42)),S20.y/S20.z_enc[-1],np.sqrt(vort_turb_s1_var_S20[-1])/S20.b_enc[-1],cmap=imola_map,levels=np.arange(0,0.105,0.005),extend='max')
cbar_ax = f.add_axes([0.15,0.1,0.8,0.03])
cbar = f.colorbar(cs1,cax=cbar_ax,orientation='horizontal')
ax1.set_ylabel(r'$z/z_\mathrm{enc}$')
ax5.set_xlabel(r'$\log_{10}(\omega_\mathrm{th}^2/\omega_0^2)$')
ax6.set_xlabel(r'$\log_{10}(\omega_\mathrm{th}^2/\omega_0^2)$')
ax1.set_title(r'(a) $z_\mathrm{enc}/L_0=15$',loc='left',fontsize=20)
ax2.set_title(r'(b) $z_\mathrm{enc}/L_0=15$',loc='left',fontsize=20)
ax3.set_title(r'(c) $z_\mathrm{enc}/L_0=18$',loc='left',fontsize=20)
ax4.set_title(r'(d) $z_\mathrm{enc}/L_0=18$',loc='left',fontsize=20)
ax5.set_title(r'(e) $z_\mathrm{enc}/L_0=21$',loc='left',fontsize=20)
ax6.set_title(r'(f) $z_\mathrm{enc}/L_0=21$',loc='left',fontsize=20)
ax1.axhline(S20.z_ig[0]/S20.z_enc[0],0.95,1,color='k',linewidth=2)
ax1.axhline(S20.z_if[0]/S20.z_enc[0],0.95,1,color='k',linewidth=2)
ax2.axhline(S20.z_ig[0]/S20.z_enc[0],0.95,1,color='k',linewidth=2)
ax2.axhline(S20.z_if[0]/S20.z_enc[0],0.95,1,color='k',linewidth=2)
ax3.axhline(S20.z_ig[3]/S20.z_enc[3],0.95,1,color='k',linewidth=2)
ax3.axhline(S20.z_if[3]/S20.z_enc[3],0.95,1,color='k',linewidth=2)
ax4.axhline(S20.z_ig[3]/S20.z_enc[3],0.95,1,color='k',linewidth=2)
ax4.axhline(S20.z_if[3]/S20.z_enc[3],0.95,1,color='k',linewidth=2)
ax5.axhline(S20.z_ig[-1]/S20.z_enc[-1],0.95,1,color='k',linewidth=2)
ax5.axhline(S20.z_if[-1]/S20.z_enc[-1],0.95,1,color='k',linewidth=2)
ax6.axhline(S20.z_ig[-1]/S20.z_enc[-1],0.95,1,color='k',linewidth=2)
ax6.axhline(S20.z_if[-1]/S20.z_enc[-1],0.95,1,color='k',linewidth=2)
plt.tight_layout(rect=[0,0.15,1,1],h_pad=2)
plt.savefig(opath+'s1rms_cond_threshold_S20_S0_vort_benc_times.pdf')
plt.show()

f, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,sharey='row',sharex='col',figsize=(10,10))
ax1.tick_params(bottom=True,top=True,left=True,right=True)
ax2.tick_params(bottom=True,top=True,left=True,right=True)
ax1.set_ylim(1,1.5)
ax3.set_ylim(1,1.5)
ax1.set_xlim(np.log10(vort_thresholds[0]**2/(ceps*B0/nu_42)),np.log10(vort_thresholds[-1]**2/(ceps*B0/nu_42)))
ax2.set_xlim(np.log10(vort_thresholds[0]**2/(ceps*B0/nu_42)),np.log10(vort_thresholds[-1]**2/(ceps*B0/nu_42)))
cs1 = ax1.contourf(np.log10(vort_thresholds**2/(ceps*B0/nu_42)),NS42.y/np.mean(NS42.z_enc),np.sqrt(vort_nonturb_s1_var_NS42_timmean)/np.mean(b_delta_NS42),cmap=imola_map,levels=np.arange(0,0.65,0.05))
cs2 = ax2.contourf(np.log10(vort_thresholds**2/(ceps*B0/nu_42)),NS42.y/np.mean(NS42.z_enc),np.sqrt(vort_turb_s1_var_NS42_timmean)/np.mean(b_delta_NS42),cmap=imola_map,levels=np.arange(0,0.65,0.05))
cs3 = ax3.contourf(np.log10(vort_thresholds**2/(ceps*B0/nu_42)),S20.y/np.mean(S20.z_enc),np.sqrt(vort_nonturb_s1_var_S20_timmean)/np.mean(b_delta_S20),cmap=imola_map,levels=np.arange(0,0.65,0.05))
cs4 = ax4.contourf(np.log10(vort_thresholds**2/(ceps*B0/nu_42)),S20.y/np.mean(S20.z_enc),np.sqrt(vort_turb_s1_var_S20_timmean)/np.mean(b_delta_S20),cmap=imola_map,levels=np.arange(0,0.65,0.05))
cbar_ax = f.add_axes([0.15,0.1,0.8,0.03])
cbar = f.colorbar(cs1,cax=cbar_ax,orientation='horizontal')
ax1.set_ylabel(r'$z/z_\mathrm{enc}$')
ax3.set_ylabel(r'$z/z_\mathrm{enc}$')
ax3.set_xlabel(r'$\log_{10}(\omega_\mathrm{th}^2/\omega_0^2)$')
ax4.set_xlabel(r'$\log_{10}(\omega_\mathrm{th}^2/\omega_0^2)$')
ax1.set_title(r'(a) $Fr_0=0$, $(b_\mathrm{rms})_\mathrm{NT}/b_\delta$',loc='left',fontsize=20)
ax2.set_title(r'(b) $Fr_0=0$, $(b_\mathrm{rms})_\mathrm{T}/b_\delta$',loc='left',fontsize=20)
ax3.set_title(r'(c) $Fr_0=20$, $(b_\mathrm{rms})_\mathrm{NT}/b_\delta$',loc='left',fontsize=20)
ax4.set_title(r'(d) $Fr_0=20$, $(b_\mathrm{rms})_\mathrm{T}/b_\delta$',loc='left',fontsize=20)
ax1.axhline(np.mean(NS42.z_ig)/np.mean(NS42.z_enc),0.95,1,color='k',linewidth=2)
ax1.axhline(np.mean(NS42.z_if)/np.mean(NS42.z_enc),0.95,1,color='k',linewidth=2)
ax2.axhline(np.mean(NS42.z_ig)/np.mean(NS42.z_enc),0.95,1,color='k',linewidth=2)
ax2.axhline(np.mean(NS42.z_if)/np.mean(NS42.z_enc),0.95,1,color='k',linewidth=2)
ax3.axhline(np.mean(S20.z_ig)/np.mean(S20.z_enc),0.95,1,color='k',linewidth=2)
ax3.axhline(np.mean(S20.z_if)/np.mean(S20.z_enc),0.95,1,color='k',linewidth=2)
ax4.axhline(np.mean(S20.z_ig)/np.mean(S20.z_enc),0.95,1,color='k',linewidth=2)
ax4.axhline(np.mean(S20.z_if)/np.mean(S20.z_enc),0.95,1,color='k',linewidth=2)
ax1.axvline(-0.83,0,NS42.y[-1]/np.mean(NS42.z_enc),color='k',linestyle='--')
ax2.axvline(-0.83,0,NS42.y[-1]/np.mean(NS42.z_enc),color='k',linestyle='--')
ax3.axvline(-0.8,0,S20.y[-1]/np.mean(S20.z_enc),color='k',linestyle='--')
ax4.axvline(-0.8,0,S20.y[-1]/np.mean(S20.z_enc),color='k',linestyle='--')
plt.tight_layout(rect=[0,0.15,1,1],h_pad=2)
plt.savefig(opath+'s1rms_cond_threshold_S20_S0_vort_timmean_b_delta.pdf')
plt.show()

f, ((ax1,ax2),(ax3,ax4),(ax5,ax6)) = plt.subplots(3,2,sharey='row',sharex='col',figsize=(10,15))
ax1.tick_params(bottom=True,top=True,left=True,right=True)
ax2.tick_params(bottom=True,top=True,left=True,right=True)
ax1.set_ylim(1,1.5)
ax3.set_ylim(1,1.5)
ax5.set_ylim(1,1.5)
ax1.set_xlim(np.log10(vort_thresholds[0]**2/(ceps*B0/nu_42)),np.log10(vort_thresholds[-1]**2/(ceps*B0/nu_42)))
ax2.set_xlim(np.log10(vort_thresholds[0]**2/(ceps*B0/nu_42)),np.log10(vort_thresholds[-1]**2/(ceps*B0/nu_42)))
cs1 = ax1.contourf(np.log10(vort_thresholds**2/(ceps*B0/nu_42)),S20.y/S20.z_enc[0],np.sqrt(vort_nonturb_s1_var_S20[0])/b_delta_S20[0],cmap=imola_map,levels=np.arange(0,0.75,0.05),extend='max')
cs2 = ax2.contourf(np.log10(vort_thresholds**2/(ceps*B0/nu_42)),S20.y/S20.z_enc[0],np.sqrt(vort_turb_s1_var_S20[0])/b_delta_S20[0],cmap=imola_map,levels=np.arange(0,0.75,0.05),extend='max')
cs3 = ax3.contourf(np.log10(vort_thresholds**2/(ceps*B0/nu_42)),S20.y/S20.z_enc[3],np.sqrt(vort_nonturb_s1_var_S20[3])/b_delta_S20[3],cmap=imola_map,levels=np.arange(0,0.75,0.05),extend='max')
cs4 = ax4.contourf(np.log10(vort_thresholds**2/(ceps*B0/nu_42)),S20.y/S20.z_enc[3],np.sqrt(vort_turb_s1_var_S20[3])/b_delta_S20[3],cmap=imola_map,levels=np.arange(0,0.75,0.05),extend='max')
cs5 = ax5.contourf(np.log10(vort_thresholds**2/(ceps*B0/nu_42)),S20.y/S20.z_enc[-1],np.sqrt(vort_nonturb_s1_var_S20[-1])/b_delta_S20[-1],cmap=imola_map,levels=np.arange(0,0.75,0.05),extend='max')
cs6 = ax6.contourf(np.log10(vort_thresholds**2/(ceps*B0/nu_42)),S20.y/S20.z_enc[-1],np.sqrt(vort_turb_s1_var_S20[-1])/b_delta_S20[-1],cmap=imola_map,levels=np.arange(0,0.75,0.05),extend='max')
cbar_ax = f.add_axes([0.15,0.1,0.8,0.03])
cbar = f.colorbar(cs1,cax=cbar_ax,orientation='horizontal')
ax1.set_ylabel(r'$z/z_\mathrm{enc}$')
ax5.set_xlabel(r'$\log_{10}(\omega_\mathrm{th}^2/\omega_0^2)$')
ax6.set_xlabel(r'$\log_{10}(\omega_\mathrm{th}^2/\omega_0^2)$')
ax1.set_title(r'(a) $z_\mathrm{enc}/L_0=15$',loc='left',fontsize=20)
ax2.set_title(r'(b) $z_\mathrm{enc}/L_0=15$',loc='left',fontsize=20)
ax3.set_title(r'(c) $z_\mathrm{enc}/L_0=18$',loc='left',fontsize=20)
ax4.set_title(r'(d) $z_\mathrm{enc}/L_0=18$',loc='left',fontsize=20)
ax5.set_title(r'(e) $z_\mathrm{enc}/L_0=21$',loc='left',fontsize=20)
ax6.set_title(r'(f) $z_\mathrm{enc}/L_0=21$',loc='left',fontsize=20)
ax1.axhline(S20.z_ig[0]/S20.z_enc[0],0.95,1,color='k',linewidth=2)
ax1.axhline(S20.z_if[0]/S20.z_enc[0],0.95,1,color='k',linewidth=2)
ax2.axhline(S20.z_ig[0]/S20.z_enc[0],0.95,1,color='k',linewidth=2)
ax2.axhline(S20.z_if[0]/S20.z_enc[0],0.95,1,color='k',linewidth=2)
ax3.axhline(S20.z_ig[3]/S20.z_enc[3],0.95,1,color='k',linewidth=2)
ax3.axhline(S20.z_if[3]/S20.z_enc[3],0.95,1,color='k',linewidth=2)
ax4.axhline(S20.z_ig[3]/S20.z_enc[3],0.95,1,color='k',linewidth=2)
ax4.axhline(S20.z_if[3]/S20.z_enc[3],0.95,1,color='k',linewidth=2)
ax5.axhline(S20.z_ig[-1]/S20.z_enc[-1],0.95,1,color='k',linewidth=2)
ax5.axhline(S20.z_if[-1]/S20.z_enc[-1],0.95,1,color='k',linewidth=2)
ax6.axhline(S20.z_ig[-1]/S20.z_enc[-1],0.95,1,color='k',linewidth=2)
ax6.axhline(S20.z_if[-1]/S20.z_enc[-1],0.95,1,color='k',linewidth=2)
plt.tight_layout(rect=[0,0.15,1,1],h_pad=2)
plt.savefig(opath+'s1rms_cond_threshold_S20_S0_vort_b_delta_times.pdf')
plt.show()
