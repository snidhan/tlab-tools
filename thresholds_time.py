#!/usr/bin/python3
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
rc('axes', linewidth=1.5)
rc('lines',linewidth=2)

colourmap_path = '/home/mpim/m300551/local/ScientificColourMaps5/'
opath = '/home/mpim/m300551/Figures/JAS2020/'

blues = matplotlib.cm.get_cmap('Blues')
######################################################################
# Grid 

path_NS117 = '/scratch/local1/m300551/ForKatherine/qCBL_3D/Re117/5120x1024x5120/'
path_NS42 = '/scratch/local1/m300551/ForKatherine/qCBL_3D/Re042/2560x576x2560/'
path_NS25 = '/scratch/local1/m300551/ForKatherine/qCBL_3D/Re025/2560x512x2560/'
path_S20 = '/scratch/local1/m300551/ForKatherine/qCBL_3D/Re042/3072x960x4608-S20/'

#######################################################################
# Constants

nu_117 = 1./70000.
nu_42 = 1./25000.
nu_25 = 1./15000.
B0 = 0.005
N = np.sqrt(3)
L0 = (B0/N**3)**0.5
ceps = 0.1

#######################################################################
# Stats

NS117 = Statistics(path_NS117+'stats/pdftimes/avg111100-149000.nc')
NS42 = Statistics(path_NS42+'stats/timmean/avg36500-47500.nc')
NS25 = Statistics(path_NS25+'stats/pdftimes/avg17000-21000.nc')
S20 = Statistics(path_S20+'stats/avg66000-84000.nc')

########################################################################
# Conditional Stats

vort_thresholds = np.array([0.16,0.32,0.48,0.64,0.8,0.96,1.12,1.28,1.44,1.6,3.2,8,16])
filenames = ['0-16','0-32','0-48','0-64','0-8','0-96','1-12','1-28','1-44','1-6','3-2','8','16']
vort_turb_area_fracs_NS117 = np.zeros((NS117.t_len,NS117.y_len,len(vort_thresholds)))
vort_turb_area_fracs_NS42 = np.zeros((NS42.t_len,NS42.y_len,len(vort_thresholds)))
vort_turb_area_fracs_NS25 = np.zeros((NS25.t_len,NS25.y_len,len(vort_thresholds)))
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
vort_turb_area_frac_zif_NS42 = np.zeros((NS42.t_len,len(vort_thresholds)))
vort_turb_area_frac_zif_S20 = np.zeros((S20.t_len,len(vort_thresholds)))
vort_turb_s1_flux_zif_NS42 = np.zeros((NS42.t_len,len(vort_thresholds)))
vort_turb_s1_flux_zif_S20 = np.zeros((S20.t_len,len(vort_thresholds)))
vort_turb_s1_mean_zif_NS42 = np.zeros((NS42.t_len,len(vort_thresholds)))
vort_turb_s1_mean_zif_S20 = np.zeros((S20.t_len,len(vort_thresholds)))
vort_turb_w_mean_zif_NS42 = np.zeros((NS42.t_len,len(vort_thresholds)))
vort_turb_w_mean_zif_S20 = np.zeros((S20.t_len,len(vort_thresholds)))
vort_turb_s1_var_zif_NS42 = np.zeros((NS42.t_len,len(vort_thresholds)))
vort_turb_s1_var_zif_S20 = np.zeros((S20.t_len,len(vort_thresholds)))
vort_turb_w_var_zif_NS42 = np.zeros((NS42.t_len,len(vort_thresholds)))
vort_turb_w_var_zif_S20 = np.zeros((S20.t_len,len(vort_thresholds)))
vort_nonturb_area_frac_zif_NS42 = np.zeros((NS42.t_len,len(vort_thresholds)))
vort_nonturb_area_frac_zif_S20 = np.zeros((S20.t_len,len(vort_thresholds)))
vort_nonturb_s1_flux_zif_NS42 = np.zeros((NS42.t_len,len(vort_thresholds)))
vort_nonturb_s1_flux_zif_S20 = np.zeros((S20.t_len,len(vort_thresholds)))
vort_nonturb_s1_mean_zif_NS42 = np.zeros((NS42.t_len,len(vort_thresholds)))
vort_nonturb_s1_mean_zif_S20 = np.zeros((S20.t_len,len(vort_thresholds)))
vort_nonturb_w_mean_zif_NS42 = np.zeros((NS42.t_len,len(vort_thresholds)))
vort_nonturb_w_mean_zif_S20 = np.zeros((S20.t_len,len(vort_thresholds)))
vort_nonturb_s1_var_zif_NS42 = np.zeros((NS42.t_len,len(vort_thresholds)))
vort_nonturb_s1_var_zif_S20 = np.zeros((S20.t_len,len(vort_thresholds)))
vort_nonturb_w_var_zif_NS42 = np.zeros((NS42.t_len,len(vort_thresholds)))
vort_nonturb_w_var_zif_S20 = np.zeros((S20.t_len,len(vort_thresholds)))
for t in range(0,NS42.t_len):
    for i in range(0,len(vort_thresholds)):
        vort_gate_file = path_NS42+'stats/gate-vorticity/gate-'+filenames[i]+'/int36500-47500.nc'
        vort_gate_data = nc.Dataset(vort_gate_file,'r')
        vort_turb_areafrac = vort_gate_data.variables['Partition2'][t,:]
        vort_turb_area_frac_zif_NS42[t,i] = vort_turb_areafrac[NS42.z_if_arg[t]]
        vort_turb_file = path_NS42+'stats/gate-vorticity/gate-'+filenames[i]+'/Partition2/cavg36500-47500.nc'
        vort_turb_data = nc.Dataset(vort_turb_file,'r')
        vort_turb_s1mean = vort_turb_data.variables['Scalar1Mom1'][t,:]
        vort_turb_s1_mean_zif_NS42[t,i] = vort_turb_s1mean[NS42.z_if_arg[t]]
        vort_turb_wmean = vort_turb_data.variables['VMom1'][t,:]
        vort_turb_w_mean_zif_NS42[t,i] = vort_turb_wmean[NS42.z_if_arg[t]]
        vort_turb_s1var = vort_turb_data.variables['Scalar1Mom2'][t,:]
        vort_turb_s1_var_zif_NS42[t,i] = vort_turb_s1var[NS42.z_if_arg[t]]
        vort_turb_wvar = vort_turb_data.variables['VMom2'][t,:]
        vort_turb_w_var_zif_NS42[t,i] = vort_turb_wvar[NS42.z_if_arg[t]]
        vort_turb_flux_file = path_NS42+'stats/gate-vorticity/gate-'+filenames[i]+'/Partition2/avgMom36500-47500.nc'
        vort_turb_flux_data = nc.Dataset(vort_turb_flux_file,'r')
        vort_turb_s1flux = vort_turb_flux_data.variables['v1Mom1'][t,:]
        vort_turb_s1_flux_zif_NS42[t,i] = vort_turb_s1flux[NS42.z_if_arg[t]]
        vort_nonturb_areafrac = vort_gate_data.variables['Partition1'][t,:]
        vort_nonturb_area_frac_zif_NS42[t,i] = vort_nonturb_areafrac[NS42.z_if_arg[t]]
        vort_nonturb_file = path_NS42+'stats/gate-vorticity/gate-'+filenames[i]+'/Partition1/cavg36500-47500.nc'
        vort_nonturb_data = nc.Dataset(vort_nonturb_file,'r')
        vort_nonturb_s1mean = vort_nonturb_data.variables['Scalar1Mom1'][t,:]
        vort_nonturb_s1_mean_zif_NS42[t,i] = vort_nonturb_s1mean[NS42.z_if_arg[t]]
        vort_nonturb_wmean = vort_nonturb_data.variables['VMom1'][t,:]
        vort_nonturb_w_mean_zif_NS42[t,i] = vort_nonturb_wmean[NS42.z_if_arg[t]]
        vort_nonturb_s1var = vort_nonturb_data.variables['Scalar1Mom2'][t,:]
        vort_nonturb_s1_var_zif_NS42[t,i] = vort_nonturb_s1var[NS42.z_if_arg[t]]
        vort_nonturb_wvar = vort_nonturb_data.variables['VMom2'][t,:]
        vort_nonturb_w_var_zif_NS42[t,i] = vort_nonturb_wvar[NS42.z_if_arg[t]]
        vort_nonturb_flux_file = path_NS42+'stats/gate-vorticity/gate-'+filenames[i]+'/Partition1/avgMom36500-47500.nc'
        vort_nonturb_flux_data = nc.Dataset(vort_nonturb_flux_file,'r')
        vort_nonturb_s1flux = vort_nonturb_flux_data.variables['v1Mom1'][t,:]
        vort_nonturb_s1_flux_zif_NS42[t,i] = vort_nonturb_s1flux[NS42.z_if_arg[t]]
vort_turb_area_frac_zif_NS42_timmean = np.mean(vort_turb_area_frac_zif_NS42,axis=0)
vort_nonturb_area_frac_zif_NS42_timmean = np.mean(vort_nonturb_area_frac_zif_NS42,axis=0)
vort_turb_s1_flux_zif_NS42_timmean = np.mean(vort_turb_s1_flux_zif_NS42-(vort_turb_s1_mean_zif_NS42*vort_turb_w_mean_zif_NS42),axis=0)
vort_nonturb_s1_flux_zif_NS42_timmean = np.mean(vort_nonturb_s1_flux_zif_NS42-(vort_nonturb_s1_mean_zif_NS42*vort_nonturb_w_mean_zif_NS42),axis=0)
vort_turb_s1_var_zif_NS42_timmean = np.mean(vort_turb_s1_var_zif_NS42,axis=0)
vort_nonturb_s1_var_zif_NS42_timmean = np.mean(vort_nonturb_s1_var_zif_NS42,axis=0)
vort_turb_w_var_zif_NS42_timmean = np.mean(vort_turb_w_var_zif_NS42,axis=0)
vort_nonturb_w_var_zif_NS42_timmean = np.mean(vort_nonturb_w_var_zif_NS42,axis=0)
for t in range(0,S20.t_len):
    for i in range(0,len(vort_thresholds)):
        vort_gate_file = path_S20+'stats/gate-vorticity/gate-'+filenames[i]+'/int66000-84000.nc'
        vort_gate_data = nc.Dataset(vort_gate_file,'r')
        vort_turb_areafrac = vort_gate_data.variables['Partition2'][t,:]
        vort_turb_area_frac_zif_S20[t,i] = vort_turb_areafrac[S20.z_if_arg[t]]
        vort_turb_file = path_S20+'stats/gate-vorticity/gate-'+filenames[i]+'/Partition2/cavg66000-84000.nc'
        vort_turb_data = nc.Dataset(vort_turb_file,'r')
        vort_turb_s1mean = vort_turb_data.variables['Scalar1Mom1'][t,:]
        vort_turb_s1_mean_zif_S20[t,i] = vort_turb_s1mean[S20.z_if_arg[t]]
        vort_turb_wmean = vort_turb_data.variables['VMom1'][t,:]
        vort_turb_w_mean_zif_S20[t,i] = vort_turb_wmean[S20.z_if_arg[t]]
        vort_turb_s1var = vort_turb_data.variables['Scalar1Mom2'][t,:]
        vort_turb_s1_var_zif_S20[t,i] = vort_turb_s1var[S20.z_if_arg[t]]
        vort_turb_wvar = vort_turb_data.variables['VMom2'][t,:]
        vort_turb_w_var_zif_S20[t,i] = vort_turb_wvar[S20.z_if_arg[t]]
        vort_turb_flux_file = path_S20+'stats/gate-vorticity/gate-'+filenames[i]+'/Partition2/avgMom66000-84000.nc'
        vort_turb_flux_data = nc.Dataset(vort_turb_flux_file,'r')
        vort_turb_s1flux = vort_turb_flux_data.variables['v1Mom1'][t,:]
        vort_turb_s1_flux_zif_S20[t,i] = vort_turb_s1flux[S20.z_if_arg[t]]
        vort_nonturb_areafrac = vort_gate_data.variables['Partition1'][t,:]
        vort_nonturb_area_frac_zif_S20[t,i] = vort_nonturb_areafrac[S20.z_if_arg[t]]
        vort_nonturb_file = path_S20+'stats/gate-vorticity/gate-'+filenames[i]+'/Partition1/cavg66000-84000.nc'
        vort_nonturb_data = nc.Dataset(vort_nonturb_file,'r')
        vort_nonturb_s1mean = vort_nonturb_data.variables['Scalar1Mom1'][t,:]
        vort_nonturb_s1_mean_zif_S20[t,i] = vort_nonturb_s1mean[S20.z_if_arg[t]]
        vort_nonturb_wmean = vort_nonturb_data.variables['VMom1'][t,:]
        vort_nonturb_w_mean_zif_S20[t,i] = vort_nonturb_wmean[S20.z_if_arg[t]]
        vort_nonturb_s1var = vort_nonturb_data.variables['Scalar1Mom2'][t,:]
        vort_nonturb_s1_var_zif_S20[t,i] = vort_nonturb_s1var[S20.z_if_arg[t]]
        vort_nonturb_wvar = vort_nonturb_data.variables['VMom2'][t,:]
        vort_nonturb_w_var_zif_S20[t,i] = vort_nonturb_wvar[S20.z_if_arg[t]]
        vort_nonturb_flux_file = path_S20+'stats/gate-vorticity/gate-'+filenames[i]+'/Partition1/avgMom66000-84000.nc'
        vort_nonturb_flux_data = nc.Dataset(vort_nonturb_flux_file,'r')
        vort_nonturb_s1flux = vort_nonturb_flux_data.variables['v1Mom1'][t,:]
        vort_nonturb_s1_flux_zif_S20[t,i] = vort_nonturb_s1flux[S20.z_if_arg[t]]
vort_turb_area_frac_zif_S20_timmean = np.mean(vort_turb_area_frac_zif_S20,axis=0)
vort_nonturb_area_frac_zif_S20_timmean = np.mean(vort_nonturb_area_frac_zif_S20,axis=0)
vort_turb_s1_flux_zif_S20_timmean = np.mean(vort_turb_s1_flux_zif_S20-(vort_turb_s1_mean_zif_S20*vort_turb_w_mean_zif_S20),axis=0)
vort_nonturb_s1_flux_zif_S20_timmean = np.mean(vort_nonturb_s1_flux_zif_S20-(vort_nonturb_s1_mean_zif_S20*vort_nonturb_w_mean_zif_S20),axis=0)
vort_turb_s1_var_zif_S20_timmean = np.mean(vort_turb_s1_var_zif_S20,axis=0)
vort_nonturb_s1_var_zif_S20_timmean = np.mean(vort_nonturb_s1_var_zif_S20,axis=0)
vort_turb_w_var_zif_S20_timmean = np.mean(vort_turb_w_var_zif_S20,axis=0)
vort_nonturb_w_var_zif_S20_timmean = np.mean(vort_nonturb_w_var_zif_S20,axis=0)
for i in range(0,len(vort_thresholds)):
    vort_gate_file = path_NS117+'stats/gate-vorticity/GATE-'+filenames[i]+'/int111100-149000.nc'
    vort_gate_data = nc.Dataset(vort_gate_file,'r')
    vort_turb_area = vort_gate_data.variables['Partition2'][:,:]
    vort_turb_area_fracs_NS117[:,:,i] = vort_turb_area
    vort_turb_area_fracs_NS117_timmean = np.mean(vort_turb_area_fracs_NS117,axis=0)
    vort_gate_file = path_NS42+'stats/gate-vorticity/gate-'+filenames[i]+'/int36500-47500.nc'
    vort_gate_data = nc.Dataset(vort_gate_file,'r')
    vort_turb_area = vort_gate_data.variables['Partition2'][:,:]
    vort_turb_area_fracs_NS42[:,:,i] = vort_turb_area
    vort_turb_area_fracs_NS42_timmean = np.mean(vort_turb_area_fracs_NS42,axis=0)
    vort_gate_file = path_NS25+'stats/gate-vorticity/gate-'+filenames[i]+'/int17000-21000.nc'
    vort_gate_data = nc.Dataset(vort_gate_file,'r')
    vort_turb_area = vort_gate_data.variables['Partition2'][:,:]
    vort_turb_area_fracs_NS25[:,:,i] = vort_turb_area
    vort_turb_area_fracs_NS25_timmean = np.mean(vort_turb_area_fracs_NS25,axis=0)
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
    vort_turb_area = vort_gate_data.variables['Partition2'][:,:]
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


########################################################################
# Colourmaps

imola_data = np.loadtxt(colourmap_path+'imola/imola.txt')
imola_map = LinearSegmentedColormap.from_list('imola',imola_data)

########################################################################
# Plot

f, (ax1,ax2) = plt.subplots(1,2,sharex='all',sharey='all',figsize=(10,5))
ax1.tick_params(bottom=True,top=True,left=True,right=True)
ax2.tick_params(bottom=True,top=True,left=True,right=True)
ax1.set_ylim(1,1.5)
ax1.set_xlim(np.log10(vort_thresholds[0]/(ceps*B0/nu_117)),np.log10(vort_thresholds[-1]/(ceps*B0/nu_25)))
cs1 = ax1.contourf(np.log10(vort_thresholds**2/(ceps*B0/nu_117)),NS117.y/np.mean(NS117.z_enc),vort_turb_area_fracs_NS117_timmean,cmap=imola_map,levels=np.arange(0,1.1,0.1))
cs2 = ax2.contourf(np.log10(vort_thresholds**2/(ceps*B0/nu_25)),NS25.y/np.mean(NS25.z_enc),vort_turb_area_fracs_NS25_timmean,cmap=imola_map,levels=np.arange(0,1.1,0.1))
cbar_ax = f.add_axes([0.3,0.15,0.5,0.03])
cbar = f.colorbar(cs1,cax=cbar_ax,orientation='horizontal')
ax1.axhline(np.mean(NS117.z_ig/NS117.z_enc),0,0.05,color='k',linewidth=2)
ax1.axhline(np.mean(NS117.z_if/NS117.z_enc),0,0.05,color='k',linewidth=2)
ax2.axhline(np.mean(NS25.z_ig/NS25.z_enc),0,0.05,color='k',linewidth=2)
ax2.axhline(np.mean(NS25.z_if/NS25.z_enc),0,0.05,color='k',linewidth=2)
cbar.ax.set_xlabel(r'$a_\mathrm{T}$')
ax1.set_ylabel(r'$z/z_\mathrm{enc}$')
ax1.set_xlabel(r'$\log_{10}(\omega_\mathrm{th}^2/\omega_0^2)$')
ax2.set_xlabel(r'$\log_{10}(\omega_\mathrm{th}^2/\omega_0^2)$')
ax1.set_title('(a) $Re_0=117$',loc='left',fontsize=20)
ax2.set_title('(b) $Re_0=25$',loc='left',fontsize=20)
ax1.axvline(-0.62,0,NS117.y[-1]/np.mean(NS117.z_enc),color='k',linestyle='--')
ax2.axvline(-0.5,0,NS117.y[-1]/np.mean(NS117.z_enc),color='k',linestyle='--')
plt.tight_layout(rect=[0,0.15,1,1])
plt.savefig(opath+'Fig12.pdf',bbox_inches='tight')
plt.show()

f, (ax1,ax2) = plt.subplots(1,2,sharex='all',sharey='all',figsize=(10,5))
ax1.tick_params(bottom=True,top=True,left=True,right=True)
ax2.tick_params(bottom=True,top=True,left=True,right=True)
ax1.set_ylim(1,1.5)
cs1 = ax1.contourf(np.log10(vort_thresholds**2/(ceps*B0/nu_42)),NS42.y/np.mean(NS42.z_enc),vort_turb_area_fracs_NS42_timmean,cmap=imola_map,levels=np.arange(0,1.1,0.1))
cs2 = ax2.contourf(np.log10(vort_thresholds**2/(ceps*B0/nu_42)),S20.y/np.mean(S20.z_enc),vort_turb_area_fracs_S20_timmean,cmap=imola_map,levels=np.arange(0,1.1,0.1))
ax1.axhline(np.mean(NS42.z_ig)/np.mean(NS42.z_enc),0,0.05,color='k',linewidth=2)
ax1.axhline(np.mean(NS42.z_if)/np.mean(NS42.z_enc),0,0.05,color='k',linewidth=2)
ax2.axhline(np.mean(S20.z_ig)/np.mean(S20.z_enc),0,0.05,color='k',linewidth=2)
ax2.axhline(np.mean(S20.z_if)/np.mean(S20.z_enc),0,0.05,color='k',linewidth=2)
cbar_ax = f.add_axes([0.3,0.15,0.5,0.03])
cbar = f.colorbar(cs1,cax=cbar_ax,orientation='horizontal')
cbar.ax.set_xlabel(r'$a_\mathrm{T}$')
ax1.set_ylabel(r'$z/z_\mathrm{enc}$')
ax1.set_xlabel(r'$\log_{10}(\omega_\mathrm{th}^2/\omega_0^2)$')
ax2.set_xlabel(r'$\log_{10}(\omega_\mathrm{th}^2/\omega_0^2)$')
ax1.set_title('(a) $Fr_0=0$',loc='left',fontsize=20)
ax2.set_title('(b) $Fr_0=20$',loc='left',fontsize=20)
ax1.axvline(-0.46,0,NS42.y[-1]/np.mean(NS42.z_enc),color='k',linestyle='--')
ax2.axvline(-0.91,0,S20.y[-1]/np.mean(S20.z_enc),color='k',linestyle='--')
plt.tight_layout(rect=[0,0.15,1,1])
plt.savefig(opath+'Fig4.pdf',bbox_inches='tight')
plt.show()

f, ((ax1,ax2),(ax3,ax4),(ax5,ax6)) = plt.subplots(3,2,sharex='all',sharey='row',figsize=(10,15))
ax1.grid(True,linewidth=1.5)
ax2.grid(True,linewidth=1.5)
ax3.grid(True,linewidth=1.5)
ax4.grid(True,linewidth=1.5)
ax5.grid(True,linewidth=1.5)
ax6.grid(True,linewidth=1.5)
ax1.tick_params(bottom=False,left=False)
ax2.tick_params(bottom=False,left=False)
ax3.tick_params(bottom=False,left=False)
ax4.tick_params(bottom=False,left=False)
ax5.tick_params(bottom=False,left=False)
ax6.tick_params(bottom=False,left=False)
ax5.set_xticks([-2,-1,0,1])
ax5.set_xlim(-2.5,1.5)
ax1.set_ylim(0,0.3)
ax3.set_ylim(0,2)
ax5.set_ylim(0,0.6)
ax1.plot(np.log10(vort_thresholds**2/(ceps*B0/nu_42)),-vort_nonturb_s1_flux_zif_NS42_timmean/B0,c=blues(0.5),label=r'$Fr_0=0$')
ax1.plot(np.log10(vort_thresholds**2/(ceps*B0/nu_42)),-vort_nonturb_s1_flux_zif_S20_timmean/B0,c=blues(0.9),label=r'$Fr_0=20$')
ax2.plot(np.log10(vort_thresholds**2/(ceps*B0/nu_42)),-vort_turb_s1_flux_zif_NS42_timmean/B0,c=blues(0.5),label=r'$Fr_0=0$')
ax2.plot(np.log10(vort_thresholds**2/(ceps*B0/nu_42)),-vort_turb_s1_flux_zif_S20_timmean/B0,c=blues(0.9),label=r'$Fr_0=20$')
ax3.plot(np.log10(vort_thresholds**2/(ceps*B0/nu_42)),np.sqrt(vort_nonturb_s1_var_zif_NS42_timmean)/(N**2*L0),c=blues(0.5),label=r'$Fr_0=0$')
ax3.plot(np.log10(vort_thresholds**2/(ceps*B0/nu_42)),np.sqrt(vort_nonturb_s1_var_zif_S20_timmean)/(N**2*L0),c=blues(0.9),label=r'$Fr_0=20$')
ax4.plot(np.log10(vort_thresholds**2/(ceps*B0/nu_42)),np.sqrt(vort_turb_s1_var_zif_NS42_timmean)/(N**2*L0),c=blues(0.5),label=r'$Fr_0=0$')
ax4.plot(np.log10(vort_thresholds**2/(ceps*B0/nu_42)),np.sqrt(vort_turb_s1_var_zif_S20_timmean)/(N**2*L0),c=blues(0.9),label=r'$Fr_0=20$')
ax5.plot(np.log10(vort_thresholds**2/(ceps*B0/nu_42)),np.sqrt(vort_nonturb_w_var_zif_NS42_timmean)/(np.mean(NS42.z_enc)*B0)**(1./3.),c=blues(0.5),label=r'$Fr_0=0$')
ax5.plot(np.log10(vort_thresholds**2/(ceps*B0/nu_42)),np.sqrt(vort_nonturb_w_var_zif_S20_timmean)/(np.mean(S20.z_enc)*B0)**(1./3.),c=blues(0.9),label=r'$Fr_0=20$')
ax6.plot(np.log10(vort_thresholds**2/(ceps*B0/nu_42)),np.sqrt(vort_turb_w_var_zif_NS42_timmean)/(np.mean(NS42.z_enc)*B0)**(1./3.),c=blues(0.5),label=r'$Fr_0=0$')
ax6.plot(np.log10(vort_thresholds**2/(ceps*B0/nu_42)),np.sqrt(vort_turb_w_var_zif_S20_timmean)/(np.mean(S20.z_enc)*B0)**(1./3.),c=blues(0.9),label=r'$Fr_0=20$')
ax1.axvline(-0.46,0,NS42.y[-1]/np.mean(NS42.z_enc),color=blues(0.5),linestyle='--')
ax2.axvline(-0.46,0,NS42.y[-1]/np.mean(NS42.z_enc),color=blues(0.5),linestyle='--')
ax1.axvline(-0.91,0,S20.y[-1]/np.mean(S20.z_enc),color=blues(0.9),linestyle='--')
ax2.axvline(-0.91,0,S20.y[-1]/np.mean(S20.z_enc),color=blues(0.9),linestyle='--')
ax3.axvline(-0.46,0,NS42.y[-1]/np.mean(NS42.z_enc),color=blues(0.5),linestyle='--')
ax4.axvline(-0.46,0,NS42.y[-1]/np.mean(NS42.z_enc),color=blues(0.5),linestyle='--')
ax3.axvline(-0.91,0,S20.y[-1]/np.mean(S20.z_enc),color=blues(0.9),linestyle='--')
ax4.axvline(-0.91,0,S20.y[-1]/np.mean(S20.z_enc),color=blues(0.9),linestyle='--')
ax5.axvline(-0.46,0,NS42.y[-1]/np.mean(NS42.z_enc),color=blues(0.5),linestyle='--')
ax6.axvline(-0.46,0,NS42.y[-1]/np.mean(NS42.z_enc),color=blues(0.5),linestyle='--')
ax5.axvline(-0.91,0,S20.y[-1]/np.mean(S20.z_enc),color=blues(0.9),linestyle='--')
ax6.axvline(-0.91,0,S20.y[-1]/np.mean(S20.z_enc),color=blues(0.9),linestyle='--')
ax5.set_xlabel(r'$\log_{10}(\omega_\mathrm{th}^2/\omega_0^2)$')
ax6.set_xlabel(r'$\log_{10}(\omega_\mathrm{th}^2/\omega_0^2)$')
ax1.set_ylabel(r'$-(\langle b^\prime w^\prime \rangle_\mathrm{NT})_{z_{i,f}}/B_0$')
ax2.set_ylabel(r'$-(\langle b^\prime w^\prime \rangle_\mathrm{T})_{z_{i,f}}/B_0$')
ax3.set_ylabel(r'$((b_\mathrm{rms})_\mathrm{NT})_{z_{i,f}}/(N^2L_0)$')
ax4.set_ylabel(r'$((b_\mathrm{rms})_\mathrm{T})_{z_{i,f}}/(N^2L_0)$')
ax5.set_ylabel(r'$((w_\mathrm{rms})_\mathrm{NT})_{z_{i,f}}/w_*$')
ax6.set_ylabel(r'$((w_\mathrm{rms})_\mathrm{T})_{z_{i,f}}/w_*$')
ax1.set_title('(a)',loc='left',fontsize=20)
ax2.set_title('(b)',loc='left',fontsize=20)
ax3.set_title('(c)',loc='left',fontsize=20)
ax4.set_title('(d)',loc='left',fontsize=20)
ax5.set_title('(e)',loc='left',fontsize=20)
ax6.set_title('(f)',loc='left',fontsize=20)
ax1.legend(loc='best',fontsize=20,handlelength=1)
plt.tight_layout()
plt.savefig(opath+'Fig10.pdf',bbox_inches='tight')
plt.show()


