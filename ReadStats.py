# Created on 16.02.2016
# @author: Armin Haghshenas
# Modified by Katherine Fodor

import netCDF4 as nc 	# reading NetCDF4 files.
import os
import struct
import numpy as np

class Statistics:
    # define one instance per statistic
    def __init__(self,data_path):
        N = np.sqrt(3)                  # stratification in free atmosphere
        
        self.data_path = data_path		# path to data 
        
        # avg statistics
        path_avg = data_path
        datafile_avg = nc.Dataset(path_avg, 'r')
        
        # Reading the variables
        self.t     = datafile_avg.variables['t'][:]		  # time				
        self.y     = datafile_avg.variables['y'][:]		  # height	
        # for 2D vars, the first index is time, the second is height
        self.rU    = datafile_avg.variables['rU'][:,:] 				
        self.rV    = datafile_avg.variables['rV'][:,:] 	  # mean velocity
        self.rW    = datafile_avg.variables['rW'][:,:] 				

        self.Tke   = datafile_avg.variables['Tke'][:,:]   # turbulent kinetic energy
        self.Rxx   = datafile_avg.variables['Rxx'][:,:]   # u variance
        self.Ryy   = datafile_avg.variables['Ryy'][:,:]   # v variance (vertical velocity)
        self.Rzz   = datafile_avg.variables['Rzz'][:,:]   # w variance
        
        self.Wx    = datafile_avg.variables['Wx'][:,:] 				
        self.Wy    = datafile_avg.variables['Wy'][:,:] 	  # mean vorticity
        self.Wz    = datafile_avg.variables['Wz'][:,:] 				
        
        self.Wx2   = datafile_avg.variables['Wx2'][:,:] 				
        self.Wy2   = datafile_avg.variables['Wy2'][:,:]   # vorticity variance
        self.Wz2   = datafile_avg.variables['Wz2'][:,:]  
        
        self.Prd   = datafile_avg.variables['Prd'  ][:,:] # shear production rate
        self.Buo   = datafile_avg.variables['Buo'  ][:,:] # buoyancy production/destruction rate 
        self.Eps   = datafile_avg.variables['Eps'  ][:,:] # viscous dissipation rate
        self.Tke_t = datafile_avg.variables['Tke_t'][:,:] # time rate of change of TKE
        
        self.Eta   = datafile_avg.variables['Eta'  ][:,:] # Kolmogorov scale
        
        # avg1s statistics 
        path_avg1s = path_avg.replace('avg','avg1s')	  # assume that name is the same except for avg -> avg1s	
        datafile_avg1s = nc.Dataset(path_avg1s, 'r')

        self.rS   = datafile_avg1s.variables['rS'  ][:,:] # mean buoyancy
        self.rS_y = datafile_avg1s.variables['rS_y'][:,:] # mean buoyancy gradient in vertical direction
        self.Rsv  = datafile_avg1s.variables['Rsv' ][:,:] # vertical buoyancy flux
        self.r2S  = datafile_avg1s.variables['rS2' ][:,:] # buoyancy variance
        # avg2s statistics 
        path_avg2s = path_avg.replace('avg','avg2s')      # assume that name is the same except for avg -> avg2s
        exists = os.path.isfile(path_avg2s)               # check that avg2s stats exist
        if exists:
            datafile_avg2s = nc.Dataset(path_avg2s, 'r')
            
            self.rS2   = datafile_avg2s.variables['rS'  ][:,:]	# mean scalar 2
            self.rS2_y = datafile_avg2s.variables['rS_y'][:,:]	# mean scalar 2 gradient in vertical direction
            self.Rs2v  = datafile_avg2s.variables['Rsv' ][:,:]	# vertical scalar 2 flux
            self.r2S2  = datafile_avg2s.variables['rS2' ][:,:]  # scalar 2 variance

        # PV statistics
        path_pv = path_avg.replace('avg','PV')             # assume that name is the same except for avg -> PV
        
        exists = os.path.isfile(path_pv)                   # check that PV stats exist
        if exists:
            datafile_pv = nc.Dataset(path_pv, 'r')
            self.PVMom1 = datafile_pv.variables['PVMom1'][:,:] # mean potential vorticity
            self.PVMom2 = datafile_pv.variables['PVMom2'][:,:] # potential vorticity variance
        
        # time and height info 
        self.t_len = len(self.t)
        self.t_ind = self.t_len-1 				# max time index	
        
        self.y_len = len(self.y)
        self.y_ind = self.y_len-1				# max y index

        # scales
        self.L_Oz  = np.zeros((self.t_len,self.y_len))     # Ozmidov length
        for n in range(0,self.t_len):
            self.L_Oz[n,:]    = (self.Eps[n,:]/self.rS_y[n,:]**(3./2.))**0.5

        # heights

        self.z_enc     = np.array([np.sqrt((2/N**2)*np.trapz(self.rS[n,:]-N**2*self.y,self.y)) for n in range(0,self.t_len)]) # encroachment height (mixed-layer depth)
        self.z_enc_arg = np.array([int(np.argmin(np.abs(self.z_enc[n]-self.y))) for n in range(0,self.t_len)])
        self.z_ig      = np.array([self.y[np.argmax(self.rS_y[n,:])] for n in range(0,self.t_len)])                           # height of maximum mean buoyancy gradient
        self.z_ig_arg  = np.array([int(np.argmin(np.abs(self.z_ig[n]-self.y))) for n in range(0,self.t_len)])
        self.z_if      = np.array([self.y[np.argmin(self.Rsv[n,:])] for n in range(0,self.t_len)])                            # height of minimum buoyancy flux
        self.z_if_arg  = np.array([int(np.argmin(np.abs(self.z_if[n]-self.y))) for n in range(0,self.t_len)])
        self.L_Oz_zif  = np.array([self.L_Oz[n,self.z_if_arg[n]] for n in range(0,self.t_len)])                               # Ozmidov length evaluated at height of min. buoyancy flux
        self.z_is      = np.array([self.z_ig[n] - 1.78*self.L_Oz_zif[n] for n in range(0,self.t_len)])                        # height of transition from lower to upper entrainment zone sublayer
        self.z_is_arg  = np.array([int(np.argmin(np.abs(self.z_is[n]-self.y))) for n in range(0,self.t_len)])

        self.b_enc = N**2*self.z_enc                                                                                          # buoyancy scale in mixed-layer


        
class Conditional_Stats:
    
    def __init__(self,data_path_int,data_path_partition1=None,data_path_partition2=None):
        self.data_path = data_path_int
        
        # intermittency factor
        path_int = data_path_int
        datafile_int = nc.Dataset(path_int,'r')
        self.int1 = datafile_int.variables['Partition1'][:,:]                      # intermittency factor in partition 1
        self.int2 = datafile_int.variables['Partition2'][:,:]                      # intermittency factor in partition 2
        
        # cavg statistics
        if data_path_partition1 is not None:
            path_partition1 = data_path_partition1
            datafile_partition1 = nc.Dataset(path_partition1,'r')
            self.P1UMom1 = datafile_partition1.variables['UMom1'][:,:]             # mean u velocity 
            self.P1UMom2 = datafile_partition1.variables['UMom2'][:,:]             # u velocity variance 
            self.P1VMom1 = datafile_partition1.variables['VMom1'][:,:]             # mean v (vertical) velocity
            self.P1VMom2 = datafile_partition1.variables['VMom2'][:,:]             # v (vertical) velocity variance
            self.P1WMom1 = datafile_partition1.variables['WMom1'][:,:]             # mean w velocity
            self.P1WMom2 = datafile_partition1.variables['WMom2'][:,:]             # w velocity variance
            self.P1S1Mom1 = datafile_partition1.variables['Scalar1Mom1'][:,:]      # mean buoyancy
            self.P1S1Mom2 = datafile_partition1.variables['Scalar1Mom2'][:,:]      # buoyancy variance
            if 'Scalar2Mom1' in datafile_partition1.variables:
                self.P1S2Mom1 = datafile_partition1.variables['Scalar2Mom1'][:,:]  # mean scalar 2
                self.P1S2Mom2 = datafile_partition1.variables['Scalar2Mom2'][:,:]  # scalar 2 variance

        if data_path_partition2 is not None:
            path_partition2 = data_path_partition2
            datafile_partition2 = nc.Dataset(path_partition2,'r')
            self.P2UMom1 = datafile_partition2.variables['UMom1'][:,:]             # as above but for partition 2...
            self.P2UMom2 = datafile_partition2.variables['UMom2'][:,:]
            self.P2VMom1 = datafile_partition2.variables['VMom1'][:,:]
            self.P2VMom2 = datafile_partition2.variables['VMom2'][:,:]
            self.P2WMom1 = datafile_partition2.variables['WMom1'][:,:]
            self.P2WMom2 = datafile_partition2.variables['WMom2'][:,:]
            self.P2S1Mom1 = datafile_partition2.variables['Scalar1Mom1'][:,:]
            self.P2S1Mom2 = datafile_partition2.variables['Scalar1Mom2'][:,:]
            if 'Scalar2Mom1' in datafile_partition2.variables:
                self.P2S2Mom1 = datafile_partition2.variables['Scalar2Mom1'][:,:]
                self.P2S2Mom2 = datafile_partition2.variables['Scalar2Mom2'][:,:]

        # avgGi statistics
        if data_path_partition1 is not None:
            path_avgGi_p1 = path_partition1.replace('cavg','avgGi') # assume that name is the same except for cavg -> avgGi
            exists = os.path.isfile(path_avgGi_p1)
            if exists:
                datafile_avgGi_p1 = nc.Dataset(path_avgGi_p1,'r')
                self.P1S1GradY = datafile_avgGi_p1.variables['GradientYMom1'][:,:] # mean buoyancy gradient in vertical direction

        if data_path_partition2 is not None:
            path_avgGi_p2 = path_partition2.replace('cavg','avgGi')
            exists = os.path.isfile(path_avgGi_p2)
            if exists:
                datafile_avgGi_p2 = nc.Dataset(path_avgGi_p2,'r')
                self.P2S1GradY = datafile_avgGi_p2.variables['GradientYMom1'][:,:] # as above but for partition 2...

        # avgEps statistics
        if data_path_partition1 is not None:
            path_avgEps_p1 = path_partition1.replace('cavg','avgEps') # assume that name is the same except for cavg -> avgEps
            exists = os.path.isfile(path_avgEps_p1)
            if exists:
                datafile_avgEps_p1 = nc.Dataset(path_avgEps_p1,'r')
                self.P1Eps = datafile_avgEps_p1.variables['EpsMom1'][:,:] # viscous dissipation rate

        if data_path_partition2 is not None:
            path_avgEps_p2 = path_partition2.replace('cavg','avgEps') 
            exists = os.path.isfile(path_avgEps_p2)
            if exists:
                datafile_avgEps_p2 = nc.Dataset(path_avgEps_p2,'r')
                self.P2Eps = datafile_avgEps_p2.variables['EpsMom1'][:,:] # as above but for partition 2...

        # avgMom statistics
        if data_path_partition1 is not None:
            path_avgMom_p1 = path_partition1.replace('cavg','avgMom') # assume that name is the same except for cavg -> avgMom
            exists = os.path.isfile(path_avgMom_p1)
            if exists:
                datafile_avgMom_p1 = nc.Dataset(path_avgMom_p1,'r')
                self.P1tauy1 = datafile_avgMom_p1.variables['tauy1Mom1'][:,:] # molecular vertical buoyancy flux
                self.P1v1 = datafile_avgMom_p1.variables['v1Mom1'][:,:] # turbulent vertical buoyancy flux

        if data_path_partition2 is not None:
            path_avgMom_p2 = path_partition2.replace('cavg','avgMom')
            exists = os.path.isfile(path_avgMom_p2)
            if exists:
                datafile_avgMom_p2 = nc.Dataset(path_avgMom_p2,'r')
                self.P2tauy1 = datafile_avgMom_p2.variables['tauy1Mom1'][:,:]
                self.P2v1 = datafile_avgMom_p2.variables['v1Mom1'][:,:] # same as above but for partiton 2...
        
    
        
class Pdfs:
    # sizeofdata in bytes: 1 for gate files, 4 otherwise
    # etype: < = little-endian, > = big-endian
    # dtype: f = floating point number, B = unsigned character (for gate files)
    # dimension: 1 for univariate, 2 for bivariate...so far only univariate
    def __init__(self,data_path_pdf_list,data_path_ydat,sizeofdata=4,etype='<',dtype='f',dimension=1):
        
        self.y = np.loadtxt(data_path_ydat)
        self.ny = len(self.y)
        
        # obtain size from first file
        # the last level contains the global Pdf
        # the last two entries per level are min/max values
        self.nb = int(os.stat(data_path_pdf_list[0]).st_size / sizeofdata / (self.ny + 1)) - 2        
        print("Files with {} bins and {} levels.".format(self.nb,self.ny))
      
        # processing data
        pdf = np.zeros((len(data_path_pdf_list),(self.ny + 1) * (self.nb + 2)),dtype=float)
        for n in range(len(data_path_pdf_list)):
            print("Processing file {} ...".format(data_path_pdf_list[n]))
            fin = open(data_path_pdf_list[n], 'rb')
            raw = fin.read()
            pdf[n,:] = np.array(struct.unpack((etype+'{}'+dtype).format(int(fin.tell()/sizeofdata)), raw))
            fin.close()
        pdf = np.reshape(pdf,(len(data_path_pdf_list),self.ny + 1,self.nb + 2))
        self.pdf = pdf

        # interpolating data onto same grid
        for n in range(len(data_path_pdf_list)):
            for j in range(self.ny):
                self.pdf[n,j,:-2] = np.interp(np.linspace(self.pdf[0,j,self.nb],self.pdf[0,j,self.nb + 1],num=self.nb),np.linspace(self.pdf[n,j,self.nb],self.pdf[n,j,self.nb + 1],num=self.nb),self.pdf[n,j,:-2])
        self.pdf_timeavg = np.mean(self.pdf,axis=0)   

        
        # normalizing histograms to obtain pdf (s.t. it integrates to 1 using midpoint rule)     
        for n in range(len(data_path_pdf_list)):
            for j in range(self.ny + 1):
                samplesize = np.sum(self.pdf[n,j,:self.nb])
                samplestep = (self.pdf[n,j,self.nb + 1] - self.pdf[n,j,self.nb]) / (self.nb - 1)
                if samplestep > 0: # otherwise the pdf is zero, by construction in Tlab
                    self.pdf[n,j,:self.nb] = self.pdf[n,j,:self.nb] / (samplesize * samplestep)
        for j in range(self.ny + 1):
            samplesize = np.sum(self.pdf_timeavg[j,:self.nb])
            samplestep = (self.pdf_timeavg[j,self.nb + 1] - self.pdf_timeavg[j,self.nb]) / (self.nb - 1)
            if samplestep > 0: 
                self.pdf_timeavg[j,:self.nb] = self.pdf_timeavg[j,:self.nb] / (samplesize * samplestep)
                
        
        # axis information
        self.xy = np.zeros((2,self.ny,self.nb),dtype=float)
        for j in range(self.ny):
            self.xy[0,j,:] = np.linspace(self.pdf[0,j,self.nb],self.pdf[0,j,self.nb + 1],num=self.nb)
            self.xy[1,j,:] = self.y[j]

        
