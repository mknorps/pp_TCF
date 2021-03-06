# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# File name: apriori_SGS_fluid.py
# Created by: mknorps 
# Creation date: 21-06-2017
# Last modified: 29-09-2017 12:35:17
# Purpose: take filtered and unfiltered fluid field 
#          in Fourier space from spectral code 
#          compute statistics of SGS fluid velocity
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import numpy as np
import matplotlib.pyplot as plt
import homstat as hs
import homfiles as hf
import homfigs as hfig
import itertools as it
from scipy.interpolate import griddata
from os.path import expanduser
import sys
import datetime

# declaration of picture attributes

LineStyle = {'DNS':'solid','LES':'dashed'}
coordinates = {1:'x',2:'y',3:'z'}
statistics = {0:("hmean_symm",),1:("hstd_symm",),
        2:("hcor_symm","Ux","Uy"),
        3:("hcor_symm","Uy","Uz"),
        4:("hcor_symm","Ux","Uz"),
        5:("hek_symm",)}
ylabels = {1:"$<U>_{x}$",2:"$<U>_{y}$",
        3:"$<U>_{z}$",4:"$<U_{x},U_{y}>$",
        5:"$<U_{y},U_z>$",6:"$<U_{x}U_{z}$",
        7:"$e_{k}$"}
ptitles = {"hmean_symm":"mean of U","hstd_symm":"Standard deviation of U","hcor_symm":"Correlation of ","hek_symm":"Kinetic energy"}


#data loading
file_path    = expanduser("~") + "/wyniki/apriori/fluid/"
fig_path     = file_path
field        = 2501

test_set_DNS = hf.ChannelFields(128,129,128,field,field,fPath = file_path, fPrefix=["upp_Ux","upp_Uy","upp_Uz"])
test_set_LES = hf.ChannelFields(32,33,64,field,field, fPath = file_path,fPrefix=["uppf_Ux","uppf_Uy","uppf_Uz"])


# The data that will be used for computations
#
DNS_data = getattr(test_set_DNS,'field_'+str(field)) 
LES_data = getattr(test_set_LES,'field_'+str(field)) 
#


print (DNS_data.K,DNS_data.N, DNS_data.M )
print (LES_data.K,LES_data.N, LES_data.M )


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# computation of SGS velocity usgs = DNS_data - LES_data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#.....................................................
#  VERSION 1: we take the difference in LES nodes only
#
#DNS nodes connected with DNS
indices_32 =  np.arange(0,127,4)
indices_33 =  np.append(np.arange(0,127,4),127)
indices_64 =  np.arange(0,128,2)

Ux_SGS = np.empty([32,33,64])
Uy_SGS = np.empty([32,33,64])
Uz_SGS = np.empty([32,33,64])


for iLES,iDNS in enumerate(indices_32):
    for jLES,jDNS in enumerate(indices_33):
        for kLES,kDNS in enumerate(indices_64):
            Ux_SGS[iLES,jLES,kLES] = DNS_data.Ux[iDNS,jDNS,kDNS] - LES_data.Ux[iLES,jLES,kLES]
            Uy_SGS[iLES,jLES,kLES] = DNS_data.Uy[iDNS,jDNS,kDNS] - LES_data.Uy[iLES,jLES,kLES]
            Uz_SGS[iLES,jLES,kLES] = DNS_data.Uz[iDNS,jDNS,kDNS] - LES_data.Uz[iLES,jLES,kLES]


SGS_LES = hs.Channel(Ux_SGS,Uy_SGS,Uz_SGS,LES_data.K,LES_data.N, LES_data.M )


#.....................................................
#  VERSION 2: we take the difference in DNS nodes
#             and interpolate LES to DNS position
#             
#             For this we need to write the data to file
#             since interpolation takes a lot of time (~ 3h for each velocity = ~9h)
#             is run only if option 'INTERPOLATE' is activated
#             in the command line run
#             e.g.  python SGS_fluid_stats.py --interpolate

INTERPOLATE = False

if '--interpolate' in sys.argv:
    INTERPOLATE = True
    print ( 'Interpolation run... it may take a while ')
    print ( 'Argument list: ',str(sys.argv))


if INTERPOLATE:
    print("interpolating!!!!")
    Ux_SGS_dense_grid = np.empty([128,129,128])
    Uy_SGS_dense_grid = np.empty([128,129,128])
    Uz_SGS_dense_grid = np.empty([128,129,128])

    #DNS grid - we will interpolate LES velocity on it
    x = np.linspace(0,4*np.pi,128)
    y = DNS_data.ynodes() 
    z = np.linspace(0,2*np.pi,128)

    grid_x,grid_y,grid_z = np.meshgrid(x,y,z)

    #LES grid
    x_LES = np.linspace(0,4*np.pi,32)
    y_LES = LES_data.ynodes() 
    z_LES = np.linspace(0,2*np.pi,64)

    print(   len(x_LES),len(y_LES),len(z_LES))

    LES_grid_x,LES_grid_y,LES_grid_z = np.meshgrid(x_LES,y_LES,z_LES)
     
    points    = zip(LES_grid_x.flatten(), LES_grid_y.flatten(),LES_grid_z.flatten())
    value_Ux  = LES_data.Ux.flatten()  
    value_Uy  = LES_data.Uy.flatten()  
    value_Uz  = LES_data.Uz.flatten()  

    print("points: ", len(points), "  \t value_Ux : ",len(value_Ux))

    #Ux
    LES_interpolated_Ux = griddata(points,value_Ux,(grid_x,grid_y,grid_z),method ='linear' )

    with open('LES_interpolated_Ux.txt','w') as f:
        f.write('# Array shape: {0}\n'.format(LES_interpolated_Ux.shape))
        for LES_slice in LES_interpolated_Ux:
            np.savetxt(f,LES_slice)
    print ('Ux interpolated to DNS grid', datetime.datetime.strftime(datetime.datetime.now(),"%Y-%m-%d %H:%M:%S"))

    #Uy
    LES_interpolated_Uy = griddata(points,value_Uy,(grid_x,grid_y,grid_z),method ='linear' )
    with open('LES_interpolated_Uy.txt','w') as f:
        f.write('# Array shape: {0}\n'.format(LES_interpolated_Uy.shape))
        for LES_slice in LES_interpolated_Uy:
            np.savetxt(f,LES_slice)
    print ('Uy interpolated to DNS grid', datetime.datetime.strftime(datetime.datetime.now(),"%Y-%m-%d %H:%M:%S"))

    #Uz
    LES_interpolated_Uz = griddata(points,value_Uz,(grid_x,grid_y,grid_z),method ='linear' )
    with open('LES_interpolated_Uz.txt','w') as f:
        f.write('# Array shape: {0}\n'.format(LES_interpolated_Uz.shape))
        for LES_slice in LES_interpolated_Uz:
            np.savetxt(f,LES_slice)
    print ('Uz interpolated to DNS grid', datetime.datetime.strftime(datetime.datetime.now(),"%Y-%m-%d %H:%M:%S"))



Ux_data = np.loadtxt('LES_interpolated_Ux.txt').reshape(128,129,128)
Uy_data = np.loadtxt('LES_interpolated_Uy.txt').reshape(128,129,128)
Uz_data = np.loadtxt('LES_interpolated_Uz.txt').reshape(128,129,128)

Ux_SGS_dense = DNS_data.Ux - Ux_data
Uy_SGS_dense = DNS_data.Uy - Uy_data
Uz_SGS_dense = DNS_data.Uz - Uz_data

SGS_DNS = hs.Channel(Ux_SGS_dense,Uy_SGS_dense,Uz_SGS_dense,DNS_data.K,DNS_data.N, DNS_data.M )



