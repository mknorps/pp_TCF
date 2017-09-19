# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# File name: time.py
# Created by: gemusia
# Creation date: 19-09-2017
# Last modified: 19-09-2017 15:02:04
# Purpose: module containing costants, helper functions
#          and declaration of atributes
#          used in all types of time scales plots
# 
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


import numpy as np
import matplotlib.pyplot as plt
import homfigs as hfig
from os.path import expanduser,isfile
from itertools import product
import csv


#--------------------------------------------------------------------#
# constants
#--------------------------------------------------------------------#

#time steps in LES and DNS
dt_DNS = 0.005
dt_LES = 0.02

#normalisation factors
utau = 0.0429
Retau = 150
ttau = 1.0/(utau*Retau)

#flow parameters
n = 32
m = 64
o = 32
length = 2.0*np.pi
width =4.0*np.pi 

#--------------------------------------------------------------------#
# helper functions
#--------------------------------------------------------------------#
def symm(lst):
    ll = len(lst)
    return 0.5*(lst[:(ll+1)/2] + np.flipud(lst)[:(ll+1)/2])

def Chebyshev_zeroes(N):
    return 0.5*(np.array(range(1,N+1))+np.array(range(N)))/float(N) * 150

def non_dimensialisation(lst):
    return Retau - lst*Retau

def y(j):
    return np.cos(float(j)*np.pi/float(n))

def filter_width(j):
    return (0.5*length*width*(y(j-1)-y(j+1))/float(m*o))**(1/3.0)

#--------------------------------------------------------------------#
# declaration of picture attributes
#--------------------------------------------------------------------#

dt = {"LES":dt_LES,"SGSles":dt_DNS,"SGSdns":dt_DNS,"DNS":dt_DNS, "apriori":dt_DNS}

#particle type
ptype = ["DNS","LES","apriori","SGSles","SGSdns"]

#list of Stokes numbers (size) of particles
St = ["fluid","St0.5","St1","St5","St25","St125"]


start_iteration = {'DNS':np.arange(2475,2525,5),  #start iteration list
                   'apriori':np.arange(2475,2525,5),
                   'LES':np.arange(6805,6855,5),
                   'SGSles':np.arange(2475,2525,5),
                   'SGSdns':np.arange(2475,2525,5)
                  }

datalines = {0:{"ls": "solid","c": "green","lbl":"y=1"},
             2:{"ls": "dashed","c": "blue","lbl":"y=3"},
             5:{"ls": "dotted","c": "black","lbl":"y=6"},
             7:{"ls": "dashdot","c": "red","lbl":"y=center"}}





file_path = expanduser("~") + "/wyniki/time_scales/"
pict_path = file_path


if __name__=='__main__':

    #--------------------------------------------------------------------#
    # Data reading, symmetrisation and averaging
    #
    # must be done once and results stored in .csv file
    #--------------------------------------------------------------------#

    #the data is written to file in Fortran program with following code
    '''
    ns   - number of slices the channel is divided

    for correlation coefficients
          write(5,'(i8,16e13.5)') t-psteps,  (Rij(j,f3),j=1,ns)

    for time scales
       do j=1,ns
           write(5,'(2i6,e13.5)') j,num_t0(j,i), Tij_cut(j,i)
       enddo

    '''

    rho = {}
    tau = {}

    for simulation, St in product (ptype, Stlist):


        print "working on ", simulation, St
        #++++++++++++++++++++++++++++++++++++++++++++++
        # reading multiple files, symmetrising and
        # averaging the data in
        # respective columns

        f_name_base = file_path + "rho_" + St + "_"  + simulation + "_"
        f_name_start = f_name_base + str(start_iteration[simulation][0]) +".dat"

        with open(f_name_start) as f:
            f_iter = [symm(np.transpose(np.loadtxt(f,skiprows=1))[1:])]
            f_iter_shape = f_iter[0].shape



        for itno in start_iteration[simulation][1:]:

            f_name_itno =  f_name_base +str(itno) +".dat"

            # first row of read data is iteration number, so it is ommited
            with open(f_name_itno) as f:
                f_new = [symm(np.transpose(np.loadtxt(f,skiprows=1))[1:])]

                if f_new[0].shape <> f_iter_shape:
                    raise AssertionError("shape of the first file dos not match shape of",  f_name_base)

                f_iter = np.append(f_iter, f_new , axis=0)



        rho[St + "_" + simulation] =  np.mean(f_iter,axis=0)

        #++++++++++++++++++++++++++++++++++++++++++++++
        # computing relaxation time tau = min {t; rho(t)<1\e}

        # print next(i for i,v in enumerate(rho[St + "_" + simulation][7]) if (v<np.exp(-1) or i==len(rho[St + "_" + simulation][7])-1))

        tau[St + "_" + simulation] = np.array(map(lambda x: ptype_dt[simulation]*
                next(i for i,v in enumerate(x) if (v<np.exp(-1) or i==len(x)-1)),rho[St + "_" + simulation]))/ttau

        n_datapoints = len(tau[St + "_" + simulation])

    x_points = Chebyshev_zeroes(n_datapoints)

    # write tau to file     
    f_name_tau = file_path + "tau.csv"

    with open(f_name_tau,'wb') as f:
        w = csv.writer(f)
        w.writerow(tau.keys())
        w.writerows(zip(*tau.values()))

        print f_name_tau, " created"



