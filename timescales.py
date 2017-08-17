# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# File name: timescales.py
# Created by: gemusia
# Creation date: 26-07-2017
# Last modified: 17-08-2017 15:32:58
# Purpose: draw plots of time scales of turbulent channel flow
#          DNS, a priori LES, LES and SGS
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


import numpy as np
import matplotlib.pyplot as plt
import homfigs as hfig
from os.path import expanduser
from itertools import product 


#--------------------------------------------------------------------#
# constants 
#--------------------------------------------------------------------#

#time steps in LES and DNS
dt_DNS = 0.005
dt_LES = 0.02


#--------------------------------------------------------------------#
# helper functions
#--------------------------------------------------------------------#
def symm(lst):
    ll = len(lst)
    return 0.5*(lst[:(ll+1)/2] + np.flipud(lst)[:(ll+1)/2])



#--------------------------------------------------------------------#
# declaration of picture attributes
#--------------------------------------------------------------------#

LineStyle = {'fterm':'solid','pterm':'dashed',
        "Ux":'solid',"Ufx":'dashed',"Vx":'dotted',
        "UxUfx":'solid',"UyUfy":'dashed',"UzUfz":'dotted',
        "fluid":'solid',"St1":'dashed',"St5":'dotted',"St25":'dashdot'}

ptype_dt = {"LES":dt_LES,"SGSles":dt_DNS,"SGSdns":dt_DNS,"DNS":dt_DNS, "apriori":dt_DNS}

#ptype = ["DNS","LES","apriori","SGSles","SGSdns"]
#Stlist = ["fluid","St0.1","St1","St5","St25"]
ptype = ["LES","DNS","apriori","SGSles","SGSdns"]

Stlist = ["fluid","St0.5","St1","St25","St125"]
#Stlist = ["fluid","St0.5","St1","St5","St25","St125"]

start_iteration = {'DNS':np.arange(2475,2525,5),  #start iteration list
                   'apriori':np.arange(2475,2525,5), 
                   'LES':np.arange(6805,6855,5), 
                   'SGSles':np.arange(2475,2525,5), 
                   'SGSdns':np.arange(2475,2525,5) 
                  }
'''
ptype = ["SGSles","SGSdns"]

Stlist = ["fluid","St1"]

start_iteration = [2555,2560,2565,2570,2573,2577,2580]  #start iteration list
'''
datalines = {0:{"ls": "solid","c": "green","lbl":"y=1"},
             2:{"ls": "dashed","c": "blue","lbl":"y=3"},
             5:{"ls": "dotted","c": "black","lbl":"y=6"},
             7:{"ls": "dashdot","c": "red","lbl":"y=center"}}



#--------------------------------------------------------------------#
# loading the data
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

file_path = expanduser("~") + "/wyniki/time_scales/"
pict_path = file_path


#--------------------------------------------------------------------#
# Data reading, symmetrisation and averaging
# 
#--------------------------------------------------------------------#

rho = {}
tau = {}

for simulation, St in product (ptype, Stlist):


    print "working on ", simulation, St
    #++++++++++++++++++++++++++++++++++++++++++++++
    # reading multiple files, symmetrising and 
    # averaging the data in
    # respective columns


    f_iter = [symm(np.transpose(np.loadtxt(file_path + "rho_" + St + "_"  + simulation + "_" 
                                  +str(start_iteration[ptype][0]) +".dat",skiprows=1))[1:])]

    for itno in start_iteration[ptype][1:]:
       
        # first row of read data is iteration number, so it is ommited
        f_iter = np.append(f_iter,[symm(np.transpose(np.loadtxt(file_path + "rho_" + St + "_" 
            + simulation + "_" +str(itno) +".dat",skiprows=1))[1:])] , axis=0)

    rho[St + "_" + simulation] =  np.mean(f_iter,axis=0)

    #++++++++++++++++++++++++++++++++++++++++++++++
    # computing relaxation time tau = min {t; rho(t)<1\e}

    print next(i for i,v in enumerate(rho[St + "_" + simulation][7]) if v<np.exp(-1))

    tau[St + "_" + simulation] = map(lambda x: ptype_dt[simulation]*
            next(i for i,v in enumerate(x) if (v<np.exp(-1) or i==len(x)-1)),rho[St + "_" + simulation]) 



#--------------------------------------------------------------------#
# PLOT DRAW
#
# we use symmetrised data
#--------------------------------------------------------------------#


for simulation, St in product (ptype, Stlist):

    corrfig = hfig.Homfig(title= St + ": Correlation coefficient for " +simulation + " simulation", 
                               ylabel="$rho (t)$", xlabel="t", xlim=[0,800])

    plotFileName = pict_path + "rho_"+ St + "_" + simulation+".eps"

    for key,arg in zip(datalines.keys(),datalines.values()):
        corrfig.add_plot(rho[St + "_" + simulation][int(key)],
                         linestyle=arg["ls"],color = arg["c"], label= arg["lbl"])
	
    corrfig.hdraw()
    corrfig.save(plotFileName)
    print "plot created: " + plotFileName
    plt.close(corrfig.fig)



for simulation in ptype:
    # correlation time plot
    taufig = hfig.Homfig(title= " Relaxation time for particles in " +simulation , 
                               ylabel="$/tau$", xlabel="y", xlim=[0,8])

    plotFileName = pict_path + "tau_"+ simulation+".eps"

    for St in Stlist:
        taufig.add_plot(tau[St + "_" + simulation],label=St)
	
    taufig.hdraw()
    taufig.save(plotFileName)
    print "plot created: " + plotFileName
    plt.close(taufig.fig)


# plots of correlation coefficient DNS vs a priori LES

for St in Stlist:
    corrfig = hfig.Homfig(title="Correlation coefficient of $u^{*}$ for " + St, ylabel="$rho (t)$",xlabel="t",xlim=[0,800])
    plotFileName = pict_path + "rho_"+St+".eps"

    
    for key,simulation in product(datalines.keys(),ptype):
        arg = datalines[key]
        corrfig.add_plot(rho[St + "_" + simulation][key],
                         linestyle=arg["ls"],color = arg["c"], label= arg["lbl"])

    corrfig.hdraw()
    corrfig.save(plotFileName)
    print "plot created: " + plotFileName
    plt.close(corrfig.fig)
# plots of relaxation time

