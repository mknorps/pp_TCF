# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# File name: timescales.py
# Created by: gemusia
# Creation date: 26-07-2017
# Last modified: 15-09-2017 20:23:41
# Purpose: draw plots of time scales of turbulent channel flow
#          DNS, a priori LES, LES and SGS
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

#--------------------------------------------------------------------#
# helper functions
#--------------------------------------------------------------------#
def symm(lst):
    ll = len(lst)
    return 0.5*(lst[:(ll+1)/2] + np.flipud(lst)[:(ll+1)/2])

def Chebyshev_zeroes(N):
    return 0.5*(np.array(range(1,N+1))+np.array(range(N)))/float(N) * 150


print "cheb zeroes: ",Chebyshev_zeroes(8)
print Chebyshev_zeroes(16)
#--------------------------------------------------------------------#
# declaration of picture attributes
#--------------------------------------------------------------------#

LineStyle = {'fterm':'solid','pterm':'dashed',
        "Ux":'solid',"Ufx":'dashed',"Vx":'dotted',
        "UxUfx":'solid',"UyUfy":'dashed',"UzUfz":'dotted',
        "fluid":'solid',"St1":'dashed',"St5":'dotted',"St25":'dashdot'}

ptype_dt = {"LES":dt_LES,"SGSles":dt_DNS,"SGSdns":dt_DNS,"DNS":dt_DNS, "apriori":dt_DNS}

ptype = ["DNS","LES","apriori","SGSles","SGSdns"]
#ptype = ["LES","DNS","SGSles","SGSdns"]

#Stlist = ["St0.5","St5","St25","St125"]
#Stlist = ["fluid","St0.5","St5","St25","St125"]
Stlist = ["fluid","St0.5","St1","St5","St25","St125"]

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
                               ylabel="$\\tau^{+}$", xlabel="$y^{+}$") 
    
    taufig_related2fluid = hfig.Homfig(title= " Relaxation time for particles in " +simulation ,
                               ylabel="$\\frac{\\tau}{\\tau_f}$", xlabel="$y^{+}$")

    plotFileName = pict_path + "tau_"+ simulation+".eps"
    plotFileName_related2fluid = pict_path + "tau_"+ simulation+"_related2fluid.eps"

    for St in Stlist:
        taufig.add_plot(x_points,tau[St + "_" + simulation],label=St)
        taufig_related2fluid.add_plot(x_points,np.divide(tau[St + "_" + simulation],tau["fluid_" + simulation]),label=St)

    taufig.hdraw()
    taufig_related2fluid.hdraw()
    taufig.save(plotFileName)
    taufig_related2fluid.save(plotFileName_related2fluid)
    print "plots created: " + plotFileName + " and " + plotFileName_related2fluid
    plt.close(taufig.fig)
    plt.close(taufig_related2fluid.fig)


for St in Stlist:
    # correlation time plot
    taufig = hfig.Homfig(title= " Relaxation time for particles " + St ,
                               ylabel="$\\tau^{+}$", xlabel="$y^{+}$")
    taufig_related2DNS = hfig.Homfig(title= " Relaxation time for particles " + St ,
                               ylabel="$\\frac{\\tau}{\\tau_{DNS}}$", xlabel="$y^{+}$")

    plotFileName = pict_path + "tau_"+ St +".eps"
    plotFileName_related2DNS = pict_path + "tau_"+ St +"_related2DNS.eps"

    for simulation in ptype:
        taufig.add_plot(x_points,tau[St + "_" + simulation],label=simulation)
        taufig_related2DNS.add_plot(x_points,np.divide(tau[St + "_" + simulation],tau[St + "_DNS"]),label=simulation)

    taufig.hdraw()
    taufig_related2DNS.hdraw()
    taufig.save(plotFileName)
    taufig_related2DNS.save(plotFileName_related2DNS)
    print "plot created: " + plotFileName + "  " + plotFileName_related2DNS
    plt.close(taufig.fig)
    plt.close(taufig_related2DNS.fig)


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


#SGS contribution related to DNS relaxation time
taufig = hfig.Homfig(title= " Relaxation time for particles in " +simulation ,
                           ylabel="$\\frac{\\tau}{\\tau_DNS}$", xlabel="$y^{+}$")

plotFileName = pict_path + "tau_SGS_normalised.eps"

for St in Stlist:
    taufig.add_plot(x_points,np.divide(tau[St + "_SGSles"],tau[St + "_DNS"]),label=St)

taufig.hdraw()
taufig.save(plotFileName)
print "plots created: " + plotFileName 
plt.close(taufig.fig)

#LES contribution related to DNS relaxation time
taufig = hfig.Homfig(title= " Relaxation time for particles in " +simulation ,
                           ylabel="$\\frac{\\tau}{\\tau_DNS}$", xlabel="$y^{+}$")

plotFileName = pict_path + "tau_LES_normalised.eps"

for St in Stlist:
    taufig.add_plot(x_points,np.divide(tau[St + "_LES"],tau[St + "_DNS"]),label=St)

taufig.hdraw()
taufig.save(plotFileName)
print "plots created: " + plotFileName 
plt.close(taufig.fig)
