# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# File name: timescales.py
# Created by: gemusia
# Creation date: 26-07-2017
# Last modified: 20-09-2017 16:16:51
# Purpose: draw plots of time scales of turbulent channel flow
#          DNS, a priori LES, LES and SGS
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import homfigs as hfig
import timeconstants as tc


#--------------------------------------------------------------------#
# Data reading
#--------------------------------------------------------------------#


f_tau     = tc.file_path + "tau.csv"
f_epsilon = tc.file_path + "diss_mean_5000.dat"
f_ksgs    = tc.file_path + "covu.dat"

#load time scales from a priori file
with open(f_tau,'r') as f:
    relaxation_time = pd.read_csv(f,sep = ',',header=0)

#load dissipation and SGS kinetic energy file
with open(f_epsilon,'r') as f:
    dissipation = np.genfromtxt(f)
    dissipation = pd.DataFrame(dissipation)
    dissipation.columns = ['y','unknown','e1','ksgs','e3']


#load  kinetic energy in LES (computed with Yoshizawa estimation) 
with open(f_ksgs,'r') as f:
    ksgsYoshizawa = np.genfromtxt(f)
    ksgsYoshizawa = pd.DataFrame(ksgsYoshizawa)
    ksgsYoshizawa.columns = ['y','u2sgsy','u2sgsz','u2sgsx','usgsyz','usgsxy','usgsxz']


#--------------------------------------------------------------------#
# Data manipulation 
#--------------------------------------------------------------------#

dissipation['tau_modell1'] = -1.0 * np.divide(dissipation['e1'],dissipation['ksgs'])/tc.ttau
dissipation['tau_modell3'] = np.divide(dissipation['e3'],dissipation['ksgs'])/tc.ttau
dissipation['tau_modell2'] = np.divide(np.array(map(tc.filter_width,range(len(dissipation)))),
                                      map(lambda x: (2.0/3.0)*x**(0.5), dissipation['ksgs']))/(tc.ttau*tc.C_epsilon)

ksgsYoshizawa['ksgs'] = 0.5* (ksgsYoshizawa['u2sgsx'] + ksgsYoshizawa['u2sgsy'] + ksgsYoshizawa['u2sgsz'])
ksgsYoshizawa['tau_Yoshizawa'] = np.divide(np.array(map(tc.filter_width,range(len(ksgsYoshizawa)))),
                                      map(lambda x: (2.0/3.0)*x**(0.5),ksgsYoshizawa['ksgs']))/tc.ttau
ksgsYoshizawa['ksgs'] = ksgsYoshizawa['ksgs']/(tc.utau**2)

# change: NaN to None
dissipation = dissipation.where(dissipation.notnull(),None)
ksgsYoshizawa = ksgsYoshizawa.where(ksgsYoshizawa.notnull(),None)

# symmetrisation of columns needed for plots
# normalisation by factor tc.ttau or tc.utau for velocity
dissipation_symm = pd.DataFrame(tc.non_dimensialisation(dissipation['y'][:(len(dissipation['y'])+1)//2]), columns = ['y'])
for col in dissipation.drop('y',axis = 1).columns:
    dissipation_symm[col] = tc.symm(dissipation[col])

dissipation_symm = dissipation_symm.where(dissipation_symm.notnull(),None)

ksgsYoshizawa_y = np.array([tc.y(j) for j in  ksgsYoshizawa['y']]) # there is enumeration, not y coordinates
ksgsYoshizawa_symm = pd.DataFrame(tc.non_dimensialisation(ksgsYoshizawa_y[:(len(ksgsYoshizawa['y'])+1)//2]), columns = ['y'])

for col in ksgsYoshizawa.drop('y',axis = 1).columns:
    ksgsYoshizawa_symm[col] = tc.symm(ksgsYoshizawa[col])#/((tc.utau)**2)

ksgsYoshizawa_symm = ksgsYoshizawa_symm.where(ksgsYoshizawa_symm.notnull(),None)

#y axis tikz
relaxation_time['y'] = tc.Chebyshev_zeroes(len(relaxation_time))


print(relaxation_time.describe())
print(dissipation_symm.describe())
print(ksgsYoshizawa_symm.describe())

#--------------------------------------------------------------------#
# PLOT DRAW
#
# we use symmetrised data
#--------------------------------------------------------------------#


if __name__ == '__main__':


    #--------------------------------------------------------------------#
    #  time scales 
    #--------------------------------------------------------------------#

    SGStime = hfig.Homfig(title= " SGS relaxation time for fluid particles" ,
                               ylabel="$\\tau$", xlabel="$y^{+}$")

    plotFileName = tc.pict_path + "SGS_aprioriVSmodel.eps"

    SGStime.add_plot(relaxation_time['y'],relaxation_time['fluid_SGSles'],label='a priori')
    SGStime.add_plot(dissipation_symm['y'],dissipation_symm['tau_modell2'],label='tau_modell2')
    SGStime.add_plot(ksgsYoshizawa_symm['y'],ksgsYoshizawa_symm['tau_Yoshizawa'],label='LES,Yoshizawa est.')

    SGStime.hdraw()
    SGStime.save(plotFileName)
    plt.close(SGStime.fig)

    #--------------------------------------------------------------------#
    #  kinetic energy 
    #--------------------------------------------------------------------#
    ksgs = hfig.Homfig(title= " SGS kinetic energy for fluid particles" ,
                               ylabel="$\\tau$", xlabel="$y^{+}$")

    plotFileName = tc.pict_path + "ksgs_aprioriVSmodel.eps"

    ksgs.add_plot(ksgsYoshizawa_symm['y'],ksgsYoshizawa_symm['ksgs'],label='Yoshizawa est.')
    ksgs.add_plot(dissipation_symm['y'],dissipation_symm['ksgs'],label='apriori')

    ksgs.hdraw()
    ksgs.save(plotFileName)
    plt.close(ksgs.fig)
