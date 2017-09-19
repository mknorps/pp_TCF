# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# File name: timescales.py
# Created by: gemusia
# Creation date: 26-07-2017
# Last modified: 19-09-2017 15:02:47
# Purpose: draw plots of time scales of turbulent channel flow
#          DNS, a priori LES, LES and SGS
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import homfigs as hfig
import timeconstants as tc

from os.path import expanduser,isfile
from itertools import product


#--------------------------------------------------------------------#
# Data reading
#--------------------------------------------------------------------#


f_tau     = tc.file_path + "tau.csv"
f_epsilon = tc.file_path + "diss_mean_5000.dat"

#load time scales from a priori file
with open(f_tau,'r') as f:
    relaxation_time = pd.read_csv(f,sep = ',',header=0)

#load dissipation and SGS kinetic energy file
with open(f_epsilon,'r') as f:
    dissipation = np.genfromtxt(f)
    dissipation = pd.DataFrame(dissipation)
    dissipation.columns = ['y','unknown','e1','ksgs','e3']
    print(dissipation.dtypes )



#--------------------------------------------------------------------#
# Data manipulation 
#--------------------------------------------------------------------#

dissipation['tau_modell1'] = -1.0 * np.divide(dissipation['e1'],dissipation['ksgs'])
dissipation['tau_modell3'] = np.divide(dissipation['e3'],dissipation['ksgs'])
dissipation['tau_modell2'] = np.divide(np.array(map(tc.filter_width,range(len(dissipation)))),
                                      map(lambda x: x**(0.5), dissipation['ksgs']))


print(dissipation)
# change: NaN to None
dissipation = dissipation.where(dissipation.notnull(),None)

#symmetrisation of columns needed for plots
dissipation_symm = pd.DataFrame(tc.non_dimensialisation(dissipation['y'][:(len(dissipation['y'])+1)//2]), columns = ['y'])
models =  ['tau_modell1','tau_modell2','tau_modell3' ]
for col in models:
    dissipation_symm[col] = tc.symm(dissipation[col])/tc.ttau

dissipation_symm = dissipation_symm.where(dissipation_symm.notnull(),None)

print(dissipation_symm)
#y axis tikz
relaxation_time['y'] = tc.Chebyshev_zeroes(len(relaxation_time))

# normalisation by factor tc.ttau
columns_without_y = [x for x in relaxation_time.columns if x != 'y']
for col in columns_without_y:
    relaxation_time[col] = relaxation_time[col]


#--------------------------------------------------------------------#
# PLOT DRAW
#
# we use symmetrised data
#--------------------------------------------------------------------#


if __name__ == '__main__':
    
    SGStime = hfig.Homfig(title= " SGS relaxation time for particles" ,
                               ylabel="$\\tau$", xlabel="$y^{+}$")

    plotFileName = tc.pict_path + "SGS_aprioriVSmodel.eps"

    SGStime.add_plot(relaxation_time['y'],relaxation_time['fluid_SGSles'],label='a priori')
    for model in models:
        SGStime.add_plot(dissipation_symm['y'],dissipation_symm[model],label=model)

    SGStime.hdraw()
    SGStime.save(plotFileName)
    plt.close(SGStime.fig)

