# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# File name: timescales.py
# Created by: gemusia
# Creation date: 26-07-2017
# Last modified: 27-07-2017 12:50:00
# Purpose: draw plots of time scales of turbulent channel flow
#          DNS, a priori LES, LES and SGS
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


import numpy as np
import matplotlib.pyplot as plt
import homfigs as hfig
from os.path import expanduser
from itertools import ifilterfalse



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


ptype = ["DNS","LES","apriori","SGSles","SGSdns"]
Stlist = ["fluid","St0.1","St1","St5","St25"]



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
# drawing plots
#--------------------------------------------------------------------#

corr_files = []
corr_files_symm = []

for i in range(4):
    corr_files.append(np.transpose(np.loadtxt(file_path + "rho_Usgs_ptype" +str(i+1) +".dat",skiprows=1)))
    # tau_files  = np.loadtxt(file_path + "Tij_" +str(i+1) +".dat")


    # plots of correlation coefficient

    corrfig = hfig.Homfig(title="Correlation coefficient for " + ptype[i], ylabel="$rho (t)$",xlabel="t",xlim=[0,600])
    plotFileName = pict_path + "rho"+str(i+1)+".eps"


    corrfig.add_plot(corr_files[i][1],linestyle='solid',color = 'green', label='y=1')
    corrfig.add_plot(corr_files[i][16],linestyle='solid',color = 'red', label='y=16')
    corrfig.add_plot(corr_files[i][3],linestyle='dotted',color = 'green', label='y=3')
    corrfig.add_plot(corr_files[i][14],linestyle='dotted',color = 'red', label='y=14')
    corrfig.add_plot(corr_files[i][6],linestyle='dashdot',color = 'green', label='y=6')
    corrfig.add_plot(corr_files[i][11],linestyle='dashdot',color = 'red', label='y=11')
    corrfig.add_plot(corr_files[i][9],linestyle='dashed', label='y=9')
	
    corrfig.hdraw()
    corrfig.save(plotFileName)
    print "plot created: " + plotFileName
    plt.close(corrfig.fig)



    #symmetric plots
    corr_files_symm.append(symm(corr_files[i][1:]))
    print corr_files_symm[i][0]
    print corr_files_symm[i][1]
    
    corrfig_symm = hfig.Homfig(title="Correlation coefficient for " + ptype[i], ylabel="$rho (t)$",xlabel="t",xlim=[0,600])
    plotFileName = pict_path + "rho"+str(i+1)+"_symm.eps"


    corrfig_symm.add_plot(corr_files[i][0],linestyle='solid',color = 'green', label='y=1')
    corrfig_symm.add_plot(corr_files[i][2],linestyle='dotted',color = 'green', label='y=3')
    corrfig_symm.add_plot(corr_files[i][5],linestyle='dashdot',color = 'green', label='y=6')
    corrfig_symm.add_plot(corr_files[i][7],linestyle='dashed', label='y=9')
	
    corrfig_symm.hdraw()
    corrfig_symm.save(plotFileName)
    print "plot created: " + plotFileName
    plt.close(corrfig_symm.fig)


    # plots of correlation coefficient DNS vs a priori LES

for i,St in zip(range(2),Stlist):
    corrfig = hfig.Homfig(title="Correlation coefficient for " + St, ylabel="$rho (t)$",xlabel="t",xlim=[0,600])
    plotFileName = pict_path + "rho"+St+"DNSvsLES.eps"


    corrfig.add_plot(corr_files[2*i][1],linestyle='solid',color = 'green', label='y=1,DNS')
    corrfig.add_plot(corr_files[2*i][16],linestyle='solid',color = 'red', label='y=16, DNS')
    corrfig.add_plot(corr_files[2*i+1][1],linestyle='solid',color = 'blue', label='y=1, LES')
    corrfig.add_plot(corr_files[2*i+1][16],linestyle='solid',color = 'orange', label='y=16, LES')
    corrfig.add_plot(corr_files[2*i][3],linestyle='dotted',color = 'green', label='y=3, DNS')
    corrfig.add_plot(corr_files[2*i][14],linestyle='dotted',color = 'red', label='y=14, DNS')
    corrfig.add_plot(corr_files[2*i+1][3],linestyle='dotted',color = 'blue', label='y=3, LES')
    corrfig.add_plot(corr_files[2*i+1][14],linestyle='dotted',color = 'orange', label='y=14, LES')
    corrfig.add_plot(corr_files[2*i][9],linestyle='dashed', color='green',label='y=9, DNS')
    corrfig.add_plot(corr_files[2*i+1][9],linestyle='dashed',color='orange', label='y=9, LES')
	
    corrfig.hdraw()
    corrfig.save(plotFileName)
    print "plot created: " + plotFileName
    plt.close(corrfig.fig)
# plots of relaxation time

