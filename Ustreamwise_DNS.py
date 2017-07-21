# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# File name: apriori_SGS_fluid.py
# Created by: mknorps 
# Creation date: 21-06-2017
# Last modified: 21-07-2017 14:45:36
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
from os.path import expanduser



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

# U streamwise, for channel flow description in PhD

DNS_5200 = np.transpose(np.loadtxt('/home/gemusia/REFERENCE_DATA/DNS_5200_Uplus_LeeMoser15.txt.csv',skiprows=6))
statfig = hfig.Homfig(title="Streamwise velocity", ylabel='$U_{x}^{+}$',xscale='log',xlim=[1,6000])

print DNS_5200

statfig.add_plot(DNS_5200[0],DNS_5200[1],linestyle=LineStyle['DNS'],label='DNS')
statfig.add_plot(DNS_5200[0][:7],DNS_5200[0][:7],linestyle='dashed',label='$y=x$')
statfig.add_plot(DNS_5200[0][6:],2.49*np.log(DNS_5200[0][6:])+5,linestyle='dashdot',label='$y=a\ln(x) + b$')
plt.axvline(x=5,ymin=0,ymax=0.6,linestyle = 'dashed',color='black',linewidth=0.2)
plt.axvline(x=30,ymin=0,ymax=0.6,linestyle = 'dashed',color='black',linewidth=0.2)
plt.axvline(x=700,ymin=0,ymax=0.8,linestyle = 'dashed',color='black',linewidth=0.2)
plt.text(2,10,"I")
plt.text(10,15,"II")
plt.text(80,10,"III")
plt.text(1100,10,"IV")

statfig.hdraw()
statfig.save(fig_path + "Uzones.eps")

#meanfig.hdraw()
plt.close(statfig.fig)
