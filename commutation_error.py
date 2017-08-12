# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# File name: commutation_error.py
# Created by: gemusia
# Creation date: 02-08-2017
# Last modified: 03-08-2017 12:38:33
# Purpose:draws commutation error coefficients 
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


import numpy as np
import matplotlib.pyplot as plt
import homfigs as hfig
from os.path import expanduser


n = 32

def f(x):
    return (np.sin((x+1)/32.0*np.pi)*np.sin((x-1)/32.0*np.pi) - (np.sin(x/32.0*np.pi))**2)/(np.sin(x/32.0*np.pi)*np.sin((x-1)/32.0*np.pi)*(np.sin(x/32.0*np.pi)+np.sin((x-1)/32.0*np.pi)))

def g(x):
    return 1/6.0*np.pi*(np.sin(np.pi/64.0))**2*(np.sin(2*np.pi*x))

def h(x):
    return np.sin(np.pi/(2*n))*np.cos(x/float(n)*np.pi)/( np.cos(np.pi/(2*n))*np.sin(x/float(n)*np.pi))

y = np.arange(2,31,1)
z = np.arange(0,1,0.1)

print z
print g(z)

file_path = expanduser("~") + "/Documents/Doktorat_current/picts/"
pict_path = file_path

asymmetry_coeff = hfig.Homfig(title= "Difference of right and left side wall-normal filter", 
                           ylabel="$c$", xlabel="y", xlim=[1,31])
                         #  ylabel="$c$", xlabel="y", xlim=[0,1])

plotFileName = pict_path + "filter_skewness.eps"

asymmetry_coeff.add_plot(h(y))
                         
	
asymmetry_coeff.hdraw()
asymmetry_coeff.save(plotFileName)
print "plot created: " + plotFileName
plt.close(asymmetry_coeff.fig)
