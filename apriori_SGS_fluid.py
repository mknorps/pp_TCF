# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# File name: apriori_SGS_fluid.py
# Created by: mknorps 
# Creation date: 21-06-2017
# Last modified: 23-06-2017 20:18:51
# Purpose: take filtered and unfiltered fluid field 
#          in Fourier space from spectral code 
#          compute statistics of SGS fluid velocity
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import numpy as np
import matplotlib.pyplot as plt
import homstat as hs



#Files names
fName  = "upp_"
fNamef = "uppf_"
start  = 2501
end    = 2532
nFiles = end-start+1

#Initialisation of variables
upp   = [[]]*(3*(end-start+1))
uppf  = [[]]*(3*(end-start+1))
sdict = {}

#Constants
n  = 128
m  = 128
o  = 128
nf = 32
mf = 64
of = 32

#Initialisation of statistics we want to compute
#stats=dict()
upp_mean  = [np.zeros(n+1),np.zeros(n+1),np.zeros(n+1)]
uppf_mean = [np.zeros(nf+1),np.zeros(nf+1),np.zeros(nf+1)]


#files with raw data
for i in range(start,end):
    j=0
    for s in ["Ux","Uy","Uz"]:
        print ("i = " + str(i) + " j = " + str(j))
        with open(fName + s + "_"+ str(i),'r') as fupp:
                upp[i-start+j]=np.loadtxt(fupp).reshape(n+1,m,o)
        with open(fNamef + s + "_"+ str(i),'r') as fuppf:
                uppf[i-start+j]=np.loadtxt(fuppf).reshape(nf+1,mf,of)
        upp_mean[j]  = upp_mean[j] + hs.hmean(upp[i-start+j])
        uppf_mean[j] = uppf_mean[j] + hs.hmean(uppf[i-start+j])
        j=j+1

for j in range(0,2):
        upp_mean[j]  = upp_mean[j]/float(nFiles)
        uppf_mean[j] = uppf_mean[j]/float(nFiles)

print upp_mean[2]
print uppf_mean[2]

