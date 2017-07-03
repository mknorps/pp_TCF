# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# File name: apriori_SGS_fluid.py
# Created by: mknorps 
# Creation date: 21-06-2017
# Last modified: 03-07-2017 16:42:58
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


test_set = hf.ChannelFields(32,33,64,2501,2505, fPrefix=["uppf_Ux","uppf_Uy","uppf_Uz"])

mean = test_set.mean_symmT()

print mean[0]
print mean[1]

meanfig = hfig.Homfig((mean[0][:],mean[1][:]),title='test plot')
meanfig.hdraw()
meanfig.save("testa.eps")
plt.close(meanfig.fig)
