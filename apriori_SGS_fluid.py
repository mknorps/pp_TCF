# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# File name: apriori_SGS_fluid.py
# Created by: mknorps 
# Creation date: 21-06-2017
# Last modified: 05-07-2017 22:55:31
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



# declaration of picture attributes

LineStyle = {'DNS':'k-^','LES':'-o'}
coordinates = {1:'x',2:'y',3:'z'}
#statistics = {0:"hmean_symm",1:"hstd_symm",2:"hcor_symm",3:"hek_symm"}
statistics = {0:"hmean_symm",1:"hstd_symm"}
ylabels = {1:"$<U>^{x}$",2:"$<U>^{y}$",3:"$<U>^{z}$"}
ptitles = {"hmean_symm":"mean of U","hstd_symm":"Standard deviation of U","hcor_symm":"Correlation of ","hek_symm":"Kinetic energy"}


#data loading
test_set_DNS = hf.ChannelFields(128,129,128,2501,2505, fPrefix=["upp_Ux","upp_Uy","upp_Uz"])
test_set_LES = hf.ChannelFields(32,33,64,2501,2505, fPrefix=["uppf_Ux","uppf_Uy","uppf_Uz"])


for skey,sval in statistics.iteritems():
    #computation of statistics
    stat_DNS = test_set_DNS.statsT(sval)
    stat_LES = test_set_LES.statsT(sval)

    #figures
    for key,val in coordinates.iteritems():
        statfig = hfig.Homfig(title=ptitles[sval]+val, ylabel=ylabels[key])
        statfig.ax.plot(stat_DNS[0],stat_DNS[key],LineStyle['DNS'],label='DNS')
        statfig.ax.plot(stat_LES[0],stat_LES[key],LineStyle['LES'],label='LES')
        statfig.save(sval+val+".eps")
        leg = statfig.ax.legend(loc=4)

    #meanfig.hdraw()
    plt.close(statfig.fig)
