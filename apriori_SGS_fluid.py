# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# File name: apriori_SGS_fluid.py
# Created by: mknorps 
# Creation date: 21-06-2017
# Last modified: 07-07-2017 14:02:12
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



# declaration of picture attributes

LineStyle = {'DNS':'solid','LES':'dashed'}
coordinates = {1:'x',2:'y',3:'z'}
Ulist = {1:'Ux',2:'Uy',3:'Uz'}
statistics = {0:"hmean_symm",1:"hstd_symm",2:"hcor_symm",3:"hek_symm"}
#statistics = {0:"hmean_symm",1:"hstd_symm"}
ylabels = {1:"$<U>^{x}$",2:"$<U>^{y}$",3:"$<U>^{z}$"}
ptitles = {"hmean_symm":"mean of U","hstd_symm":"Standard deviation of U","hcor_symm":"Correlation of ","hek_symm":"Kinetic energy"}


#data loading
test_set_DNS = hf.ChannelFields(128,129,128,2501,2505, fPrefix=["upp_Ux","upp_Uy","upp_Uz"])
test_set_LES = hf.ChannelFields(32,33,64,2501,2505, fPrefix=["uppf_Ux","uppf_Uy","uppf_Uz"])


for skey,sval in statistics.iteritems():
    #computation of statistics
    statchoice = sval
    if skey == 2:  # for correlation we need to give correlting directions
        for x in it.combinations(Ulist.values(),2):
            stat_DNS = test_set_DNS.statsT(statchoice,*x)
            stat_LES = test_set_LES.statsT(statchoice,*x)
    else:
        stat_DNS = test_set_DNS.statsT(statchoice)
        stat_LES = test_set_LES.statsT(statchoice)

    #figures
    for a in range(1,len(stat_DNS)) : #len(stat_DNS) = len(stat_LES)
        statfig = hfig.Homfig(title=ptitles[statistics[a]]+str(a), ylabel=ylabels[a])
        statfig.add_plot(stat_DNS[0],stat_DNS[a],linestyle=LineStyle['DNS'],label='DNS')
        statfig.add_plot(stat_LES[0],stat_LES[a],linestyle=LineStyle['LES'],label='LES')
        statfig.hdraw()
        statfig.save(sval+coordinates[a]+".eps")

    #meanfig.hdraw()
    plt.close(statfig.fig)
