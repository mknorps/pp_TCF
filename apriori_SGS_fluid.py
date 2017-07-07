# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# File name: apriori_SGS_fluid.py
# Created by: mknorps 
# Creation date: 21-06-2017
# Last modified: 07-07-2017 15:31:19
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
test_set_DNS = hf.ChannelFields(128,129,128,2501,2502, fPrefix=["upp_Ux","upp_Uy","upp_Uz"])
test_set_LES = hf.ChannelFields(32,33,64,2501,2502, fPrefix=["uppf_Ux","uppf_Uy","uppf_Uz"])


for skey,sval in statistics.iteritems():

    print skey,sval,type(sval[0])
    #computation of statistics
    stat_DNS = test_set_DNS.statsT(*sval)
    stat_LES = test_set_LES.statsT(*sval)

    #figures
    for a in range(1,len(stat_DNS)) : #len(stat_DNS) = len(stat_LES)
        #handle for different plot title
        if skey in range(2,5) :
            statfig = hfig.Homfig(title=ptitles[sval[0]] + sval[1] 
                    + " " + sval[2], ylabel=ylabels[skey+1]) #correlations
        elif skey == 5:
            statfig = hfig.Homfig(title=ptitles[sval[0]], ylabel=ylabels[6]) #kinetic energy
        else:
            statfig = hfig.Homfig(title=ptitles[sval[0]]+coordinates[a], ylabel=ylabels[a])

        statfig.add_plot(stat_DNS[0],stat_DNS[a],linestyle=LineStyle['DNS'],label='DNS')
        statfig.add_plot(stat_LES[0],stat_LES[a],linestyle=LineStyle['LES'],label='LES')
        statfig.hdraw()
        statfig.save(''.join(sval) + ".eps")

    #meanfig.hdraw()
    plt.close(statfig.fig)
