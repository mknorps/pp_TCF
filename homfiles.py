# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# File name: apriori_SGS_fluid.py
# Created by: mknorps 
# Creation date: 21-06-2017
# Last modified: 06-07-2017 09:52:37
# Purpose: take filtered and unfiltered fluid field 
#          in Fourier space from spectral code 
#          compute statistics of SGS fluid velocity
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import numpy as np
import homstat as hs



# class for loading data files
# naming convension is: 'name'_number
class ChannelFields:

    # (K,N,M) are number of nodes in streamwise, wall-normal 
    #                            and spanwise directions respectively
    def __init__(self,K,N,M,fstart,fend,fPrefix=["upp_Ux","upp_Uy","upp_Uz"]):

        self.nFiles = fend-fstart + 1
        self.fStart = fstart
        self.fEnd   = fend
        self.fNames = []

        for x in range(fstart,fend+1):
            self.fNames.append(map(lambda y: y+"_"+str(x),fPrefix))

            # original input has permuted directions 
            # (left - input, right - standard naming convention, applied here)
            # (Ux,Uy,Uz) --> (Uy,Uz,Ux)
            with open(self.fNames[x-fstart][0],'r') as Uy, \
                    open(self.fNames[x-fstart][1],'r') as Uz, \
                    open(self.fNames[x-fstart][2],'r') as Ux :

	    # input data is in a matrix, and we first have to reshape 
	    # it to be a 3D field
	    # then directions must be permuted (np.transpose)
            # attribute of type Channel is added

		 setattr(self,"field_"+str(x),
                        hs.Channel(np.transpose(np.loadtxt(Ux).reshape(N,M,K),axes=(2,0,1)),
                        np.transpose(np.loadtxt(Uy).reshape(N,M,K),axes=(2,0,1)),
                        np.transpose(np.loadtxt(Uz).reshape(N,M,K),axes=(2,0,1)),K,N,M))
                
    # statistics averaged over time (T)
    def statsT(self,stats,*args):
       
        # tuple is immutable, list is mutable
        # mean_T is gathering mean values from consecutive timesteps, 
        #        so it has to be a list
        stats_T = list(getattr(getattr(self,"field_"+str(self.fStart)),stats)(*args))

        for attrNo in range(self.fStart+1,self.fEnd+1): 
            stats_tmp =getattr(getattr(self,"field_"+str(attrNo)),stats)(*args)
            for y in range(1,len(stats_T)):
                stats_T[y] = stats_T[y] + stats_tmp[y]

	for y in range(1,len(stats_T)):
	    stats_T[y] = stats_T[y]/float(self.nFiles)
        
        return stats_T 



