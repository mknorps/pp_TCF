# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# File name: apriori_SGS_fluid.py
# Created by: mknorps 
# Creation date: 21-06-2017
# Last modified: 10-07-2017 22:39:39
# Purpose: take filtered and unfiltered fluid field 
#          in Fourier space from spectral code 
#          compute statistics of SGS fluid velocity
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import numpy as np
import particlestat as ps
from itertools import ifilter



# class for loading data files
# naming convension is: 'name'_number
class ParticleFields:

    def __init__(self,fstart,fend,fCoreName="fede_terms_",fSuffix=["St1","St5","St25","fluid"],**kwargs):

        self.nFiles = fend-fstart + 1
        self.fStart = fstart
        self.fEnd   = fend
        self.fNames = []
        self.fSuffix = fSuffix
        self.kwargs = kwargs # dictionary keeping variable name and column in which it is located

        for x in range(fstart,fend+1):
            self.fNames = self.fNames+ map(lambda y: fCoreName +str(x)+"_"+y,fSuffix)

                
    # statistics averaged over time (T)
    # symmetrised
    # it computes statistic of name 'stattype' with arguments *statargs
    # and symmetrisation type 'symm'
    # for all type of particles

    def statsP(self,StNo,*args):
     
        # initialise vector containing particle time statistics
        stats_PT = {}
        for arg in args:
            stats_PT[arg[0]+''.join(arg[2])] = list(np.zeros(ps.Particles.Nnodes/2+1))


        for name in ifilter(lambda x: x.endswith(StNo),self.fNames): #loop over time steps


            with open(name,'r') as current_file:

                raw_data = np.transpose(np.loadtxt(current_file))

                # setting up data structure
                arglist  = [self.kwargs["x"],self.kwargs["y"],self.kwargs["z"]]
                kwarglist= {}
                for key,val in ifilter(lambda x: x[0] not in ["x","y","z"] ,self.kwargs.iteritems()):
                    kwarglist[key]=raw_data[val]

                # tuple is immutable, list is mutable
                # mean_T is gathering mean values from consecutive timesteps, 
                #        so it has to be a list
                Pdata = ps.Particles(*map(lambda i: raw_data[i],arglist),**kwarglist) #object of class Particles

                stats_PT["yplus"] = Pdata.y_nondim()  #nodes in y direction


                for cntr,arg in enumerate(args): #loop over choosen statistics and variables

                    stattype = arg[0] #string
                    symmtype = arg[1] #string
                    statargs = arg[2] #list

                    dictarg =[arg[0]+''.join(arg[2])] 
                    stats_PT[dictarg] = stats_PT[dictarg] + Pdata.stat_symm(stattype,symmtype,*statargs)[1]

        for i in stats_PT.keys():

            stats_PT[i] = stats_PT[i]/float(self.nFiles)
            

        return stats_PT 



