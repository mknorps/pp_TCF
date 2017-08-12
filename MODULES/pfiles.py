# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# File name: apriori_SGS_fluid.py
# Created by: mknorps 
# Creation date: 21-06-2017
# Last modified: 12-08-2017 18:40:00
# Purpose: take filtered and unfiltered fluid field 
#          in Fourier space from spectral code 
#          compute statistics of SGS fluid velocity
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import numpy as np
import particlestat as ps
from itertools import ifilter,ifilterfalse



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
            stats_PT[''.join(arg[2])+arg[0]] = list(np.zeros(ps.Particles.Nnodes/2+1))


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

                    dictarg =''.join(arg[2]) +arg[0]
                    stats_PT[dictarg] = stats_PT[dictarg] + Pdata.stat_symm(stattype,symmtype,*statargs)[1]

        for i in ifilterfalse(lambda x: x=="yplus", stats_PT.keys()):

            stats_PT[i] = stats_PT[i]/float(self.nFiles)
            

        return stats_PT 

#TODO - reduce the spaghetti code
    #method for performing some computations (function f) on the input (*args)
    # args is a list of list of arguments of ndarray type, 
    # we want to compare several 
    # different input for the same function
    # for example f(x,y,z) = (x-y)*z
    def equationP(self,StNo,f,stattype,symmetryType='none',*args): 
        
        # initialise vector containing particle time statistics
        stats_PT = {}


        for argl in args:

            if symmetryType <> 'none':
                stats_PT[''.join(argl)] = list(np.zeros(ps.Particles.Nnodes/2+1))
            else:
                stats_PT[''.join(argl)] = list(np.zeros(ps.Particles.Nnodes))


        for name in ifilter(lambda x: x.endswith(StNo),self.fNames): #loop over time steps

            with open(name,'r') as current_file:

                raw_data = np.transpose(np.loadtxt(current_file))

                # setting up data structure
                spaceCoordinates  = [self.kwargs["x"],self.kwargs["y"],self.kwargs["z"]]
                kwarglist= {}

                #computing function of input raw data arguments
                for arglist in args:
                    kwarglist[''.join(arglist)]= f(*map(lambda x: raw_data[self.kwargs[x]],arglist)) 

                # tuple is immutable, list is mutable
                # mean_T is gathering mean values from consecutive timesteps, 
                #        so it has to be a list
                Pdata = ps.Particles(*map(lambda xx: raw_data[xx],spaceCoordinates),**kwarglist) #object of class Particles

                if symmetryType <> 'none':
                    stats_PT["yplus"] = Pdata.y_nondim()  #nodes in y direction
                else:
                    stats_PT["yplus"] = Pdata.ynodes()  #nodes in y direction


                for kwlistkey,kwlistval in kwarglist.iteritems(): #loop over choosen statistics and variables

                    if symmetryType <> 'none':
                        stats_PT[kwlistkey] = stats_PT[kwlistkey] + Pdata.stat_symm(stattype,symmetryType,kwlistkey)[1]
                    else:
                        stats_PT[kwlistkey] = stats_PT[kwlistkey] + Pdata.statistics(stattype,kwlistkey)[1]

#                    print kwarglist.keys(), kwlistkey, type(stats_PT[kwlistkey])

        for i in ifilterfalse(lambda x: x=="yplus", stats_PT.keys()):


            stats_PT[i] = stats_PT[i]/float(self.nFiles)
            

        return stats_PT 

