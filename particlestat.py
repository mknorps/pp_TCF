# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# File name: particlestat.py
# Created by: mknorps
# Creation date: 07-07-2017
# Last modified: 08-07-2017 21:43:09
# Purpose: module for computing statistics of 
#   particles in turbulent channel flow.
#
#   Statistics are taken by averaging over homogeneity directions x and z
#   and are presented a functions of y
#
#   Statistics may be presented in non-dimensional form as functions of y^+
#
#   Input np.array is a 2D array.
#   for every particle we have following values:
#    (x,y,z,Vx,Vy,Vz,
#     Ux,Uy,Uz,Ufx,Ufy,Ufz,
#     dUx/dx,dUx/dy,dUx/dz,dUy/dx,dUy/dy,dUz/dz,dUz/dx,dUz/dy,dUz/dz,
#     dUfx/dx,dUfx/dy,dUfx/dz,dUfy/dx,dUfy/dy,dUfz/dz,dUfz/dx,dUfz/dy,dUfz/dz)
#
#    where "f" denote filtered (LES) velocties
#   
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


import numpy as np

#~~~~~~~~~~ Particles ~~~~~~~~~~
#   - class for paritcles in wall-bounded turbulent channel flow
#
# Directions of the flow:
# x - streamwise 
# y - wall-normal
# z - spanwise
#
# Number of paricles - N
# Number of bins in wall-normal direction - Nbins


# zeroes of Chebyshev polynomials
# usefull for creating bins in wall-normal direction
def ChebZeros(N):
    def ChebCurr(y):
        return np.cos(y*np.pi/N)
    return ChebCurr


#symmetrisation of 1D nparray
#useful for presenting statistics in most readable format
#input may be symmetric or asymmetric function
def symm(mode,lst):
    opt = {"symm":1.0,"asymm":-1.0}
    ll = len(lst)

    return 0.5*(opt[mode] * lst[:(ll+1)/2] + np.flipud(lst)[:(ll+1)/2])



#class Particles:
class Particles:

    # homogeneous directions axis
    ha=(0,2)

    #constats - change is rarely neede
    Retau=150
    nuinv=3500
    Nbins = 32
    Nnodes = Nbins+1

    # dictionary of possible statistics
    stat_dict = {"mean":(lambda x: x),
            "mean_sqr":(lambda x: x**2),
            "cov":(lambda x,y:x*y),
            "ke":(lambda x,y,z: x**2+y**2+z**2)}



    #friction velocity used for non-dimensionalisation
    def utau(self):
        return float(self.Retau)/float(self.nuinv)


    # for object initiation we need velocities in three directions (Ux,Uy,Uz)
    # it should be given as three numpy arrays of size (K,N,M)
    def __init__(self,x,y,z,Vx,Vy,Vz,**kwargs):

        self.x     = np.array(x)
        self.y     = np.array(y) #wall-normal direction
        self.z     = np.array(z)
        self.N     = len(self.x) #nuber of particles
        self.Vx    = np.array(Vx)
        self.Vy    = np.array(Vy)
        self.Vz    = np.array(Vz)
        self.ibin  = map(lambda xx: int((Particles.Nbins/np.pi)*np.arccos(xx)),self.y) #number of bin in which a particle is located


        self.pbin = [list() for _ in range(Particles.Nbins+1)]

        for particle,yp in enumerate(self.y):
            self.pbin[self.ibin[particle]].append(particle)    #list of particles located in a separate bin


        for key,val in kwargs.iteritems():
            setattr(self,key,val)

    # nodes are taken as zeroes of Chebyshev pomynomials
    def ynodes(self):
        return np.array(map(ChebZeros(Particles.Nnodes-1),range(Particles.Nnodes)))

    # symmetrised nodes in nondimensional notation
    # if N is even, we have N/2 nodes
    # if N is odd, we have N/2+1 nodes
    def y_nondim(self):
        y = self.ynodes()[:(Particles.Nnodes+1)/2]
        return np.array(map(lambda x: (1.0-x)*self.Retau,y))

    # cell centers in y direction, non-dimensionalized
    def y_centers(self):
        y = self.y_nondim()
        return 0.5*(y[:len(y)-1]+y[1:])

#...............................................................
#      STATISTICS
#...............................................................

# Symmetrised statistics ar presented in non-dimensional form

    #one-point  statistics
    def stat1P (self,stat_type,*args):

        var_list = []
        statInBin = np.zeros(Particles.Nbins+1)

        for arg in args:
            var_list.append(getattr(self,arg))

        val = zip(*var_list)
        

        for particle,_ in enumerate(self.y):
            ibinp = self.ibin[particle]  # bin where particle 'particle' is located
            statInBin[ibinp] = statInBin[ibinp] + stat_type(*val[particle])


        # to reduce number of divisions, we divide statistics by particle number in bin in a separate loop

        for ibin in range(Particles.Nbins+1):
            np_ibin = len(self.pbin[ibin]) #number of particles in ith bin
            if np_ibin == 0:
                statInBin[ibin] = 0.0
            else:
                statInBin[ibin] = statInBin[ibin]/float(np_ibin)

        return statInBin



    def pmean (self,arg):
        return (self.ynodes(),self.stat1P(Particles.stat_dict["mean"],arg))

    def pvar (self,arg):
        return (self.ynodes(), self.stat1P(Particles.stat_dict["mean_sqr"],arg) - (self.stat1P(Particles.stat_dict["mean"],arg))**2)
                
    def pstd (self,arg):
        return  (self.ynodes(),np.sqrt( self.pvar(arg)[1] ))

    def pcov (self,arg1,arg2):
        return  (self.ynodes(),self.stat1P(Particles.stat_dict["cov"],arg1,arg2)
                 - self.stat1P(Particles.stat_dict["mean"],arg1)*self.stat1P(Particles.stat_dict["mean"],arg2))

    def pke (self,arg1,arg2,arg3):
        return  (self.ynodes(),self.stat1P(Particles.stat_dict["ke"],arg1,arg2,arg3) )



    #symmetrisation of chosen statistic (stat: pmean, pvar, pstd, pcov,pke) 
    # for given argument (for example velocity)
    def stat_symm (self,stat,mode,*args):
        return (self.y_nondim(),symm(mode,getattr(self,stat)(*args)[1]))
