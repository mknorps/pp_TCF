# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# File name: homstat.py
# Created by: gemusia
# Creation date: 22-06-2017
# Last modified: 07-07-2017 22:31:13
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
    Retau=150
    nuinv=3500
    Nbins = 32

    # for object initiation we need velocities in three directions (Ux,Uy,Uz)
    # it should be given as three numpy arrays of size (K,N,M)
    def __init__(self,x,y,z,Vx,Vy,Vz,**kwargs):

        self.x     = np.array(x)
        self.y     = np.array(y)
        self.z     = np.array(z)
        self.N     = len(self.x) #nuber of particles
        self.Vx    = np.array(Vx)
        self.Vy    = np.array(Vy)
        self.Vz    = np.array(Vz)
        self.ibin  = np.array() #number of bin in which a particle is located

#TODO - finish bin assignment




        for key,val in kwargs.iteritems():
            setattr(self,key,val)

    # nodes are taken as zeroes of Chebyshev pomynomials
    def ynodes(self):
        return np.array(map(ChebZeros(self.N-1),range(self.N)))

    # symmetrised nodes in nondimensional notation
    # if N is even, we have N/2 nodes
    # if N is odd, we have N/2+1 nodes
    def y_nondim(self):
        y = self.ynodes()[:(self.N+1)/2]
        return np.array(map(lambda x: (1.0-x)*self.Retau,y))

    # cell centers in y direction, non-dimensionalized
    def y_centers(self):
        y = self.y_nondim()
        return 0.5*(y[:len(y)-1]+y[1:])

    #friction velocity used for non-dimensionalisation
    # here giben as parameter
    # TODO - compute it from the data
    def utau(self):
        return float(self.Retau)/float(self.nuinv)
#...............................................................
#      STATISTICS
#...............................................................

# Symmetrised statistics ar presented in non-dimensional form
    #mean
    def hmean (self):
        return (self.ynodes(),
                np.mean(self.Ux,axis=self.ha), 
                np.mean(self.Uy,axis=self.ha), 
                np.mean(self.Uz,axis=self.ha))

    #mean symmetrised
    def hmean_symm (self):
        return (self.y_nondim(),
                symm("symm",np.mean(self.Ux,axis=self.ha)/self.utau()), 
                symm("asymm",np.mean(self.Uy,axis=self.ha)/self.utau()), 
                symm("symm",np.mean(self.Uz,axis=self.ha)/self.utau()))

    #standard deviation
    def hstd (self):
        return (self.ynodes(),
                np.std(self.Ux,axis=self.ha), 
                np.std(self.Uy,axis=self.ha), 
                np.std(self.Uz,axis=self.ha))

    #standard deviation symmetrised
    def hstd_symm (self):
        return (self.y_nondim(),
                symm("symm",np.std(self.Ux,axis=self.ha)/self.utau()), 
                symm("symm",np.std(self.Uy,axis=self.ha)/self.utau()), 
                symm("symm",np.std(self.Uz,axis=self.ha)/self.utau()))

        
    #correlation of velocties
    #E(U1*U2)-E(U1)*E(U2)
    def hcor (self,vel1,vel2):
        A1 = getattr(self,vel1)
        A2 = getattr(self,vel2)
        return (self.ynodes(),
                np.mean(A1*A2,axis=self.ha)
                -np.mean(A1,axis=self.ha)*np.mean(A2,axis=self.ha))


    #correlation of velocties symmetrised
    #E(U1*U2)-E(U1)*E(U2)
    def hcor_symm (self,vel1,vel2):
        A1 = getattr(self,vel1)
        A2 = getattr(self,vel2)
        symm_dict = {"Ux":1.0,"Uy":-1.0,"Uz":1.0}
        def symm_kind(a,b):
            if a*b==1:
                return "symm"
            else:
                return "asymm"

        return (self.y_nondim(),
                symm(symm_kind(symm_dict[vel1],symm_dict[vel2]), np.mean(A1*A2,axis=self.ha)
                -np.mean(A1,axis=self.ha)*np.mean(A2,axis=self.ha))/(self.utau()**2))

    #kinetic energy
    def hek (self):
        return  (self.ynodes(),
                (np.var(self.Ux,axis=self.ha) +  np.var(self.Uy,axis=self.ha) +  np.var(self.Uz,axis=self.ha))/3.0)


    #kinetic energy symmetrised
    def hek_symm (self):
        return  (self.y_nondim(),
                symm("symm",np.var(self.Ux,axis=self.ha) +  np.var(self.Uy,axis=self.ha) +  np.var(self.Uz,axis=self.ha))/(3.0*(self.utau()**2))) 



