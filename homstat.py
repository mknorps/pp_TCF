# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# File name: homstat.py
# Created by: gemusia
# Creation date: 22-06-2017
# Last modified: 30-06-2017 19:00:33
# Purpose: module for computing statistics of 
#   turbulent channel flow.
#
#   Statistics are taken by averaging over homogeneity directions x and z
#   and are presented a functions of y
#
#   Statistics may be presented in non-dimensional form as functions of y^+
#
#   Input np.array is a nested 3D list of 3-elemnt lists upp - (velocities in x,y,z)
#    upp (0:n,-2:m+2,-2:o+2,3)
#    uppf (0:nf,-1:mf+1,-1:of+1,3)
#    for example:
#   [[[upp(0,-2,-2),...,upp(0,-2,o+2)],...,[upp(0,m+2,-2),...,upp(0,m+2,o+2)]],
#    ...  
#    [[upp(n,-2,-2),...,upp(n,-2,o+2)],...,[upp(n,m+2,-2),...,upp(n,m+2,o+2)]]]
#   
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


import numpy as np

#~~~~~~~~~~ channel ~~~~~~~~~~
#   - class for velocity field of wall bounded turbulent channel flow
#
# Directions of the flow:
# x - streamwise 
# y - wall-normal
# z - spanwise
#
# Mesh size: K x N x M
#
# We assume that stremwise and spanwise directions have uniform mesh
# and that wal-normal direction has nodes in zeroes of Chebyshev poynomials
# y(j)=cos(j*pi/N)



# zeroes of Chebyshev polynomials
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


#class Channel(np.ndarray):
class Channel:

    # homogeneous directions axis
    # note that input data  from files is in format (Uy,Uz,Ux)
    ha=(0,2)

    # for object initiation we need velocities in three directions (Ux,Uy,Uz)
    # it should be given as three numpy arrays of size (K,N,M)
    def __init__(self,Ux,Uy,Uz,K,N,M,Retau=150):

        #TODO - think if it can be coded simpler
        wrongSize = []
        UxShape = np.array(Ux).shape
        UyShape = np.array(Uy).shape
        UzShape = np.array(Uz).shape
        if UxShape <> (K,N,M):
            wrongSize.append(("Ux has wrong shape: ",UxShape))
        if UyShape <> (K,N,M):
            wrongSize.append(("Uy has wrong shape: ",UyShape))
        if UzShape <> (K,N,M):
            wrongSize.append(("Uz has wrong shape: ",UzShape))

        if  wrongSize <> []:
            raise ValueError("Incorrect mesh size; ", wrongSize, 'correct is: ( %d, %d, %d) ' % (K,N,M))
        else: 
            self.Ux    = np.array(Ux)
            self.Uy    = np.array(Uy)
            self.Uz    = np.array(Uz)
            self.K     = K
            self.N     = N
            self.M     = M
            self.Retau = Retau

    def displayU(self):
        print "Ux = ",self.Ux
        print "Uy = ",self.Uy
        print "Uz = ",self.Uz

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

#...............................................................
#      STATISTICS
#...............................................................

    #mean
    def hmean (self):
        return (self.ynodes(),
                np.mean(self.Ux,axis=self.ha), 
                np.mean(self.Uy,axis=self.ha), 
                np.mean(self.Uz,axis=self.ha))

    #mean symmetrised
    def hmean_symm (self):
        return (self.y_nondim(),
                symm("symm",np.mean(self.Ux,axis=self.ha)), 
                symm("asymm",np.mean(self.Uy,axis=self.ha)), 
                symm("symm",np.mean(self.Uz,axis=self.ha)))

    #standard deviation
    def hstd (self):
        return (self.ynodes(),
                np.std(self.Ux,axis=self.ha), 
                np.std(self.Uy,axis=self.ha), 
                np.std(self.Uz,axis=self.ha))

    #standard deviation symmetrised
    def hstd_symm (self):
        return (self.y_nondim(),
                symm("symm",np.std(self.Ux,axis=self.ha)), 
                symm("symm",np.std(self.Uy,axis=self.ha)), 
                symm("symm",np.std(self.Uz,axis=self.ha)))

        
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
                -np.mean(A1,axis=self.ha)*np.mean(A2,axis=self.ha)))

    #kinetic energy
    def hek (self):
        return  (self.ynodes(),
                np.var(self.Ux,axis=ha) +  np.var(self.Uy,axis=ha) +  np.var(self.Uz,axis=ha) )


    #kinetic energy symmetrised
    def hek_symm (self):
        return  (self.y_nondim(),
                symm("symm",np.var(self.Ux,axis=ha) +  np.var(self.Uy,axis=ha) +  np.var(self.Uz,axis=ha))) 



