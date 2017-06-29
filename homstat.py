# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# File name: homstat.py
# Created by: gemusia
# Creation date: 22-06-2017
# Last modified: 28-06-2017 18:46:52
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

#class Channel(np.ndarray):
class Channel:

    # homogeneous directions axis
    # note that input data  from files is in format (Uy,Uz,Ux)
    ha=(1,2)
    # for object initiation we need velocities in three directions (Ux,Uy,Uz)
    # it should be given as three numpy arrays of size (K,N,M)
    def __init__(self,Ux,Uy,Uz,K,N,M):
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
            self.Ux = np.array(Ux)
            self.Uy = np.array(Uy)
            self.Uz = np.array(Uz)
            self.shape = (K,N,M)

    def displayU(self):
        print "Ux = ",self.Ux
        print "Uy = ",self.Uy
        print "Uz = ",self.Uz

    # zeroes of Chebyshev polynomials
    def ChebZeros(N):
        def ChebCurr(y):
            return np.cos(y*np.pi/N)
        return ChebCurr

    # nodes are taken as zeroes of Chebyshev pomynomials
    def ynodes(self):
        return map(ChebZeros(self.shape[1]),range(self.shape[1]+1))

    #mean
    def hmean (self):
        return (self.ynodes,np.mean(self.Ux,axis=ha), np.mean(self.Uy,axis=ha), np.mean(self.Uz,axis=ha))

'''
    #standard deviation
    def hstd (lst):
        return np.std(lst,axis=ha)

    #correlation of velocties
    #E(U1*U2)-E(U1)*E(U2)
    def hcor (lst1,lst2):
        return hmean(lst1*lst2)-hmean(lst1)*hmean(lst2)


    #kinetic energy
    def hek (lstUx,lstUy,lstUz):
        return np.var(lstUx,axis=ha) +  np.var(lstUy,axis=ha) +  np.var(lstUz,axis=ha) 



'''

