# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# File name: homstat.py
# Created by: gemusia
# Creation date: 22-06-2017
# Last modified: 25-06-2017 14:35:43
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

class Channel:

    # for object initiation we need velocities in three directions (Ux,Uy,Uz)
    # it should be given as three numpy arrays of size (K,N,M)
    def __init__(self,Ux,Uy,Uz,K,N,M):
        wrongSize = []
        UxShape = Ux.shape()
        UyShape = Uy.shape()
        UzShape = Uz.shape()
        if UxShape <> (K,N,M):
            wrongSize.append(("Ux has wrong shape: ",UxShape))
        if UyShape <> (K,N,M):
            wrongSize.append(("Uy has wrong shape: ",UyShape))
        if UzShape <> (K,N,M):
            wrongSize.append(("Uz has wrong shape: ",UzShape))

        if  wrongSize <> []:
            raise ValueError("Incorrect mesh size; ", wrongSize, "correct is: ( %d, %d, %d) ") % (K,N,M)
        else: 
            self.Ux = Ux
            self.Uy = Uy
            self.Uz = Uz

    def displayU(self):
        print "Ux = ",Ux
        print "Uy = ",Uy
        print "Uz = ",Uz







#axis for computing means over homogeneity directions
ha = (1,2)


#mean
def hmean (lst):
    return np.mean(lst,axis=ha)

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


