# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# File name: homstat.py
# Created by: gemusia
# Creation date: 22-06-2017
# Last modified: 23-06-2017 16:05:18
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


