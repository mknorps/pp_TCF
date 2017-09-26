
#
# File name: griddata_test.py
# Created by: gemusia
# Creation date: 24-09-2017
# Last modified: 26-09-2017 11:42:12
# Purpose: 
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


import numpy as np
import matplotlib.pyplot as plt
import itertools as it
from scipy.interpolate import griddata


#DNS grid - we will interpolate LES velocity on it
x = np.linspace(0,4*np.pi,8)
y = np.linspace(-1,1,8)
z = np.linspace(0,2*np.pi,8)

grid_x,grid_y,grid_z = np.meshgrid(x,y,z)


#LES grid
x_LES = np.linspace(0,4*np.pi,3)
y_LES = np.linspace(-1,1,3)
z_LES = np.linspace(0,2*np.pi,3)

print(   len(x_LES),len(y_LES),len(z_LES))

LES_grid_x,LES_grid_y,LES_grid_z = np.meshgrid(x_LES,y_LES,z_LES)
 
points    = zip(LES_grid_x.flatten(), LES_grid_y.flatten(),LES_grid_z.flatten())
value_Ux  = np.arange(27)
value_Uy  = np.ones(27)
value_Uz  = np.ones(27)

print("points: ", len(points), "  \t value_Ux : ",len(value_Ux))

LES_interpolated_Ux = griddata(points,value_Ux,(grid_x,grid_y,grid_z),method ='linear' )
print( LES_interpolated_Ux)
print( len(grid_x),len(grid_y),len(grid_z))

with open('griddata_test.txt','w') as f:
    f.write('# Array shape: {0}\n'.format(LES_interpolated_Ux.shape))

    for LES_slice in LES_interpolated_Ux:
        np.savetxt(f,LES_slice)


with open('griddata_test.txt','r') as f:
    new_data = np.loadtxt('griddata_test.txt')
    print new_data.shape

    new_data = new_data.reshape(8,8,8)

assert np.all(new_data == LES_interpolated_Ux)
