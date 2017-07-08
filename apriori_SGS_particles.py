# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# File name: apriori_SGS_particles.py
# Created by: gemusia
# Creation date: 08-07-2017
# Last modified: 08-07-2017 22:18:28
# Purpose:computation of apriori statistics of particles,
#         deterministic terms of equation are compared
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


import numpy as np
import matplotlib.pyplot as plt
import particlestat as ps
import homfigs as hfig



# declaration of picture attributes

LineStyle = {'DNS':'solid','LES':'dashed'}
coordinates = {1:'x',2:'y',3:'z'}
statistics = {0:("hmean_symm",),1:("hstd_symm",),
        2:("hcor_symm","Ux","Uy"),
        3:("hcor_symm","Uy","Uz"),
        4:("hcor_symm","Ux","Uz"),
        5:("hek_symm",)}
ylabels = {1:"$<U>_{x}$",2:"$<U>_{y}$",
        3:"$<U>_{z}$",4:"$<U_{x},U_{y}>$",
        5:"$<U_{y},U_z>$",6:"$<U_{x}U_{z}$",
        7:"$e_{k}$"}
ptitles = {"hmean_symm":"mean of U","hstd_symm":"Standard deviation of U","hcor_symm":"Correlation of ","hek_symm":"Kinetic energy"}


#data loading

#initialising Particle class instance
#the data is written to file in Fortran program with following code
'''
       do j=1,npar
            write(4,'(30e13.5)')(pos(j,i),i=1,3),(vpar(j,i),i=1,3),
     $     (upar(j,i),i=1,3), (uparf(j,i),i=1,3),
     $     (dupar(j,i,1),i=1,3), (dupar(j,i,2),i=1,3), (dupar(j,i,3),i=1,3),
     $     (duparf(j,i,1),i=1,3), (duparf(j,i,2),i=1,3), (duparf(j,i,3),i=1,3)
        enddo
'''
#the code has different direction naming and we have to permute directions
# fortran code     |    current convention (widely used in presenting channel flow data)
#---------------------------------------------------------------------------
# x - wall-normal  |   x - streamwise
# y - spanwise     |   y - wall-normal
# z - streamwise   |   z - spanwise
#
#we permute (x,y,z)->(z,x,y)
part_2501 = ps.Particles(tf[2],tf[0],tf[1],tf[5],tf[3],tf[4],
        Ux=tf[8],Uy=tf[6],Uz=tf[7],Ufx=tf[11],Ufy=tf[9],Ufz=tf[10],
        dUxdx=tf[20],dUxdy=tf[18],dUxdz=tf[19],
        dUydx=tf[14],dUydy=tf[12],dUydz=tf[13],
        dUzdx=tf[17],dUzdy=tf[15],dUzdz=tf[16],
        dUfxdx=tf[29],dUfxdy=tf[27],dUfxdz=tf[28],
        dUfydx=tf[23],dUfydy=tf[21],dUfydz=tf[22],
        dUfzdx=tf[26],dUfzdy=tf[24],dUfzdz=tf[25])

