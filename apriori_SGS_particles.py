# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# File name: apriori_SGS_particles.py
# Created by: gemusia
# Creation date: 08-07-2017
# Last modified: 09-07-2017 17:11:49
# Purpose:computation of apriori statistics of particles,
#         deterministic terms of equation are compared
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


import numpy as np
import matplotlib.pyplot as plt
import particlestat as ps
import homfigs as hfig
from os.path import expanduser



# declaration of picture attributes

LineStyle = {'fterm':'solid','pterm':'dashed'}
coordinates = {0:'x',1:'y',2:'z'}
terms = {0:"termx",1:"termy",2:"termz"}
ptype = {"St1","St5","St25","fluid"}


statistics = {0:("pmean",),1:("pstd",),
        2:("pcor","Ux","Uy"),
        3:("pcor","Uy","Uz"),
        4:("pcor","Ux","Uz"),
        5:("pek",)}



#data loading and
#initialising Particle class instance
#the data is written to file in Fortran program with following code
'''
       do j=1,npar
            write(4,'(9e12.4)')pos(j,1),pos(j,2),pos(j,3),       FTERM
     $    ((upar(j,1)-uparf(j,1))*duparf(j,i,1)+
     $     (upar(j,2)-uparf(j,2))*duparf(j,i,2)+(upar(j,3)-uparf(j,3))*duparf(j,i,3), i=1,3),
     $     ((vpar(j,1)-upar(j,1))*(dupar(j,i,1)-duparf(j,i,1))+   PTERM
     $      (vpar(j,2)-upar(j,2))*(dupar(j,i,2)-duparf(j,i,2))+
     $      (vpar(j,3)-upar(j,3))*(dupar(j,i,3)-duparf(j,i,3)), i=1,3)

'''
#the code has different direction naming and we have to permute directions
# fortran code     |    current convention (widely used in presenting channel flow data)
#---------------------------------------------------------------------------
# x - wall-normal  |   x - streamwise
# y - spanwise     |   y - wall-normal
# z - streamwise   |   z - spanwise
#
#we permute (x,y,z)->(z,x,y)

path = expanduser("~") + "/wyniki/apriori/fede_terms"

for val in ptype:
    tf=np.transpose(np.loadtxt(path+"/fede_terms_2502_"+val))

    part = ps.Particles(tf[2],tf[0],tf[1],
          ftermx=tf[5],ftermy=tf[3],ftermz=tf[4],
          ptermx=tf[8],ptermy=tf[6],ptermz=tf[7])

    for stattype in ("pmean","pstd"):
        for key,val2 in terms.iteritems():
            stat_fterm = part.stat_symm(stattype,"symm","f"+val2)
            stat_pterm = part.stat_symm(stattype,"symm","p"+val2)

            #figures

            statfig = hfig.Homfig(title=key, ylabel=key)

            statfig.add_plot(*stat_fterm,linestyle=LineStyle['fterm'],label='fterm')
            statfig.add_plot(*stat_pterm,linestyle=LineStyle['pterm'],label='pterm')
            
            statfig.hdraw()
            statfig.save(val2 +"_"+stattype+"_"+val+".eps")

            plt.close(statfig.fig)

'''
tf=np.transpose(np.loadtxt("200_particles"))

part = ps.Particles(tf[2],tf[0],tf[1],
      Vx=tf[5],Vy=tf[3],Vz=tf[4],
      Ux=tf[8],Uy=tf[6],Uz=tf[7],
      Ufx=tf[11],Ufy=tf[9],Ufz=tf[10])

for stattype in ("pmean","pstd"):
    for key,val2 in coordinates.iteritems():
        stat_V = part.stat_symm(stattype,"symm","V"+val2)
        stat_U = part.stat_symm(stattype,"symm","U"+val2)
        stat_Uf = part.stat_symm(stattype,"symm","Uf"+val2)

        #figures

        statfig = hfig.Homfig(title=stattype, ylabel=val2)

        statfig.add_plot(*stat_V,linestyle="solid",label='V')
        statfig.add_plot(*stat_U,linestyle="dashed",label='U')
        statfig.add_plot(*stat_Uf,linestyle="dotted",label='Uf')
        
        statfig.hdraw()
        statfig.save(stattype+"_"+val2+".eps")

        plt.close(statfig.fig)
'''
