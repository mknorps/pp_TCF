# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# File name: apriori_SGS_particles.py
# Created by: gemusia
# Creation date: 08-07-2017
# Last modified: 11-07-2017 13:58:47
# Purpose:computation of apriori statistics of particles,
#         deterministic terms of equation are compared
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


import numpy as np
import matplotlib.pyplot as plt
import particlestat as ps
import pfiles as pf
import homfigs as hfig
from os.path import expanduser
from itertools import ifilterfalse


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

file_path = expanduser("~") + "/wyniki/apriori/fede_terms_pDNS/"
pict_path = file_path

pfields=pf.ParticleFields(2501,2519,fCoreName=file_path+"fede_terms_",x=2,y=0,z=1,ftermx=5,ftermy=3,ftermz=4, ptermx=8,ptermy=6,ptermz=7)


for StNo in ptype:

    statArgList = [] #list of statistics that will be computed

    for stattype in ("pmean","pstd"):
        for key,val2 in terms.iteritems():
            statArgList.append([stattype,"symm",["f"+val2]])
            statArgList.append([stattype,"symm",["p"+val2]])
   
   # print statArgList

    pstat = pfields.statsP(StNo,*statArgList) #computation of all required statistic
                                            #to improve efficiency and not opening
                                            #big data files too many times
    print map(lambda x: x[1:],pstat.keys())

            #figures
    for arg in ifilterfalse(lambda x: x=="plus", set(map(lambda y: y[1:],pstat.keys()))): #custom , for this case

        statfig = hfig.Homfig(title=arg, ylabel=arg)
        print arg

        #TODO
        statfig.add_plot(pstat["yplus"],pstat["f"+arg],linestyle=LineStyle['fterm'],label='fterm')
        statfig.add_plot(pstat["yplus"],pstat["p"+arg],linestyle=LineStyle['pterm'],label='pterm')
        
        statfig.hdraw()
        plotFileName = pict_path + arg +"_"+StNo+".eps"
        statfig.save(plotFileName)
        print "plot created: " + plotFileName

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
