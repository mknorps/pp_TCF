# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# File name: apriori_SGS_fede_timestat.py
# Created by: gemusia
# Creation date: 11-07-2017
# Last modified: 12-08-2017 08:48:12
# Purpose:computation of apriori statistics of particles,
#        statistic derived from scratch 
#      - test of possible substitution of (V-U)du^*/dx term 
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

LineStyle = {'fterm':'solid','pterm':'dashed',
        "Ux":'solid',"Ufx":'dashed',"Vx":'dotted',
        "UxUfx":'solid',"UyUfy":'dashed',"UzUfz":'dotted',
        "fluid":'solid',"St1":'dashed',"St5":'dotted',"St25":'dashdot'}
coordinates = {0:'x',1:'y',2:'z'}
terms = {0:"termx",1:"termy",2:"termz"}
ptype = {"St1","St5","St25","fluid"}




#data loading and
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

file_path = expanduser("~") + "/wyniki/apriori/fede_terms_pDNS/"
pict_path = file_path


pfields= pf.ParticleFields(2501,2508,fCoreName=file_path+"SGS_terms_",x=2,y=0,z=1,Vx=5,Vy=3,Vz=4,
     Ux=8,Uy=6,Uz=7,Ufx=11,Ufy=9,Ufz=10, 
     dUxdx=20,dUxdy=14,dUxdz=17, dUydx=18,dUydy=12,dUydz=15, dUzdx=19,dUzdy=13,dUzdz=16, 
     dUfxdx=29,dUfxdy=23,dUfxdz=26, dUfydx=27,dUfydy=21,dUfydz=24, dUfzdx=28,dUfzdy=22,dUfzdz=25)
#pfields= pf.ParticleFields(2501,2508,fCoreName=file_path+"SGS_terms_",x=2,y=0,z=1,Vx=5,Vy=3,Vz=4,
#     Ux=8,Uy=6,Uz=7,Ufx=11,Ufy=9,Ufz=10, 
#     dUxdx=20,dUxdy=18,dUxdz=19, dUydx=14,dUydy=12,dUydz=13, dUzdx=17,dUzdy=15,dUzdz=16, 
#     dUfxdx=29,dUfxdy=27,dUfxdz=28, dUfydx=23,dUfydy=21,dUfydz=22, dUfzdx=26,dUfzdy=24,dUfzdz=25)

def pterm(V1,V2,V3,U1,U2,U3,dUdx,dUdy,dUdz,dUfdx,dUfdy,dUfdz):
    return (V1-U1)*(dUdx-dUfdx) + (V2-U2)*(dUdy-dUfdy) + (V3-U3)*(dUdz-dUfdz)

ptermArglist = [['Vx','Vy','Vz','Ux','Uy','Uz','dUxdx','dUxdy','dUxdz','dUfxdx','dUfxdy','dUfxdz'],
    ['Vx','Vy','Vz','Ux','Uy','Uz','dUydx','dUydy','dUydz','dUfydx','dUfydy','dUfydz'],
    ['Vx','Vy','Vz','Ux','Uy','Uz','dUzdx','dUzdy','dUzdz','dUfzdx','dUfzdy','dUfzdz']]


# separate plots for each stokes number
for StNo in ptype:

    for stattype in ("pmean","pstd"):
      
        # STATISTICS
        pstat = pfields.equationP(StNo,(lambda x:x),stattype,"symm",["Ux"],["Ufx"],["Vx"])  
        sgsStat = pfields.equationP(StNo,(lambda x,y: x-y),stattype,"symm",["Ux","Ufx"],["Uy","Ufy"],["Uz","Ufz"])  


        # FIGURES

        # velocity statistics
        statfig = hfig.Homfig(title="Velocities streamwise", ylabel="U")
        plotFileName = pict_path +stattype +"_"+StNo+".eps"

        for arg in ifilterfalse(lambda x: x=="yplus", set(pstat.keys())): #custom , for this case
            statfig.add_plot(pstat["yplus"],pstat[arg],linestyle=LineStyle[arg],label=arg)
            
        statfig.hdraw()
        statfig.save(plotFileName)
        print "plot created: " + plotFileName
        plt.close(statfig.fig)


        # SGS statistics
        sgsfig = hfig.Homfig(title="SGS velocity", ylabel="$u^*$")
        plotFileNameSGS = pict_path + "Usgs_"+stattype +"_"+StNo+".eps"

        for arg in ifilterfalse(lambda x: x=="yplus", set(sgsStat.keys())): #custom , for this case
            sgsfig.add_plot(sgsStat["yplus"],sgsStat[arg],linestyle=LineStyle[arg],label=arg)
            
        sgsfig.hdraw()
        sgsfig.save(plotFileNameSGS)
        print "plot created: " + plotFileNameSGS
        plt.close(sgsfig.fig)



#several stokes number on one plot
for stattype in ("pmean","pstd"):

    ptermStat = {}

    for StNo in ptype:
            #(V-U)_j*du/dx_j
        ptermStat[StNo] = pfields.equationP(StNo,pterm,stattype,"symm",*ptermArglist)  

    print "ptermStat: ", type(ptermStat['St1']), ptermStat['St1'].keys()
    # pterm statistics : (V-U)_j*du/dx_j
    for direction,ptermKey in zip(range(3),ifilterfalse(lambda x: x=="yplus", set(ptermStat['St1'].keys()))):
        ptermfig = hfig.Homfig(title="pterm ", ylabel="$(V-U)_j*du/dx_j$")
        plotFileNamePterm = pict_path + "pterm_"+stattype +coordinates[direction]+".eps"

        for StNo in ptype: 
            ptermfig.add_plot(ptermStat[StNo]["yplus"],ptermStat[StNo][ptermKey],linestyle=LineStyle[StNo],label=StNo)
            
        ptermfig.hdraw()
        ptermfig.save(plotFileNamePterm)
        print "plot created: " + plotFileNamePterm
        plt.close(ptermfig.fig)
