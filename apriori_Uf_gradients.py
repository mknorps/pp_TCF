# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# File name: apriori_tests_timestat.py
# Created by: gemusia
# Creation date: 21-07-2017
# Last modified: 19-08-2017 08:57:33
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
from itertools import ifilter


# declaration of picture attributes

LineStyle = {'fterm':'solid','pterm':'dashed',
        "Ux":'solid',"Ufx":'dashed',"Vx":'dotted',
        "UxUfx":'solid',"UyUfy":'dashed',"UzUfz":'dotted',
        "fluid":'solid',"St1":'dashed',"St5":'dotted',"St25":'dashdot'}
coordinates = {0:'x',1:'y',2:'z'}
terms = {0:"termx",1:"termy",2:"termz"}
ptype = {"St1","St5","St25","fluid"}


# normalisation constants
utau = 0.0429
Retau = 150
termplus = utau * Retau #parameter for non-dimensialisation

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

#DATA STRUCTURE
pfields= pf.ParticleFields(2501,2508,fCoreName=file_path+"SGS_terms_",x=2,y=0,z=1,Vx=5,Vy=3,Vz=4,
     Ux=8,Uy=6,Uz=7,Ufx=11,Ufy=9,Ufz=10, 
     dUxdx=20,dUxdy=14,dUxdz=17, dUydx=18,dUydy=12,dUydz=15, dUzdx=19,dUzdy=13,dUzdz=16, 
     dUfxdx=29,dUfxdy=23,dUfxdz=26, dUfydx=27,dUfydy=21,dUfydz=24, dUfzdx=28,dUfzdy=22,dUfzdz=25)
#pfields= pf.ParticleFields(2501,2508,fCoreName=file_path+"SGS_terms_",x=2,y=0,z=1,Vx=5,Vy=3,Vz=4,
#     Ux=8,Uy=6,Uz=7,Ufx=11,Ufy=9,Ufz=10, 
#     dUxdx=20,dUxdy=18,dUxdz=19, dUydx=14,dUydy=12,dUydz=13, dUzdx=17,dUzdy=15,dUzdz=16, 
#     dUfxdx=29,dUfxdy=27,dUfxdz=28, dUfydx=23,dUfydy=21,dUfydz=22, dUfzdx=26,dUfzdy=24,dUfzdz=25)

#functions for testing:
#$(V-U)_j*du/dx_j$  and  $(V-Uf)_j*du/dx_j$
def pterm(V1,V2,V3,U1,U2,U3,dUdx,dUdy,dUdz,dUfdx,dUfdy,dUfdz):
    return (V1-U1)*(dUdx-dUfdx) + (V2-U2)*(dUdy-dUfdy) + (V3-U3)*(dUdz-dUfdz)

#$(V-Uf)_j*dUf/dx_j$
def pterm_test2(V1,V2,V3,U1,U2,U3,dUdx,dUdy,dUdz):
    return (V1-U1)*dUdx + (V2-U2)*dUdy + (V3-U3)*dUdz

#$(V-U)_j*du/dx_j$
ptermArglist = [['Vx','Vy','Vz','Ux','Uy','Uz','dUxdx','dUxdy','dUxdz','dUfxdx','dUfxdy','dUfxdz'],
    ['Vx','Vy','Vz','Ux','Uy','Uz','dUydx','dUydy','dUydz','dUfydx','dUfydy','dUfydz'],
    ['Vx','Vy','Vz','Ux','Uy','Uz','dUzdx','dUzdy','dUzdz','dUfzdx','dUfzdy','dUfzdz']]

#$(V-Uf)_j*du/dx_j$
ptermArglist_test1 = [['Vx','Vy','Vz','Ufx','Ufy','Ufz','dUxdx','dUxdy','dUxdz','dUfxdx','dUfxdy','dUfxdz'],
    ['Vx','Vy','Vz','Ufx','Ufy','Ufz','dUydx','dUydy','dUydz','dUfydx','dUfydy','dUfydz'],
    ['Vx','Vy','Vz','Ufx','Ufy','Ufz','dUzdx','dUzdy','dUzdz','dUfzdx','dUfzdy','dUfzdz']]

#$(V-Uf)_j*dUf/dx_j$
ptermArglist_test2 = [['Vx','Vy','Vz','Ufx','Ufy','Ufz','dUfxdx','dUfxdy','dUfxdz'],
    ['Vx','Vy','Vz','Ufx','Ufy','Ufz','dUfydx','dUfydy','dUfydz'],
    ['Vx','Vy','Vz','Ufx','Ufy','Ufz','dUfzdx','dUfzdy','dUfzdz']]

#$(V-U)_j*dU/dx_j$
ptermArglist_test3 = [['Vx','Vy','Vz','Ux','Uy','Uz','dUxdx','dUxdy','dUxdz'],
    ['Vx','Vy','Vz','Ux','Uy','Uz','dUydx','dUydy','dUydz'],
    ['Vx','Vy','Vz','Ux','Uy','Uz','dUzdx','dUzdy','dUzdz']]

#$(V-U)_j*dUf/dx_j$
ptermArglist_test4 = [['Vx','Vy','Vz','Ux','Uy','Uz','dUfxdx','dUfxdy','dUfxdz'],
    ['Vx','Vy','Vz','Ux','Uy','Uz','dUfydx','dUfydy','dUfydz'],
    ['Vx','Vy','Vz','Ux','Uy','Uz','dUfzdx','dUfzdy','dUfzdz']]

gradients = {}
gradients['dx'] = [['dUfxdx'],['dUfydx'],['dUfzdx']]
gradients['dy'] = [['dUfxdy'],['dUfydy'],['dUfzdy']]
gradients['dz'] = [['dUfxdz'],['dUfydz'],['dUfzdz']]

stype={'pmean':'asymm','pstd':'symm'}
labels = {'dUfxdx':'$d\overline{U}_x/dx$',
          'dUfydx':'$d\overline{U}_y/dx$',
          'dUfzdx':'$d\overline{U}_z/dx$',
          'dUfxdy':'$d\overline{U}_x/dy$',
          'dUfydy':'$d\overline{U}_y/dy$',
          'dUfzdy':'$d\overline{U}_z/dy$',
          'dUfxdz':'$d\overline{U}_x/dz$',
          'dUfydz':'$d\overline{U}_y/dz$',
          'dUfzdz':'$d\overline{U}_z/dz$'
        }
# separate plots for each stokes number
for StNo in ptype:

    for stattype in ("pmean","pstd"):

        def keys_no_yplus(arg):   #keys except "yplus"
            #return ifilterfalse(lambda x: x=="yplus", set(arg))
            arg_sorted = sorted(arg)
            return ifilter(lambda x: x<>"yplus", arg_sorted)
        

        # STATISTICS
#        pstat = pfields.equationP(StNo,pterm,stattype,"none",*ptermArglist)  
#        pstat_test1 = pfields.equationP(StNo,pterm,stattype,"none",*ptermArglist_test1)  
#        pstat_test2 = pfields.equationP(StNo,pterm_test2,stattype,"none",*ptermArglist_test2)  
#        pstat_test3 = pfields.equationP(StNo,pterm_test2,stattype,"none",*ptermArglist_test3)  
#        pstat_test4 = pfields.equationP(StNo,pterm_test2,stattype,"none",*ptermArglist_test4)  

        # FIGURES




        pstat_gradient = {}
        gradfig = hfig.Homfig(title="gradients  of $\overline{U}$" , ylabel="$d\overline{U}/dx_j$", xscale='log',xlim=[1,160])
        plotFileNamePterm = pict_path + "gradients_"+stattype +"_"+StNo+".eps"

        for component in ['dx','dz']:
            pstat_gradient[component] = pfields.equationP(StNo,lambda x : x,stattype,'symm',*gradients[component])  
    
            iterable1 =  keys_no_yplus(pstat_gradient[component].keys())
            for pKey in iterable1:
                gradfig.add_plot(pstat_gradient[component]["yplus"],pstat_gradient[component][pKey]/termplus,linestyle='dotted',label=labels[pKey[1:]])
        
        pstat_gradient['dy'] = pfields.equationP(StNo,lambda x : x,stattype,stype[stattype],*gradients['dy'])  
       

        iterable1 =  keys_no_yplus(pstat_gradient['dy'].keys())
        for pKey in iterable1:
            gradfig.add_plot(pstat_gradient['dy']["yplus"],pstat_gradient['dy'][pKey]/termplus,label=labels[pKey[1:]])

        gradfig.hdraw()
        gradfig.save(plotFileNamePterm)
        print "plot created: " + plotFileNamePterm
        plt.close(gradfig.fig)

        '''
        # velocity statistics
        iterable =  zip(range(3),keys_no_yplus(pstat.keys()),
                keys_no_yplus(pstat_test1.keys()),
                keys_no_yplus(pstat_test2.keys()),
                keys_no_yplus(pstat_test3.keys()),keys_no_yplus(pstat_test4.keys()))

        print iterable

        for direction,pKey,pKey_test1,pKey_test2,pKey_test3,pKey_test4 in iterable:
           # ptermfig = hfig.Homfig(title="pterm ", ylabel="$(V-U)_j*du/dx_j$",xlim=[-1,1])
           # plotFileNamePterm = pict_path + "test_pterm_nosymm_"+stattype +coordinates[direction]+"_"+StNo+".eps"

            #ptermfig.add_plot(pstat["yplus"],pstat[pKey],linestyle='solid',label='excact term, $(V-U)_j*du/dx_j$')
            #ptermfig.add_plot(pstat_test1["yplus"],pstat_test1[pKey_test1],linestyle='dotted',label='$(V-Uf)_j*du/dx_j$')
            #ptermfig.add_plot(pstat_test2["yplus"],pstat_test2[pKey_test2],linestyle='dashed',label='$(V-Uf)_j*dUf/dx_j$')
                
            #ptermfig.hdraw()
            #ptermfig.save(plotFileNamePterm)
            #print "plot created: " + plotFileNamePterm
            #plt.close(ptermfig.fig)


            ptermfullvel = hfig.Homfig(title="pterm ", ylabel="$(V-U)_j*dU/dx_j$",xlim=[-1,1])
            plotFileNamePterm = pict_path + "test_VUdUdx_nosymm"+stattype +coordinates[direction]+"_"+StNo+".eps"

            ptermfullvel.add_plot(pstat_test2["yplus"],pstat_test2[pKey_test2],linestyle='solid',label='$(V-Uf)_j*dUf/dx_j$')
            ptermfullvel.add_plot(pstat_test3["yplus"],pstat_test3[pKey_test3],linestyle='dotted',label='$(V-U)_j*dU/dx_j$')
            ptermfullvel.add_plot(pstat_test4["yplus"],pstat_test4[pKey_test4],linestyle='dashed',label='$(V-U)_j*dUf/dx_j$')
                
            ptermfullvel.hdraw()
            ptermfullvel.save(plotFileNamePterm)
            print "plot created: " + plotFileNamePterm
            plt.close(ptermfullvel.fig)
        '''
