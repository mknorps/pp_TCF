# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# File name: homfigs.py
# Created by: gemusia
# Creation date: 30-06-2017
# Last modified: 18-07-2017 11:59:53
# Purpose: module for creating matplotlib figures of
#          statistics created with module Channel
#          from 'homstat.py'
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


import numpy as np
import matplotlib.pyplot as plt



# class of figures for channel flow
# with one subfigure

class Homfig:
   

    # **kwargs_axes: list of axes features for ax1.set_$(STH)
    #     possible keys:
    #     title,xlabel,ylabel,xlim,ylim,xscale,yscale


    def __init__(self,**kwargs_axes):
        self.kwargs_axes = kwargs_axes
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111)

        self.plt_data = []

        #default values
        self.ax.set_xlabel('$y^{+}$')
        self.ax.set_xlim([0,160])

        #parameters read from the input
        for key, val in self.kwargs_axes.iteritems():
            getattr(self.ax,'set_'+key)(val)

    def add_plot(self,*plt_args,**plt_kwargs):
        self.plt_data.append((plt_args,plt_kwargs))

    def hdraw(self,leg_loc=0):
        for args,kwargs in self.plt_data:
            # *args - unpacked as positional arguments
            self.ax.plot(*args,**kwargs)
            leg = self.ax.legend(loc=leg_loc)

    def save(self,name):
        self.fig.savefig(name)



'''
    ax1.plot(DataFiles[cntr][:,0],DataFiles[cntr][:,1], 'k-^', label='model')
    ax1.plot(DataFiles_LES[cntr][:,0],DataFiles_LES[cntr][:,1], '-', label='LES')
    ax1.plot(DataFiles_DNS[cntr][:,0],DataFiles_DNS[cntr][:,1], '--', label='DNS')


'''
