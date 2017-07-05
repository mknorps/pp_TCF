# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# File name: homfigs.py
# Created by: gemusia
# Creation date: 30-06-2017
# Last modified: 05-07-2017 15:12:17
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
   

    # *args: list of plot features for ax1.plot
    #     (xdata,ydata,str linetype,str label)
    # **kwargs: list of axes features for ax1.set_$(STH)
    #     possible keys:
    #     title,xlabel,ylabel,xlim,ylim,xscale,yscale


    def __init__(self,*args,**kwargs):
        self.args = args
        self.kwargs = kwargs
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111)

        #default values
        self.ax.set_xlabel('$y^{+}$')
        self.ax.set_xlim([0,160])

        #parameters read from the input
        for key, val in self.kwargs.iteritems():
            getattr(self.ax,'set_'+key)(val)

    def save(self,name):
        self.fig.savefig(name)


#ff = Homfig (title='very important title',xlabel="xxx",ylabel="yyy")
#ff.hdraw()
#ff.save()

'''
    ax1.plot(DataFiles[cntr][:,0],DataFiles[cntr][:,1], 'k-^', label='model')
    ax1.plot(DataFiles_LES[cntr][:,0],DataFiles_LES[cntr][:,1], '-', label='LES')
    ax1.plot(DataFiles_DNS[cntr][:,0],DataFiles_DNS[cntr][:,1], '--', label='DNS')

    def hdraw(self):
        for arg in self.args:
            print arg
            # *args - unpacked as positional arguments
            self.ax.plot(arg[0],arg[1],**arg[2])
        #self.ax.plot([0, 1, 2, 3, 4], [0, 1, 2, 3, 4], label="Test", color='g')
        leg = self.ax.legend(loc=4)

'''
