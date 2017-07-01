# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# File name: homfigs.py
# Created by: gemusia
# Creation date: 30-06-2017
# Last modified: 01-07-2017 21:03:06
# Purpose: module for creating matplotlib figures of
#          statistics created with module Channel
#          from 'homstat.py'
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


import numpy as np
import matplotlib.pyplot as plt



class Homfig(Figure):
    
    #tuple of possible subplot names and configurations
    subplot_dict = [{ax1:111}, {ax1:121,ax2:122},{ax1:221,ax2:222,ax3:223},{ax1:221,ax2:222,ax3:223,ax4:224}]

    def __init__(self,n):
        self=plt.figure()
        for x in range(n):
            self. = self.add_subplot(value)

    def get_axes_params(self):
        return self.axes_params

    def set_axes_params(self, **kwargs)
        if kwargs is not None:
            for key,value in kwargs.iteritems():
                self.
        title,xlabel,ylabel,xlim,xscale,yscale):
'''
    fig2=plt.figure()
    ax1 = fig2.add_subplot(111)  
    ax1.set_title('St'+ x)
    ax1.set_xlabel('$y^{+}$')
    ax1.set_ylabel('$concentration$')
    ax1.set_xlim([0,160])
    ax1.set_xscale('log')
    ax1.set_yscale('log')

    ax1.plot(DataFiles[cntr][:,0],DataFiles[cntr][:,1], 'k-^', label='model')
    ax1.plot(DataFiles_LES[cntr][:,0],DataFiles_LES[cntr][:,1], '-', label='LES')
    ax1.plot(DataFiles_DNS[cntr][:,0],DataFiles_DNS[cntr][:,1], '--', label='DNS')


    leg = ax1.legend(loc=4)
    plt.savefig(directory + 'Conc_St'+x+'.eps')
    plt.close(fig2)
'''
