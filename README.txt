###########################################################
#
#  Documentation of pp_TCF:
#       postprocessing of Turbulent Channel flow
#
#  Created by: mknorps
#  Creation date: 10.07.2017
#
#  Purpose: 
#     - Computation of statistics of fluid and particles
#      from two-phase wall-bounded channel flow
#      with heavy particles.
#      Statistics are averaged over homogeneous directions.
#      There is also a helper module for computing time-space averages.
#     - Generating plots from computed statistics
#     
#
###########################################################

Project consists of modules:

	homfigs.py      -  plotting the data
	homfiles.py     -  loading of fluid data from multiple time points
	homstat.py      -  computing statistics of fluid
	particlestat.py -  computing statistics of particles

and their test suites:

	homfigs_tests.py
	homfiles_tests.py
	homstat_tests.py
	particlestat_tests.py


Examples of usage are:
	apriori_SGS_fluid.py
	apriori_SGS_particles.py


###########################################################
#
#  General description
#
###########################################################


#~~~~~~~~~~~~~~~~~  INPUT  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

The input may be:

(1)	- a text file of fluid field (velocity or other)
	  on a mesh. 
	- mesh nodes number (N,M,O)  in 
	  streamwise, wall-normal and spanwise directions
	  respectively.
	  The mesh of channel flow has uniform grid points 
	  distribution in periodic directions and Chebyshev 
	  polynomial zeros in wall-normal direction 
	  y(j)=cos(j*pi/N)

(2)     - a text file of particle related quantities
