'''

	Example 1 on using WindFLO library

  Purpose:
  	
  	This example demonstrate how to use the Python WindFLO API for serial optimization 
  	of wind farm layout. It uses WindFLO to analyze each layout configuration
  	and a the particle swarm optimization (PSO) algorithm from the pyswarm package
  	to perform the optimization.  The layout is optimized for maximum power generation
  	and incorporates constraints on the minimum allowable clearance between turbines. 
  
  	IMPORTANT: The PSO minimizes the fitness function, therefore the negative of the
  			   farm power is returned by the EvaluateFarm function
  
  Licensing:
  
	This code is distributed under the Apache License 2.0 
	
  Author:
  	Sohail R. Reddy
  	sredd001@fiu.edu
  	
'''


import sys
# Append the path to the API directory where WindFLO.py is
sys.path.append('../../API/')
# Append the path to the Optimizers directory where pso.py is
sys.path.append('../../Optimizers/')


import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg

from pso import pso
from WindFLO import WindFLO

###############################################################################
#	WindFLO Settings and Params
nTurbines = 25			# Number of turbines
libPath = '../../release/'		# Path to the shared library libWindFLO.so
inputFile = 'WindFLO.dat'	# Input file to read
turbineFile = 'V90-3MW.dat'	# Turbine parameters
terrainfile = 'terrain.dat'	# Terrain file
diameter = 90.0			# Diameter to compute clearance constraint

windFLO = WindFLO(inputFile = inputFile, nTurbines = nTurbines, libDir = libPath, 
		  turbineFile = turbineFile,
		  terrainfile = terrainfile)
windFLO.terrainmodel = 'IDW'	# change the default terrain model from RBF to IDW




# Function to evaluate the farm's performance
def EvaluateFarm(x):

	k = 0
	for i in range(0, nTurbines):
		for j in range(0, 2):
	   		# unroll the variable vector 'x' and assign it to turbine positions
			windFLO.turbines[i].position[j] = x[k]
			k = k + 1
	# Run WindFLO analysis
	windFLO.run(clean = True)	
	
	# Return the farm power or any other farm output
  	# NOTE: The negative value is returns since PSO minimizes the fitness value
	return -windFLO.farmPower



# Compute the minimum clearance between turbines in farm must be greater 
# than turbine diameter. The constraint to satisfy is g(x) >= 0
def ComputeClearance(x):

	position = np.zeros((nTurbines,2))
	k = 0
	for i in range(0, windFLO.nTurbines):
		for j in range(0, 2):
	   		# unroll the variable vector 'x' and assign it to turbine positions		
			position[i,j] = x[k]
			k = k + 1

	minClearence = 10000000
	for i in range(0, nTurbines):
		for j in range(i+1, nTurbines):
			# get the minimum clearance between turbines in the farm
			minClearence = min(minClearence, linalg.norm( position[i,0:2] - position[j,0:2] ))
			
	# if minClearence < diameter, then g(x) < 0 and the constraint is violated
	# if minClearence >= 0 then g(x) >= 0 and the constraint is satisfied			
	return (minClearence - diameter ) # g(x) = minClearence - diameter



# Main function
if __name__ == "__main__":

	# Two variable per turbines (its x and y coordinates)
	lbound = np.zeros(nTurbines*2)	#lower bounds of x and y of turbines
	ubound = np.ones(nTurbines*2)*2000	#upper bounds of x and y of turbines

	# Solve the optimization problem 
	xBest, bestPower = pso(lbound, ubound, func = EvaluateFarm, ieqcons=[ComputeClearance], args=(), kwargs={},
			swarmsize=50, omega=0.5, phip=0.5, phig=0.5, maxiter=10, minstep=1e-8,
			minfunc=1e-8, debug=True)


	# Evaluate the best farm layout
	EvaluateFarm(xBest)
	windFLO.run(clean = True, resFile = 'WindFLO.res')
	
	# Plot the optimum configuration	
	fig = plt.figure(figsize=(8,5), edgecolor = 'gray', linewidth = 2)
	ax = windFLO.plotWindFLO2D(fig, plotVariable = 'P', scale = 1.0e-3, title = 'P [kW]')
	windFLO.annotatePlot(ax)
	plt.show()
	
	# save the optimum to a file
	np.savetxt('optimum.dat', xBest)

