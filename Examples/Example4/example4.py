'''

	Example 4 on using WindFLO library

  Purpose:
  	
  	This example demonstrate how to use the Python WindFLO API to perform 
  	wind farm layout optimization in parallel
  
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


import os
import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Pool
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
MaxSimulRun = 10		# maximum allowable parallel runs (depends on core availability)

windFLO = WindFLO(inputFile = inputFile, nTurbines = nTurbines, libDir = libPath, 
		  turbineFile = turbineFile, terrainfile = terrainfile, runDir = '')
windFLO.terrainmodel = 'IDW'	# change the default terrain model from RBF to IDW



# Function to evaluate a single farm's performance
def EvaluateFarm(x, runDir = ''):

#	windFLO = WindFLO(inputFile = inputFile, nTurbines = nTurbines, libDir = libPath, 
#		  turbineFile = turbineFile, terrainfile = terrainfile, runDir = runDir)

	x = np.array(x).flatten()
	
	k = 0
	for i in range(0, windFLO.nTurbines):
		for j in range(0, 2):
	   		# unroll the variable vector 'x' and assign it to turbine positions		
			windFLO.turbines[i].position[j] = x[k]
			k = k + 1
			
	# Run WindFLO analysis in a particular director runDir
	windFLO.run(clean = True, runDir = runDir)	
	
	# Return the farm power or any other farm output
	return -windFLO.farmPower

# Function called in parallel to run EvaluateFarm
def EvaluateFarm_Helper(x):
	
	# Unzip the name of the directory and the location of turbines
	dir, x = zip(*x)
	x = np.array(x).flatten()

	# Make the directory to run this particular configuration
	os.system('mkdir '+dir[0])
	
	# Evaluate this particular wind farm configuration
	power = EvaluateFarm(x,runDir = dir[0])
	
	# Clean up and remove the created directory
	os.system('rm -rf '+dir[0])

	# Return the power or any other farm output
	return power


# Compute the minimum clearance between turbines in farm must be greater 
# than turbine diameter. The constraint to satisfy is g(x) >= 0
def ComputeClearance(x):

	position = np.zeros((nTurbines,2))
	k = 0
	for i in range(0, nTurbines):
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
	return ( minClearence - diameter ) # g(x) = minClearence - diameter


# Given a set of configurations (i.e. a matrix of turbines x), compute the
# performance of each farm configuration in parallel 
def RunParallel(x):

	# Number of individual farm configurations
	n = x.shape[0]

	# Create pool of worker process... taken as minimum
	# of number of allowable parallel process and number of configurations
	casePool = Pool(min(n,MaxSimulRun))
	
	# Run the cases in parallel... first argument is the function to unzip arguments
	# prep the directories and call the EvaluateFarm function... the second argument
	# zips the name of the directory and the turbine locations in each configuration
	outputs = casePool.map( EvaluateFarm_Helper , zip( ('Run'+str(i)+'/', x[i,:]) for i in range(0,n) ) )
	casePool.close()

	return np.array(outputs)


if __name__ == "__main__":

	# Two variable per turbines (its x and y coordinates)
	lbound = np.zeros(nTurbines*2)	#lower bounds of x and y of turbines
	ubound = np.ones(nTurbines*2)*2000	#upper bounds of x and y of turbines


	# Solve the optimization problem in parallel
	xBest, bestPower = pso(lbound, ubound, parallelFunc = RunParallel, ieqcons=[ComputeClearance],
			 args=(), kwargs={}, swarmsize=50, omega=0.5, phip=0.5, phig=0.5, maxiter=100, 
			 minstep=1e-8, minfunc=1e-8, debug=True)


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

	

	

