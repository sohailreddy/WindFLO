'''

	Example 3 on using WindFLO library

  Purpose:
  	
  	This example demonstrate how to use the Python WindFLO API to compute several
  	wind farm layouts in serial and in parallel.
  
  Licensing:
  
	This code is distributed under the Apache License 2.0 
	
  Author:
  	Sohail R. Reddy
  	sredd001@fiu.edu
  	
'''


import sys
# Append the path to the API directory where WindFLO.py is
sys.path.append('../../API/')

import os
import numpy as np
import matplotlib.pyplot as plt
import datetime

from multiprocessing import Pool
from WindFLO import WindFLO

###############################################################################
#	WindFLO Settings and Params
nTurbines = 25			# Number of turbines
libPath = '../../release/'		# Path to the library
inputFile = 'WindFLO.dat'	# Input file to read
turbineFile = 'V90-3MW.dat'	# Turbine parameters
terrainfile = 'terrain.dat'	# Terrain file
MaxSimulRun = 10		# maximum allowable parallel runs (depends on core availability)

windFLO = WindFLO(inputFile = inputFile, nTurbines = nTurbines, libDir = libPath, 
	  turbineFile = turbineFile, terrainfile = terrainfile, runDir = '')


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


# Runs each configuration in serial where 'x' is now a matrix
def RunSerial(x):

	# Number of individual farm configurations
	n = x.shape[0]
	
	# Outputs for each farm	
	outputs = np.zeros(n)	
	
	# Start timing
	start = datetime.datetime.now()
	
	for i in range(0,n):
		print('Running config: ',i)
		# Evaluate one farm configuration at a time
		outputs[i] = EvaluateFarm(x[i,:], runDir = 'Serial/')

	# Stop timing
	end = datetime.datetime.now()
	print('Time for serial: ', end-start)
	
	# Save outputs
	np.savetxt('SerialOutputs.dat',outputs)


# Runs all configurations in parallel where 'x' is now a matrix
def RunParallel(x):

	# Number of individual farm configurations
	n = x.shape[0]

	# Start timing	
	start = datetime.datetime.now()
	
	# Create pool of worker process... taken as minimum
	# of number of allowable parallel process and number of configurations
	casePool = Pool(min(n,MaxSimulRun))
	
	# Run the cases in parallel... first argument is the function to unzip arguments
	# prep the directories and call the EvaluateFarm function... the second argument
	# zips the name of the directory and the turbine locations in each configuration
	outputs = casePool.map( EvaluateFarm_Helper , zip( ('Run'+str(i)+'/', x[i,:]) for i in range(0,n) ) )
	casePool.close()

	# Stop timing
	end = datetime.datetime.now()
	print('Time for parallel: ', end-start)

	# Save outputs
	np.savetxt('ParallelOutputs.dat',outputs)



if __name__ == "__main__":

	# Total number of configurations to run
	nRuns = 300

	
	lb = 0.0	# Lower bound of the farm
	ub = 2000.0	# Upper bound of the farm

	# Generate matrix of configurations
	x = np.zeros((nRuns, nTurbines*2))
	for i in range(0,nRuns):
		x[i,:] = lb + (ub-lb)*np.random.rand(nTurbines*2) 

	# Runs all configurations in serial
	RunSerial(x)
	
	# Runs all configurations in serial	
	RunParallel(x)

	

