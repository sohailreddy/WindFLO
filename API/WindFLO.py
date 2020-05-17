#!/usr/bin/python3

from ctypes import *
import numpy as np
import f90nml
import tempfile
import os

import matplotlib.pyplot as plt
from matplotlib  import cm
from numpy import linalg as LA
import mpl_toolkits.mplot3d.art3d as art3d
import matplotlib.ticker as ticker


class TurbineData:

	def __init__(self, **kwargs):
	
		params = kwargs.get('params', np.zeros(30))	
		self.turbineNum = params[0]

		self.position = np.zeros(3)
		self.orientation = np.zeros(3)
		self.velocity = np.zeros(3)
		
		self.position[0] = params[1]
		self.position[1] = params[2]
		self.position[2] = params[3]

		self.orientation[0] = params[4]
		self.orientation[1] = params[5]
		self.orientation[2] = params[6]
		
		self.velocity[0] = params[7]
		self.velocity[1] = params[8]
		self.velocity[2] = params[9]
		
		self.height = params[10]
		self.radius = params[11]
		self.diameter = 2.0 * self.radius
		self.area = params[12]
		self.ratedPower = params[13] 
		self.power = params[14] 
		self.fictitious = bool(params[15])
		self.yaw = True
		self.cpCurve = np.zeros((10,2))

		self.variables = {'position': self.position, 'V': self.velocity, 'H': self.height, 'R': self.radius, \
		'D': self.diameter, 'A': self.area, 'P_r': self.ratedPower, 'P': self.power}

	def UpdateDict(self):
		self.variables = {'position': self.position, 'V': self.velocity, 'H': self.height, 'R': self.radius, \
		'D': self.diameter, 'A': self.area, 'P_r': self.ratedPower, 'P': self.power}

	def WriteInputFile(self, filename):


		tmpFile = tempfile.NamedTemporaryFile(mode = 'w+', delete = True)
		try:
			tmpFile.write('&turbine_data \n /')
			tmpFile.seek(0)			
			namelist = f90nml.read(tmpFile.name)			
		finally:
			tmpFile.close()
						
		namelist['turbine_data']['turbinenum'] = self.turbineNum
		namelist['turbine_data']['position'] = self.position.tolist()
		namelist['turbine_data']['orientation'] = self.orientation.tolist()
		namelist['turbine_data']['yaw'] = self.yaw
		namelist['turbine_data']['fictitious'] = self.fictitious
		namelist['turbine_data']['height'] = self.height
		namelist['turbine_data']['diameter'] = self.radius * 2.0
		namelist['turbine_data']['cpcurve'] = self.cpCurve.tolist()		
		namelist['turbine_data']['ratedpower'] = self.ratedPower
		
		namelist.write(filename, force=True)
		
		
	def SetFromNamelist(self, filename):
		
		namelist = f90nml.read(filename)			
	
		if 'turbinenum' in namelist['turbine_data']:					
			self.turbineNum = namelist['turbine_data']['turbinenum']
			
		if 'position' in namelist['turbine_data']:			
			self.position = np.array(namelist['turbine_data']['position'])
			
		if 'orientation' in namelist['turbine_data']:			
			self.orientation = np.array(namelist['turbine_data']['orientation'])
			
		if 'yaw' in namelist['turbine_data']:			
			self.yaw = namelist['turbine_data']['yaw']

		if 'fictitious' in namelist['turbine_data']:			
			self.fictitious = namelist['turbine_data']['fictitious']
			
		if 'height' in namelist['turbine_data']:			
			self.height = namelist['turbine_data']['height']
			
		if 'diameter' in namelist['turbine_data']:		
			self.diameter	= namelist['turbine_data']['diameter']
		self.radius = self.diameter / 2.0
		self.area = np.pi * self.radius * self.radius
		
		if 'cpcurve' in namelist['turbine_data']:		
			self.cpCurve = np.array(namelist['turbine_data']['cpcurve'])
		
		if 'ratedpower' in namelist['turbine_data']:
			self.ratedPower	= namelist['turbine_data']['ratedpower']
			
		
def fmt(x, pos):
	return str('{0:.2f}'.format(x))

class WindFLO:

	def __init__(self,  **kwargs):
	
		self.runDir = kwargs.get('runDir', '')	
	
		self.rho = kwargs.get('rho', 1.2)
		self.turbulenceintensity = kwargs.get('turbulenceIntensity', 0.0)
		self.windmodel = kwargs.get('windModel', 'constant')
		self.modelvelocity = kwargs.get('modelVelocity', np.zeros(3))
		self.surfaceroughness = kwargs.get('surfaceRoughness', 0.0)
		self.referenceheight = kwargs.get('surfaceRoughness', 1.0)
		self.wakemodel = kwargs.get('wakeModel', 'Frandsen')
		self.wakemergemodel = kwargs.get('wakeMergeModel', 'Quadratic')
		self.wakeexpansioncoeff = kwargs.get('wakeExpansionCoeff', np.zeros(2))
		self.gaussorder = kwargs.get('gaussOrder', 4)
		self.montecarlopts = kwargs.get('monteCarloPts', 100)
		self.coe = kwargs.get('coe', 0.0)
		self.terrainmodel = kwargs.get('terrainModel', 'RBF')
		self.poweridw = kwargs.get('powerIDW', 4)
		self.rbfkernel = kwargs.get('rbfKernel', 1)
		self.shapefactor = kwargs.get('shapeFactor', 5)
		self.octreemaxpts = kwargs.get('octreeMaxPts', 100)
		self.octreedepth = kwargs.get('octreeDepth', 1)
		self.terrainfile = kwargs.get('terrainFile', '')
		self.windrosefile = kwargs.get('windRoseFile', '')


		self.AEP = 0.0
		self.farmPower = 0.0
		self.farmEfficiency = 0.0
		self.farmCost = 0.0
		self.landUsed = 0.0
		self.COE = 0.0
		self.totalRatedPower = 0.0
		self.normalizedAEP = 0.0


		self.nTurbines = kwargs.get('nTurbines', 0)
		self.turbines = []
		for i in range(0, self.nTurbines):
			self.turbines.append( TurbineData( **kwargs ) )


		inputFile = kwargs.get('inputFile', '')
		if inputFile == '':
			tmpFile = tempfile.NamedTemporaryFile(mode = 'w+', delete = True)
			try:
				tmpFile.write('&windflo_data \n /')
				tmpFile.seek(0)			
				self.namelist = f90nml.read(tmpFile.name)			
			finally:
				tmpFile.close()
		else:	
			self.namelist = f90nml.read(inputFile)
			self.SetFromNamelist(**kwargs);


		resFile = kwargs.get('resFile', '')
		if resFile != '':
			self.ReadResultsFile(resFile)
			

		libPath = kwargs.get('libDir', '')
		if libPath != '':
			self.libWindFLO = cdll.LoadLibrary('./'+libPath+'libWindFLO.so')


	def __del__(self):
	
		try:
			self.libWindFLO.CleanAll_()
			dlclose_func = CDLL(None).dlclose
			dlclose_func.argtypes = [c_void_p]
			dlclose_func.restype = c_int
			dlclose_func(self.libWindFLO._handle)
		except:
			pass

	def ParseKwargsForAnalysisParams(self, **kwargs):
	
		self.rho = kwargs.get('rho', self.rho)
		self.turbulenceintensity = kwargs.get('turbulenceIntensity', self.turbulenceintensity)
		self.windmodel = kwargs.get('windModel', self.windmodel)
		self.modelvelocity = kwargs.get('modelVelocity', self.modelvelocity)
		self.surfaceroughness = kwargs.get('surfaceRoughness', self.surfaceroughness)
		self.referenceheight = kwargs.get('surfaceRoughness', self.referenceheight)
		self.wakemodel = kwargs.get('wakeModel', self.wakemodel)
		self.wakemergemodel = kwargs.get('wakeMergeModel', self.wakemergemodel)
		self.wakeexpansioncoeff = kwargs.get('wakeExpansionCoeff', self.wakeexpansioncoeff)
		self.gaussorder = kwargs.get('gaussOrder', self.gaussorder)
		self.montecarlopts = kwargs.get('monteCarloPts', self.montecarlopts)
		self.coe = kwargs.get('coe', self.coe)
		self.terrainmodel = kwargs.get('terrainModel', self.terrainmodel)
		self.poweridw = kwargs.get('powerIDW', self.poweridw)
		self.rbfkernel = kwargs.get('rbfKernel', self.rbfkernel)
		self.shapefactor = kwargs.get('shapeFactor', self.shapefactor)
		self.octreemaxpts = kwargs.get('octreeMaxPts', self.octreemaxpts)
		self.octreedepth = kwargs.get('octreeDepth', self.octreedepth)
		self.terrainfile = kwargs.get('terrainFile', self.terrainfile)
		self.windrosefile = kwargs.get('windRoseFile', self.windrosefile)


	def SetFromNamelist(self, **kwargs):
		
		if 'rho' in self.namelist['windflo_data']:		
			self.rho = self.namelist['windflo_data']['rho']
		
		if 'modelvelocity' in self.namelist['windflo_data']:		
			self.modelvelocity = np.array(self.namelist['windflo_data']['modelvelocity'])
			
		if 'windmodel' in self.namelist['windflo_data']:		
			self.windmodel = self.namelist['windflo_data']['windmodel']
				
		if 'gaussorder' in self.namelist['windflo_data']:		
			self.gaussorder = self.namelist['windflo_data']['gaussorder']
				
		if 'montecarlopts' in self.namelist['windflo_data']:		
			self.montecarlopts = self.namelist['windflo_data']['montecarlopts']
				
		if 'referenceheight' in self.namelist['windflo_data']:		
			self.referenceheight = self.namelist['windflo_data']['referenceheight']
				
		if 'surfaceroughness' in self.namelist['windflo_data']:		
			self.surfaceroughness = self.namelist['windflo_data']['surfaceroughness']
				
		if 'turbulenceintensity' in self.namelist['windflo_data']:		
			self.turbulenceintensity = self.namelist['windflo_data']['turbulenceintensity']
				
		if 'wakemodel' in self.namelist['windflo_data']:		
			self.wakemodel = self.namelist['windflo_data']['wakemodel']
				
		if 'wakemergemodel' in self.namelist['windflo_data']:		
			self.wakemergemodel = self.namelist['windflo_data']['wakemergemodel']
				
		if 'wakeexpansioncoeff' in self.namelist['windflo_data']:		
			self.wakeexpansioncoeff = np.array(self.namelist['windflo_data']['wakeexpansioncoeff'])
				
		if 'coe' in self.namelist['windflo_data']:		
			self.coe = self.namelist['windflo_data']['coe']
				
		if 'terrainmodel' in self.namelist['windflo_data']:		
			self.terrainmodel = self.namelist['windflo_data']['terrainmodel']
				
		if 'terrainfile' in self.namelist['windflo_data']:		
			self.terrainfile = self.namelist['windflo_data']['terrainfile']
				
		if 'poweridw' in self.namelist['windflo_data']:		
			self.poweridw = self.namelist['windflo_data']['poweridw']
				
		if 'rbfkernel' in self.namelist['windflo_data']:		
			self.rbfkernel = self.namelist['windflo_data']['rbfkernel']
				
		if 'shapefactor' in self.namelist['windflo_data']:		
			self.shapefactor = self.namelist['windflo_data']['shapefactor']
				
		if 'octreemaxpts' in self.namelist['windflo_data']:		
			self.octreemaxpts = self.namelist['windflo_data']['octreemaxpts']
				
		if 'octreedepth' in self.namelist['windflo_data']:		
			self.octreedepth = self.namelist['windflo_data']['octreedepth']
				
		if 'windrosefile' in self.namelist['windflo_data']:		
			self.windrosefile = self.namelist['windflo_data']['windrosefile']
			
		nTurbines = kwargs.get('nTurbines', 0)			
		if 'turbinefiles' in self.namelist['windflo_data'] and nTurbines == 0:
			self.nTurbines = len(self.namelist['windflo_data']['turbinefiles'])
		elif 'turbinefiles' in self.namelist['windflo_data']:
			nTurbines = len(self.namelist['windflo_data']['turbinefiles'])
			
		turbineFile = kwargs.get('turbineFile', '')
		self.turbines = []
		for i in range(0, self.nTurbines):
			self.turbines.append( TurbineData( params = np.zeros(30) ) )
			if turbineFile != '':
				self.turbines[i].SetFromNamelist(turbineFile)
			elif i < nTurbines:
				if os.path.exists(self.namelist['windflo_data']['turbinefiles'][i]):
					self.turbines[i].SetFromNamelist(self.namelist['windflo_data']['turbinefiles'][i])
		



	def WriteInputFile(self, **kwargs):
	
		runDir = kwargs.get('runDir', self.runDir)		
	
		self.namelist['windflo_data']['rho'] = self.rho
		self.namelist['windflo_data']['modelvelocity'] = self.modelvelocity.tolist()
		self.namelist['windflo_data']['windmodel'] = self.windmodel
		self.namelist['windflo_data']['gaussorder'] = self.gaussorder
		self.namelist['windflo_data']['montecarlopts'] = self.montecarlopts
		self.namelist['windflo_data']['referenceheight'] = self.referenceheight
		self.namelist['windflo_data']['surfaceroughness'] = self.surfaceroughness
		self.namelist['windflo_data']['turbulenceintensity'] = self.turbulenceintensity
		self.namelist['windflo_data']['wakemodel'] = self.wakemodel
		self.namelist['windflo_data']['wakemergemodel'] = self.wakemergemodel
		self.namelist['windflo_data']['wakeexpansioncoeff'] = self.wakeexpansioncoeff.tolist()
		self.namelist['windflo_data']['coe'] = self.coe
		self.namelist['windflo_data']['terrainmodel'] = self.terrainmodel
		self.namelist['windflo_data']['terrainfile'] = self.terrainfile
		self.namelist['windflo_data']['poweridw'] = self.poweridw
		self.namelist['windflo_data']['rbfkernel'] = self.rbfkernel
		self.namelist['windflo_data']['shapefactor'] = self.shapefactor
		self.namelist['windflo_data']['octreemaxpts'] = self.octreemaxpts
		self.namelist['windflo_data']['octreedepth'] = self.octreedepth
		self.namelist['windflo_data']['windrosefile'] = self.windrosefile

	
		self.namelist['windflo_data']['turbinefiles'] = []
		for i in range(0, self.nTurbines):
			self.namelist['windflo_data']['turbinefiles'].append( runDir +'turbine'+str(i+1)+'.inp' )
			self.turbines[i].WriteInputFile( runDir + 'turbine'+str(i+1)+'.inp' )

		filename = kwargs.get('inFile', 'WindFLO.inp')
		self.namelist.write(runDir + filename, force=True)

		


	def run(self,**kwargs):
		
		self.ParseKwargsForAnalysisParams(**kwargs)
		
		self.runDir = kwargs.get('runDir', self.runDir)			
		self.WriteInputFile(**kwargs)
			
		self.libWindFLO.Python_WindFLO_API.argtypes = [POINTER(c_char),		# infilename
												 	  POINTER(c_char),		# infilename
												  	  (c_int) ,				# nturbines
												  	  POINTER(c_double),	# velocities
												  	  POINTER(c_double),	# power
												  	  POINTER(c_double),	# rated power
												  	  POINTER(c_double)]	# outputs
		
		
		
		inFile = self.runDir + kwargs.get('inFile', 'WindFLO.inp')
		CinFile = np.asarray(inFile + '\0', dtype = c_char).ctypes.data_as(POINTER(c_char))

		resFile = self.runDir + kwargs.get('resFile', '')
		CoutFile = np.asarray(resFile + '\0', dtype = c_char).ctypes.data_as(POINTER(c_char))

		
		velocities = np.zeros(self.nTurbines*3, dtype = c_double)
		Cvelocities = velocities.ctypes.data_as(POINTER(c_double))
		
		power = np.zeros(self.nTurbines, dtype = c_double)		
		Cpower = power.ctypes.data_as(POINTER(c_double))

		ratedpower = np.zeros(self.nTurbines, dtype = c_double)		
		Cratedpower = ratedpower.ctypes.data_as(POINTER(c_double))


		outputs = np.zeros(100, dtype = c_double)				
		Coutputs = outputs.ctypes.data_as(POINTER(c_double))

		self.libWindFLO.Python_WindFLO_API(CinFile, CoutFile, c_int(self.nTurbines), Cvelocities, Cpower,Cratedpower, Coutputs)

		self.AEP = outputs[0]
		self.farmPower = outputs[1]
		self.farmEfficiency = outputs[2]
		self.farmCost = outputs[3]
		self.landUsed = outputs[4]
		self.COE = self.farmCost / self.AEP		
		self.totalRatedPower = np.sum(ratedpower)
		self.normalizedAEP = self.AEP / (self.totalRatedPower * 365.0 * 24.0)
		
		self.UpdateDict()
		
		clean = kwargs.get('clean', True)		
		k = 0
		for i in range(0, self.nTurbines):
			for j in range(0, 3):
				self.turbines[i].velocity[j] = velocities[k]
				k = k + 1
			self.turbines[i].power = power[i]
			self.turbines[i].ratedPower = ratedpower[i]
			self.turbines[i].UpdateDict()			
			if clean:			
				os.remove(self.namelist['windflo_data']['turbinefiles'][i])
		if clean:		
			os.remove(inFile)

			
	def clean(self):
		self.libWindFLO.Clean_()
	def cleanall(self):
		self.libWindFLO.CleanAll_()
	

	def ReadResultsFile(self, inFile):
	
		f = open(inFile,'r')
		line = f.readline()
		
		
		# Get the farm performance
		for i in range(0,6):
			line = f.readline()

			varName = line.split('=')[0].strip()
			varVal = float(line.split('=')[1])

			if(varName == 'nTurbines'):
				self.nTurbines = int(varVal)
			if(varName == 'AEP'):
				self.AEP = varVal
			if(varName == 'Power'):
				self.farmPower = varVal
			if(varName == 'Efficiency'):
				self.farmEfficiency = varVal
			if(varName == 'Cost'):
				self.farmCost = varVal
			if(varName == 'Land Used'):
				self.landUsed = varVal


		# Get turbines in farm
		line = f.readline()
		line = f.readline()
		self.turbines = []
		for i in range(0, self.nTurbines):
			line = f.readline()
			
			if(line.strip() == '$ConvexHull'):
				break
			
			turbineParams = [float(i) for i in line.split(',')]
			self.turbines.append( TurbineData( params = turbineParams ) )
		
		# Get convex hull params		
		line = f.readline()
		line = f.readline()

		self.convexHull = np.zeros((self.nTurbines,2))
		j = 0
		for i in range(0,self.nTurbines):
			line = f.readline()
			if(line.strip() == '$End'):
				break
		
			points = line.split(',')
			self.convexHull[j,0] = points[1]
			self.convexHull[j,1] = points[2]
			j = j + 1
	
		self.convexHull.resize((j,2))
		
		
		ratedPower = np.array([i.ratedPower for i in self.turbines])
		self.totalRatedPower = np.sum(ratedPower)
		self.farmPower = np.sum(np.array([i.power for i in self.turbines]))
		self.farmEfficiency = self.farmPower / self.totalRatedPower
		self.normalizedAEP = self.AEP / (self.totalRatedPower * 365.0 * 24.0)


#		Take Care of Units
		if self.AEP > 1e7:
			aepUnit = 'MWh'
		elif self.AEP > 1e4:
			aepUnit = 'kWh'
		else:
			aepUnit = 'Wh'

		if self.farmPower > 1e7:
			powerUnit = 'MW'
		elif self.farmPower > 1e4:
			powerUnit = 'kW'
		else:
			powerUnit = 'W'
		
		if self.landUsed > 1e6:
			areaUnit = 'km^2'
		else:
			areaUnit = 'm^2'

		self.COE = self.farmCost / self.AEP
		coeUnits = '$/'+aepUnit

		self.units = {'P': powerUnit,'Pr':powerUnit, 'A': areaUnit, 'AEP': aepUnit, 'COE': coeUnits}

	def UpdateDict(self):

#		Take Care of Units
		if self.AEP > 1e7:
			aepUnit = 'MWh'
		elif self.AEP > 1e4:
			aepUnit = 'kWh'
		else:
			aepUnit = 'Wh'

		if self.farmPower > 1e7:
			powerUnit = 'MW'
		elif self.farmPower > 1e4:
			powerUnit = 'kW'
		else:
			powerUnit = 'W'
		
		if self.landUsed > 1e6:
			areaUnit = 'km^2'
		else:
			areaUnit = 'm^2'

		coeUnits = '$/'+aepUnit

		self.units = {'P': powerUnit,'Pr':powerUnit, 'A': areaUnit, 'AEP': aepUnit, 'COE': coeUnits}


		self.scaling = {'MWh':1e-6, 'MW':1e-6, 'kWh':1e-3, 'kW':1e-3, 'Wh':1, 'W':1, 'km^2':1e-6, 'm^2':1,\
						'$/MWh':1e6, '$/kWh':1e3, '$/Wh':1}


		self.props = {'P': self.farmPower, 'A': self.landUsed, 'AEP': self.AEP, 'COE': self.COE, \
					  'nAEP': self.normalizedAEP, 'eff': self.farmEfficiency, 'cost':self.farmCost,\
					  'Pr': self.totalRatedPower}


	

	def getVar(self, variable):
	
		self.UpdateDict()
		value = self.props.get(variable, 0.0)
		unit = self.units.get( variable , '')	
		scale = self.scaling.get(unit,1)
		
		return value * scale, unit
		
		

	def plotWindFLO2D(self, fig, plotVariable = 'V', scale = 1.0, title = ''):
		
		ax = plt.subplot(1,1,1)
				
		ax.set_xlabel('x [m]', fontsize=16, labelpad = 5)
		ax.set_ylabel('y [m]', fontsize=16, labelpad = 5)	

		ax.tick_params(axis = 'x', which = 'major', labelsize = 15, pad=0)
		ax.tick_params(axis = 'y', which = 'major', labelsize = 15, pad=0)		

		x = np.array([i.position[0] for i in self.turbines])
		y = np.array([i.position[1] for i in self.turbines])
		var = np.array([LA.norm(i.variables[plotVariable]) for i in self.turbines])	* scale	

		jet_map = plt.get_cmap('jet')			
		scatterPlot = ax.scatter(x,y, c=var, marker = '^',  s = 100, cmap = jet_map, alpha = 1)

		if( (max(var) - min(var)) > 0):

			colorTicks = np.linspace(min(var), max(var), 7, endpoint = True)
			colorBar = plt.colorbar(scatterPlot,pad = 0.06,shrink = 0.8, format = ticker.FuncFormatter(fmt), ticks = colorTicks)	
		
			colorBar.ax.tick_params(labelsize=16)
			colorBar.ax.set_title(title, fontsize = 16,ha='left', pad = 15)
		#	tick_locator = ticker.MaxNLocator(nbins=7)
		#	colorBar.locator = tick_locator
			colorBar.update_ticks()


		plt.locator_params(axis='x', nbins=6)
		plt.locator_params(axis='y', nbins=6)

		plt.tight_layout()
		

		return ax



	def plotWindFLO3D(self, fig, plotVariable = ['V'], scale = [1.0], title = ['']):
		
		
		ax = plt.subplot(1,1,1, projection='3d')
				
		ax.set_xlabel('x [m]', fontsize=16, labelpad = 5)
		ax.set_ylabel('y [m]', fontsize=16, labelpad = 5)	
		ax.set_zlabel('z [m]', fontsize=16, labelpad = 10)	

		ax.tick_params(axis = 'x', which = 'major', labelsize = 15, pad=0)
		ax.tick_params(axis = 'y', which = 'major', labelsize = 15, pad=0)		
		ax.tick_params(axis = 'z', which = 'major', labelsize = 15, pad=5)


		nPlot = len(plotVariable)
	
		x = np.array([i.position[0] for i in self.turbines])
		y = np.array([i.position[1] for i in self.turbines])
		z = np.array([i.position[2] for i in self.turbines])
		H = np.array([i.height for i in self.turbines])
		var1 = np.array([LA.norm(i.variables[plotVariable[0]]) for i in self.turbines])	* scale[0]
	
		for i in range(0,len(x)):
			ax.plot( [x[i],x[i]], [y[i],y[i]], [z[i],z[i]+H[i]], markersize = 8, markerfacecolor = 'none', color = 'gray' )

		jet_map = plt.get_cmap('jet')


		if nPlot > 1:
			var2 = np.array([LA.norm(i.variables[plotVariable[1]]) for i in self.turbines])	* scale[1]
			minMarkerSize = 50
			maxMarkerSize = 150
			markerSize = minMarkerSize + (var2 - min(var2)) * (maxMarkerSize / (max(var2) - min(var2)))
		else:
			markerSize = 100


		scatterPlot = ax.scatter(x,y,z+H, c=var1, marker = 'o',  s = markerSize, cmap = jet_map, alpha = 1)

		colorTicks = np.linspace(min(var1), max(var1), 7, endpoint = True)
		colorBar = plt.colorbar(scatterPlot,pad = 0.0,shrink = 0.8, format = ticker.FuncFormatter(fmt), ticks = colorTicks)	
		
		colorBar.ax.tick_params(labelsize=15)
		colorBar.ax.set_title(title[0], fontsize = 16,ha='left', pad = 10)
#		tick_locator = ticker.MaxNLocator(nbins=7)
#		colorBar.locator = tick_locator
		colorBar.update_ticks()


		if nPlot > 1:

			# produce a legend with a cross section of sizes from the scatter
			tmp = np.linspace(min(var2),max(var2), 5, endpoint = True)			
#			tmp = np.arange(min(var2),max(var2), (max(var2) - min(var2)) / 5)
			tmpMs =  minMarkerSize + (tmp - min(tmp)) * (maxMarkerSize / (max(tmp) - min(tmp)))	
			labels = ['{0:.1f}'.format(s) for s in tmp]
			points = [ax.scatter([], [], [], s=s, c='gray') for s in tmpMs]
			sizeLegend = plt.legend(reversed(points), reversed(labels), scatterpoints=1,loc=(-.35,0.2), fontsize = 15,frameon=True)	
			sizeLegend.set_title(title[1], prop = {'size':'16'})
			sizeLegend.get_frame().set_linewidth(2.0)


		
		plt.locator_params(axis='x', nbins=6)
		plt.locator_params(axis='y', nbins=6)

		ax.view_init(azim=-150,elev=30.)	

		plt.subplots_adjust(left=0.25, right=1., top=0.9, bottom=0.01)

		
		return ax
		

		
	def getScaledPerformance(self, scale = False):
	
		AEP = 0.0
		farmPower = 0.0
		landUsed = 0.0
		
		if scale == False:
			AEP = self.AEP 
			farmPower = self.farmPower 
			landUsed = self.landUsed

			return AEP, farmPower, landUsed
			
		
		if self.units['AEP'] == 'MW':
			AEP = self.AEP / 1e6
		elif self.units['AEP'] == 'kW':
			AEP = self.AEP / 1e3
		else:
			AEP = self.AEP 

		if self.units['P'] == 'MW':
			farmPower = self.farmPower / 1e6
		elif self.units['P'] == 'kW':
			farmPower = self.farmPower / 1e3
		else:
			farmPower = self.farmPower 
		
		if 	self.units['A'] == 'km^2':
			landUsed = self.landUsed / 1e6
		else:
			landUsed = self.landUsed
		

		return AEP, farmPower, landUsed

	def annotatePlot(self, ax):

		AEP, power, landArea = self.getScaledPerformance(scale = True)

		COE, _ = self.getVar('COE')
		AEP, _ = self.getVar('AEP')
		power, _ = self.getVar('P')
		landArea, _ = self.getVar('A')

		ax.annotate( r'CoE = '+str('{0:.2f} '.format(COE)) + self.units['COE'], (0.0,0.0), (-.1,1.08), fontsize = 16, xycoords='axes fraction',va='top')
#		ax.annotate( r'$AEP = $'+str('{0:.2f}'.format(AEP)) + self.units['AEP'], (0.0,0.0), (0.,1.08), fontsize = 16, xycoords='axes fraction',va='top')
		ax.annotate( r'$\overline{\eta}_{Farm} = $'+str('{0:.2f}'.format(self.farmEfficiency)) , (0.0,0.0), (0.5,1.08), fontsize = 16, xycoords='axes fraction',va='top')   
		ax.annotate( r'$A_{Farm} = $'+str('{0:.2f}'.format(landArea)) + r' $'+ self.units['A']+'$', (0.0,0.0), (0.85,1.08), fontsize = 16, xycoords='axes fraction',va='top')   
#		ax.annotate( r'$AEP = $'+str('{0:.0f}'.format(self.AEP)) + ' kWH', (0.0,0.0), (-0.15,0.07), fontsize = 16)
#		ax.annotate( r'$\overline{Cost}_{Farm} = $'+str('{0:.0f}'.format(self.farmCost/self.farmPower)) + ' USD/kW', (0.0,0.0), (-0.15,0.07), fontsize = 16)

		plt.tight_layout()


	def plotConvexHull(self, ax):
		ax.plot(self.convexHull[:,0],self.convexHull[:,1], 'k-')

