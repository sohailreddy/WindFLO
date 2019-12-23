#!/usr/bin/python3

import matplotlib
import matplotlib.pyplot as plt
import scipy as sp
import math as mth
import os
import numpy as np
import sys
from matplotlib  import cm
from numpy import linalg as LA
import mpl_toolkits.mplot3d.art3d as art3d
import matplotlib.ticker as ticker


class TurbineData:

	turbineNum = 0
#	position = sp.zeros(3)
#	velocity = sp.zeros(3)
	height = 0.0
	radius = 0.0
	area = 0.0
	ratedPower = 0.0
	power = 0.0
	fictitious = False

	def __init__(self, params):
		self.turbineNum = params[0]

		self.position = sp.zeros(3)
		self.velocity = sp.zeros(3)
		
		self.position[0] = params[1]
		self.position[1] = params[2]
		self.position[2] = params[3]
		self.velocity[0] = params[4]
		self.velocity[1] = params[5]
		self.velocity[2] = params[6]
		self.height = params[7]
		self.radius = params[8]
		self.area = params[9]
		self.ratedPower = params[10]
		self.power = params[11]
		self.fictitious = params[12]


def fmt(x, pos):
	return str('{0:.2f}'.format(x))

def plotTurbine(ax, turbineData):
		
	ax.set_xlabel('x [m]', fontsize=16, labelpad = 5)
	ax.set_ylabel('y [m]', fontsize=16, labelpad = 5)	
	ax.set_zlabel('z [m]', fontsize=16, labelpad = 10)	

	ax.tick_params(axis = 'x', which = 'major', labelsize = 15, pad=0)
	ax.tick_params(axis = 'y', which = 'major', labelsize = 15, pad=0)		
	ax.tick_params(axis = 'z', which = 'major', labelsize = 15, pad=5)

	x = turbineData[:,1]
	y = turbineData[:,2]
	z = turbineData[:,3]
	V = LA.norm(turbineData[:,4:7], axis = 1)
	H = turbineData[:,7]
	R = turbineData[:,8]
	P = turbineData[:,11]
	


	for i in range(0,len(x)):
		ax.plot( [x[i],x[i]], [y[i],y[i]], [z[i],z[i]+H[i]], markersize = 8, markerfacecolor = 'none', color = 'gray' )
	

	minMarkerSize = 50
	maxMarkerSize = 150
	markerSize = minMarkerSize + (P - min(P)) * (maxMarkerSize / (max(P) - min(P)))


	jet_map = plt.get_cmap('jet')
	scatterPlot = ax.scatter(x,y,z+H, c=V, marker = 'o', s = markerSize, cmap = jet_map, alpha = 1)
	vBarTicks = sp.arange(min(V), max(V), (max(V) - min(V))/6)
	velocityBar = plt.colorbar(scatterPlot,pad = -0.05,shrink = 0.8, format = ticker.FuncFormatter(fmt), ticks = vBarTicks)	
	velocityBar.ax.tick_params(labelsize=16)
	velocityBar.ax.set_title('V [m/s]', fontsize = 16)
	tick_locator = ticker.MaxNLocator(nbins=7)
	velocityBar.locator = tick_locator
	velocityBar.update_ticks()
	
	
	# produce a legend with a cross section of sizes from the scatter
	tmp = np.arange(min(P),max(P), (max(P) - min(P)) / 5)
	tmpMs =  minMarkerSize + (tmp - min(tmp)) * (maxMarkerSize / (max(tmp) - min(tmp)))	
	labels = ['{0:.3f}'.format(s) for s in tmp]
	points = [ax.scatter([], [], [], s=s, c='gray') for s in tmpMs]
	sizeLegend = plt.legend(reversed(points), reversed(labels), scatterpoints=1,loc=(-.35,0.2), fontsize = 15,frameon=True)	
	sizeLegend.set_title("P [W]", prop = {'size':'16'})
	sizeLegend.get_frame().set_linewidth(2.0)
#	ax.w_zaxis.line.set_lw(0.)
#	ax.set_zticks([])
#	ax.grid(False)

	plt.locator_params(axis='x', nbins=6)
	plt.locator_params(axis='y', nbins=6)
#	plt.locator_params(axis='z', nbins=6)		

#	ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
#	ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
#	ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))

	ax.view_init(azim=-150,elev=30.)	

	plt.subplots_adjust(left=0.25, right=1., top=0.9, bottom=0.01)

def plotConvexHull(ax, convexHull):

	ax.plot(convexHull[:,1],convexHull[:,2], 'k-')

def annotatePlot(ax, farmPower, farmEfficiency, farmCost, landUsed):

    ax.annotate( r'$P_{Farm} = $'+str('{0:.2f}'.format(farmPower)) + ' W', (0.0,0.0), (-0.15,0.1), fontsize = 16)
    ax.annotate( r'$\overline{\eta}_{Farm} = $'+str('{0:.2f}'.format(farmEfficiency)) , (0.0,0.0), (-0.15,0.085), fontsize = 16)   
    ax.annotate( r'$\overline{Cost}_{Farm} = $'+str('{0:.0f}'.format(farmCost/farmPower)) + ' USD/kW', (0.0,0.0), (-0.15,0.07), fontsize = 16)
    ax.annotate( r'$A_{Farm} = $'+str('{0:.2f}'.format(landUsed)) + r' $m^2$', (0.0,0.0), (-0.15,0.055), fontsize = 16)

	

if __name__ == "__main__":

	Case = 'Case0'
	Model = 'Frandsen'
	Dir = ''#Case+'/'+Model+'/'

	WindFLO_res = sp.loadtxt(Dir+"WindFLO.res", delimiter = '=', skiprows = 1, usecols = 1, max_rows = 4)
	turbineData = np.genfromtxt(Dir+"WindFLO.res", delimiter = ',',skip_header=7, skip_footer=1,invalid_raise=False)
	convexHull = np.genfromtxt(Dir+"WindFLO.res", delimiter = ',',skip_header=9+len(turbineData[:,0]), skip_footer=1,loose=True)
	
	turbineData = sp.array(turbineData[turbineData[:,12] == 0,:])
	
	farmPower = WindFLO_res[0]
	farmEfficiency = WindFLO_res[1]	
	farmCost = WindFLO_res[2]
	landUsed = WindFLO_res[3]


	fig = plt.figure(figsize=(8,5), edgecolor = 'gray', linewidth = 2)
	ax = plt.subplot(1,1,1, projection='3d')
	
	windTurbines = []
	for i in range(0, len(turbineData[:,0])):
		windTurbines.append( TurbineData( turbineData[i,:] ) )
	
	plotTurbine(ax, turbineData)
	plotConvexHull(ax, convexHull)
	annotatePlot(ax, farmPower, farmEfficiency, farmCost, landUsed)
	
	plt.savefig(Case+'_'+Model +'.png',  edgecolor=fig.get_edgecolor())
#	plt.show()			
		
		
