import sys
sys.path.append('../../API/')

import numpy as np
from WindFLO import WindFLO
import matplotlib.pyplot as plt


if __name__ == "__main__":


	windFLO = WindFLO(resFile = 'WindFLO.res')	
	fig = plt.figure(figsize=(8,5), edgecolor = 'gray', linewidth = 2)
	ax = windFLO.plotWindFLO2D(fig, plotVariable = 'P', scale = 1.0e-3, title = 'P [kW]')
	windFLO.annotatePlot(ax)
	plt.show()
	


