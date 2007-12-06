#! /usr/bin/env python

import os
import sys
import numpy as np
from scipy.optimize import leastsq
#import pylab

usage = """FitStrains(datafile, tensor, cellvolume)
First arg is the name of the data file
(two column matrix of strain amplitude and energy). 
Second arg is deformation tensor, third arg is cellvolume.
"""

def residuals(p, y, x): 
	err = y-peval(x,p) 
	return err

def peval(x, p): 
	return p[0] + p[1] * x**2

def numlistequal(x,y,eps):
	if len(x) <> len(y): return False
	for i in range(len(x)):
		if abs(x[i]-y[i]) > eps: return False
	return True


def FitParabola(data):
	"""Fits a parabola y = a + b*x^2 to two-column data."""
    
	data = np.array(data)
    
	xdata = data[:,0]
	ydata = data[:,1]

	#    polycoeffs = scipy.polyfit(xdata, ydata, 2)
	# We do not want to have any linear coefficient to fit with 
	# as it should be fixed to zero

	# instead we use the peval function above with only zero-th and second order terms:
	pname = (['a','b'])
	p0 = [0.0, 0.0]
	plsq = leastsq(residuals, p0, args=(ydata, xdata), maxfev=1000)
	yfit=peval(xdata, plsq[0])
	#pylab.plot(xdata, ydata, 'k.')
	#pylab.plot(xdata, yfit, 'r-')
	b = plsq[0][1]
	return b

def FitCubicElasticConstant(data, tensor, cellvolume):
	"""Fits the elastic constant for a given strain type,
	corresponding to the deformation tensor.
	Returns the elastic constant in eV/A^2.
	The data must be in eV vs A. Cellvolume is in A^3."""
	
	b = FitParabola(data)
	tensor = np.array(tensor)
	eps = 0.001

	# Definitions of possible strain deformation tensors:
	import CubicStandardStrainTensors as cubic    

	if numlistequal(tensor.flatten(), cubic.C11.flatten(), eps):
		return (2 * b / cellvolume)
	elif numlistequal(tensor.flatten(), cubic.C11C12_1.flatten(), eps):
		return (b / cellvolume) 
	elif numlistequal(tensor.flatten(), cubic.C11C12_2.flatten(), eps):
		return (b / (3 * cellvolume))
	elif numlistequal(tensor.flatten(), cubic.C44_1.flatten(), eps):
		return (b / (2 * cellvolume))
	elif numlistequal(tensor.flatten(), cubic.C44_2.flatten(), eps):
		return (b / (6 * cellvolume))
	else:
		return 0

def FitCubicStrain(datafilename, straintype, cellvolume):
	"""Fits the strain data from datafile to extract corresponding elastic constants.
	The straintype is a string corresponding to standard cubic strain types."""
    
	import CubicStandardStrainTensors as cubic
	from pylab import load
	try:
		data = load(datafilename)
	except: print "Could not load strain data."

	tensor = cubic.__dict__[straintype]
	elastConst = FitCubicElasticConstant(data, tensor, cellvolume)
	return elastConst

        
if __name__ == "__main__":
	try:
		datafile=sys.argv[1] ; tensor=sys.argv[2] ; cellvolume=sys.argv[3]
		FitCubicStrain(sys.argv[1], sys.argv[2], sys.argv[3])
	except:
		print "Usage:", usage
		
	
    
