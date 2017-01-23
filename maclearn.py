#Python
from __future__ import division
import os sys itertools math numpy pickle

#Math functions
def variance(alist,type):
	alist_adjust = []
	var = 1
	mean_alist = float(sum(alist))/len(alist)
	for i in range(len(alist)):
		alist_adjust.append(alist[i] - mean_alist)
		var = var*alist_adjust[i]
	if type == 1:
		return alist_adjust
	elif type == 2:
		return var/(len(alist) - 1)

def covariance(alist, blist):
	covar = 1
	if len(alist) == len(blist):
		a = variance(alist,1)
		b = variance(blist,1)
		for i in range(len(alist)):
			covar = covar*a[i]*b[i]
		return float(covar/(len(alist) - 1))
	else:
		return 0

def update_covariance(value,beta,n,m):
	if n == 0:
		n = 0.1
	e = 2.71828
	exponent = (((n - m)/n) - 1)
	return float(value + (beta * (e**exponent) * value))

def gDimIndex(covar, value):
	if covar > 0:
		if value > 0:
			RDI = -1 #This is taken as negative as we are calculating dissimilarity
			NDI = 1
		elif value < 0:
			RDI = 0.5
			NDI = -0.5
	elif covar < 0:
		if value > 0:
			RDI = -1
			NDI = 1
		elif value < 0:
			RDI = -0.5
			NDI = 0.5
	return RDI, NDI

class PathwayVariable:
	def __init__(self):
		self.tuple = []
		self.mean = None
		self.mean_adjust = []
		self.covar = None
		self.match = None #n
		self.mismatch = None #m

	def extract_variance_values(alist, blist):
		if len(alist) == len(blist):
			for i in range(len(alist)):
				couple = (alist[i], blist[i])
				self.tuple.append(couple)
				couple = None

	def calculate_covariance:
		pathway_value = []
		response = []
		for i in range(len(self.tuple)):
			pathway_value.append(self.tuple[i][0])
			response.append(self.tuple[i][1])
		self.covar = covariance(pathway_value,response)

	def update_covariance(beta):
		self.covar = update_covariance(self.covar,beta,self.match,self.mismatch)




