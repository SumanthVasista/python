from __future__ import division
import math as mt
import pickle

class Patient:
	def __init__(self):
		self.response = None
		self.x_score = None
		self.y_score = None
		self.ktag = None

patients_feature_create = []

def source_cluster(a_file):
		patientlist = [patientlist.rsplit(' ')[0] for patientlist in open(a_file)]
		x_score = [x_score.rsplit(' ')[1] for x_score in open(a_file)]
		y_score = [y_score.rsplit(' ')[2] for y_score in open(a_file)]
		try:
			response = [response.rsplit(' ')[3] for response in open(a_file)]
		except:
			print "non labelled data detected. Please enter labelled data to create features-separator"
			exit(-1)
		for i in range(len(patientlist)):
			patientlist[i] = Patient()
			if response[i] == 1:
				patientlist[i].response = 'R'
			elif response[i] == 0:
				patientlist[i].response = 'N'
			patientlist[i].x_score = x_score[i]
			patientlist[i].y_score = y_score[i]
			patients_feature_create.append(patientlist[i])


def create_feature_line():
	'''calculate y-intercept'''
	y_rlist = []
	y_nlist = []
	for element in patients_feature_create:
		if element.response == 'R':
			y_rlist.append(element.y_score)
		if element.response == 'N':
			y_nlist.append(element.y_score)
	c = (max(y_nlist)+min(y_rlist))/2
	


if __name__=="__main__":
