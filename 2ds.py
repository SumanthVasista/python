from __future__ import division
import math as mt
import numpy as np
import pickle
import os
import sys

cwl = os.path.dirname(os.path.realpath(__file__))

class Patient:
	def __init__(self):
		self.id = None
		self.response = None
		self.position = None
		self.ktag = None

class color:
   PURPLE = '\033[95m'
   CYAN = '\033[96m'
   DARKCYAN = '\033[36m'
   BLUE = '\033[94m'
   GREEN = '\033[92m'
   YELLOW = '\033[93m'
   RED = '\033[91m'
   BOLD = '\033[1m'
   UNDERLINE = '\033[4m'
   END = '\033[0m'

patient_features_create = []
new_patient_features_create = []

def qualify(a_file):
	source_cluster(a_file)
	print "File found. Checking Quality Metrics..."
	m1, m2, m, c = create_calibration_line()
	#This module qualifies a(the) base cluster.dat file that represents the model produced results for the drug
	#Checks the 1D-PPV, 1D-NPV, 2D-alpha, oneDCS, sensitivity and judges usability of conditional-usability of the model
	p_match = 0
	p_list = 0
	n_match = 0
	n_list = 0
	outlier = 0
	two_d_match = 0
	usage_signal = 1
	r_xlist = []
	n_xlist = []
	warnings = []
	conditions = []
	for element in patient_features_create:
	#1. Checking 1D-PPV - This is 1D
		if element.response == 'R':
			r_xlist.append(element.position[0]) #4. Sensitivity
			p_list += 1
			if element.position[0] > 0:
				p_match += 1
			if element.position[1] < (element.position[0]*m + c): #3. 2D-alpha
				outlier += 1
	#2. Checking 1D-NPV - This is 1D
		if element.response == 'N':
			n_xlist.append(element.position[0]) #4. Sensitivity
			n_list += 1
			if element.position[0] < 0:
				n_match += 1
			if element.position[1] > (element.position[0]*m + c): #3. 2D-alpha
				outlier += 1
		predicted_response = return_discreet_response(m1, m2, m, c, element.position[0], element.position[1])
		if predicted_response == element.response:
			two_d_match += 1
	PPV = (p_match/p_list)*100
	NPV = (n_match/n_list)*100
	alpha = outlier / len(patient_features_create) * 100
	spread = float(np.std(r_xlist) + np.std(n_xlist)*0.5*100)
	sensitivity = float((np.mean(r_xlist) - np.mean(n_xlist))*100)
	s_by_p = sensitivity / spread
	oneDCS = ((p_match + n_match)/(p_list + n_list)) * 100
	twoDCS = (two_d_match/len(patient_features_create)) * 100
	#Conditions for non-usability of model
	if PPV < 60 or NPV < 60:
		usage_signal = 0
	if 55 < oneDCS < 65:
		warnings.append(str('Low 1D-cluster score, consider re-modeling - proceed with caution'))
	if oneDCS < 55:
		usage_signal = 0 
	if twoDCS < 75:
		usage_signal = 0
	if alpha > 15:
		warnings.append(str('High overlap of R/N - proceed with caution'))
	if usage_signal == 0:
		if PPV > 70:
			conditions.append(str('Usage for Responder probability testing'))
		if NPV > 70:
			conditions.append(str('Usage for Non-Responder probability testing'))
	print "Writing report"
	print "\n"
	print "1-D QUALIFICATION METRICS"
	print "1D-CLUSTERING CORREALTION: ", oneDCS, "%", " (CUT-OFF - 65%)"
	print "1D-PPV: ", PPV, "%", " (CUT-OFF - 60%)"
	print "1D-NPV: ", NPV, "%", " CUT-OFF - 60%)"
	print "SENSITIVITY: ", sensitivity, " (CUT-OFF - ~)"
	print "SPREAD: ", spread, " (CUT-OFF - ~)"
	print "s/p RATIO: ", s_by_p, " (CUT-OFF - ~)"
	print "\n"
	print " 2-D QUALIICATION METRICS"
	print "2D-CLUSTERING CORREALTION", twoDCS, "%", " (CUT-OFF - 75%)"
	print "alpha (p-value): ", alpha, "%", " (CUT-OFF - 15%)"
	print "\n"
	if usage_signal == 1 and len(warnings) == 0:
		print "USAGE_SIGNAL: " + color.GREEN + color.BOLD + "GREEN" + color.END
	elif usage_signal == 1 and len(warnings) > 0:
		print "USAGE_SIGNAL: " + color.YELLOW + color.BOLD + "YELLOW" + color.END
	else:
		print "USAGE_SIGNAL: " + color.RED + color.BOLD + "RED" + color.END
	if usage_signal == 0:
		print "CONDITIONAL_USAGE:"
		print "\n"
		for message in conditions:
			print message
	print "OTHER WARNINGS:"
	for message in warnings:
		print message
	print "\n"
	print "~~END~~"	

def source_cluster(a_file):
		patientlist = [patientlist.rsplit(' ')[0] for patientlist in open(a_file)]
		x_score = [x_score.rsplit(' ')[1] for x_score in open(a_file)]
		x_score = map(float, x_score)
		y_score = [y_score.rsplit(' ')[2] for y_score in open(a_file)]
		y_score = map(float, y_score)
		try:
			response = [response.rsplit(' ')[3].rstrip('\n') for response in open(a_file)]
		except:
			print "non labelled data detected. Please enter labelled data to create features-separator"
			exit(-1)
		for i in range(len(patientlist)):
			name = str(patientlist[i])
			patientlist[i] = Patient()
			patientlist[i].id = name
			if response[i] == '1':
				patientlist[i].response = 'R'
			elif response[i] == '0':
				patientlist[i].response = 'N'
			patientlist[i].position = x_score[i], y_score[i]
			patient_features_create.append(patientlist[i])

def source_newdata(a_file):
		new_patientlist = [new_patientlist.rsplit(' ')[0] for new_patientlist in open(a_file)]
		new_x_score = [new_x_score.rsplit(' ')[1] for new_x_score in open(a_file)]
		new_x_score = map(float, new_x_score)
		new_y_score = [new_y_score.rsplit(' ')[2] for new_y_score in open(a_file)]
		new_y_score = map(float, new_y_score)
		print "Finished sourcing new data"
		for i in range(len(new_patientlist)):
			name = str(new_patientlist[i])
			new_patientlist[i] = Patient()
			new_patientlist[i].id = name
			new_patientlist[i].position = new_x_score[i], new_y_score[i]
			new_patient_features_create.append(new_patientlist[i])

def create_calibration_line():
	'''calculate y-intercept'''
	y_rlist = []
	y_nlist = []
	for element in patient_features_create:
		if element.response == 'R':
			y_rlist.append(element.position[1])
		elif element.response == 'N':
			y_nlist.append(element.position[1])
	mean_y_rlist = np.mean(y_rlist)
	std_y_rlist = np.std(y_rlist)
	mean_y_nlist = np.mean(y_rlist)
	std_y_nlist = np.std(y_rlist)
	y_rlist_abs = []
	y_nlist_abs = []
	for element in y_rlist:
		if element > (mean_y_rlist - 1.2*std_y_rlist):
			y_rlist_abs.append(element)
	for element in y_nlist:
		if element < (mean_y_nlist + 1.2*std_y_nlist):
			y_nlist_abs.append(element)
	c = (max(y_nlist_abs)+min(y_rlist_abs))/2
	'''calculate m (slope)'''
	x_rlist = []
	x_nlist = []
	for element in patient_features_create:
		if element.response == 'R':
			x_rlist.append(element.position[0])
		elif element.response == 'N':
			x_nlist.append(element.position[0])
	mean_x_rlist = np.mean(x_rlist)
	std_x_rlist = np.std(x_rlist)
	mean_x_nlist = np.mean(x_rlist)
	std_x_nlist = np.std(x_rlist)
	x_rlist_abs = []
	x_nlist_abs = []
	for element in x_rlist:
		if element > (mean_x_rlist - 1.2*std_x_rlist):
			x_rlist_abs.append(element)
	for element in x_nlist:
		if element < (mean_y_nlist + 1.2*std_y_nlist):
			x_nlist_abs.append(element)
	for element in patient_features_create:
		if element.position[0] == min(x_nlist_abs):
			y_high = element.position[1]
		if element.position[0] == max(x_rlist_abs):
			y_low = element.position[1]
	p1 = min(x_nlist_abs), float(y_low-0.1)
	p2 = 0, c
	p3 = max(x_rlist_abs), float(y_high+0.1)

	m1 = (p2[1] - p1[1])/(p2[0] - p1[0])
	m2 = (p3[1] - p2[1])/(p3[0] - p2[0])

	m = (m1+m2)/2
	return m1, m2, m, c

def export_co_efficients():
	if os.path.isfile(cwl+'/cluster.dat') == True:
		source_cluster(cwl+'/cluster.dat')
		m1, m2, m, c = create_calibration_line()
		print m1, m2, m, c

def return_discreet_response(m1, m2, m, c, patient_x, patient_y):
	if float(patient_y) > float(m*patient_x) + c:
		n_response = 'R'
	elif float(patient_y) < float(m*patient_x) + c:
		n_response = 'N'
	else:
		n_response = 'A'
	return n_response

def export_reports(a_file):
	if os.path.isfile(cwl+'/cluster.dat') == True:
		source_cluster(cwl+'/cluster.dat')
		m1, m2, m, c = create_calibration_line()
		print m1, m2, m, c
		print "Writing reports..."
		source_newdata(a_file)
		#Writing reports
		target = open(cwl+'/report.csv','w')
		str0 = str('PatientID,Relative_Efficacy,Predicted_Response')
		target.write(str0)
		target.write('\n')
		for element in new_patient_features_create:
			dist = (abs(element.position[0]*m + element.position[1]*-1 + c)/(mt.sqrt(m*m + 1)))
			if element.position[1] > m*element.position[0] + c:
				n_response = 'R'
			elif element.position[1] < m*element.position[0] + c:
				n_response = 'N'
			else:
				n_response = 'A'
			if n_response == 'R':
				beta = 1
			elif n_response == 'N':
				beta = -1
			relative_efficacy = (beta*mt.tanh(dist)*100) + 30
			string = str(str(element.id)+','+str(relative_efficacy)+','+str(n_response))
			target.write(string)
			target.write('\n')
		target.close()
		print "Finished printing report"


if __name__=="__main__":
	script, mode, newpath = sys.argv
	if sys.argv[1] == 'viz':
		export_co_efficients()

	if sys.argv[1] == 'report':
		export_reports(sys.argv[2])
				
