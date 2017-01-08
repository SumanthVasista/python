#Python program to cluster mutations and pathways with drug response
from __future__ import division
import sys
import pickle
import os
import math

cwl = os.path.dirname(os.path.realpath(__file__))

#Defining global variales
training_list = []
pathways_initialized = {}
pathways_selected_for_patients = {}

def UpdatePathwayScores(default_location):
	'''Extracting information from Universal Pathways Repository.Generating a pickle file for pathways score once'''
	if os.path.isfile(default_location+'/list') == False:
		print "Pathway Relationsip Score data not found. Creating from scratch"
		loc_pathways = raw_input("Please enter location of pathways score file: ")
		if os.path.isfile(loc_pathways+'/list') == True:
			pathways = [pathways.split(',')[0] for pathways in open(loc_pathways+'/list')]
			target2 = open(default_location+'/list','w')
			for listed_paths in pathways:
				pathname = str(listed_paths)
				print pathname
				genes_in_pathways = [genes_in_pathways.split(',')[0] for genes_in_pathways in open(loc_pathways+'/'+listed_paths)]
				score_in_pathways = [score_in_pathways.split(',')[1] for score_in_pathways in open(loc_pathways+'/'+listed_paths)]
				score_in_pathways_strip = []
				for elements in score_in_pathways:
					score_in_pathways_strip.append(float(elements))
				listed_paths = dict(zip(genes_in_pathways, score_in_pathways_strip))
				target = open(default_location+'/'+pathname,'w')
				pickle.dump(listed_paths,target)
				target2.write(pathname)
				target2.write('\n')
				target.close()
			target2.close()
			print "Finished creating Pathway Relationship Score data files. Updated in location", default_location
		else:
			print "No list found in", loc_pathways, "to update pathways"
	else:
		print "Pathway Relationship Score data found. Proceeding..."

class Patient:
	def __init__(self):
		self.type = None
		self.gender = None
		self.disease = None
		self.response = None
		self.mutations = {}
		self.pathways = {}
		list_of_pathways = [list_of_pathways.rstrip('\n') for list_of_pathways in open(cwl+'/UpdatedPathwaysScores/list')]
		for component in list_of_pathways:
			pathway_name = str(component)
			self.pathways[pathway_name] = 0

	def source_mutation(self,a_list):
		gene = [gene.split(',')[1] for gene in open(a_list)]
		mut = [mut.split(',')[2] for mut in open(a_list)]
		disease = dict(zip(gene, mut))
		for component in gene:
			if disease[component] == 'OE' or disease[component] == 'OE\n':
				self.mutations[component] = 'GOF'
			elif disease[component] == 'KD' or disease[component] == 'KD\n':
				self.mutations[component] = 'LOF'
		print "Finished sourcing mutations for patient"

	def generate_pathways_signature(self,id,source,dest):
		for mutation in self.mutations:
			for pathway in self.pathways:
				target = open(source+'/'+pathway,'r')
				expanded_pathway = pickle.load(target)
				target.close()
				if mutation in expanded_pathway:
					if self.mutations[mutation] == 'GOF':
						self.pathways[pathway] += expanded_pathway[mutation]
					elif self.mutations[mutation] == 'LOF':
						self.pathways[pathway] += -1*expanded_pathway[mutation]
		target = open(dest+'/'+str(id), 'w')
		pathways_selected_for_patients = dict((k, v) for k, v in self.pathways.items() if v >= 1 or v <= -1)
		if len(pathways_selected_for_patients) < 1:
			pickle.dump(self.pathways,target)
			print "Patient does not express any strong pathways signatures... Proceeding with weak signal..."
		else:
			pickle.dump(pathways_selected_for_patients,target)
		pathways_selected_for_patients = {}
		target.close()
		target2 = open(dest+'/list','w')
		target2.write(str(id))
		target2.write('\n')
		target2.close()
		print "Finished generating patient specific pathways data for", id

def Categorize_Drug_Response(givenlist):
	drug_specific_pathways = {}
	list_of_pathways = [list_of_pathways.rstrip('\n') for list_of_pathways in open(cwl+'/UpdatedPathwaysScores/list')]
	for component in list_of_pathways:
		pathway_name = str(component)
		drug_specific_pathways[pathway_name] = 0
	for patient in givenlist:
		patient_name = str(patient)
		patient = Patient()
		patient.response = str(patient_response[patient_name])
		patient.source_mutation(patient_msource[patient_name])
		patient.generate_pathways_signature(patient_name, cwl+'/UpdatedPathwaysScores', cwl+'/PatientPathwaysTraining')
		target = open(cwl+'/PatientPathwaysTraining/'+patient_name,'r')
		patient_pathways = {}
		patient_pathways = pickle.load(target)
		target.close()
		for pathway in patient_pathways:
			if patient.response == 'R':
				drug_specific_pathways[pathway] += 0.2*float(patient_pathways[pathway])
			elif patient.response == 'N':
				drug_specific_pathways[pathway] += -0.2*(float(patient_pathways[pathway]))
	pathways_selected_for_patients = dict((k, v) for k, v in drug_specific_pathways.items() if v >= 0.5 or v <= -0.5)
	print pathways_selected_for_patients
	target = open(cwl+'/PatientPathwaysTraining/DrugPathways','w')
	pickle.dump(pathways_selected_for_patients,target)
	target.close()
	print "Finished clustering pathways to drug response"

def ForceWeightage():
	if os.path.isfile(cwl+'/PatientPathwaysTraining/DrugPathways') == True:
		print "Begining to modify pathways weightages... Please type !x when finished"
		target = open(cwl+'/PatientPathwaysTraining/DrugPathways', 'r')
		drug_specific_pathways = pickle.load(target)
		target.close()
		print "Selected pathways:"
		for component in drug_specific_pathways:
			print component, drug_specific_pathways[component]
		while True:
			key = raw_input('Please enter key to be modified: ')
			if key == "!x":
				break
			else:
				value = raw_input('Please enter new value to be appended: ')
				drug_specific_pathways[key] = value
		target = open(cwl+'/PatientPathwaysTraining/DrugPathways', 'w')
		target.truncate()
		pickle.dump(drug_specific_pathways,target)
		target.close()
	else:
		print "pathway score for drugs not created to modify weightages"
	print "Modification completed and stored"


def Calibrate_Response_Coefficients(givenlist):
	calib_patient_pathways = {}
	calib_patient_response_probability = {}
	predicted_response = {}
	actual_response = {}
	accuracy_by_pvalue = {}
	if os.path.isfile(cwl+'/PatientPathwaysTraining/DrugPathways') == True:
		if not os.path.exists(cwl+'/CalibrationData'):
			os.makedirs(cwl+'/CalibrationData')
		target0 = open(cwl+'/PatientPathwaysTraining/DrugPathways','r')
		drug_specific_pathways = pickle.load(target0)
		target0.close()
		for calib_patient in givenlist:
			calib_patient_name = str(calib_patient)
			calib_patient = Patient()
			calib_patient.response = str(calib_patient_response[calib_patient_name])
			actual_response[calib_patient_name] = calib_patient.response
			calib_patient.source_mutation(calib_patient_msource[calib_patient_name])
			calib_patient.generate_pathways_signature(calib_patient_name, cwl+'/UpdatedPathwaysScores', cwl+'/CalibrationData')
			target = open(cwl+'/CalibrationData/'+calib_patient_name,'r')
			calib_patient_pathways = pickle.load(target)
			target.close()
			calib_patient_response_probability[calib_patient_name] = 0.0
			n_calib_patient = len(calib_patient_pathways)
			n_drug = len(drug_specific_pathways)
			for paths in calib_patient_pathways:
				for selected_paths in drug_specific_pathways:
					if paths == selected_paths:
						calib_patient_response_probability[calib_patient_name] += float((2/n_calib_patient)*(2/n_drug)*(calib_patient_pathways[paths])*(drug_specific_pathways[selected_paths]))
			print calib_patient_response_probability[calib_patient_name]
		for pvalue in range(-10,991):
			positive_corelation_score = 0
			for calib_patient_name in calib_patient_response_probability:
				if calib_patient_response_probability[calib_patient_name] >= float(pvalue*0.001):
					predicted_response[calib_patient_name] = 'R'
				elif calib_patient_response_probability[calib_patient_name] < float(pvalue*0.001):
					predicted_response[calib_patient_name] = 'N'
			for calib_patient in predicted_response:
				if predicted_response[calib_patient] == actual_response[calib_patient]:
						positive_corelation_score += 1
			accuracy_by_pvalue[pvalue] = float(positive_corelation_score*100/len(calib_patient_response_probability))
		sorted_pvalue =  sorted(accuracy_by_pvalue, key=accuracy_by_pvalue.get)
		pout = sorted_pvalue[1000] #If the pvalue list length is modified, please update the last index number here
		print "Calibration co-relation: ", accuracy_by_pvalue[pout]
		return (pout*0.001)

def Test_Patient_Response(givenlist,giventh):
	test_patient_pathways = {}
	test_patient_response_probability = {}
	predicted_response = {}
	if os.path.isfile(cwl+'/PatientPathwaysTraining/DrugPathways') == True:
		if not os.path.exists(cwl+'/PredictionData'):
			os.makedirs(cwl+'/PredictionData')
		for test_patient in givenlist:
			# Opening drug pathways dict from a pickle with every iteration because, OptimizePathwayWeightage will update values with every patient
			target0 = open(cwl+'/PatientPathwaysTraining/DrugPathways','r')
			drug_specific_pathways = pickle.load(target0)
			target0.close()
			test_patient_name = str(test_patient)
			test_patient = Patient()
			test_patient.source_mutation(test_patient_msource[test_patient_name])
			test_patient.generate_pathways_signature(test_patient_name, cwl+'/UpdatedPathwaysScores', cwl+'/PredictionData')
			target = open(cwl+'/PredictionData/'+test_patient_name,'r')
			test_patient_pathways = pickle.load(target)
			target.close()
			test_patient_response_probability[test_patient_name] = 0.0
			print "Identified pathways in patient:", test_patient_pathways
			n_patient = len(test_patient_pathways)
			n_drug = len(drug_specific_pathways)
			for paths in test_patient_pathways:
				for selected_paths in drug_specific_pathways:
					if paths == selected_paths:
						test_patient_response_probability[test_patient_name] += float((2/n_patient)*(2/n_drug)*(test_patient_pathways[paths])*(drug_specific_pathways[selected_paths]))
			print test_patient_response_probability[test_patient_name]
			if test_patient_response_probability[test_patient_name] >= float(giventh):
				predicted_response[test_patient_name] = 'R'
			elif test_patient_response_probability[test_patient_name] < float(giventh):
				predicted_response[test_patient_name] = 'N'
			print "Predicted Response for Patient", test_patient_name, "is", predicted_response[test_patient_name]
			if update_mode == 0:
				continue
			else:
				OptimizePathwayWeightage(test_patient_pathways,drug_specific_pathways)
		targetn = open(cwl+'/PredictionData/prediction.log','w')
		pickle.dump(predicted_response,targetn)
		targetn.close()
		

def OptimizePathwayWeightage(adict, bdict):
	expected_result = int(raw_input('press [1] for match | [0] for mismatch '))
	if expected_result == 1:
		print "modifying values..."
		for component in adict:
			if component in bdict:
				bdict[component] = bdict[component] * 1.15
				print "Updated weightage for", component, "from", bdict[component] / 1.15, "to", bdict[component]
	elif expected_result == 0:
		print "modifying values..."
		for component in adict:
			if component in bdict:
				bdict[component] = bdict[component] / 1.15
				print "Updated weightage for", component, "from", bdict[component] * 1.15, "to", bdict[component]
	target0 = open(cwl+'/PatientPathwaysTraining/DrugPathways','w')
	pickle.dump(bdict,target0)
	target0.close()

def ValidatePredictions(adict):
	if os.path.isfile(cwl+'/PredictionData/prediction.log') == True:
		target0 = open(cwl+'/PredictionData/prediction.log','r')
		bdict = pickle.load(target0)
		target0.close()
		total_positive = 0
		total_negative = 0
		correct_positive = 0
		correct_negative = 0
		for component in adict:
			if adict[component] == 'R':
				total_positive += 1
			elif adict[component] == 'N':
				total_negative += 1
		target1 = open(cwl+'/PredictionData/prediction_match.log','w')
		for component in bdict:
			status = 'MisMatch'
			entry = None
			if bdict[component] == adict[component]:
				status = 'Match'
				if bdict[component] == 'R':
					correct_positive += 1
				elif bdict[component] == 'N':
					correct_negative += 1
			entry = str(component)+',Predicted :'+str(bdict[component])+',Actual :'+str(adict[component])+','+str(status)
			target1.write(entry)
			target1.write('\n')
		target1.close()
		print "Tabulated results..."
		print ""
		print "ACCURACY: ", ((correct_positive + correct_negative)/(total_positive + total_negative)*100), " % (", correct_positive + correct_negative, "/", total_positive + total_negative, ")"
		if total_positive > 0:
			print "PPV: ", (correct_positive/total_positive)*100, " % (", correct_positive, "/", total_positive, ")" 
		if total_negative > 0:
			print "NPV: ", (correct_negative/total_negative)*100, " % (", correct_negative, "/", total_negative, ")" 
	else:
		print "Prediction log not found. Please run prediction first..."
	
	

if __name__ == "__main__":
	if len(sys.argv) != 4:
		print "Usage: Pred.py <training-list> <calibration-list> <test.list>"
		# NZEC for chaining
		sys.exit(-1)
	script, trainlist, caliblist, testlist = sys.argv

	if not os.path.exists(cwl+'/UpdatedPathwaysScores'):
		os.makedirs(cwl+'/UpdatedPathwaysScores')
	if not os.path.exists(cwl+'/PatientPathwaysTraining'):
		os.makedirs(cwl+'/PatientPathwaysTraining')

	UpdatePathwayScores(cwl+'/UpdatedPathwaysScores')
	
	update_mode = int(raw_input('Update mode? [1/0] '))

	print "Reading training data..."

	patients = [patients.split(',')[0] for patients in open(trainlist)]
	response = [response.split(',')[1] for response in open(trainlist)]
	msource = [msource.split(',')[2] for msource in open(trainlist)]
	mosurce_split = []
	for elements in msource:
		mosurce_split.append(elements.rstrip('\n'))
	patient_response = dict(zip(patients, response))
	patient_msource = dict(zip(patients, mosurce_split))

	print "Finished reading training data"
	if os.path.isfile(cwl+'/PatientPathwaysTraining/DrugPathways') == False:

		Categorize_Drug_Response(patients)
	
		print "Would you like to modify Pathway weightages?"
		modify_status = raw_input('Y/N: ')
		if modify_status == 'Y' or modify_status == 'y':
			ForceWeightage()

		print "Reading data for Calibration..."

		calib_patients = [calib_patients.split(',')[0] for calib_patients in open(caliblist)]
		calib_response = [calib_response.split(',')[1] for calib_response in open(caliblist)]
		calib_msource = [calib_msource.split(',')[2] for calib_msource in open(caliblist)]
		calib_msource_split = []
		for elements in calib_msource:
			calib_msource_split.append(elements.rstrip('\n'))
		calib_patient_response = dict(zip(calib_patients, calib_response))
		calib_patient_msource = dict(zip(calib_patients, calib_msource_split))
		print "Finished reading calibration data"
		pout =	Calibrate_Response_Coefficients(calib_patients)
	else:
		print "Drug Specific Pathways data found. Proceeding to calibration..."
		print "Reading data for Calibration..."

		calib_patients = [calib_patients.split(',')[0] for calib_patients in open(caliblist)]
		calib_response = [calib_response.split(',')[1] for calib_response in open(caliblist)]
		calib_msource = [calib_msource.split(',')[2] for calib_msource in open(caliblist)]
		calib_msource_split = []
		for elements in calib_msource:
			calib_msource_split.append(elements.rstrip('\n'))
		calib_patient_response = dict(zip(calib_patients, calib_response))
		calib_patient_msource = dict(zip(calib_patients, calib_msource_split))
		print "Finished reading calibration data"

		pout =	Calibrate_Response_Coefficients(calib_patients)
	print "Initiating Prediction..."
	print "Reading test data..."
	test_patients = [test_patients.split(',')[0] for test_patients in open(testlist)]
	test_msource = [test_msource.split(',')[2] for test_msource in open(testlist)]
	test_msource_split = []
	for elements in test_msource:
		test_msource_split.append(elements.rstrip('\n'))
	test_patient_msource = dict(zip(test_patients, test_msource_split))
	print "Finished reading test data"
	print "Current calibrated Pout value is", pout
	y = raw_input("Would you like to update a new threshold value?[Y/N]: ")
	if y == 'Y' or y =='y':
		th = raw_input('Enter threshold value: ')
	else:
		th = pout
	Test_Patient_Response(test_patients,th)
	print "Finished Prediction for test set. Do you want to validate predictions?"
	nod = raw_input('[Y/N]: ')
	if nod == 'Y' or nod == 'y':
		vpatients = [vpatients.split(',')[0] for vpatients in open(testlist)]
		vresponse = [vresponse.split(',')[1] for vresponse in open(testlist)]
		vpatient_response = dict(zip(vpatients, vresponse))

		ValidatePredictions(vpatient_response)
	
