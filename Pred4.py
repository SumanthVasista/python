#Python program to cluster mutations and pathways with drug response
from __future__ import division
import sys
import pickle
import os
import math as mt
import numbers
import numpy

cwl = os.path.dirname(os.path.realpath(__file__))

#Defining global variales
pathways_initialized = {}
pathways_causal_depth = {}
pathways_selected_for_patients = {}
pathway_causal_depth = {}
holder_patient = {}

class Patient:
	def __init__(self):
		self.type = None
		self.id = None
		self.gender = None
		self.disease = None
		self.response = None
		self.mutations = {}
		self.pathways = {}
		self.pathwaysmap = {}
		list_of_pathways = [list_of_pathways.rstrip('\n') for list_of_pathways in open(cwl+'/UpdatedPathwaysScores/list')]
		for component in list_of_pathways:
			pathway_name = str(component)
			self.pathways[pathway_name] = 0

	def source_mutation(self,a_list):
		gene = [gene.split(',')[0] for gene in open(a_list)]
		mut = [mut.split(',')[1] for mut in open(a_list)]
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
				target_gof = open(str(source)+'/'+str(pathway)+'-GOF','r')
				target_lof = open(str(source)+'/'+str(pathway)+'-LOF','r')
				expanded_pathway_gof = pickle.load(target_gof)
				expanded_pathway_lof = pickle.load(target_lof)
				target_gof.close()
				target_lof.close()
				if mutation in expanded_pathway_gof:
					if self.mutations[mutation] == 'GOF':
						self.pathways[pathway] += expanded_pathway_gof[mutation]
					elif self.mutations[mutation] == 'LOF':
						self.pathways[pathway] += -1*expanded_pathway_lof[mutation]
		target = open(dest+'/'+str(id), 'w')
		pathways_selected_for_patients = dict((k, v) for k, v in self.pathways.items() if v >= 0.7 or v <= -0.7)
		if len(pathways_selected_for_patients) < 1:
			pickle.dump(self.pathways,target)
			print "Patient does not express any strong pathways signatures... Proceeding..."
		else:
			if len(pathways_selected_for_patients) > 10:
				print "Trimming..."
				pathways_reselected_for_patients = dict((k, v) for k, v in pathways_selected_for_patients.items() if v >= 1.25 or v <= -1.25)
				if len(pathways_reselected_for_patients) > 0:
					pickle.dump(pathways_reselected_for_patients,target)
				else:
					pickle.dump(pathways_selected_for_patients,target)
			else:
				pickle.dump(pathways_selected_for_patients,target)
		self.pathwaysmap = pathways_selected_for_patients
		pickle.dump(self.pathways,target)
		target.close()
		target2 = open(dest+'/list','w')
		target2.write(str(id))
		target2.write('\n')
		target2.close()
		print "Finished generating patient specific pathways data for", id


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
				genes_in_pathways = [genes_in_pathways.split(',')[0] for genes_in_pathways in open(loc_pathways+'/'+listed_paths)]
				gof_score_in_pathways = [gof_score_in_pathways.split(',')[1] for gof_score_in_pathways in open(loc_pathways+'/'+listed_paths)]
				try:
					lof_score_in_pathways = [lof_score_in_pathways.split(',')[2] for lof_score_in_pathways in open(loc_pathways+'/'+listed_paths)]
				except:
					print sys.exc_info()
					lof_score_in_pathways = [lof_score_in_pathways.split(',')[1] for lof_score_in_pathways in open(loc_pathways+'/'+listed_paths)]
				gof_score_in_pathways_strip = []
				lof_score_in_pathways_strip = []
				for elements in gof_score_in_pathways:
					gof_score_in_pathways_strip.append(float(elements))
				for elements in lof_score_in_pathways:
					lof_score_in_pathways_strip.append(float(elements))
				gof_listed_paths = dict(zip(genes_in_pathways, gof_score_in_pathways_strip))
				lof_listed_paths = dict(zip(genes_in_pathways, lof_score_in_pathways_strip))
				target = open(default_location+'/'+pathname+'-GOF','w')
				pickle.dump(gof_listed_paths,target)
				target.close()
				lof_listed_paths = dict(zip(genes_in_pathways, lof_score_in_pathways_strip))
				target = open(default_location+'/'+pathname+'-LOF','w')
				pickle.dump(lof_listed_paths,target)
				target.close()
				target2.write(pathname)
				target2.write('\n')
			target2.close()
			print "Finished creating Pathway Relationship Score data files. Updated in location", default_location
		else:
			print "No list found in", loc_pathways, "to update pathways"
	else:
		print "Pathway Relationship Score data found. Proceeding..."

def Categorize_Drug_Response(givenlist):
	adaptive_threshold = []
	drug_specific_pathways = {}
	list_of_pathways = [list_of_pathways.rstrip('\n') for list_of_pathways in open(cwl+'/UpdatedPathwaysScores/list')]
	for component in list_of_pathways:
		pathway_name = str(component)
		drug_specific_pathways[pathway_name] = 0
	target1 = open(cwl+'/PatientPathwaysTraining/.kbplist','w')
	for patient in givenlist:
		patient_name = str(patient)
		patient = Patient()
		patient.id = patient_name
		patient.response = str(patient_response[patient_name])
		patient.source_mutation(patient_msource[patient_name])
		patient.generate_pathways_signature(patient_name, cwl+'/UpdatedPathwaysScores', cwl+'/PatientPathwaysTraining')
		target = open(cwl+'/PatientPathwaysTraining/'+patient_name,'r')
		patient_pathways = {}
		patient_pathways = pickle.load(target)
		target.close()
		holder_patient[patient_name] = patient
		target0 = open(cwl+'/PatientPathwaysTraining/.'+patient_name+'.clf','w')
		target1.write(patient_name)
		target1.write('\n')
		pickle.dump(patient,target0)
		target0.close()
		adaptive_threshold.append(len(patient_pathways))
		for pathway in patient_pathways:
			if patient.response == 'R':
				drug_specific_pathways[pathway] += 0.2*patient_pathways[pathway]
			elif patient.response == 'N':
				drug_specific_pathways[pathway] -= 0.2*patient_pathways[pathway]
	target1.close()
	print "ADAPTING"
	noise = sum(adaptive_threshold)/len(adaptive_threshold)
	a = float(0.035*noise)
	b = float(-0.035*noise)
	print a, b
	pathways_selected_for_drug = dict((k, v) for k, v in drug_specific_pathways.items() if v >= a or v <= b)
	print pathways_selected_for_drug
	target = open(cwl+'/PatientPathwaysTraining/DrugPathways','w')
	pickle.dump(pathways_selected_for_drug,target)
	target.close()
	print "Finished clustering pathways to drug response"

def cluster_patient(apatient):
	search_space = []
	search_space_reduced = []
	if len(holder_patient) > 5:
#	step 1
#		for query in apatient.pathways:
#			for patient in holder_patient:
#				if query in patient.pathways and patient.pathways[query]*apatient.pathways[query] > 0:
#					search_space.append(patient)
#			print len(search_space), "similar patients found to begin"
#			marker = str(query)
#			break
		for element1 in holder_patient:
			search_space.append(holder_patient[element1])
#	step 2
		while len(search_space) != 0:
			for query in apatient.pathwaysmap:
				for patient in search_space:
					if query not in patient.pathwaysmap or patient.pathwaysmap[query]*apatient.pathwaysmap[query] <= 0:
						search_space_reduced = []
						for element2 in search_space:
							search_space_reduced.append(element2)
						search_space.pop(search_space.index(patient))
					if len(search_space) == 0:
						break
			break
	print "Closest relatives found for: ", apatient.id
	similarity = []
	for element in search_space_reduced:
		if element.response == 'R':
			score = +1
		elif element.response == 'N':
			score = -1
		elif element.response == None:
			score = 0
		if abs(len(element.pathwaysmap)/len(apatient.pathwaysmap)) >= 1:
			update_score = float(score*len(apatient.pathwaysmap)/(0.01+len(element.pathwaysmap)))
			if abs(update_score) > 0.3:
				similarity.append(update_score)
		else:
			update_score = float(score*len(element.pathwaysmap)/(0.01+len(apatient.pathwaysmap)))
			if abs(update_score) > 0.3:
				similarity.append(update_score)
		print element.id, update_score, element.response
	response_score = sum(similarity)/(0.01+len(similarity))
	return response_score


def Calibrate_Response_Coefficients(givenlist):
	calib_patient_pathways = {}
	calib_patient_response_probability = {}
	calib_patient_causal_depth = {}
	predicted_response = {}
	actual_response = {}
	accuracy_by_pvalue = {}
	if os.path.isfile(cwl+'/PatientPathwaysTraining/DrugPathways') == True:
		if not os.path.exists(cwl+'/CalibrationData'):
			os.makedirs(cwl+'/CalibrationData')
		target0 = open(cwl+'/PatientPathwaysTraining/DrugPathways','r')
		drug_specific_pathways = pickle.load(target0)
		target0.close()
		#target0 = open(cwl+'/PatientPathwaysTraining/PathwayAntagony', 'r')
		#pathways_causal_depth = pickle.load(target0)
		#target0.close()
		kbplist = [kbplist.rstrip('\n') for kbplist in open(cwl+'/PatientPathwaysTraining/.kbplist')]
		if len(holder_patient) == 0:
			for element in kbplist:
				target0 = open(cwl+'/PatientPathwaysTraining/.'+element+'.clf','r')
				p = pickle.load(target0)
				holder_patient[element] = p
				p = None
				target0.close()
			print "captured knowledge of", len(holder_patient), "patients"
		for calib_patient in givenlist:
			calib_patient_name = str(calib_patient)
			calib_patient = Patient()
			calib_patient.id = calib_patient_name
			calib_patient.response = str(calib_patient_response[calib_patient_name])
			actual_response[calib_patient_name] = calib_patient.response
			calib_patient.source_mutation(calib_patient_msource[calib_patient_name])
			calib_patient.generate_pathways_signature(calib_patient_name, cwl+'/UpdatedPathwaysScores', cwl+'/CalibrationData')
			target = open(cwl+'/CalibrationData/'+calib_patient_name,'r')
			calib_patient_pathways = pickle.load(target)
			target.close()
			calib_patient_response_probability[calib_patient_name] = 0.0
			calib_patient_causal_depth[calib_patient_name] = 0
			n_calib_patient = len(calib_patient_pathways)
			n_drug = len(drug_specific_pathways)
			for paths in calib_patient_pathways:
				for selected_paths in drug_specific_pathways:
					if paths == selected_paths:
						if calib_patient_pathways[selected_paths] != 0:
							calib_patient_response_probability[calib_patient_name] += float((5/(n_calib_patient*n_calib_patient))*(5/n_drug)*(drug_specific_pathways[selected_paths])*(calib_patient_pathways[selected_paths]))
						else:
							calib_patient_response_probability[calib_patient_name] += 0
			calib_patient_causal_depth[calib_patient_name] = cluster_patient(calib_patient)
			print calib_patient_response_probability[calib_patient_name], calib_patient_causal_depth[calib_patient_name]
			if calib_patient_name not in holder_patient:
				holder_patient[calib_patient_name] = calib_patient
			if calib_patient_name not in kbplist:
				target0 = open(cwl+'/PatientPathwaysTraining/.'+calib_patient_name+'.clf','w')
				target1 = open(cwl+'/PatientPathwaysTraining/.kbplist','a')
				pickle.dump(calib_patient,target0)
				target1.write(calib_patient_name)
				target1.write('\n')
				target1.close()
				target0.close()
		for pvalue in range(-50,50):
			positive_corelation_score = 0
			for calib_patient_name in calib_patient_response_probability:
				if calib_patient_response_probability[calib_patient_name] >= float(pvalue*0.05):
					predicted_response[calib_patient_name] = 'R'
				elif calib_patient_response_probability[calib_patient_name] < float(pvalue*0.05):
					predicted_response[calib_patient_name] = 'N'
			for calib_patient in predicted_response:
				if predicted_response[calib_patient] == actual_response[calib_patient]:
						positive_corelation_score += 1
			accuracy_by_pvalue[pvalue] = float(positive_corelation_score*100/len(calib_patient_response_probability))
		sorted_pvalue =  sorted(accuracy_by_pvalue, key=accuracy_by_pvalue.get)
#		for elem in sorted_pvalue:
#			print elem, ":", accuracy_by_pvalue[elem]
		pout = sorted_pvalue[99] #If the pvalue list length is modified, please update the last index number here
		print "Calibration co-relation: ", accuracy_by_pvalue[pout]
		PCR_list = []
		NCR_list = []
		PAR_list = []
		NAR_list = []
		for calib_patient in actual_response:
			if actual_response[calib_patient] == 'R':
				PCR_list.append(float(calib_patient_response_probability[calib_patient]))
				PAR_list.append(float(calib_patient_causal_depth[calib_patient]))
			elif actual_response[calib_patient] == 'N':
				NCR_list.append(float(calib_patient_response_probability[calib_patient]))
				NAR_list.append(float(calib_patient_causal_depth[calib_patient]))
		PCR = numpy.mean(PCR_list)
		NCR = numpy.mean(NCR_list)
		PCL = numpy.std(PCR_list)
		NCL = numpy.std(NCR_list)
		PCR_list_update = []
		NCR_list_update = []
		for component in PCR_list:
			if component > float(PCR - PCL):
				PCR_list_update.append(component)
		for component in NCR_list:
			if component < float(NCR + NCL):
				NCR_list_update.append(component)
		PAR = sum(PAR_list)/len(PAR_list)
		NAR = sum(NAR_list)/len(NAR_list)
		o_threshold = (pout*0.05)
		x_threshold = float((min(PCR_list_update) + max(NCR_list_update))/2)
		y_threshold = float((PAR+NAR)/2)
		print "Positive Clustered region(x):", PCR
		print "Negative Clustered region(x):", NCR
		print "Clustering distance(x):", abs(min(PCR_list) - max(NCR_list))
		print "Optimal threshold(x):", x_threshold
		target = open(cwl+'/CalibrationData/pout_static','w')
		Pout_to_store = [o_threshold, x_threshold, y_threshold]
		pickle.dump(Pout_to_store, target)
		target.close()
		return o_threshold ,x_threshold, y_threshold

def Test_Patient_Response(givenlist,g_x,g_y):
	test_patient_pathways = {}
	test_patient_response_probability = {}
	test_patient_causal_depth = {}
	predicted_response = {}
	if os.path.isfile(cwl+'/PatientPathwaysTraining/DrugPathways') == True:
		if not os.path.exists(cwl+'/PredictionData'):
			os.makedirs(cwl+'/PredictionData')
		tx = open(cwl+'/PredictionData/cluster.dat','w')
		#target0 = open(cwl+'/PatientPathwaysTraining/PathwayAntagony', 'r')
		#pathways_causal_depth = pickle.load(target0)
		#target0.close()
		kbplist = [kbplist.rstrip('\n') for kbplist in open(cwl+'/PatientPathwaysTraining/.kbplist')]
		if len(holder_patient) == 0:
			for element in kbplist:
				target0 = open(cwl+'/PatientPathwaysTraining/.'+element+'.clf','r')
				p = pickle.load(target0)
				holder_patient[element] = p
				p = None
				target0.close()
			print "captured knowledge of", len(holder_patient), "patients"
		for test_patient in givenlist:
			# Opening drug pathways dict from a pickle with every iteration because, OptimizePathwayWeightage will update values with every patient
			target0 = open(cwl+'/PatientPathwaysTraining/DrugPathways','r')
			drug_specific_pathways = pickle.load(target0)
			target0.close()
			test_patient_name = str(test_patient)
			test_patient = Patient()
			test_patient.id = test_patient_name
			test_patient.source_mutation(test_patient_msource[test_patient_name])
			test_patient.generate_pathways_signature(test_patient_name, cwl+'/UpdatedPathwaysScores', cwl+'/PredictionData')
			test_patient.response = test_patient_response[test_patient_name]
			if test_patient.response == '':
				 test_patient.response = None
			target = open(cwl+'/PredictionData/'+test_patient_name,'r')
			test_patient_pathways = pickle.load(target)
			target.close()
			test_patient_response_probability[test_patient_name] = 0.0
			test_patient_causal_depth[test_patient_name] = 0
			print "Identified pathways in patient:", test_patient_pathways
			n_patient = len(test_patient_pathways)
			n_drug = len(drug_specific_pathways)
			for paths in test_patient_pathways:
				for selected_paths in drug_specific_pathways:
					if paths == selected_paths:
						test_patient_response_probability[test_patient_name] += float((5/(n_patient*n_patient))*(5/n_drug)*(drug_specific_pathways[selected_paths])*(test_patient_pathways[selected_paths]))
			print test_patient_response_probability[test_patient_name], test_patient_causal_depth[test_patient_name]
			if test_patient_response_probability[test_patient_name] > float(g_x):
				predicted_response[test_patient_name] = 'R'
			elif test_patient_response_probability[test_patient_name] <= float(g_x):
				predicted_response[test_patient_name] = 'N'
			test_patient_causal_depth[test_patient_name] = cluster_patient(test_patient)
			if test_patient_name not in holder_patient:
				holder_patient[test_patient_name] = test_patient
			if test_patient.response != None:
				if test_patient.response == 'R':
					j = 1
				else:
					j = 0
			else:
				j = ''
			string = str(test_patient_name)+' '+str(float((test_patient_response_probability[test_patient_name])) - g_x)+' '+str(float(test_patient_causal_depth[test_patient_name]))+' '+str(j)
			tx.write(string)
			tx.write('\n')
			print "Predicted Response for Patient", test_patient_name, "is", predicted_response[test_patient_name]
			if update_mode == 0:
				continue
			elif update_mode == 1:
				diff = test_patient_response_probability[test_patient_name] - g_x
				if test_patient.response != None:
					if predicted_response[test_patient_name] == test_patient.response:
						OptimizePathwayWeightage_pa(test_patient_pathways,drug_specific_pathways,1,diff,test_patient.response)
					else:
						OptimizePathwayWeightage_pa(test_patient_pathways,drug_specific_pathways,-1,diff,test_patient.response)
				else:
					pass
			elif update_mode == 2:
				diff = test_patient_response_probability[test_patient_name] - g_x
				if test_patient.response != None:
					if predicted_response[test_patient_name] == test_patient.response:
						OptimizePathwayWeightage_hyp(test_patient_pathways,drug_specific_pathways,1,diff,test_patient.response)
					else:
						OptimizePathwayWeightage_hyp(test_patient_pathways,drug_specific_pathways,-1,diff,test_patient.response)
				else:
					pass
			if test_patient_name not in kbplist:
				target0 = open(cwl+'/PatientPathwaysTraining/.'+test_patient_name+'.clf','w')
				target1 = open(cwl+'/PatientPathwaysTraining/.kbplist','a')
				pickle.dump(test_patient,target0)
				target1.write(test_patient_name)
				target1.write('\n')
				target1.close()
				target0.close()
			print "\n"
		tx.close()
		targetn = open(cwl+'/PredictionData/prediction.log','w')
		pickle.dump(predicted_response,targetn)
		targetn.close()


def OptimizePathwayWeightage_pa(adict, bdict, beta, skew, actual):
	if beta == 1:
		pass

	elif beta == -1:
		print "modifying values (P/A)"
		for component in adict:
			if component in bdict:
				if bdict[component] >= -0.1 and bdict[component] < 0:
					bdict[component] = 0.1
				elif bdict[component] <= 0.1 and bdict[component] > 0:
					bdict[component] = -0.1
			if component in bdict:
				if skew < 0:
					if actual == 'R':
						if adict[component] > 0 and bdict[component] > 0:
							bdict[component] = bdict[component] * (1.25* float(1 + abs(skew)))
						elif adict[component] > 0 and bdict[component] < 0:
							bdict[component] = bdict[component] / (1.1 * float(1 + abs(skew)))
						elif adict[component] < 0 and bdict[component] > 0:
							bdict[component] = bdict[component] / (1.1 * float(1 + abs(skew)))
						elif adict[component] < 0 and bdict[component] < 0:
							bdict[component] = bdict[component] * (1.25* float(1 + abs(skew)))
					elif actual == 'N':
						pass
				elif skew > 0:
					if actual == 'R':
						pass
					elif actual == 'N':
						if adict[component] > 0 and bdict[component] > 0:
							bdict[component] = bdict[component] / (1.1 * float(1 + abs(skew)))
						elif adict[component] > 0 and bdict[component] < 0:
							bdict[component] = bdict[component] * (1.25 * float(1 + abs(skew)))
						elif adict[component] < 0 and bdict[component] > 0:
							bdict[component] = bdict[component] * (1.25 * float(1 + abs(skew)))
						elif adict[component] < 0 and bdict[component] < 0:
							bdict[component] = bdict[component] / (1.1 * float(1 + abs(skew)))						
	target0 = open(cwl+'/PatientPathwaysTraining/DrugPathways','w')
	pickle.dump(bdict,target0)
	target0.close()

def OptimizePathwayWeightage_hyp(adict, bdict, beta, skew, actual):
	if beta == 1:
		print "modifying values(Hyperbolic - P)"
		for component in adict:
			if component in bdict:
				if skew < 0:
					if actual == 'R':
						pass
					elif actual == 'N':
						if adict[component] > 0 and bdict[component] > 0:
							bdict[component] = bdict[component] * (1.25 * float(1 + abs(skew)))
						elif adict[component] > 0 and bdict[component] < 0:
							bdict[component] = bdict[component] / (1.25 * float(1 + abs(skew)))
						elif adict[component] < 0 and bdict[component] > 0:
							bdict[component] = bdict[component] / (1.25 * float(1 + abs(skew)))
						elif adict[component] < 0 and bdict[component] < 0:
							bdict[component] = bdict[component] * (1.25 * float(1 + abs(skew)))
				elif skew > 0:
					if actual == 'R':
						if adict[component] > 0 and bdict[component] > 0:
							bdict[component] = bdict[component] / (1.25 * float(1 + abs(skew)))
						elif adict[component] > 0 and bdict[component] < 0:
							bdict[component] = bdict[component] * (1.25 * float(1 + abs(skew)))
						elif adict[component] < 0 and bdict[component] > 0:
							bdict[component] = bdict[component] * (1.25 * float(1 + abs(skew)))
						elif adict[component] < 0 and bdict[component] < 0:
							bdict[component] = bdict[component] / (1.25 * float(1 + abs(skew)))	

					elif actual == 'N':
						pass
	elif beta == -1:
		print "modifying values(Hyperbolic - N)"
		for component in adict:
			if component in bdict:
				if bdict[component] > -0.1 and bdict[component] < 0:
					bdict[component] = 0.1
				elif bdict[component] < 0.1 and bdict[component] > 0:
					bdict[component] = -0.1
			if component in bdict:
				if skew < 0:
					if actual == 'R':
						if adict[component] > 0 and bdict[component] > 0:
							bdict[component] = bdict[component] * (1.25 * float(1 + abs(skew)))
						elif adict[component] > 0 and bdict[component] < 0:
							bdict[component] = bdict[component] / (1.1 * float(1 + abs(skew)))
						elif adict[component] < 0 and bdict[component] > 0:
							bdict[component] = bdict[component] / (1.1 * float(1 + abs(skew)))
						elif adict[component] < 0 and bdict[component] < 0:
							bdict[component] = bdict[component] * (1.25 * float(1 + abs(skew)))
					elif actual == 'N':
						pass
				elif skew > 0:
					if actual == 'R':
						pass
					elif actual == 'N':
						if adict[component] > 0 and bdict[component] > 0:
							bdict[component] = bdict[component] / (1.1 * float(1 + abs(skew)))
						elif adict[component] > 0 and bdict[component] < 0:
							bdict[component] = bdict[component] * (1.25 * float(1 + abs(skew)))
						elif adict[component] < 0 and bdict[component] > 0:
							bdict[component] = bdict[component] * (1.25 * float(1 + abs(skew)))
						elif adict[component] < 0 and bdict[component] < 0:
							bdict[component] = bdict[component] / (1.1 * float(1 + abs(skew)))						
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
	if len(sys.argv) != 5:
		print "Usage: Pred.py <training-list> <calibration-list> <test.list> <update-mode[0/1/2]>"
		print "[ Update 0 : Prediction Mode ] [ Update 1 : P/A Classification Mode ] [ Update 2 : Hyperbolic Classification Mode ]"
		# NZEC for chaining
		sys.exit(-1)
	script, trainlist, caliblist, testlist, upmode = sys.argv

	if not os.path.exists(cwl+'/UpdatedPathwaysScores'):
		os.makedirs(cwl+'/UpdatedPathwaysScores')
	if not os.path.exists(cwl+'/PatientPathwaysTraining'):
		os.makedirs(cwl+'/PatientPathwaysTraining')

	UpdatePathwayScores(cwl+'/UpdatedPathwaysScores')
	
	update_mode = int(upmode)
	if isinstance(update_mode, int) == False:
		sys.exit(-1)
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
	
#		print "Would you like to modify Pathway weightages?"
#		modify_status = raw_input('Y/N: ')
#		if modify_status == 'Y' or modify_status == 'y':
#			ForceWeightage()

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
		Pout =	Calibrate_Response_Coefficients(calib_patients)
	else:
		if os.path.isfile(cwl+'/CalibrationData/pout_static') == False:
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

			Pout = Calibrate_Response_Coefficients(calib_patients)
		else:
			print "Calibration point found. Proceeding to prediction"
			target = open(cwl+'/CalibrationData/pout_static','r')
			Pout = pickle.load(target)
			target.close()
	print "Initiating Prediction..."
	print "Reading test data..."
	test_patients = [test_patients.split(',')[0] for test_patients in open(testlist)]
	test_response = [test_patients_repsonse.split(',')[1] for test_patients_repsonse in open(testlist)]
	test_msource = [test_msource.split(',')[2] for test_msource in open(testlist)]
	test_msource_split = []
	for elements in test_msource:
		test_msource_split.append(elements.rstrip('\n'))
	test_patient_msource = dict(zip(test_patients, test_msource_split))
	test_patient_response = dict(zip(test_patients,test_response))
	print "Finished reading test data"
	print "Current calibrated Pout value is", Pout[0]
	print "Current optimal Popt value is", Pout[1]
#	y = raw_input("Would you like to update a new threshold value?[Y/N]: ")
	if update_mode == 1:
		th_x = Pout[1]
		th_y = Pout[2]
		print "Selected OCV"
	else:
		th_x = Pout[1]
		th_y = Pout[2]
		print "Selected Highest co-relation value"
	Test_Patient_Response(test_patients,th_x,th_y)
	print "Finished Prediction for test set. Do you want to validate predictions?"
	print "Results present in - ./PredictionData"
	print "~~END~~"
