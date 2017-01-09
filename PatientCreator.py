#Python
from __future__ import division
import os
import sys
import math as mt
import random as rn
import numpy as np
import pickle

cwl = os.path.dirname(os.path.realpath(__file__))
patients_in_pipeline = []

class Mutation:
	def __init__(self):
		self.tag = None
		self.nature = None
		self.indication_score = None
		self.sis = {}
		self.collateral = {}

	def source_species_interations_score(self, a_picklefile):
		target = open(a_picklefile,'r')
		self.sis = pickle.load(target)
		target.close()

	def pick_collateral_components(self):
		for entries in self.sis:
			if self.sis[entries] >= 0.1 or self.sis[entries] < 0.1:
				self.collateral[entries] = self.sis[entries]


	def dominance(self):
		if self.tag == 'TP' or self.tag == 'NEUTRAL':
			if self.nature == 'GOF':
				return 1
		elif self.tag == 'TS' or self.tag == 'NEUTRAL':
			if self.nature == 'LOF':
				return 1
		else:
			return -1

class Patient:
	def __init__(self):
		self.patientid = None
		self.indication = None
		self.mutations = {}
		self.biomarkers = {}
		self.disease = {}

	def source_mutation_components(self, mutation_file):
		gene = [gene.split(',')[0] for gene in open(mutation_file)]
		mut = [mut.split(',')[1] for mut in open(mutation_file)]
		disease = dict(zip(gene, mut))
		for component in gene:
			if disease[component] == 'OE' or disease[component] == 'OE\n':
				self.mutations[component] = 'GOF'
			elif disease[component] == 'KD' or disease[component] == 'KD\n':
				self.mutations[component] = 'LOF'
		print "Finished sourcing mutations for patient"

	def create_pseudo_id(self):
		value = rn.randint(1000,99999)
		while value in patients_in_pipeline:
			value = rn.randint(1000,99999)
		self.patientid = value
		patients_in_pipeline.append(value)

'''Global Functions'''

def source_global_indices(a_listfile): #F1
	target = open(a_listfile,'r')
	lt = []
	for entries in target.read().splitlines():
		lt.append(str(entries))
	if os.path.isfile(cwl+'/dependecies') == True:
		target2 = open(cwl+'/dependecies/global_indicies.p','w')
		pickle.dump(lt,target2)
		target2.close()
	else:
		print "/dependencies directory not found, redirecting output to", cwl
		target2 = open(cwl+'/global_indicies.p','w')
		pickle.dump(lt,target2)
		target2.close()
	print "Finished sourcing Indices"
	target.close()

def source_critical_components(a_listfile): #F2
	target = open(a_listfile,'r')
	lt = []
	for entries in target.read().splitlines():
		lt.append(str(entries))
	if os.path.isfile(cwl+'/dependecies') == True:
		target2 = open(cwl+'/dependecies/global_critical_components.p','w')
		pickle.dump(lt,target2)
		target2.close()
	else:
		print "/dependencies directory not found, redirecting output to", cwl
		target2 = open(cwl+'/global_critical_components.p','w')
		pickle.dump(lt,target2)
		target2.close()
	print "Finished sourcing critical components"
	target.close()

def check_dependecies(): #F0
	print "Current Working Directory :", cwl
	print "Checking Dependecies..."
	proceed_signal = 1
	if os.path.isfile(cwl+'/dependecies') == False:
		print "Dependecies directory not found. Check for dependecies directory, update the requirements and try again.."
		proceed_signal = 0
	else:
		print "Dependecies directory found. Checking files..."
		if os.path.isfile(cwl+'/dependecies/global_indicies.p') == True:
			print "global Indicies list found"
		else:
			print "global Indicies list not found. Creating new list from scratch..."
			if os.path.isfile(cwl+'/dependecies/backup/global_indicies.txt') == True:
				source_global_indicies(cwl+'/dependecies/global_indicies.txt') #f1
			else:
				print "No files to source from. Exiting..."
				proceed_signal = 0
				exit(0)
		if os.path.isfile(cwl+'/dependecies/global_critical_components.p') == True:
			print "global critical_components list found"
		else:
			print "critical_components list not found. Creating new list from scratch..."
			if os.path.isfile(cwl+'/dependecies/backup/global_critical_components.txt') == True:
				source_critical_components(cwl+'/dependecies/global_critical_components.txt') #f2
			else:
				print "No files to source from. Exiting..."
				proceed_signal = 0
				exit(0)
						#Update dependecies list hereafter...
		if os.path.isfile(cwl+'/dependecies/gene_matrix/list') == True:
			print "gene matrix found"
		else:
			print "gene matrix not found. Please create mutation curve before initiation patients"
			proceed_signal = 0
			exit(0)
		if os.path.isfile(cwl+'/dependecies/targetnodedata.p') == True:
			print "target-node data found"
		else:
			print "target-node not found. Creating new list from scratch..."
			if os.path.isfile(cwl+'/dependecies/backup/targetnodedata.txt') == True:
				source_targetnodes(cwl+'/dependecies/backup/targetnodedata.txt') #f2
			else:
				print "No files to source from. Exiting..."
				proceed_signal = 0
				exit(0)
		return proceed_signal
