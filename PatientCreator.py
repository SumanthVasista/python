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
		self.mutation_type = None
		self.mutation_effect = None
		self.indication_score = None
		self.anchor = None
		self.sis = {}
		self.collateral = {}

	def source_indication_score(self,gene_name,code):
		if os.path.isfile(cwl+'/dependencies/'+code+'_indication_matrix.p') == True:
			target = open(cwl+'/dependencies/'+code+'_indication_matrix.p','r')
			matrix = pickle.load(target)
			target.close()
			return matrix[gene_name]
		else:
			return float(1)

	def source_species_interactive_score(self, a_picklefile):
		target = open(a_picklefile,'r')
		self.sis = pickle.load(target)
		target.close()

	def pick_collateral_components(self):
		for entries in self.sis:
			if self.sis[entries] >= 0.1 or self.sis[entries] < 0.1:
				self.collateral[entries] = self.sis[entries]
class Patient:
	def __init__(self):
		self.patientid = None
		self.indication = None 
		self.mutations = {}
		self.biomarkers = {}
		self.disease = {}
		self.mutation_mtype = {}
		self.gene_anchor = {}

	def source_mutation_components(self, mutation_file):
		gene = [gene.split(',')[0].rstrip('\n') for gene in open(mutation_file)]
		mut = [mut.split(',')[1].rstrip('\n') for mut in open(mutation_file)]
		mtype = [mtype.split(',')[2].rstrip('\n') for mtype in open(mutation_file)]
		anchor = [anchor.split(',')[3].rstrip('\n') for anchor in open(mutation_file)]
		disease = dict(zip(gene, mut))
		self.mutation_mtype = dict(zip(gene, mtype))
		self.gene_anchor = dict(zip(gene, anchor))
		for element in self.gene_anchor:
			if self.gene_anchor[element] == 'auto':
				self.gene_anchor[element] = None
		for component in disease:
			if disease[component] == 'OE' and self.mutation_mtype[component] == 'CNV':
				self.mutations[component] = 'CNVGOF'
			elif disease[component] == 'OE' and self.mutation_mtype[component] == 'MUT':
				self.mutations[component] = 'MUTGOF'
			elif disease[component] == 'KD' and self.mutation_mtype[component] == 'CNV':
				self.mutations[component] = 'CNVLOF'
			elif disease[component] == 'KD' and self.mutation_mtype[component] == 'MUT':
				self.mutations[component] = 'MUTLOF'
		print "Finished sourcing mutations for patient"

	def create_pseudo_id(self):
		value = rn.randint(1000,99999)
		while value in patients_in_pipeline:
			value = rn.randint(1000,99999)
		self.patientid = value
		patients_in_pipeline.append(value)

## Global Functions
def source_global_indicies(a_listfile): #F1
	target = open(a_listfile,'r')
	lt = []
	for entries in target.read().splitlines():
		lt.append(str(entries))
	if os.path.exists(cwl+'/dependencies') == True:
		target2 = open(cwl+'/dependencies/global_indicies.p','w')
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
	if os.path.exists(cwl+'/dependencies') == True:
		target2 = open(cwl+'/dependencies/global_critical_components.p','w')
		pickle.dump(lt,target2)
		target2.close()
	else:
		print "/dependencies directory not found, redirecting output to", cwl
		target2 = open(cwl+'/global_critical_components.p','w')
		pickle.dump(lt,target2)
		target2.close()
	print "Finished sourcing Critical Components"
	target.close()

def source_targetnodes(a_listfile): #F3
	gene = [gene.split(',')[0] for gene in open(a_listfile)]
	tg = [tg.split(',')[1].rstrip('\n') for tg in open(a_listfile)]
	key_list = dict(zip(tg, gene))
	if os.path.exists(cwl+'/dependencies') == True:
		target = open(cwl+'/dependencies/targetnodedata.p','w')
		pickle.dump(key_list,target)
	else:
		print "/dependencies directory not found, redirecting output to", cwl
		target = open(cwl+'/targetnodedata.p','w')
		pickle.dump(key_list,target)
	print "Finished sourcing Target Nodes data"
	target.close()

def source_tstptag(a_listfile): #F3
	gene = [gene.split(',')[0] for gene in open(a_listfile)]
	tag = [tag.split(',')[1].rstrip('\n') for tag in open(a_listfile)]
	key_list = dict(zip(gene, tag))
	if os.path.exists(cwl+'/dependencies') == True:
		target = open(cwl+'/dependencies/tstptag.p','w')
		pickle.dump(key_list,target)
	else:
		print "/dependencies directory not found, redirecting output to", cwl
		target = open(cwl+'/tstptag.p','w')
		pickle.dump(key_list,target)
	print "Finished sourcing TS-TP Tag data"
	target.close()

## Dependencies Function
def check_dependencies(): #F0
	print "Current Working Directory :", cwl
	print "Checking dependencies..."
	proceed_signal = 1
	if not os.path.exists(cwl+'/dependencies'):
		print "dependencies directory not found. Check for dependencies directory, update the requirements and try again.."
		proceed_signal = 0
		exit(0)
	else:
		print "dependencies directory found. Checking files..."
		if os.path.isfile(cwl+'/dependencies/global_indicies.p') == True:
			print "global Indicies list found"
		else:
			print "global Indicies list not found. Creating new list from scratch..."
			if os.path.isfile(cwl+'/dependencies/backup/global_indicies.txt') == True:
				source_global_indicies(cwl+'/dependencies/backup/global_indicies.txt') #f1
			else:
				print "No files to source from. Exiting..."
				proceed_signal = 0
				exit(0)
		if os.path.isfile(cwl+'/dependencies/global_critical_components.p') == True:
			print "global critical_components list found"
		else:
			print "critical_components list not found. Creating new list from scratch..."
			if os.path.isfile(cwl+'/dependencies/backup/global_critical_components.txt') == True:
				source_critical_components(cwl+'/dependencies/backup/global_critical_components.txt') #f2
			else:
				print "No files to source from. Exiting..."
				proceed_signal = 0
				exit(0)
		if os.path.isfile(cwl+'/dependencies/gene_matrix/list') == True:
			print "gene matrix found"
		else:
			print "gene matrix not found. Please create mutation curve before initiation patients"
			proceed_signal = 0
			exit(0)
		if os.path.isfile(cwl+'/dependencies/targetnodedata.p') == True:
			print "target-node data found"
		else:
			print "target-node not found. Creating new list from scratch..."
			if os.path.isfile(cwl+'/dependencies/backup/GS_Target_Nodes.csv') == True:
				source_targetnodes(cwl+'/dependencies/backup/GS_Target_Nodes.csv') #f3
			else:
				print "No files to source from. Exiting..."
				proceed_signal = 0
				exit(0)
		if os.path.isfile(cwl+'/dependencies/tstptag.p') == True:
			print "TS-TP data found"
		else:
			print "TS TP data not found. Creating new list from scratch..."
			if os.path.isfile(cwl+'/dependencies/backup/tstptag.txt') == True:
				source_tstptag(cwl+'/dependencies/backup/tstptag.txt') #f4
			else:
				print "No files to source from. Exiting..."
				proceed_signal = 0
				exit(0)
	return proceed_signal

def Build_Cancer_for_Patient(patient, pid, indc, mutation_file):
	print "Initiated patient creation for:", pid
	global_collateral = {}
	patient = Patient()
	patient.patientid = pid
	patient.indication = indc
	patient.source_mutation_components(mutation_file)
	critical_component_override_signal = 0
	anchored_mutations = []
	if len(patient.mutations) < 30:
		critical_component_override_signal = 1
		print "Overriding Critical Crosstalks check due to less mutations/CNV"
	for mutation in patient.mutations:
		gene_name = str(mutation)
		print "Gene: ", gene_name
		mutation = Mutation()
		if gene_name in tstptag:
			mutation.tag = tstptag[gene_name]
		else:
			mutation.tag = 'NEUTRAL'
			print "Mutation tag not found for", gene_name
		mutation.mutation_type = patient.mutation_mtype[gene_name]
		mutation.mutation_effect = patient.mutations[gene_name]
		mutation.anchor = patient.gene_anchor[gene_name]
		if mutation.anchor < 0 and mutation.mutation_type == 'CNVGOF' or mutation.mutation_type == 'MUTGOF':
			print "Illogical anchor value. Ignoring"
			mutation.anchor = None
		elif mutation.anchor > 0 and mutation.mutation_type == 'CNVLOF' or mutation.mutation_type == 'MUTLOF':
			print "Illogical anchor value. Ignoring"
			mutation.anchor = None
		if mutation.anchor != None:
			anchored_mutations.append(gene_name)
			print gene_name, "anchored to value", mutation.anchor
		i_score = mutation.source_indication_score(gene_name, indc)
		mutation.source_species_interactive_score(cwl+'/dependencies/gene_matrix/'+gene_name+'_curated_matrix')
#		mutation.source_mutation_elements(cwl+'/dependencies/gene_matrix/'+gene_name+'_mutation_elements')
		index_score = {}
		index_score_sorted = []
		for index in global_indicies:
			index_score[index] = abs(mutation.sis[index])
		for key, value in sorted(index_score.iteritems(), key=lambda (k,v): (v,k)):
    			index_score_sorted.append(key)
		index_selected = index_score_sorted[len(index_score)-1]
		#Starting with n% on selected index
		direc = 1
		if mutation.mutation_effect == 'CNVLOF' and mutation.tag == 'TSG':
			direc = float(-6)
		elif mutation.mutation_effect == 'MUTLOF' and mutation.tag == 'TSG':
			direc = float(-3)
		elif mutation.mutation_effect == 'CNVGOF' and mutation.tag == 'OG':
			direc = float(2.5)
		elif mutation.mutation_effect == 'MUTGOF' and mutation.tag == 'OG':
			direc = float(1.5)
		elif mutation.mutation_effect == 'CNVGOF' and mutation.tag == 'TSG':
			direc = float(0.5)
		elif mutation.mutation_effect == 'MUTGOF' and mutation.tag == 'TSG':
			direc = float(0.5)
		elif mutation.mutation_effect == 'CNVLOF' and mutation.tag == 'OG':
			direc = float(-0.05)
		elif mutation.mutation_effect == 'MUTLOF' and mutation.tag == 'OG':
			direc = float(-0.05)
		elif mutation.mutation_effect == 'CNVGOF' and mutation.tag == 'NEUTRAL':
			direc = float(0.5)
		elif mutation.mutation_effect == 'MUTGOF' and mutation.tag == 'NEUTRAL':
			direc = float(0.5)
		elif mutation.mutation_effect == 'CNVLOF' and mutation.tag == 'NEUTRAL':
			direc = float(-0.05)
		elif mutation.mutation_effect == 'MUTLOF' and mutation.tag == 'NEUTRAL':
			direc = float(-0.05)
		else:
			direc = float(0.1)
		mutation.pick_collateral_components()
		for component in mutation.collateral:
			if component in global_collateral:
				global_collateral[component] += mutation.sis[component]*(direc/abs(direc))*0.1
			else:
				global_collateral[component] = 1
		m_perb = 0
		a = len(patient.mutations)
		if abs(mutation.sis[index_selected]) >= 1 and a < 30:
			mutation.sis[index_selected] = mutation.sis[index_selected]/15
			print "Increasing contribution bias of", gene_name, ": L1"
		elif abs(mutation.sis[index_selected]) >= 1 and a > 30:
			mutation.sis[index_selected] = mutation.sis[index_selected]/10
			print "Increasing contribution bias of", gene_name, ": L2"
#		min(base_mca_contri) ~ 4 varies between 24% to 4% 
		base_mca_contri = float(2.5/(a/(50+a)))
		if abs(mutation.sis[index_selected]) >= 0.00001:
			if mutation.anchor == None:
##########################################################################################################

				m_perb = float(direc*base_mca_contri/mutation.sis[index_selected])*i_score

##########################################################################################################
				if m_perb <= -99:
					m_perb = -75.00
				if m_perb > 250:
					m_perb = 250.00
			else:
				m_perb = float(mutation.anchor)
			patient.disease[gene_name] = float(m_perb)
		if critical_component_override_signal == 0 and mutation.anchor == None:
			for component in global_critical_components:
				if m_perb*mutation.sis[component] >= 35:
					print gene_name, "identified as regulator of critical component: ", component, ". Reducing effect"
					patient.disease[gene_name] = patient.disease[gene_name]/1.5
				elif m_perb*mutation.sis[component] >= 20 and m_perb*mutation.sis[component] < 35:
					print gene_name, "identified as regulator of critical component: ", component, ". Reducing effect"
					patient.disease[gene_name] = patient.disease[gene_name]/1.2
			else:
				pass
		if len(patient.disease) <= 1:
			pass
		else:
			for component in patient.disease:
				temp = {}
				if component != gene_name:
					target = open(cwl+'/dependencies/gene_matrix/'+component+'_curated_matrix','r')
					temp = pickle.load(target)
					query = str('_'+component+'_')
					for entry in mutation.sis:
						if query in entry:
							if mutation.sis[entry] >= 0.01:
								offset = m_perb*mutation.sis[entry]
								change_by_component = offset*temp[index_selected]
								if patient.disease[component] > change_by_component and change_by_component >= 0.25:
									if offset > 0:
										patient.disease[component] = patient.disease[component] - offset
										print "Modified gene: ", component
#									elif offset < 0:
#										patient.disease[component] = patient.disease[component] + offset/4
#										print "Modified gene: ", component
							elif mutation.sis[entry] < -0.01:
								offset = m_perb*mutation.sis[entry]
								change_by_component = offset*temp[index_selected]
								if patient.disease[component] > change_by_component and change_by_component >= 0.25:
									if offset > 0:
										patient.disease[component] = patient.disease[component] - offset
										print "Modified gene: ", component
#									elif offset < 0:
#										patient.disease[component] = patient.disease[component] - offset/4
#										print "Modified gene: ", component
					target.close()
#	Write_Cancer_for_Patient
	if not os.path.exists(cwl+'/'+pid):
		try:
			os.makedirs(cwl+'/'+pid)
			print "Created directory ", pid
		except:
			print "Cannot create directory. Either directory name already exists in path or Unknown error"
			exit(-1)
	tx = open(cwl+'/'+pid+'/global_collateral.csv','w')	
	for keys in global_collateral:
		if global_collateral[keys] >= 1.5:
			string = str(keys)+','+str(global_collateral[keys])
			tx.write(string)
			tx.write('\n')
	tx.close()
	target0 = open(cwl+'/'+pid+'/'+pid+'.p','w')
	pickle.dump(patient.disease,target0)
	target0.close()
	target = open(cwl+'/'+pid+'/'+pid+'.cwd','w')
	for mutation in patient.mutations:
		header = '#MUT-START-'+str(mutation)
		target.write(header)
		target.write('\n')
		for entry in targetnodedata:
			if targetnodedata[entry] == mutation:
				if '_Txn' in entry or '_Src' in entry or '_Frm' in entry:
					if patient.mutations[mutation] == 'CNVLOF':# and patient.disease[mutation] <= 0:
						if abs(patient.disease[mutation]) > 75:
							patient.disease[mutation] = float(75.00)
						string = str('ADD KD REACTION '+entry+' '+str(abs(patient.disease[mutation])))
						target.write(string)
						target.write('\n')
					elif patient.mutations[mutation] == 'CNVGOF':# and patient.disease[mutation] >= 0:
						string = str('ADD OE REACTION '+entry+' '+str(abs(patient.disease[mutation])))
						target.write(string)
						target.write('\n')
					else:
						pass
				else:
					if patient.mutations[mutation] == 'CNVLOF':# and patient.disease[mutation] <= 0:
						if abs(patient.disease[mutation]) > 75:
							patient.disease[mutation] = float(75.00)
						string = str('ADD KD SPECIES '+entry+' '+str(abs(patient.disease[mutation])))
						target.write(string)
						target.write('\n')
					elif patient.mutations[mutation] == 'CNVGOF':# and patient.disease[mutation] >= 0:
						string = str('ADD OE SPECIES '+entry+' '+str(abs(patient.disease[mutation])))
						target.write(string)
						target.write('\n')
					else:
						pass
		target.write('\n')
	target.close()		

#	Push_Envelope_for_Patient

	envs = [2, 1.42, 1.25, 1.11, 0.91, 0.8, 0.74, 0.66, 0.5]
	for env in envs:
		d = {}
		for component in patient.disease:
			if component not in anchored_mutations:
				d[component] = float(patient.disease[component]*env)
			else:
				d[component] = float(patient.disease[component])
		target = open(cwl+'/'+pid+'/'+pid+'_f_'+str(env)+'.cwd','w')
		for mutation in patient.mutations:
			header = '#MUT-START-'+str(mutation)
			target.write(header)
			target.write('\n')
			for entry in targetnodedata:
				if targetnodedata[entry] == mutation:
					if '_Txn' in entry or '_Src' in entry or '_Frm' in entry:
						if patient.mutations[mutation] == 'CNVLOF':# and patient.disease[mutation] <= 0:
							if abs(d[mutation]) > 99:
								d[mutation] = 99.00
							string = str('ADD KD REACTION '+entry+' '+str(abs(d[mutation])))
							target.write(string)
							target.write('\n')
						elif patient.mutations[mutation] == 'CNVGOF':# and patient.disease[mutation] >= 0:
							string = str('ADD OE REACTION '+entry+' '+str(abs(d[mutation])))
							target.write(string)
							target.write('\n')
						else:
							pass
					else:
						if patient.mutations[mutation] == 'CNVLOF':# and patient.disease[mutation] <= 0:
							if abs(d[mutation]) > 99:
								d[mutation] = 99.00
							string = str('ADD KD SPECIES '+entry+' '+str(abs(d[mutation])))
							target.write(string)
							target.write('\n')
						elif patient.mutations[mutation] == 'CNVGOF':# and patient.disease[mutation] >= 0:
							string = str('ADD OE SPECIES '+entry+' '+str(abs(d[mutation])))
							target.write(string)
							target.write('\n')
						else:
							pass
			target.write('\n')
		target.close()
	print "Finished generation patient with id", pid		

def Initiate(a_list):
	for ids in a_list:
		patients_in_pipeline.append(ids)
	print "Queued patients from setup"

if __name__ == "__main__":
	if len(sys.argv) != 2 or sys.argv[1] == '--help':
		print "Usage: PatientCreator.py <setup-file>"
		print "Usage: Setup File Format : <patient-id>,<indication-id>,<path-to-mutation-file>"
		# NZEC for chaining
		sys.exit(-1)
	script, setupfile = sys.argv
	print "~~~Welcome to Patient Creator Suite~~~"
	print ""
	patientids = [patientids.split(',')[0] for patientids in open(setupfile)]
	indicationids = [indicationids.split(',')[1] for indicationids in open(setupfile)]
	mutationinfo_file = [mutationinfo_file.split(',')[2].rstrip('\n') for mutationinfo_file in open(setupfile)]

	patientids_indications = dict(zip(patientids, indicationids))
	patientids_mutationifo_file = dict(zip(patientids, mutationinfo_file))
	y = check_dependencies()
	if y == 1:
		Initiate(patientids)
		# Expanding all dependent lists
		target = open(cwl+'/dependencies/global_indicies.p','r')
		global_indicies = pickle.load(target)
		target.close()
		target = open(cwl+'/dependencies/global_critical_components.p','r')
		global_critical_components = pickle.load(target)
		target.close()
		target = open(cwl+'/dependencies/tstptag.p','r')
		tstptag = pickle.load(target)
		target.close()
		target = open(cwl+'/dependencies/targetnodedata.p','r')
		targetnodedata = pickle.load(target)
		target.close()

		for entry in patients_in_pipeline:
			###
			Build_Cancer_for_Patient(entry,str(entry),patientids_indications[entry],patientids_mutationifo_file[entry])
			###
