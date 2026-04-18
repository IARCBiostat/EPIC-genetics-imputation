#!/usr/bin/env python2.7

# -f or --fam : Fam file (mandatory)
# -i or --id : an ID file (mandatory, i.e. Subj_Id_2015.txt)
# -q : search for questionnaire IDs and replace them by sample IDs
# -l or --link : a file to link Laboratory IDs to Epic IDs
# -s : specify the field with the "Fam file ID" (starting from 0)
# -e : specify the field with the "EPIC ID" (starting from 0)
# -nlh : Number of header lines in the link file
# -k : keep IDs that are not in the link file
# -c : specify 8, 9 or B. Usefull only if one of the countries 8, 9 or B have IDs without the country code.


#output:
# - <project>_goodID.txt
# - <project>_removeID.txt

#written by Manon Knuchel, September 2019

##MODULES
import sys
import argparse #module to make it easy to write user-friendly command-line interfaces
import re #module that allows to do regular expressions (= regex) search

##############################Input file options##############################

#you can access to the help section with '-h' option

parser = argparse.ArgumentParser()

parser.add_argument('--fam', '-f', type=str, help='Path to the Fam file.', required=True)
parser.add_argument('--id', '-i', type=str, help='Path to the EPIC ID file.', required=True)
parser.add_argument('--questID','-q', action='store_true', help='If this option is present, the script will search Fam file IDs within Questionnaire IDs instead of Sample IDs')
parser.add_argument('--link', '-l', type=str, help='Path to the link file.')
parser.add_argument('--sample', '-s', type=int, help='Number of the field in the link file with the Fam file ID.')
parser.add_argument('--epicID', '-e', type=int, help='Number of the field in the link file with the EPIC ID.')
parser.add_argument('--nbLineHead', '-nlh', type=int, help='Number of header lines in the link file.')
parser.add_argument('--keep','-k', action='store_true', help='if we want to keep IDs that are not in the link file')
parser.add_argument('--country', '-c', choices=['8', '9', 'B', '0'], help='Country (8, 9 or B) for IDs without country code. 0 to delete these IDs.')
parser.add_argument('--summary_file', '-sum_f', type=str, help='Path to the Normalization summary files.')
args = parser.parse_args()

fam_file = args.fam # option -f
fam_file_name = fam_file.split('/')[-1].split('.')[0] #get the 'project name'

summary_file = args.summary_file

epic_file = args.id # otion -i

id_quest = args.questID # option -q
if id_quest == None:
	id_quest = False

link_file = args.link # option -l

field_sample = args.sample # option -s

field_epic = args.epicID # option -e

nblhead = args.nbLineHead # option -nlh
if nblhead == None:
	nblhead = 0

keep_ids = args.keep #option -k
if keep_ids == None:
	keep_ids = False

missing_country = args.country # option -c

# if a file is specified, we need both fields to do the ID linkage
if (link_file != None and (field_sample == None or field_epic == None)) or (link_file == None and (field_sample != None or field_epic != None)):
	sys.stderr.write("\n\n**********************************************************\n\nError : If you specify a link file, you have to specify the file name (-l) and the fields with the Fam file ID (-s) and with EPIC ID (-e).\n\n**********************************************************\n\n")
	exit(0)


##############################Storage of the link file##############################

dico_link = {} # dictionnary to save the link file
flag_sep = 0

if link_file != None:
	with open(link_file) as file_link:
		for nb, line in enumerate(file_link):
			if nb >= nblhead : # skip header lines
				if flag_sep == 0:
					if re.search(r'(\t)+',line): #try if the separators are tabulations
						flag_sep = 1
					elif re.search(r'(,)+',line): #try if the separators are commas
						flag_sep = 2
					elif re.search(r'( )+',line): #try if the separators are spaces
						flag_sep = 3
					else:
						sys.stderr.write("\n\n**********************************************************\n\nError : Separators of the link file are neither tabulations nor commas or spaces !\n\n**********************************************************\n\n")
						exit(0)

				if flag_sep == 1:
					data = line.split('\t')
				elif flag_sep == 2:
					data = line.split(',')
				else:
					data = line.split(' ')

				lab_ID = data[field_sample].strip() # Get the ID which should be in the Fam file
				epic_ID = data[field_epic].strip() # Get the ID which should correspond to an EPIC ID
				
				if epic_ID :
					dico_link[lab_ID] = epic_ID # creation of the link dictionnary


##############################Storage of the EPIC ID file##############################

dico_epic = {} # dictionnary to save the database

with open(epic_file) as file_id: # file Subj_Id_2015.txt
	for line in file_id:
		data = line.split(',') #DON'T use the csv module because the id '00123' will become '123'......
		country = data[0].strip() # Get all data needed
		center = data[1].strip()
		questID = data[2].strip()
		sex = data[3].strip()
		sampleID1 = data[5].strip()
		sampleID2 = data[7].strip()

		if sampleID1 != '': # we keep only subjects (one line in the database) who have at least one sample ID
			place = country + center
			#duplicate1 = 0
			#duplicate2 = 0

			if sampleID2 != '': # if they have two sample IDs
				subject = [place, sex, sampleID1, questID, sampleID2]
				#subject = [place, duplicate1, sampleID1, questID, duplicate2, sampleID2]
			else: # if they have one sample
				subject = [place, sex, sampleID1, questID]

			if place not in dico_epic: # dictionnary per "country center" with an array of subjects of that center
				dico_epic[place]=[]

			dico_epic[place].append(subject)


##############################Function to treat IDs from Fam file##############################

def test_ID_equality(detectedID, subject, ID_found): # function that test the equality of the ID detected in the Fam file (or link file) with the sampleID (IDBio or IDBio2 if it exists) and store it in an array
	if detectedID == subject[2]: # test equality with IDBio
		ID_found.append(subject) # if found, add the result to an array
	elif len(subject) == 5 and detectedID == subject[4]: # test equality with IDBio2
		ID_found.append(subject)

def compute_last_digit(ID): # for some center (25, 42, 51, 52), the last digit may have been replaced by 0.
	som_mod = (int(ID[0]) + (int(ID[1])*2)%10 + int(ID[2]) + (int(ID[3])*2)%10 + int(ID[4]) + (int(ID[5])*2)%10 + int(ID[6]))%10 # compute the last digit
	last_digit = (10-som_mod)%10

	compute_ID = ID[0] + ID[1] + ID[2] + ID[3] + ID[4] + ID[5] + ID[6] + str(last_digit) # the function sent back the entire ID with the good last digit
	return compute_ID


##############################Read the Fam file##############################

#declaring counter
nb_IDFam = 0 # number of subjects in the initial Fam file
nb_IDfound = 0 # number of Fam file IDs linked to an EPIC ID
nb_IDremove = 0 # number of Fam file IDs not linked to an EPIC ID
nb_IDbio2 = 0  # number of IDbio2 found in the Fam file
nb_IDdiffResult = 0 # number of Fam file IDs linked to more than one EPIC ID
nb_sampleQC = 0 # number of Fam file IDs marked with the sign of QC

flag_sep = 0 # to know which separator is used in the Fam file

dico_result = {} # dico to save results to be able to count duplicate and write them in the *_goodID.txt file

with open(fam_file, 'r') as file_fam, open(fam_file_name + '_removeID.txt', 'w') as remove_file:#, open(fam_file_name + '_duplicateID.txt', 'w') as duplicate_file:
	for line in file_fam:
		nb_IDFam +=1

		if flag_sep == 0:
			sep = re.search(r'(\t)+',line) #try if the separators are tabulations
			if sep: # sep!= None : separators are tabulations
				flag_sep = 1
			else: #sep == None : separators are spaces
				flag_sep = 2

		if flag_sep == 1: # separators are tabulations
			data = line.split('\t')
		else: # separators are spaces
			data = line.split(' ')

		id_famille = data[0].strip() # Get Familly ID and Subject ID from the Fam file
		id_lab = data[1].strip()

		id_work = ''
		if link_file != None: # if we have a link file
			if id_lab in dico_link:# and dico_link[id_lab] != '':
				id_work = dico_link[id_lab] # we will work with the linked ID
			elif keep_ids == True:
				id_work = id_lab # if we want to keep IDs that are not in the link file, we will work with the ID of the Fam file
			else:
				remove_file.write(str(id_famille) + "\t" + str(id_lab) + "\n") # if the ID is not in the file, we delete it
				nb_IDremove +=1
		else:
			id_work = id_lab # if we don't have a link file, we will work with the ID of the Fam file
		
		#start checking id_work
		flag_QCsample = False
		if re.search(r".+<qc>$",id_work):
			flag_QCsample = True

		ID_found = [] # the array will contain all results found in the EPIC database for the current Fam file ID
		detectedID = ''

		completeID = re.search(r"([^a-zA-Z0-9]|^)([B1-9][0-9])[-_]+(((P[0-9]{4,})|([0-9]{6,}A[0-9]{3,})|([0-9]+))[A-Z]?)([^a-zA-Z0-9]|$)", id_work)
		if completeID: # if the Fam file ID has two numbers (or B and a number), one or more separators and a string that can be an ID (ex : B1_______P0483 or 82__248296A482 or 52____52526275 or 81___________1)
			place = completeID.group(2) # we get the potential "country center"
			detectedID = completeID.group(3) # and the potential ID

			if place[1] == '0': # if only the country is known (case happens only if we modify the link file)
				for i in range(1,7): # we test all possible centers (1 to 6 maximum)
					index = place[0] + str(i)
					if index in dico_epic:	
						for value in dico_epic[index]: # test ID equality with subjects of the database
							if id_quest: # if we know that the Fam file ID corresponds to a questionnaire ID
								if detectedID == value[3]: # we test the equality with the questionnaire ID
									ID_found.append(value)
							else: # else we test the equality with IdBio and IdBio2
								test_ID_equality(detectedID, value, ID_found)
								
			else:	
				if place in dico_epic: # we search within subjects of the "center country"
					for value in dico_epic[place]: # test ID equality with subjects of the database
						if id_quest : # search questionnaire ID
							if detectedID == value[3]:
								ID_found.append(value)
						else: # search IdBio and IdBio2 ID
							if len(detectedID) == 8 and detectedID[7] == '0': # if the ID ends by 0, recompute the last digit
								detectedID = compute_last_digit(detectedID)
							test_ID_equality(detectedID, value, ID_found) 
							
		else: # Fam file ID does not have an EPIC ID format
			partID = re.search(r"([^a-zA-Z0-9]|^)(((P[0-9]{4,})|([0-9]{6,}A[0-9]{3,})|([0-9]+))[A-Z]?)([^a-zA-Z0-9]|$)", id_work)
			if partID: #if we found a string that may correspond to an Idbio

				detectedID = partID.group(2)
				length_id = len(detectedID) # depending of the length, we can determine more information

				if id_quest: 
					if length_id == 10 and detectedID[0] == '2' and detectedID[1] == '4': # Fam file ID from Umea
						for value in dico_epic["82"]:
							if detectedID == value[3]:
								ID_found.append(value)
					
					else : #Other Fam file IDs are rejected (without "country center", the questionnaire IDs are not unique)
						sys.stderr.write("\n\n**********************************************************\n\nError : We have a questionnaire ID without country/center code: " + str(id_work) + "\n\n**********************************************************\n\n")
						remove_file.write(str(id_famille) + "\t" + str(id_lab) + "\n")
						id_work = ''
						nb_IDremove += 1
						#exit(0)

				elif length_id == 10: # IDs with a length of 10 are from Umea (82)
					for value in dico_epic["82"]:  
						test_ID_equality(detectedID, value, ID_found)	# within subjects of this place we search IdBio and IdBio2

				elif length_id == 8: # can be country from 1 to 7. IDs with a length of 8 have : [country] [center] [5 numbers] [1 number computed with the 7 others] 
					if detectedID[7] == '0':
						detectedID = compute_last_digit(detectedID) # we compute the last digit

					place = detectedID[0] + detectedID[1] # we extract "country" and "center" from ID
					
					if detectedID[0] == '1': # the IDs for France do not contain the center code (the first two numbers are 11)
						for i in range(1,7): # for each French center
							index = place[0] + str(i)
							if index in dico_epic:
								for value in dico_epic[index]: #test equality
									test_ID_equality(detectedID, value, ID_found)
					else: # for country from 2 to 7
						if place in dico_epic:
							for value in dico_epic[place]:  #test equality
								test_ID_equality(detectedID, value, ID_found)
							
				elif length_id <= 5:# three countries/centers can have an ID length <= 5 : 81, 9 and B
					if missing_country == None:# IDs from these three countries might be similar.
						sys.stderr.write("\n\n**********************************************************\n\nError : we can't define the country of some IDs, please specify a country with the option -c or do a pretreatment on the file. "+ id_work +"\n\n**********************************************************\n\n")
						exit(0)

					elif missing_country == '0':
						remove_file.write(str(id_famille) + "\t" + str(id_lab) + "\n") # we remove it
						id_work = ''

					elif missing_country == '8':
						for value in dico_epic["81"]:
							test_ID_equality(detectedID, value, ID_found) #test equality
					else :
						while len(detectedID) < 5:
							detectedID = '0' + detectedID #we increase the size of the ID to get a length of 5 (only for country 9 and B)
						for i in range(1,3): # for each center
							index = missing_country + str(i)
							if index in dico_epic:
								for value in dico_epic[index]:
									test_ID_equality(detectedID, value, ID_found) #test equality

	
		if len(ID_found) == 0: # ID not found in the database
			if id_work != '':
				remove_file.write(str(id_famille) + "\t" + str(id_lab) + "\n") # we remove it
				nb_IDremove +=1
		
		elif len(ID_found) == 1:
			nb_IDfound += 1
			lab_info = (id_famille,id_lab)
			
			if id_quest == False and len(ID_found[0]) == 5 and detectedID != ID_found[0][2]:# and detectedID == ID_found[0][4]:
				new_ID = ID_found[0][4] # if IdBio2 is detected, we write Idepic_Bio2 as new EPIC ID
				nb_IDbio2 += 1
			else: # if IdBio or questionnaire ID is detected, we write Idepic Bio as new EPIC ID
				new_ID = ID_found[0][2]

				
			while len(new_ID) < 12: # we put as many '_' as necessary to have 14 characters in the ID
				new_ID = '_' + new_ID
			new_ID = ID_found[0][0] + new_ID

			# we also get the questionnaire ID as a familly ID
			fam_quest = ID_found[0][3]
			while len(fam_quest) < 12: # we put as many '_' as necessary to have 14 characters in the ID questionnaire
				fam_quest = '_' + fam_quest
			fam_quest = ID_found[0][0] + fam_quest



			if new_ID not in dico_result: # if it is the first time we have this ID, we have to create the structure
				tab_replicate = []
				tab_qc = []
				dico_result[new_ID]=[ID_found[0][1],tab_replicate,tab_qc,fam_quest]

			if flag_QCsample: # if the ID is a QC, we write it in the QC array else we write it in the normal array (where we can have duplicate)
				dico_result[new_ID][2].append(lab_info)
				nb_sampleQC += 1
			else :
				dico_result[new_ID][1].append(lab_info)
			
		else : # more than one Epic ID found
			sys.stderr.write("\n\n**********************************************************\n\nError with this Fam file ID that match with " + str(len(ID_found)) + " different EPIC IDs : " + str(id_lab) + " ("+str(id_work)+")\n\n**********************************************************\n\n")
			remove_file.write(str(id_famille) + "\t" + str(id_lab) + "\n") # we remove it
			nb_IDdiffResult += 1
			#nb_IDremove += 1



nb_ID_duplicate = 0 # number of IDs found more than once in the Fam file
nb_times_duplicate = 0 # number of IDs found more than once in the Fam file

with open(fam_file_name + '_goodID.txt', 'w') as res_file :
	for key, value in dico_result.items():
		if len(value[1]) == 1:
			# format : old ID_Fam, old ID_ind, ID questionnaire, new_ID_ind, sex
			res_file.write(str(value[1][0][0]) + "\t" + str(value[1][0][1]) + "\t" + str(value[3]) + "\t" + str(key) + "\t" + str(value[0]) + "\n")

		elif len(value[1]) > 1:
			nb_ID_duplicate += 1
			nb_times_duplicate += len(value[1])
			for i in range(len(value[1])):
				# we write new_ID_ind_R1 _R2 ...
				res_file.write(str(value[1][i][0]) + "\t" + str(value[1][i][1]) + "\t" + str(value[3]) + "\t" + str(key) + "_R" + str(i+1) + "\t" + str(value[0]) + "\n")

		if len(value[2]) == 1:
			res_file.write(str(value[2][0][0]) + "\t" + str(value[2][0][1]) + "\t" + str(value[3]) + "\t" + str(key) + "_QC\t" + str(value[0]) + "\n")

		elif len(value[2]) > 1:
			for i in range(len(value[2])):
				res_file.write(str(value[2][i][0]) + "\t" + str(value[2][i][1]) + "\t" + str(value[3]) + "\t" + str(key) + "_QC" + str(i+1) + "\t" + str(value[0]) + "\n")



sys.stdout.write("\n\n**********************************************************\n\nNumber of subjects in the Fam file: " + str(nb_IDFam) +"\n- Number of EPIC IDs identified: " + str(nb_IDfound) +"\n- Number of IDs that are not (or not identified as) EPIC IDs: " + str(nb_IDremove) + "\n- Number of IDs that match more than one EPIC ID: " + str(nb_IDdiffResult))
sys.stdout.write("\n\nAbout the EPIC IDs identified:\n- Number of QC sample IDs: " + str(nb_sampleQC) + "\n- Number of duplicates: " + str(nb_times_duplicate) + "\n- Number of distinct IDs with duplicates: " + str(nb_ID_duplicate) + "\n- Number of 2nd blood sample IDs: " + str(nb_IDbio2) +"\n\n**********************************************************\n\n")

if summary_file != None:
	with open(summary_file, 'a') as out_file :
		out_file.write("\n*************************** ID Part *******************************\n\nNumber of subjects in the Fam file: " + str(nb_IDFam) +"\n- Number of EPIC IDs identified: " + str(nb_IDfound) +"\n- Number of IDs that are not (or not identified as) EPIC IDs: " + str(nb_IDremove) + "\n- Number of IDs that match more than one EPIC ID: " + str(nb_IDdiffResult))
		out_file.write("\n\nAbout the EPIC IDs identified:\n- Number of QC sample IDs: " + str(nb_sampleQC) + "\n- Number of duplicates: " + str(nb_times_duplicate) + "\n- Number of distinct IDs with duplicates: " + str(nb_ID_duplicate) + "\n- Number of 2nd blood sample IDs: " + str(nb_IDbio2) +"\n\n**********************************************************\n\n")

