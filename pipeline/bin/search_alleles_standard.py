#!/usr/bin/env python2.7

#script that finds the alleles for a SNP, if the manifest file in csv is available
#input: 3 parameters:
# a bim file (mandatory: -b or --bim)
# a manifest file (Illumina csv) (mandatory: -m or --man)
# a loci to rsID file (Illumina txt) (optionnal: -t or --txt)
#output:
# - [file]_goodAlleles_[genomeVersion].txt
# - [file]_goodAlleles_[genome_version]_results.txt
# - [file]_goodAlleles_strand.txt

#bsub -eo search_alleles_csv_err.txt -oo search_alleles_csv_out.txt 'source /opt/rh/python27/enable ; /data/epic-nmb/Scripts/Normalisation/Completion_Lift/Alleles/search_alleles_csv.py -b [bim_file] -m [csv_file_from_Illumina] -t [loci_to_rsID_file]'

#script written by Manon Knuchel in May 2019 based on the first version developped in May 2018 by Emilie Gerard-Marchant & Benjamin Bourgeois

##MODULES
import argparse #module to make it easy to write user-friendly command-line interfaces
import csv #module to read csv
#import pprint #module to print 'pretty' dictionnary
import re #module that allows to do regular expressions (= regex) search
#import time #module to see how long the script takes
import sys # module to write in the error output

#declaring complementary allele dictionary
complementary = {"A":"T", "T":"A", "C":"G", "G":"C"}


##############################Input file options##############################

#you can access to the help section with '-h' option

parser = argparse.ArgumentParser()

parser.add_argument('--bim', '-b', type=str, help='Path to the bim file.')
parser.add_argument('--man', '-m', type=str, help='Path to the manifest file (Illumina.csv).')
parser.add_argument('--txt', '-t', type=str, help='Path to the loci to rsID file (Illumina.txt).')
parser.add_argument('--build', '-bd', type=str, help='build version of the bim/manifest file')

args = parser.parse_args()

bim_file = args.bim
if bim_file == None: #check if the bim file is given
	sys.stderr.write("\n\n**********************************************************\n\nNo bim file given. Argument '-b' is mandatory.\n Exiting...\n\n**********************************************************\n\n")
	exit(0)
bim_file_name = bim_file.split('/')[-1]

manifest_file = args.man
if manifest_file == None: #check if the manifest file is given
	sys.stderr.write("\n\n**********************************************************\n\nNo manifest file (Illumina csv). Argument '-m' is mandatory.\n Exiting...\n\n**********************************************************\n\n")
	exit(0)

genome_version = args.build #variable that will contain the genome version
if genome_version == None:  #check if we know th genome version
	sys.stderr.write("\n\n**********************************************************\n\nNo version build specified.  Argument '-bd' is mandatory.\n Exiting...\n\n**********************************************************\n\n")
	exit(0)

file_rsID = args.txt

##############################Storage of the bim file##############################

bim_dico = {} # dictionnary where all the data from the bim file will be stored, it will be like : {line_nb: ('snp', 'allele_1', 'allele_2')}

file_bim_length = 0 #variable that will contain the number of lines of the old bim file = number of line of the new bim file

allele1_test = True #variable that will tell us about the presence of the 1st allele info in the bim file: True if present, False if missing
allele2_test = True #variable that will tell us about the presence of the 2nd allele info in the bim file: True if present, False if missing

with open(bim_file) as file_bim: #bim_file  is the bim file given in the terminal
	for line_nb, line in enumerate(file_bim):
		sep = re.search(r'(\t)+',line) #try if the separators are tabulations

		if sep: # separator are tabulations, sep != None
			data = line.split('\t')
		else: #sep == None : separators are spaces
			data = line.split(' ')

		#get all the datas we need
		snp = data[1]
		allele_1 = data[4] #minor allele
		allele_2 = data[5].strip() #major allele

		if allele_1 == '0':
			allele1_test = False #we will need to change the allele info later
		if allele_2 == '0':
			allele2_test = False #we will need to change the allele info later

		bim_dico[line_nb] = (snp, allele_1, allele_2) #creation of the dictionnary, key = int, value = tuple
		file_bim_length = line_nb + 1 # = number of SNPs in the bim file

#pprint.pprint(bim_dico)

##############################Storage of the Illumina file##############################

#col_names = ['Chr','MapInfo','Name','Allele1','Allele2','Strand'] #the name of the columns in the csv file

#file_illumina_length = 0  #variable that will contain the number of lines in the Illumina file

alleles = {} #dictionnary where all the necessary datas from Illumina will be stored, it will be like {'snp': ('allele_1', 'allele_2')} all will be in + strand

#with open(manifest_file) as illumina_file: #sys.argv[2] is the Illumina file (csv) given in the terminal
#	for line in illumina_file:
#		file_illumina_length += 1 #we count the number of lines

#	print(file_illumina_length)

with open(manifest_file) as illumina_file: 
	reader = csv.reader(illumina_file, delimiter='\t')
	for line_nb, row in enumerate(reader):
		#we get all the datas we want
		snp = row[2]
		allele_1 = row[3].strip()
		allele_2 = row[4].strip()
		strand_not = row[5].strip()
		
		# verifier contenu du bim
		# pas de AB dans manifest (pas besoin de le traiter)
		
		#we check if the nomenclature is '-' or 'N' and we put it missing (0)
		if allele_1 == '-':
			allele_1 = '0'
		if allele_2 == '-':
			allele_2 = '0'
		if allele_1 == 'N'or allele_2 == 'N' or allele_1 == allele_2:
			allele_1 = '0'
			allele_2 = '0'
			
		# ID si on ne remplace pas par les valeurs, ne rien faire	
		#if allele_1 == 'I' or allele_1 == 'D' or allele_2 == 'I' or allele_2 == 'D':
			#reccuperer un autre champs : determiner avec la longueur le D et le I (si egal, chercher si on a un "-" qui correspond au D): mettre toujours dans le meme sens pour faciliter la completion ?
		#	if strand_not == "BOT" or strand_not == "R" or strand_not == "MINUS" or strand_not == "-":
				#!!! le complementaire de ATG est CAT

		# if we have regular alleles (ATCG, 0), if they are on the negative strand, we change it : (the dictionnary have all his SNP on the same strand)	
		if (allele_1 in complementary or allele_1 =='0') and (allele_2 in complementary or allele_2 =='0'):
			if strand_not == "BOT" or strand_not == "R" or strand_not == "MINUS" or strand_not == "-":
				if allele_1 != '0':
					allele_1 = complementary.get(allele_1)
				if allele_2 != '0':
					allele_2 = complementary.get(allele_2)
			
		
		alleles[snp] = (allele_1, allele_2) #creation of the dictionnary, key = str, value = tuple

# pprint.pprint(alleles)

##############################Storage of the loci to rsID file datas##############################

snp_name = {} #dictionnary that will be like: {'kgp...': 'rs...'}

if file_rsID:
	with open(file_rsID) as loci_file: #file_rsID is the loci to rsID file given in the terminal
		for line_nb, line in enumerate(loci_file):
			
			if line_nb >= 1: #because the 1st line is the name of the columns 
				data = line.split('\t')
				#we get all the data we need
				new_snp = data[0] #the SNP that is not necessarily in the 'rs[number]' nomenclature
				good_snp = data[1].strip() #the SNP that is in the 'rs[number]' nomenclature
				if good_snp != '.': #there is a name for the SNP
					if ',' in good_snp : #there was a merge between several SNPs
						snps = []
						snps = good_snp.split(',')
						for var in snps:
							snp_name[var] = new_snp
					else:
						snp_name[good_snp] = new_snp ##construction of the dictionnary,key = str, value = str
						

#pprint.pprint(snp_name)

##############################Creation of the new bim file##############################

def find_other_allele(man_A1, man_A2, snp, a_know, liste_allele_error): 
	"""Function that choose the 2nd allele
	"""
	if a_know == man_A1: #the allele we know is = to allele_1
		other_allele = man_A2 #the other allele is = to allele_2
	elif a_know == man_A2: #the allele we know is = to allele_2
		other_allele = man_A1 #the other allele is = to allele_1
	else: # allele known is not one of the two alleles in the manifest, it should not happen
		add_liste_pb(liste_allele_error, snp)
		other_allele = '0' # keep the second allele missing
	return other_allele, liste_allele_error

def add_liste_pb(liste, snp):
	"""Function that adds the 10 1st SNP name in the problem list: if the chr information is missing
	"""
	if len(liste) < 10:
		if snp not in liste:
			liste.append(snp)

def write_info(fichier, snp, old_allele1, old_allele2, new_allele1, new_allele2):
	"""Function that write all the info in the file in order to allow modification in the bim file with Plink 
	"""
	fichier.write(snp + '\t' + old_allele1 + '\t' + old_allele2 + '\t' + new_allele1 + '\t' + new_allele2 + '\n') 

def replace_alleles_PM(allele_know, allele_2, rs, liste_allele_error): # remplace the allele nomenclature, +/- by I/D
	if allele_know == '+':
		allele1_new = 'I' # allele1_new correspond to allele_know
		allele2_new = 'D'
	elif allele_know == '-':
		allele1_new = 'D'
		allele2_new = 'I'
	else:
		add_liste_pb(liste_allele_error, bim_dico[nb][0])
		allele1_new = allele_know
		allele2_new = allele_2

	return allele1_new, allele2_new, liste_allele_error

cancer_bim = bim_file_name.split('.')[0] #the name of the bim file without the extension (.bim)
file_name = cancer_bim + '_goodAlleles_' + str(genome_version) + '.txt'

allele1_buf = 0 #variable that will contain the number of allele1 != 0 at the end of this script
allele2_buf = 0 #variable that will contain the number of allele2 != 0 at the end of this script

#ambiguous_SNP = 0 #variable to count new ambiguous SNP from the manifest

liste_allele_pb = [] #list that will contain the 10 1st SNPs which have a nomenclature which need to be change but were not found in the manifest file 
liste_allele_error = [] #list that will contain the 10 1st SNPs which have incoherent allele with the manifest

with open(file_name, 'w') as info_file:
	if allele1_test == False or allele2_test == False:
		for nb in range(file_bim_length): #to keep the order of the SNP the same as in the former bim file
			snp = ""
			if bim_dico[nb][0] in alleles: #the SNP is in the manifest file
				snp = bim_dico[nb][0]
			elif file_rsID != None and bim_dico[nb][0] in snp_name:#the SNP is in the loci to rsID file
				if snp_name[bim_dico[nb][0]] in alleles: #the old name is in the manifest file
					snp = snp_name[bim_dico[nb][0]]


			if snp != "" : # if we find the snp directly in the manifest or if we find it with the rsID file
				if (bim_dico[nb][1] not in complementary and bim_dico[nb][1] != '0') or (bim_dico[nb][2] not in complementary and bim_dico[nb][2] != '0'): # when allele are neither ATCG nor 0
					# !!! "snp" variable a utiliser comme clef pour trouver le SNP dans le dico "alleles"
					
					# to change the allele nomenclature I/D
					if bim_dico[nb][1] == 'I' or bim_dico[nb][2] == 'I' or bim_dico[nb][1] == 'D' or bim_dico[nb][2]== 'D': 
						# just to complete the other letter
						if bim_dico[nb][1] == '0' :
							allele1_buf += 1
							if bim_dico[nb][2] == 'I':
								allele1_new = 'D'
							else:
								allele1_new = 'I'
							allele2_buf += 1
							write_info(info_file, bim_dico[nb][0], bim_dico[nb][1], bim_dico[nb][2], allele1_new, bim_dico[nb][2])
						elif bim_dico[nb][2] == '0':
							allele2_buf += 1
							if bim_dico[nb][1] == 'I':
								allele2_new = 'D'
							else:
								allele2_new = 'I'
							allele1_buf += 1
							write_info(info_file, bim_dico[nb][0], bim_dico[nb][1], bim_dico[nb][2], bim_dico[nb][1], allele2_new)
						else :
							allele1_buf += 1
							allele2_buf += 1
					# to change the allele nomenclature +/-
					elif bim_dico[nb][1] == '+' or bim_dico[nb][2] == '+' or bim_dico[nb][1] == '-' or bim_dico[nb][2]== '-':
						allele_know = bim_dico[nb][1]
						if allele_know == '0':
							allele_know = bim_dico[nb][2]
							allele2_new, allele1_new, liste_allele_error = replace_alleles_PM(allele_know, bim_dico[nb][1], bim_dico[nb][0], liste_allele_error)
						else:
							allele1_new, allele2_new, liste_allele_error = replace_alleles_PM(allele_know, bim_dico[nb][1], bim_dico[nb][0], liste_allele_error)
						allele1_buf += 1
						allele2_buf += 1
						write_info(info_file, bim_dico[nb][0], bim_dico[nb][1], bim_dico[nb][2], allele1_new, allele2_new)

					# if allele B, nomenclature A/B
					elif  bim_dico[nb][1] == 'B': # !!! utile que lorsque le fichier n'est pas completement en notation A/B -> si on peut avoir B0 on peut avoir A0, les A0 ne seront pas detecte comme appartenant a la notation A/B, mais comme le A de ATCG.
						allele1_buf += 1
						allele2_buf += 1
						write_info(info_file, bim_dico[nb][0], bim_dico[nb][1], bim_dico[nb][2], alleles[snp][1], alleles[snp][0])
					elif  bim_dico[nb][2] == 'B':
						allele1_buf += 1
						allele2_buf += 1
						write_info(info_file, bim_dico[nb][0], bim_dico[nb][1], bim_dico[nb][2], alleles[snp][0], alleles[snp][1])
					else: # if the nomenclature is not known, the SNP is add to the list and we don't modify it
						add_liste_pb(liste_allele_pb, bim_dico[nb][0])
				
				# if both alleles are missing
				elif bim_dico[nb][1] == '0' and bim_dico[nb][2] == '0':
					allele1_new = alleles[snp][0]
					allele2_new = alleles[snp][1]
					if allele1_new != '0':
						allele1_buf += 1
					if allele2_new != '0':
						allele2_buf += 1
					
					if allele1_new != '0' or allele2_new != '0':
						write_info(info_file, bim_dico[nb][0] , 'Z', bim_dico[nb][2], allele1_new, allele2_new)

				# if one allele is missing
				elif bim_dico[nb][1] == '0' or bim_dico[nb][2] == '0':
					allele_know = bim_dico[nb][1]
					if allele_know == '0':
						allele_know = bim_dico[nb][2]
						allele2_buf += 1
						allele1_new, liste_allele_error = find_other_allele(alleles[snp][0], alleles[snp][1], bim_dico[nb][0], allele_know, liste_allele_error)
						if allele1_new != '0': #the allele we got back is known
							write_info(info_file, bim_dico[nb][0], bim_dico[nb][1], bim_dico[nb][2], allele1_new, bim_dico[nb][2])
							allele1_buf += 1
					else:
						allele1_buf += 1
						allele2_new, liste_allele_error = find_other_allele(alleles[snp][0], alleles[snp][1], bim_dico[nb][0], allele_know, liste_allele_error)
						if allele2_new != '0':
							write_info(info_file, bim_dico[nb][0], bim_dico[nb][1], bim_dico[nb][2], bim_dico[nb][1], allele2_new)
							allele2_buf += 1
				else: # both alleles are known in the bim file
					allele1_buf += 1
					allele2_buf += 1
	

#we calculate the % of 1st alleles that are known (= allele != 0)
if allele1_test == False: #there were some missing 1st alleles in the bim file
	allele1_found = allele1_buf / float(file_bim_length) * 100
	nb_allele1 = round(allele1_found, 4) #rounded number of 1st alleles find to the nearest 10th
else: #	allele1_test == True, all the 1st alleles were present in the bim file
	nb_allele1 = 100

#we calculate the % of 2nd alleles that are known (= allele != 0)
if allele2_test == False: #there were some missing 2nd alleles in the bim file
	allele2_found = allele2_buf / float(file_bim_length) * 100
	nb_allele2 = round(allele2_found, 4) #rounded number of 2nd alleles find to the nearest 10th
else: #allele2_test == True, all the 2nd alleles were present in the bim file
	nb_allele2 = 100


sys.stdout.write('\n\n**********************************************************\n\nResult for the allele research of the bim file: ' + bim_file + '\n\n\n')
	
#sys.stdout.write("length : "+str(file_bim_length)+"\n")
#sys.stdout.write("Allele 1 : "+str(allele1_buf)+"\n")
#sys.stdout.write("Allele 2 : "+str(allele2_buf)+"\n\n")

sys.stdout.write(str(nb_allele1) + '% of first alleles are known.\n')
sys.stdout.write(str(nb_allele2) + '% of second alleles are known.\n')

#sys.stdout.write(str(ambiguous_SNP) + ' new ambiguous SNP after completed the bim file alleles.\n')

sys.stdout.write('\n\nThe first SNPs which have an allele nomenclature not handle by the script:')
if liste_allele_pb != [] :
	for index in range(len(liste_allele_pb)):
		sys.stdout.write('\n' + liste_allele_pb[index])
else:
	sys.stdout.write('\nNone')

sys.stdout.write('\n\nThe first SNPs which have incoherent alleles with the manifest:')
if liste_allele_error != [] :
	for index in range(len(liste_allele_error)):
		sys.stdout.write('\n' + liste_allele_error[index])
else:
	sys.stdout.write('\nNone')

sys.stdout.write('\n\n**********************************************************\n\n')
