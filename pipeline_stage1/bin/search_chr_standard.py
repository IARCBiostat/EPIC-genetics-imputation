#!/usr/bin/env python2.7

#script that finds the chromosome for a SNP, if the manifest file in csv is available
#input: 3 parameters:
# a bim file (mandatory: -b or --bim)
# a manifest file (Illumina csv) (mandatory: -m or --man)
# a loci to rsID file (Illumina txt) (optionnal: -t or --txt)
#output:
# - [file]_goodChr_[genomeVersion].txt
# - [file]_goodChr_[genome_version]_results.txt'

#script written by Manon Knuchel in May 2019 based on the first version developped in May 2018 by Emilie Gerard-Marchant & Benjamin Bourgeois

##MODULES
import argparse #module to make it easy to write user-friendly command-line interfaces
import csv #module to read csv
#import pprint #module to print 'pretty' dictionnary
import re #module that allows to do regular expressions (= regex) search
#import time #module to see how long the script takes
import sys # module to write in the error output

#start = time.time()

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

bim_dico = {} #dictionnary where all the datas from the bim file will be stored, it will be like : {line_nb: (chr, 'snp')}

file_bim_length = 0 #variable that will contain the number of lines of the old bim file

chr_test = True #variable that will tell us about the presence of the chr info in the bim file: True if present, False if missing


with open(bim_file) as file_bim: #bim_file  is the bim file given in the terminal
	for line_nb, line in enumerate(file_bim):
		sep = re.search(r'(\t)+',line) #try if the separators are tabulations

		if sep: #separator are tabulations, sep != None
			data = line.split('\t')
		else: #sep == None : separators are spaces
			data = line.split(' ')

		#get all the datas we need
		snp = data[1]
		chromo = int(data[0])

		if chromo == 0: #the chr info is missing
			chr_test = False

		bim_dico[line_nb] = (chromo, snp) #creation of the dictionnary, key = int, value = tuple
		file_bim_length = line_nb + 1 # = number of SNPs in the bim file

#pprint.pprint(bim_dico) #to have a 'pretty' print of the dictionary

##############################Storage of the manifest file datas##############################

#col_names = ['Chr','MapInfo''Name','Allele1','Allele2','strand'] #the name of the columns in the csv file

file_illumina_length = 0  #variable that will contain the number of lines in the Illumina file (csv)

chr_MapInfo = {} #dictionnary that will be like {'snp': chr}

with open(manifest_file) as illumina_csv: #manifest_file is the manifest file given in the terminal
	for line in illumina_csv:
		file_illumina_length += 1 #we count the nb of lines

#       print(file_illumina_length)

with open(manifest_file) as illumina_csv:
	reader = csv.reader(illumina_csv, delimiter='\t')
	for line_nb, row in enumerate(reader):
		#we get all the datas we want
		snp = row[2]
		chromo = row[0]

		if chromo == 'X': #the X chr
			chromo = 23 #we change the nomenclature
		elif chromo == 'Y': #the Y chr
			chromo = 24
		elif chromo == 'XY': #the pseudo-autosomal region (PAR)
			chromo = 25
		elif chromo == 'MT' or chromo == 'M' or chromo == 'chrM': #the chr mitochondrial
			chromo = 26

		chr_MapInfo[snp] = int(chromo) #construction of the dictionnary, key = str, value = int

#pprint.pprint(chr_MapInfo)

##############################Storage of the loci to rsID file datas##############################

snp_name = {} #dictionnary that will be like: {'kgp...': 'rs...'}

if file_rsID: #the loci to rsId file was given in the terminal
	with open(file_rsID) as loci_file: #file_rsID is the loci to rsID file given in the terminal
		for line_nb, line in enumerate(loci_file):

			if line_nb >= 1: #because the 1st line is the name of the columns
				data = line.split('\t')
				#we get all the datas we need
				new_snp = data[0] #the SNP that is not necessarily in the 'rs[number]' nomenclature
				good_snp = data[1].strip() #the SNP that is in the 'rs[number]' nomenclature
                                if good_snp != '.': #there is a name for the SNP
					if ',' in good_snp: #there was a merge between several SNPs
						snps = []
						snps = good_snp.split(',')
						for var in snps:
							snp_name[var] = new_snp
					else:
						snp_name[good_snp] = new_snp ##construction of the dictionnary,key = str, value = str


#pprint.pprint(snp_name)

##############################Creation of the result file##############################

def change_name(dico, snp):
	"""Function that gets the old name of a SNP if the new one is missing in the manifest file
	"""
	new_name = dico[snp]
	return new_name

def get_chr(dico, snp):
	"""Function that gets the chr if the SNP is present in the manifest file""
	"""
	chromo = dico[snp]
	return chromo

def add_liste_pb(liste, snp):
	"""Function that adds the SNP name in the problem list: if the chr information is missing
	"""
	if snp not in liste:
		liste.append(snp)


cancer_bim = bim_file_name.split('.')[0] #the name of the bim file without the extension (.bim)
file_name = cancer_bim + '_goodChr_' + str(genome_version) + '.txt' #the name of the file that will contain the new chr info

chr_buf = 0 #variable that will contain the number of chr != 0 at the end of this script

liste_pb_chr = [] #list that will contain all the SNPs which have missing chr info


with open(file_name, 'w') as info_file:
	if chr_test == False: #there were missing chr in the bim file
		for nb in range(file_bim_length):

			#the chr is missing
			if bim_dico[nb][0] == 0: 

				if bim_dico[nb][1] in chr_MapInfo: #the SNP is in the manifest file
					chromo = get_chr(chr_MapInfo, bim_dico[nb][1]) #we get the new chr
					if chromo != 0:
						chr_buf += 1
						info_file.write(bim_dico[nb][1] + "\t" + str(chromo) + "\n") #we write the info in the file containing the new chr info
				elif file_rsID != None: #the SNP is not in the manifest file and the loci to rsID was given in the terminal
					if bim_dico[nb][1] in snp_name: #the SNP is in the loci to rsID file
						snp = change_name(snp_name,bim_dico[nb][1]) #we get the new name

						if snp in chr_MapInfo: #the new name is in the manifest file
							chromo = get_chr(chr_MapInfo, snp) #we get the new chr
							if chromo != 0:
								chr_buf += 1
								info_file.write(bim_dico[nb][1] + "\t" + str(chromo) + "\n") #we write the info in the file containing the new chr info
						else: #the new name is not in the manifest file
							add_liste_pb(liste_pb_chr, bim_dico[nb][1])
					else: #the SNP is not in the loci to rsID file
						add_liste_pb(liste_pb_chr, bim_dico[nb][1])

				else: #the SNP is neither in the manifest file nor in the loci to rsID
					add_liste_pb(liste_pb_chr, bim_dico[nb][1])

			#the chr is not missing
			else: 
				chr_buf += 1


#we calculate the % of chr that are known (= chr != 0)
if chr_test == False: #there were some missing chr in the bim file
	chr_found = chr_buf / float(file_bim_length) * 100
	nb_chr = round(chr_found, 4) #rounded number of chr find to the nearest 10th
else: #chr_test == True, all the chr were present in the bim file
	nb_chr = 100


sys.stdout.write('\n\n**********************************************************\n\nResult for the chromosome research of the bim file: ' + bim_file_name + '\n\n\n')
sys.stdout.write(str(nb_chr) + '% of chromosomes are known.\n')


if liste_pb_chr != []: #there are some SNPs with missing chr info
	sys.stdout.write('\n\nThe first SNPs which still have missing chromosomes:')

for nb_line, snp in enumerate(liste_pb_chr):
	sys.stdout.write('\n' + snp)
	if nb_line == 10:
		break
		
sys.stdout.write('\n\n**********************************************************\n\n')
#end = time.time()

#print(end - start)
