#!/usr/bin/env python2.7

#this script rename the SNPs in the bim file with the name of SNPs in the Manifest file

#input: 2 parameters:
# the bim file (without any extension)
# the standard Manifest file (.csv)

#output: [project_file]_ManifestNames.txt

#script written by Manon Knuchel in August 2019

##MODULES
import sys, os, re, datetime, argparse, csv

############################ getting arguments given to the script #######################################

parser = argparse.ArgumentParser()

parser.add_argument('--file', '-f', type=str, help='Path to the data files.', required = True)
parser.add_argument('--man', '-m', type=str, help='Path to the manifest file (Illumina.csv).', required = True)
parser.add_argument('--txt', '-t', type=str, help='Path to the loci to rsID file (Illumina.txt).')

args = parser.parse_args()

bim = args.file+'.bim'
cancer = args.file.split('/')[-1]

manifest_file = args.man
file_rsID = args.txt

########### we read the Manifest file and keep information of its SNPs

dico_man = {} # keep data from the Manifest file

with open(manifest_file) as man_file:
	reader = csv.reader(man_file, delimiter='\t')
	for line_nb, row in enumerate(reader):
		#we get all the data we need
		snpName = row[2].strip()
		allele_1 = row[3].strip()
		allele_2 = row[4].strip()

		if (allele_1 == 'N'or allele_2 == 'N') and (allele_1 == allele_2): # if we have a SNP with the Alleles NA, we put them as unknown
			allele_1 = '0'
			allele_2 = '0'

		dico_man[snpName]=(allele_1, allele_2)


snp_name = {} #dictionnary that will be like: {'name Manifest file': 'other name of the SNP(might be used in the Bim file'}

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


test1=0
test2=0
test3=0
with open(bim) as bim_file, open(cancer+"_alleles_tradAB.txt", 'w') as out_file:
	for line in bim_file:
		sep = re.search(r'(\t)+',line) #try if the separators are tabulations
		if sep: # separators are tabulations, sep != None
			fields=line.split("\t")
		else:
			fields=line.split(" ")

		#get all the information needed
		name=fields[1].strip() #the SNP name
		allele1=fields[4].strip() #the first allele
		allele2=fields[5].strip() #the second allele

		snp = ""
		if name in dico_man:
			snp = name
		elif file_rsID != None and name in snp_name :
			if snp_name[name] in dico_man :
				snp = snp_name[name]
		
		if snp != "" :
			if allele1 != "0" or allele2 != "0": # if one of the allele is known
				if allele1 == "A" or allele2 == "B" :
					test1 += 1
					out_file.write(name + "\t" + allele1 + "\t" + allele2 + "\t" + dico_man[snp][0] + "\t" + dico_man[snp][1] + "\n")
				elif allele1 == "B" or allele2 == "A" :
					test1 += 1
					out_file.write(name + "\t" + allele1 + "\t" + allele2 + "\t" + dico_man[snp][1] + "\t" + dico_man[snp][0] + "\n")
				else : # if allele are not A or B
					test3 += 1
					sys.stderr("\n\n**********************************************************\n\nwhat allele have you ? '"+allele1+" "+allele2+"'\n\n**********************************************************\n\n")
			else : # if both allele are unknown
				test2 += 1
				out_file.write(name + "\tZ\t" + allele2 + "\t" + dico_man[snp][0] + "\t" + dico_man[snp][1] + "\n")

		

sys.stdout.write("\n\n**********************************************************\n\nNumber of SNPs with alleles translated in A, C, G, T nomenclature: "+str(test1)+"\n")
sys.stdout.write("Number of allele with both alleles missing that have been completed: "+str(test2)+"\n")
sys.stdout.write("!!!!!!!!!!!! Warning : Number of SNPs without A/B nomenclature: "+str(test3)+"\n\n**********************************************************\n\n")
