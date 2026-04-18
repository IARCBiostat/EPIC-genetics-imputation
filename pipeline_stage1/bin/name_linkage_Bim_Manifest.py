#!/usr/bin/env python2.7

#this script rename the SNPs in the bim file with the name of SNPs in the Manifest file

#input: 2 parameters:
# the bim file (without any extension)
# the standard Manifest file (.csv)

#output: [project_file]_ManifestNames.txt

#script written by Manon Knuchel in August 2019

##MODULES
import sys, os, re, datetime, argparse, csv

#declaring complementary allele dictionary, it's also used to verify if the alleles are ATCG
complementary = {"A":"T", "T":"A", "C":"G", "G":"C"}

############################ getting arguments given to the script #######################################

parser = argparse.ArgumentParser()

parser.add_argument('--file', '-f', type=str, help='Path to the data files.', required = True)
parser.add_argument('--man', '-m', type=str, help='Path to the manifest file (Illumina.csv).', required = True)

args = parser.parse_args()

bim = args.file+'.bim'
cancer = args.file.split('/')[-1]

manifest_file = args.man

########### we read the Manifest file and keep information of its SNPs

dico_chrpos_man = {} # keep data from the Manifest file

with open(manifest_file) as man_file:
	reader = csv.reader(man_file, delimiter='\t')
	for line_nb, row in enumerate(reader):
		#we get all the data we need
		chro = row[0].strip()

		if chro == "M" or chro == "MT": # convertion of the letter in the number notation
			chro = "26"
		elif chro == "XY":
			chro = "25"
		elif chro == "Y":
			chro = "24"
		elif chro == "X":
			chro = "23"
			
		pos = row[1].strip()

		chrpos = chro+":"+pos # the position of the SNP will be the key of the dictionary

		snpName = row[2].strip()
		allele_1 = row[3].strip()
		allele_2 = row[4].strip()

		if allele_1 == 'N'or allele_2 == 'N': # if we have a SNP with the Alleles NA, we put them as unknown
			allele_1 = '0'
			allele_2 = '0'

		alleles = allele_1 + allele_2

		if chro != '0' and pos != '0': # if one of the position is missing we might do a wrong link between the two files
			if chrpos not in dico_chrpos_man :
				dico_chrpos_man[chrpos] = []

			curent_snp = (alleles, snpName)
			dico_chrpos_man[chrpos].append(curent_snp) # we have for each position a list of (alleles, name) in case there is more than one SNP at this position
		#sys.stdout.write(chrpos + "\n")

		

########### read the bim file and write result in the output file #############
link_found = 0

#test1 = 0

#test2 = 0

with open(bim) as bim_file, open(cancer+"_linkManifestNames.txt", 'w') as out_file:
	for line in bim_file:
		sep = re.search(r'(\t)+',line) #try if the separators are tabulations
		if sep: # separators are tabulations, sep != None
			fields=line.split("\t")
		else:
			fields=line.split(" ")
		
		#get all the information needed
		chro=fields[0].strip() #the chromosome
		name=fields[1].strip() #the SNP name
		pos=fields[3].strip() #the position
		allele1=fields[4].strip() #the first allele
		allele2=fields[5].strip() #the second allele
		chrpos = chro +":"+ pos

		
		if allele1 != "0" and allele2 != "0":
			#test1 += 1
			alleles12 = allele1 + allele2
			alleles21 = allele2 + allele1
			Calleles12 = complementary.get(allele1) + complementary.get(allele2)
			Calleles21 = complementary.get(allele2) + complementary.get(allele1)
			
			good_alleles = 0

			if chrpos in dico_chrpos_man :
				#test2 += 1
				for snp in dico_chrpos_man[chrpos]:
					if snp[0] == alleles12 or snp[0] == alleles21 :
						if good_alleles == 0 :
							out_file.write(name + "\t" + snp[1] + "\n")
							good_alleles = 1
							link_found += 1
					elif snp[0] == Calleles12 or snp[0] == Calleles21 :
						if good_alleles == 0 :
							out_file.write(name + "\t" + snp[1] + "\n")
							good_alleles = 1
							link_found += 1
					#else : # no good alleles

				if good_alleles == 0 :
					sys.stderr.write("\n\n**********************************************************\n\n" + name + "don't correspond to :\n")
					for snp in dico_chrpos_man[chrpos]:
						sys.stderr.write(snp[1] + " ")
					sys.stderr.write("\n\n**********************************************************\n\n")
			#else :
				#sys.stdout.write(chrpos+"\n")





sys.stdout.write("\n\n**********************************************************\n\nNumber of link found between SNPs of the Bim file and those from the Manifest file : "+str(link_found)+"\n\n**********************************************************\n\n")

#sys.stdout.write("\n\ntest1 : " + str(test1) + "\n")
#sys.stdout.write("\n\ntest2 : " + str(test2)+ "\n")


