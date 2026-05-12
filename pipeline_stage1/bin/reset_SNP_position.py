#!/usr/bin/env python2.7

# this script delete all position information of SNPs in the Bim file

#input: 1 parameter:
# the bim file (without any extension)

#output: Bim file without position : file + _reset_all_pos

#script written by Manon Knuchel in August 2019

##MODULES
import sys, os, argparse

############################ getting arguments given to the script #######################################

parser = argparse.ArgumentParser()

parser.add_argument('--file', '-f', type=str, help='Path to the data files.', required = True)
args = parser.parse_args()

file_plink = args.file
plink = os.environ.get("PLINK_BIN", "plink")

with open(file_plink+".bim") as bim_file, open(file_plink+"_reset_chr.txt", 'w') as out_chr, open(file_plink+"_reset_pos.txt", 'w') as out_pos:
	for line in bim_file:
		fields=line.split("\t")

		#get all the information needed
		name=fields[1].strip() #the SNP name

		out_chr.write(name + "\t0\n") # we write "0" to have missing information
		out_pos.write(name + "\t0\n")

# we use plink to update chromosome and position with the missing data
os.system(plink + ' --bfile ' + file_plink + ' --update-chr ' + file_plink + '_reset_chr.txt --make-bed --out ' + file_plink + '_reset_chr')
os.system(plink + ' --bfile ' + file_plink + '_reset_chr --update-map ' + file_plink + '_reset_pos.txt --make-bed --out ' + file_plink + '_reset_all_pos')

os.system("rm " + file_plink + "_reset_chr.* " + file_plink + "_reset_pos.txt")
