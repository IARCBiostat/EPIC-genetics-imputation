#!/usr/bin/env python2.7

#script that modify the position of a list of SNP by adding +1 at their position.

#input: 4 parameters:
# a bim file (mandatory: -b or --bim)
# the list of SNP name you need to modify. <cancer>_SNP_build_....txt (mandatory: -l or --list)
# a loci to rsID file (Illumina txt) (optionnal: -t or --txt)
# the number you want to add at your SNP position -d or --diff

#output:
# - [file]_goodSubversionPos.txt

#script written by Manon Knuchel in August 2019

##MODULES
import argparse #module to make it easy to write user-friendly command-line interfaces
#import pprint #module to print 'pretty' dictionnary
import re #module that allows to do regular expressions (= regex) search
#import time #module to see how long the script takes
import sys # module to write in the error output

#start = time.time()

##############################Input file options##############################

#you can access to the help section with '-h' option

parser = argparse.ArgumentParser()

parser.add_argument('--bim', '-b', type=str, help='Path to the bim file.', required = True)
parser.add_argument('--txt', '-t', type=str, help='Path to the file with the SNP list you want to modify the positions (<cancer>_SNP_build....txt).', required = True)
parser.add_argument('--diff', '-d', type=str, help='build version of the bim/manifest file', required = True)

args = parser.parse_args()

bim_file = args.bim

fic_SNP_list = args.txt
cancer = fic_SNP_list.split('/')[-1]
cancer = cancer.replace(".txt", "", 1) #the name of the file without the extension (.txt)


pos_modiff = int(args.diff)

##############################Storage of the SNP list##############################

SNP_list =[] # list of SNP which we want to modify the position

with open(fic_SNP_list) as list_file :
	for line in list_file :
		SNP_list.append(line.strip())


##############################Reading the Bim file and Creation of the output file##############################

with open(bim_file) as bim, open(cancer+"_goodSubversionPos.txt", 'w') as out_file:
	for line in bim:
		fields=line.split("\t")

		name=fields[1].strip() #the SNP name
		pos=int(fields[3].strip()) #the position

		if name in SNP_list :
			new_pos = pos + pos_modiff
			out_file.write(name +"\t" + str(new_pos) + "\n") #we write the info in the file containing the new pos info


