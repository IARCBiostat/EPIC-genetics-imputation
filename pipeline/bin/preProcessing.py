#!/usr/bin/env python2.7

#this script  do an inventory of the data from the bim and fam files
#input: 2 parameters:
# the project file (without any extension)
# the manifest file (.csv)
#output: [project_file]_results.txt


#script written by Manon Knuchel in May 2019 based on the first version developed in May 2018 by Emilie Gerard-Marchant & Benjamin Bourgeois

plink = "plink"


##MODULES
import sys, os, re, datetime, argparse

#My module:
#from preProcessing_functions.SNP_research_BM import * #importing all function of the module SNP_research_BM who is in the folder pre-processing_functions
from preProcessing_functions.checks_listing import * #importing all function of the module checks_listing
from preProcessing_functions.strand_determination import * #research_strand()

############################ getting arguments given to the script #######################################

parser = argparse.ArgumentParser()

parser.add_argument('--file', '-f', type=str, help='Path to the data files.', required = True)
parser.add_argument('--man', '-m', type=str, help='Path to the manifest file (Illumina.csv).', required = True)
parser.add_argument('--txt', '-t', type=str, help='Path to the loci to rsID file (Illumina.txt).')
parser.add_argument('--chr', '-c', type=int, help='Field number (start at 0) of SNPs Chromosome in the manifest file.', required = True)
parser.add_argument('--pos', '-p', type=int, help='Field number (start at 0) of SNPs position in the manifest file.', required = True)
parser.add_argument('--name', '-n', type=int, help='Field number (start at 0) of SNPs name in the manifest file.', required = True)
parser.add_argument('--allele', '-a', type=int, help='Field number (start at 0) of SNPs alleles in the manifest file.', required = True)

parser.add_argument('--strandPlus', '-sp', type=int, help='Field number (start at 0) of SNPs strand (+/- nomenclature) in the manifest file.')
parser.add_argument('--strandForward', '-sf', type=int, help='Field number (start at 0) of the IlmnID field (contain For/Rev information) in the manifest file.')
parser.add_argument('--strandTop', '-st', type=int, help='Field number (start at 0) of SNPs strand (Top/Bot nomenclature) in the manifest file.')

parser.add_argument('--build', '-b', type=int, help='Field number (start at 0) of SNPs build in the manifest file.')

parser.add_argument('--nlhead', '-nlh', type=int, help='Number of line without SNP at the begining of the manifest file.')
parser.add_argument('--nltail', '-nlt', type=int, help='Number of line without SNP at the end of the manifest file.')

parser.add_argument('--summary_file', '-sum_f', type=str, help='Path to the Normalization summary files.')
parser.add_argument('--summary_number', '-sum_n', choices=['a', 'n', '1', '2', '3'], help='QC number : input data = 1 ; after completing data = 2 ; after expluding part2 = 3 ; if you do not want to print in the summery file = n ; if you want to print every thing = a.')

args = parser.parse_args()

bed = args.file+'.bed'
bim = args.file+'.bim'
fam = args.file+'.fam'
cancer = args.file.split('/')[-1]

summary_file = args.summary_file

print_summary = args.summary_number
if print_summary == None:
	print_summary = 'a'

if summary_file == None:
	print_summary = 'n'

manifest_file = args.man
manifest_name = manifest_file.split('/')[-1]

file_rsID = args.txt

field_chr = args.chr
field_pos = args.pos
field_name = args.name
field_allele = args.allele

field_strandP = args.strandPlus
field_strandF = args.strandForward
field_strandT = args.strandTop

if field_strandP == None and field_strandF == None and field_strandT == None:
	sys.stderr.write("\n\n**********************************************************\n\nERROR: you need to specify at least one of these parameters:\n-sp : the field number (start at 0) of SNPs strand (+/- nomenclature) in the manifest file.\n-sf : Field number (start at 0) of the IlmnID field (contain For/Rev information) in the manifest file.\n-st : Field number (start at 0) of SNPs strand (Top/Bot nomenclature) in the manifest file.\n\n**********************************************************\n\n")
	exit(0)
field_build = args.build

nlhead = args.nlhead
if nlhead == None:
	nlhead = 0

nltail = args.nltail
if nltail == None:
	nltail = 0

################## creating output files ############################
date = datetime.datetime.now()

#we create a new file with all results
FSOR=open(cancer+"_results.txt","w")

#we create a new file for the SNPs in the strand -
if field_strandP != None:
	FSSTRANDp=open(cancer+"_negative_strand_manifest_PM.txt", 'w')
	FSMINUS=open(cancer+"_negative_strand_PM.txt", 'w')
	FSUNSTRANDp=open(cancer+"_unknown_strand_PM.txt", 'w') # file with SNP with unknown strand in the Manifest
#we create a new file for the SNPs in the strand Rev
if field_strandF != None:
	FSSTRANDf=open(cancer+"_negative_strand_manifest_FR.txt", 'w')
	FSREV=open(cancer+"_negative_strand_FR.txt", 'w')
	FSUNSTRANDf=open(cancer+"_unknown_strand_FR.txt", 'w') # file with SNP with unknown strand in the Manifest
#we create a new file for the SNPs in the strand Bot
if field_strandT != None:
	FSSTRANDt=open(cancer+"_negative_strand_manifest_TB.txt", 'w')
	FSBOT=open(cancer+"_negative_strand_TB.txt", 'w')
	FSUNSTRANDt=open(cancer+"_unknown_strand_TB.txt", 'w') # file with SNP with unknown strand in the Manifest
#we create a new file with only needed fields of the manifest
FSMAN=open(cancer+"_manifest.csv", 'w')
#we create a new file with SNP which have not the same chr:pos than the manifest
FSDIFPOS=open(cancer+"_diff_chrpos.txt", 'w')
#we create a new file with SNP which have not the same alleles than the manifest
FSDIFALL=open(cancer+"_diff_alleles.txt", 'w')
#file with SNP which are not found in the manifest file
FSNOTFOUND=open(cancer+"_not_found.txt", 'w')
#file with misplaced SNP with the name Mito
FSMITO=open(cancer+"_missplaced_mito.txt", 'w')
#file with duplicate SNP... (at least same position but can also have the same name and they have the same alleles)
FSDUPLSAME=open(cancer+"_duplicate_T1_same_allele.txt", 'w')
#file with duplicate SNP, tri-allelic (same name and position but differents alleles)
FSDUPLNAME=open(cancer+"_duplicate_T2_name_pos.txt", 'w')
#file with duplicate SNP, tri-allelic (different names, same positions and differents alleles)
FSDUPLPOS=open(cancer+"_duplicate_T3_pos.txt", 'w')
#file with duplicate SNP ? (same names but different positions. They can have the same alleles or differents alleles)
FSDUPLERROR=open(cancer+"_duplicate_T4_error.txt", 'w')


#FSOR.write(str(date)+"\n")
FSOR.write("bed/bim/fam data used: "+cancer+"\n")
FSOR.write("Manifest file used: "+manifest_name+"\n")


#declaring complementary allele dictionary, it's also used to verify if the allele are ATCG
complementary = {"A":"T", "T":"A", "C":"G", "G":"C"}

# nb1 et nb2 correspond to the length of error lists (if you want to modify it see in th list declaration below, wich variable use your list)
nb1=10
nb2=30


#################### declaration of all error list ############################

liste_wrong_mito_bim = [] # nb1 : list 10 first SNPs which are not on chr 0 nor chr 26
liste_wrong_mito_manifest = [] # nb1 : list 10 first SNPs from the manifest file which are not on chr 0 nor chr 26
liste_notrs = [] # nb2 : list 30 first name which do not start with Mito nor rs
liste_diff_alleles = [] # nb1 : list of SNPs with at least one alleles different from A T C G 0
liste_diff_alleles_manifest = [] # nb1 : list of SNPs from the manifest file with at least one alleles different from A T C G 0
liste_diff_bim_manifest_p = [] # nb1 : list incoherent alleles between the bim and manifest file
liste_diff_bim_manifest_f = [] # nb1 : list incoherent alleles between the bim and manifest file
liste_diff_bim_manifest_t = [] # nb1 : list incoherent alleles between the bim and manifest file

liste_double_rs_bim = [] # nb1 : list first SNP with same name in bim file
liste_double_chrpos_bim = [] # nb1 : list first SNP with same chr and pos in bim file
liste_double_rschrpos_bim = [] # nb1 : list first SNP with same name, chr and pos in bim file

liste_double_rs_man_p = []# nb1 : list first SNP with same name in manifest file
liste_double_chrpos_man_p = [] # nb1 : list first SNP with same chr and pos in manifest file
liste_double_rschrpos_man_p = []# nb1 : list first SNP with same name, chr and pos in manifest file

liste_double_rs_man_f = []# nb1 : list first SNP with same name in manifest file
liste_double_chrpos_man_f = [] # nb1 : list first SNP with same chr and pos in manifest file
liste_double_rschrpos_man_f = []# nb1 : list first SNP with same name, chr and pos in manifest file

liste_double_rs_man_t = []# nb1 : list first SNP with same name in manifest file
liste_double_chrpos_man_t = [] # nb1 : list first SNP with same chr and pos in manifest file
liste_double_rschrpos_man_t = []# nb1 : list first SNP with same name, chr and pos in manifest file

liste_diff_pos_BM = [] # nb1 : list first SNP with different position in the bim and manifest file

##########################################################################
######################## Start the sample checks #########################
##########################################################################
FSOR.write(str("\n\n********************************************************** \n"))
FSOR.write(str("******************** Sample Checks *********************** \n"))
FSOR.write(str("********************************************************** \n \n"))

#calculation number of subjects
nb_people=file_len(fam)
FSOR.write(str("These file contains ")+str(nb_people)+str(" subjects")+"\n") 

#calculating number of SNPs
nb_SNP=file_len(bim)
FSOR.write(str("These file contains ")+str(nb_SNP)+str(" SNPs")+"\n \n")



##########################################################################
####################### Start the file Check #############################
##########################################################################
FSOR.write(str("********************************************************** \n"))
FSOR.write(str("********************** Check bim file ******************** \n"))
FSOR.write(str("********************************************************** \n \n"))
FSOR.write(str("---------------- Check of chromosome information ---------------------- \n"))
f_bim=open(bim)#check the bim file

###################### declaring all variable use to count ############################

compteur_chro_0=0 #nb of chr = 0 (= unknown chr)
compteur_chro_23=0 #nb of chr = 23 (= X chr)
compteur_chro_24=0 #nb of chr = 24 (= Y chr)
compteur_chro_25=0  #nb of chr = 25 (= pseudo-autosomal chr)
compteur_chro_26 = 0 #nb of chr = 26 (= mitochondrial chr)
compteur_chro_sup_27 = 0 #nb of chr >= 27 (= anormal chr)
compteur_goodrs=0 #rs[number]
compteur_mito=0 #count name of SNP if they start with Mito (should be on chr26, but can be on 0)
compteur_mito_pos=0 #count number of Mitochondrial SNP which need to be mooved on chr26
compteur_mito_wrong=0 #count number of Mitochondrial SNP which are not placed on the good chr
compteur_badrs=0 #not rs[number]
compteur_goodpos=0 #pos != 0
compteur_badpos=0 #pos == 0
compteur_bad_chrpos=0 #chro=0 and pos=0
chro_pasbeau=0
compteur_missing_pos = 0 # chro=0 or pos=0

couple_good_allele = 0 # both allele are A,T,C or G
couple_unusual_allele = 0 # at least one allele is not A, T ,C, G or 0
couple_missing_allele = 0 # at least one allele is 0, if only one missing the other is A,T,C or G

couple_good_allele_manifest = 0 # in the manifest, both allele are A,T,C or G
couple_unusual_allele_manifest = 0 # in the manifest, at least one allele is not A, T ,C, G or 0
couple_missing_allele_manifest = 0 # in the manifest, at least one allele is 0, if only one missing the other is A,T,C or G

man_minus_strand_p = 0 # Number of SNP on the strand - in the manifest
man_minus_strand_f = 0 # Number of SNP on the strand Rev in the manifest
man_minus_strand_t = 0 # Number of SNP on the strand Bot in the manifest

compteur_mito_wrong_manifest = 0 #count number of Mitochondrial SNP which are not placed on the good chr in the manifest file

double_rs_bim = 0 # count number of SNP (name) witch are at least twice in the bim file
double_chrpos_bim = 0 # count number of SNP (chr and pos) witch are at least twice in the bim file
double_rschrpos_bim = 0 # count number of SNP (name, chr and pos) witch are at least twice in the bim file

double_rschrpos_man_p = 0 # count number of SNP (name, chr and pos) witch are at least twice in the manifest
double_rs_man_p = 0 # count number of SNP (name) witch are at least twice in the manifest
double_chrpos_man_p = 0 # count number of SNP (chr and pos) witch are at least twice in the manifest

double_rschrpos_man_f = 0 # count number of SNP (name, chr and pos) witch are at least twice in the manifest
double_rs_man_f = 0 # count number of SNP (name) witch are at least twice in the manifest
double_chrpos_man_f = 0 # count number of SNP (chr and pos) witch are at least twice in the manifest

double_rschrpos_man_t = 0 # count number of SNP (name, chr and pos) witch are at least twice in the manifest
double_rs_man_t = 0 # count number of SNP (name) witch are at least twice in the manifest
double_chrpos_man_t = 0 # count number of SNP (chr and pos) witch are at least twice in the manifest

double_same_allele_bim = 0 # count the number of SNP (name, chr and pos) witch are at least twice in the bim file and have identical alleles

double_same_allele_manifest_p=0 # count the number of SNP (name, chr and pos) witch are at least twice in the manifest file and have identical alleles
double_same_allele_manifest_f=0# count the number of SNP (name, chr and pos) witch are at least twice in the manifest file and have identical alleles
double_same_allele_manifest_t=0# count the number of SNP (name, chr and pos) witch are at least twice in the manifest file and have identical alleles


############### dictionary to count SNP in double in the bim file #######################

dico_bim_rs = {} # keep name of SNP in the bim file to count the number of SNP in double
dico_bim_chrpos = {} # keep chr and pos of SNP in the bim file to count the number of SNP in double

dico_count_lvl_name_duplicate = {} # count how many name duplicate of 2 SNP we have, 3 SNP ....
dico_count_lvl_pos_duplicate = {} # count how many position duplicate of 2 SNP we have, 3 SNP ....

############### declaring dictionnary which count the number of different alleles found in files ###############
compteur_allele1 = {"A":0,"C":0,"G":0,"T":0,"0":0}
compteur_allele2 = {"A":0,"C":0,"G":0,"T":0,"0":0}

compteur_allele1_manifest = {"A":0,"C":0,"G":0,"T":0,"0":0}
compteur_allele2_manifest = {"A":0,"C":0,"G":0,"T":0,"0":0}


##### function which count the different alleles found, and the SNP kind #####
def check_allele(rs, allele1, allele2, compteur_allele1, compteur_allele2, liste, couple_good_allele, couple_unusual_allele, couple_missing_allele):
	# count the number found and add new one in the dictionnary
	if allele1 in compteur_allele1:
		compteur_allele1[allele1] += 1
	else:
		compteur_allele1[allele1] = 1
	if allele2 in compteur_allele2:
		compteur_allele2[allele2] += 1
	else:
		compteur_allele2[allele2] = 1
	
	# count the number of SNP with different nomenclature
	if allele1 in complementary and allele2 in complementary: # in both alleles are ATCG
		couple_good_allele +=1
	elif (allele1 not in complementary and allele1!="0") or (allele2 not in complementary and allele2!="0"): # if at least one allele is different from ATCG and 0 (ex : ID, NA, ...)
		snp_unusual_alleles = (rs, allele1, allele2)
		add_liste_pb(liste, snp_unusual_alleles, nb1)
		couple_unusual_allele +=1
	elif allele1 == "0" or allele2 == "0": # if it have missing alleles (if one allele is know, it is in ATCG)
		couple_missing_allele +=1
	else : # this should never happen
		sys.stderr.write("\n\n**********************************************************\n\nA case is missing in the script pre-processing.py (line 288): "+ allele1 +" "+ allele2 +"\n\n**********************************************************\n\n")
		exit(0)
	return compteur_allele1, compteur_allele2, couple_good_allele, couple_unusual_allele, couple_missing_allele



####################################### reading the bim file (first verification)######################################
for line in f_bim:
	sep = re.search(r'(\t)+',line) #try if the separators are tabulations
	if sep: # separator are tabulations, sep != None
		fields=line.split("\t")
	else: #sep == None : separators are spaces
		fields=line.split(" ") 
	#get all the information needed
	chro=int(fields[0].strip()) #the chromosome
	rs=fields[1].strip() #the SNP name
	pos=fields[3].strip() #the position
	allele1=fields[4].strip() #the first allele
	allele2=fields[5].strip() #the second allele
	chrpos=str(chro) +":"+ pos
	alleles= allele1 + allele2

      	##we check the chromosome##
	if chro<1: #if the file contains unknow chromosome
		compteur_chro_0 += 1
	elif chro==23:
		compteur_chro_23 += 1
	elif chro==24:
		compteur_chro_24 += 1
	elif chro==25:
		compteur_chro_25 += 1
	elif chro==26:
		compteur_chro_26 += 1
	elif chro > 26:
		compteur_chro_sup_27 += 1
		chro_pasbeau=chro

        ##we check the rs##
	if rs.startswith("rs"): #check if we get rs+number
		compteur_goodrs +=1
	elif rs.startswith("Mito"):
		compteur_mito +=1
		if chro != 26:
			if chro==0:
				compteur_mito_pos +=1
			else:
				compteur_mito_wrong +=1
				add_liste_pb(liste_wrong_mito_bim, (rs,chro), nb1)
				FSMITO.write(rs+"\n")
	else: #if we don't get rs+number or Mito
		compteur_badrs +=1
		add_liste_pb(liste_notrs, rs, nb2)

	##we check the position##
	if pos =="0":
		compteur_badpos +=1
	else:
		compteur_goodpos +=1

	#if the SNP doesn't have position and chromosome
	if chro<1 and pos =="0":
		compteur_bad_chrpos +=1
	if chro=="0" or pos =="0":
		compteur_missing_pos +=1
	##we count the different alleles found and SNP kind##
	compteur_allele1, compteur_allele2, couple_good_allele, couple_unusual_allele, couple_missing_allele = check_allele(rs, allele1, allele2, compteur_allele1, compteur_allele2, liste_diff_alleles, couple_good_allele, couple_unusual_allele, couple_missing_allele)
	
	#count number of SNP in double in the bim file
	bim_duplicate_result = count_duplicate(dico_bim_rs, dico_bim_chrpos, rs, chrpos, allele1, allele2)

	if bim_duplicate_result != 0:
		if bim_duplicate_result == 1: # Duplicate Type 1
			double_same_allele_bim += 1
			if len(dico_bim_chrpos[chrpos])==1:# writting the first duplicate
				FSDUPLSAME.write(dico_bim_chrpos[chrpos][0]+"\n")
			FSDUPLSAME.write(rs+"\n")
		elif bim_duplicate_result == 2: # Duplicate Type 2
			double_rschrpos_bim += 1
			snp = (rs, chrpos, alleles, dico_bim_rs[rs][len(dico_bim_rs[rs])-1][1])
			add_liste_pb(liste_double_rschrpos_bim,snp,nb1)
			FSDUPLNAME.write(rs+"\n") # do we need to write the first duplicate ? they have the same name ....
		elif bim_duplicate_result == 3: # Duplicate Type 3
			double_chrpos_bim += 1
			snp = (rs, dico_bim_chrpos[chrpos][len(dico_bim_chrpos[chrpos])-1], chrpos, alleles, dico_bim_rs[dico_bim_chrpos[chrpos][len(dico_bim_chrpos[chrpos])-1]][len(dico_bim_rs[dico_bim_chrpos[chrpos][len(dico_bim_chrpos[chrpos])-1]])-1][1])
			add_liste_pb(liste_double_chrpos_bim,snp,nb1)
			if len(dico_bim_chrpos[chrpos])==1:# writting the first duplicate
				FSDUPLPOS.write(dico_bim_chrpos[chrpos][0]+"\n")
			FSDUPLPOS.write(rs+"\n")
		elif bim_duplicate_result == 4: # Duplicate Type 4
			double_rs_bim += 1
			snp = (rs, chrpos, dico_bim_rs[rs][len(dico_bim_rs[rs])-1][0], alleles, dico_bim_rs[rs][len(dico_bim_rs[rs])-1][1])
			add_liste_pb(liste_double_rs_bim,snp,nb1)
			FSDUPLERROR.write(rs+"\n") # do we need to write the others duplicates ? they have the same name ....
		elif bim_duplicate_result == -1 : # programmeur mistake : forget a case
			exit(0)

		# we count the number of SNP in each duplicate
		if rs in dico_bim_rs:
			nb_name_duplicate = len(dico_bim_rs[rs])
			
			if nb_name_duplicate > 1:
				dico_count_lvl_name_duplicate[nb_name_duplicate] -= 1
			
			if nb_name_duplicate+1 not in dico_count_lvl_name_duplicate:
				dico_count_lvl_name_duplicate[nb_name_duplicate+1] = 1
			else:
				dico_count_lvl_name_duplicate[nb_name_duplicate+1] += 1

		if chrpos in dico_bim_chrpos :
			nb_pos_duplicate = len(dico_bim_chrpos[chrpos])
			
			if nb_pos_duplicate > 1:
				dico_count_lvl_pos_duplicate[nb_pos_duplicate] -= 1
			
			if nb_pos_duplicate+1 not in dico_count_lvl_pos_duplicate:
				dico_count_lvl_pos_duplicate[nb_pos_duplicate+1] = 1
			else:
				dico_count_lvl_pos_duplicate[nb_pos_duplicate+1] += 1


	# we put the SNP in the two dictionnaries
	if rs not in dico_bim_rs:
		dico_bim_rs[rs] = []
	if chrpos not in dico_bim_chrpos:
		dico_bim_chrpos[chrpos] = []
	
	curent_snp = (chrpos, alleles)

	dico_bim_rs[rs].append(curent_snp)
	dico_bim_chrpos[chrpos].append(rs)


#we write the results on the file

FSOR.write("Number of SNPs with unknown chromosomes: "+str(compteur_chro_0)+" \n")

FSOR.write("Number of SNPs on X chromosomes (=23): "+str(compteur_chro_23)+" \n")

FSOR.write("Number of SNPs on Y chromosomes (=24): "+str(compteur_chro_24)+" \n")

FSOR.write("Number of SNPs on pseudo-autosomal chromosomes (=25): "+str(compteur_chro_25)+" \n")

FSOR.write("Number of SNPs on mitochondrial chromosomes (=26): "+str(compteur_chro_26)+" \n")

FSOR.write("Number of SNPs on anormal chromosomes (>=27): "+str(compteur_chro_sup_27)+" \n")


FSOR.write(str("\n---------------- Check of SNP positions ---------------------- \n"))
FSOR.write("Number of SNPs with known positions: "+str(compteur_goodpos)+"\n")

FSOR.write("Number of SNPs with unknown positions: "+str(compteur_badpos)+"\n")


FSOR.write(str("\n---------------- SNPs without chromosome and position ------------------ \n"))
FSOR.write("Number of SNPs without positions and chromosomes: "+str(compteur_bad_chrpos)+"\n")
FSOR.write("Number of SNPs with missing information on position or chromosome: "+str(compteur_missing_pos)+"\n")

FSOR.write(str("\n---------------- Check of SNP names ---------------------- \n"))

FSOR.write("Number of SNPs nomenclature in rs[number]: "+str(compteur_goodrs)+"\n")

FSOR.write("Number of mitochondrials SNPs (name starts by 'Mito'): "+str(compteur_mito)+"\n")

FSOR.write("Number of mitochondrials SNPs (name starts by 'Mito') with unknown position: "+str(compteur_mito_pos)+"\n")

FSOR.write("Number of mitochondrials SNPs (name starts by 'Mito') placed on a wrong chromosome: "+str(compteur_mito_wrong)+"\n")

FSOR.write("Number of SNPs with a nomenclature different from rs[number]: "+str(compteur_badrs)+"\n")

FSOR.write("\nList of the first mitochondrial SNPs (name starting with 'Mito') on a wrong chromosome (not on 26 or 0):")
if liste_wrong_mito_bim != [] :
	for index in range(len(liste_wrong_mito_bim)):
		FSOR.write("\n" + liste_wrong_mito_bim[index][0] + " on chr : " + str(liste_wrong_mito_bim[index][1]))
else:
	FSOR.write("\nNone")

FSOR.write("\n\nList of the first SNP with a name different from rs[number] or 'Mito':")
if liste_notrs != [] :
	for index in range(len(liste_notrs)):
		FSOR.write("\n" + liste_notrs[index])
	FSOR.write("\n")
else:
	FSOR.write("\nNone\n")

FSOR.write(str("\n---------------- Check of allele codes ---------------------- \n"))

FSOR.write("\nNumber of SNPs with both allele codes as A, T, C, G: "+str(couple_good_allele) +"\n")
FSOR.write("Number of SNPs with at least one missing allele code (if only one is missing, the other is A, T, C or G): " + str(couple_missing_allele)+"\n")
FSOR.write("Number of SNPs with at least one allele code different from A, T, C, G, 0: " + str(couple_unusual_allele) +"\n")

FSOR.write("\nNumber of allele code A, for the first allele: "+str(compteur_allele1["A"]) + ", for the second allele: " + str(compteur_allele2["A"]) + "\n")
FSOR.write("Number of allele code T, for the first allele: " + str(compteur_allele1["T"]) + ", for the second allele: " + str(compteur_allele2["T"]) + "\n")
FSOR.write("Number of allele code G, for the first allele: " + str(compteur_allele1["G"]) + ", for the second allele: " + str(compteur_allele2["G"]) + "\n")
FSOR.write("Number of allele code C, for the first allele: " + str(compteur_allele1["C"]) + ", for the second allele: " + str(compteur_allele2["C"]) + "\n")
FSOR.write("Number of SNPs with 0 as allele code (missing), for the first allele: " + str(compteur_allele1["0"]) + ", for the second allele: " + str(compteur_allele2["0"]) + "\n\n")

#compteur_bad_allele1=0
for key, value in compteur_allele1.items():
	if key != 'A' and key != 'T' and key != 'G' and key != 'C' and key != '0':
		FSOR.write("Number of allele code "+str(key)+" for the first allele: "+ str(value) +"\n")
		#compteur_bad_allele1 += value

#compteur_bad_allele2=0
for key, value in compteur_allele2.items():
	if key!="A" and key!="T" and key!="G" and key!="C" and key!="0":
		FSOR.write("Number of allele code "+str(key)+" for the second allele: "+ str(value) +"\n")
		#compteur_bad_allele2 += value


FSOR.write("\nList of the first SNPs with an allele code different from A, T, C, G or 0:")
if liste_diff_alleles != []:
	for index in range(len(liste_diff_alleles)):
		FSOR.write("\n" + liste_diff_alleles[index][0] + "\t" + liste_diff_alleles[index][1] + "\t" + liste_diff_alleles[index][2])

else:
	 FSOR.write("\nNone")

FSOR.write("\n\n**************Identification of duplicates*********\nThese SNP numbers should be taken with caution when information on SNP are not completed (missing information on chromosomal position or allele codes).\n")

FSOR.write("\n\n***TYPE 1:\nNumber of SNPs that appear at least twice in the Bim file with same chr:pos and same allele codes (whatever their SNP name): " + str(double_same_allele_bim) +"\n")

FSOR.write("\n***TYPE 2:\nNumber of SNPs (same name, chr and pos) that appear at least twice with different allele codes: " + str(double_rschrpos_bim) +"\n")
FSOR.write("\nList of the first SNPs Type 2:\n")
if liste_double_rschrpos_bim != [] :
	for index in range(len(liste_double_rschrpos_bim)):
		for i in range(len(liste_double_rschrpos_bim[index])):
			FSOR.write(liste_double_rschrpos_bim[index][i] +"\t")
		FSOR.write("\n")
else :
	FSOR.write("None\n")

FSOR.write("\n***TYPE 3:\nNumber of SNPs (same chr and pos but different names) that appear at least twice and with different allele codes: " + str(double_chrpos_bim) +"\n")
FSOR.write("\nList of the first SNPs Type 3:\n")
if liste_double_chrpos_bim != [] :
	for index in range(len(liste_double_chrpos_bim)):
		for i in range(len(liste_double_chrpos_bim[index])):
			FSOR.write(liste_double_chrpos_bim[index][i] +"\t")
		FSOR.write("\n")
else :
	FSOR.write("None\n")

FSOR.write("\n***TYPE 4:\nNumber of SNPs that appear at least twice with the same name but with different chr pos (whatever their allele codes are): " + str(double_rs_bim) +"\n")
FSOR.write("\nList of the first SNPs Type 4:\n")
if liste_double_rs_bim != [] :
	for index in range(len(liste_double_rs_bim)):
		for i in range(len(liste_double_rs_bim[index])):
			FSOR.write(liste_double_rs_bim[index][i] +"\t")
		FSOR.write("\n")
else:
	FSOR.write("None\n")

FSOR.write("\nNumber of times a SNP appears in the bim file with the same name:\n")
flag_none = 0
for key,value in dico_count_lvl_name_duplicate.items():
	if value != 0:
		FSOR.write(str(value)+" SNP duplicates that have "+str(key)+" occurences in the bim file\n")
		flag_none = 1
if flag_none == 0:
	FSOR.write("None\n")

flag_none = 0
FSOR.write("\nNumber of times a SNP appears in the bim file with the same position:\n")
for key,value in dico_count_lvl_pos_duplicate.items():
	if value != 0:
		FSOR.write(str(value)+" SNP duplicates that have "+str(key)+" occurences in the bim file\n")
		flag_none = 1
if flag_none == 0:
	FSOR.write("None\n")


f_bim.close() #we close the bim file


FSOR.write(str("\n\n---------------- Check manifest file ---------------------- \n"))

### function which verify if SNP with a name starting by Mito is not on Mitochondrial chromosome or unknow chromosome
def pos_mito(rs, chropos, liste, nb):
	chro = chropos.split(":")[0]
	if rs.startswith("Mito") and (chro != "26" and chro != "MT" and chro != "M" and chro != "chrM" and chro !="0"):
		mito_error = (rs, chro)
		add_liste_pb(liste, mito_error, nb1)
		nb +=1
	return nb

############### Declaring dictionnary for each strand notation ##############################

# we want to identify the different chromosome of the Manifest
manifest_chr = {}

# for the +/- notation
dicostrand_p = {} #use a dict to save the informations of the SNP and their alleles in the manifest file
dico_chropos_strand_p = {} #use an another dict if the SNP names is chr:pos in the bim file

# for the For/Rev notation
dicostrand_f = {} #use a dict to save the informations of the SNP and their alleles in the manifest file
dico_chropos_strand_f = {} #use an another dict if the SNP names is chr:pos in the bim file

# for the Top/Bot notation
dicostrand_t = {} #use a dict to save the informations of the SNP and their alleles in the manifest file
dico_chropos_strand_t = {} #use an another dict if the SNP names is chr:pos in the bim file


line_manifest_control=len(open(manifest_file).readlines())-nltail #remove the tail of the manifest file
m_line = 0 # use to know the number of SNP in the manifest file


###### function which count the SNP in double in the manifest file and write SNP information in the two dictionnary where we are going to search SNP from the bim file
def write_dico_strand(rs, rs_chropos, alleles, strand, dicostrand, dico_chropos_strand, double_rs_man, double_chrpos_man, double_rschrpos_man, double_same_allele, build, liste_double_rs_man, liste_double_chrpos_man, liste_double_rschrpos_man):
	# search for SNP which appear twice in the manifest
	allele1 = str(alleles[0])
	allele2 = str(alleles[1])
	
	permute_alleles = allele2+allele1
	if (rs in dicostrand):
		if (dicostrand[rs][0]==rs_chropos):# and (chro != "0" and pos != "0"): # name and chr:pos identical
			if (dicostrand[rs][1] != alleles and dicostrand[rs][1] != permute_alleles)or dicostrand[rs][0] != rs_chropos:
				double_rschrpos_man +=1
				snp = (rs, rs_chropos, alleles, dicostrand[rs][1])
				add_liste_pb(liste_double_rschrpos_man,snp,nb1)
			else:
				double_same_allele +=1
		else:# same name but different position
			double_rs_man += 1
			snp = (rs, rs_chropos, dicostrand[rs][0], alleles, dicostrand[rs][1])
			add_liste_pb(liste_double_rs_man,snp,nb1)

	elif rs_chropos in dico_chropos_strand and chro != "0" and pos != "0": # if only chr:pos are identical
		if dicostrand[dico_chropos_strand[rs_chropos]][1] != alleles and dicostrand[dico_chropos_strand[rs_chropos]][1] != permute_alleles:
			double_chrpos_man +=1
			snp = (rs, dico_chropos_strand[rs_chropos], rs_chropos, alleles, dicostrand[dico_chropos_strand[rs_chropos]][1])
			add_liste_pb(liste_double_chrpos_man,snp,nb1)
		else:
			double_same_allele += 1

	# else : it's not duplicate
	# keep information of the SNP in the manifest to compare it with the bim file
	dicostrand[rs] = (rs_chropos, alleles, strand, build) # rs_chropos have the notation chr:pos, alleles are on the +/For/Top strand, strand is the original strand of the SNP in the manifest
	dico_chropos_strand[rs_chropos] = rs

	return double_rs_man, double_chrpos_man, double_rschrpos_man, double_same_allele 


missing_chro = 0
missing_pos = 0
missing_chropos = 0
chro_wrong_value = 0
pos_wrong_value = 0

dico_char_PM = {}
dico_char_FR = {}
dico_char_TB = {}

########################################################## we read the manifest file ##########################################################
with open(manifest_file) as file:
	for nb_line,line in enumerate(file):
		if nb_line >= nlhead and nb_line < line_manifest_control:#without the head and the tail (without control ...) if they were specified
			if re.search(r'(\t)+',line) != None: # we verify what kind of separator we have
				fields=line.split('\t')
			elif re.search(r'(,)+',line) != None: # separator are commas
				fields=line.split(",")
			else:
				sys.stderr.write("\n\n**********************************************************\n\nSeparator is neither tabulation nor comma !\nVerify your file, if necessary add it in the code (pre-processing.py line 699)\n\n**********************************************************\n\n")
				exit(0)
			# getting all informations needed
			rs = fields[field_name].strip() #name of the SNP
			alleles = fields[field_allele].strip() # field with allele 
			res = re.search(r'^\[(\D+)/(\D+)\]',alleles) # we test if alleles have the notation [A/B]
			if res != None:
				allele1 = res.group(1)
				allele2 = res.group(2)
				alleles = allele1 + allele2
			elif len(alleles) == 2: # else they should have the notation AB
				allele1 = alleles[0]
				allele2 = alleles[1]
			elif len(alleles) == 1: # else we should have given the manifest normalise, the second allele is in the next field
				allele1 = alleles
				allele2 = fields[(field_allele+1)].strip()
				alleles = str(allele1) + str(allele2)
			else:
				sys.stderr.write("\n\n**********************************************************\n\nAlleles "+ alleles +" different from [A/B] and AB notation !\n(code in pre-processing.py line 717)\n\n**********************************************************\n\n")
				exit(0)
			
			#we count the different alleles found and SNP kind
			compteur_allele1_manifest, compteur_allele2_manifest, couple_good_allele_manifest, couple_unusual_allele_manifest, couple_missing_allele_manifest = check_allele(rs, allele1, allele2, compteur_allele1_manifest, compteur_allele2_manifest, liste_diff_alleles_manifest, couple_good_allele_manifest, couple_unusual_allele_manifest, couple_missing_allele_manifest)
	
			
			chro = fields[field_chr].strip()
			pos = fields[field_pos].strip()
			rs_chropos=chro+":"+pos # we put chr and pos information in the notation chr:pos
			
			if chro == "0": # verify if there is missing position in the Manifest file
				missing_chro += 1
			if pos == "0" :
				missing_pos += 1
			if chro == "0" and pos == "0":
				missing_chropos += 1
			
			if chro.isdigit():
				if int(chro) > 22:
					if chro in manifest_chr:
						manifest_chr[chro] += 1
					else:
						manifest_chr[chro] = 1
			else:
				if chro in manifest_chr:
					manifest_chr[chro] += 1
				else:
					manifest_chr[chro] = 1

			if chro != "X" and chro != "Y"  and chro != "XY" and chro != "M" and chro != "MT" and chro != "chrM": # if the chromosome value is not one of these letters, it has to be a number between 0 and 26
				valchro = re.search(r'^\d+$',chro) # verify if the chromosome value is a number
				if valchro != None:
					intchro = int(chro)
					if intchro < 0 or intchro > 26 : # verify that the number correspond to chromome value
						chro_wrong_value += 1
						sys.stderr.write("\n\n**********************************************************\n\nIn the Manifest : Strange Chromosome value : '"+chro+"'\n\n**********************************************************\n\n")
				else :
					chro_wrong_value += 1
					sys.stderr.write("\n\n**********************************************************\n\nIn the Manifest : Strange Chromosome value : '"+chro+"'\n\n**********************************************************\n\n")

			valpos = re.search(r'^\d+$',pos) # verify if the position value is a number
			if valpos == None:
				pos_wrong_value += 1
				sys.stderr.write("\n\n**********************************************************\n\nStrange Chr:Position value : '"+chro+":"+pos+"'\n\n**********************************************************\n\n")
			

			# we verify if there is a SNP with the name Mito on a chr different from 0, 26, MT, M
			compteur_mito_wrong_manifest = pos_mito(rs, rs_chropos, liste_wrong_mito_manifest, compteur_mito_wrong_manifest)
	
			if field_strandP != None: # if the strand notation +/- is specified, we verify if we have only + - characters (else we can not treat them latter, so we stop the execution)
				PMstrand=fields[field_strandP].strip() # strand of the SNP
				if PMstrand != "+" and PMstrand != "-":					
					if PMstrand not in dico_char_PM :
						dico_char_PM[PMstrand] = 1
					else :
						dico_char_PM[PMstrand] += 1

			if field_strandF != None: # if For / Rev notation specified
				FRstrand = fields[field_strandF].strip()
				if FRstrand != 'F' and FRstrand != 'R':
					IlmnID = fields[field_strandF].split("_") # read only the notation of the field IlmnID			
					FRstrand = IlmnID[len(IlmnID)-2].strip() # it's in the string, before the last field (tab delimited by '_')
					
					if FRstrand != "F" and FRstrand != "R":
						if FRstrand not in dico_char_FR :
							dico_char_FR[FRstrand] = 1
						else :
							dico_char_FR[FRstrand] += 1

			
			if field_strandT != None: # if TOP / Bot notation specified
				TBstrand = fields[field_strandT].strip()
				TBstrand = TBstrand.upper() # in case characters are in minuscule top or Top instead of TOP. by know we made the change, I treat only upper case notation (TOP) 
				if TBstrand != "TOP" and TBstrand != "BOT" and TBstrand != "PLUS"  and TBstrand != "MINUS": # when alleles ID, notation PLUS and MINUS
					if TBstrand not in dico_char_TB :
						dico_char_TB[TBstrand] = 1
					else :
						dico_char_TB[TBstrand] += 1
	
			if field_build != None :
				build = fields[field_build].strip()
			else :
				build = "0"

			#If we don't want to write exceptions like NA, put the next lines (653 to 695) in this condition : (dans ces lignes il y a, l'ecriture dans le fichier manifest normalise, dans le dico, le contage des doublons) On peut vouloir se contente de ne pas ecrire dans le fichier manifest normaliser. Dans ce cas decaler seulement les lignes de 654 a 661.
			#if allele1 != 'N' or allele2 != 'N':
				
			# Write in a file with the chosen information of the manifest to help treatment for the next step the Strand is specify in the order : +/- For/Rev Top/Bot
			FSMAN.write(chro +"\t"+ pos +"\t"+ rs +"\t"+ allele1 +"\t"+ allele2)
			if field_strandP != None:
				FSMAN.write("\t"+ PMstrand)
			if field_strandF != None:
				FSMAN.write("\t"+ FRstrand)
			if field_strandT != None:
				FSMAN.write("\t"+ TBstrand)
			FSMAN.write("\t"+ build)
			FSMAN.write("\n")
			
			# For each nomenclature specified, the dicionary of the nomenclature is filled with allele in +/For/Top, but we need to know which was the original strand in the manifest !!!!
			if field_strandP != None:
				alleles_p = alleles
				strand_kind_p = 0 # it's 0 if you don't need to flip the SNP in this nomenclature
				if PMstrand == "-": #we need the change to +
					alleles_p = complementary_strand_manifest(allele1, allele2) # if the strand is -, we get the complementary allele
					man_minus_strand_p += 1
					strand_kind_p = 1 # it's 1 if you need to flip the SNP in this nomenclature
					#PMstrand = "+"
				if PMstrand != "+" and PMstrand != "-":
					strand_kind_p = 2
				double_rs_man_p, double_chrpos_man_p, double_rschrpos_man_p, double_same_allele_manifest_p = write_dico_strand(rs, rs_chropos, alleles_p, strand_kind_p, dicostrand_p, dico_chropos_strand_p, double_rs_man_p, double_chrpos_man_p, double_rschrpos_man_p, double_same_allele_manifest_p, build, liste_double_rs_man_p, liste_double_chrpos_man_p, liste_double_rschrpos_man_p)

			if field_strandF != None: # For/Rev : we do the same step than for the +/- nomenclature
				alleles_f = alleles
				strand_kind_f = 0
				if FRstrand == "R":
					alleles_f = complementary_strand_manifest(allele1, allele2)
					man_minus_strand_f += 1
					strand_kind_f = 1
					#FRstrand = "F"
				if FRstrand != "F" and FRstrand != "R":
					strand_kind_f = 2
				double_rs_man_f, double_chrpos_man_f, double_rschrpos_man_f, double_same_allele_manifest_f = write_dico_strand(rs, rs_chropos, alleles_f, strand_kind_f, dicostrand_f, dico_chropos_strand_f, double_rs_man_f, double_chrpos_man_f, double_rschrpos_man_f, double_same_allele_manifest_f, build, liste_double_rs_man_f, liste_double_chrpos_man_f, liste_double_rschrpos_man_f)

			if field_strandT != None:# Top/Bot : we do the same step than for the +/- nomenclature
				alleles_t = alleles
				strand_kind_t = 0
				if TBstrand == "BOT" or TBstrand == "MINUS":
					alleles_t = complementary_strand_manifest(allele1, allele2)
					man_minus_strand_t += 1
					strand_kind_t = 1
					#if TBstrand == "MINUS":
					#       TBstrand = "PLUS"
					#else:
					#       TBstrand = "TOP"
				if TBstrand != "TOP" and TBstrand != "BOT" and TBstrand != "PLUS" and TBstrand != "MINUS":
					strand_kind_t = 2
				double_rs_man_t, double_chrpos_man_t, double_rschrpos_man_t, double_same_allele_manifest_t = write_dico_strand(rs, rs_chropos, alleles_t, strand_kind_t, dicostrand_t , dico_chropos_strand_t, double_rs_man_t, double_chrpos_man_t, double_rschrpos_man_t, double_same_allele_manifest_t, build, liste_double_rs_man_t, liste_double_chrpos_man_t, liste_double_rschrpos_man_t)

		
			m_line = (nb_line + 1)-nlhead
				
flag_sep = 0
#if the rsID is given, we create a dictionary with the SNPs names
snp_link = {} #dictionnary that will be like: {'kgp...': 'rs...'}
if file_rsID != None:
	with open(file_rsID) as loci_file: #file_rsID is the loci to rsID file given in the terminal
		for line_nb, line in enumerate(loci_file):
			#if line_nb >= 1: #because the 1st line is the name of the columns
			
			if flag_sep == 0:
				if re.search(r'(\t)+',line): #try if the separators are tabulations
					flag_sep = 1
				elif re.search(r'(,)+',line): #try if the separators are commas
					flag_sep = 2
				elif re.search(r'( )+',line): #try if the separators are spaces
					flag_sep = 3
				else :
					sys.stderr.write("\n\n**********************************************************\n\nError : Separators of the link file are neither tabulations nor commas or spaces !\n\n**********************************************************\n\n")
					exit(0)
			if flag_sep == 1:
				data = line.split('\t')
			elif flag_sep == 2:
				data = line.split(',')
			else :
				data = line.split(' ')
			
			new_snp = data[0] #the SNP that is not necessarily in the 'rs[number]' nomenclature
			good_snp = data[1].strip() #the SNP that is in the 'rs[number]' nomenclature
			
			if good_snp != '.' and good_snp != '' and good_snp != ' ': #there is a name for the SNP
				if ',' in good_snp: #there was a merge between several SNPs
					snp_name = []
					snp_name = good_snp.split(',') #selection of the good one
					for var in snp_name:
						snp_link[var] = new_snp
				else:
					snp_link[good_snp] = new_snp ##construction of the dictionnary, key = str, value = str

######### Declaring all variable to count the strand in the bim file ##################
nb_line=len(open(bim).readlines())#calculate the line number in the bim file
number_find_p=0
number_find_f=0
number_find_t=0

number_pos=0
number_neg=0
number_for=0
number_rev=0
number_top=0
number_bot=0

number_imp_p=0
number_imp_f=0
number_imp_t=0

number_imp_pos=0
number_imp_neg=0
number_imp_for=0
number_imp_rev=0
number_imp_top=0
number_imp_bot=0

number_missing_p=0
number_missing_f=0
number_missing_t=0

number_ambiguous_p=0
number_ambiguous_f=0
number_ambiguous_t=0

number_notation_p=0
number_notation_f=0
number_notation_t=0

nb_diff_strand_BM_p=0
nb_diff_strand_BM_f=0
nb_diff_strand_BM_t=0

man_minus_strand_found_p=0
man_minus_strand_found_f=0
man_minus_strand_found_t=0

man_unknown_strand_p=0
man_unknown_strand_f=0
man_unknown_strand_t=0

number_diff_pos_BM=0

dico_build = {}

########### function listing SNP by build in a file #################
def compare_SNPbuild(rs, build, dico_build):
	if build != "0" and print_summary == '1' : # we check the build only if it is the first QC
		if build in dico_build :
			dico_build[build].write(rs+"\n")
		else :
			dico_build[build] = open(cancer+"_SNP_build_"+build+".txt", 'w')
			dico_build[build].write(rs+"\n")



########### function that compare the SNP position in the bim file and the manifest file ###########
def compare_chrpos_BM(rs, chr_pos, man_chrpos, number_diff_pos_BM, liste_diff_pos_BM, fichier):
	res = re.search(r'([A-Z]+):(\d+)',man_chrpos) # verify if the position in the manifest contein a letter instead of a number for the chromosome
	if res != None: 
		man_chr = res.group(1)
		man_pos = res.group(2)
		if man_chr == "M" or man_chr == "MT" or man_chr == "chrM": # convertion of the letter in the number notation
			man_chr = "26"
		elif man_chr == "XY":
			man_chr = "25"
		elif man_chr == "Y":
			man_chr = "24"
		elif man_chr == "X":
			man_chr = "23"
		else:
			sys.stderr.write("\n\n**********************************************************\n\nThere is a strange notation of chr (line 945): '"+man_chr+"' !\n\n**********************************************************\n\n")
		man_chrpos = man_chr +":"+man_pos
	bimresult = re.search(r'(\d+):(\d+)',chr_pos)
	manresult = man_chrpos.split(":")
	if bimresult.group(1) == "0": # if chr is missing but we have the position, we compare the position
		if bimresult.group(2)!= "0" and manresult[1] != "0" and bimresult.group(2) != manresult[1]:
			number_diff_pos_BM += 1 # we count it
			snp = (rs, chr_pos, man_chrpos)
			add_liste_pb(liste_diff_pos_BM, snp, nb1) # put the first one in the list
			fichier.write(rs+"\n") # write them in a file
	elif bimresult.group(2) == "0": # if we have the chr but the position is missing, we compare the chr
		if bimresult.group(1) != manresult[0] and manresult[0] != 0 :
			number_diff_pos_BM += 1 # we count it
			snp = (rs, chr_pos, man_chrpos)
			add_liste_pb(liste_diff_pos_BM, snp, nb1) # put the first one in the list
			fichier.write(rs+"\n") # write them in a file
	elif manresult[0] == "0" : # in the bim file there are the two positions information, we test if it's the case for the Manifest file
		if manresult[1] != "0" and bimresult.group(2) != manresult[1]:
			number_diff_pos_BM += 1
			snp = (rs, chr_pos, man_chrpos)
			add_liste_pb(liste_diff_pos_BM, snp, nb1)
			fichier.write(rs+"\n")
	elif manresult[1] == "0": # in the Manifest we have the chr but the position is missing, we compare the chr
		if bimresult.group(1) != manresult[0]:
			number_diff_pos_BM += 1
			snp = (rs, chr_pos, man_chrpos)
			add_liste_pb(liste_diff_pos_BM, snp, nb1)
			fichier.write(rs+"\n")
	else:
		if chr_pos != man_chrpos : # if position in the bim and the manifest is different
			number_diff_pos_BM += 1 # we count it
			snp = (rs, chr_pos, man_chrpos)
			add_liste_pb(liste_diff_pos_BM, snp, nb1) # put the first one in the list
			fichier.write(rs+"\n") # write them in a file
	return number_diff_pos_BM


################################# read the bim file to compare it with the manifest file ########################################## 
with open(bim) as file:
	for line in file:
		sep = re.search(r'(\t)+',line) #try if the separators are tabulations
		if sep: # separator are tabulations, sep != None
			fields=line.split("\t")
		else: #sep == None : separators are spaces
			fields=line.split(" ")
		# get information of the SNP of the bim file
		rs=fields[1].strip()
		allele1=fields[4].strip()
		allele2=fields[5].strip()
		alleles=allele1+allele2
		chro=fields[0].strip()
		pos=fields[3].strip()
		chr_pos=chro+":"+pos
		
		################## if we have the +/- strand notation specified #################
		if field_strandP != None:
			flag = 1
			if rs in dicostrand_p: #if the name of the SNP in the bim file is in the manifest file
				#sys.stderr.write("\ndicostrand_P "+rs)
				compare_SNPbuild(rs, dicostrand_p[rs][3], dico_build)
				# we compare the position
				number_diff_pos_BM = compare_chrpos_BM(rs, chr_pos, dicostrand_p[rs][0], number_diff_pos_BM, liste_diff_pos_BM, FSDIFPOS)
				# we verify his strand
				number_find_p, number_pos, number_neg, number_ambiguous_p, number_missing_p, number_imp_p, number_imp_pos, number_imp_neg, nb_diff_strand_BM_p , number_notation_p, man_minus_strand_found_p,man_unknown_strand_p = research_strand(rs, dicostrand_p[rs], allele1, allele2, number_find_p, number_pos, number_neg, number_ambiguous_p, number_missing_p, number_imp_p, number_imp_pos, number_imp_neg, nb_diff_strand_BM_p, number_notation_p, man_minus_strand_found_p, man_unknown_strand_p, liste_diff_bim_manifest_p, flag, FSMINUS, FSDIFALL, FSUNSTRANDp)
				# we put all SNP in the negative strand in the manifest in a file (in case we need to flip them) -> in a raw bim file, the strand of a SNP should be the same as the manifest
				if dicostrand_p[rs][2] == 1: # so we base our decision on the manifest strand and not on the bim strand. (we suppose that if the bim file is not raw, it is already on one strand of one of the notation)
					FSSTRANDp.write(rs+"\n")
				
				""" # this two research might be a problem with duplicate SNP, to complete positions or alleles we don't use them anyway... can bring flase count in the QC results
			elif rs in dico_chropos_strand_p: # if the name of the SNP in the bim file is a chr:pos, we look for it in manifest dictionnary with all positions
				#sys.stderr.write("\ndico_chropos_P "+rs)
				# we compare the position
				number_diff_pos_BM = compare_chrpos_BM(rs, chr_pos, dicostrand_p[dico_chropos_strand_p[rs]][0], number_diff_pos_BM, liste_diff_pos_BM, FSDIFPOS)
				# we verify his strand
				number_find_p, number_pos, number_neg, number_ambiguous_p, number_missing_p, number_imp_p, number_imp_pos, number_imp_neg, nb_diff_strand_BM_p, number_notation_p, man_minus_strand_found_p = research_strand(rs, dicostrand_p[dico_chropos_strand_p[rs]], allele1, allele2, number_find_p, number_pos, number_neg, number_ambiguous_p, number_missing_p, number_imp_p, number_imp_pos, number_imp_neg, nb_diff_strand_BM_p, number_notation_p, man_minus_strand_found_p, liste_diff_bim_manifest_p, flag, FSMINUS, FSDIFALL)
				# we put all SNP in the negative strand in the manifest in a file
				if dicostrand_p[dico_chropos_strand_p[rs]][2] == 1:
					FSSTRANDp.write(rs+"\n")
	
			elif chr_pos in dicostrand_p: # if the name of the SNP in the manifest file is a chr:pos, we take the chr:pos of the SNP from the bim file and search it in the manifest dictionnary with the name of the SNP
				#sys.stderr.write("\ndico_gall_P "+chr_pos)
				# we compare the position
				number_diff_pos_BM = compare_chrpos_BM(rs, chr_pos, dicostrand_p[chr_pos][0], number_diff_pos_BM, liste_diff_pos_BM, FSDIFPOS)
				# we verify his strand
				number_find_p, number_pos, number_neg, number_ambiguous_p, number_missing_p, number_imp_p, number_imp_pos, number_imp_neg, nb_diff_strand_BM_p, number_notation_p, man_minus_strand_found_p = research_strand(rs, dicostrand_p[chr_pos], allele1, allele2, number_find_p, number_pos, number_neg, number_ambiguous_p, number_missing_p, number_imp_p, number_imp_pos, number_imp_neg, nb_diff_strand_BM_p, number_notation_p, man_minus_strand_found_p, liste_diff_bim_manifest_p, flag, FSMINUS, FSDIFALL)
				# we put all SNP in the negative strand in the manifest in a file
				if dicostrand_p[chr_pos][2] == 1:
					FSSTRANDp.write(rs+"\n")
				"""
			elif file_rsID != None and rs in snp_link: # if the file with the link betwen new and old SNP name have been specified, and if we found the SNP name of the bim file in the rsID file (second column with rs, in case the bim file have already been annotate)
				#sys.stderr.write("\ndico_rsID_P "+rs)
				new_rs=snp_link[rs] #we take the old SNPs name (first column. It should correspond to the manifest name)
				if new_rs in dicostrand_p: # if the old SNPs name is found in the manifest file
					compare_SNPbuild(rs, dicostrand_p[new_rs][3], dico_build)
					# we compare the position
					number_diff_pos_BM = compare_chrpos_BM(rs, chr_pos, dicostrand_p[new_rs][0], number_diff_pos_BM, liste_diff_pos_BM, FSDIFPOS)
					# we verify his strand
					#print("\n"+str(dicostrand_p[new_rs])+" !!!!\n")
					number_find_p, number_pos, number_neg, number_ambiguous_p, number_missing_p, number_imp_p, number_imp_pos, number_imp_neg, nb_diff_strand_BM_p, number_notation_p, man_minus_strand_found_p, man_unknown_strand_p = research_strand(rs, dicostrand_p[new_rs], allele1, allele2, number_find_p, number_pos, number_neg, number_ambiguous_p, number_missing_p, number_imp_p, number_imp_pos, number_imp_neg, nb_diff_strand_BM_p, number_notation_p, man_minus_strand_found_p, man_unknown_strand_p, liste_diff_bim_manifest_p, flag, FSMINUS, FSDIFALL, FSUNSTRANDp)
					# we put all SNP in the negative strand in the manifest in a file
					if dicostrand_p[new_rs][2] == 1:
						FSSTRANDp.write(rs+"\n")
				else:
					FSNOTFOUND.write(rs+"\n")
			else : #the SNP is not find
				FSNOTFOUND.write(rs+"\n")
			

		################## if we have the Forward/Reverse strand notation specified (same behaviour as the +/- nomenclature) #################
		if field_strandF!= None:
			flag = 0
			if rs in dicostrand_f: # if the name of the SNP in the bim file is in the manifest file
				#sys.stderr.write("\ndicostrand_F "+rs)
				if field_strandP == None:
					compare_SNPbuild(rs, dicostrand_f[rs][3], dico_build)
					number_diff_pos_BM = compare_chrpos_BM(rs, chr_pos, dicostrand_f[rs][0], number_diff_pos_BM, liste_diff_pos_BM, FSDIFPOS)
					flag = 1
				number_find_f, number_for, number_rev, number_ambiguous_f, number_missing_f, number_imp_f, number_imp_for, number_imp_rev, nb_diff_strand_BM_f, number_notation_f, man_minus_strand_found_f, man_unknown_strand_f = research_strand(rs, dicostrand_f[rs], allele1, allele2, number_find_f, number_for, number_rev, number_ambiguous_f, number_missing_f, number_imp_f, number_imp_for, number_imp_rev, nb_diff_strand_BM_f, number_notation_f, man_minus_strand_found_f, man_unknown_strand_f, liste_diff_bim_manifest_f, flag, FSREV, FSDIFALL, FSUNSTRANDf)
				if dicostrand_f[rs][2] == 1:
					FSSTRANDf.write(rs+"\n")
				""" # this two research might be a problem with duplicate SNP, to complete positions or alleles we don't use them anyway... can bring flase count in the QC results
			elif rs in dico_chropos_strand_f: # if the name of the SNP in the bim file is a chr:pos, we look for it in manifest dictionnary with all positions
				#sys.stderr.write("\ndico_chropos_F "+rs)
				if field_strandP == None:
					number_diff_pos_BM = compare_chrpos_BM(rs, chr_pos, dicostrand_f[dico_chropos_strand_f[rs]][0], number_diff_pos_BM, liste_diff_pos_BM, FSDIFPOS)
					flag = 1
				number_find_f, number_for, number_rev, number_ambiguous_f, number_missing_f, number_imp_f, number_imp_for, number_imp_rev, nb_diff_strand_BM_f, number_notation_f, man_minus_strand_found_f = research_strand(rs, dicostrand_f[dico_chropos_strand_f[rs]], allele1, allele2, number_find_f, number_for, number_rev, number_ambiguous_f, number_missing_f, number_imp_f, number_imp_for, number_imp_rev, nb_diff_strand_BM_f, number_notation_f, man_minus_strand_found_f, liste_diff_bim_manifest_f, flag, FSREV, FSDIFALL)
				if dicostrand_f[dico_chropos_strand_f[rs]][2] == 1:
					FSSTRANDf.write(rs+"\n")
	
			elif chr_pos in dicostrand_f: # if the name of the SNP in the manifest file is a chr:pos, we take the chr:pos of the SNP from the bim file and search it in the manifest dictionnary with the name of the SNP
				#sys.stderr.write("\ndico_gall_F "+chr_pos)
				if field_strandP == None:
					number_diff_pos_BM = compare_chrpos_BM(rs, chr_pos, dicostrand_f[chr_pos][0], number_diff_pos_BM, liste_diff_pos_BM, FSDIFPOS)
					flag = 1
				number_find_f, number_for, number_rev, number_ambiguous_f, number_missing_f, number_imp_f, number_imp_for, number_imp_rev, nb_diff_strand_BM_f, number_notation_f, man_minus_strand_found_f = research_strand(rs, dicostrand_f[chr_pos], allele1, allele2, number_find_f, number_for, number_rev, number_ambiguous_f, number_missing_f, number_imp_f, number_imp_for, number_imp_rev, nb_diff_strand_BM_f, number_notation_f, man_minus_strand_found_f, liste_diff_bim_manifest_f, flag, FSREV, FSDIFALL)
				if dicostrand_f[chr_pos][2] == 1:
					FSSTRANDf.write(rs+"\n")
				"""
			elif file_rsID != None and rs in snp_link: # if we found the SNP name of the bim file in the rsID file (second column with rs, in case the bim file have already been annotate)
				#sys.stderr.write("\ndico_rsID_F "+rs)
				new_rs=snp_link[rs] #we take the old SNPs name (manifest name)
				if new_rs in dicostrand_f: #and if the old SNPs name is found in the manifest file
					if field_strandP == None:
						compare_SNPbuild(rs, dicostrand_f[new_rs][3], dico_build)
						number_diff_pos_BM = compare_chrpos_BM(rs, chr_pos, dicostrand_f[new_rs][0], number_diff_pos_BM, liste_diff_pos_BM, FSDIFPOS)
						flag = 1
					number_find_f, number_for, number_rev, number_ambiguous_f, number_missing_f, number_imp_f, number_imp_for, number_imp_rev, nb_diff_strand_BM_f, number_notation_f, man_minus_strand_found_f, man_unknown_strand_f = research_strand(rs, dicostrand_f[new_rs], allele1, allele2, number_find_f, number_for, number_rev, number_ambiguous_f, number_missing_f, number_imp_f, number_imp_for, number_imp_rev, nb_diff_strand_BM_f, number_notation_f, man_minus_strand_found_f, man_unknown_strand_f, liste_diff_bim_manifest_f, flag, FSREV, FSDIFALL, FSUNSTRANDf)
					if dicostrand_f[new_rs][2] == 1:
						FSSTRANDf.write(rs+"\n")
				else:
					if field_strandP == None:
						FSNOTFOUND.write(rs+"\n")
			else : #the SNP is not find
				if field_strandP == None:
					FSNOTFOUND.write(rs+"\n")
	
		################## if we have the TOP/BOT strand notation specified (same behaviour as the +/- nomenclature) #################
		if field_strandT!= None:
			flag = 0
			if rs in dicostrand_t: #if the name of the SNP in the bim file is in the manifest file
				#sys.stderr.write("\ndicostrand_T "+rs)
				if field_strandP == None and field_strandF == None:
					compare_SNPbuild(rs, dicostrand_t[rs][3], dico_build)
					number_diff_pos_BM = compare_chrpos_BM(rs, chr_pos, dicostrand_t[rs][0], number_diff_pos_BM, liste_diff_pos_BM, FSDIFPOS)
					flag = 1
				number_find_t, number_top, number_bot, number_ambiguous_t, number_missing_t, number_imp_t, number_imp_top, number_imp_bot, nb_diff_strand_BM_t, number_notation_t, man_minus_strand_found_t, man_unknown_strand_t = research_strand(rs, dicostrand_t[rs], allele1, allele2, number_find_t, number_top, number_bot, number_ambiguous_t, number_missing_t, number_imp_t, number_imp_top, number_imp_bot, nb_diff_strand_BM_t, number_notation_t, man_minus_strand_found_t, man_unknown_strand_t, liste_diff_bim_manifest_t, flag, FSBOT, FSDIFALL, FSUNSTRANDt)
				if dicostrand_t[rs][2] == 1:
					FSSTRANDt.write(rs+"\n")
				""" # this two research might be a problem with duplicate SNP, to complete positions or alleles we don't use them anyway... can bring flase count in the QC results
			elif rs in dico_chropos_strand_t: # if the name of the SNP in the bim file is a chr:pos, we look for it in manifest dictionnary with all positions
				#sys.stderr.write("\ndico_chropos_T "+rs)
				if field_strandP == None and field_strandF == None:
					number_diff_pos_BM = compare_chrpos_BM(rs, chr_pos, dicostrand_t[dico_chropos_strand_t[rs]][0], number_diff_pos_BM, liste_diff_pos_BM, FSDIFPOS)
					flag = 1
				number_find_t, number_top, number_bot, number_ambiguous_t, number_missing_t, number_imp_t, number_imp_top, number_imp_bot, nb_diff_strand_BM_t, number_notation_t, man_minus_strand_found_t = research_strand(rs, dicostrand_t[dico_chropos_strand_t[rs]], allele1, allele2, number_find_t, number_top, number_bot, number_ambiguous_t, number_missing_t, number_imp_t, number_imp_top, number_imp_bot, nb_diff_strand_BM_t, number_notation_t, man_minus_strand_found_t, liste_diff_bim_manifest_t, flag, FSBOT, FSDIFALL)
				if dicostrand_t[dico_chropos_strand_t[rs]][2] == 1:
					FSSTRANDt.write(rs+"\n")
	
			elif chr_pos in dicostrand_t: # if the name of the SNP in the manifest file is a chr:pos, we take the chr:pos of the SNP from the bim file and search it in the manifest dictionnary with the name of the SNP
				#sys.stderr.write("\ndico_gall_T "+chr_pos)
				if field_strandP == None and field_strandF == None:
					number_diff_pos_BM = compare_chrpos_BM(rs, chr_pos, dicostrand_t[chr_pos][0], number_diff_pos_BM, liste_diff_pos_BM, FSDIFPOS)
					flag = 1
				number_find_t, number_top, number_bot, number_ambiguous_t, number_missing_t, number_imp_t, number_imp_top, number_imp_bot, nb_diff_strand_BM_t, number_notation_t, man_minus_strand_found_t = research_strand(rs, dicostrand_t[chr_pos], allele1, allele2, number_find_t, number_top, number_bot, number_ambiguous_t, number_missing_t, number_imp_t, number_imp_top, number_imp_bot, nb_diff_strand_BM_t, number_notation_t, man_minus_strand_found_t, liste_diff_bim_manifest_t, flag, FSBOT, FSDIFALL)	
				if dicostrand_t[chr_pos][2] == 1:
					FSSTRANDt.write(rs+"\n")
				"""
			elif file_rsID != None and rs in snp_link: # if we found the SNP name of the bim file in the rsID file 
				#sys.stderr.write("\ndico_rsID_T "+rs)
				new_rs=snp_link[rs] #we take the old SNPs name (manifest name)
				if new_rs in dicostrand_t: #and if the old SNPs name is found in the manifest file
					if field_strandP == None and field_strandF == None:
						compare_SNPbuild(rs, dicostrand_t[new_rs][3], dico_build)
						number_diff_pos_BM = compare_chrpos_BM(rs, chr_pos, dicostrand_t[new_rs][0], number_diff_pos_BM, liste_diff_pos_BM, FSDIFPOS)
						flag = 1
					number_find_t, number_top, number_bot, number_ambiguous_t, number_missing_t, number_imp_t, number_imp_top, number_imp_bot, nb_diff_strand_BM_t, number_notation_t, man_minus_strand_found_t, man_unknown_strand_t = research_strand(rs, dicostrand_t[new_rs], allele1, allele2, number_find_t, number_top, number_bot, number_ambiguous_t, number_missing_t, number_imp_t, number_imp_top, number_imp_bot, nb_diff_strand_BM_t, number_notation_t, man_minus_strand_found_t, man_unknown_strand_t, liste_diff_bim_manifest_t, flag, FSBOT, FSDIFALL, FSUNSTRANDt)
					if dicostrand_t[new_rs][2] == 1:
						FSSTRANDt.write(rs+"\n")
				else :
					if field_strandP == None and field_strandF == None:
						FSNOTFOUND.write(rs+"\n")
			else : #the SNP is not find
				if field_strandP == None and field_strandF == None:
					FSNOTFOUND.write(rs+"\n")



#for key, value in dicostrand_t.items():
#	sys.stderr.write(key + " : " + value+"\n")


#m_line = m_line - (nlhead + nltail)

################################################ printing the result of the verification of the manifest file #####################################################

FSOR.write(str(m_line) + " lines (SNP) in the manifest \n")

if field_build != None:
	FSOR.write("\nThere are "+str(len(dico_build))+" build versions specified in the Manifest file (Checked only for the SNP presents in the bim file)\n")
	for key, value in dico_build.items():
		FSOR.write(key + "\n")
else :
	FSOR.write("\nThere is no build version specified in the Manifest file.\n")

FSOR.write("\nNumber of SNPs with unknown chromosome : "+str(missing_chro)+"\n")
FSOR.write("Number of SNPs with unknown position : "+str(missing_pos)+"\n")
FSOR.write("Number of SNPs with unknown chromosome and position: "+str(missing_chropos)+"\n")

FSOR.write("\nNumber of strange chromosome value (if not null, take a look to the error output) : "+str(chro_wrong_value)+"\n")
FSOR.write("Number of strange position value (if not null, take a look to the error output) : "+str(pos_wrong_value)+"\n")

FSOR.write("\nChromosomes that are not autosomal (>22):\n")
for key, value in manifest_chr.items():
	FSOR.write("Number of SNPs with the chromosome code '"+str(key)+"': "+ str(value) +"\n")


FSOR.write("\nNumber of SNPs with both allele codes as A, T, C, G: "+str(couple_good_allele_manifest) +"\n")
FSOR.write("Number of SNPs with at least one missing allele (if only one is missing, the other is A, T, C or G: " + str(couple_missing_allele_manifest)+"\n")
FSOR.write("Number of SNPs with at least one allele code different of from A, T, C, G, 0:" + str(couple_unusual_allele_manifest) +"\n")

FSOR.write("\nNumber of allele codes A, for the first allele: "+str(compteur_allele1_manifest["A"]) + ", for the second allele: " + str(compteur_allele2_manifest["A"]) + "\n")
FSOR.write("Number of allele codes T, for the first allele: " + str(compteur_allele1_manifest["T"]) + ", for the second allele: " + str(compteur_allele2_manifest["T"]) + "\n")
FSOR.write("Number of allele codes G, for the first allele: " + str(compteur_allele1_manifest["G"]) + ", for the second allele: " + str(compteur_allele2_manifest["G"]) + "\n")
FSOR.write("Number of allele codes C, for the first allele: " + str(compteur_allele1_manifest["C"]) + ", for the second allele: " + str(compteur_allele2_manifest["C"]) + "\n")
FSOR.write("Number of SNPs with 0 as allele code (missing), for the first allele: " + str(compteur_allele1_manifest["0"]) + ", for the second allele: " + str(compteur_allele2_manifest["0"]) + "\n\n")


#compteur_bad_allele1=0
for key, value in compteur_allele1_manifest.items():
	if key != 'A' and key != 'T' and key != 'G' and key != 'C' and key != '0':
		FSOR.write("Number of allele "+str(key)+" for the first allele: "+ str(value) +"\n")
		#compteur_bad_allele1 += value

#compteur_bad_allele2=0
for key, value in compteur_allele2_manifest.items():
	if key!="A" and key!="T" and key!="G" and key!="C" and key!="0":
		FSOR.write("Number of allele "+str(key)+" for the second allele: "+ str(value) +"\n")
		#compteur_bad_allele2 += value


FSOR.write("\n\nList of the first SNPs from the Manifest file with an allele code different from A, T, C, G or 0:")
if liste_diff_alleles_manifest != []:
	for index in range(len(liste_diff_alleles_manifest)):
		FSOR.write("\n"+ liste_diff_alleles_manifest[index][0] + "\t" + liste_diff_alleles_manifest[index][1] + "\t" + liste_diff_alleles_manifest[index][2])
else :
	FSOR.write("\nNone\n")

FSOR.write("\n\nNumber of mitochondrial SNPs (starting with 'Mito') that are not on the chr 26 or 0 in the Manifest file: " + str(compteur_mito_wrong_manifest) + "\n")
FSOR.write("\nList of the first mitochondrial SNPs (starting with 'Mito') that are not on the chr 26 or 0 in the Manifest file:")
if liste_wrong_mito_manifest != []:
	for index in range(len(liste_wrong_mito_manifest)):
		FSOR.write("\n"+ liste_wrong_mito_manifest[index][0] + " on chr : " +  liste_wrong_mito_manifest[index][1])
else :
	FSOR.write("\nNone")

def print_double(double_rs_man, liste_double_rs_man, double_chrpos_man, double_same_allele_manifest, liste_double_chrpos_man, double_rschrpos_man, liste_double_rschrpos_man):
	FSOR.write("\n\n**************Identification of duplicates*********\n\n***TYPE 1:\nNumber of SNPs that appear at least twice in the bim file with same chr:pos and same allele codes (whatever their SNP name): " + str(double_same_allele_manifest) +"\n")
	
	FSOR.write("\n***TYPE 2:\nNumber of SNPs (same name, chr and pos) that appear at least twice with different allele codes: " + str(double_rs_man) +"\n")
	FSOR.write("\nList of the first SNPs Type 2:\n")
	if liste_double_rs_man != [] :
		for index in range(len(liste_double_rs_man)):
			for i in range(len(liste_double_rs_man[index])):
				FSOR.write(liste_double_rs_man[index][i] +"\t")
			FSOR.write("\n")
	else:
		FSOR.write("None\n")

	FSOR.write("\n***TYPE 3:\nNumber of SNPs (same chr and pos but different names) that appear at least twice and with different allele codes: " + str(double_chrpos_man) +"\n")
	FSOR.write("\nList of the first SNPs Type 3:\n")
	if liste_double_chrpos_man != [] :
		for index in range(len(liste_double_chrpos_man)):
			for i in range(len(liste_double_chrpos_man[index])):
				FSOR.write(liste_double_chrpos_man[index][i] +"\t")
			FSOR.write("\n")
	else:
		FSOR.write("None\n")

	FSOR.write("\n***TYPE 4:\nNumber of SNPs that appear at least twice with the same name but with different chr pos (whatever their allele codes are): " + str(double_rschrpos_man) +"\n")
	FSOR.write("\nList of the first SNPs Type 4:\n")
	if liste_double_rschrpos_man != [] :
		for index in range(len(liste_double_rschrpos_man)):
			for i in range(len(liste_double_rschrpos_man[index])):
				FSOR.write(liste_double_rschrpos_man[index][i] +"\t")
			FSOR.write("\n")
	else:
		FSOR.write("None\n")

if field_strandP != None:
	print_double(double_rs_man_p, liste_double_rs_man_p, double_chrpos_man_p, double_same_allele_manifest_p, liste_double_chrpos_man_p, double_rschrpos_man_p, liste_double_rschrpos_man_p)
elif field_strandF!=None:
	print_double(double_rs_man_f, liste_double_rs_man_f, double_chrpos_man_f, double_same_allele_manifest_f, liste_double_chrpos_man_f, double_rschrpos_man_f, liste_double_rschrpos_man_f)
elif field_strandT!=None:
	print_double(double_rs_man_t, liste_double_rs_man_t, double_chrpos_man_t, double_same_allele_manifest_t, liste_double_chrpos_man_t, double_rschrpos_man_t, liste_double_rschrpos_man_t)
else:
	sys.stderr.write("\n\n**********************************************************\n\nThere isn't any strand nomenclature specified !!!!\nIt should not happened (this error message is line 1256, the problem should come from the parameters at the begining of script pre-processing.py\n\n**********************************************************\n\n")
	exit(0)


if field_strandP!=None:
	FSOR.write("\nNumber of SNPs on the strand - in the manifest (Nomenclature +/-): " +str(man_minus_strand_p))
else:
	FSOR.write("\nWe don't have the information on the strand in the nomenclature +/- in this Manifest.")

if field_strandF!=None:
	FSOR.write("\nNumber of SNPs on the strand Rev in the manifest (Nomenclature For/Rev): " +str(man_minus_strand_f))
else:
	FSOR.write("\nWe don't have the information on the strand in the nomenclature For/Rev in this Manifest.")

if field_strandT!=None:
	FSOR.write("\nNumber of SNPs on the strand Bot in the manifest (Nomenclature Top/Bot): " +str(man_minus_strand_t))
else:
	FSOR.write("\nWe don't have the information on the strand in the nomenclature Top/Bot in this Manifest.")


#FSOR.write("\n\n\n !!!!! Be careful : the next result are important !!!!!\n")
FSOR.write("\n\n\nNumber of Characters other than + or - in the field that contains the strand in nomenclature +/- :\n")
if field_strandP != None and dico_char_PM != {} :
	for key,value in dico_char_PM.items():
		FSOR.write("Charater : '"+str(key)+"' found "+str(value)+" time\n")
else :
	FSOR.write("None\n")

FSOR.write("\nNumber of Characters other than F or R in the field that contains the strand in nomenclature Forward/Reverse :\n")
if field_strandF != None and dico_char_FR != {} :
	for key,value in dico_char_FR.items():
		FSOR.write("Charater : '"+str(key)+"' found "+str(value)+" time\n")
else :
	FSOR.write("None\n")

FSOR.write("\nNumber of Characters other than TOP, BOT, PLUS or MINUS in the field that contains the strand in nomenclature Top/Bot :\n")
if field_strandT != None and dico_char_TB != {} :
	for key,value in dico_char_TB.items():
		FSOR.write("Charater : '"+str(key)+"' found "+str(value)+" time\n")
else :
	FSOR.write("None\n")



#################################### Printing result of the comparaison betwen the bim and the manifest file #############################################

def print_once(number_find, number_missing, number_ambiguous, number_imp, number_notation):
	number_find_print = number_find
	number_missing_print = number_missing
	number_ambiguous_print = number_ambiguous 
	number_imp_print = number_imp
	number_notation_print = number_notation
	return number_find_print, number_missing_print, number_ambiguous_print, number_imp_print, number_notation_print

if field_strandP!=None:
	number_find_print, number_missing_print, number_ambiguous_print, number_imp_print, number_notation_print = print_once(number_find_p ,number_missing_p ,number_ambiguous_p ,number_imp_p, number_notation_p)
elif field_strandF!=None:
	number_find_print, number_missing_print, number_ambiguous_print, number_imp_print, number_notation_print = print_once(number_find_f ,number_missing_f ,number_ambiguous_f ,number_imp_f, number_notation_f)
else:
	number_find_print, number_missing_print, number_ambiguous_print, number_imp_print, number_notation_print = print_once(number_find_t ,number_missing_t ,number_ambiguous_t ,number_imp_t, number_notation_t)

FSOR.write(str("\n\n---------------- Strand check ---------------------- \n\n"))
# name field in the bim match whith the name field in the manifest
FSOR.write(str(number_find_print)+" of the SNPs in the bim file (on "+str(nb_line)+") are found in the Manifest file (linkage using SNP name)\n")

if field_strandP!=None:
	FSOR.write("\nNumber of SNPs from the Bim file that are defined on the strand '-' in the Manifest file (whatever the allele codes reported in the Bim file): " +str(man_minus_strand_found_p))
else:
	FSOR.write("\nWe don't have the information on the strand in the notation +/-.")

if field_strandF!=None:
	FSOR.write("\nNumber of SNPs from the Bim file that are defined on the strand 'Rev' in the Manifest file (whatever the allele codes reported in the Bim file): " +str(man_minus_strand_found_f))
else:
	FSOR.write("\nWe don't have the information on the strand in the notation For/Rev.")

if field_strandT!=None:
	FSOR.write("\nNumber of SNPs from the Bim file that are defined on the strand 'Bot' in the Manifest file (whatever the allele codes reported in the Bim file): " +str(man_minus_strand_found_t))
else:
	FSOR.write("\nWe don't have the information on the strand in the notation Top/Bot.")


FSOR.write("\n\nDetection of strand for the SNPs in the bim file:")
if field_strandP != None:
	FSOR.write("\n\nSNPs in strand +: "+str(number_pos)+"\tSNP with undetermined strand but reported in + in the Manifest: "+str(number_imp_pos)+"\n")
	FSOR.write("SNPs in strand -: "+str(number_neg)+"\tSNP with undetermined strand but reported in - in the Manifest: "+str(number_imp_neg)+"\n")
	FSOR.write("SNPs in a different strand between the Bim and the Manifest file: "+str(nb_diff_strand_BM_p)+"\n")
	FSOR.write("SNPs with unknown charactere instead of nomenclature +/- in the Manifest file: "+str(man_unknown_strand_p)+"\n")
else:
	FSOR.write("\n\nWe don't have the information on the strand in the nomenclature +/-.\n")

if field_strandF != None:
	FSOR.write("\nSNPs in strand For: "+str(number_for)+"\tSNP with undetermined strand but reported in For in the Manifest: "+str(number_imp_for)+"\n")
	FSOR.write("SNPs in strand Rev: "+str(number_rev)+"\tSNP with undetermined strand but reported in Rev in the Manifest: "+str(number_imp_rev)+"\n")
	FSOR.write("SNPs in a different strand between the Bim and the Manifest file: "+str(nb_diff_strand_BM_f)+"\n")
	FSOR.write("SNPs with unknown charactere instead of nomenclature For/Rev in the Manifest file: "+str(man_unknown_strand_f)+"\n")
else:
	FSOR.write("\nWe don't have the information on the strand in the nomenclature For/Rev.\n")
if field_strandT != None:
	FSOR.write("\nSNPs in strand Top: "+str(number_top)+"\tSNP with undetermined strand but reported in Top in the Manifest: "+str(number_imp_top)+"\n")
	FSOR.write("SNPs in strand Bot: "+str(number_bot)+"\tSNP with undetermined strand but reported in Bot in the Manifest: "+str(number_imp_bot)+"\n")
	FSOR.write("SNPs in a different strand between the Bim and the Manifest file: "+str(nb_diff_strand_BM_t)+"\n")
	FSOR.write("SNPs with unknown charactere instead of nomenclature Top/Bot in the Manifest file: "+str(man_unknown_strand_t)+"\n")
else:
	FSOR.write("\nWe don't have the information on the strand in the nomenclature Top/Bot.\n")


FSOR.write("\n\nDifferent cases where we are not able to define SNP strand:\n")
FSOR.write("\nNumber of SNPs with missing information on both allele codes (in the Bim or Manifest file): "+str(number_missing_print)+"\n")
FSOR.write("\nNumber of ambiguous SNPs in the bim file and in the Manifest file: " + str(number_ambiguous_print) + "\n")

FSOR.write("\nNumber of SNPs from the Bim file found in the Manifest file with allele codes different from A, T, C, G, 0: " + str(number_notation_print) + "\n")


FSOR.write("\nNumber of SNPs with different allele codes in the Bim file and in the Manifest file: " + str(number_imp_print) + "\n")
FSOR.write("\nList of the first SNPs with different allele codes in the Bim file and in the Manifest file:")
if liste_diff_bim_manifest_p != []:
	for index in range(len(liste_diff_bim_manifest_p)):
		FSOR.write("\n"+ liste_diff_bim_manifest_p[index])
else :
	FSOR.write("\nNone")


FSOR.write("\n\nNumber of SNPs from the Bim file found in the Manifest with a chr:pos different from the information indicated in the Manifest file: " + str(number_diff_pos_BM) +"\n")
FSOR.write("\nList of the first SNPs with different chr:pos in the Bim and Manifest file:\n")
if liste_diff_pos_BM != []:
	for index in range(len(liste_diff_pos_BM)):
		for i in range(len(liste_diff_pos_BM[index])):
			FSOR.write(liste_diff_pos_BM[index][i]+"\t")
		FSOR.write("\n")
else:
	FSOR.write("None\n")


################################ verification of the fam file ################################

FSOR.write(str("\n\n********************************************************** \n"))
FSOR.write(str("********************** Check fam file ******************** \n"))
FSOR.write(str("********************************************************** \n \n"))
FSOR.write(str("---------------- Check the availability of family ID ---------------------- \n"))
f_fam=open(fam)# we open the fam file

compteur_ID_family=0 #nb of family's ID
compteur_ID_indi=0 #nb of subject's ID
compteur_ID_father=0 #nb of father's ID
compteur_ID_mother=0 #nb of mother's ID
compteur_men=0 #nb of men
compteur_women=0 #nb of women
compteur_nosex=0 #if the sex is not informed
compteur_control=0 #controls
compteur_case=0 #cases
compteur_miss=0 #missing infos
compteur_pb = 0 #if the info is not either 1, 2, 0 or -9

for line in f_fam: #when we read the file
	try :
		fields=line.split(" ") #get
		ID_family=fields[0] #we get the family ID
		ID_indi=fields[1] #we get the individual ID
	except:
		fields=line.split("\t") #get:
		ID_family=fields[0] #we get the family ID
		ID_indi=fields[1] #we get the individual ID        
        sexe=fields[4]#we get the sexe
        pheno=fields[5].strip()#we look if informations is missing

        #we check the family id
        if ID_family =="0":
                compteur_ID_family +=1

        #we check the individual id
        if ID_indi =="0":
                compteur_ID_indi +=1

        #we check the sex
        if sexe =="1":
                compteur_men +=1
        elif sexe =="2":
                compteur_women +=1
        else:
                compteur_nosex +=1

        #we check the missing informations or cases/controls
        if pheno =="1":
                compteur_control +=1
        elif pheno =="2":
                compteur_case +=1
	elif pheno == "0" or pheno == "-9":
		compteur_miss += 1
        else: #if there are some pbs
                compteur_pb +=1


#we write the results in the output file
if compteur_ID_family !=0:
        FSOR.write("This file has unknown family ID, number: "+str(compteur_ID_family)+"\n")
else:
        FSOR.write("All family IDs are completed \n")


FSOR.write(str("\n---------------- Check the availability of individual ID ---------------------- \n"))
if compteur_ID_indi !=0:
        FSOR.write("This file has unknown individual ID, number: "+str(compteur_ID_indi)+"\n")
else:
        FSOR.write("All individual IDs are completed \n")


FSOR.write(str("\n------------------ Check the availability of sex -------------------------- \n"))
if compteur_men !=0:
        FSOR.write("This file has "+str(compteur_men)+" men\n")

if compteur_women !=0:
        FSOR.write("This file has "+str(compteur_women)+" women \n")

if compteur_nosex !=0:
        FSOR.write("This file has unknown sex, number: "+str(compteur_nosex)+"\n")


FSOR.write(str("\n----------- Check the availability of phenotype information ---------------- \n"))
if compteur_control !=0:
        FSOR.write("This file has "+str(compteur_control)+" controls \n")

if compteur_case !=0:
        FSOR.write("This file has "+str(compteur_case)+" cases \n")

if compteur_miss !=0:
        FSOR.write("This file contains missing phenotype information, number: "+str(compteur_miss)+" \n")

if compteur_pb !=0:
        FSOR.write("This file contains problematic phenotype information, number: "+str(compteur_pb)+" \n")

f_fam.close()



##########################################################################
######################## Start the success rate ##########################
##########################################################################
FSOR.write(str("\n********************************************************** \n"))
FSOR.write(str("********************* Success rate *********************** \n"))
FSOR.write(str("********************************************************** \n \n"))

limit = 15

###Test with 95%###
FSOR.write(str("---------------- success rate at 95% ---------------------------- \n"))
if compteur_chro_sup_27==0: #if chro <27
	os.system(plink+" --bfile "+args.file+" --mind 0.05 --make-bed --out outQC") #we use plink with a shell command
else: #if we have some chro >= 27
	os.system(plink+" --bfile "+args.file+" --mind 0.05 --chr-set "+str(chro_pasbeau)+" no-x no-y no-xy no-mt --make-bed --out outQC") #we use plink with a shell command

if os.path.exists("./outQC.log"):
        log=open("outQC.log")
        for line in log:
                if "loaded" in line or "rate" in line or "removed" in line or "pass" in line:
			FSOR.write(line)

liste_ID=""
if os.path.exists("./outQC.irem"): #if plink creates a file for people removed
	irem=open("outQC.irem")
	display = 0
	for line in irem:
		display += 1
		if display == limit:
			break
		ID=line.split("\t")[0] #get the id and write them in the file
		liste_ID=liste_ID+"\n"+ID
	FSOR.write("List of the first subjects ID who don't pass the test: \n"+liste_ID)
	liste_ID=""

FSOR.write("\n \n")
os.system('rm outQC*') #and we delete all the files than plink created

###Test with 90%###
FSOR.write(str("---------------- success rate at 90% ---------------------------- \n"))
if compteur_chro_sup_27==0:
	os.system(plink+" --bfile "+args.file+" --mind 0.1 --make-bed --out outQC")
else:
	os.system(plink+" --bfile "+args.file+" --mind 0.1 --chr-set "+str(chro_pasbeau)+" no-x no-y no-xy no-mt --make-bed --out outQC") #we use plink with a shell command

if os.path.exists("./outQC.log"):
        log=open("outQC.log")
        for line in log:
                if "loaded" in line or "rate" in line or "removed" in line or "pass" in line:
			FSOR.write(line)

if os.path.exists("./outQC.irem"):
        irem=open("outQC.irem")
	display = 0
        for line in irem:
		display += 1
		if display == limit:
			break
                ID=line.split("\t")[0]
                liste_ID=liste_ID+"\n"+ID
	FSOR.write("List of the first subjects ID who don't pass the test: \n"+liste_ID)
	liste_ID=""

FSOR.write("\n \n")
os.system('rm outQC*')

###Test with 80%###
FSOR.write(str("---------------- success rate at 80% ---------------------------- \n"))
if compteur_chro_sup_27==0:
	os.system(plink+" --bfile "+args.file+" --mind 0.2 --make-bed --out outQC")
else:
	os.system(plink+" --bfile "+args.file+" --mind 0.2 --chr-set "+str(chro_pasbeau)+" no-x no-y no-xy no-mt --make-bed --out outQC") #we use plink with a shell command

if os.path.exists("./outQC.log"):
	log=open("outQC.log")
	for line in log:
		if "loaded" in line or "rate" in line or "removed" in line or "pass" in line:
			FSOR.write(line)


if os.path.exists("./outQC.irem"):
        irem=open("outQC.irem")
	display = 0
        for line in irem:
		display += 1
		if display == limit:
			break
                ID=line.split("\t")[0]
                liste_ID=liste_ID+"\n"+ID
	FSOR.write("List of subjects ID who don't pass the test: \n"+liste_ID)
	liste_id=""

FSOR.write("\n \n")
os.system('rm outQC*')



##########################################################################
######################## Start the Sex Check##############################
##########################################################################
FSOR.write(str("********************************************************** \n"))
FSOR.write(str("*********************** Sex Check ************************ \n"))
FSOR.write(str("********************************************************** \n \n"))

def impute_sex(compteur_chro_sup_27,chro_pasbeau, nb_people, FSOR):
	if compteur_chro_sup_27==0:
		os.system(plink+" --bfile "+args.file+" --check-sex --out outQC") #we use plink with a shell command
	else:
		os.system(plink+" --bfile "+args.file+" --chr-set "+str(chro_pasbeau)+" no-x no-y no-xy no-mt --check-sex --out outQC")

	if os.path.exists("./outQC.sexcheck"): #if plink create a file when the file has a problem with the sex
		os.system("grep 'PROBLEM' outQC.sexcheck > outQC.problem") #i filtre the result
		sexcheck=open("outQC.problem")
		FSOR.write(str("FID")+"\t"+str("IID")+"\t"+str("PEDSEX")+"\t"+str("SNPSEX")+"\t"+str("STATUS")+"\t"+str("F")+"\n")
		display = 0
		for line in sexcheck:
			display += 1
			if display < limit:
				FSOR.write(str(line)) #and i just write the people in the file who have a problem
		FSOR.write("\nNumber of subjects with sex problem: "+str(display)+" on "+str(nb_people)+"\n")

	if os.path.exists("./outQC.log"): #if plink create a file when we have an error
		log=open("outQC.log")
		for line in log:
			if "Error" in line:
				FSOR.write(line)

	FSOR.write("\n\n")
	os.system('rm outQC*') #and i delete all files create by plink

if (compteur_chro_23 + compteur_chro_24 + compteur_chro_25 + compteur_chro_26) != 0:
	impute_sex(compteur_chro_sup_27,chro_pasbeau, nb_people, FSOR)
else:
	FSOR.write("We can not check sex because there is no SNPs on sexual chromosomes.\n\n")

##########################################################################
###################### Start the Genome Version ##########################
##########################################################################
FSOR.write(str("\n********************************************************** \n"))
FSOR.write(str("********************** Genome Version ******************** \n"))
FSOR.write(str("********************************************************** \n \n"))

os.system("/data/Epic_central_genetics/work/Scripts/Normalisation/Preprocessing/check_RefBuild.sh "+bim+" > temp.txt")

with open("temp.txt",'r') as genome:
	for line in genome:
		FSOR.write(line+"\n")


# We close all file

keepOne_key = ""
for key, value in dico_build.items():
	value.close()
	keepOne_key = key

if len(dico_build) == 1:
	os.system('rm '+cancer+'_SNP_build_'+keepOne_key+'.txt')


FSOR.close()
FSMAN.close()
if field_strandP != None:
	FSSTRANDp.close()
	FSMINUS.close()
	FSUNSTRANDp.close()
if field_strandF != None:
	FSSTRANDf.close()
	FSREV.close()
	FSUNSTRANDf.close()
if field_strandT != None:
	FSSTRANDt.close()
	FSBOT.close()
	FSUNSTRANDt.close()
FSDIFPOS.close()
FSDIFALL.close()
FSNOTFOUND.close()
FSMITO.close()
FSDUPLSAME.close()
FSDUPLNAME.close()
FSDUPLPOS.close()
FSDUPLERROR.close()


##########################################################################
########################## print QC summary #############################
##########################################################################

if print_summary != 'n':
	with open(summary_file, 'a') as sum_file :
		sum_file.write("\nNumber of subjects in the Fam file: " + str(nb_people) + "\n")
		sum_file.write("Number of SNPs in the Bim file: " + str(nb_SNP) + "\n")
		
		if print_summary == '1' or print_summary == 'a':
			#sum_file.write("\n*** Usefull information from the first QC ***\n")
	
			if field_build != None:# print build found in the Manifest
				sum_file.write("\nManifest build: "+str(len(dico_build))+" build version(s) specified in the Manifest file (Checked only for the SNP present in the Bim file)\n")
				for key, value in dico_build.items():
					sum_file.write(key + "\n")
				if len(dico_build) > 1:
					sum_file.write("\nWarning !!!\nMaybe SNP from some of those build need to be deleted.\n")
			else :
				sum_file.write("\nThere is no build version specified in the Manifest file.\n")
		
			sum_file.write("\nDetected build in the Bim file:\n")
			with open("temp.txt",'r') as genome:
				for line in genome: # print build detected in the Bim file
					sum_file.write(line)

			sum_file.write("\nNumber of SNPs from the Bim file not found in the Manifest file: " + str(nb_SNP - number_find_print)+"\n")
	
			sum_file.write("\nNumber of SNPs with different chr:pos information between the Bim file and the Manifest file: " + str(number_diff_pos_BM) +"\n")
	
			sum_file.write("Number of SNPs with different allele codes between the Bim file and the Manifest file: " + str(number_imp_print) + "\n")
	
			sum_file.write("\n\nDetection of strand for the SNPs in the Bim file:\n")

			if field_strandP != None:
				nb_diff_strand_BM = nb_diff_strand_BM_p
			elif field_strandF != None:
				nb_diff_strand_BM = nb_diff_strand_BM_f
			else :
				nb_diff_strand_BM = nb_diff_strand_BM_t
			sum_file.write("\n- SNPs in a different strand between the Bim and the Manifest file: "+str(nb_diff_strand_BM)+"\n")
	
			if field_strandP != None:
				sum_file.write("\n- SNPs in strand +: "+str(number_pos)+"\n")
				sum_file.write("- SNPs in strand -: "+str(number_neg)+"\n")
				sum_file.write("- SNPs with characters other than +/- in the Manifest file: "+str(man_unknown_strand_p)+"\n")
			else:
				sum_file.write("\n- We don't have the information on the strand in the nomenclature +/-.\n")
	
			if field_strandF != None:
				sum_file.write("\n- SNPs in strand For: "+str(number_for)+"\n")
				sum_file.write("- SNPs in strand Rev: "+str(number_rev)+"\n")
				sum_file.write("- SNPs with characters other than F/R in the Manifest file: "+str(man_unknown_strand_f)+"\n")
			else:
				sum_file.write("\n- We don't have the information on the strand in the nomenclature For/Rev.\n")
			if field_strandT != None:
				sum_file.write("\n- SNPs in strand Top: "+str(number_top)+"\n")
				sum_file.write("- SNPs in strand Bot: "+str(number_bot)+"\n")
				sum_file.write("- SNPs with characters other than Top/Bot in the Manifest file: "+str(man_unknown_strand_t)+"\n\n")
			else:
				sum_file.write("\n- We don't have the information on the strand in the nomenclature Top/Bot.\n\n")
	
		if print_summary == '1' or print_summary == '3' or print_summary == 'a':
			#sum_file.write("\n\n*** Information to compare between the first and the third QC ***\n\n")
	
			sum_file.write("\nNumber of SNPs with missing information on position or chromosome: "+str(compteur_missing_pos)+"\n")

			sum_file.write("Number of SNPs with at least one missing allele code (if only one is missing, the other is A, T, C or G): " + str(couple_missing_allele)+"\n")

	
		if print_summary == '2' or print_summary == 'a':
			#sum_file.write("\n\n*** Usefull information from the second QC ***\n\n")
			
			sum_file.write("\nNumber of women: " + str(compteur_women) + "\nNumber of men: " + str(compteur_men) + "\n\n")

			if print_summary == 'a':
				if (compteur_chro_23 + compteur_chro_24 + compteur_chro_25 + compteur_chro_26) != 0:
					impute_sex(compteur_chro_sup_27,chro_pasbeau, nb_people, sum_file)
				else:
					sum_file.write("We can not impute sex because there is no SNPs on sexual chromosomes.\n\n")

			sum_file.write("Number of mitochondrials SNPs (name starts by 'Mito') not located on chromosome 26: "+str(compteur_mito_wrong)+"\n")

			nb_dupl_T1 = file_len(cancer+"_duplicate_T1_same_allele.txt")
			sum_file.write("Number of SNPs duplicates with the same position and the same alleles, whatever their names (T1): " + str(nb_dupl_T1) + "\n")

			nb_dupl_T3 = file_len(cancer+"_duplicate_T3_pos.txt")
			sum_file.write("Number of SNPs duplicates with different names and alleles, but with the same position (T3): " + str(nb_dupl_T3) + "\n")

	
		if print_summary == '3' or print_summary == 'a':
			#sum_file.write("\n\n*** Usefull information from the third QC ***\n\n")

			sum_file.write("\nNumber of SNPs on chromosomes 23, 24, 25 and 26: " + str(compteur_chro_23 + compteur_chro_24 + compteur_chro_25 + compteur_chro_26) + "\n")	

			#sum_file.write("Number of SNPs with unusual allele (they will not be annotate), they usually are DIP: " + str(couple_unusual_allele) + "\n")

		sum_file.write("\n**********************************************************\n\n")

os.system("rm temp.txt")
