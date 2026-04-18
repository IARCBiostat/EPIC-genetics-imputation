#!/usr/bin/env python2.7

#this script has to be use as a module for the script : pre-processing
#it contain function who do verifications done on bim and manifest files, and some of them write in exclusion files.

# this module contain several functions:

################### add_liste_pb(liste, snp, n):################### 
#Fonction that adds the n 1st SNP name in the problem list

################### def file_len(fname):###################
#count the number of line in a file

################### complementary_strand_manifest(allele1, allele2): ######################
#take two parameters : allele1 and allele2 and send back their complementary # exemple : A- give T- or A0 give T0  ....      we think it's not correct for NA to be send back NT

################### count_duplicate(): ###################
#count duplicates the 4 types and number of SNP in one duplicate. It send back the type of the duplicate

###################  #####################

#script written by Manon Knuchel in June 2019

##MODULES
import sys, os, re, datetime, argparse

#declaring complementary allele dictionary, it's also used to verify if the allele are ATCG
complementary = {"A":"T", "T":"A", "C":"G", "G":"C"}

# nb1 et nb2 correspond to the length of error lists (if you want to modify it see in th list declaration below, wich variable use your list)
nb1=10
nb2=30

#Fonction that adds the n 1st SNP name in the problem list
def add_liste_pb(liste, snp, n):
	if len(liste) < n:
		if snp not in liste:
			liste.append(snp)

#Fonction that count the number of line in a file
def file_len(fname):
	with open(fname) as f:
		for i, l in enumerate(f):
			pass
		try:
			i
		except NameError:
			i = -1
		return i+1

#### fonction that send back the complementary alleles, used for the manifest file : (alleles in dico are keep in + / For / Top)######
# !!!!!!!!!!!  be carefull of special case in entry : actually if we have one strange allele and one with ATCG in entry it send back the same strange allele but the complementary for the ATCG one, except for NA !!!!!!!!!!!!!!
# exemple : A- give T- or A0 give T0  ....      we think it's not correct for NA to be send back NT (no allele AB in the Manifest)
def complementary_strand_manifest(allele1, allele2):
	alleleA =""
	alleleB =""
	if allele1 in complementary:
		alleleA = complementary.get(allele1)
	else:
		alleleA = allele1
	if allele2 in complementary:
		alleleB = complementary.get(allele2)
	else:
		alleleB = allele2
	if allele1 == 'N' or allele2 == 'N' : # we don't want to have the complement of the allele A who is with an N
		alleleA = allele1
		alleleB = allele2
	return alleleA+alleleB


#count different type of duplicate :
# 1 : position and allele are equal (whatever the name is equal or not) 			-> real duplicate
# 2 : name and position are equal but alleles are different 					-> tri-allelic with the same name
# 3 : position is equal but alleles are different 						-> tri-allelic with different name
# 4 : same name but positions are different (whatever the alleles are different or equal)	-> they might be error
def count_duplicate(dico_name, dico_pos, name, pos, allele1, allele2):
	alleles = allele1+allele2
	permute = allele2+allele1

	flag_duplicate = 1
	flag_nom_pos = 0
	flag_nom_allele = 0
	flag_pos = 0
	flag_pos_allele = 0
	result = pos.split(":")

	if name in dico_name: # if the name already exist
		#for i in range(len(dico[name])): # i here is an index
		#	if pos == dico[name][i][0]: # we compare the SNP position with all his duplicate
		for snp in dico_name[name]: # for each snp in with the same name
			if pos == snp[0]: # we compare position
				flag_nom_pos = 1
				if alleles == snp[1] or permute == snp[1]: # we compare alleles
					flag_nom_allele = 1
	elif pos in dico_pos and (result[0] != "0" and result[1] != "0"): # if their is already a snp at this position
		flag_pos = 1
		for pos_name in dico_pos[pos]: # I take each name of snp who are at this position
			for snp in dico_name[pos_name]: # for each name, I verify allele of each duplicate
				if snp[1] == alleles or snp[1] == permute: 
					flag_pos_allele = 1 
	else:# it's not a duplicate
		flag_duplicate = 0
	
	if flag_duplicate == 1: # it's a duplicate
		if flag_nom_allele == 1 or flag_pos_allele == 1: # the snp pass all control and they it have the same allele with at least one other duplicate
			return 1 # Duplicate Type 1
		elif flag_nom_pos == 1 and flag_nom_allele == 0: # same name and pos, but alleles different
			return 2 # Duplicate Type 2
		elif flag_pos == 1 and flag_pos_allele == 0: # no identical name but same pos and different alleles
			return 3 # Duplicate Type 3
		elif flag_nom_pos == 0 and flag_pos == 0: # same name but different pos (whatever alleles are)
			return 4 # Duplicate Type 4
		else : 
			sys.stderr.write("\n\n**********************************************************\n\nA case is missing in preProcessing_functions/checks_listing.py (line 101) in count_duplicate !\n\n**********************************************************\n\n")
			return -1
	else: # it's not a duplicate
		return 0 


