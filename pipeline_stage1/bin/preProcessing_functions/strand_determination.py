##!/usr/bin/env python2.7

#this script has to be use as a module for the script : pre-processing
#it do the research strand of SNP from the bim file by comparing it with the manifest file.

# this module contain several functions:
# research_strand()
#
#

#script written by Manon Knuchel in June 2019

##MODULES
import sys, os, re, datetime, argparse

#My module:
from .checks_listing import * ##importing all function of the module cheks_listing who is in the same folder as the curent module

######## fonction who count the strand in the bim file by comparing it to the manifest file. (it need manifest allele on the +/For/Top strand, and the information on the original strand of the SNP in the manifest file) ########
def research_strand(rs, snp, allele1, allele2, nb_find, nb_p, nb_n, nb_amb, nb_mis, nb_imp, nb_imp_pos, nb_imp_neg, nb_diff_strand_BM, nb_not, nb_neg_man, nb_unknown_strand, liste_diff_bim_manifest, flag, fichier, fic_diff_allele, fic_unknown_strand):
	nb_find += 1  #SNP find, +1
	if snp[2] == 2:
		nb_unknown_strand += 1
		fic_unknown_strand.write(rs+"\n")

		return nb_find, nb_p, nb_n, nb_amb, nb_mis, nb_imp, nb_imp_pos, nb_imp_neg, nb_diff_strand_BM, nb_not, nb_neg_man, nb_unknown_strand

	if snp[2] == 1:
		nb_neg_man += 1
	 
	allele1 = allele1.strip()
	allele2 = allele2.strip()
	if (allele1 == "0" and allele2 == "0") or (snp[1][0] == "0" and snp[1][1] == "0"): # we don't know on wich strand this SNP because both alleles are missing.
		nb_mis += 1
		if snp[2] == 0: # verify in witch strand the SNP is in the manifest file
			nb_imp_pos += 1
		else:
			nb_imp_neg += 1
	elif (allele1 not in complementary and allele1 != "0") or (allele2 not in complementary and allele2 != "0"): # at least one of the bim file alleles are not ATCG or 0
		nb_not +=1
		if allele1 != "B" and allele2 != "B" and snp[1][0] != "N" and snp[1][1] != "N" :
			if((allele1 != "0") and (allele1 != snp[1][0] and allele1 != snp[1][1])): # if allele are not the same in the bim and the manifest file
				add_liste_pb(liste_diff_bim_manifest, rs, nb1)
				nb_imp += 1
				if flag == 1:
					fic_diff_allele.write(rs+"\n")
			elif ((allele2 != "0") and (allele2 != snp[1][0] and allele2 != snp[1][1])):
				add_liste_pb(liste_diff_bim_manifest, rs, nb1)
				nb_imp += 1	
				if flag == 1:
					fic_diff_allele.write(rs+"\n")
			#else:
			#	FSOR.write("\nNothing to do !!!\n")
		if snp[2] == 0: # we count the strand of the SNP in the manifest
			nb_imp_pos += 1
		else:
			nb_imp_neg += 1
	
	elif allele1 == "0" or allele2 == "0": # one of the alleles is missing
		know = allele1
		if allele2 != "0":
			know = allele2
		if snp[1][0] == complementary.get(snp[1][1]): # if the alleles in the manifest file are ambiguous.
			if know == snp[1][0] or know == snp[1][1]: # if the known allele correspond at an allele of the manifest file
				nb_amb +=1
			else: # the SNP in the manifest is ambiguous, but the allele know in the manifest is different (AT and G0 : alleles are different between the two files
				add_liste_pb(liste_diff_bim_manifest, rs, nb1)
				nb_imp += 1
				if flag == 1:
					fic_diff_allele.write(rs+"\n")
			if snp[2] == 0:# we count the strand of the SNP in the manifest
				nb_imp_pos += 1
			else:
				nb_imp_neg += 1
	
		elif know == snp[1][0] or know == snp[1][1]: # if the allele is equal to one of the manifest, it is one the + strand
			nb_p += 1
			if snp[2] == 1: # count if the strand is not the same in the bim and manifest file
				nb_diff_strand_BM += 1
		elif know == complementary.get(snp[1][0]) or know == complementary.get(snp[1][1]): # if the allele is the complementary to one of the manifest, it is one the - strand
			nb_n += 1
			fichier.write(rs+"\n")
			if snp[2] == 0:  # count if the strand is not the same in the bim and manifest file
				nb_diff_strand_BM += 1	
		else: # the allele of the bim file doesn't mach with one of the manifest, their is an incoherence
			add_liste_pb(liste_diff_bim_manifest, rs, nb1)
			nb_imp += 1
			if flag == 1:
				fic_diff_allele.write(rs+"\n")
			if snp[2] == 0:# we count the strand of the SNP in the manifest
				nb_imp_pos += 1
			else:
				nb_imp_neg += 1
	else: # both allele are present
		if snp[1][0] == complementary.get(snp[1][1]): # alleles in the manifest are ambiguous
			if allele1 == complementary.get(allele2) and (allele1 == snp[1][0] or allele1 == snp[1][1]): # alleles in the bim file are ambiguous and are the same as the manifest
				nb_amb +=1
			else : # alleles in the bin and manifest file are different
				add_liste_pb(liste_diff_bim_manifest, rs, nb1)
				nb_imp += 1
				if flag == 1:
					fic_diff_allele.write(rs+"\n")
			if snp[2] == 0:# we count the strand of the SNP in the manifest
				nb_imp_pos += 1
			else:
				nb_imp_neg += 1

		elif allele1 == snp[1][0]: # if the first allele in bim is equal to the first allele in the manifest
			if allele2 == snp[1][1]: #the second allele is equal, so the strand is +
				nb_p += 1
				if snp[2] == 1: # count if the strand is not the same in the bim and manifest file
					nb_diff_strand_BM += 1
			else : # there is an incoherence between the alleles of the files
				add_liste_pb(liste_diff_bim_manifest, rs, nb1)
				nb_imp += 1
				if flag == 1:
					fic_diff_allele.write(rs+"\n")
				if snp[2] == 0: # we count the strand of the SNP in the manifest
					nb_imp_pos += 1
				else:
					nb_imp_neg += 1
		elif allele1 == snp[1][1]: # if the first allele in bim is equal to the second allele in the manifest
			if allele2 == snp[1][0]: #the second allele is equal, so the strand is +
				nb_p += 1
				if snp[2] == 1: # count if the strand is not the same in the bim and manifest file
					nb_diff_strand_BM += 1
			else : # there is an incoherence between the alleles of the files
				add_liste_pb(liste_diff_bim_manifest, rs, nb1)
				nb_imp += 1
				if flag == 1:
					fic_diff_allele.write(rs+"\n")
				if snp[2] == 0: # we count the strand of the SNP in the manifest
					nb_imp_pos += 1
				else:
					nb_imp_neg += 1

		elif allele1 == complementary.get(snp[1][0]):  # if the first allele in bim is the complementary of the first allele in the manifest
			if allele2 == complementary.get(snp[1][1]): #the second allele is the complementary, so the strand is -
				nb_n += 1
				fichier.write(rs+"\n")
				if snp[2] == 0: # count if the strand is not the same in the bim and manifest file 
					nb_diff_strand_BM += 1	
			else: # the second allele is not the same as the manifest
				add_liste_pb(liste_diff_bim_manifest, rs, nb1)
				nb_imp += 1
				if flag == 1:
					fic_diff_allele.write(rs+"\n")
				if snp[2] == 0: # we count the strand of the SNP in the manifest
					nb_imp_pos += 1
				else:
					nb_imp_neg += 1

		elif allele1 == complementary.get(snp[1][1]):  # if the first allele in bim is the complementary of the second allele in the manifest
			if allele2 == complementary.get(snp[1][0]): #the second allele is the complementary, so the strand is -
				nb_n += 1
				fichier.write(rs+"\n")
				if snp[2] == 0: # count if the strand is not the same in the bim and manifest file
					nb_diff_strand_BM += 1	
			else: # the second allele is not the same as the manifest
				add_liste_pb(liste_diff_bim_manifest, rs, nb1)
				nb_imp += 1
				if flag == 1:
					fic_diff_allele.write(rs+"\n")
				if snp[2] == 0: # we count the strand of the SNP in the manifest
					nb_imp_pos += 1
				else:
					nb_imp_neg += 1

		else: # allele 1 is not 1 of the 2 allele from the manifest file, neither 1 of the 2 complmentary allele of the manifest.
			add_liste_pb(liste_diff_bim_manifest, rs, nb1)
			nb_imp += 1
			if flag == 1:
				fic_diff_allele.write(rs+"\n")
			if snp[2] == 0:  # we count the strand of the SNP in the manifest
				nb_imp_pos += 1
			else:
				nb_imp_neg += 1

				
	return nb_find, nb_p, nb_n, nb_amb, nb_mis, nb_imp, nb_imp_pos, nb_imp_neg, nb_diff_strand_BM, nb_not, nb_neg_man, nb_unknown_strand

