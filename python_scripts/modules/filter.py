#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#   LOncoG : a software for Longitudinal OncoGenomics analysis
#   Authors: S. Besseau, G. Siekaniec, W. Gouraud
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

from pathlib import Path
import copy
import logging
import math
import os
import re
import subprocess
import sys
from copy import deepcopy
from itertools import takewhile
from tqdm import tqdm
from python_scripts.manage_parameters import list_colors
from python_scripts.reusable_functions.check_files import verif_input_vcf, verif_output
from python_scripts.manage_parameters import manage_parameters

def line_without_fail(fail,record, to_suppr : str, intIndex : dict, strLotusFilterCode : str, strFilterPassKey : str, strFilterPassOk : str):
	'''
	Takes a record and deletes all variants that do not pass the filters, if no variants pass the filters returns False otherwise returns True
	Input : record and list of variants to suppressed
	Output : True (Record modified) or False (No record left)
	'''

	# Get the number of alternative variants
	nb_alt = len(record[intIndex['Alt']].split(','))

	# Suppression of variants that don't pass LOncoG filters
	alternative = [alt for i, alt in enumerate(record[intIndex['Alt']].split(',')) if i not in to_suppr]
	if alternative == []:
		return False
	else:
		record[intIndex['Alt']] = ','.join(alternative)

	# Filter modification
	if not fail:
		record[intIndex['Filter']] = 'PASS'


	# Suppresion of the useless modification
	informations = [[info.split('=')[0], info.split('=')[1]] if (len(info.split('=')) == 2) else tuple([info.split('=')[0],'']) for info in [info for info in record[intIndex['Info']].split(';')]]

	for j, information in enumerate(informations):
		id, info = information
		if id == 'OTHER_FILTER':
			informations[j][1]=strFilterPassOk
		
		elif id == 'AS_FilterStatus':
			informations[j][1]='|'.join([status for i, status in enumerate(info.split('|')) if i not in to_suppr])

		elif id == 'AS_SB_TABLE':
			informations[j][1]='|'.join([fr for i, fr in enumerate(info.split('|')) if i-1 not in to_suppr])
	
		else:
			if len(info) > 1 and not info.isnumeric():
				if len(info.split(',')) == nb_alt:
					informations[j][1] = ','.join([elmt for i, elmt in enumerate(info.split(',')) if i not in to_suppr])
				elif len(info.split(',')) == nb_alt+1:
					informations[j][1] = ','.join([elmt for i, elmt in enumerate(info.split(',')) if i-1 not in to_suppr])
	record[intIndex['Info']] = ';'.join(['='.join([id, info]) for id, info in informations])

	names = record[intIndex['Format']].split(':')
	
	values = [value for value in record[intIndex['Values']].split(':')]
	for i, v in enumerate(values):
		if names[i] == 'GT':
			values[i] = '/'.join([str(i) for i in range(nb_alt-len(to_suppr)+1)]) 
		elif len(v.split(',')) == nb_alt:
			values[i] = '/'.join([char for i, char in enumerate(v.split(',')) if i not in to_suppr])
		elif len(v.split(',')) == nb_alt+1:
			values[i] = '/'.join([char for i, char in enumerate(v.split(',')) if i-1 not in to_suppr])
	record[intIndex['Values']] = ':'.join(values)

	return True
	

def fail_filters(record, dict_para, info: {}, format_values: {}, nb_alt: int, AD: [], filter_qual = 50, filter_mqsbz = 0.5, filter_dp: int = 10, filter_af_population: float = 0.0001):
	'''
	Do variant(s) at a position passes the filters ?
	Input : info, AD and AF fields and the number of variants
	Output : False (don't fail) or a dictionnary containing the failed filters
	'''
	#check if known good variant is not failing

	filter_on_mutation_position = dict_para['filter_on_mutation_location'].upper()
	if 'NO' in filter_on_mutation_position or 'NONE' in filter_on_mutation_position:
		filter_on_mutation_position = False
		mutation_to_save = {'IGR', 'THREE_PRIME_UTR', 'FIVE_PRIME_UTR', 'INTRON', 'FIVE_PRIME_FLANK', 'MISSENSE', 'NONSENSE', 'NONSTOP', 'RNA', 'LINCRNA', 'START_CODON_SNP',
							'DE_NOVO_START_IN_FRAME', 'DE_NOVO_START_OUT_FRAME', 'IN_FRAME_DEL', 'IN_FRAME_INS', 'FRAME_SHIFT_INS', 'FRAME_SHIFT_DEL', 'START_CODON_INS',
							'START_CODON_DEL', 'DE_NOVO_START_IN_FRAME', 'DE_NOVO_START_OUT_FRAME', 'exonic', 'splicing', 'ncRNA', 'ncRNA_intronic', 'ncRNA_exonic', 'UTR5', 'UTR3',
							'intronic', 'upstream', 'downstream', 'intergenic', 'SPLICE_SITE', 'START_CODON_SNP', 'START_CODON_INS', 'START_CODON_DEL', 'IGR'}

	elif 'YES' in filter_on_mutation_position or 'TRUE' in filter_on_mutation_position:
		filter_on_mutation_position = True
		mutation_to_save = {'MISSENSE', 'NONSENSE', 'NONSTOP', 'SILENT', 'RNA', 'LINCRNA', 'START_CODON_SNP', 'DE_NOVO_START_IN_FRAME', 'DE_NOVO_START_OUT_FRAME', 'IN_FRAME_DEL',
							'IN_FRAME_INS', 'FRAME_SHIFT_INS', 'FRAME_SHIFT_DEL', 'START_CODON_INS', 'START_CODON_DEL', 'DE_NOVO_START_IN_FRAME', 'DE_NOVO_START_OUT_FRAME',
							'exonic', 'splicing', 'SPLICE_SITE', 'START_CODON_SNP', 'START_CODON_INS', 'START_CODON_DEL'}
	else:
		filter_on_mutation_position = True
		mutation_to_save = set(filter_on_mutation_position.replace(" ","").split(','))

	fail = {'all': set()}
	for i in range(nb_alt):
		fail[i] = set()

	filter_alt_MBQ = dict_para['min_alt_MBQ']
	filter_ad = dict_para['min_alt_AD']
	max_SOR = float(dict_para['max_SOR'])

	alts = False
	if ',' in record[4]:  # IF SEVERAL VARIANTS ARE FOUND IN THE SAME LINE IN VCF
		failed_variants = []
		passed_variants = []
		n_alts = len(record[4].split(','))
		alts = True
		anno = True
		try:
			beginning = str(record).split('ANNO')[0]
		except:
			print('not annotated with ANNOVAR, developer should find another split in previous lines, to split alleles into 2 new variants!')
			anno = False
		chr = record[0]
		pos = record[1]
		id = record[2]
		ref = record[3]
		if '1016963' in record:
			print('1016963')
		if pos == '1016963' and ref == 'GC':
			print('1016963')
		alt = record[4]
		qual = record[5]
		filter_pass = record[6]
		dp = int(re.search(r"DP=(.*?);", str(record)).group(1))
		mbqs = re.search(r"MBQ=(.*?);", str(record)).group(1).split(',')
		ads = format_values['AD'].split(',')
		try:
			popafs = re.search(r"POPAF=(.*?);", str(record)).group(1)
		except:
			print('not called with Mutect2 : no POPAF')
		alt_fails_ori = ['fails: ']
		if dp < filter_dp:
			alt_fails_ori.append('low_DP')
		l_variants = []
		dict_fails = {}
		o = 0
		for v in alt.split(','):
			alt_fails = alt_fails_ori.copy()
			if l_variants == []:
				variant = str(record).split('ALLELE_END;ANNOVAR_DATE')[0].split('ANNO')[1]
				l_variants.append(variant)
			else:
				variant = str(record).split(l_variants[o - 1])[1].split('ALLELE_END')[1].split('ALLELE_END')[0]
				l_variants.append(variant)
				if len(variant) > 10:
					pass
				else:
					variant = l_variants[o - 1]
			alt_mbq = mbqs[o + 1]
			if float(alt_mbq) < float(filter_alt_MBQ):
				alt_fails_ori.append('low_MBQ')
			try:
				popaf = popafs.split(',')[o]
				if float(popaf) >= float(filter_af_population):
					alt_fails.append('high_POPAF')
			except:
				print('not called with Mutect2 : no POPAF')
			ad = ads[o]
			if float(ad) < float(filter_ad):
				alt_fails.append('low_AD')
			try:
				strandqs = re.search(r"STRANDQ=(.*?);", str(record)).group(1)
				strandq = strandqs.split(',')[o]
				if float(strandq) < float(dict_para['min_STRANDQ']):
					alt_fails.append('low_STRANDQ')
			except:
				try:
					SB_values = format_values['SB']
					SB_values = [int(number) for number in SB_values.split(',')]
					refFw = SB_values[0] + 1
					refRv = SB_values[1] + 1
					altFw = SB_values[2] + 1
					altRv = SB_values[3] + 1
					symmetricalRatio = (refFw * altRv) / (altFw * refRv)
					refRatio = refRv / refFw
					altRatio = altFw / altRv
					SOR = round(math.log(symmetricalRatio) + refRatio - altRatio, 3)
					if SOR > max_SOR:
						if 'highSOR' not in alt_fails:
							alt_fails.append('high_SOR')
				except:
					print('no strand bias related information found in vcf')
			try:
				try:
					vaf_samples = format_values['VAF'].split(',')
					vaf_sample = vaf_samples[o]
					if float(vaf_sample) > float(dict_para['max_VAF_sample'].split("#")[0]):
						alt_fails.append('high_VAF_sample')

				except:
					vaf_samples = format_values['AF'].split(',')
					vaf_sample = vaf_samples[o]
					if float(vaf_sample) > float(dict_para['max_VAF_sample'].split("#")[0]):
						alt_fails.append('high_VAF_sample')
			except:
				print('no VAF sample found in formats field')

			try:
				vaf_pop_match = re.search(r"gnomad40_exome_AF=(.*?);", variant).group(1)
				if vaf_pop_match == '.':
					if dict_para['keep_variant_if_no_VAF_pop'].upper() == 'FALSE' or dict_para['keep_variant_if_no_VAF_pop'].upper() == 'NO':
						alt_fails.append('no_VAF_pop')
				elif float(vaf_pop_match) >= float(filter_af_population):
					alt_fails.append('high_VAF_pop')
			except:
				try:
					vaf_pop_raw_match1 = re.search(r"gnomad40_exome_AF_raw=(.*?);", variant).group(1)
					if vaf_pop_raw_match1 == '.':
						if dict_para['keep_variant_if_no_VAF_pop'].upper() == 'FALSE' or dict_para['keep_variant_if_no_VAF_pop'].upper() == 'NO':
							alt_fails.append('no_VAF_pop')
					elif float(vaf_pop_raw_match1) >= float(filter_af_population):
						alt_fails.append('high_VAF_pop')
				except:
					print('no gnomad40_exome_AF in variant')
			try:
				if dict_para['filter_on_mutation_location'].upper() == 'TRUE' or dict_para['filter_on_mutation_location'].upper() == 'YES':
					location = re.search(r"Func.refGene=(.*?);", variant).group(1)
					if location == '.':
						pass
					elif location not in mutation_to_save:
						alt_fails.append('mutation_location')
			except:
				print("can't filter on mutation location since VCF is not annotated by ANNOVAR or SnpEff")
				anno = False
			try:
				poly_score_match = re.search(r"Polyphen2_HDIV_score=(.*?);", variant).group(1)
				poly_pred_match = re.search(r"Polyphen2_HDIV_pred=(.*?);", variant).group(1)
				if poly_score_match == '.':
					if 'NO' in dict_para['keep_variant_if_no_Polyphen2_info'].upper() or 'FALSE' in dict_para['keep_variant_if_no_Polyphen2_info'].upper():
						alt_fails.append('no_Polyphen2_info')
				elif float(poly_score_match) < float(dict_para['min_PolyPhen2_score']):
					alt_fails.append('low_Polyphen2_score')
				if poly_pred_match == '.':
					if 'NO' in dict_para['keep_variant_if_no_Polyphen2_info'].upper() or 'FALSE' in dict_para['keep_variant_if_no_Polyphen2_info'].upper():
						alt_fails.append('no_Polyphen2_info')
				elif poly_pred_match not in dict_para['PolyPhen2_preds_to_keep'].upper().split(',') and 'ALL' not in dict_para['PolyPhen2_preds_to_keep'].upper():
					if 'NO' in dict_para['keep_variant_if_no_Polyphen2_info'].upper() or 'FALSE' in dict_para['keep_variant_if_no_Polyphen2_info'].upper():
						alt_fails.append('polyphen2_pred_not_allowed')
			except:
				print('no Poplyphen2 info in variant')
			try:
				sift_score_match = re.search(r"SIFT_score=(.*?);", variant).group(1)
				sift_pred_match = re.search(r"SIFT_pred=(.*?);", variant).group(1)
				if sift_score_match == '.':
					if 'NO' in dict_para['keep_variant_if_no_SIFT_info'].upper() or 'FALSE' in dict_para['keep_variant_if_no_SIFT_info'].upper():
						alt_fails.append('no_SIFT_info')
				elif 'NO' not in dict_para['filter_on_SIFT_score'].upper() and float(sift_score_match) > float(dict_para['max_SIFT_score']):
					alt_fails.append('high_SIFT_score')
				if sift_pred_match == '.':
					if 'NO' in dict_para['keep_variant_if_no_SIFT_info'].upper() or 'FALSE' in dict_para['keep_variant_if_no_SIFT_info'].upper():
						alt_fails.append('no_SIFT_info')
				elif sift_pred_match not in dict_para['SIFT_preds_to_keep'].upper().split(',') and 'ALL' not in dict_para['SIFT_preds_to_keep'].upper():
					if 'NO' in dict_para['keep_variant_if_no_SIFT_info'].upper() or 'FALSE' in dict_para['keep_variant_if_no_SIFT_info'].upper():
						alt_fails.append('sift_pred_not_allowed')
			except:
				print('no SIFT info in variant')

			fails = (str(alt_fails).replace('[',"").replace(']',"").replace("'","").replace('"',"").
					 replace(": ,",":").replace(" ",""))
			alt_v = alt.split(',')[o]
			if anno:
				variant_line = (chr + '\t' + pos + '\t' + id + '\t' + ref + '\t' + alt_v + '\t' + qual + '\t' + filter_pass + '\t' + record[7].split('ANNO')[0] + 'ANNO' + variant +
						   fails + '\t' + record[-2] + '\t' + record[-1] + '\n')
			dict_fails[variant_line] = alt_fails
			o += 1

		for key in dict_fails.keys():
			if dict_fails[key] == ['fails: ']:
				dict_fails[key] = ['PASS']
				return False, dict_fails
			else:
				dict_fails[key] = ''.join(dict_fails[key])
				return True, dict_fails

	else:
		fail = {'all': set()}
		for i in range(nb_alt):
			fail[i] = set()

		# Minimum allelic depth
		filter_alt_ad = dict_para['min_alt_AD']
		alt_ad = float(AD[1])  # allelic depth of alternative allele
		if float(alt_ad) < float(filter_alt_ad):
			fail[i].add('low_AD')

		try:
			max_vaf_sample = float(dict_para['max_VAF_sample'].split("#")[0])
		except:
			pass

		vaf_sample_ok = False
		for format in format_values.keys():
			if format == 'VAF':
				elmt = format_values[format]
				if ',' in elmt:
					elmt = elmt.split(',')
					elmt = float(elmt[1])
				else:
					elmt = float(elmt)
				if elmt > max_vaf_sample:
					fail['all'].add('high_VAF_sample')
				else:
					vaf_sample_ok = True

			if format == 'AF':
				elmt = format_values[format]
				if ',' in elmt:
					elmt = elmt.split(',')
					elmt = float(elmt[1])
				else:
					elmt = float(elmt)
				if elmt > float(max_vaf_sample):
					fail['all'].add('high_VAF_sample')
				else:
					vaf_sample_ok = True

		for id, elmt in info.items():
			if id == 'SIFT_score':
				if elmt == '.':
					if 'NO' in dict_para['keep_variant_if_no_SIFT_info'].upper() or 'FALSE' in dict_para['keep_variant_if_no_SIFT_info'].upper():
						fail['all'].add('no_SIFT_info')
				elif float(elmt) > float(dict_para['max_SIFT_score']):
					fail['all'].add('high_SIFT_score')

			if id == 'Polyphen2_HDIV_score':
				if elmt == '.':
					if 'NO' in dict_para['keep_variant_if_no_Polyphen2_info'].upper() or 'FALSE' in dict_para['keep_variant_if_no_Polyphen2_info'].upper():
						fail['all'].add('no_Polyphen2_info')
				elif 'NO' not in dict_para['filter_on_Polyphen2_score'].upper() and float(elmt) < float(dict_para['min_PolyPhen2_score']):
					fail['all'].add('low_Polyphen2_score')

			if id == 'SIFT_pred':
				if elmt == '.':
					if 'NO' in dict_para['keep_variant_if_no_SIFT_info'].upper() or 'FALSE' in dict_para['keep_variant_if_no_SIFT_info'].upper():
						fail['all'].add('no_SIFT_info')
				elif 'ALL' not in dict_para['SIFT_preds_to_keep']:
					allowed_sift_list = dict_para['SIFT_preds_to_keep'].upper().split(',')
					if elmt not in allowed_sift_list and 'YES' not in dict_para['keep_variant_if_no_SIFT_info'].upper() and 'TRUE' not in dict_para['keep_variant_if_no_SIFT_info'].upper():
						fail['all'].add('sift_pred_not_allowed')

			if id == 'Polyphen2_HVAR_pred':
				if elmt == '.':
					if 'NO' in dict_para['keep_variant_if_no_Polyphen2_info'].upper() or 'FALSE' in dict_para['keep_variant_if_no_Polyphen2_info'].upper():
						fail['all'].add('no_Polyphen2_info')
				if 'ALL' not in dict_para['PolyPhen2_preds_to_keep']:
					allowed_polyphen2_list = dict_para['PolyPhen2_preds_to_keep'].upper().split(',')
					if elmt not in allowed_polyphen2_list and 'YES' not in dict_para['keep_variant_if_no_Polyphen2_info'].upper() and 'TRUE' not in dict_para['keep_variant_if_no_Polyphen2_info'].upper():
						fail['all'].add('polyphen2_pred_not_allowed')

			if id == 'STRANDQ':
				if float(elmt) < float(dict_para['min_STRANDQ']):
					fail['all'].add('low_STRANDQ')

			if id == 'ExonicFunc.refGene':
				if elmt == '.' or elmt == 'unknown':
					if dict_para['remove_unknown_mutations'].upper() == 'YES' or dict_para['remove_unknown_mutations'].upper() == 'TRUE':
						fail['all'].add('unknown_mutation')
					elif 'unknown' in elmt:
						fail['all'].add('unknown_mutation')
				elif dict_para['remove_non_driver_mutations'].upper() == 'YES' or dict_para['remove_non_driver_mutations'].upper() == 'TRUE':
					if 'syno' in elmt and 'nons' not in elmt:
						fail['all'].add('non_driver_mutation')
					elif 'nonf' in elmt:
						fail['all'].add('non_driver_mutation')

			if id == 'DP':
				try:
					if float(elmt[0]) < float(filter_dp):
						fail['all'].add('low_DP')
				except:
					print('DP problem')

			if id == 'MQSBZ':
				if float(elmt) < -1*float(filter_mqsbz) or float(elmt) > float(filter_mqsbz):
					fail['all'].add('high_strand_bias')
				else:
					#print('MQSBZ OK')
					pass

			# median base quality by allele
			if id == 'MBQ':
				alt_mbq = re.search(r"MBQ=(.*?);", str(record)).group(1).split(',')[1]
				filter_alt_MBQ = dict_para['min_alt_MBQ']
				if float(alt_mbq) < float(filter_alt_MBQ):
					fail[i].add('low_MBQ')

			if id == 'POPAF':
				try:
					popaf1 = re.search(r"POPAF=(.*?);", str(record)).group(1)
					if float(popaf1) >= float(filter_af_population):
						fail[i].add('high_POPAF')
				except:
					print('not called with Mutect2 : no POPAF')

			gnomad_af_not_found = False
			# population allele frequencies
			if id == 'gnomad40_exome_AF':
				if elmt == '.':
					# if 'SIFT_score=.' in str(record):
					if dict_para['keep_variant_if_no_VAF_pop'].upper() == 'FALSE' or dict_para['keep_variant_if_no_VAF_pop'].upper() == 'NO':
						fail['all'].add('no_VAF_pop')
				elif 'gnomad40_exome_AF=.' in str(record):
					# if 'SIFT_score=.' in str(record):
					if dict_para['keep_variant_if_no_VAF_pop'].upper() == 'FALSE' or dict_para['keep_variant_if_no_VAF_pop'].upper() == 'NO':
						fail['all'].add('no_VAF_pop')
				else:
					af = float(elmt)
					gnomad_af_not_found = True
					try:
						if af >= float(filter_af_population):
							fail[i].add('high_VAF_pop')
					except:
						# if 'SIFT_score=.' in str(record):
						if dict_para['keep_variant_if_no_VAF_pop'].upper() == 'FALSE' or dict_para['keep_variant_if_no_VAF_pop'].upper() == 'NO':
							fail['all'].add('no_VAF_pop')

			# only allele with a potential functional effect (not silent)
			# todo : ADD FUNCOTATOR FIELDS ANALYSIS FOR stats.txt
			# todo : ADD FUNCOTATOR DRIVER MUTATIONS SELECTIONS
			# mutation_to_save = {'MISSENSE', 'NONSENSE', 'NONSTOP', 'START_CODON_SNP', 'DE_NOVO_START_IN_FRAME', 'DE_NOVO_START_OUT_FRAME',
			# 					'IN_FRAME_DEL', 'IN_FRAME_INS', 'FRAME_SHIFT_INS', 'FRAME_SHIFT_DEL', 'START_CODON_INS', 'START_CODON_DEL',
			# 					'DE_NOVO_START_IN_FRAME', 'DE_NOVO_START_OUT_FRAME', 'START_CODON_SNP', 'START_CODON_INS', 'START_CODON_DEL'}
			if id == 'FUNCOTATION':
				for i, funco in enumerate(elmt):
					t = funco.split('|')[5]  # mutation_type
					t2 = funco.split('|')[6]
					if not (t in mutation_to_save or (t == 'SPLICE_SITE' and t2 in mutation_to_save)):
						fail[i].add('not_functional')  # examples : "silent", "COULD_NOT_DETERMINE") "5'UTR" if WES (IGR = intergenic region)

			if id == 'Func.refGene':
				t = info['Func.refGene']
				if '\\' in t:
					t = t.split('\\')[0]
				if not (t in mutation_to_save):
					fail[i].add('mutation_location')
					# print('non exonic')

			if id == 'QUAL':
				try:
					if float(elmt) < float(filter_qual):
						fail['all'].add('low_QUAL')
				except:
					pass  # no QUAL because no ANNOVAR

		if all([True if v == set() else False for v in fail.values()]):
			return False, {}
		else:
			return fail, {}


def tumor_sample_header_logging(x, logger):
	'''
	Is the tumor sample ID exists in the header
        Input : a record (x) and logger 
	Output : False if tumor sample ID if does not exist or True otherwise
	'''
	# if x.startswith('##source=Mutect2'):
	# 	print(x)
	if x.startswith('##tumor_sample='):
		strSampleCode = x.split('=')[1].rstrip()
		logger.info(f'Sample id : {strSampleCode}')
		return strSampleCode
	return False


def get_values(dict_para, record : str, strFieldSplit : dict, intIndex : dict, InfoFieldsToProcess:set):
	'''
        Extract useful informations too filter a line of the VCF
        Input : the record and dictionnary needed to find and split line fields
        Output : information needed to filter the line (AD, AF, Info fields such as Funcotation...)
        '''

	# Extract data from the record to feed further functions
	# 1. the Info field -> from string or list to Dictionary
	recInfo = {}

	# if dict_para['vcf_annotation_method'].upper() == 'ANNOVAR':
	# 	str_record = str(record)
	#
	# 	# AC (Allele Count)
	# 	ac_start_index = str_record.find("AC=") + 3
	# 	ac_end_index = str_record.find(";", ac_start_index)
	# 	ac_value_string = str_record[ac_start_index:ac_end_index]
	# 	ac_values = list(map(int, ac_value_string.split(",")))
	#
	# 	# AN (Allele Number)
	# 	an_start_index = str_record.find("AN=") + 3
	# 	an_end_index = str_record.find(";", an_start_index)
	# 	an_value_string = str_record[an_start_index:an_end_index]
	# 	an_value_string = an_value_string.replace("'", "")
	# 	an_value = int(an_value_string)
	#
	# 	# AF (Allele Frequency)
	# 	ac_sum = sum(ac_values)
	# 	af_value = ac_sum / an_value
	#
	# 	record[7] = f'AF={af_value};{record[7]}'
	# 	recInfo['AF'] = af_value

	recInfo['QUAL'] = record[5]
	recInfoOrig = record[intIndex['Info']].split(strFieldSplit['Info1'])
	for r in recInfoOrig:
		# To avoid empty values:
		if strFieldSplit['Info2'] in r:
			info = r.split(strFieldSplit['Info2'])
			try:
				if info[0] in InfoFieldsToProcess:
					info[1] = info[1].split(strFieldSplit['Funcota'])
					if (info[0] == 'FUNCOTATION'):
						for k, v in enumerate(info[1]):
							info[1][k] = v.strip("[]")
				recInfo[info[0]] = info[1]
			except:
				recInfo[info[0]] = info[1]
		else:
			# Key without a value (eg: PON)
			recInfo[r] = ''
	else:
		recFormat = record[intIndex['Format']].split(strFieldSplit['Format1'])

	try:
		recValues = record[intIndex['Values']].split(strFieldSplit['Format1'])
		recAD = recValues[recFormat.index('AD')].split(strFieldSplit['Format2'])
	except:
		sys.exit("VCF file is not correctly formatted. Please check the format of the file (choose tabulations, not spaces.")

	# if dict_para['vcf_annotation_method'].upper() != 'ANNOVAR':
	# 	recAF = recValues[recFormat.index('AF')].split(strFieldSplit['Format2'])
		#bizarre pcq 'DP=6;VDB=5.074656e-02;AF1=1;AC1=2;DP4=0,0,3,3;MQ=60;FQ=-45' dans info mais ni AD ni AF en tout cas pas dans recValues

	# 3. get the filter field in case record doesn't pass FUNCOTATOR filter
	try:
		if record[intIndex['Filter']] == 'PASS':
			recOriginalFilter = False
		else:
			recOriginalFilter = record[intIndex['Filter']]
	except:
		if record[intIndex['Filter']] == '.':
			recOriginalFilter = ''
	return (recInfo, recFormat, recValues, recAD, recOriginalFilter)


def create_new_records(dict_para, record, strFieldSplit, intIndex, InfoFieldsToProcess, strCRLF, filter_ad=1, filter_mbq=20, filter_dp=10, filter_qual=50,
					   filter_mqsbz=0.5,filter_af_population=0.00001):
	'''
	Creation of new records containing new LOncoG filter or without variants that dont pass filters  
        Input : record and dictionnary needed to find and split line fields
        Output : boolean (is the record pass filters) and the two new records
	'''

	recInfo, recFormat, recValues, recAD, recOriginalFilter = get_values(dict_para, record, strFieldSplit, intIndex, InfoFieldsToProcess)

	# #####################
	# Treatment of failure
	if len(record[intIndex['Alt']].split(',')) != 1:
		pass #to debug and understand alt variants

	format_values = {}
	j = 0
	for format in recFormat:
		if not format in format_values.keys():
			format_values[format] = recValues[j]
		j += 1

	fail, dict_alts = fail_filters(record, dict_para, recInfo, format_values, len(record[intIndex['Alt']].split(',')), recAD , filter_qual, filter_mqsbz, filter_dp, filter_af_population) #False if the variant does not pass the filters and the name of the filter that does not pass otherwise

	# ####################
	# Rebuilt and storage of filtered data
	# Variable
	strLotusFilterCode = 'LOncoG_filter'
	strFilterPassKey = 'OTHER_FILTER='  # LOncoG filter
	strFilterPassOk = 'PASS'
	blnPass = False                                                 # set as default : the filter was not successful

	# Treatment
	if dict_alts == {}:
		if fail:                                                        # suppress variant that don't pass filters
			to_suppr = set()
			filters_failed = ''
			try:
				if 'all' in fail.keys() and fail['all'] != set():
					filters_failed += 'all'
			except:
				print('a')

			for i in range(len(fail)-1):
				if fail[i] != set():
					filters_failed += str(i)+':'+':'.join([f for f in fail[i]])+','
					to_suppr.add(i)
				elif fail['all'] != {''}:
					pass
				else:
					filters_failed += str(i)+':PASS'+','

			filters_failed += str(fail['all']).replace('{', '').replace('}', '').replace("'", '').replace(' ', '').replace(":all", ':')
			filters_failed = filters_failed.rstrip(',')

			# Add a new item (new filters) to Info field
			record[intIndex['Info']] = record[intIndex['Info']]+str(';')+strFilterPassKey+str(filters_failed)

			if filters_failed != '':
				record[intIndex['Filter']] = strLotusFilterCode

			try:
				if not recOriginalFilter:
					cleanRecord = deepcopy(record) #Using of deep_copy to avoid modification of the original record
					try:
						if line_without_fail(fail, cleanRecord, to_suppr, intIndex, strLotusFilterCode, strFilterPassKey, strFilterPassOk): # Modification of the current record to keep only variants that pass filter
							# blnPass = True
							pass
					except:
						blnPass = False
				else:
					cleanRecord = record

			except:
				if not recOriginalFilter and not 'DP' in filters_failed:
					cleanRecord = deepcopy(record) #Using of deep_copy to avoid modification of the original record
					if line_without_fail(fail, cleanRecord, to_suppr, intIndex, strLotusFilterCode, strFilterPassKey, strFilterPassOk): # Modification of the current record to keep only variants that pass filter
						# blnPass = True
						pass
				else:
					cleanRecord = record
		else:
			#Add a new item to Info field
			record[intIndex['Info']] = record[intIndex['Info']]+str(';')+strFilterPassKey+strFilterPassOk
			cleanRecord = record
			if not recOriginalFilter or recOriginalFilter == '.': #if there is no filter in the original vcf file (no filter in the FILTER field)
				blnPass = True

		# Rebuilt the new full record
		newRecord = strFieldSplit[1].join(record)+strCRLF #for filtered file
		newRecord2 = strFieldSplit[1].join(cleanRecord)+strCRLF #for passed file
		return blnPass, newRecord, newRecord2, dict_alts  #blnPass is True if the variant passes the filters,
		# newRecord is the record with the new filter and newRecord2 is the record without variants that don't pass filters
	else:
		return 'x', 'x', 'x', dict_alts


def merge_and_write_alt_variants(l_variants):
	if len(l_variants) == 1:
		return l_variants[0]
	else:
		common_start = l_variants[0].split('ANNOVAR_DATE')[0]
		alternative = common_start.split('\t')[4]
		alts = []
		common_end = l_variants[0].split('\t')[-2] + '\t' + l_variants[0].split('\t')[-1]
		mid_part = []
		for variant in l_variants:
			alt = variant.split('\t')[4]
			alts.append(alt)
			mid_part.append("ANNOVAR_DATE" + variant.split('ANNOVAR_DATE')[1].split('\t')[-3])
		mid_part = ';ALLELE_END;'.join(mid_part)
		common_start = common_start.replace(alternative, ",".join(alts))
		return common_start + mid_part + '\t' + common_end


def write_filter_stats(dict_para, dict_colors, write_filtered_file, chr_variants, passed_count):
	stats = {'mutation_location':0, 'no_VAF_pop':0, 'high_POPAF':0, 'high_VAF_pop':0, 'high_VAF_sample':0, 'high_SOR':0, 'low_STRANDQ':0, 'low_MBQ':0, 'low_AD':0, 'high_POPAF':0, \
			 'low_DP':0, 'low_QUAL':0, 'high_strand_bias':0, 'VAF_sample':0, 'low_Polyphen2_score':0, 'high_SIFT_score':0, 'polyphen2_pred_not_allowed':0, 'sift_pred_not_allowed':0,
			 'no_SIFT_info':0, 'no_Polyphen2_info':0, 'non_driver_mutation':0, 'unknown_mutation':0, 'no_SIFT_info':0, 'no_Polyphen2_info':0, 'not_functional':0}

	if dict_para['colored_execution'].upper() == 'FALSE' or dict_para['colored_execution'].upper() == 'NONE':
		no_tqdm_bar = True
		color = ''
	else:
		no_tqdm_bar = False
		color = dict_colors['yellow2']
	with tqdm(total=chr_variants, disable=no_tqdm_bar, bar_format='{l_bar}{bar:30}{r_bar}', ncols=150, smoothing=1) as pbar:
		pbar.set_description(color + f'-> reporting failed filters counts in stats.txt: ')  # progress bar
		with open(write_filtered_file, 'r') as file:
			for line in file:
				if not line.startswith('#'):
					fails = re.findall(r'fails:.*?(?=\t)', line)[0]
					for key in stats.keys():
						if key.upper() in fails.upper():
							stats[key] += 1
					pbar.update(1)
	keys_to_remove = [key for key, value in stats.items() if value == 0]
	for key in keys_to_remove:
		del stats[key]

	out_stats = dict_para['output_path_sample'].replace(".vcf", "") + 'filtered_stats.txt'
	if dict_para['verbose_prints'].upper() == 'TRUE':
		print('Write stats file...')
	with open(out_stats, 'w') as f:
		f.write('Number of variants: ' + str(chr_variants) + '\n')
		f.write('Number of variants that passed the filter: ' + str(passed_count) + '\n')
		f.write('Number of variants that did not pass the filters: ' + str(chr_variants - passed_count) + '\n')
		f.write('\n----- FILTER -----\n')
		f.write('number of variants that failed to pass each filter\n')
		for key in stats.keys():
			f.write(key + ': ' + str(stats[key]) + '\n')

def filter(dict_colors, dict_para, output1, output2, vcf_file: str, logger: str, working_method: str, filter_ad: int, filter_qual: int, filter_mbq: int, filter_dp: int,
		   filter_mqsbz: int, filter_af_population: float):
	'''
	Filters out variants:
	Input : vcf file
	Output : 2 vcf files : one containing variant with new filter column and another one without variants that don't pass filters
	logfile : the full path with file name for the log file to write in
	working_method : 2 possibilities : 'InMemory' (more speed but higher memory consumption) or 'Direct' (slow speed but low memory consumption)
	'''

	filter_mbq = int(filter_mbq)
	filter_dp = int(filter_dp)
	filter_qual = int(filter_qual)
	filter_mqsbz = float(filter_mqsbz)
	filter_af_population = float(filter_af_population)

	passed_count = 0

	if dict_para['verbose_prints'].upper() == 'TRUE':
		print(f'Read file :\t{vcf_file}\nWrite in {output1} and {output2}\n')

	# log data
	logger.info(f'Working method chosen by user : {str(working_method)}')
	logger.info(f'Read input file : {str(vcf_file)}')
	logger.info(f'Write to output files : {str(output1)} and {str(output2)}')

	# ######################################################################################
	# Opening the file, and read it line by line, store it in a list for further treatments

	mutations_file = dict_para['output_path_sample'].split('.')[0] + '/strong_mutations.txt'
	if os.path.exists(mutations_file):
		# Remove the existing mutations file if there is one
		os.remove(mutations_file)

	with open(vcf_file, mode='r', encoding='latin-1', errors='backslashreplace') as obFileInput:		 #read in latin-1 because some char from funcotator anotation are not in utf-8

		n = 0				# total number of lines in the VCF
		n_header = 0			# number of header lines from the VCF

		content = obFileInput.read()
		if 'chrX' not in content:
			sys.exit(f"ERROR: {vcf_file} looks cut (no X chromosome). Please check your input file.")

		if working_method.upper() == 'INMEMORY':
			lstInput1=[]	# The original data
			lstInput1 = obFileInput.readlines()
			n = len(lstInput1)
			n_header = len(list(takewhile(lambda x: x[0]=='#', lstInput1)))

		elif working_method.upper() == 'DIRECT':
			# read the open file in r mode for counting lines, then close and reopen it
			n = sum(1 for _ in obFileInput)
			obFileInput.seek(0,0)
			n_header = len(list(takewhile(lambda x: x[0]=='#', (line for line in obFileInput))))
			obFileInput.seek(0,0)

		if dict_para['verbose_prints'].upper() == 'TRUE':
			print(f'Number of lines to read : {n}\nNumber of header lines : {n_header}')
		logger.info(f'Input file contains : {str(n)} lines ({n_header} header and {n-n_header} data)')


		#####################################################
		# Opening for writing outputs
		filtered_file = open(output1, 'w', encoding='latin-1')  # filtered.vcf -> all variants with a new complementary filter
		passed_file = open(output2, 'w', encoding='latin-1')	# passed.vcf -> only variants that passed the filter

		# Criteria for header (only 6 (intlenghtCrit) char long for searching frame)
		intlenghtCrit = 6
		strHeaderFilter = '##FILT'				# To locate filter fields
		strHeaderInfo = '##INFO'				# To locate info fields

		# New lines to add in header
		strNewFilter = '##FILTER=<ID=LOncoG_filter,Description="Mutation does not pass LOncoG filters">'
		strNewInfo = '##INFO=<ID=OTHER_FILTER,Number=A,Type=String,Description="Other filter that dont pass.">'

		# Line ending
		strCRLF = '\n'

		# Counters
		intCountHeaderIns = 0					# Inserted header lines counter
		intCountHeaderMod = 0					# Modified header lines counter
		count = 0						# Total lines counter

		id_column = ''						# Saved the id of colunm for later

		# Presetting variable for reading
		sample_id = []
		if working_method.upper() == 'INMEMORY':
			strLinePrevious = lstInput1[0]
		elif working_method.upper() == 'DIRECT':
			strLinePrevious = '#'

		# searching through the header
		if working_method.upper() == 'INMEMORY':
			try:
				command = f"grep -c -v '^#' {vcf_file}"
				result = subprocess.run(command, shell=True, capture_output=True, text=True)
				count_variants = int(result.stdout.strip())
			except:
				count_variants = 0
				with open(vcf_file, "r") as file:
					for line in file:
						if not line.startswith("#"):
							count_variants += 1

			if dict_para['verbose_prints'].upper() == 'TRUE':
				print('Header treatment')
			id_column = lstInput1[n_header-1].lstrip('#')
			lstSave = deepcopy(lstInput1)
			dict_colors = list_colors()
			if dict_para['verbose_prints'].upper() == 'TRUE':
				print('Header treatment')

			if dict_para['colored_execution'].upper() == 'TRUE' or dict_para['colored_execution'].upper() == 'YES':
				blue_color = '\033[34m'
				no_tqdm_bar = False
			else:
				blue_color = ''
				no_tqdm_bar = True
			for i in tqdm(range(n_header), disable=no_tqdm_bar, bar_format=f"{blue_color}{{l_bar}}{{bar}}{blue_color}{{r_bar}}"):
				x = lstSave[i]
				# 1. Search for a specific information without affecting data integrity
				sample_id.append(tumor_sample_header_logging(x, logger)) # if ##tumor_sample in header, return the sample id
				# 2. Search for the end of a theme : to insert a line at the end of the theme
				strLinePreviousStart = strLinePrevious[0:intlenghtCrit]
				if x[0:intlenghtCrit] != strLinePreviousStart:
					if strLinePreviousStart == strHeaderFilter:
						lstInput1 = lstInput1[0:i+intCountHeaderIns]+[strNewFilter+strCRLF]+lstInput1[i+intCountHeaderIns:]
						intCountHeaderIns +=1
					elif strLinePreviousStart == strHeaderInfo:
						lstInput1 = lstInput1[0:i+intCountHeaderIns]+[strNewInfo+strCRLF]+lstInput1[i+intCountHeaderIns:]
						intCountHeaderIns +=1
				strLinePrevious = x
			del lstSave

		elif working_method.upper() == 'DIRECT':
			try:
				command = f"grep -c -v '^#' {vcf_file}"
				result = subprocess.run(command, shell=True, capture_output=True, text=True)
				count_variants = int(result.stdout.strip())
			except:
				count_variants = 0
				with open(vcf_file, "r") as file:
					for line in file:
						if not line.startswith("#"):
							count_variants += 1
			x = obFileInput.readline().strip()
			if dict_para['verbose_prints'].upper() == 'TRUE':
				print('Header treatment')
			
			if dict_para['colored_execution'].upper() == 'TRUE' or dict_para['colored_execution'].upper() == 'YES':
				color = dict_colors['yellow2']
				no_tqdm_bar = False
			else:
				color = ''
				no_tqdm_bar = True
			with tqdm(total=n_header, disable=no_tqdm_bar, bar_format='{l_bar}{bar:30}{r_bar}', ncols=100, smoothing=1) as pbar:
				pbar.set_description(color + f'-> parsing vcf header')
				while x[0] == '#':
					# 1 Write unmodified lines
					if strLinePrevious != '#':
						filtered_file.write(str(strLinePrevious+strCRLF)) # Simple copy of the actual value
						passed_file.write(str(strLinePrevious+strCRLF))
					# 2. Search for a specific information without affecting data integrity
					sample_id.append(tumor_sample_header_logging(x, logger)) # if ##tumor_sample in header, return the sample id containing this line
					# 3. Search for the end of a theme : to insert a line at the end of the theme
					strLinePreviousStart = strLinePrevious[0:intlenghtCrit]
					if x[0:intlenghtCrit] != strLinePreviousStart:
						if strLinePreviousStart == strHeaderFilter:
							filtered_file.write(str(strNewFilter+strCRLF))
							passed_file.write(str(strNewFilter+strCRLF))
							intCountHeaderIns += 1
						elif strLinePreviousStart == strHeaderInfo:
							filtered_file.write(str(strNewInfo+strCRLF))
							passed_file.write(str(strNewInfo+strCRLF))
							intCountHeaderIns +=1
						if x[0] == '#' and x[1] != '#': #if we get to #CHROM line, just before informations
							intDataStartPos=obFileInput.tell()
							id_column = x.lstrip('#')
							filtered_file.write(str(x+strCRLF))
							passed_file.write(str(x+strCRLF))
					strLinePrevious = x
					count += 1
					pbar.update(1)
					x = obFileInput.readline().strip()
			pbar.close()
		if not any(sample_id):
			logger.warning('No sample id !')
		logger.info(f'{intCountHeaderIns} line(s) added to header')
		logger.info(f'Actual number of header line is {n_header+intCountHeaderIns} ({n_header} before)')
		n_header = n_header+intCountHeaderIns


		# Criteria for data
		# To know how to split each sub-fields of datas

		strFieldSplit = {1: '\t', 'Info1': ';', 'Info2' : '=', 'Funcota' : ',', 'Format1' : ':', 'Format2' : ','}

		# index for fields
		intIndex = {'Chr' : id_column.split(strFieldSplit[1]).index('CHROM'), 'Alt' : id_column.split(strFieldSplit[1]).index('ALT'),
					'Filter' : id_column.split(strFieldSplit[1]).index('FILTER'), 'Info' : id_column.split(strFieldSplit[1]).index('INFO'),
					'Format' : id_column.split(strFieldSplit[1]).index('FORMAT')}
		intIndex['Values'] = intIndex['Format']+1

		InfoFieldsToProcess = {'DP', 'MFRL', 'MBQ', 'POPAF', 'FUNCOTATION'}

		# Counter for new files
		intCountData = 0
		intCountPass = 0

		# ###################################################
		# Analyze and extract genotype results, line by line

		if working_method.upper() == 'INMEMORY':
			lstInput2 = deepcopy(lstInput1)
			for id_line, x in enumerate(lstInput1[n_header:]):
				record = x_line.split(strFieldSplit[1])
				blnPass, newRecord, newRecord2, dict_alts = create_new_records(dict_para, record, strFieldSplit, intIndex,
																			   InfoFieldsToProcess, strCRLF, filter_ad, filter_mbq, filter_dp, filter_qual, filter_mqsbz,
																			   filter_af_population)
				if dict_alts == {}:
					count += 1
					if newRecord.startswith('chr'):
						newRecord = re.sub(r'OTHER_FILTER=', 'fails:', newRecord)
						newRecord = re.sub(":all", ":", newRecord)
						l_newRecord = newRecord.split('\t')
						l_newRecord[6] = 'LOncoG_filter'
						newRecord = '\t'.join(l_newRecord)
						if blnPass:
							newRecord = re.sub(r'fails:PASS', 'fails:NONE', newRecord)
							if newRecord.startswith('chr'):
								newRecord = re.sub(r'OTHER_FILTER=', 'fails:', newRecord)
								newRecord = re.sub(":all", ":", newRecord)
								l_newRecord = newRecord.split('\t')
								l_newRecord[6] = 'LOncoG_filter'
								newRecord = '\t'.join(l_newRecord)
								passed_count += 1
								passed_file.write(newRecord)
						filtered_file.write(newRecord)
				else:
					count += len(dict_alts.keys())
					l_multiple_passed = []
					for key in dict_alts:
						l_key = key.split('\t')
						l_key[6] = 'LOncoG_filter'
						new_key = '\t'.join(l_key)
						if dict_alts[key] == ['PASS'] or dict_alts[key] == ['fails: ']:
							new_key = re.sub(r';fails:', ';fails:NONE', key)
							variant_without_fail = re.sub(r';fails:.*?\t', '\t', key)
							l_variant_without_fail = variant_without_fail.split('\t')
							l_variant_without_fail[6] = 'LOncoG_filter'
							variant_without_fail = '\t'.join(l_variant_without_fail)
							if variant_without_fail.startswith('chr'):
								passed_count += 1
								l_multiple_passed.append(variant_without_fail)
						if 'LOncoG_filter' not in new_key:
							l_new_key = new_key.split('\t')
							l_new_key[6] = 'LOncoG_filter'
							new_key = '\t'.join(l_new_key)
						if new_key.startswith('chr'):
							filtered_file.write(new_key)
					try:
						merged_line = merge_and_write_alt_variants(l_multiple_passed)
						if merged_line.startswith('chr'):
							passed_file.write(merged_line)
					except:
						pass
				x_line = obFileInput.readline().strip()
				count += 1
				pbar.update(1)

		elif working_method.upper() == 'DIRECT':
			obFileInput.seek(intDataStartPos)
			x_line = obFileInput.readline().strip()
			count = 1
			if dict_para['verbose_prints'].upper() == 'TRUE':
				print(f'\nSaving {str(output1)} and {str(output2)} on disk', flush=True)

			if dict_para['verbose_prints'].upper() == 'TRUE':
				print('Header treatment')

			if dict_para['colored_execution'].upper() == 'TRUE' or dict_para['colored_execution'].upper() == 'YES':
				color = dict_colors['yellow2']
				no_tqdm_bar = False
			else:
				color = ''
				no_tqdm_bar = True

			chr_variants = 0

			with tqdm(total=count_variants, disable=no_tqdm_bar, bar_format='{l_bar}{bar:30}{r_bar}', ncols=140, smoothing=1) as pbar:
				pbar.set_description(color + f'-> parsing vcf content (filtering variants)')
				while x_line != '':
					record = x_line.split(strFieldSplit[1])
					blnPass, newRecord, newRecord2, dict_alts = create_new_records(dict_para, record, strFieldSplit, intIndex,
															InfoFieldsToProcess, strCRLF, filter_ad, filter_mbq, filter_dp, filter_qual, filter_mqsbz, filter_af_population)
					if dict_alts == {}:
						count += 1
						if newRecord.startswith('chr'):
							newRecord = re.sub(r'OTHER_FILTER=', 'fails:', newRecord)
							newRecord = re.sub(":all",":", newRecord)
							l_newRecord = newRecord.split('\t')
							l_newRecord[6] = 'LOncoG_filter'
							newRecord = '\t'.join(l_newRecord)
							if blnPass:
								newRecord = re.sub(r'fails:PASS', 'fails:NONE', newRecord)
								if newRecord.startswith('chr'):
									newRecord = re.sub(r'OTHER_FILTER=', 'fails:', newRecord)
									newRecord = re.sub(":all", ":", newRecord)
									l_newRecord = newRecord.split('\t')
									l_newRecord[6] = 'LOncoG_filter'
									newRecord = '\t'.join(l_newRecord)
									passed_count += 1
									passed_file.write(newRecord)
							chr_variants += 1
							filtered_file.write(newRecord)
					else:
						count += len(dict_alts.keys())
						l_multiple_passed = []
						for key in dict_alts:
							l_key = key.split('\t')
							l_key[6] = 'LOncoG_filter'
							new_key = '\t'.join(l_key)
							if dict_alts[key] == ['PASS'] or dict_alts[key] == ['fails: ']:
								new_key = re.sub(r';fails:', ';fails:NONE', key)
								variant_without_fail = re.sub(r';fails:.*?\t', '\t', key)
								l_variant_without_fail = variant_without_fail.split('\t')
								l_variant_without_fail[6] = 'LOncoG_filter'
								variant_without_fail = '\t'.join(l_variant_without_fail)
								if variant_without_fail.startswith('chr'):
									passed_count += 1
									l_multiple_passed.append(variant_without_fail)
							if 'LOncoG_filter' not in new_key:
								l_new_key = new_key.split('\t')
								l_new_key[6] = 'LOncoG_filter'
								new_key = '\t'.join(l_new_key)
							if new_key.startswith('chr'):
								chr_variants += 1
								filtered_file.write(new_key)
						try:
							merged_line = merge_and_write_alt_variants(l_multiple_passed)
							if merged_line.startswith('chr'):
								passed_file.write(merged_line)
						except:
							pass
					x_line = obFileInput.readline().strip()
					pbar.update(1)
				pbar.close()

			if passed_count == 0:
				word = 'variant'
			else:
				word = 'variants'
			if dict_para['colored_execution'].upper() == 'TRUE' or dict_para['colored_execution'].upper() == 'YES':
					print(f'\n\033[1m{color}{passed_count}/{count-2}\033[0m{color} {word} have passed the filter!')
			else:
				print(f'\n{passed_count}/{count-1} {word} have passed the filter!')
		if working_method == 'InMemory':
			print(f'Saving {str(output2)} and {str(output2)} on disk', flush=True)
			for count, record in enumerate(tqdm(lstInput1)):
				if record.starswith('chr'):
					filtered_file.write(record)
				if lstInput2[count]:
					passed_file.write(lstInput2[count])
				count += 1

		filtered_file.close()

		filtered_file_path = dict_para['output_path_sample'].split('.vcf')[0] + '/' + dict_para['output_path_sample'].split('samples/')[1].replace("/","").replace(".vcf","_filtered.vcf")
		write_filter_stats(dict_para, dict_colors, filtered_file_path, chr_variants, passed_count-1)

		if dict_para['verbose_prints'].upper() == 'TRUE':
			print(f'Number dof lines containing variants : {str(intCountData)}')
			print(f'Number of lines flagged as PASS : {str(intCountPass)}')

		logger.info(f'Number of lines for input header :\t{str(n_header)}')
		logger.info(f'Number of lines for output header :\t{str(n_header)} of which {str(n_header-intCountHeaderIns)} lines added and {str(intCountHeaderIns)} line modified')
		logger.info(f'Number of variants flagged as PASS :\t{str(intCountPass)}')

		'''# ##############################################
		# Close the 2 destination files
		filtered_file.close()
		passed_file.close()
		del filtered_file
		del passed_file
		if working_method == 'InMemory':
			del lstInput1
			del lstInput2

		logger.info(f' -> file {str(output1)} :\t{str(n_header+intCountData)} total lines of which {str(intCountData)} genotype data lines')
		logger.info(f' -> file {str(output2)} :\t{str(n_header+intCountPass)} total lines of which {str(intCountPass)} genotype data lines')
		logger.info(f'* End of filtering *')
		logger.info('**************************************************************************************************************')'''


def main(dict_colors, args, vcf_file, recover):

	if recover:
		config_path = args['config_path']
		args = manage_parameters(config_path, True)
	dict_para = copy.copy(args)

	output_path = dict_para['output_path_sample'].replace(".vcf","")
	output_path = output_path + output_path.split('samples/')[1].replace('/', '_')
	output1 = output_path + 'filtered.vcf'
	output2 = output_path + 'passed.vcf'
	# print(output1, output2)

	filter_mbq = args['min_alt_MBQ']
	filter_dp = args['min_DP']
	filter_mqsbz = args['max_SB']
	filter_qual = args['min_QUAL']
	filter_ad = args['min_alt_AD']
	filter_af_sample = args['max_VAF_sample']
	filter_af_population = args['max_VAF_pop']

	working_method = args['working_method']	# working_method = 'InMemory' (more speed but higher memory consumption) or 'Direct' (slow speed but low memory consumption)
	working_directory = Path(vcf_file).parent.absolute()

	# Logger configuration
	logger = logging.getLogger('LOncoG filter')
	logger.setLevel(logging.DEBUG)
	fh = logging.FileHandler(args['log'])
	fh.setLevel(logging.DEBUG)
	formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
	fh.setFormatter(formatter)
	logger.addHandler(fh)

	# Verification of given arguments
	try:
		logger.info('Verification of vcf file')
		verif_input_vcf(vcf_file)
		logger.info('- Vcf file ok -')
	except ValueError:
		print(f'Problem with vcf file {vcf_file}:', sys.exc_info())
		logger.error('- Problem with vcf file -')
		exit(1)
	try:
		logger.info('Verification of output files')
		verif_output(output1)
		logger.info('- Output file filtered ok -')
		verif_output(output2)
		logger.info('- Output file passed ok -')
	except ValueError:
		print(f'Problem with output file {output1}:', sys.exc_info())
		logger.error('- Problem with output file filtered -')
		print(f'Problem with output file {output2}:', sys.exc_info())
		logger.error('- Problem with output file passed -')
		exit(1)

	# Start
	logger.info('**************************************************************************************************************')
	logger.info('*** LOncoG filtering module ***')
	no_argument = ''
	logger.info(f'** cmd line : python loncog.py filter -v {str(vcf_file)} -o {output1} -wm {working_method} --DP {filter_dp} --MBQ {filter_mbq} '+str(no_argument)+' **')
	logger.info('* Start filtering *')
	logger.info(f'Working directory (vcf files folder) : {working_directory}')
	logger.info(f'Current directory : {Path().absolute()}')

	filter(dict_colors, dict_para, output1, output2, vcf_file, logger, working_method, filter_ad, filter_qual, filter_mbq, filter_dp, filter_mqsbz, filter_af_population)
	# End
