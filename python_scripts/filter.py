#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
import time
import logging
import os
import sys
from copy import deepcopy
from itertools import takewhile
from tqdm import tqdm
from python_scripts.manage_parameters import list_colors
from python_scripts.check_files import verif_input_vcf, verif_output
from python_scripts.path_modification import true_stem

def line_without_fail(record, to_suppr : str, intIndex : dict, strLotusFilterCode : str, strFilterPassKey : str, strFilterPassOk : str):
	'''
	Takes a record and deletes all variants that do not pass the filters, if no variants pass the filters returns False otherwise returns True
	Input : record and list of variants to suppressed
	Output : True (Record modified) or False (No record left)
	'''

	# Get the number of alternative variants
	nb_alt = len(record[intIndex['Alt']].split(','))

	# Suppression of variants that don't pass LOTUS filters
	alternative=[alt for i, alt in enumerate(record[intIndex['Alt']].split(',')) if i not in to_suppr]
	if alternative == []:
		return False
	else:
		record[intIndex['Alt']] = ','.join(alternative)

	# Filter modification
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
			values[i] = '/'.join( [char for i, char in enumerate(v.split(',')) if i not in to_suppr] )
		elif len(v.split(',')) == nb_alt+1:
			values[i] = '/'.join( [char for i, char in enumerate(v.split(',')) if i-1 not in to_suppr] )
	record[intIndex['Values']] = ':'.join(values)

	return True
	

def fail_filters(dict_para, info: {}, nb_alt: int, AD : [], filter_ad: int = 1, filter_qual=50, filter_mbq: int = 20, filter_mqsbz = 0.5, filter_dp: int = 10, filter_popaf: float = 0.0001, unpaired: bool = False):
	'''
	Do variant(s) at a position passes the filters ?
	Input : info, AD and AF fields and the number of variants
	Output : False (don't fail) or a dictionnary containing the failed filters
	'''

	#valeurs de filtres redéfinies ici si oubliées --> éviter cela!
	if dict_para['vcf_annotation_method'].upper() == 'ANNOVAR':
		if dict_para['sequencing'].upper() == 'WGS':
			mutation_to_save = set({'exonic', 'splicing', 'ncRNA', 'ncRNA_intronic', 'ncRNA_exonic', 'UTR5', 'UTR3', 'intronic', 'upstream', 'downstream', 'intergenic'})
		else :
			mutation_to_save = set({'exonic', 'splicing'})

	if dict_para['vcf_annotation_method'].upper() == 'FUNCOTATOR':
		if dict_para['sequencing'].upper() == 'WGS':
			mutation_to_save = set({'IGR', 'THREE_PRIME_UTR', 'FIVE_PRIME_UTR', 'INTRON', 'FIVE_PRIME_FLANK', 'MISSENSE', 'NONSENSE',
									'NONSTOP', 'RNA', 'LINCRNA', 'START_CODON_SNP', 'DE_NOVO_START_IN_FRAME', 'DE_NOVO_START_OUT_FRAME', 'IN_FRAME_DEL', 'IN_FRAME_INS',
									'FRAME_SHIFT_INS', 'FRAME_SHIFT_DEL', 'START_CODON_INS', 'START_CODON_DEL', 'DE_NOVO_START_IN_FRAME', 'DE_NOVO_START_OUT_FRAME'})
		elif dict_para['sequencing'].upper() == 'WES':
			mutation_to_save = set({'MISSENSE', 'NONSENSE', 'NONSTOP', 'RNA', 'LINCRNA', 'START_CODON_SNP', 'DE_NOVO_START_IN_FRAME',
									'DE_NOVO_START_OUT_FRAME', 'IN_FRAME_DEL', 'IN_FRAME_INS', 'FRAME_SHIFT_INS', 'FRAME_SHIFT_DEL', 'START_CODON_INS', 'START_CODON_DEL',
									'DE_NOVO_START_IN_FRAME', 'DE_NOVO_START_OUT_FRAME'})
	splice_site = 'SPLICE_SITE'

	fail = {'all': set()}
	for i in range(nb_alt):
		fail[i] = set()
	# if dict_para['vcf_annotation_method'].upper() == 'FUNCOTATOR':
	# 	if type(AF) != type([]):
	# 		AF = [AF]

	# Minimum allele sequencing depth
	for i, ad in enumerate(AD[1:]):
		if float(ad[0]) < float(filter_ad):
			fail[i].add('AD')

	# Minimum allele frequency
	# if dict_para['vcf_annotation_method'].upper() == 'FUNCOTATOR':
	# 	for i, af in enumerate(AF):
	# 		if float(af) < float(filter_af):
	# 			fail[i].add('AF')

	for id, elmt in info.items():
		if id == 'DP':
			if dict_para['vcf_annotation_method'].upper() == 'FUNCOTATOR':
				elmt = float(elmt[0])
			elif dict_para['vcf_annotation_method'].upper() == 'ANNOVAR':
				elmt = float(elmt)
			if elmt < int(filter_dp):
				fail['all'].add('DP')
			else:
				#print('DP OK')
				pass

		if id == 'MQSBZ':
			if dict_para['vcf_annotation_method'].upper() == 'FUNCOTATOR':
				print('check strand bias filter')
			elif dict_para['vcf_annotation_method'].upper() == 'ANNOVAR':
				elmt = float(elmt)
			if elmt < -1*float(filter_mqsbz) or elmt > float(filter_mqsbz):
				fail['all'].add('MQSBZ')
			else:
				#print('MQSBZ OK')
				pass

		# median fragment length by allele
		if not unpaired:
			if id == 'MFRL':
				for i, mfrl in enumerate(elmt[1:]):
					if float(mfrl) == 0:
						fail[i].add('MFRL')
		# median base quality by allele
		if id == 'MBQ':
			for i, mbq in enumerate(elmt[1:]):
				if float(mbq) < int(filter_mbq):
					fail[i].add('MBQ')
		# population allele frequencies
		if id == 'POPAF':
			for i, popaf in enumerate(elmt):
				popaf = 10**(-float(popaf)) #make it a probability by taking the negative log with base 10
				if float(popaf) >= float(filter_popaf):
					fail[i].add('POPAF')
		# only allele with a potential functional effect (not silent)
		if dict_para['vcf_annotation_method'].upper() == 'FUNCOTATOR':
			if id == 'FUNCOTATION':
				for i, funco in enumerate(elmt):
					t = funco.split('|')[5] #mutation_type
					t2 = funco.split('|')[6]
					if not (t in mutation_to_save or (t == 'SPLICE_SITE' and t2 in mutation_to_save)):
						fail[i].add('NOT_FUNCTIONAL') # examples : "silent", "COULD_NOT_DETERMINE") "5'UTR" if WES (IGR = intergenic region)
		elif dict_para['vcf_annotation_method'].upper() == 'ANNOVAR':
			if id == 'Func.refGene':
				t = info['Func.refGene']
				if '\\' in t:
					t = t.split('\\')[0]
				if not (t in mutation_to_save):
					fail[i].add('NOT_FUNCTIONAL')
					#print('non exonic')
				elif t == 'exonic':
					#print(t)
					pass
			if id == 'QUAL':
				qual = elmt
				if float(qual) < float(filter_qual):
					fail['all'].add('QUAL')

	if all([True if v == set() else False for v in fail.values()]):
		return(False)
	else:
		return(fail)


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


def get_values (dict_para, record : str, strFieldSplit : dict, intIndex : dict, InfoFieldsToProcess : set):
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
			if dict_para['vcf_annotation_method'].upper() == 'FUNCOTATOR':
			# Special case in Info field that needs treatment
				if info[0] in InfoFieldsToProcess:
					info[1] = info[1].split(strFieldSplit['Funcota'])
					if (info[0] == 'FUNCOTATION'):
						for k, v in enumerate(info[1]):
							info[1][k] = v.strip("[]")
				recInfo[info[0]]=info[1]
			elif dict_para['vcf_annotation_method'].upper() == 'ANNOVAR':
				recInfo[info[0]]=info[1]
		else:
			# Key without a value (eg: PON)
			recInfo[r]=''

	# 2. the AD and AF fields -> into List of float
	# if dict_para['vcf_annotation_method'].upper() == 'ANNOVAR':
	# 	try:
	# 		recFormat = (record[intIndex['Format']].split(strFieldSplit['Format1']))
	# 	except:
	# 		print('a')
	# 	recFormat.append('AF')
	# 	recAF = [str(af_value)]

	else:
		recFormat = record[intIndex['Format']].split(strFieldSplit['Format1'])

	recValues = record[intIndex['Values']].split(strFieldSplit['Format1'])
	recAD = recValues[recFormat.index('AD')].split(strFieldSplit['Format2'])

	# if dict_para['vcf_annotation_method'].upper() != 'ANNOVAR':
	# 	recAF = recValues[recFormat.index('AF')].split(strFieldSplit['Format2'])
		#bizarre pcq 'DP=6;VDB=5.074656e-02;AF1=1;AC1=2;DP4=0,0,3,3;MQ=60;FQ=-45' dans info mais ni AD ni AF en tout cas pas dans recValues

	# 3. get the filter field in case record doesn't pass FUNCOTATOR filter
	if dict_para['vcf_annotation_method'].upper() == 'FUNCOTATOR':
		if record[intIndex['Filter']] == 'PASS':
			recOriginalFilter = False
		else :
			recOriginalFilter = record[intIndex['Filter']]
	elif dict_para['vcf_annotation_method'].upper() == 'ANNOVAR':
		if record[intIndex['Filter']] == '.':
			recOriginalFilter = ''
	return (recInfo, recFormat, recValues, recAD, recOriginalFilter)


def create_new_records(dict_para, record, strFieldSplit, intIndex, InfoFieldsToProcess, strCRLF, filter_ad=1, filter_mbq=20, filter_dp=10, filter_qual=50, filter_mqsbz=0.5,filter_popaf=0.00001, unpaired: bool = False):
	'''
	Creation of new records containing new Lotus filter or without variants that dont pass filters  
        Input : record and dictionnary needed to find and split line fields
        Output : boolean (is the record pass filters) and the two new records
	'''

	recInfo, recFormat, recValues, recAD, recOriginalFilter = get_values(dict_para, record, strFieldSplit, intIndex, InfoFieldsToProcess)

	# #####################
	# Treatment of failure
	if len(record[intIndex['Alt']].split(',')) != 1:
		pass #to debug and understand alt variants

	fail = fail_filters(dict_para, recInfo, len(record[intIndex['Alt']].split(',')), recAD, filter_ad, filter_qual, filter_mbq, filter_mqsbz, filter_dp, filter_popaf, unpaired) #False if the variant does not pass the filters and the name of the filter that does not pass otherwise

	# ####################
	# Rebuilt and storage of filtered data

	# Variable
	strLotusFilterCode = 'LOTUS_filter'
	strFilterPassKey = 'OTHER_FILTER=' # lotus filter
	strFilterPassOk = 'PASS'
	blnPass = False                                                 # set as default : the filter was not successful

	#Treatment
	if fail:                                                        # suppress variant that don't pass filters
		to_suppr = set()
		filters_failed = ''
		if 'all' in fail.keys() and fail['all'] != set():
			filters_failed += 'all'
		for i in range(len(fail)-1):
			if fail[i] != set():
				filters_failed += str(i)+':'+':'.join([f for f in fail[i]])+','
				to_suppr.add(i)
			else:
				filters_failed += str(i)+':PASS'+','
		filters_failed = filters_failed.rstrip(',')
		#print(filters_failed)

		# Add a new item (new filters) to Info field
		record[intIndex['Info']] = record[intIndex['Info']]+str(';')+strFilterPassKey+str(filters_failed)

		if filters_failed != '':
			if not recOriginalFilter: #if there is no filter in the original vcf file (no filter in the FILTER field)
				# Add the 'Lotus' filter item in Filter field
				record[intIndex['Filter']] = strLotusFilterCode
			else:
				record[intIndex['Filter']] = record[intIndex['Filter']]+';'+strLotusFilterCode

		if dict_para['vcf_annotation_method'].upper() == 'ANNOVAR':
			if not recOriginalFilter:
				cleanRecord = deepcopy(record) #Using of deep_copy to avoid modification of the original record
				if dict_para['vcf_annotation_method'].upper() == 'FUNCOTATOR':
					if line_without_fail(cleanRecord, to_suppr, intIndex, strLotusFilterCode, strFilterPassKey, strFilterPassOk): # Modification of the current record to keep only variants that pass filter
						blnPass = True
				elif dict_para['vcf_annotation_method'].upper() == 'ANNOVAR':
					blnPass = False
			else:
				cleanRecord = record

		elif dict_para['vcf_annotation_method'].upper() == 'FUNCOTATOR':
			if not recOriginalFilter and not 'DP' in filters_failed:
				cleanRecord = deepcopy(record) #Using of deep_copy to avoid modification of the original record
				if line_without_fail(cleanRecord, to_suppr, intIndex, strLotusFilterCode, strFilterPassKey, strFilterPassOk): # Modification of the current record to keep only variants that pass filter
					blnPass = True
			else:
				cleanRecord = record
	else:
		#Add a new item to Info field
		record[intIndex['Info']] = record[intIndex['Info']]+str(';')+strFilterPassKey+strFilterPassOk
		cleanRecord = record
		if not recOriginalFilter: #if there is no filter in the original vcf file (no filter in the FILTER field)
			blnPass = True

	# Rebuilt the new full record
	newRecord = strFieldSplit[1].join(record)+strCRLF #for filtered file
	newRecord2 = strFieldSplit[1].join(cleanRecord)+strCRLF #for passed file
	return blnPass, newRecord, newRecord2 #blnPass is True if the variant passes the filters,
	# newRecord is the record with the new filter and newRecord2 is the record without variants that don't pass filters


def filter(dict_colors, dict_para, output1, output2, vcf_file: str, logger: str, working_method : str, filter_ad: int, filter_qual: int, filter_mbq: int, filter_dp: int,
			filter_mqsbz: int, filter_popaf: float, unpaired: bool):
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
	filter_popaf = float(filter_popaf)

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


		# ####################################################
		# Opening for writing outputs
		filtered_file = open(output1,'w', encoding='latin-1')		# filtered.vcf -> all variants with a new complementary filter
		passed_file = open(output2,'w', encoding='latin-1')		# passed.vcf -> only variants that passed the filter

		# Criteria for header (only 6 (intlenghtCrit) char long for searching frame)
		intlenghtCrit = 6
		strHeaderFilter = '##FILT'				# To locate filter fields
		strHeaderInfo = '##INFO'				# To locate info fields

		# New lines to add in header
		strNewFilter = '##FILTER=<ID=LOTUS_filter,Description="Mutation does not pass LOTUS filters">'
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
			if dict_para['verbose_prints'].upper() == 'TRUE':
				print('Header treatment')
			id_column = lstInput1[n_header-1].lstrip('#')
			lstSave = deepcopy(lstInput1)
			dict_colors = list_colors()
			blue_color = '\033[34m'
			if dict_para['colored_execution'].upper() == 'FALSE' or dict_para['colored_execution'].upper() == 'NONE':
				blue_color = ''
			for i in tqdm(range(n_header), bar_format=f"{blue_color}{{l_bar}}{{bar}}{blue_color}{{r_bar}}"):
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
			x = obFileInput.readline().strip()
			if dict_para['verbose_prints'].upper() == 'TRUE':
				print('Header treatment')

			color = dict_colors['yellow2']
			if dict_para['colored_execution'].upper() == 'FALSE' or dict_para['colored_execution'].upper() == 'NONE':
				color = ''
			with tqdm(total=n_header, bar_format='{l_bar}{bar:30}{r_bar}', ncols=100, smoothing=1) as pbar:
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
							intCountHeaderIns +=1
						elif strLinePreviousStart == strHeaderInfo :
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

		strFieldSplit = {1 : '\t', 'Info1' : ';', 'Info2' : '=', 'Funcota' : ',', 'Format1' : ':', 'Format2' : ','}

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
				record = x.strip().split(strFieldSplit[1])

				blnPass, newRecord, newRecord2 = create_new_records(dict_para, record, strFieldSplit, intIndex,
																		InfoFieldsToProcess, strCRLF, filter_ad, filter_mbq, filter_dp, filter_qual, filter_mqsbz, filter_popaf, unpaired)

				lstInput1[n_header+id_line]=str(newRecord)
				intCountData +=1
				if (blnPass):
					lstInput2[n_header+id_line]=str(newRecord2)
					blnPass = False		# reset to default
					intCountPass +=1
				else:
					lstInput2[n_header+id_line]=None

		elif working_method.upper() == 'DIRECT':
			obFileInput.seek(intDataStartPos)
			x_line = obFileInput.readline().strip()
			count = 1
			if dict_para['verbose_prints'].upper() == 'TRUE':
				print(f'\nSaving {str(output1)} and {str(output2)} on disk', flush=True)

			color = dict_colors['yellow2']
			with tqdm(total=n - (n_header - intCountHeaderIns), bar_format='{l_bar}{bar:30}{r_bar}', ncols=140, smoothing=1) as pbar:
				pbar.set_description(color + f'-> parsing vcf content (filtering variants)')
				while x_line != '':
					record = x_line.split(strFieldSplit[1])
					blnPass, newRecord, newRecord2 = create_new_records(dict_para, record, strFieldSplit, intIndex,
																		InfoFieldsToProcess, strCRLF, filter_ad, filter_mbq, filter_dp, filter_qual, filter_mqsbz, filter_popaf, unpaired)
					filtered_file.write(newRecord)
					intCountData += 1
					if blnPass:
						passed_count+=1
						passed_file.write(newRecord2)
						txt_file_path = dict_para['output_path_sample'].split('.')[0] + '/'
						if dict_para['discard_weak_variants'].upper() == 'FALSE' or 'NO' in dict_para['discard_weak_variants'].upper():
							txt_file = open(txt_file_path + "strong_mutations.txt", "a")
							variant = newRecord2.split('\t')[0:5]
							variant[4] = variant[4].replace(',','/')
							variant = str(variant).replace("'", "").replace("[", "").replace("]", "").replace('"', '').replace(" ", "")
							#print(variant)
							txt_file.write(variant + '\n')
						blnPass = False  # reset to default
					elif dict_para['discard_weak_variants'].upper() == 'FALSE' or dict_para['discard_weak_variants'].upper() == 'NO':
						#print('NOT PASSED')
						if newRecord.startswith("chr"):
							passed_file.write(newRecord)
					x_line = obFileInput.readline().strip()
					count += 1
					pbar.update(1)
				pbar.close()
			if dict_para['colored_execution'].upper() == 'TRUE' or dict_para['colored_execution'].upper() == 'YES':
					print(f'\033[1m{passed_count}\033[0m{color} variants have passed the filter!')
			else:
				print(f'{passed_count} variants have passed the filter!')
		if working_method == 'InMemory':
			print(f'Saving {str(output2)} and {str(output2)} on disk', flush=True)
			for count, record in enumerate(tqdm(lstInput1)):
				filtered_file.write(record)
				if lstInput2[count]:
					passed_file.write(lstInput2[count])
				count += 1

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


def main(dict_colors, dict_para, vcf_file, args, output_path):

<<<<<<< HEAD
<<<<<<< HEAD
	output1 = output_path + 'filtered.vcf'
	output2 = output_path + 'passed.vcf'
=======
	output1 = output_path + args['output_filtered']
	output2 = output_path + args['output_passed']
>>>>>>> 9c8a90e37f5f0a6717b01efd2a2d4ef069c38abb
=======
	output1 = output_path + args['output_filtered']
	output2 = output_path + args['output_passed']
>>>>>>> 9c8a90e37f5f0a6717b01efd2a2d4ef069c38abb
	filter_mbq = args['MBQ']
	filter_dp = args['DP']
	filter_mqsbz = args['MQSBZ']
	if dict_para['vcf_annotation_method'].upper() == 'ANNOVAR':
		filter_qual = args['QUALITY']
	filter_ad = args['AD']
	# filter_af = args['AF']
	filter_popaf = args['POPAF']
	unpaired = args['unpaired']
	if unpaired.upper() == 'TRUE':
		unpaired = True
	else :
		unpaired = False

	working_method = args['working-method']	# working_method = 'InMemory' (more speed but higher memory consumption) or 'Direct' (slow speed but low memory consumption)
	working_directory = Path(vcf_file).parent.absolute()

	# Logger configuration
	logger = logging.getLogger('LOTUS filter')
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
		print (f'Problem with vcf file {vcf_file}:', sys.exc_info())
		logger.error('- Problem with vcf file -')
		exit(1)
	try:
		logger.info('Verification of output files')
		verif_output(output1)
		logger.info('- Output file filtered ok -')
		verif_output(output2)
		logger.info('- Output file passed ok -')
	except ValueError:
		print (f'Problem with output file {output1}:', sys.exc_info())
		logger.error('- Problem with output file filtered -')
		print(f'Problem with output file {output2}:', sys.exc_info())
		logger.error('- Problem with output file passed -')
		exit(1)

	# Start
	logger.info('**************************************************************************************************************')
	logger.info('*** LOTUS filtering module ***')
	no_argument=''
	if unpaired:
		no_argument+=' --unpaired'
	logger.info(f'** cmd line : python lotus.py filter -v {str(vcf_file)} -o {output1} -wm {working_method} --DP {filter_dp} --MBQ {filter_mbq} --POPAF {filter_popaf}'+str(no_argument)+' **')
	logger.info('* Start filtering *')
	logger.info(f'Working directory (vcf files folder) : {working_directory}')
	logger.info(f'Current directory : {Path().absolute()}')

	filter(dict_colors, dict_para, output1, output2, vcf_file, logger, working_method, filter_ad, filter_qual, filter_mbq, filter_dp, filter_mqsbz, filter_popaf, unpaired)
	# End
