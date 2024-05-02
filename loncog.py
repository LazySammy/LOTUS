#!/usr/bin/env python3
# coding: utf-8

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#   LOncoG : a software for Longitudinal OncoGenomics analysis
#   Authors: S. Besseau, G. Siekaniec, W. Gouraud
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import os
import logging
from playsound import playsound
import pandas as pd
import sys
import time
import warnings
from python_scripts.modules import filter, summarise, compare, merge
from python_scripts.manage_parameters import calculate_time
from python_scripts.manage_parameters import list_colors
from python_scripts.manage_parameters import compare_parameters
from python_scripts.manage_parameters import manage_parameters

warnings.filterwarnings("ignore")

# real run or toy run?
toy_dataset = 'FALSE'
with open('config.txt', "r") as file:
	lines = file.readlines()
	for line in lines:
		if "use_toy_dataset" in line:
			toy_dataset = line.split("use_toy_dataset =")[1].strip().upper()
if toy_dataset == 'YES' or toy_dataset == 'TRUE':
	toy_dataset = True
	config_path = 'toy_dataset/toy_config.txt'
else:
	toy_dataset = False
	config_path = 'config.txt'

dict_para = manage_parameters(config_path, False)
if toy_dataset:
	dict_para['config_path'] = 'toy_dataset/toy_config.txt'
else:
	dict_para['config_path'] = 'config.txt'

if dict_para['colored_execution'].upper() == 'TRUE' or dict_para['colored_execution'].upper() == 'YES':
	colored = True
else:
	colored = False
dict_colors = list_colors()
if not colored:
	yellow = ''
	yellow1 = ''
	yellow2 = ''
	yellow3 = ''
	yellow4 = ''
	yellow5 = ''
	purple = ''
	dark_blue = ''
	dark_red = ''
	pink = ''
	bold = ''
	reset_bold = ''
else:
	yellow = dict_colors['yellow']
	yellow1 = dict_colors['yellow1']
	yellow2 = dict_colors['yellow2']
	yellow3 = dict_colors['yellow3']
	yellow4 = dict_colors['yellow4']
	yellow5 = dict_colors['yellow5']
	purple = dict_colors['purple']
	dark_blue = dict_colors['dark_blue']
	dark_red = dict_colors['dark_red']
	pink = dict_colors['pink']
	bold = dict_colors['bold']
	reset_bold = dict_colors['reset_bold']

print(f"{yellow1}{bold}Reading parameters...{reset_bold}")
start_time = time.time()
parameters_time = time.time()
previous = dict_para['output_folder_path']
modules = dict_para['module(s)'].upper()
output_path = dict_para['output_path']
# if previous.upper() != 'FALSE' and not os.path.isdir(output_path + 'samples'):
#     sys.exit(f"{dark_red}{bold}ERROR: {reset_bold}Without filter and summarise files, you can't run compare or merge modules :"
#              f" \n-> please run LOTUS with filter and summarise modules in this folder, or choose another previous_folder name with a samples folder."
#              f"\n{pink}-> the samples folder must contain every output from filter and summarise modules, for each sample folder.")

print(f"{yellow1}Working in {bold}{output_path}{reset_bold} directory...")

if not os.path.exists(output_path) and (previous.upper() != 'FALSE' or previous.upper() != 'NO' or previous.upper() != 'NONE'):
	os.makedirs(output_path)

if dict_para['module(s)'].replace(" ", "").upper() != 'COMPARE':
	i = 0
	for file in dict_para['vcf_files']:
		if file.endswith(".vcf"):
			sample_id = file.split('.funco')[0]
			output_path_sample = output_path + 'samples/' + sample_id + '/'
			dict_para['output_path_sample'] = output_path_sample
			if '.vcf' in output_path_sample:
				output_path_sample = output_path_sample.replace('.vcf', '')
			if not os.path.exists(output_path_sample):
				os.makedirs(output_path_sample)
			index = dict_para['vcf_files'].index(file) + 1
			pre_filter_time = time.time()

			if 'FILTER' in modules:
				file_extension = dict_para['dataset_path'].split('.')[-1]
				if file_extension == 'csv':
					try:
						df = pd.read_csv(dict_para['dataset_path'], delimiter=',')
					except:
						df = pd.read_excel(dict_para['dataset_path'])
				elif file_extension == 'xlsx':
					df = pd.read_excel(dict_para['dataset_path'])
				elif file_extension == 'tsv':
					df = pd.read_csv(dict_para['dataset_path'], delimiter='\t')
				else:
					raise ValueError("Unsupported file format.")

				vcf_folder_path = dict_para['vcf_folder_path']
				file_count = 0

				vcf_names = []
				for file_name in os.listdir(vcf_folder_path):
					if file_name.endswith('.vcf'):
						file_count += 1
						vcf_names.append(file_name)

				print(f"\n{yellow2}{bold}Filtering vcf file {index}/{file_count} ({sample_id})...{reset_bold}")

				# check if you have weird folders to delete in samples/ folder
				vcf_names = [name.replace('.vcf', '') if name.endswith('.vcf') else name for name in vcf_names]
				if i == 0:
					output_samples_path = dict_para['output_path_sample'].split('samples/')[0] + '/samples/'
					for root, dirs, files in os.walk(output_samples_path, topdown=False):
						for dir_name in dirs:
							dir_path = os.path.join(root, dir_name)
							if dir_name not in vcf_names:
								if os.path.exists(dir_path):
									if os.path.isdir(dir_path):
										file_list = os.listdir(dir_path)
										for file_name in file_list:
											file_path = os.path.join(dir_path, file_name)
											if os.path.isfile(file_path):
												os.remove(file_path)
										os.rmdir(dir_path)

				if previous.upper() == 'FALSE':
					input_vcf = (dict_para['vcf_folder_path'] + file)
					filter.main(dict_colors, dict_para, input_vcf, False)
					filter_time = calculate_time(dict_para, pre_filter_time, 'Filter')
				elif previous.upper() != 'FALSE':
					try:
						for filename in os.listdir(output_path_sample):
							file_path = os.path.join(output_path_sample, filename)
							if filename.endswith(".vcf"):
								os.remove(file_path)
					except:
						pass
					input_vcf = dict_para['vcf_folder_path'] + file
					if not os.path.exists(dict_para['output_path_sample'].replace(".vcf", "")) and (
							previous.upper() != 'FALSE' or previous.upper() != 'NO' or previous.upper() != 'NONE'):
						os.makedirs(dict_para['output_path_sample'].replace(".vcf", ""))
					filter.main(dict_colors, dict_para, input_vcf, False)
					filter_time = calculate_time(dict_para, pre_filter_time, 'Filter')
			if dict_para['keep_filtered_vcf_after_run'].upper() == 'FALSE' or dict_para['keep_filtered_vcf_after_run'].upper() == 'NO':
				if os.path.exists(output_path_sample + output_path_sample.split('/')[-2] + '_filtered.vcf'):
					os.remove(output_path_sample + output_path_sample.split('/')[-2] + '_filtered.vcf')
			filter_time = time.time()

			if 'SUMMARISE' in modules and previous.upper() == 'FALSE' and 'FILTER' not in modules:
				sys.exit("We can't use summarise module without filter before, or without using a previus folder containing filter files.")
			elif 'SUMMARISE' in modules:
				file_extension = dict_para['dataset_path'].split('.')[-1]
				if file_extension == 'csv':
					try:
						df = pd.read_csv(dict_para['dataset_path'], delimiter=',')
					except:
						df = pd.read_excel(dict_para['dataset_path'])
				elif file_extension == 'xlsx':
					df = pd.read_excel(dict_para['dataset_path'])
				elif file_extension == 'tsv':
					df = pd.read_csv(dict_para['dataset_path'], delimiter='\t')
				else:
					raise ValueError("Unsupported file format.")
				vcf_folder_path = dict_para['vcf_folder_path']
				file_count = 0

				for file_name in os.listdir(vcf_folder_path):
					if file_name.endswith('.vcf'):
						file_count += 1

				if dict_para['output_folder_path'].upper() != 'FALSE':
					for filename in os.listdir(output_path):
						if filename.endswith("xlsx") or filename.endswith("tsv") or filename.endswith("csv"):
							os.remove(output_path + filename)
					for filename in os.listdir(output_path_sample):
						file_path = os.path.join(output_path_sample, filename)
						if not filename.endswith(".vcf") and not filename.endswith("stats.txt"):
							os.remove(file_path)
				print(f"\n{yellow3}{bold}Summarising vcf file {index}/{file_count} ({sample_id}):{reset_bold}")
				input_vcf_filtered = output_path_sample + output_path_sample.split('samples/')[1].replace('/', "_") + 'filtered.vcf'
				input_vcf_passed = output_path_sample + output_path_sample.split('samples/')[1].replace('/', "_") + 'passed.vcf'
				if i != len(dict_para['vcf_files']) - 1:
					last = False
				else:
					last = True
				summarise.main(dict_para, dict_colors, dict_para, input_vcf_filtered, input_vcf_passed, output_path_sample, last)
				summarise_time = calculate_time(dict_para, filter_time, 'Summarise')
			i += 1

if 'COMPARE' in modules:
	# y a-t-il les fichiers filtered?
	if os.path.exists(output_path + 'samples'):
		for file_name in os.listdir(output_path + 'samples/'):
			if os.path.isdir(output_path + 'samples/' + file_name + '/'):
				for file in os.listdir(output_path + 'samples/' + file_name + '/'):
					if file.endswith(".vcf"):
						if file.endswith("filtered.vcf"):
							break
						elif dict_para['keep_filtered_vcf_after_run'].upper() == 'TRUE' or dict_para['keep_filtered_vcf_after_run'].upper() == 'YES':
							input_vcf = dict_para['vcf_folder_path'] + file_name + '.vcf'
							print(f"\n{yellow2}{bold}Filtering to recover {file_name}.vcf ...{reset_bold}")
							output_path_sample = output_path + 'samples/' + file_name + '/'
							filter.main(dict_colors, dict_para, input_vcf, output_path_sample, True)
	if 'XLSX' in dict_para['dataset_path'].upper():
		df = pd.read_excel(dict_para['dataset_path'])
	elif 'CSV' in dict_para['dataset_path'].upper():
		try:
			df = pd.read_csv(dict_para['dataset_path'])
		except:
			df = pd.read_excel(dict_para['dataset_path'])
	elif 'TSV' in dict_para['dataset_path'].upper():
		try:
			df = pd.read_csv(dict_para['dataset_path'], delimiter='\t')
		except:
			df = pd.read_excel(dict_para['dataset_path'])
	n_comp = 0
	counter = 0
	for value in df[dict_para['time1_column_name']]:
		if str(value) != 'nan' and not str(value).isspace():
			n_comp += 1
		counter += 1
	if previous.upper() != 'FALSE' and previous.upper() != 'NO' and previous.upper() != 'NONE':
		if os.path.exists(output_path + 'comparisons'):
			for file_name in os.listdir(output_path + 'comparisons/'):
				if os.path.isdir(output_path + 'comparisons/' + file_name):
					for file in os.listdir(output_path + 'comparisons/' + file_name + '/'):
						os.remove(output_path + 'comparisons/' + file_name + '/' + file)
					os.rmdir(output_path + 'comparisons/' + file_name)
				else:
					os.remove(output_path + 'comparisons/' + file_name)

	for i in range(n_comp):
		before_compare = time.time()
		dict_para = compare_parameters(dict_para, i)
		print(f"\n{yellow4}{bold}Comparing vcf pair {i + 1}/{n_comp} ({dict_para['file1']} & {dict_para['file2']}):{reset_bold}")
		compare.main(dict_para, df)
		compare_time = calculate_time(dict_para, before_compare, 'Compare')

if previous.upper() != 'FALSE':
	compare_time = time.time()

if 'MERGE' in modules:
	if dict_para['keep_filtered_vcf_after_run'].upper() == 'FALSE' or dict_para['keep_filtered_vcf_after_run'].upper() == 'NO':
		folder_path = dict_para['output_folder_path'] + 'samples/'
		for folder_name in os.listdir(folder_path):
			folder = os.path.join(folder_path, folder_name)
			if os.path.isdir(folder):
				for file_name in os.listdir(folder):
					file = os.path.join(folder, file_name)
					if file.endswith('filtered.vcf'):
						os.remove(file)
	n_pairs = sum(os.path.isdir(os.path.join(output_path + "comparisons/", item)) for item in os.listdir(output_path + "comparisons/"))
	pair = 'pair' if n_pairs == 1 else 'pairs'
	print(f"\n{yellow5}{bold}Merging {n_pairs} {pair}:{reset_bold}")
	output_path_merge = output_path + 'merge/'
	if os.path.exists(output_path_merge) and previous.upper() == 'FALSE':
		os.rmdir(output_path_merge)
	elif os.path.exists(output_path_merge) and previous.upper() != 'FALSE':
		for file_name in os.listdir(output_path_merge):
			file_path = os.path.join(output_path_merge, file_name)
			os.remove(file_path)
		os.rmdir(output_path_merge)
	os.makedirs(output_path_merge)

	dict_para['M_MutatedGenes_name'] = output_path_merge + 'Mutated_Genes'
	dict_para['chromosomes_plot_name'] = output_path_merge + 'Chromosomes_map'
	merge.main(dict_para, output_path)
	merge_time = calculate_time(dict_para, compare_time, 'Merge')

audio_file = "input/resources/job_done.wav"
playsound(audio_file)
final_time = calculate_time(dict_para, start_time, 'End')

__version__ = '0.0.1'
ascii_summarise = r'''                                                      _
                                                     (_)         
             ___ _   _ _ __ ___  _ __ ___   __ _ _ __ _ ___  ___ 
            / __| | | | '_ ` _ \| '_ ` _ \ / _` | '__| / __|/ _ \
            \__ \ |_| | | | | | | | | | | | (_| | |  | \__ \  __/
            |___/\__,_|_| |_| |_|_| |_| |_|\__,_|_|  |_|___/\___|                                                         

-----------------------------------------------------------------------------------

-----------------------------------------------------------------------------------
'''

ascii_filter = r'''                              __ _ _ _            
                             / _(_) | |           
                            | |_ _| | |_ ___ _ __ 
                            |  _| | | __/ _ \ '__|
                            | | | | | ||  __/ |   
                            |_| |_|_|\__\___|_|

-----------------------------------------------------------------------------------

-----------------------------------------------------------------------------------
'''

ascii_compare = r'''                  ___ ___  _ __ ___  _ __   __ _ _ __ ___ 
                 / __/ _ \| '_ ` _ \| '_ \ / _` | '__/ _ \
                | (_| (_) | | | | | | |_) | (_| | | |  __/
                 \___\___/|_| |_| |_| .__/ \__,_|_|  \___|
                                    | |                   
                                    |_|                  

-----------------------------------------------------------------------------------

-----------------------------------------------------------------------------------
'''

ascii_merge = r'''                        _ __ ___   ___ _ __ __ _  ___ 
                       | '_ ` _ \ / _ \ '__/ _` |/ _ \
                       | | | | | |  __/ | | (_| |  __/
                       |_| |_| |_|\___|_|  \__, |\___|
                                            __/ |     
                                           |___/

-----------------------------------------------------------------------------------

-----------------------------------------------------------------------------------
'''

# Logger configuration
if dict_para["log"].upper() != 'NONE':
	logger = logging.getLogger('LOncoG main')
	logger.setLevel(logging.DEBUG)
	logger.info(f'---------------- LOncoG v{__version__} ----------------')
