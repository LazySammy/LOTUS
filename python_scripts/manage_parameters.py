import os
import pandas as pd
import time
from datetime import datetime
import shutil
import sys
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#   LOncoG : a software for Longitudinal OncoGenomics analysis
#   Authors: S. Besseau, G. Siekaniec, W. Gouraud
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def manage_parameters(config_file, recover):

	current_time = datetime.now().strftime('%d-%m-%y(%Hh%M)')  # Utilisez '_' au lieu de '|'
	dict_parameters = {}

	with open(config_file, 'r') as file:
		# getting module name
		lines = file.readlines()
		for line in lines:
			line = line.strip()
			if line.startswith("module"):
				parameter, module_value = line.split("=")
				module_value = module_value.upper()
				if recover:
					module_value = module_value + 'FILTER'

		# getting parameters
		in_global_section = False
		filtered_lines = []
		for line in lines:
			if line.startswith("### Global") or section_line == 'global':
				section_line = "global"
			if (line.startswith("### Filter") or section_line == 'filter'):
				section_line = "filter"
			if (line.startswith("### Summarise") or section_line == 'summarise'):
				section_line = "summarise"
			if (line.startswith("### Compare") or section_line == 'compare'):
				section_line = "compare"
			if (line.startswith("### Merge") or section_line == "merge"):
				section_line = "merge"
			if line.startswith("# "):
				section_line = ''
			if section_line == 'global' or section_line in module_value.lower():
				filtered_lines.append(line)
			if line.startswith('variants_selection_approach'):
				dict_parameters['variants_selection_approach'] = line.split('=')[1].strip()

		# filling parameters dictionary with parameters
		for line in filtered_lines:
			line.strip().replace(" ", "")
			if not line.startswith("#"):
				if "=" in line:
					try:
						parameter, value = line.split("=")
					except:
						print('a')
					parameter = parameter.replace(" ", "")
					value = value.strip()
					if parameter == 'reference_genome' or parameter == 'gff3' or parameter == 'cytoband_file':
						value = value.replace("\\", "/")
						dict_parameters[parameter] = value
					elif parameter == 'vcf_S':
						value = 'input/vcf/' + value
						dict_parameters[parameter] = value
					else:
						dict_parameters[parameter] = value
	output_folder = dict_parameters['output_folder_path']
	if output_folder.upper() != 'FALSE':
		if not output_folder.endswith('/'):
			dict_parameters['output_folder_path'] = output_folder + '/'
			output_folder = dict_parameters['output_folder_path']

	if 'S_enrichment' in dict_parameters.keys():
		enrichment = dict_parameters['S_enrichment']
		if 'YES' in enrichment or 'TRUE' in enrichment or 'PANTHER' in enrichment or 'TOPPGENE' in enrichment:
			dict_parameters['S_enrichment'] = True

	if output_folder.upper() == 'FALSE':
		output_path = 'output/' + current_time + '/'
		i = 0
		while os.path.exists(output_path):
			if os.path.exists(output_path):
				i += 1
				output_path = 'output/' + current_time + i * "'" + '/'
		os.mkdir(output_path)

		# writing parameters in a file
		with open(output_path + "parameters.txt", 'w') as output_file:
			for line in filtered_lines:
				if not line.startswith("#"):
					cleaned_line = line.replace('\n', '')
				elif 'GLOBAL' in line.upper():
					cleaned_line = line[:-1]
				else:
					cleaned_line = '\n' + line[:-1]
				if cleaned_line.strip():
					if cleaned_line.startswith("# " + module_value):
						output_file.write('\n')
					output_file.write(cleaned_line + '\n')

		input_vcf = dict_parameters['vcf_folder_path']
		if not input_vcf.endswith('/'):
			dict_parameters['vcf_folder_path'] = input_vcf + '/'
			input_vcf += '/'
		vcf_files = [filename for filename in os.listdir(input_vcf) if filename.endswith('.vcf')]
		os.makedirs(output_path + 'samples/')
		dict_parameters['vcf_files'] = vcf_files
		dict_parameters['output_path'] = output_path


	else:
		output_path = output_folder
		# writing parameters in another file
		with open(output_path + "parameters.txt", 'w') as output_file:
			for line in filtered_lines:
				if not line.startswith("#"):
					cleaned_line = line.replace('\n', '')
				elif 'GLOBAL' in line.upper():
					cleaned_line = line[:-1]
				else:
					cleaned_line = '\n' + line[:-1]
				if cleaned_line.strip():
					if cleaned_line.startswith("# " + module_value):
						output_file.write('\n')
					output_file.write(cleaned_line + '\n')

		input_vcf = dict_parameters['vcf_folder_path']
		if not input_vcf.endswith('/'):
			dict_parameters['vcf_folder_path'] = input_vcf + '/'
			input_vcf += '/'
		vcf_files = [filename for filename in os.listdir(input_vcf) if ".vcf" in filename]
		dict_parameters['vcf_files'] = vcf_files

		dict_parameters['output_path'] = output_path

	if 'COMPARE' in module_value:
		shutil.copy(dict_parameters['dataset_path'], output_path)

	# writing parameters in a file
	with open(output_path + "parameters.txt", 'w') as output_file:
		for line in filtered_lines:
			if not line.startswith("#"):
				cleaned_line = line.replace('\n', '')
			elif 'GLOBAL' in line.upper():
				cleaned_line = line[:-1]
			else:
				cleaned_line = '\n' + line[:-1]
			if cleaned_line.strip():
				if cleaned_line.startswith("# " + module_value):
					output_file.write('\n')
				output_file.write(cleaned_line + '\n')

	return dict_parameters




def compare_parameters(dict_para, pair_index):
	df_dataset = pd.read_excel(dict_para['dataset_path'])
	output_path = dict_para['output_path']
	l_col1_names = []
	l_col2_names = []
	l_col_pair_names = []

	counter = 0
	for value in df_dataset[dict_para['time1_column_name']]:
		if str(value) != 'nan' and not str(value).isspace():
			l_col1_names.append(value)
		counter += 1
	for value in df_dataset[dict_para['time2_column_name']]:
		if str(value) != 'nan' and not str(value).isspace():
			l_col2_names.append(value)
		counter += 1
	for value in df_dataset[dict_para['pair_names_column_name']]:
		if str(value) != 'nan' and not str(value).isspace():
			l_col_pair_names.append(value)

	for i in range(len(l_col1_names)):
		l_col1_names[i] = l_col1_names[i].split('.funco')[0]
	for i in range(len(l_col2_names)):
		l_col2_names[i] = l_col2_names[i].split('.funco')[0]
	n_lines1 = len(l_col1_names)
	n_lines2 = len(l_col2_names)
	if n_lines1 != n_lines2:  # changer pour vérifier si bien en face de l'autre
		print('The number of samples in the two comparison columns is not the same. Please make sure to have pairs.')
		exit(2)

	output_folder = dict_para['output_folder_path']
	if pair_index == 0 and output_folder.upper() == 'FALSE' and not os.path.exists(output_path + 'comparisons/'):
		os.makedirs(output_path + 'comparisons/')
	elif pair_index == 0 and output_folder.upper() != 'FALSE' and not os.path.exists(output_path + 'comparisons/'):
		os.makedirs(output_path + 'comparisons/')
	elif pair_index == 0 and output_folder.upper() != 'FALSE' and os.path.exists(output_path + 'comparisons/'):
		comp_path = output_path + 'comparisons/'
		for folder_name in os.listdir(comp_path):
			folder_path = os.path.join(comp_path, folder_name) + '/'
			file_list = os.listdir(folder_path)
			for file in file_list:
				file_path = os.path.join(folder_path, file)
				os.remove(file_path)
			os.rmdir(folder_path)
		os.rmdir(comp_path)
		os.makedirs(output_path + 'comparisons/')

	if not l_col_pair_names == []:
		comp_id = l_col_pair_names[pair_index]
		comp_id2 = l_col1_names[pair_index] + '___' + l_col2_names[pair_index]
		l_comp = comp_id2.split('___')
	else:
		comp_id = l_col1_names[pair_index] + '___' + l_col2_names[pair_index]
		comp_id = comp_id.replace(".vcf", "")
		l_comp = comp_id.split('___')
		l_comp[0] = l_comp[0].lstrip().rstrip()
		l_comp[1] = l_comp[1].lstrip().rstrip()

	output_path_comparison = output_path + 'comparisons/' + comp_id + '/'
	os.makedirs(output_path_comparison)

	# input_vcf_passed1 = 'output/greg/samples/' + l_comp[0].split('.fun')[0] + '/' + l_comp[0].split('ESM')[0] + '.passed.vcf'
	# input_vcf_passed2 = 'output/greg/samples/' + l_comp[1].split('.fun')[0] + '/' + l_comp[1].split('ESM')[0] + '.passed.vcf'
	input_vcf_passed1 = output_path + 'samples/' + l_comp[0] + '/passed.vcf'
	input_vcf_passed2 = output_path + 'samples/' + l_comp[1] + '/passed.vcf'
	# input_vcf_filtered1 = 'output/greg/samples/' + l_comp[0].split('.fun')[0] + '/' + l_comp[0].split('ESM')[0] + '.filtered.vcf'
	# input_vcf_filtered2 = 'output/greg/samples/' + l_comp[1].split('.fun')[0] + '/' + l_comp[1].split('ESM')[0] + '.filtered.vcf'
	input_vcf_filtered1 = output_path + 'samples/' + l_comp[0] + '/filtered.vcf'
	input_vcf_filtered2 = output_path + 'samples/' + l_comp[1] + '/filtered.vcf'
	# input_vcf_indel1 = 'output/greg/samples/' + l_comp[0].split('.fun')[0] + '/' + l_comp[0].split('ESM')[0] + '_indel.deletion.tsv'
	# input_vcf_indel2 = 'output/greg/samples/' + l_comp[1].split('.fun')[0] + '/' + l_comp[1].split('ESM')[0] + '_indel.deletion.tsv'
	input_vcf_indel1 = output_path + 'samples/' + l_comp[0] + '/indel.deletion.tsv'
	input_vcf_indel2 = output_path + 'samples/' + l_comp[1] + '/indel.deletion.tsv'
	# input_vcf_insert1 = 'output/greg/samples/' + l_comp[0].split('.fun')[0] + '/' + l_comp[0].split('ESM')[0] + '_indel.insertion.tsv'
	# input_vcf_insert2 = 'output/greg/samples/' + l_comp[1].split('.fun')[0] + '/' + l_comp[1].split('ESM')[0] + '_indel.insertion.tsv'
	input_vcf_insert1 = output_path + 'samples/' + l_comp[0] + '/indel.insertion.tsv'
	input_vcf_insert2 = output_path + 'samples/' + l_comp[1] + '/indel.insertion.tsv'
	# input_vcf_snp1 = 'output/greg/samples/' + l_comp[0].split('.fun')[0] + '/' + l_comp[0].split('ESM')[0] + '_profile.tsv'
	# input_vcf_snp2 = 'output/greg/samples/' + l_comp[1].split('.fun')[0] + '/' + l_comp[1].split('ESM')[0] + '_profile.tsv'sv'
	out_indel = output_path_comparison + 'indel.svg'

	dict_comp = {'file1': l_comp[0], 'file2': l_comp[1], 'output_path_comparison': output_path_comparison,
				 'passed1': input_vcf_passed1, 'passed2': input_vcf_passed2,
				 'filtered1': input_vcf_filtered1, 'filtered2': input_vcf_filtered2,
				 'indel1': input_vcf_indel1, 'indel2': input_vcf_indel2,
				 'insert1': input_vcf_insert1, 'insert2': input_vcf_insert2, 'indel': out_indel}
	dict_para.update(dict_comp)

	return dict_para


def list_colors():
	colored_execution = True
	dict_colors = {}
	forest_green = "\033[32m"
	purple = "\033[35m"
	yellow = "\033[33m"
	red = "\033[31m"
	dark_blue = "\033[34m"
	dark_red = "\033[31m"
	pink = "\033[95m"
	reset_color = "\033[0m"
	bold = "\033[1m"
	reset_bold = "\033[22m"
	italic = "\033[3m"
	reset_italic = "\033[23m"

	yellow1 = '\033[38;2;255;244;213m'
	yellow2 = '\033[38;2;255;231;104m'
	yellow3 = '\033[38;2;255;215;0m'
	yellow4 = '\033[38;2;255;193;7m'
	yellow5 = '\033[38;2;255;152;0m'

	dict_colors['forest_green'] = forest_green
	dict_colors['purple'] = purple
	dict_colors['yellow'] = yellow
	dict_colors['red'] = red
	dict_colors['dark_blue'] = dark_blue
	dict_colors['dark_red'] = dark_red
	dict_colors['pink'] = pink
	dict_colors['reset_color'] = reset_color
	dict_colors['bold'] = bold
	dict_colors['reset_bold'] = reset_bold
	dict_colors['italic'] = italic
	dict_colors['reset_italic'] = reset_italic

	dict_colors['yellow1'] = yellow1
	dict_colors['yellow2'] = yellow2
	dict_colors['yellow3'] = yellow3
	dict_colors['yellow4'] = yellow4
	dict_colors['yellow5'] = yellow5

	return dict_colors


def calculate_time(dict_para, time1, module):
	time_now = time.time()
	minutes = int((time_now - time1) // 60)
	seconds = int((time_now - time1) % 60)

	if dict_para['colored_execution'].upper() == 'TRUE' or dict_para['colored_execution'].upper() == 'YES' or dict_para['verbose_prints'].upper() == 'YES':
		dict_colors = list_colors()
		color_reset = dict_colors['reset_color']
		module_format = dict_colors['italic']
		format_reset = dict_colors['reset_italic']
		bold = dict_colors['bold']
		module_color = dict_colors['forest_green']
	else:
		module_color = ''
		color_reset = ''
		module_format = ''
		format_reset = ''
		bold = ''

	if 'PARAMETERS' in module.upper():
		# module_color = dict_colors['yellow']
		module = module_color + 'Parameters preparation'
	elif 'FILTER' in module.upper():
		# module_color = dict_colors['purple']
		module = module_color + 'Filtering'
	elif 'SUMMARISE' in module.upper():
		# module_color = dict_colors['dark_blue']
		module = module_color + 'Summarising'
	elif 'COMPARE' in module.upper() :
		# module_color = dict_colors['dark_red']
		module = module_color + 'Comparing'
	elif 'MERGE' in module.upper() :
		# module_color = dict_colors['pink']
		module = module_color + 'Merging'
	elif 'END' in module.upper():
		# module_color = dict_colors['pink']
		module = module_color + '\nFull program'


	if dict_para['colored_execution'].upper() == 'TRUE' or dict_para['colored_execution'].upper() == 'YES' or dict_para['verbose_prints'].upper() == 'YES':
		if minutes == 0 and seconds == 0:
			print(f"{module_format}{module} "
				  f"done in {bold}{minutes}m and {seconds}s!{color_reset}{format_reset}")
		elif minutes == 0:
			print(f"{module_format}{module} "
				  f"done in {bold}{seconds}s!{color_reset}{format_reset}")
		elif seconds == 0:
			print(f"{module_format}{module} "
				  f"done in {bold}{minutes}m!{color_reset}{format_reset}")
		else:
			print(f"{module_format}{module} "
				  f"done in {bold}{minutes}m and {seconds}s!{color_reset}{format_reset}")
	elif dict_para['colored_execution'].upper() != 'TRUE' and dict_para['colored_execution'].upper() != 'YES' and dict_para['verbose_prints'].upper() != 'YES':
		if minutes == 0 and seconds == 0:
			print(f"{module} "
				  f"done in {minutes}m and {seconds}s!")
		elif minutes == 0:
			print(f"{module} "
				  f"done in {seconds}s!")
		elif seconds == 0:
			print(f"{module} "
				  f"done in {minutes}m!")
		else:
			print(f"{module} "
				  f"done in {minutes}m and {seconds}s!")

	time_now = time.time()
	return time_now
