#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pickle as pk
from math import isnan
from pathlib import Path
import sys
import os
import numpy as np
import logging
import re
import warnings
from tqdm import tqdm
import pandas as pd
import matplotlib.pyplot as plt
from copy import deepcopy
from more_itertools import powerset, unique_everseen
from upsetplot import from_memberships
from upsetplot import UpSet
from python_scripts.check_files import verif_output, verif_input_config_merge, verif_input_xlsx, verif_input_tsv, \
	verif_input, verif_supplementary_information_file
from python_scripts.toppgene_api import ToppGene_GOEA
from python_scripts.panther_api import Panther_GOEA
from python_scripts.chromosomes_plot import create_chromosomes_plot
import matplotlib.pyplot as plt


def get_informations_for_genes(args, info_file, logger):
	'''
	Extract informations from the external cancer databases file.
	'''
	df = pd.read_excel(info_file, index_col=1)
	df = df.drop(['Ordre'], axis=1)
	df.set_axis([source.split(' Info')[0] for source in df.columns], axis="columns")
	if args['verbose_prints'].upper() == 'TRUE':
		print(f'Extracting informations from {len(list(df.columns))} sources: {", ".join(list(df.columns))}')
	logger.info(f'Extract informations from {len(list(df.columns))} sources: {", ".join(list(df.columns))}')
	return (df)


def read_config(args, config_file):
	'''
	Read the merge config file (composed of one gene.tsv compare file per line).
	'''
	with open(config_file, 'r') as f:
		for line in f:
			line = line.strip()
			if line != '':
				line = line.split(',')
				if verif_input_xlsx(line[0]):
					data = pd.read_excel(line[0], index_col=0)
				elif verif_input_tsv(line[0]):
					data = pd.read_csv(line[0], sep='\t', index_col=0)
				if len(line) > 1:
					yield line[1], data
				else:
					yield line[0], data


def add_variants_to_dictionnary(df_all_genes_info: {}, gene_name: str, pos: int, field: str, row):
	'''
	Add a variant in count dictionnary.
	'''
	if field == 'gb1' or field == 'gb2' or field == 'cb1' or field == 'cb2' or field == 'pb1' or field == 'pb2':
		if df_all_genes_info[gene_name][field] == {}:
			df_all_genes_info[gene_name][field] = []
		if '|' in str(row.iloc[pos]):
			l_variants = row.iloc[pos].split('|')
			for variant in l_variants:
				df_all_genes_info[gene_name][field].append(variant)
		elif str(row.iloc[pos]).startswith('g.'):
			df_all_genes_info[gene_name][field].append(row.iloc[pos])
		elif 'nan' in str(row.iloc[pos]):
			df_all_genes_info[gene_name][field] = []
	elif field != 'mutation types':
		if not isnan(row.iloc[pos]):
			for i in row.iloc[pos]:
				try:
					if not i in df_all_genes_info[gene_name][field].keys():
						df_all_genes_info[gene_name][field][i] = 1
					else:
						df_all_genes_info[gene_name][field][i] += 1
				except:
					if not i in df_all_genes_info[gene_name][field]:
						df_all_genes_info[gene_name][field].append(row.iloc[pos])


def count_mutation_types(args, mutation_types_counting_method):
	if args['vcf_annotation_method'].upper() == 'ANNOVAR':
		annovar = True

	# Which file belongs do what time?
	df_path = args['dataset_path']
	if df_path.endswith('xlsx'):
		df = pd.read_excel(df_path)
	elif df_path.endswith('csv'):
		df = pd.read_csv(df_path)
	time1_column_values = df[args['time1_column_name']].tolist()
	time2_column_values = df[args['time2_column_name']].tolist()
	time1_column_values = [value for value in time1_column_values if not pd.isna(value)]
	time2_column_values = [value for value in time2_column_values if not pd.isna(value)]
	new_values_t1 = []
	for value in time1_column_values:
		if '.vcf' in value:
			new_values_t1.append(value.split('.vcf')[0])
		else:
			new_values_t1.append(value)
	time1_column_values = new_values_t1
	new_values_t2 = []
	for value in time2_column_values:
		if '.vcf' in value:
			new_values_t2.append(value.split('.vcf')[0])
		else:
			new_values_t2.append(value)
	time2_column_values = new_values_t2

	# Get passed files paths
	samples_path = args['output_path'] + 'samples/'
	vcf_passed_files = []
	for root, dirs, files in os.walk(samples_path):
		for file in files:
			if file.endswith('passed.vcf'):
				folder = os.path.basename(root)
				vcf_passed_files.append(samples_path + folder + '/' + file)

	t1_files = []
	t2_files = []
	for file in vcf_passed_files:
		if file.split('samples/')[1].split('/')[0] in time1_column_values:
			t1_files.append(file)
		for name in time2_column_values:
			if file.split('samples/')[1].split('/')[0] in name:
				t2_files.append(file)

	# Mutation types and subtypes counting
	l_unique_variants_t1 = []
	l_unique_variants_t2 = []
	dic_unique_mutation_types_t1 = {}
	dic_unique_mutation_types_t2 = {}
	dic_unique_mutation_subtypes_t1 = {}
	dic_unique_mutation_subtypes_t2 = {}
	dic_total_mutation_types_t1 = {}
	dic_total_mutation_types_t2 = {}
	dic_total_mutation_subtypes_t1 = {}
	dic_total_mutation_subtypes_t2 = {}

	for file in t1_files:
		with open(file, "r") as vcf_file:
			for line in vcf_file:
				if not line.startswith("#"):
					line = line.split("\t")
					chr = line[0]
					pos = line[1]
					ref = line[3]
					alt = line[4]
					variant = chr + ':' + pos + ':' + ref + ':' + alt

					if mutation_types_counting_method == 'UNIQUE':
						if not variant in l_unique_variants_t1:
							l_unique_variants_t1.append(variant)
							if len(ref) == len(alt):
								if len(ref) == 1:
									mutation_type = 'SNP'
								elif len(ref) == 2:
									mutation_type = 'DNP'
								elif len(ref) == 3:
									mutation_type = 'TNP'
								elif len(ref) > 3:
									mutation_type = 'ONP'
							elif len(ref) > len(alt):
								mutation_type = 'Deletion'
							elif len(ref) < len(alt):
								mutation_type = 'Insertion'
							if not mutation_type in dic_unique_mutation_types_t1.keys():
								dic_unique_mutation_types_t1[mutation_type] = 0
							dic_unique_mutation_types_t1[mutation_type] += 1

							if annovar:
								match = re.search(r'ExonicFunc\.refGene=(.*?);', str(line))
								if match:
									mutation_subtype = match.group(1)
									if not mutation_subtype in dic_unique_mutation_subtypes_t1.keys():
										dic_unique_mutation_subtypes_t1[mutation_subtype] = 0
									dic_unique_mutation_subtypes_t1[mutation_subtype] += 1

					elif mutation_types_counting_method == 'TOTAL':
						if len(ref) == len(alt):
							if len(ref) == 1:
								mutation_type = 'SNP'
							elif len(ref) == 2:
								mutation_type = 'DNP'
							elif len(ref) == 3:
								mutation_type = 'TNP'
							elif len(ref) > 3:
								mutation_type = 'ONP'
						elif len(ref) > len(alt):
							mutation_type = 'Deletion'
						elif len(ref) < len(alt):
							mutation_type = 'Insertion'
						if not mutation_type in dic_total_mutation_types_t1.keys():
							dic_total_mutation_types_t1[mutation_type] = 0
						dic_total_mutation_types_t1[mutation_type] += 1

						if annovar:
							match = re.search(r'ExonicFunc\.refGene=(.*?);', str(line))
							if match:
								mutation_subtype = match.group(1)
								if not mutation_subtype in dic_total_mutation_subtypes_t1.keys():
									dic_total_mutation_subtypes_t1[mutation_subtype] = 0
								dic_total_mutation_subtypes_t1[mutation_subtype] += 1

	for file in t2_files:
		with open(file, "r") as vcf_file:
			for line in vcf_file:
				if not line.startswith("#"):
					line = line.split("\t")
					chr = line[0]
					pos = line[1]
					ref = line[3]
					alt = line[4]
					variant = chr + ':' + pos + ':' + ref + ':' + alt

					if mutation_types_counting_method == 'UNIQUE':
						if not variant in l_unique_variants_t2:
							l_unique_variants_t2.append(variant)
							if len(ref) == len(alt):
								if len(ref) == 1:
									mutation_type = 'SNP'
								elif len(ref) == 2:
									mutation_type = 'DNP'
								elif len(ref) == 3:
									mutation_type = 'TNP'
								elif len(ref) > 3:
									mutation_type = 'ONP'
							elif len(ref) > len(alt):
								mutation_type = 'Deletion'
							elif len(ref) < len(alt):
								mutation_type = 'Insertion'
							if not mutation_type in dic_unique_mutation_types_t2.keys():
								dic_unique_mutation_types_t2[mutation_type] = 0
							dic_unique_mutation_types_t2[mutation_type] += 1

							if annovar:
								match = re.search(r'ExonicFunc\.refGene=(.*?);', str(line))
								if match:
									mutation_subtype = match.group(1)
									if not mutation_subtype in dic_unique_mutation_subtypes_t2.keys():
										dic_unique_mutation_subtypes_t2[mutation_subtype] = 0
									dic_unique_mutation_subtypes_t2[mutation_subtype] += 1

					elif mutation_types_counting_method == 'TOTAL':
						if len(ref) == len(alt):
							if len(ref) == 1:
								mutation_type = 'SNP'
							elif len(ref) == 2:
								mutation_type = 'DNP'
							elif len(ref) == 3:
								mutation_type = 'TNP'
							elif len(ref) > 3:
								mutation_type = 'ONP'
						elif len(ref) > len(alt):
							mutation_type = 'Deletion'
						elif len(ref) < len(alt):
							mutation_type = 'Insertion'
						if not mutation_type in dic_total_mutation_types_t2.keys():
							dic_total_mutation_types_t2[mutation_type] = 0
						dic_total_mutation_types_t2[mutation_type] += 1

						if annovar:
							match = re.search(r'ExonicFunc\.refGene=(.*?);', str(line))
							if match:
								mutation_subtype = match.group(1)
								if not mutation_subtype in dic_total_mutation_subtypes_t2.keys():
									dic_total_mutation_subtypes_t2[mutation_subtype] = 0
								dic_total_mutation_subtypes_t2[mutation_subtype] += 1

	if mutation_types_counting_method == 'UNIQUE' and annovar:
		return dic_unique_mutation_types_t1, dic_unique_mutation_types_t2, dic_unique_mutation_subtypes_t1, dic_unique_mutation_subtypes_t2
	elif mutation_types_counting_method == 'UNIQUE' and not annovar:
		return dic_unique_mutation_types_t1, dic_unique_mutation_types_t2
	elif mutation_types_counting_method == 'TOTAL' and annovar:
		return dic_total_mutation_types_t1, dic_total_mutation_types_t2, dic_total_mutation_subtypes_t1, dic_total_mutation_subtypes_t2
	elif mutation_types_counting_method == 'TOTAL' and not annovar:
		return dic_total_mutation_types_t1, dic_total_mutation_types_t2


# print('a')


def create_mutation_types_barplot(args, dic_mutation_types_t1, dic_mutation_types_t2, counting_method):
	plt.clf()
	labels = [label for label in dic_mutation_types_t1 if dic_mutation_types_t1[label] != 0 or dic_mutation_types_t2[label] != 0]
	values_t1 = [dic_mutation_types_t1[label] for label in labels]
	values_t2 = [dic_mutation_types_t2[label] for label in labels]
	x = np.arange(len(labels))

	bar_width = 0.2
	bar_gap = 0.1

	fig, ax = plt.subplots(figsize=(12, 8))

	rects1 = ax.bar(x - bar_width - bar_gap / 2, values_t1, bar_width, label='time 1', color=(0, 0.5, 0), edgecolor='black', linewidth=1, align="center")
	rects2 = ax.bar(x + bar_gap / 2, values_t2, bar_width, label='time 2', color=(0, 0, 0.5), edgecolor='black', linewidth=1, align="center")

	for rect in rects1:
		height = rect.get_height()
		ax.annotate('{}'.format(height),
					xy=(rect.get_x() + rect.get_width() / 2, height),
					xytext=(0, 3),
					textcoords="offset points",
					ha='center', va='bottom',
					fontsize=15)

	for rect in rects2:
		height = rect.get_height()
		ax.annotate('{}'.format(height),
					xy=(rect.get_x() + rect.get_width() / 2, height),
					xytext=(0, 3),
					textcoords="offset points",
					ha='center', va='bottom',
					fontsize=15)

	ax.set_ylabel('Count', fontsize=14.5, labelpad=10)
	mode = counting_method.lower()
	mode = mode.capitalize()
	ax.set_title(mode + ' mutation types comparison between time 1 and time 2\n(crossover of all patients from the cohort)', fontsize=16, pad=10)
	ax.set_xticks(x)
	ax.set_xticklabels(labels, fontsize=15)
	ax.set_ylim(0, max(max(values_t1), max(values_t2)) * 1.1)

	offset = bar_width / 2
	ax.set_xticks(x - offset, minor=False)
	ax.set_xticklabels(labels, minor=False)
	plt.tick_params(bottom=False)

	ax.legend(fontsize=14, edgecolor='white')
	plt.tight_layout()

	file_formats = args['M_types_barplot_format(s)'].upper()
	path = args['output_path'] + 'merge/' + args['M_types_barplot_name']
	if 'PNG' in file_formats:
		plt.savefig(path + '.png', dpi=400)
	if 'PDF' in file_formats:
		plt.savefig(path + '.pdf', dpi=400)
	if 'SVG' in file_formats:
		plt.savefig(path + '.svg', dpi=400)
	if 'JPG' in file_formats:
		plt.savefig(path + '.jpg', dpi=400)

	plt.close()


def create_mutation_subtypes_barplot(args, dic_mutation_subtypes_t1, dic_mutation_subtypes_t2, counting_method):
	if '.' in dic_mutation_subtypes_t1.keys():
		dic_mutation_subtypes_t1['unknown'] += dic_mutation_subtypes_t1['.']
		del dic_mutation_subtypes_t1['.']
	if '.' in dic_mutation_subtypes_t2.keys():
		if 'unknown' not in dic_mutation_subtypes_t2.keys():
			dic_mutation_subtypes_t2['unknown'] = 0
		dic_mutation_subtypes_t2['unknown'] += dic_mutation_subtypes_t2['.']
		del dic_mutation_subtypes_t2['.']

	plt.clf()
	for key in dic_mutation_subtypes_t1.keys():
		if not key in dic_mutation_subtypes_t2.keys():
			dic_mutation_subtypes_t2[key] = 0
	labels = [label for label in dic_mutation_subtypes_t1 if dic_mutation_subtypes_t1[label] != 0 and dic_mutation_subtypes_t2[label] != 0]
	values_t1 = [dic_mutation_subtypes_t1[label] for label in labels]
	values_t2 = [dic_mutation_subtypes_t2[label] for label in labels]
	x = np.arange(len(labels))

	bar_width = 0.3
	bar_gap = 0.2

	fig, ax = plt.subplots(figsize=(12, 12))

	x_indexes = np.arange(len(labels))
	x1 = x_indexes - bar_width / 2 - bar_gap / 2
	x2 = x_indexes + bar_width / 2 + bar_gap / 2

	rects1 = ax.bar(x1, values_t1, bar_width, label='time 1', color=(0, 0.5, 0), edgecolor='black', linewidth=1, align="center")
	rects2 = ax.bar(x2, values_t2, bar_width, label='time 2', color=(0, 0, 0.5), edgecolor='black', linewidth=1, align="center")

	for rect in rects1:
		height = rect.get_height()
		ax.annotate('{}'.format(height),
					xy=(rect.get_x() + rect.get_width() / 2, height),
					xytext=(0, 3),
					textcoords="offset points",
					ha='center', va='bottom',
					fontsize=13)

	for rect in rects2:
		height = rect.get_height()
		ax.annotate('{}'.format(height),
					xy=(rect.get_x() + rect.get_width() / 2, height),
					xytext=(0, 3),
					textcoords="offset points",
					ha='center', va='bottom',
					fontsize=13)

	ax.set_ylabel('Count', fontsize=14.5, labelpad=14)

	mode = counting_method.lower()
	mode = mode.capitalize()
	ax.set_title(mode + ' mutation subtypes comparison between time 1 and time 2\n(crossover of all patients from the cohort)', fontsize=16, pad=10)

	ax.set_xticks(x_indexes)
	ax.set_xticklabels(labels, rotation=90, fontsize=13, ha='right')  # Adjust ha='right'
	ax.set_ylim(0, max(max(values_t1), max(values_t2)) * 1.1)

	ax.tick_params(axis='y', labelsize=13)
	ax.tick_params(axis='x', pad=10)  # Add padding to x-axis ticks
	plt.tick_params(bottom=False)

	ax.legend(fontsize=14, edgecolor='white')
	plt.tight_layout()

	fig.autofmt_xdate()
	file_formats = args['M_subtypes_barplot_format(s)'].upper()
	path = args['output_path'] + 'merge/' + args['M_subtypes_barplot_name']
	if 'PNG' in file_formats:
		plt.savefig(path + '.png', dpi=400)
	if 'PDF' in file_formats:
		plt.savefig(path + '.pdf', dpi=400)
	if 'SVG' in file_formats:
		plt.savefig(path + '.svg', dpi=400)
	if 'JPG' in file_formats:
		plt.savefig(path + '.jpg', dpi=400)

	plt.close()


def merge_results(args, file_paths, category, output, upset_output, infos, cytoband_file, chromosomes_output, step, nb_files, enrichment, logger):
	# Count mutation types and subtypes (and plot them)
	if args['vcf_annotation_method'].upper() == 'ANNOVAR':
		annovar = True
	mutation_types_counting_method = args['mutations_types_counting'].upper()

	if mutation_types_counting_method == 'UNIQUE' and annovar:
		print('Counting and plotting unique mutations types and subtypes...')
		dic_unique_mutation_types_t1, dic_unique_mutation_types_t2, dic_unique_mutation_subtypes_t1, dic_unique_mutation_subtypes_t2 = count_mutation_types(args, 'UNIQUE')
		create_mutation_types_barplot(args, dic_unique_mutation_types_t1, dic_unique_mutation_types_t2, 'UNIQUE')
		create_mutation_subtypes_barplot(args, dic_unique_mutation_subtypes_t1, dic_unique_mutation_subtypes_t2, 'UNIQUE')
	elif mutation_types_counting_method == 'UNIQUE' and not annovar:
		print('Counting and plotting unique mutations types...')
		dic_unique_mutation_types_t1, dic_unique_mutation_types_t2 = count_mutation_types(args, 'UNIQUE')
		create_mutation_types_barplot(args, dic_unique_mutation_types_t1, dic_unique_mutation_types_t2, 'UNIQUE')
	elif mutation_types_counting_method == 'TOTAL' and annovar:
		print('Counting and plotting unique mutations types and subtypes...')
		dic_total_mutation_types_t1, dic_total_mutation_types_t2, dic_total_mutation_subtypes_t1, dic_total_mutation_subtypes_t2 = count_mutation_types(args, 'TOTAL')
		create_mutation_types_barplot(args, dic_total_mutation_types_t1, dic_total_mutation_types_t2, 'TOTAL')
		create_mutation_subtypes_barplot(args, dic_total_mutation_subtypes_t1, dic_total_mutation_subtypes_t2, 'TOTAL')
	elif mutation_types_counting_method == 'TOTAL' and not annovar:
		print('Counting and plotting unique mutations types...')
		dic_total_mutation_types_t1, dic_total_mutation_types_t2 = count_mutation_types(args, 'TOTAL')
		create_mutation_types_barplot(args, dic_total_mutation_types_t1, dic_total_mutation_types_t2, 'TOTAL', )

	# Get genes infos if not None from the 6 databases
	if infos:
		df_databases_info = get_informations_for_genes(args, infos, logger)

	names = set()
	df_all_genes_info = {}  # gene dictionnary to create tsv union file (list format)

	if args['vcf_annotation_method'].upper() == 'FUNCOTATOR':
		dic_gene_fields = {'chr': '', 'start': '', 'end': '', 'weakness': [], 'gb1': {}, 'cb1': {}, 'pb1': {}, 'gb2': {}, 'cb2': {},
						   'pb2': {}, 'samples': []}
	elif args['vcf_annotation_method'].upper() == 'ANNOVAR':
		dic_gene_fields = {'chr': '', 'start': '', 'end': '', 'weakness': [], 'mutation types': [], 'gb1': {}, 'cb1': {}, 'pb1': {}, 'gb2': {},
						   'cb2': {}, 'pb2': {}, 'samples': []}
	if args['verbose_prints'].upper() == 'TRUE':
		print('Reading files...')

	if args['colored_execution'].upper() == 'TRUE' or args['colored_execution'].upper() == 'YES':
		color = '\033[38;2;255;152;0m'
		no_tqdm_bar = False
	else:
		color = ''
		no_tqdm_bar = True

	pbar_files = tqdm(total=int(nb_files / 2), disable=no_tqdm_bar, bar_format='{l_bar}{bar:30}{r_bar}', ncols=150, smoothing=1)
	pbar_files.set_description(color + f' -> Extracting information from \033[3mCompare\033[0m' + color + ' results (' + str(round(nb_files / 2)) + ' patients)')
	dic_unique_variants = {}
	dic_genes_in_pairs = {}

	for name in file_paths:
		if 'genes' in name:
			element = name.split("/")[-1]
			folder = name.split(element)[0]
			files = [f for f in os.listdir(folder) if os.path.isfile(os.path.join(folder, f))]
			for file in files:
				if 'genes' in file and 'xlsx' in file:
					name = folder + file
			id = name.split('sons/')[1].split('/')[0]
			dic_genes_in_pairs[id] = []
			if 'tsv' in name:
				df_from_compare = pd.read_csv(name, sep='\t', index_col=0)
			elif 'xlsx' in name:
				df_from_compare = pd.read_excel(name, index_col=0)
			elif 'csv' in name:
				df_from_compare = pd.read_csv(name, index_col=0)
			# if args['verbose_prints'].upper() == 'TRUE':
			#     print(name)
			if not name in names:
				names.add(name)
		# df_from_compare.to_excel('test.xlsx', index=False)
		# if args['vcf_annotation_method'].upper() == 'FUNCOTATOR':
		#     df_from_compare = df_from_compare.drop(df_from_compare.columns[0], axis=0) # drop the first column

		col_passed1 = df_from_compare.columns.values.tolist()[5]
		col_passed2 = df_from_compare.columns.values.tolist()[9]
		logger.info(f'Processing of samples B1: {col_passed1} and B2: {col_passed2}')  # b1 and b2 are the names of the samples
		# row is an object: doesn't have gene column (used as index)
		selection = args['C_variants_selection_approach'].upper()
		union = False
		change = False
		common = False
		if selection == 'UNION':
			union = True
		elif selection == 'CHANGE':
			change = True
		elif selection == 'COMMON':
			common = True

		if args['vcf_annotation_method'].upper() == 'FUNCOTATOR':
			passed_index = 5  # because we removed index column
			gb_index = 6
			cb_index = 7
			pb_index = 8
			passed1_index = 9
			gb1_index = 10
			cb1_index = 11
			pb1_index = 12

		elif args['vcf_annotation_method'].upper() == 'ANNOVAR':
			if union:
				passed_index = 6
				mutation_type = 7
				gb_index = 8
				cb_index = 9
				pb_index = 10
				passed1_index = 11
				gb1_index = 12
				cb1_index = 13
				pb1_index = 14
			elif change:
				passed_index = 5
				mutation_type = 6
				gb_index = 7
				cb_index = 8
				pb_index = 9
				passed1_index = 10
				gb1_index = 11
				cb1_index = 12
				pb1_index = 13
			elif common:
				passed_index = 4
				mutation_type = 5
				gb_index = 6
				cb_index = 7
				pb_index = 8
				passed1_index = 9
				gb1_index = 10
				cb1_index = 11
				pb1_index = 12

		variants_count = []

		for index, row in df_from_compare.iterrows():
			if not index in dic_genes_in_pairs[id]:
				dic_genes_in_pairs[id].append(index)

			if index not in df_all_genes_info.keys():
				df_all_genes_info[index] = deepcopy(dic_gene_fields)  # index is the gene name

			if df_all_genes_info[index]['chr'] == '':
				df_all_genes_info[index]['chr'] = row['Chromosome']

			if df_all_genes_info[index]['start'] == '':
				try:
					df_all_genes_info[index]['start'] = int(row['Gene start position'])
				except:
					df_all_genes_info[index]['start'] = 'error'
			if df_all_genes_info[index]['end'] == '':
				try:
					df_all_genes_info[index]['end'] = int(row['Gene end position'])
				except:
					df_all_genes_info[index]['start'] = 'error'

			# df_all_genes_info[index]['weakness'].append(float(row['Gene weakness (in %)']))

			if args['vcf_annotation_method'].upper() == 'ANNOVAR':
				add_variants_to_dictionnary(df_all_genes_info, index, mutation_type, 'mutation types', row)

			add_variants_to_dictionnary(df_all_genes_info, index, gb_index, 'gb1', row)
			add_variants_to_dictionnary(df_all_genes_info, index, cb_index, 'cb1', row)
			add_variants_to_dictionnary(df_all_genes_info, index, pb_index, 'pb1', row)
			add_variants_to_dictionnary(df_all_genes_info, index, gb1_index, 'gb2', row)
			add_variants_to_dictionnary(df_all_genes_info, index, cb1_index, 'cb2', row)
			add_variants_to_dictionnary(df_all_genes_info, index, pb1_index, 'pb2', row)
			df_all_genes_info[index]['samples'].append(name)

		pbar_files.update(1)
	pbar_files.close()

	dic_intermediate_genes = {}  # gene dictionary to create tsv union file (mean computation for each list)
	dic_patients_genes_lists = {}  # dictionary containing genes list for each sample - use to create the upset plot
	genes_pos_for_chromosomes = {}  # dictionary containing genes position for each chromosome - use to create the chromosomes plot
	genes_pos_for_chromosomes_t1 = {}  # s1 = sample 1
	genes_pos_for_chromosomes_t2 = {}  # s2 = sample 2

	if args['colored_execution'].upper() == 'TRUE' or args['colored_execution'].upper() == 'YES':
		color = '\033[38;2;255;152;0m'
		no_tqdm_bar = False
	else:
		color = ''
		no_tqdm_bar = True

	pbar_files = tqdm(total=len(df_all_genes_info.items()), disable=no_tqdm_bar, bar_format='{l_bar}{bar:30}{r_bar}', ncols=150, smoothing=1)
	pbar_files.set_description(color + f' -> Preparing information for \033[3mgenes_union \033[0m' + color + 'dataframe')

	dic_mutation_types = {}
	dic_mutation_subtypes = {}

	counting_method = args['mutations_types_counting'].upper()
	if counting_method == 'UNIQUE':
		unique = True
	elif counting_method == 'TOTAL':
		total = True

	patients_threshold = args['table_patients_threshold']
	df_names = pd.read_excel(args['dataset_path'])
	ids_column = args['pairs_ids_column_name']
	is_id_column_filled = False
	if ids_column.upper() != 'NONE' and ids_column != 'NO' and ids_column != '' and ids_column != 'FALSE':
		is_id_column_filled = True
		ids = df_names[ids_column].dropna().tolist()
		names1 = df_names[args['time1_column_name']].dropna().tolist()
		names2 = df_names[args['time2_column_name']].dropna().tolist()
		names1 = [name.replace(".vcf", "") for name in names1]
		names2 = [name.replace(".vcf", "") for name in names2]
		dic_names = {}
		j = 0
		for i in range(len(ids)):
			key = f"{names1[j]}___{names2[j]}"
			dic_names[key] = str(ids[j])
			j += 1

	for k, v in df_all_genes_info.items():
		if not v['chr'] in genes_pos_for_chromosomes_t1.keys():
			genes_pos_for_chromosomes_t1[v['chr']] = []
		if not v['chr'] in genes_pos_for_chromosomes_t2.keys():
			genes_pos_for_chromosomes_t2[v['chr']] = []
		genes_pos_for_chromosomes[v['chr']] = [(v['start'], v['end'])] * int(len(v['samples']))

		genes_pos_for_chromosomes_t1[v['chr']] += [(v['start'], v['end'])] * len(v['gb1'])
		genes_pos_for_chromosomes_t2[v['chr']] += [(v['start'], v['end'])] * len(v['gb2'])
		# A CORRIGER : remettre gb1 en dic et si ca change rien en tournant, vérifier positions
		# * int(sum(list(v['gb1'].values())))

		if not k in dic_intermediate_genes.keys():
			dic_intermediate_genes[k] = {}

		in_pairs = []
		for id in dic_genes_in_pairs.keys():
			if k in dic_genes_in_pairs[id]:
				in_pairs.append(id)
		if is_id_column_filled:
			final_in_pairs = []
			for pair in in_pairs:
				final_in_pairs.append(dic_names[pair])
		else:
			final_in_pairs = in_pairs
		dic_intermediate_genes[k]['Samples pairs'] = str(';'.join(final_in_pairs))

		for sample in v['samples']:
			if not sample in dic_patients_genes_lists.keys():  # i
				dic_patients_genes_lists[sample] = []
			dic_patients_genes_lists[sample].append(k)  # add only adds an element to the set if the element is not already present

		dic_intermediate_genes[k]['Chromosome'] = v['chr']
		dic_intermediate_genes[k]['Gene start position'] = v['start']
		dic_intermediate_genes[k]['Gene end position'] = v['end']
		dic_intermediate_genes[k]['Nb samples'] = int(len(v['samples']))

		# unique variants counting for each gene found in all samples from the whole dataset provided
		try:
			if k not in dic_unique_variants.keys():
				dic_unique_variants[k] = {}
				dic_unique_variants[k]['gb1'] = []
			for variant in df_all_genes_info[k]['gb1']:
				if variant not in dic_unique_variants[k]:
					dic_unique_variants[k]['gb1'].append(variant)
		except:
			if k not in dic_unique_variants.keys():
				dic_unique_variants[k] = {}
				dic_unique_variants[k]['gb1'] = []
		try:
			if 'gb2' not in dic_unique_variants[k].keys():
				dic_unique_variants[k]['gb2'] = []
			elif k not in dic_unique_variants.keys():
				dic_unique_variants[k] = {}
				dic_unique_variants[k]['gb2'] = []
			for variant in df_all_genes_info[k]['gb2']:
				if variant not in dic_unique_variants[k]:
					dic_unique_variants[k]['gb2'].append(variant)
		except:
			print(k)
			if k not in dic_unique_variants.keys():
				dic_unique_variants[k] = {}
				dic_unique_variants[k]['gb2'] = []

		# dic_intermediate_genes[k]['Total variants'] = len(set(v['gb1'].keys()).union(set(v['gb2'].keys())))  # gb1 and gb2 are dict with the number of mutation for each gene
		# dic_intermediate_genes[k]['Mean gene weakness'] = np.mean(v['weakness'])
		dic_intermediate_genes[k]['Unique variants'] = len(set(dic_unique_variants[k]['gb1']).union(set(dic_unique_variants[k]['gb2'])))
		dic_intermediate_genes[k]['Total variants'] = len(dic_unique_variants[k]['gb1']) + len(dic_unique_variants[k]['gb2'])
		max_times_variants_number = max(len(set(dic_unique_variants[k]['gb1'])), len(set(dic_unique_variants[k]['gb2'])))
		if not common and not change:
			dic_intermediate_genes[k]['Unique variation (%)'] = round(
				(len(set(dic_unique_variants[k]['gb2'])) - len(set(dic_unique_variants[k]['gb1']))) / max_times_variants_number * 100)
			dic_intermediate_genes[k]['Total variation (%)'] = round((len(dic_unique_variants[k]['gb2']) - len(dic_unique_variants[k]['gb1'])) / max_times_variants_number * 100)
		if not common:
			dic_intermediate_genes[k]['Δ Unique variants'] = len(set(dic_unique_variants[k]['gb2'])) - len(set(dic_unique_variants[k]['gb1']))
			dic_intermediate_genes[k]['Δ total variants'] = len(dic_unique_variants[k]['gb2']) - len(dic_unique_variants[k]['gb1'])
		dic_intermediate_genes[k]['Unique mutations (t1)'] = len(set(dic_unique_variants[k]['gb1']))
		dic_intermediate_genes[k]['Total mutations (t1)'] = len(dic_unique_variants[k]['gb1'])
		dic_intermediate_genes[k]['g.(t1)'] = '|'.join([str(variant) for variant in v['gb1']])
		dic_intermediate_genes[k]['c.(t1)'] = '|'.join([str(variant) for variant in v['cb1']])
		dic_intermediate_genes[k]['p.(t1)'] = '|'.join([str(variant) for variant in v['pb1']])
		dic_intermediate_genes[k]['Unique mutations (t2)'] = len(set(dic_unique_variants[k]['gb2']))
		dic_intermediate_genes[k]['Total mutations (t2)'] = len(dic_unique_variants[k]['gb2'])
		dic_intermediate_genes[k]['g.(t2)'] = '|'.join([str(variant) for variant in v['gb2']])
		dic_intermediate_genes[k]['c.(t2)'] = '|'.join([str(variant) for variant in v['cb2']])
		dic_intermediate_genes[k]['p.(t2)'] = '|'.join([str(variant) for variant in v['pb2']])
		# dic_intermediate_genes[k]['g.TPn union'] = '|'.join([str(k) + '(' + str(v) + ')' for k, v in v['gb1'].items()])
		# dic_intermediate_genes[k]['c.TPn union'] = '|'.join([str(k) + '(' + str(v) + ')' for k, v in v['cb1'].items()])
		# dic_intermediate_genes[k]['p.TPn union'] = '|'.join([str(k) + '(' + str(v) + ')' for k, v in v['pb1'].items()])
		# dic_intermediate_genes[k]['Nb Mutation TPn+1'] = len(set(v['gb2'].keys()))
		# dic_intermediate_genes[k]['g.TPn+1 union'] = '|'.join([str(k) + '(' + str(v) + ')' for k, v in v['gb2'].items()])
		# dic_intermediate_genes[k]['c.TPn+1 union'] = '|'.join([str(k) + '(' + str(v) + ')' for k, v in v['cb2'].items()])
		# dic_intermediate_genes[k]['p.TPn+1 union'] = '|'.join([str(k) + '(' + str(v) + ')' for k, v in v['pb2'].items()])

		pbar_files.update(1)
	pbar_files.close()

	for gene in dic_unique_variants.keys():
		dic_unique_variants[gene]['set'] = set(dic_unique_variants[gene]['gb1']).union(dic_unique_variants[gene]['gb2'])

	print('Creating genes_union.tsv dataframe...')
	# del df_all_genes_info
	df_final_genes = pd.DataFrame(dic_intermediate_genes).T  # get pandas Dataframe from dic_intermediate_genes
	# del dic_intermediate_genes
	df_final_genes.index.name = 'Gene'  # add a name to the index column

	df_final_genes = df_final_genes.sort_values(by=['Nb samples', 'Unique variants', 'Gene'],
												ascending=[False, False, True])  # sort columns
	df_final_genes['Nb samples'] = df_final_genes['Nb samples'].astype('int')
	# if args['vcf_annotation_method'].upper() == 'FUNCOTATOR':
	#     df_final_genes[['Unique variants', 'Total variants', 'g.(t1)', 'c.(t1)', 'p.(t1)',
	#         'Mutations (t1)', 'g.(t2)', 'c.(t2)', 'p.(t2)', 'Mutations (t2)']] = (
	#         np.round(df_final_genes[['Total variants', 'g.(t1)', 'c.(t1)', 'p.(t1)',
	#             'Mutations (t1)', 'g.(t2)', 'c.(t2)', 'p.(t2)', 'Mutations (t2)']], 2))
	# elif args['vcf_annotation_method'].upper() == 'ANNOVAR':
	#     df_final_genes[['Unique variants', 'Total variants', 'g.(t1)', 'c.(t1)', 'p.(t1)',
	#         'Mutations (t1)', 'g.TPn+1 union', 'c.(t2)', 'p.(t2)', 'Mutations (t2)']] = np.round(
	#         df_final_genes[['Unique variants', 'Total variants', 'g.(t1)', 'c.(t1)', 'p.(t1)',
	#             'Mutations (t1)', 'g.(t2)', 'c.(t2)', 'p.(t2)', 'Mutations (t2)']], 2)
	patients_threshold = args['table_patients_threshold']
	if infos:
		# Add additional cancer centric genes information
		df_final_genes = df_final_genes.join(df_databases_info)

		rows_to_remove = []
		for index, row in df_final_genes.iterrows():
			sample_pairs = row['Nb samples']
			if int(sample_pairs) < int(patients_threshold):
				rows_to_remove.append(index)

		# Drop empty informational columns
		df_final_genes = df_final_genes.drop(rows_to_remove)
		empty_cols = [col for col in df_final_genes.columns if df_final_genes[col].isnull().all()]
		df_final_genes.drop(empty_cols, axis=1, inplace=True)

	genes_list = list(df_final_genes.index)

	##########################################
	# Save genes list files (.xlsx and .tsv) #
	logger.info(f'Save genes list files in {Path(output).with_suffix(".MutatedGenes.xlsx")} and {Path(output).with_suffix(".MutatedGenes.tsv")}')

	if 'XLSX' in args['M_MutatedGenes_format(s)'].upper():
		df_final_genes = pd.DataFrame(df_final_genes)
		df_final_genes.to_excel(args['output_path'] + 'merge/genes_union.xlsx')
	elif 'TSV' in args['M_MutatedGenes_format(s)'].upper():
		df_final_genes.to_csv(output, sep='\t')
	elif 'CSV' in args['M_MutatedGenes_format(s)'].upper():
		df_final_genes.to_csv(output, sep=',')

	# Biological process enrichment using the genes list with the ToppGene and Panther API
	if enrichment:
		discard = args['discard_weak_variants'].upper()
		if enrichment and discard == 'TRUE' or discard == 'YES':
			print(f"\033[1m{len(genes_list)}\033[0m{color} genes are concerned.")
			if len(genes_list) < 1500:
				print('Computing ToppGene GOEA analysis...')
				toppgene_name = args['M_ToppGene_name']
				ToppGene_GOEA('summarise', args, genes_list, toppgene_name, logger)
			else:
				print("The VCF is heavy, too many genes are concerned for ToppGene GO to be run.")
			if len(genes_list) < 2000:
				print('Computing Panther GOEA analysis...')
				panther_name = args['M_Panther_name']
				Panther_GOEA('summarise', args, genes_list, panther_name, logger)
			else:
				print("The VCF is too heavy, too many genes are concerned for Panther GO to be run.")

	print('Creating mutations cartography plot...')
	if cytoband_file.upper() != 'NONE' or cytoband_file == '':
		##### Create the Chromosome plot
		create_chromosomes_plot(args, genes_pos_for_chromosomes_t1, genes_pos_for_chromosomes_t2,
								genes_pos_for_chromosomes, cytoband_file, chromosomes_output, step, logger)

	print('Creating upset plot...')
	if dic_patients_genes_lists == {}:
		print('a')
	if upset_output:
		if len(names) < 16:  # Actually upset plot can not be calculated for more than 15 samples
			#####create_upset Create the UpSetPlot
			create_upsetplot(args, dic_patients_genes_lists, category, upset_output, names, logger)


def create_upsetplot(args, data, category, upset_output, names, logger):
	'''
	Create an Upsetplot showing the number of common genes to the different sample set
	'''
	# data = list of genes for each pair
	plt.clf()
	old_names = names
	new_names = set()
	for file_name in names:
		start_index = file_name.find('sons/') + len('sons/')
		end_index = file_name.find('/compared')
		if start_index != -1 and end_index != -1 and start_index < end_index:
			extracted_name = file_name[start_index:end_index].replace('___', '-')
			new_names.add(extracted_name)

	i = 0
	for cat in category:
		j = 0
		for file_name in cat:
			start_index = file_name.find('sons/') + len('sons/')
			end_index = file_name.find('/compared')
			if start_index != -1 and end_index != -1 and start_index < end_index:
				category[i][j] = file_name[start_index:end_index].replace('___', '-')
			j += 1
		i += 1

	mutations_upset_threshold = args['mutations_upset_threshold']
	if mutations_upset_threshold == 0:
		mutations_upset_threshold = max([len(v) for v in data.values()])

	no_gene = old_names - set([k for k in data.keys()])
	if no_gene != set():
		for sample in no_gene:
			data[sample] = set()
			print(f'Warning : sample {sample} don\'t have gene !')
			logger.warning(f'Sample {sample} don\'t have gene !')

	category.sort(key=len)

	all_genes = set()
	for v in data.values():
		if not all_genes:
			all_genes = set(v)
		else:
			all_genes = all_genes.union(set(v))

	data1 = {}
	for key in list(data.keys()):
		start_index = key.find('sons/') + len('sons/')
		element = key.split("/")[-1]
		end_index = key.find(element)
		if start_index != -1 and end_index != -1 and start_index < end_index:
			new_name = key[start_index:end_index].replace('___', '-')
			data1[new_name] = data[key]

	del data
	data = data1

	data2 = []
	cat = []
	to_suppr = set()
	to_suppr2 = set()
	# remove genes that are not found in all samples comparisons

	for sublist in category:
		sublist[0] = sublist[0].rsplit('/', 1)[0]
	updated_data = {}
	for key in data.keys():
		updated_key = key.rstrip('/')
		updated_data[updated_key] = data[key]
	data = updated_data

	cat_count = 0
	for cat in category:
		#print(cat)
		value_count = 0
		for value in cat:
			new_value = value.split('sons/')[1]
			if '___' in new_value:
				new_value = new_value.replace('___', '-')
			if 'genes' in new_value:
				new_value = new_value.split('/genes')[0]
			#print(value + ' -> ' + new_value)
			category[cat_count][value_count] = new_value
			value_count += 1
		cat_count += 1

	for c in category[::-1]:
		save_set = None
		for v in c:
			if save_set is None:
				save_set = set(data[v])
			else:
				save_set = save_set.intersection(set(data[v]))
		data2.append(len(save_set))

	# remove categories that we don't want according to parameters
	category.reverse()
	patients_upset_threshold = int(args['patients_upset_threshold'])
	l_len_cat = []
	for cat in category:
		j = 0
		for x in cat:
			for o in x.split(','):
				j += 1
		l_len_cat.append(j)
	k = 0
	for i in l_len_cat:
		if i < patients_upset_threshold:
			del category[k]
			del data2[k]
			k = k - 1
		k += 1

	mutations_upset_threshold = int(args['mutations_upset_threshold'])
	j = 0
	new_data = []
	new_list = []
	for i in data2:
		if int(i) >= mutations_upset_threshold:
			new_data.append(i)
			new_list.append(category[j])
		j += 1

	input_table = args['dataset_path']
	if 'CSV' in input_table.upper():
		input_table = pd.read_csv(input_table, sep=',', index_col=0)
	elif 'TSV' in input_table.upper():
		input_table = pd.read_csv(input_table, sep='\t', index_col=0)
	elif 'XLSX' in input_table.upper():
		input_table = pd.read_excel(input_table, index_col=0)
	col1_name = args['time1_column_name']
	col2_name = args['time2_column_name']
	col_id_name = args['pair_names_column_name']

	dict_names_ids = {}
	for i in range(len(input_table)):
		try:
			row = input_table.iloc[i]
			key = row[col1_name] + '-' + row[col2_name]
			value = row[col_id_name]
			dict_names_ids[key] = value
		except:
			pass

	for sublist in new_list:
		for i in range(len(sublist)):
			if sublist[i] in dict_names_ids:
				sublist[i] = dict_names_ids[sublist[i]]

	df = from_memberships(new_list, data=new_data)

	# UPSET PLOT
	with warnings.catch_warnings():
		warnings.filterwarnings("ignore")
		UpSet(df, show_counts=True, element_size=40).plot()

	plt.ylabel('Number of mutated genes')
	info = f"min patients in subset: {patients_upset_threshold} | min genes in subset: {mutations_upset_threshold}"
	if len(data2) < 4:
		plt.title('Number of mutated genes ', y=1.05)
		plt.figtext(0.55, 0.05, info, ha="center", fontsize=7,
					bbox={"facecolor": "lightgray", "alpha": 0.5, "pad": 3})  # Adjust fontsize parameter

	if 'SVG' in args['upset_plot_format(s)'].upper():
		upset_output = args['output_path'] + 'merge/' + args['upset_plot_name'] + '.svg'
		if args['verbose_prints'].upper() == 'TRUE':
			print(f'Create UpSetPlot in {upset_output} !')
		logger.info(f'Create UpSetPlot in {upset_output} !')
		try:
			plt.savefig(upset_output, format='svg', dpi=500)
		except:
			plt.savefig(upset_output, format='svg', dpi=200)

	if 'PNG' in args['upset_plot_format(s)'].upper():
		upset_output = args['output_path'] + 'merge/' + args['upset_plot_name'] + '.png'
		if args['verbose_prints'].upper() == 'TRUE':
			print(f'Create UpSetPlot in {upset_output} !')
		logger.info(f'Create UpSetPlot in {upset_output} !')
		try:
			plt.savefig(upset_output, format='png', dpi=500)
		except:
			plt.savefig(upset_output, format='png', dpi=200)
	if 'PDF' in args['upset_plot_format(s)'].upper():
		upset_output = args['output_path'] + 'merge/' + args['upset_plot_name'] + '.pdf'
		if args['verbose_prints'].upper() == 'TRUE':
			print(f'Create UpSetPlot in {upset_output} !')
		logger.info(f'Create UpSetPlot in {upset_output} !')
		try:
			plt.savefig(upset_output, format='pdf', dpi=500)
		except:
			plt.savefig(upset_output, format='pdf', dpi=200)


def main(args, output_path):
	output = args['M_MutatedGenes_name']
	upset_output = args['upset_plot_name']
	enrichment = args['M_enrichment']

	# chromosomes plot variables
	cytoband_file = args['cytoband']
	chromosomes_output = args['chromosomes_plot_name']
	step = args['chromosome-step']
	current_directory = os.getcwd()

	# Logger configuration
	logger = logging.getLogger('LOTUS merge')
	logger.setLevel(logging.DEBUG)
	fh = logging.FileHandler(args['log'])
	fh.setLevel(logging.DEBUG)
	formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
	fh.setFormatter(formatter)
	logger.addHandler(fh)

	if cytoband_file.upper() != 'NONE':
		try:
			logger.info('Verification of {cytoband_file}')
			verif_input(cytoband_file)
			try:
				logger.info('Verification of {chromosomes_output}')
				verif_output(chromosomes_output)
				logger.info('- Output file ok -')
			except ValueError:
				print(f'Problem with {chromosomes_output}: ', sys.exc_info()[0])
				logger.error('- Problem with {chromosomes_output} -')
				exit(1)
			logger.info('- Input file ok -')
		except ValueError:
			print(f'Problem with {cytoband_file}: ', sys.exc_info()[0])
			logger.error('- Problem with {cytoband_file} -')
			exit(1)

	if int(step) < 500000:
		warnings.warn("A step value below 500000 may cause a long calculation time !", RuntimeWarning)
		logger.warning('A step value below 500000 may cause a long calculation time !')

	if args['verbose_prints'].upper() == 'TRUE':
		print('Search for file Lotus_ExternalBases_202301.xlsx ...')
	infos = 'Lotus_ExternalBases_202301.xlsx'
	logger.info('Verification of {infos}')
	try:
		infos = verif_supplementary_information_file(infos, current_directory)
	except ValueError:
		print(f'Problem with {infos}: ', sys.exc_info()[0])
		logger.error('- Problem with {infos} -')
		exit(1)

	logger.info(
		'**************************************************************************************************************')
	logger.info('*** LOTUS merging module ***')
	no_argument = ''
	if enrichment != '' or enrichment.upper() != 'FALSE':
		no_argument += ' --enrichment'
	logger.info('* Start merging *')
	logger.info(f'Current directory : {Path().absolute()}')

	df_dataset = pd.read_excel(args['dataset_path'])
	l_col1_names = df_dataset[args['time1_column_name']].dropna().str.strip().replace('', np.nan).dropna().tolist()
	l_col1_new = [name.split('.funco')[0] for name in l_col1_names]
	l_col2_names = df_dataset[args['time2_column_name']].dropna().str.strip().replace('', np.nan).dropna().tolist()
	l_col2_new = [name.split('.funco')[0] for name in l_col2_names]
	nb_files = len(l_col1_new) + len(l_col2_new)
	names = [col1 + "___" + col2 for col1, col2 in zip(l_col1_new, l_col2_new)]
	path = "config.txt"
	with open(path, 'r') as file:
		for line in file:
			if line.startswith('#') or line == '\n' or line.startswith('------'):
				pass
			else:
				try:
					key, value = line.strip().replace(" ", "").split('=')
				except:
					pass
				if key == "C_MutatedGenes_name":
					C_MutatedGenes_name = value.lower()
				elif key == "C_MutatedGenes":
					C_MutatedGenes_format = value.lower()
				elif key == 'C_variants_selection_approach':
					C_variants_selection_approach = value.lower()
	file_paths = [output_path + 'comparisons/' + name + '/' + C_MutatedGenes_name + '_' + C_variants_selection_approach + "." + C_MutatedGenes_format for name in names]
	file_paths = [path.replace(".vcf", "") for path in file_paths]

	if args['keep_filtered_vcf_after_run'].upper() == 'FALSE' or args['keep_filtered_vcf_after_run'].upper() == 'NO':
		folder_path = args['output_folder_path'] + 'samples/'
		for folder_name in os.listdir(folder_path):
			folder = os.path.join(folder_path, folder_name)
			if os.path.isdir(folder):
				for file_name in os.listdir(folder):
					file = os.path.join(folder, file_name)
					if file.endswith('.vcf') and 'filter' in file_name:
						os.remove(file)

	if upset_output != '' and upset_output.upper() != 'NONE':
		category = [list(i) for i in list(powerset(unique_everseen(file_paths))) if
					i != ()]  # create all possible combinations of samples
	else:
		warnings.warn("Upset plot not computed because of the combination explosion !", RuntimeWarning)
		logger.warning('Upset plot not computed because of the combination explosion !')

	logger.info(f'Merging {nb_files} files !')
	if enrichment.upper() == 'TRUE':
		enrichment = True
	else:
		enrichment = False
	if args['verbose_prints'].upper() == 'TRUE':
		print('Start merging...')

	merge_results(args, file_paths, category, output, upset_output, infos, cytoband_file, chromosomes_output, int(step), nb_files, enrichment, logger)

	logger.info('* End merging *')
	logger.info(
		'**************************************************************************************************************')

# End
