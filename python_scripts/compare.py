#!/usr/bin/env python 3
# -*- coding: utf-8 -*-

import logging
import sys
import os
import pandas as pd
import re
import numpy as np
import warnings
import matplotlib.pyplot as plt
from pathlib import Path
from python_scripts.check_files import verif_input_vcf, verif_output, verif_input_config, verif_input, verif_supplementary_information_file
from python_scripts.toppgene_api import ToppGene_GOEA
from python_scripts.panther_api import Panther_GOEA
from python_scripts.read_vcf import read_vcf, get_vcf_header
from python_scripts.path_modification import true_stem
from python_scripts.read_gff3 import read_gff3


def get_informations_for_genes(info_file, logger):
	df = pd.read_excel(info_file, index_col=1)
	df = df.drop(['Ordre'], axis=1)
	df.set_axis([source.split(' Info')[0] for source in df.columns], axis="columns")
	print(f'Extracting information from {len(list(df.columns))} sources: {", ".join(list(df.columns))}')
	logger.info(f'Extract information from {len(list(df.columns))} sources: {", ".join(list(df.columns))}')
	#print('DF',df)
	return (df)


def get_variants(file, variants_save):
	'''
	Take a vcf file and return all its variants in a set (or a dictionnary with the index fields)
	Input : vcf file and empty set to save variants
	Output : set containing the variants
	'''
	for line in read_vcf(file):
		if type(line) == type({}):
			intfield = line
		else:
			chr = line[intfield['idx_chr']]
			pos = line[intfield['idx_pos']]
			ref = line[intfield['idx_ref']]
			alts = line[intfield['idx_alts']]
			for alt in alts.split(','):
				variants_save.add((chr,pos,ref,alt))
	return variants_save


def add_to_genes_anno(genes : dict, gene : str, type : str, chr : str, gene_position : (), mg : str, mc : str, mp : str, anno_info : str):
	'''
	Add 1 to a gene in the genes dictionary according to its type
        Input : genes dictionary, gene name, type of the gene ('weak' or 'strong') + add chromosome, gene position, mutation in genomic sequence, mutation in coding sequence, mutation in proteic sequence
	'''
	if not gene in genes.keys():
		genes[gene]={}
		genes[gene]['weak']=0
		genes[gene]['strong']=0
		genes[gene]['common']=0
		genes[gene]['mg']=[]
		genes[gene]['mc']=[]
		genes[gene]['mp']=[]
		genes[gene]['anno_info']=[]
	genes[gene][type]+=1
	genes[gene]['chr']=chr
	genes[gene]['gene_pos']=gene_position
	genes[gene]['mg'].append(mg)
	genes[gene]['mc'].append(mc)
	genes[gene]['mp'].append(mp)
	genes[gene]['anno_info'].append(anno_info)


def add_to_genes_funco(genes : dict, gene : str, type : str, chr : str, gene_position : (), mg : str, mc : str, mp : str):
	'''
	Add 1 to a gene in the genes dictionary according to its type
        Input : genes dictionary, gene name, type of the gene ('weak' or 'strong') + add chromosome, gene position, mutation in genomic sequence, mutation in coding sequence, mutation in proteic sequence
	'''
	if not gene in genes.keys():
		genes[gene]={}
		genes[gene]['weak']=0
		genes[gene]['strong']=0
		genes[gene]['mg']=[]
		genes[gene]['mc']=[]
		genes[gene]['mp']=[]
	genes[gene][type]+=1
	genes[gene]['chr']=chr
	genes[gene]['gene_pos']=gene_position
	genes[gene]['mg'].append(mg)
	genes[gene]['mc'].append(mc)
	genes[gene]['mp'].append(mp)


def modify_variants_pass_and_get_genes(dict_para, file1, file2, variants, weak, strong, gene_name_dico : {}, gene_id_dico : {}, transcript_dico : {}, parent_name : str, logger):
	'''
	Get the pass variants lines from a vcf and add the variant type ('weak'/'strong'/'common') + get the count of types weak and strong for every gene 
	Input : vcf file 1 to modify, vcf file 2 (just used for the name), set of specific variants from file 1, set of weak variants from file 1, set of strong variants from file 1 and the logger
	Output : dictionary containing genes and their number of weak and strong variants (specific to vcf file 1)
	'''

	outfile = parent_name + '/' + 'common_variants.passed.vcf'
	genes = {}
	with open(outfile, 'w') as o:
		for line in get_vcf_header(file1):
			o.write(line+'\n')
		for line in read_vcf(file1):
			if type(line) == type({}):
				intfield = line
			else:
				chromosome = line[intfield['idx_chr']]
				ID_ensemble = []
				anno_info = []
				mg = [] #mutation in genomic sequence
				mc = [] #mutation in coding sequence
				mp = [] #mutation in proteic sequence
				infos = [tuple(infos.split('=')) if len(infos.split('='))>1 else (infos,'') for infos in line[intfield['idx_info']].split(';')]
				nb_alt = len(line[intfield['idx_alts']].split(','))

				if dict_para['vcf_annotation_method'].upper() == "ANNOVAR":
					dictionary = {key: value for key, value in infos}
					gene = dictionary['Gene.refGene']
					if '\\' in gene:
						gene = gene.split('\\')[0]
					#print(gene)
					try:
						gene_position = gene_name_dico[gene][0]
					except:
						#print('a')
						pass
					g_mutation = "g." + chromosome + ":" + line[intfield['idx_pos']] + line[intfield['idx_ref']] + ">" + line[intfield['idx_alts']]
					for i in range(nb_alt):
						if nb_alt > 1:
							g_mutation = "g." + chromosome + ":" + line[intfield['idx_pos']] + line[intfield['idx_ref']] + ">" + line[intfield['idx_alts']]
							mutation1 = g_mutation.split(',')[0]
							if '>' in mutation1:
								mutation2 = mutation1[:-1] + mutation1.split('>')[1]
							elif 'ins' in mutation1:
								mutation2 = mutation1.split['ins'][0] + 'ins' + mutation1.split('ins')[1]
							elif 'del' in mutation1:
								mutation2 = mutation1.split['del'][0] + 'del' + mutation1.split('del')[1]
							mg.append(mutation1)
							mg.append(mutation2)
						else:
							mg.append(g_mutation)

						if dictionary['ExonicFunc.refGene'] == ".":
							anno_info.append('NO INFO')
							mc.append('NO INFO')
							mp.append('NO INFO')
						if dictionary['ExonicFunc.refGene'] == "unknown":
							anno_info.append('unknown')
							mc.append('unknown')
							mp.append('unknown')
						else:
							anno_info.append(dictionary['ExonicFunc.refGene'])
						i=0
						for info in anno_info:
							if info == '.':
								anno_info[i] = 'NO INFO'
							i+=1
						i=0
						pattern1 = r'NM(.*?:.*?):'
						matches1 = re.findall(pattern1, dictionary['AAChange.refGene'])
						if len(matches1) == 1:
							dictionary['AAChange.refGene'] + ','
						modified_matches = ['NM' + match for match in matches1]
						pattern2 = r'c\.(.*?):'
						matches2 = re.findall(pattern2, dictionary['AAChange.refGene'])
						k = 0
						pattern3 = r'p\.(.*?)(?=,|$)'
						matches3 = re.findall(pattern3, dictionary['AAChange.refGene'])
						dic_mc = {}
						dic_mp = {}
						list1 = []
						list2 = []

						for transcript in modified_matches:
							if 'NMT:' in transcript:
								transcript = transcript.replace('NMT:', '')
							list1.append(matches2[k])
							dic_mc[transcript] = list1
							list2.append(matches3[k])
							dic_mp[transcript] = list2
							if k == 0: # NOUS MANQUE SUREMENT L AUTRE MOITIE SI NB_ALTS > 1, A VERIFIER
								if k != len(modified_matches) -1:
									mc.append(transcript + ' -> ' + matches2[k] + ' | ')
									mp.append(transcript + ' -> ' + matches3[k] + ' | ')
								elif k == len(modified_matches) - 1:
									mc.append(transcript + ' -> ' + matches2[k])
									mp.append(transcript + ' -> ' + matches3[k])
							elif k >= 1:
								if k != len(modified_matches):
									mc[i] = mc[i] + transcript + ' -> ' + matches2[k] + ' | '
									mp[i] = mp[i] + transcript + ' -> ' + matches3[k] + ' | '
								elif k == len(modified_matches) - 1:
									mc[i] = mc[i] + transcript + ' -> ' + matches2[k]
									mp[i] = mp[i] + transcript + ' -> ' + matches3[k]
							k += 1

				elif dict_para['vcf_annotation_method'].upper() == "FUNCOTATOR":
					idx_funcotation = list(zip(*infos))[0].index('FUNCOTATION')
					info_funcotation = list(zip(*infos))[1][idx_funcotation].split(',')
					gene = info_funcotation[0].split('|')[0].lstrip('[')
					id_ensembl = info_funcotation[0].split('|')[12]
					try:
						if len(gene_name_dico[gene]) > 1:
							gene_position = transcript_dico[id_ensembl.split('.')[0]][0]
						else:
							gene_position = gene_name_dico[gene][0]
					except KeyError:
						try:
							gene_position = transcript_dico[id_ensembl.split('.')[0]][0]
						except KeyError:
							gene_position = ('',chromosome,'')
					for i in range(nb_alt):
						mg.append(info_funcotation[i].split('|')[11])
						mc.append(info_funcotation[i].split('|')[16])
						mp.append(info_funcotation[i].split('|')[18])
						ID_ensemble.append(info_funcotation[i].split('|')[12])

				try:
					if mc[-1].endswith(' | '):
						mc[-1] = mc[-1].rstrip(' | ')
					if mp[-1].endswith(' | '):
						mp[-1] = mp[-1].rstrip(' | ')
				except:
					print('a')

				variant_type = []
				for i, alt in enumerate(line[intfield['idx_alts']].split(',')):							#For each variants
					variant = (line[intfield['idx_chr']], line[intfield['idx_pos']], line[intfield['idx_ref']], alt)
					# Add variant type to the info field

					if dict_para['vcf_annotation_method'].upper() == "ANNOVAR":
						if variant in weak:
							variant_type.append('weak')
							add_to_genes_anno(genes, gene, 'weak', chromosome, gene_position, mg[i], mc[i], mp[i], anno_info[i])
						elif variant in strong:
							variant_type.append('strong')
							add_to_genes_anno(genes, gene, 'strong', chromosome, gene_position, mg[i], mc[i], mp[i], anno_info[i])
						else:
							variant_type.append('common')
							try:
								add_to_genes_anno(genes, gene, 'common', chromosome, gene_position, mg[i], mc[i], mp[i], anno_info[i])
							except:
								print('b')
					if genes[gene]['anno_info'] is list:
						if len(genes[gene]['anno_info']) >= 2:
							print(genes[gene])
					elif dict_para['vcf_annotation_method'].upper() == "FUNCOTATOR":
						variant = (line[intfield['idx_chr']], line[intfield['idx_pos']], line[intfield['idx_ref']], alt)
						# Add variant type to the info field
						if variant in weak:
							variant_type.append('weak')
							add_to_genes_funco(genes, gene, 'weak', chromosome, gene_position, mg[i], mc[i], mp[i])
						elif variant in strong:
							variant_type.append('strong')
							add_to_genes_funco(genes, gene, 'strong', chromosome, gene_position, mg[i], mc[i], mp[i])
						else:
							variant_type.append('common')

				# Recreate and write the vcf line
				line[intfield['idx_info']] = line[intfield['idx_info']]+';VARIANT_TYPE='+','.join(variant_type)		#Add the variant type to the vcf Info field
				o.write("\t".join(line)+'\n')
	return genes


def create_graph_snp(df_snp, df_snp2, outname, logger):
	'''
	Snp count plot creation
	Input : snp profile of file 1, snp profile of file 2, output name for the plot and the logger
	'''

	name1 = true_stem(df_snp.columns[2])
	name2 = true_stem(df_snp2.columns[2])

	# Colors
	white_96 = ['white']*96
	color_6 = ['darkblue', 'blue', 'lightblue', 'darkgreen', 'green', 'lightgreen']
	color_96 = []
	for i in color_6:
		color_96 += [i]*16
	########


	bars = [row[1] for index, row in df_snp.iterrows()]
	height = [float(i) for i in df_snp.iloc[:,2]-df_snp2.iloc[:,2]]
	height2 = [float(i) for i in df_snp.iloc[:,2]]
	height3 = [float(-i) for i in df_snp2.iloc[:,2]]
	height4 = [i if i > 0 else 0 for i in height]
	height5 = [i if i < 0 else 0 for i in height]
	del height
	group = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']
	x_pos = np.arange(len(bars))

	# Create bars in two subplot

	fig = plt.figure(figsize=(15, 10))
	plt.gcf().subplots_adjust(left = 0.1, bottom = 0.1, right = 0.9, top = 0.9, wspace = 0, hspace = 0.12)
	ax1 = plt.subplot(211)
	ax1.bar(x_pos, height2, color=white_96, edgecolor=color_96, linestyle="--")
	ax1.bar(x_pos, height4, color=color_96)
	ax1.set_xticks(x_pos)
	ax1.set_xticklabels(bars, rotation=90, fontsize=9)
	ax2 = plt.subplot(212)
	ax2.bar(x_pos, height3, color=white_96, edgecolor=color_96, linestyle="--")
	ax2.bar(x_pos, height5, color=color_96)

	# x-axis and y-axis creation

	ax2.set_xticks(x_pos)
	ax2.set_xticklabels([])
	ax3 = ax2.twiny()
	newpos = [8, 24, 40, 56, 72, 88]
	ax3.set_xticks(newpos)
	ax3.set_xticklabels(group)
	ax2.xaxis.set_ticks_position('top')
	ax3.xaxis.set_ticks_position('bottom')
	ax3.xaxis.set_label_position('bottom')
	ax3.spines['bottom'].set_position(('outward', 5))

	ax3.set_xlabel('Mutation types')
	ax3.set_xlim(ax2.get_xlim())

	plt.tick_params(
	axis='x',          # changes apply to the x-axis
	which='both',      # both major and minor ticks are affected
	bottom=False,      # ticks along the bottom edge are off
	labelbottom=True)  # labels along the bottom edge are on

	for xtick, color in zip(ax3.get_xticklabels(), color_6):
		xtick.set_color(color)
		xtick.set_size(12)
	ax3.spines['bottom'].set_visible(False)

	plt.draw()		# populate the yticklabels

	ylabels = [round(float(ytick.get_text()), 2) if (ytick.get_text()[0] == '0' or ytick.get_text()[0] == '1') else round(float(ytick.get_text()[1:]),2) for ytick in ax1.get_yticklabels()]
	ylabels2 = [round(float(ytick.get_text()), 2) if (ytick.get_text()[0] == '0' or ytick.get_text()[0] == '1') else round(float(ytick.get_text()[1:]),2) for ytick in ax2.get_yticklabels()]
	max_ylabel = max(ylabels+ylabels2)

	ylab = [x for x in np.arange(0, max_ylabel+0.01, 0.01)]
	ax1.set_yticks(ylab)
	ax1.set_yticklabels([str(round(i*100,0)) for i in  ylab])
	ax2.set_yticks([-i for i in ylab[::-1]])
	ax2.set_yticklabels([str(round(i*100,0)) for i in  ylab[::-1]])
	
	# Sample name
	ax1.annotate(name1, xy=(0, ax1.get_ylim()[1] - 0.08 * ax1.get_ylim()[1]))
	ax2.annotate(name2, xy=(0, ax2.get_ylim()[0] - 0.05 * ax2.get_ylim()[0]))

	# y-axis title
	fig.text(0.05, 0.5, 'Difference between the percentages of each mutation-type', va='center', rotation='vertical')

	logger.info(f'Save profile comparison graph in {outname}')
	plt.savefig(outname) #.svg
	outname = Path(outname).with_suffix('.png')
	logger.info(f'Save profile comparison graph in {outname}')
	plt.savefig(outname) #.png
	
	plt.close()


def create_merged_profile_tsv(dict_para, df1, df2, out_C_SNP_profile):
	dftot = pd.merge(df1, df2)
	dftot = dftot.rename({'Unnamed: 0' : 'Mutation types', 'Unnamed: 1' : 'DNA context'}, axis=1)
	if 'TSV' in dict_para['C_SNP_profile'].upper():
		dftot.to_csv(Path(out_C_SNP_profile).with_suffix(".tsv"), sep='\t', index=False)
	if 'CSV' in dict_para['C_SNP_profile'].upper():
		dftot.to_csv(Path(out_C_SNP_profile).with_suffix(".csv"), sep=',', index=False)
	if 'XLSX' in dict_para['C_SNP_profile'].upper():
		dftot.to_excel(Path(out_C_SNP_profile).with_suffix(".xlsx"), index=False)

def graph_snp(dict_para, snp_profile_files, out_C_SNP_profile, logger):
	'''
	Create the snp profile plot for each snp profile file (if exist)
	Input : snp profile files, snp profile files, output name for the plot and the logger	
	'''
	if dict_para['verbose_prints'].upper() == 'TRUE':
		print('Create profile comparison graph...')
	for num in range(1,len(snp_profile_files)):
		df1 = pd.read_csv(snp_profile_files[num-1], sep='\t')
		df2 = pd.read_csv(snp_profile_files[num], sep='\t')
		create_graph_snp(df1,df2,out_C_SNP_profile, logger)
		create_merged_profile_tsv(dict_para, df1, df2, out_C_SNP_profile)


def create_graph_indel(deletion1, deletion2, insertion1, insertion2, outname, logger):
	'''
	Create the indel count comparison plot
        Input : Insertion count file 1, Insertion count file 2, Deletion count file 1, Deletion count file 2, output name for the plot and the logger
	'''
	
	insert = True
	delet = True
	if (deletion1 == None or deletion2 == None) and (insertion1 == None or insertion2 == None):  # no insertion nor deltion counts file
		print('Warning ! No indel files ! No graphs will be created !')
		logger.warning(f'No indel files ! No graphs will be created !')
		return None
	elif insertion1 == None or insertion2 == None:  # insertion counts file does not exist
		print('Warning ! No insertion file !')
		logger.warning(f'No insertion file !')
		insert = False
	elif deletion1 == None or deletion2 == None:  # deltion counts file does not exist
		print('Warning ! No deletion file !')
		logger.warning(f'No deletion file !')
		delet = False

	# Get dataframe from tsv file

	if insert:
		ins1 = pd.read_csv(insertion1, sep='\t', header=0, index_col=0).iloc[:, 0]
		sum_ins1 = sum(ins1)
		ins2 = pd.read_csv(insertion2, sep='\t', header=0, index_col=0).iloc[:, 0]
		sum_ins2 = sum(ins2)
		height_ins1 = list([float(i) / sum_ins1 for i in ins1])
		height_ins2 = list([-(float(i) / sum_ins2) for i in ins2])
		bars_ins1 = list(ins1.index)
		bars_ins2 = list(ins2.index)
		name_ins1 = ins1.name
		name_ins2 = ins2.name
	if delet:
		del1 = pd.read_csv(deletion1, sep='\t', header=0, index_col=0).iloc[:, 0]
		del2 = pd.read_csv(deletion2, sep='\t', header=0, index_col=0).iloc[:, 0]
		sum_del1 = sum(del1)
		sum_del2 = sum(del2)
		height_del1 = list([float(i) / sum_del1 for i in del1])
		height_del2 = list([-(float(i) / sum_del2) for i in del2])
		bars_del1 = list(del1.index)
		bars_del2 = list(del2.index)
		name_del1 = del1.name
		name_del2 = del2.name

	if delet and insert:
		maximum = max([max(bars_ins1)+1, max(bars_ins2)+1, max(bars_del1)+1, max(bars_del2)+1])
		max_y = max(height_del1+height_del2+height_ins1+height_ins2)+0.05	
	elif delet:
		maximum = max([max(bars_del1)+1, max(bars_del2)+1])
		max_y = max(height_del1+height_del2)+0.05
	elif not delet:
		maximum = max([max(bars_ins1)+1, max(bars_ins2)+1])
		max_y = max(height_ins1+height_ins2)+0.05

	x1 = [0]+[i+1 for i in range(maximum)]
	x2 = [0]+[i+1 for i in range(maximum)]

	width = 0.25

	fig = plt.figure(figsize=(15, 10))
	ax1 = plt.subplot(1,1,1)

	# create plot according to deletion and insertion counts

	if delet and insert:
		plt.bar([float(i)+(width*0.65) for i in bars_ins1], height_ins1, color = 'r', width = width, edgecolor = 'r', label='Insertion_'+true_stem(name_ins1))
		plt.bar([float(i)-(width*0.65) for i in bars_del1], height_del1, color = 'darkred', width = width, edgecolor = 'darkred', label='Deletion_'+true_stem(name_del1))
		plt.bar([float(i)+(width*0.65) for i in bars_ins2], height_ins2, color = 'b', width = width, edgecolor = 'b', label='Insertion_'+true_stem(name_ins2))
		plt.bar([float(i)-(width*0.65) for i in bars_del2], height_del2, color = 'darkblue', width = width, edgecolor = 'darkblue', label='Deletion_'+true_stem(name_del2))
		x1 = sorted(list(set(bars_del1).union(set(bars_ins1))))
		x2 = sorted(list(set(bars_ins2).union(set(bars_del2))))
	elif not insert:
		plt.bar([float(i) for i in bars_del1], height_del1, color = 'darkred', width = width, edgecolor = 'darkred', label='Deletion_'+true_stem(name_del1))
		plt.bar([float(i) for i in bars_del2], height_del2, color = 'darkblue', width = width, edgecolor = 'darkblue', label='Deletion_'+true_stem(name_del2))
		x1 = bars_del1
		x2 = bars_del2
	elif not delet:
		plt.bar([float(i) for i in bars_ins1], height_ins1, color = 'r', width = width, edgecolor = 'r', label='Insertion_'+true_stem(name_ins1))
		plt.bar([float(i) for i in bars_ins2], height_ins2, color = 'b', width = width, edgecolor = 'b', label='Insertion_'+true_stem(name_ins2))
		x1 = bars_ins1
		x2 = bars_ins2

	plt.legend()
	ax1.grid(axis='x', color='lightgrey', linestyle='-', linewidth=0.5)
	max_x = max(x1+x2)

	ax1.set_xlim(0, max_x+1)
	ax1.set_xlabel("Indel size (bp)")
	ax1.set_ylabel("Indel percentage")
	ax1.set_xticks(x2)
	ax1.set_xticklabels(x2, fontsize=10)

	ax2 = ax1.twiny()
	ax2.set_xlim(0, max_x+1)
	ax2.set_xticks(x1)
	ax2.set_xticklabels(x1, fontsize=10)
	ax2.grid(axis='x', color='lightgrey', linestyle='-', linewidth=0.5)

	fig.set_size_inches(max_x/3, 10)


	y_pos = list(np.arange(-1, 1.1, 0.1))
	y_value = [round(i,1) for i in (list(np.arange(1, 0, -0.1))+list(np.arange(0, 1.1, 0.1)))]
	ax1.set_yticks(y_pos)
	ax1.set_yticklabels(y_value)

	ylabels = [str(round(float(ytick.get_text()), 2)) if (ytick.get_text()[0] == '0' or ytick.get_text()[0] == '1') else str(round(float(ytick.get_text()[1:]),2)) for ytick in ax1.get_yticklabels()]
	ax1.set_yticks(ax1.get_yticks())
	ax1.set_yticklabels([str(round(float(i)*100,1)) for i in ylabels])
	ax1.set_ylim(-max_y, max_y)

	plt.hlines(0, 0, maximum+1, color='black', linewidth=1)
	

	logger.info(f'Draw indel size barplot in {outname}')

	plt.savefig(outname)

	outname = Path(outname).with_suffix('.png')
	plt.savefig(outname)

	plt.close()


def graph_indel(dict_para, insertions_count_files, deletions_count_files, out_indel, logger):
	'''
	Create the indel count comparison plot for each indel count file (if exist) 
	Input : Insertion count files, Deletion count files, output name for the plot and the logger 
	'''
	if len(insertions_count_files) != len(deletions_count_files):
		logger.error(f'Not the same number of insertions and deletions count files. No graph will be produced !')
		raise ValueError(f'Not the same number of insertions and deletions count files. No graph will be produced !')
		exit(1)
	if dict_para['verbose_prints'].upper() == 'TRUE':
		print('Draw indel size barplot...')
	for num in range(1, len(insertions_count_files)):
		create_graph_indel(insertions_count_files[num-1], insertions_count_files[num], deletions_count_files[num-1], deletions_count_files[num], out_indel, logger)

def add_mutation_type(dic_mutation_types, dic_mutation_subtypes, genes, gene):
	# Mutation type counting
	for mutation in genes[gene]['mg']:
		ref = mutation.split('>')[0].split(':')[1]
		ref = re.sub(r'[^A-Za-z]', '', ref)
		alt = mutation.split('>')[1]
		if len(ref) == len(alt) == 1:
			dic_mutation_types['SNP'] += 1
		elif len(ref) == len(alt) == 2:
			dic_mutation_types['DNP'] += 1
		elif len(ref) == len(alt) == 3:
			dic_mutation_types['TNP'] += 1
		elif len(ref) == len(alt) > 3:
			dic_mutation_types['ONP'] += 1
		elif len(ref) > len(alt):
			dic_mutation_types['DELETION'] += 1
			dic_mutation_types['INDEL'] += 1
		elif len(ref) < len(alt):
			dic_mutation_types['INSERTION'] += 1
			dic_mutation_types['INDEL'] += 1

	# Mutation subtype counting
	for mutation in genes[gene]['anno_info']:
		if mutation == 'synonymous_SNV':
			dic_mutation_subtypes['synonymous SNV'] += 1
		elif mutation == 'nonsynonymous_SNV':
			dic_mutation_subtypes['nonsynonymous SNV'] += 1
		elif mutation == 'frameshift_substitution':
			dic_mutation_subtypes['frameshift substitution'] += 1
		elif mutation == 'nonframeshift_substitution':
			dic_mutation_subtypes['non frameshift substitution'] += 1
		elif mutation == 'stopgain':
			dic_mutation_subtypes['stopgain'] += 1
		elif mutation == 'stoploss':
			dic_mutation_subtypes['stoploss'] += 1
		elif mutation == 'startloss':
			dic_mutation_subtypes['startloss'] += 1
		elif mutation == 'frameshift_insertion':
			dic_mutation_subtypes['frameshift insertion'] += 1
		elif mutation == 'nonframeshift_insertion':
			dic_mutation_subtypes['non frameshift insertion'] += 1
		elif mutation == 'frameshift_deletion':
			dic_mutation_subtypes['frameshift deletion'] += 1
		elif mutation == 'nonframeshift_deletion':
			dic_mutation_subtypes['non frameshift deletion'] += 1
		elif mutation == 'unknown':
			dic_mutation_subtypes['unknown'] += 1

	return dic_mutation_types, dic_mutation_subtypes

def create_mutation_types_barplot(dict_para, dic_mutation_types_t1, dic_mutation_types_t2):
	plt.clf()

	labels = [label for label in dic_mutation_types_t1 if dic_mutation_types_t1[label] != 0 and dic_mutation_types_t2[label] != 0]
	values_t1 = [dic_mutation_types_t1[label] for label in labels]
	values_t2 = [dic_mutation_types_t2[label] for label in labels]
	x = np.arange(len(labels))

	bar_width = 0.3
	bar_gap = 0.1

	fig, ax = plt.subplots(figsize=(12, 8))

	rects1 = ax.bar(x - bar_width - bar_gap / 2, values_t1, bar_width, label='time 1', color=(0, 0.5, 0), edgecolor='black', linewidth=1)
	rects2 = ax.bar(x + bar_gap / 2, values_t2, bar_width, label='time 2', color=(0, 0, 0.5), edgecolor='black', linewidth=1)

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
	patient = dict_para['output_path_comparison'].split('isons/')[1].replace("/", "").replace("___", "_")
	ax.set_title('Mutation types comparison between time 1 and time 2\n(patient ' + patient + ')', fontsize=16, pad=10)
	ax.set_xticks(x)
	ax.set_xticklabels(labels, ha='right', fontsize=15)
	ax.set_ylim(0, max(max(values_t1), max(values_t2)) * 1.1)

	offset = bar_width / 2
	ax.set_xticks(x - offset, minor=False)
	ax.set_xticklabels(labels, minor=False)
	plt.tick_params(bottom=False)

	with warnings.catch_warnings():
		warnings.simplefilter("ignore")
		ax.set_yticklabels(ax.get_yticks().astype(int),fontsize=12)
	plt.tick_params(bottom=False)

	ax.legend(fontsize=14, edgecolor='white')
	plt.tight_layout()

	file_formats = dict_para['M_types_barplot_format(s)'].upper()
	path = dict_para['output_path_comparison'] + dict_para['M_types_barplot_name']
	if 'PNG' in file_formats:
		plt.savefig(path + '.png', dpi=400)
	if 'PDF' in file_formats:
		plt.savefig(path + '.pdf', dpi=400)
	if 'SVG' in file_formats:
		plt.savefig(path + '.svg', dpi=400)
	if 'JPG' in file_formats:
		plt.savefig(path + '.jpg', dpi=400)
def create_mutation_subtypes_barplot(dict_para, dic_mutation_subtypes_t1, dic_mutation_subtypes_t2):
	plt.clf()

	labels = [label for label in dic_mutation_subtypes_t1 if dic_mutation_subtypes_t1[label] != 0 and dic_mutation_subtypes_t2[label] != 0]
	values_t1 = [dic_mutation_subtypes_t1[label] for label in labels]
	values_t2 = [dic_mutation_subtypes_t2[label] for label in labels]

	x = np.arange(len(labels))
	bar_width = 0.3
	bar_gap = 0.1

	fig, ax = plt.subplots(figsize=(16, 9))

	rects1 = ax.bar(x - bar_width - bar_gap / 2, values_t1, bar_width, label='time 1', color=(0, 0.5, 0), edgecolor='black', linewidth=1)
	rects2 = ax.bar(x + bar_gap / 2, values_t2, bar_width, label='time 2', color=(0, 0, 0.5), edgecolor='black', linewidth=1)

	for rect in rects1:
		height = rect.get_height()
		ax.annotate('{}'.format(height),
					xy=(rect.get_x() + rect.get_width() / 2, height),
					xytext=(0, 3),
					textcoords="offset points",
					ha='center', va='bottom',
					fontsize=14)

	for rect in rects2:
		height = rect.get_height()
		ax.annotate('{}'.format(height),
					xy=(rect.get_x() + rect.get_width() / 2, height),
					xytext=(0, 3),
					textcoords="offset points",
					ha='center', va='bottom',
					fontsize=14)

	ax.set_ylabel('Count', fontsize=14.5, labelpad=10)
	patient = dict_para['output_path_comparison'].split('isons/')[1].replace("/", "").replace("___", "_")
	ax.set_title('Mutation subtypes comparison between time 1 and time 2 \n(patient ' + patient + ')', fontsize=16, pad=10)
	ax.set_xticks(x)
	ax.set_xticklabels(labels, rotation=90, ha='center', fontsize=14)
	ax.set_ylim(0, max(max(values_t1), max(values_t2)) * 1.1)

	offset = bar_width / 2
	ax.set_xticks(x - offset, minor=False)
	ax.set_xticklabels(labels, minor=False)
	plt.tick_params(bottom=False)

	ax.legend(fontsize=14, edgecolor='white')
	plt.tight_layout()

	with warnings.catch_warnings():
		warnings.simplefilter("ignore")
		ax.set_yticklabels(ax.get_yticks().astype(int), fontsize=12)
	plt.tick_params(bottom=False)

	plt.subplots_adjust(bottom=0.3)

	file_formats = dict_para['C_subtypes_barplot_format(s)'].upper()
	path = dict_para['output_path_comparison'] + dict_para['C_subtypes_barplot_name']
	if 'PNG' in file_formats:
		plt.savefig(path + '.png', dpi=400)
	if 'PDF' in file_formats:
		plt.savefig(path + '.pdf', dpi=400)
	if 'SVG' in file_formats:
		plt.savefig(path + '.svg', dpi=400)
	if 'JPG' in file_formats:
		plt.savefig(path + '.jpg', dpi=400)

def compare_vcf(out_gene, dict_para, gene_name_dico : {}, gene_id_dico : {}, transcript_dico : {}, infos : str, enrichment : bool, logger):
	'''
	Compare each vcf files with the n-1 vcf file (starting from the second file)
	Input : vcf files (pass and filtered from the filter module), boolean to know if GOEA need to be done and the logger
	'''
	vcf_pass_files = []
	vcf_pass_files.append(dict_para['passed1'])
	vcf_pass_files.append(dict_para['passed2'])
	vcf_filtered_files = []
	vcf_filtered_files.append(dict_para['filtered1'])
	vcf_filtered_files.append(dict_para['filtered2'])
	if dict_para['vcf_annotation_method'].upper() == "ANNOVAR":
		annotation = 'ANNOVAR'
	elif dict_para['vcf_annotation_method'].upper() == "FUNCOTATOR":
		annotation = 'FUNCOTATOR'

	# Get genes infos if not None
	if infos:
		infos_df = get_informations_for_genes(infos, logger)

	print('Creating comparison genes mutations table...')

	if dict_para['verbose_prints'].upper() == 'TRUE':
		print('Save genes lists from comparisons...')

	color = '\033[38;2;255;215;0m'

	for num in range(1, len(vcf_pass_files)):							# Comparison file by file starting from the second with is n-1
		logger.info(f'Start comparing {true_stem(vcf_pass_files[num-1])} and {true_stem(vcf_pass_files[num])} !')
		variants_pass_save = set()
		variants_pass_save2 = set()
		variants_filtered_save = set()
		variants_filtered_save2 = set()
		variants_pass_save = get_variants(vcf_pass_files[num-1], variants_pass_save)
		variants_pass_save2 = get_variants(vcf_pass_files[num], variants_pass_save2)
		variants_filtered_save = get_variants(vcf_filtered_files[num-1], variants_filtered_save)
		variants_filtered_save2 = get_variants(vcf_filtered_files[num], variants_filtered_save2)

		# Get variants sets specific to each file and separate this sets in weak and strong variants using the filtered file from the other vcf
		croissant1 = variants_pass_save-variants_pass_save2				# variants from pass 1 without pass 2
		croissant2 = variants_pass_save2-variants_pass_save				# variants from pass 2 without pass 1
		strong1 = croissant1-variants_filtered_save2					# variants from pass 1 without pass 2 and filter 2
		strong2 = croissant2-variants_filtered_save						# variants from pass 2 without pass 1 and filter 1
		weak1 = croissant1.intersection(variants_filtered_save2)		# variants from pass 1 without pass 2 but present in filter 2
		weak2 = croissant2.intersection(variants_filtered_save)			# variants from pass 2 without pass 1 but present in filter 1

		# Creation of the dataframe containing impacted gene and creation of the new vcf files containing the genes type ('strong'/'weak'/'common')
		dfGenes = pd.DataFrame(columns=['Gene', 'Chromosome', 'Gene start position', 'Gene end position', 'Unique variants',  'Total variants', 'Variation (%)',
										'|Variants delta|', 'Mutations (t1)', 'g.(t1)', 'c.(t1)', 'p.(t1)', 'Mutations (t2)', 'g.(time2 union)',
										'c.(time2 union)', 'p.(time2 union)'])
		chromosomes_list = []
		gene_position_list = []
		weakness_list = []
		charge_list = []
		dic_unique_variants = {}
		dic_variants_delta = {}
		tumor1_list = []
		tumor2_list = []
		genomic_variant_annotation_list = []
		coding_variant_annotation_list = []
		proteomic_variant_annotation_list = []
		genomic_variant2_annotation_list = []
		coding_variant2_annotation_list = []
		proteomic_variant2_annotation_list = []
		dic_anno_info = {}
		genes1 = modify_variants_pass_and_get_genes(dict_para, vcf_pass_files[num-1], vcf_pass_files[num], croissant1, weak1, strong1, gene_name_dico, gene_id_dico, transcript_dico, str(Path(out_gene).resolve().parent), logger)
		genes2 = modify_variants_pass_and_get_genes(dict_para, vcf_pass_files[num], vcf_pass_files[num-1], croissant2, weak2, strong2, gene_name_dico, gene_id_dico, transcript_dico, str(Path(out_gene).resolve().parent), logger)
		genes_list = list(set(genes1.keys()).union(set(genes2.keys())))

		dic_mutation_types_t1 = {'SNP': 0, 'DNP': 0, 'TNP': 0, 'ONP': 0, 'INDEL': 0, 'INSERTION': 0, 'DELETION': 0}
		dic_mutation_types_t2 = {'SNP': 0, 'DNP': 0, 'TNP': 0, 'ONP': 0, 'INDEL': 0, 'INSERTION': 0, 'DELETION': 0}
		anno = False
		if dict_para['vcf_annotation_method'].upper() == 'ANNOVAR':
			anno = True
			dic_mutation_subtypes_t1 = {'synonymous SNV': 0, 'nonsynonymous SNV': 0, 'frameshift substitution': 0,
									 'non frameshift substitution': 0, 'stopgain': 0, 'stoploss': 0,
									 'startloss': 0, 'frameshift insertion': 0,
									 'non frameshift insertion': 0, 'frameshift deletion': 0,
									 'non frameshift deletion': 0, 'unknown': 0}
			dic_mutation_subtypes_t2 = {'synonymous SNV': 0, 'nonsynonymous SNV': 0, 'frameshift substitution': 0,
									  'non frameshift substitution': 0, 'stopgain': 0, 'stoploss': 0,
									  'startloss': 0, 'frameshift insertion': 0,
									  'non frameshift insertion': 0, 'frameshift deletion': 0,
									  'non frameshift deletion': 0, 'unknown': 0}

		# unique variants counting for each gene found in at least one of the two samples of current pair
		for gene in genes1:
			# Mutation types and subtypes counting
			if anno:
				dic_mutation_types_t1, dic_mutation_subtypes_t1 = add_mutation_type(dic_mutation_types_t1, dic_mutation_subtypes_t1, genes1, gene)

			try:
				if gene not in dic_unique_variants.keys():
					dic_unique_variants[gene] = {}
					dic_unique_variants[gene]['t1'] = []
				for variant in genes1[gene]['mg']:
					if variant not in dic_unique_variants[gene]['t1']:
						dic_unique_variants[gene]['t1'].append(variant)
			except:
				print('a')
		for gene in genes2:
			if anno:
				dic_mutation_types_t2, dic_mutation_subtypes_t2 = add_mutation_type(dic_mutation_types_t2, dic_mutation_subtypes_t2, genes2, gene)
			try:
				if gene not in dic_unique_variants.keys():
					dic_unique_variants[gene] = {}
					dic_unique_variants[gene]['t2'] = []
				if 't2' not in dic_unique_variants[gene].keys():
					dic_unique_variants[gene]['t2'] = []
				for variant in genes2[gene]['mg']:
					if variant not in dic_unique_variants[gene]['t2']:
						dic_unique_variants[gene]['t2'].append(variant)
			except:
				print('b')


		for gene in genes_list:

			# |Variants delta| counting for each gene found in at least one of the two samples of current pair
			if gene in genes1.keys() and gene in genes2.keys():
				unique_set = set(genes1[gene]['mg']).union(set(genes2[gene]['mg']))
				unique_set1 = set(genes1[gene]['mg'])
				unique_set2 = set(genes2[gene]['mg'])
				dic_variants_delta[gene] = len(unique_set2) - len(unique_set1)
				dic_unique_variants[gene]['set'] = unique_set

			elif gene in genes1.keys() and not gene in genes2.keys():
				try:
					dic_variants_delta[gene] = -len(genes1[gene]['mg'])
					dic_unique_variants[gene]['set'] = set(genes1[gene]['mg'])
				except:
					print('a')
			elif gene in genes2.keys() and not gene in genes1.keys():
				try:
					dic_variants_delta[gene] = len(genes2[gene]['mg'])
					dic_unique_variants[gene]['set'] = set(genes2[gene]['mg'])
				except:
					print('a')

			if annotation == 'ANNOVAR':
				if gene in genes1.keys() and gene in genes2.keys():
					chromosomes_list.append(genes1[gene]['chr'])
					gene_position_list.append(genes1[gene]['gene_pos'])
					try:
						weakness_list.append(round(((genes1[gene]['weak']+genes2[gene]['weak'])/(genes1[gene]['weak']+genes2[gene]['weak']+genes1[gene]['strong']+genes2[gene]['strong']))*100,2))
					except:
						weakness_list.append(666)

					charge_list.append(genes1[gene]['weak'] + genes2[gene]['weak'] + genes1[gene]['strong'] + genes2[gene][
							'strong'] + genes1[gene]['common'] + genes2[gene]['common'])
					tumor1_list.append(genes1[gene]['weak'] + genes1[gene]['strong'] + genes1[gene]['common'])
					tumor2_list.append(genes2[gene]['weak'] + genes2[gene]['strong'] + genes2[gene]['common'])
					dic_anno_info[gene] = []
					try:
						# if len(genes1[gene]['anno_info']) > 15:
						# 	print(genes1[gene]['anno_info'])
						for info in genes1[gene]['anno_info']:
							dic_anno_info[gene].append(info)
						for info in genes2[gene]['anno_info']:
							dic_anno_info[gene].append(info)
					except:
						print('error anno info')

					genomic_variant_annotation_list.append('|'.join(genes1[gene]['mg']))
					coding_variant_annotation_list.append('|'.join(genes1[gene]['mc']))
					proteomic_variant_annotation_list.append('|'.join(genes1[gene]['mp']))
					genomic_variant2_annotation_list.append('|'.join(genes2[gene]['mg']))
					coding_variant2_annotation_list.append('|'.join(genes2[gene]['mc']))
					proteomic_variant2_annotation_list.append('|'.join(genes2[gene]['mp']))

				elif gene in genes1.keys():
					chromosomes_list.append(genes1[gene]['chr'])
					gene_position_list.append(genes1[gene]['gene_pos'])
					weakness_list.append(round((genes1[gene]['weak']/(genes1[gene]['weak']+genes1[gene]['strong']))*100,2))
					charge_list.append(genes1[gene]['weak']+genes1[gene]['strong']+genes1[gene]['common'])
					tumor1_list.append(genes1[gene]['weak'] + genes1[gene]['strong'] + genes1[gene]['common'])
					dic_anno_info[gene] = []
					for info in genes1[gene]['anno_info']:
						dic_anno_info[gene].append(info)
					tumor2_list.append(0)
					genomic_variant_annotation_list.append('|'.join(genes1[gene]['mg']))
					coding_variant_annotation_list.append('|'.join(genes1[gene]['mc']))
					proteomic_variant_annotation_list.append('|'.join(genes1[gene]['mp']))
					genomic_variant2_annotation_list.append('')
					coding_variant2_annotation_list.append('')
					proteomic_variant2_annotation_list.append('')

				elif gene in genes2.keys():
					chromosomes_list.append(genes2[gene]['chr'])
					gene_position_list.append(genes2[gene]['gene_pos'])
					weakness_list.append(round((genes2[gene]['weak']/(genes2[gene]['weak']+genes2[gene]['strong']))*100,2))
					tumor1_list.append(0)
					dic_anno_info[gene] = []
					for info in genes2[gene]['anno_info']:
						dic_anno_info[gene].append(info)
					charge_list.append(genes2[gene]['weak'] + genes2[gene]['strong'] + genes2[gene]['common'])
					tumor2_list.append(genes2[gene]['weak'] + genes2[gene]['strong'] + genes2[gene]['common'])
					genomic_variant_annotation_list.append('')
					coding_variant_annotation_list.append('')
					proteomic_variant_annotation_list.append('')
					genomic_variant2_annotation_list.append('|'.join(genes2[gene]['mg']))
					coding_variant2_annotation_list.append('|'.join(genes2[gene]['mc']))
					proteomic_variant2_annotation_list.append('|'.join(genes2[gene]['mp']))
			elif annotation == 'FUNCOTATOR':
				if gene in genes1.keys() and gene in genes2.keys():
					chromosomes_list.append(genes1[gene]['chr'])
					gene_position_list.append(genes1[gene]['gene_pos'])
					try :
						weakness_list.append(round(((genes1[gene]['weak']+genes2[gene]['weak'])/(genes1[gene]['weak']+genes2[gene]['weak']+genes1[gene]['strong']+genes2[gene]['strong']))*100,2))
					except:
						weakness_list.append(666)
					charge_list.append(genes1[gene]['weak']+genes2[gene]['weak']+genes1[gene]['strong']+genes2[gene]['strong'])
					tumor1_list.append(genes1[gene]['weak']+genes1[gene]['strong'])
					tumor2_list.append(genes2[gene]['weak']+genes2[gene]['strong'])
					genomic_variant_annotation_list.append('|'.join(genes1[gene]['mg']))
					coding_variant_annotation_list.append('|'.join(genes1[gene]['mc']))
					proteomic_variant_annotation_list.append('|'.join(genes1[gene]['mp']))
					genomic_variant2_annotation_list.append('|'.join(genes2[gene]['mg']))
					coding_variant2_annotation_list.append('|'.join(genes2[gene]['mc']))
					proteomic_variant2_annotation_list.append('|'.join(genes2[gene]['mp']))
				elif gene in genes1.keys():
					chromosomes_list.append(genes1[gene]['chr'])
					gene_position_list.append(genes1[gene]['gene_pos'])
					weakness_list.append(round((genes1[gene]['weak']/(genes1[gene]['weak']+genes1[gene]['strong']))*100,2))
					charge_list.append(genes1[gene]['weak']+genes1[gene]['strong'])
					tumor1_list.append(genes1[gene]['weak']+genes1[gene]['strong'])
					tumor2_list.append(0)
					genomic_variant_annotation_list.append('|'.join(genes1[gene]['mg']))
					coding_variant_annotation_list.append('|'.join(genes1[gene]['mc']))
					proteomic_variant_annotation_list.append('|'.join(genes1[gene]['mp']))
					genomic_variant2_annotation_list.append('')
					coding_variant2_annotation_list.append('')
					proteomic_variant2_annotation_list.append('')
				elif gene in genes2.keys():
					chromosomes_list.append(genes2[gene]['chr'])
					gene_position_list.append(genes2[gene]['gene_pos'])
					weakness_list.append(round((genes2[gene]['weak']/(genes2[gene]['weak']+genes2[gene]['strong']))*100,2))
					charge_list.append(genes2[gene]['weak']+genes2[gene]['strong'])
					tumor1_list.append(0)
					genomic_variant_annotation_list.append('')
					coding_variant_annotation_list.append('')
					proteomic_variant_annotation_list.append('')
					tumor2_list.append(genes2[gene]['weak']+genes2[gene]['strong'])
					genomic_variant2_annotation_list.append('|'.join(genes2[gene]['mg']))
					coding_variant2_annotation_list.append('|'.join(genes2[gene]['mc']))
					proteomic_variant2_annotation_list.append('|'.join(genes2[gene]['mp']))

		dfGenes['Gene'] = genes_list
		dfGenes['Chromosome'] = chromosomes_list
		# some gene names are empty, error out of range -> try/except to replace empty gene name by empty string
		for i in range(len(gene_position_list)):
			if not all(gene_position_list[i]):
				gene_position_list[i] = tuple('empty' if not element else element for element in gene_position_list[i])
		gene_position_start_list = [element[2][0] for element in gene_position_list]
		dfGenes['Gene start position'] = gene_position_start_list
		gene_position_end_list = [element[2][1] for element in gene_position_list]
		dfGenes['Gene end position'] = gene_position_end_list
		for gene in dic_unique_variants.keys():
			try:
				len(dic_unique_variants[gene]['t1'])
			except:
				dic_unique_variants[gene]['t1'] = []
			try:
				len(dic_unique_variants[gene]['t2'])
			except:
				dic_unique_variants[gene]['t2'] = []
			dfGenes.loc[dfGenes['Gene'] == gene, 'Unique variants'] = len(dic_unique_variants[gene]['set'])
			if dic_variants_delta[gene] == 0:
				dfGenes.loc[dfGenes['Gene'] == gene, '|Variants delta|'] = 0
			else:
				dfGenes.loc[dfGenes['Gene'] == gene, '|Variants delta|'] = abs(dic_variants_delta[gene])  # | variants delta | (helps sorting dataframe)
			if gene in genes1.keys() and gene in genes2.keys():
				dfGenes.loc[dfGenes['Gene'] == gene, 'Total variants'] = len(genes1[gene]['mg']) + len(genes2[gene]['mg'])
				gene_times_max_variants = max(len(genes2[gene]['mg']), len(genes1[gene]['mg']))
				dfGenes.loc[dfGenes['Gene'] == gene, 'Variation (%)'] = round((len(genes2[gene]['mg']) - len(genes1[gene]['mg'])) / gene_times_max_variants * 100)
			elif gene in genes1.keys():
				dfGenes.loc[dfGenes['Gene'] == gene, 'Total variants'] = len(genes1[gene]['mg'])
				dfGenes.loc[dfGenes['Gene'] == gene, 'Variation (%)'] = -100
			elif gene in genes2.keys():
				dfGenes.loc[dfGenes['Gene'] == gene, 'Total variants'] = len(genes2[gene]['mg'])
				dfGenes.loc[dfGenes['Gene'] == gene, 'Variat ion (%)'] = 100
			try:
				dfGenes.loc[dfGenes['Gene'] == gene, 'Mutations (t1)'] = len(genes1[gene]['mg'])
			except:
				dfGenes.loc[dfGenes['Gene'] == gene, 'Mutations (t1)'] = 0
			try:
				dfGenes.loc[dfGenes['Gene'] == gene, 'Mutations (t2)'] = len(genes2[gene]['mg'])
			except:
				dfGenes.loc[dfGenes['Gene'] == gene, 'Mutations (t2)'] = 0
		# dfGenes['Total variants'] = charge_list
		# dfGenes['Gene weakness (in %)'] = weakness_list
		i = 0
		# if annotation == 'ANNOVAR':
		# 	for index, row in dfGenes.iterrows():
		# 		gene = row['Gene']
		# 		dfGenes['Mutation Types'] = dfGenes['Mutation Types'].astype(str)
		# 		for i, info in enumerate(dic_anno_info[gene]):
		# 			dfGenes.at[index, 'Mutation Types'] += info
		# 			if i < len(dic_anno_info[gene]) - 1:
		# 				dfGenes.at[index, 'Mutation Types'] += '|'
		# 		if 'nan' in dfGenes.at[index, 'Mutation Types']:
		# 			dfGenes.at[index, 'Mutation Types'] = dfGenes.at[index, 'Mutation Types'].replace('nan', '')

		# genomic_variant_annotation_list = [x for x in genomic_variant_annotation_list if x]
		# coding_variant_annotation_list = [x for x in coding_variant_annotation_list if x]
		# proteomic_variant_annotation_list = [x for x in proteomic_variant_annotation_list if x]
		# genomic_variant2_annotation_list = [x for x in genomic_variant2_annotation_list if x]
		# coding_variant2_annotation_list = [x for x in coding_variant2_annotation_list if x]
		# proteomic_variant2_annotation_list = [x for x in proteomic_variant2_annotation_list if x]

		dfGenes['g.(t1)'] = genomic_variant_annotation_list
		dfGenes['c.(t1)'] = coding_variant_annotation_list
		dfGenes['p.(t1)'] = proteomic_variant_annotation_list
		dfGenes['g.(time2 union)'] = genomic_variant2_annotation_list
		dfGenes['c.(time2 union)'] = coding_variant2_annotation_list
		dfGenes['p.(time2 union)'] = proteomic_variant2_annotation_list

		dfGenes = dfGenes.sort_values(by=['|Variants delta|', 'Unique variants', 'Variation (%)', 'Gene'], ascending=[False, False, False, True])
		dfGenes = dfGenes[(dfGenes['Gene'] != 'NONE') & (dfGenes['Gene'] != 'Unknown')] #skipping variants with no corresponding gene name by the annotator
		#dfGenes = dfGenes.set_index('Gene')

		# Save impacted genes information
		if infos:
			dfGenes = dfGenes.join(infos_df)
			# Drop empty informational columns
			dfGenes = dfGenes.dropna(axis=1, how='all')
			#dfGenes.to_excel('comparison_test.xlsx', index=False)

		logger.info(f'Save genes list in {out_gene}')
		if dict_para['verbose_prints'].upper() == 'TRUE':
			print(f'Save genes list in {out_gene}')
		if 'TSV' in dict_para['C_MutatedGenes'].upper():
			dfGenes.to_csv(out_gene, sep='\t')
		elif 'CSV' in dict_para['C_MutatedGenes'].upper():
			dfGenes.to_csv(out_gene, sep=',')
		elif 'XLSX' in dict_para['C_MutatedGenes'].upper():
			dfGenes.to_excel(out_gene, index=False)
		logger.info(f'Save genes list in {Path(out_gene).with_suffix(".xlsx")}')
		if dict_para['verbose_prints'].upper() == 'TRUE':
			print(f'Save genes list in {Path(out_gene).with_suffix(".xlsx")}')
		#dfGenes.to_excel(Path(out_gene).with_suffix('.xlsx'))


		print('Creating mutation types and subtypes plots...')
		if not 'NO' in dict_para['C_types_barplot_name'] and not 'FALSE' in dict_para['C_types_barplot_name']:
			create_mutation_types_barplot(dict_para, dic_mutation_types_t1, dic_mutation_types_t2)
		if not 'NO' in dict_para['C_subtypes_barplot_name'] and not 'FALSE' in dict_para['C_subtypes_barplot_name']:
			create_mutation_subtypes_barplot(dict_para, dic_mutation_subtypes_t1, dic_mutation_subtypes_t2)

		################################################
		# Biological process enrichment using the genes list with the ToppGene and Panther API
		if enrichment:
			discard = dict_para['discard_weak_variants'].upper()
			if enrichment.upper() == 'TRUE' and discard == 'TRUE' or discard == 'YES':
				print(f"\033[1m{len(genes_list)}\033[0m\033[38;2;255;193;7m genes are concerned.")
				if len(genes_list) < 1500:
					print('Computing ToppGene GOEA analysis...')
					toppgene_name = dict_para['C_ToppGene_name']
					ToppGene_GOEA('summarise', dict_para, genes_list, toppgene_name, logger)
				else:
					print("The VCF is heavy, too many genes are concerned for ToppGene GO to be run.")
				if len(genes_list) < 2000:
					print('Computing Panther GOEA analysis...')
					panther_name = dict_para['C_Panther_name']
					Panther_GOEA('summarise', dict_para, genes_list, panther_name, logger)
				else:
					print("The VCF is heavy, too many genes are concerned for Panther GO to be run.")

def main(args):
	
	output_path = args['output_path_comparison']
	if args['C_enrichment'].upper() == 'FALSE':
		enrichment = False
	else :
		enrichment = args['C_enrichment']
	out_gene = output_path + args['C_MutatedGenes_name'] + "." + args['C_MutatedGenes'].lower()
	if "." in out_gene:
		out_gene = re.sub(r"\..*", "", out_gene) + "." + args['C_MutatedGenes'].lower()
	gff3 = args['gff3']
	current_directory = os.getcwd()

	# Logger configuration
	logger = logging.getLogger('LOTUS compare')
	logger.setLevel(logging.DEBUG)
	fh = logging.FileHandler(args['log'])
	fh.setLevel(logging.DEBUG)
	formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
	fh.setFormatter(formatter)
	logger.addHandler(fh)
	try:
		logger.info('Verification of input in config file')
		#vcf_filtered_files, vcf_pass_files, snp_profile_files, insertions_count_files, deletions_count_files = verif_input_config(args.config)
		logger.info('- Inputs file ok -')
	except ValueError:
		print('Problem with config file : ', sys.exc_info()[1])
		logger.error('- Problem with config file -')
		exit(1)

	try:
		verif_input(gff3)
		genes_positions = read_gff3(gff3)
		transcript_dico = genes_positions[2]
		gene_id_dico = genes_positions[1]
		gene_name_dico = genes_positions[0]
	except UnicodeDecodeError:
		print (f'{gff3} can not be read as a classic gff3 file : ', sys.exc_info()[1])
		logger.error(f'- {gff3} can not be read as a classic gff3 file ! -')
		exit(1)
	except ValueError:
		print (f'{gff3} is not a gff3 file (or the first # line missing) : ', sys.exc_info()[1])
		logger.error(f'- {gff3} is not a gff3 file ! -')
		exit(1)

	# Verification of given arguments
	vcf_filtered1 = args['filtered1']
	vcf_filtered2 = args['filtered2']
	vcf_strong_weak1 = args['passed1']
	vcf_strong_weak2 = args['passed2']

	out_indel = args['C_indel_profile_name']
	if "." in out_indel:
		out_indel = re.sub(r"\..*", "", out_indel)

	out_gene_test = args['genes_test']

	try:
		logger.info('Verification of outputs file')
		verif_output(out_gene_test)
		verif_output(vcf_strong_weak1)
		verif_output(vcf_strong_weak2)

		'''if enrichment:
			# verif_output(str(Path(out_gene).resolve().parent)+'/'+str(Path(file_one+'_'+file_two+'_ToppGene_enrichment.xlsx')))
			# verif_output(str(Path(out_gene).resolve().parent)+'/'+str(Path(file_one+'_'+file_two+'_ToppGene_enrichment.tsv')))
			# verif_output(str(Path(out_gene).resolve().parent)+'/'+str(Path(file_one+'_'+file_two+'_Panther_enrichment.xlsx')))
			# verif_output(str(Path(out_gene).resolve().parent)+'/'+str(Path(file_one+'_'+file_two+'_Panther_enrichment.tsv')))
			print('')

		logger.info('- Output files ok -')'''

	except ValueError:
		print ('Problem with one or more output files: ', sys.exc_info()[1])
		logger.error('- Problem with output files -')
		exit(1)


	infos = None
	if args['verbose_prints'].upper() == 'TRUE':
		print('Search for file Lotus_ExternalBases_202301.xlsx ...')
	infos = 'input/ressources/Lotus_ExternalBases_202301.xlsx'
	logger.info('Verification of {infos}')
	try:
		infos = verif_supplementary_information_file(infos, current_directory)
	except ValueError:
		print(f'Problem with {infos}: ', sys.exc_info()[0])
		logger.error('- Problem with {infos} -')
		exit(1)

	# Start
	logger.info('********************************************************************************************************************')
	logger.info('*** LOTUS compare module ***')
	no_argument = ''
	if enrichment:
		no_argument+=' --enrichment'
	configuration = args['dataset']
	logger.info(f'** cmd line : python lotus.py compare -c {str(configuration)} '+str(no_argument)+' **')
	logger.info('* Start comparing *')
	logger.info(f'Current directory : {Path().absolute()}')
	
	##############################
	# Comparison of both vcf files

	compare_vcf(out_gene, args, gene_name_dico, gene_id_dico, transcript_dico, infos, enrichment, logger)

	####################################################
	# Create snp mutation types plot and indel size plot

	snp_profile_files = []
	snp1 = args['output_path'] + 'samples/' + args ['file1'] + '/SNP_profile.tsv'
	snp_profile_files.append(snp1)
	snp2 = args['output_path'] + 'samples/' + args['file2'] + '/SNP_profile.tsv'
	snp_profile_files.append(snp2)
	insertions_count_files = []
	if os.path.exists(args['insert1']):
		insertions_count_files.append(args['insert1'])
		insertions_count_files.append(args['insert2'])
	else :
		print("No common insertion mutation found between both files.")
	deletions_count_files = []
	deletions_count_files.append(args['indel1'])
	deletions_count_files.append(args['indel2'])

	#kiwi
	if args['verbose_prints'].upper() == 'FALSE':
		print('Creating plots...')
	C_SNP_profile = args['C_SNP_profile_name']
	if not "NONE" in args['C_SNP_profile'].upper() and not "FALSE" in args['C_SNP_profile'].upper():
		if "." in C_SNP_profile:
			C_SNP_profile = re.sub(r"\..*", "", C_SNP_profile)
	out_C_SNP_profile = output_path + C_SNP_profile
	graph_snp(args, snp_profile_files, out_C_SNP_profile, logger)
	if os.path.exists(args['insert1']):
		graph_indel(args, insertions_count_files, deletions_count_files, out_indel, logger)
	else:
		print('Indel plot not created because no common insertion mutation found between both files.')

	logger.info('* End comparing *')
	logger.info('********************************************************************************************************************')

	# End








