#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#   LOncoG : a software for Longitudinal OncoGenomics analysis
#   Authors: S. Besseau, G. Siekaniec, W. Gouraud
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def read_vcf(vcf_file, get_variant_idx=False):
	'''
	Create a generator to read the vcf file
	Input : path to the vcf file
	Output : generator returning a dictionary containing the index position or the strip and split line
	'''
	idx = {}
	a = 0
	with open(vcf_file, 'r', encoding='latin-1') as vcf:
		if not get_variant_idx:
			for line in vcf:
				if line != '':
					if line[0] != '#':
						a += 1  # if line is not a header
						yield (line.strip().split('\t'))  # add the line as a list of strings in the object
					elif line.startswith('#') and not line.startswith('##'):
						# if we are on the header line
						line = line.lstrip('#').split('\t')  # make it a list removing \t
						idx['idx_chr'] = line.index('CHROM')  # get the index for each field and put them in index dict
						idx['idx_pos'] = line.index('POS')
						idx['idx_ref'] = line.index('REF')
						idx['idx_alts'] = line.index('ALT')
						idx['idx_filt'] = line.index('FILTER')
						idx['idx_info'] = line.index('INFO')
						idx['idx_format'] = line.index('FORMAT')
						idx['idx_values'] = idx['idx_format'] + 1
						yield (idx)  # add the index as a dictionary in the object
		else:
			print('error')
			pos = vcf.tell()
			line = vcf.readline()
			if line.startswith('#') and not line.startswith('##'):
				line = line.lstrip('#').split('\t')
				idx['idx_chr'] = line.index('CHROM')
				idx['idx_pos'] = line.index('POS')
				idx['idx_ref'] = line.index('REF')
				idx['idx_alts'] = line.index('ALT')
				idx['idx_filt'] = line.index('FILTER')
				idx['idx_info'] = line.index('INFO')
				idx['idx_format'] = line.index('FORMAT')
				idx['idx_values'] = idx['idx_format'] + 1
				yield (pos, idx)
			while line != '':
				pos = vcf.tell()
				line = vcf.readline()
				if line != '':
					if line[0] != '#':
						yield (pos, line.strip().split('\t'))
					elif line.startswith('#') and not line.startswith('##'):
						line = line.lstrip('#').split('\t')
						idx['idx_chr'] = line.index('CHROM')
						idx['idx_pos'] = line.index('POS')
						idx['idx_ref'] = line.index('REF')
						idx['idx_alts'] = line.index('ALT')
						idx['idx_filt'] = line.index('FILTER')
						idx['idx_info'] = line.index('INFO')
						idx['idx_format'] = line.index('FORMAT')
						idx['idx_values'] = idx['idx_format'] + 1
						yield (pos, idx)
def get_vcf_header(vcf_file):
	'''
	Create a generator for the headers line of a vcf
	'''
	with open(vcf_file, 'r', encoding='latin-1') as vcf:
		for line in vcf:
			if line != '':
				if line.startswith('#'):
					yield (line.strip())
