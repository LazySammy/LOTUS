#!/usr/bin/env python3
# coding: utf-8

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#   LOTUS : a program for LOngiTUdinal comparative genomic Study
#   Authors: Samuel Besseau, G. Siekaniec, W. Gouraud
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import os
import logging
from playsound import playsound
import sys
import time
from python_scripts import filter, summarise, compare, merge
from python_scripts.manage_parameters import calculate_time
from python_scripts.manage_parameters import list_colors
from python_scripts.manage_parameters import compare_parameters
from python_scripts.manage_parameters import manage_parameters

dict_para = manage_parameters('config.txt')
if dict_para['colored_execution'].upper() == 'TRUE' or dict_para['colored_execution'].upper() == 'YES':
    colored = True
else:
    colored = False
dict_colors = list_colors()
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
if colored:
    print(f"{yellow1}{bold}Reading parameters...{reset_bold}")
else:
    print("Reading parameters...")
start_time = time.time()
parameters_time = time.time()
previous = dict_para['previous_folder']
modules = dict_para['module(s)'].upper()
output_path = dict_para['output_path']
# if previous.upper() != 'FALSE' and not os.path.isdir(output_path + 'samples'):
#     sys.exit(f"{dark_red}{bold}ERROR: {reset_bold}Without filter and summarise files, you can't run compare or merge modules :"
#              f" \n-> please run LOTUS with filter and summarise modules in this folder, or choose another previous_folder name with a samples folder."
#              f"\n{pink}-> the samples folder must contain every output from filter and summarise modules, for each sample folder.")

if colored:
    print(f"\n{yellow1}Working in {bold}{output_path}{reset_bold} directory...")
else:
    print(f"\nWorking in {output_path} directory...")

i = 0
for file in dict_para['vcf_files']:
    sample_id = file.split('.funco')[0]
    output_path_sample = output_path + 'samples/' + sample_id + '/'
    dict_para['output_path_sample'] = output_path_sample
    if '.vcf' in output_path_sample:
        output_path_sample = output_path_sample.replace('.vcf', '')
    if not os.path.exists(output_path_sample):
        os.makedirs(output_path_sample)
    index = dict_para['vcf_files'].index(file) + 1
    pre_filter_time = time.time()

    if 'FILTER' in modules and previous.upper() == 'FALSE':
        if colored:
            print(f"\n{yellow2}{bold}Filtering vcf file n°{index} ({sample_id}):{reset_bold}")
        else:
            print(f"\nFiltering vcf file n°{index} ({sample_id}):")
        input_vcf = 'input/vcf/' + file
        filter.main(dict_colors,dict_para,input_vcf, dict_para, output_path_sample)
        filter_time = calculate_time(dict_para, pre_filter_time, 'Filter')
    elif 'FILTER' in modules and previous.upper() != 'FALSE':
        #Replacing filter files in your previous_folder (if you want to avoid that, put False to create a new folder)
        #sys.exit("We can't use filter module using previous files, since it is the first module to produce files")
        for filename in os.listdir(output_path_sample):
            file_path = os.path.join(output_path_sample, filename)
            if filename.endswith(".vcf"):
                os.remove(file_path)
        if colored:
            print(f"\n{yellow2}{bold}Filtering vcf file n°{index} ({sample_id})...{reset_bold}")
        else:
            print(f"\nFiltering vcf file n°{index} ({sample_id})...")
        input_vcf = 'input/vcf/' + file
        filter.main(dict_colors, dict_para, input_vcf, dict_para, output_path_sample)
        filter_time = calculate_time(dict_para, pre_filter_time, 'Filter')

    filter_time = time.time()

    if 'SUMMARISE' in modules and previous.upper() == 'FALSE' and 'FILTER' not in modules:
        sys.exit("We can't use summarise module without filter before, or without using a previus folder containing filter files.")
    elif 'SUMMARISE' in modules:
        if dict_para["previous_folder"].upper() != 'FALSE':
            for filename in os.listdir(output_path_sample):
                file_path = os.path.join(output_path_sample, filename)
                if not filename.endswith(".vcf") and not filename.endswith("mutations.txt"):
                    os.remove(file_path)
        if colored:
            print(f"\n{yellow3}{bold}Summarising vcf file n°{index} ({sample_id}):{reset_bold}")
        else:
            print(f"\nSummarising vcf file n°{index} ({sample_id}):")
        input_vcf_filtered = output_path_sample + 'filtered.vcf'
        input_vcf_passed = output_path_sample + 'passed.vcf'
        if i != len(dict_para['vcf_files']) - 1:
            last = False
        else:
            last = True
        summarise.main(dict_para, dict_colors, dict_para, input_vcf_filtered, input_vcf_passed, output_path_sample, last)
        summarise_time = calculate_time(dict_para, filter_time, 'Summarise')
    i+=1

if 'COMPARE' in modules:
    n_comp = int(len([f for f in os.listdir('input/vcf/') if f.endswith('.vcf')])/2)
    for i in range(n_comp):
        before_compare = time.time()
        dict_para = compare_parameters(dict_para,i)
        if colored:
            print(f"\n{yellow4}{bold}Comparing vcf pair n°{i+1} ({dict_para['file1']} & {dict_para['file2']}):{reset_bold}")
        else:
            print(f"\nComparing vcf pair n°{i+1} ({dict_para['file1']} & {dict_para['file2']}):")
        compare.main(dict_para)
        compare_time = calculate_time(dict_para, before_compare, 'Compare')

if previous.upper() != 'FALSE':
    compare_time = time.time()

if 'MERGE' in modules:
    n_pairs = sum(os.path.isdir(os.path.join(output_path + "comparisons/", item)) for item in os.listdir(output_path + "comparisons/"))
    pair = 'pair' if n_pairs == 1 else 'pairs'
    if colored:
        print(f"\n{yellow5}{bold}Merging {n_pairs} {pair}:{reset_bold}")
    else:
        print(f"\nMerging {n_pairs} {pair}:")
    output_path_merge = output_path + 'merge/'
    if os.path.exists(output_path_merge) and previous.upper() == 'FALSE':
        os.rmdir(output_path_merge)
    elif os.path.exists(output_path_merge) and previous.upper() != 'FALSE':
        for file_name in os.listdir(output_path_merge):
            file_path = os.path.join(output_path_merge, file_name)
            os.remove(file_path)
        os.rmdir(output_path_merge)
    os.makedirs(output_path_merge)

    dict_para['M_MutatedGenes_name'] = output_path_merge + dict_para['M_MutatedGenes_name']
    dict_para['chromosomes_plot_name'] = output_path_merge + dict_para['chromosomes_plot_name']
    merge.main(dict_para, output_path)
    merge_time = calculate_time(dict_para, compare_time, 'Merge')
playsound("input/ressources/job_done.wav")
final_time = calculate_time(dict_para, start_time, 'End')

__version__ = '0.0.1'
lotus_ascii = r'''
-----------------------------------------------------------------------------------

-----------------------------------------------------------------------------------

                                        .
                                      .. ..
                                     .     .
                                    .       .
          ...           . ..       .         .       ....           ...
          .    ..      ..   ..     .          .    .    .      ...   ..
          .       ..   .      ..  .           .. .      .   ..       .
           .         ...        . .           ...       ....         .
           ..          .         ..           .         ..          .
 .    ..... ..         .           .         .          ..         .  ....   ..
  .          ...       .            .       .           .        ..          .
   ..             ..   ..            .     .            .    ..             .
     .               .. .            ..   .            ....               ..
      .                 ..            .  .            ..                 ..
       ..                  ..         .. .          ..                  .
         .                   ..        ...       ..                   ..
           .                    .      ..      ..                   ..
             ..                   .    ..    ..                  ..
                ...                 .  ..  ..                ...
                      ....            . ...           ....
                            ......   .......   .....
                            _       _
                           | |     | |
                           | | ___ | |_ _   _ ___
                           | |/ _ \| __| | | / __|
                           | | (_) | |_| |_| \__ \
                           |_|\___/ \__|\__,_|___/
                         
                         
'''

ascii_help= r'''

-----------------------------------------------------------------------------------

-----------------------------------------------------------------------------------
'''

ascii_summarise = r'''                                                      _
                                                     (_)         
             ___ _   _ _ __ ___  _ __ ___   __ _ _ __ _ ___  ___ 
            / __| | | | '_ ` _ \| '_ ` _ \ / _` | '__| / __|/ _ \
            \__ \ |_| | | | | | | | | | | | (_| | |  | \__ \  __/
            |___/\__,_|_| |_| |_|_| |_| |_|\__,_|_|  |_|___/\___|                                                         

-----------------------------------------------------------------------------------

-----------------------------------------------------------------------------------
'''

ascii_filter= r'''                              __ _ _ _            
                             / _(_) | |           
                            | |_ _| | |_ ___ _ __ 
                            |  _| | | __/ _ \ '__|
                            | | | | | ||  __/ |   
                            |_| |_|_|\__\___|_|

-----------------------------------------------------------------------------------

-----------------------------------------------------------------------------------
'''

ascii_compare= r'''                  ___ ___  _ __ ___  _ __   __ _ _ __ ___ 
                 / __/ _ \| '_ ` _ \| '_ \ / _` | '__/ _ \
                | (_| (_) | | | | | | |_) | (_| | | |  __/
                 \___\___/|_| |_| |_| .__/ \__,_|_|  \___|
                                    | |                   
                                    |_|                  

-----------------------------------------------------------------------------------

-----------------------------------------------------------------------------------
'''

ascii_merge= r'''                        _ __ ___   ___ _ __ __ _  ___ 
                       | '_ ` _ \ / _ \ '__/ _` |/ _ \
                       | | | | | |  __/ | | (_| |  __/
                       |_| |_| |_|\___|_|  \__, |\___|
                                            __/ |     
                                           |___/

-----------------------------------------------------------------------------------

-----------------------------------------------------------------------------------
'''

# Logger configuration
logger = logging.getLogger('LOTUS main')
logger.setLevel(logging.DEBUG)
logger.info(f'---------------- LOTUS v{__version__} ----------------')