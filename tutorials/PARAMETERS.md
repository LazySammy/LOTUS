# Parameters

This is an exhaustive description of all possible parameters of [SLOOP](../).\
Each parameter must be filled in the configuration file before running the pipeline.\

```<Name of parameter> = <required option>```\
module(s) = Filter, Summarise\
min_DP = 40

☑️ Default value\
🔵 Key parameter\
🟢 Secondary parameter

## Global arguments
Here are some key arguments about the program in general, and also about input and output folders.

### LOTUS
#### 🔵 ```module(s)```
- Choose the module(s) you want to run:
    - **Filter** the variants among their sequencing quality and allelic frequency
    - **Summarise** the filtered variants (statistics, mutation types, mutated genes, etc)
    - **Compare** the filtered variants from time 1 and time 2 for each patient
    - **Merge** the comparisons results together, to get a cross analysis of all patients
- Options: 
  - ```Filter, Summarise, Compare, Merge``` ☑️ 
  - ```Filter```
  - ```Summarise```
  - ```Compare```
  - ```Merge```

#### 🔵 ```discard_weak_variants```
- Choose if you want to keep the variants that don't pass the filter you chose:
    - **Yes** > the variants that don't pass the filter you chose are removed (strong variants)
    - **Summarise** > the variants that don't pass the filter you chose are kept (weak variants)
- Options: 
  - ```yes``` ☑️ 
  - ```no```

#### 🔵 ```output_folder_path```
- Choose the path where you want to save the results of the analysis:
    - **new** > creates a new folder in "output" folder with date as a name (format : YYYY-MM-DD)
    - **<your path>** > creates a new folder located at the path specified
    > Note: if you already ran Filter and Summarise in a folder, you can use its path to complete it with Compare and Merge analysis
    (if you choose to run Filter or Summarise again, it will overwrite the previous results if you don't change the folderpath)
- Options: 
  - ```new``` ☑️ 
  - ```<you path>```

#### 🔵 ```keep_filtered_vcf_after_run```
- After FILTER execution, filtered vcf files are created (indicate which variants did pass the filter), but
  you can save a lot of space by removing them when no more required:
    - **yes** > the filtered vcf files are kept in the "samples" subfolder (in your output folder)
    - **no** > the filtered vcf files are removed after the execution of Compare module
- Options: 
  - ```no``` ☑️ 
  - ```yes```
    > Caution! Summarise and Compare modules need the filtered vcf files to be correctly run!
    (if you removed the filtered vcf files by error, and you run Summarise or Compare, it will re-run Filter module)

### Inputs
#### 🔵 ```vcf_folder_path```
- We need to know the path of the folder containing all of your annotated vcf, as LOTUS can find them.
- Options: 
  - ```input/vcf/``` ☑️ 
  - ```<your_path>```

#### 🔵 ```dataset_path```
- We need to know which files are paired within your dataset, so you need to fill a table with the filenames of each pair.
- Options: 
  - ```input/dataset.xlsx``` ☑️ 
  - ```<your_path>```

#### 🔵 ```time1_column_name```
- In your table, all filenames from time1 must be in the same column, so we need to know its name.
- Options: 
  - ```time1``` ☑️ 
  - ```<your_column_name_for_time1>```
  
#### 🔵 ```time2_column_name```
- In your table, all filenames from time2 must be in the same column, so we need to know its name.
- Options: 
  - ```time2``` ☑️ 
  - ```<your_column_name_for_time2>```

#### 🟢 ```pair_names_column_name```
- In your table, you can add an optional column of pair names, that LOTUS will display in analysis results (plots, etc):
  - **<your_column_name_for_pairs_names>** > you can put your patients ids there, the outputs will be clearer (short names)
  - **none** > if you don't have a column with pair names, we will use ```file_time1_name___file_time2_name``` as pair id (very long names)
- Options: 
  - ```patients``` ☑️ 
  - ```<your_column_name_for_pairs_names>```
  - ```None```

### Execution
#### 🟢 ```log```
- For people running LOTUS on linux, you can specify a log file to save the execution logs:
  - **<your_log_file_name>** > the log file will be created in the output folder with this name (include .log)
  - **<none>** > no log file will be created
- Options: 
  - ```LOTUS.log``` ☑️ 
  - ```<your_log_file_name>```
  - ```none```

#### 🟢 ```verbose_prints```
- During LOTUS running, you can choose to display synthetic or detailed information in the console:
    - **yes** > detailed information will be displayed in the console
    - **no** > synthetic information will be displayed in the console
- Options: 
  - ```no``` ☑️ 
  - ```yes```

#### 🟢 ```colored_execution```
- During LOTUS running, the information appearing in the console can be colored or not:
    - **yes** > the information will be colored in the console (recommanded for recent consoles such as VSCode, PyCharm, etc)
    - **no** > the information will not be colored in the console (recommanded for old consoles such as PyScripter, etc)
- Options: 
  - ```yes``` ☑️ 
  - ```no```
  

---------------------------

---------------------------
## FILTER
Here are the parameters you can set to filter the variants with LOTUS first module.
If you don't want a filter to be applied, just set the parameter to 0 if min threshold, and 1000 if max threshold.

#### 🟢 ```colored_execution```
- During LOTUS running, you can choose the execution method (efficiency, memory).
    - **Direct** > the execution will be faster but will use more memory
    - **In memory** > the execution will be slower but will use less memory
- Options:
  - ```Direct``` ☑️ 
  - ```In memory```

#### 🔵 ```filter_on_mutation_location```
- The annotators (such as ANNOVAR) can provide mutations localisations, so you can choose to filter the variants on this basis:
    - **yes** > only mutations from exons and splicing (2 bp adjacent to the splicing site) will be conserved, recommended for WES data (Whole Exome Sequencing)
    - **no** > all mutations types will be conserved, recommended for WGS data (Whole Genome Sequencing)
    - **<your_mutations>** > you can specify the mutations you want to keep (ex: MISSENSE, NONSENSE, etc)
    > Example: 'exonic', 'splicing', 'ncRNA', 'ncRNA_intronic', 'ncRNA_exonic', 'UTR5', 'UTR3', 'intronic', 'upstream', 'downstream', 'intergenic'
    are all possible values for ANNOVAR, you can tell LOTUS the ones you want to keep, separated by a comma
- Options:
  - ```yes``` ☑️ 
  - ```no```
  - ```<your_mutations>```

#### 🟢 ```min_AD```
- The minimum allelic depth required to keep the variant. If the variant has an allelic depth > min_AD, it will be kept.
- Options:
  - ```20``` ☑️ 
  - ```<your_value>```

#### 🔵 ```min_DP```
- The minimum depth of sequencing required to keep the variant. If the variant has a depth > min_DP, it will be kept.
- Options:
  - ```80``` ☑️ 
  - ```<your_value>```

#### 🟢 ```min_MBQ```
- The minimum base quality required to keep the variant. If the variant has a median base quality > min_MBQ, it will be kept.
- Options:
  - ```20``` ☑️ 
  - ```<your_value>```

#### 🟢 ```min_QUAL```
- QUAL is a quality score, proportional to the confidence you can have in the variant call. If the variant has a quality > min_QUAL, it will be kept.
- Options:
  - ```50``` ☑️ 
  - ```<your_value>```
  > provided by ANNOVAR

#### 🔵 ```max_SB```
- The strand bias is a measure of the bias in the allelic depth between the forward and reverse strands (between -1 and 1). If the variant has a strand bias < max_SB, it will be kept.
- Options:
  - ```0.5``` ☑️ 
  - ```<your_value>```
  > provided by ANNOVAR, Funcotator, SnpEff

#### 🔵 ```max_VAF_pop```
- The minimum allelic frequency in the population (number of times the mutation has been observed in the population). 
If the variant has a frequency in population > min_VAF_pop, it will be kept.
VAF_sample_max_threshold = 0.3
#POPAF = 0.00001
> provided by ANNOVAR (using gnomAD database for example)

#### 🔵 ```max_VAF_sample```
- The maximum allelic frequency in the sample (number of times the mutation has been observed in sample). 
If the variant has a frequency < max_VAF_sample, it will be kept.
- Options:
  - ```0.002``` ☑️ 
  - ```<your_value>```
  > provided by bcftools (```bcftools +fill-tags $input -Ov -o $output -- -t FORMAT/VAF```)

## SUMMARISE
Here are the parameters you can set to summarise the filtered variants (statistics, plots, tables).

### Parameters



genome = reference_genome.pk
enrichment=False