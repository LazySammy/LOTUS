# Parameters

This is an exhaustive description of all possible parameters of [LOncoG](../).\
Each parameter must be filled in the configuration file before running the pipeline.

```<Name of parameter> = <required option>```
> module(s) = Filter, Summarise

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

#### 🔵 ```variants_selection_approach```
- Choose the method you want to use to select the variants to keep:
    - **change** > the variants that appeared or disappeared between times are selected
    - **common** > the variants that are common between times are selected
    - **all** > all the variants are selected (found in time 1, time 2, or both)
- Options: 
  - ```change``` ☑️ 
  - ```common```
  - ```all```
  > ```change``` is advised, so you can study variants that may have a relation with the disease evolution or treatment response.

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
  - ```none```

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

#### 🔵 ```filter_on_SIFT_score```
- The SIFT score is a measure of the impact of a mutation on the protein function. The lower it is, the better. It goes from 0 to 1.
- You can choose to filter the variants on this basis:
    - **yes** > only mutations with a SIFT score < max_SIFT_score will be conserved
    - **no** > all mutations will be conserved
- Options:
  - ```yes``` ☑️ 
  - ```no```
  > ```yes``` is recommended, as it can help to remove false positive mutations and focus on the most impactful ones (maybe driver mutations)

#### 🔵 ```filter_on_PolyPhen2_score```
- The PolyPhen2 score is a more recent measure of the impact of a mutation on protein function. The higher it is, the better. It goes from 0 to 1.
- You can choose to filter the variants on this basis:
    - **yes** > only mutations with a PolyPhen2 score > min_PolyPhen_score will be conserved
    - **no** > all mutations will be conserved
- Options:
    - ```yes``` ☑️
    - ```no```
    > ```yes``` is recommended, as it can help to remove false positive mutations and focus on the most impactful ones (maybe driver mutations)

#### 🟢 ```keep_variant_if_no_VAF_pop ```
- You can choose to keep the variants where the allelic frequency in the population is not available (not provided by annotators).
- Options:
  - ```yes``` ☑️ 
  - ```no```

#### 🟢 ```keep_variant_if_no_SIFT_info```
- You can choose to keep the variants where the SIFT score is not available (not provided by annotators).
- Options:
  - ```yes``` ☑️ 
  - ```no```

#### 🟢 ```keep_variant_if_no_Polyphen2_info ```
- You can choose to keep the variants where the PolyPhen2 score is not available (not provided by annotators).
- Options:
  - ```yes``` ☑️ 
  - ```no```

#### 🟢 ```min_alt_AD```
- The minimum allelic depth required to keep the variant. If the variant has an allelic depth > min_AD, it will be kept.
- Options:
  - ```20``` ☑️ 
  - ```<your_value>```

#### 🔵 ```min_DP```
- The minimum depth of sequencing required to keep the variant. If the variant has a depth > min_DP, it will be kept.
- Options:
  - ```80``` ☑️ 
  - ```<your_value>```
  
#### 🟢 ```min_alt_AD ```
- The minimum allelic depth required to keep the variant. If the variant has an allelic depth > min_alt_AD, it will be kept.
- Options : 
- ```20``` ☑️
- ```<your_value>```


#### 🟢 ```min_alt_MBQ```
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

#### 🔵 ```min_STRANDQ ```
- Phred-scaled minimum strand quality required to keep the variant. If the variant has a strand quality > min_STRANDQ, it will be kept.
- Options:
  - ```20``` ☑️ 
  - ```<your_value>```
  > provided by FilterMutectCalls from GATK, helps filtering out more germline-like variants

#### 🔵 ```max_SB```
- The strand bias is a measure of the bias in the allelic depth between the forward and reverse strands (between -1 and 1). If the variant has a strand bias < max_SB, it will be kept.
- Options:
  - ```0.5``` ☑️ 
  - ```<your_value>```
  > provided by ANNOVAR, Funcotator, SnpEff

#### 🔵 ```max_SOR ```
- Sometimes, a single SB value is not provided by annotators. In these cases, callers like Mutect2 are able to produce an SB table with 4 values.
From these values, we can compute a SOR value (Strand Odds Ratio). If the variant has a SOR < max_SOR, it will be kept.
- Options:
  - ```10``` ☑️ 
  - ```<your_value>```
  > calculation method available here : [GATK website](https://gatk.broadinstitute.org/hc/en-us/articles/360036361772-StrandOddsRatio)

#### 🔵 ```max_VAF_pop```
- The minimum allelic frequency in the population (number of times the mutation has been observed in the population). 
If the variant has a frequency in population > min_VAF_pop, it will be kept.
- Options:
  - ```0.002``` ☑️ 
  - ```<your_value>```
> provided by ANNOVAR (using gnomad40_exome database for example)

#### 🔵 ```max_VAF_sample```
- The maximum allelic frequency in the sample (number of times the mutation has been observed in sample). 
If the variant has a frequency < max_VAF_sample, it will be kept.
- Options:
  - ```0.3``` ☑️ 
  - ```<your_value>```
> provided by bcftools (```bcftools +fill-tags $input -Ov -o $output -- -t FORMAT/VAF```) or Mutect2 caller

#### 🔵 ```max_SIFT_score```
- The maximum SIFT score allowed to keep the variant. If the variant has a SIFT score < max_SIFT_score, it will be kept.
- If you chose not to filter on SIFT_score (```filter_on_SIFT_score = no```), it will just add SIFT threshold as a red line in protein impacts plots.
- Options:
  - ```0.05``` ☑️ 
  - ```<your_value>``
> provided by ANNOVAR (using dbnsfp41a for example), advised to keep only the most impactful mutations (potentially driver mutations)
[advised threshold](https://ionreporter.thermofisher.com/ionreporter/help/GUID-2097F236-C8A2-4E67-862D-0FB5875979AC.html)


#### 🔵 ```SIFT_preds_to_keep```
- You can choose to keep only the mutations predicted as deleterious by SIFT (D) or tolerated (T) or both.
- Options:
  - ```D``` ☑️ 
  - ```T```
  - ```all```
> provided by ANNOVAR (using dbnsfp41a for example), advised to keep only the most impactful mutations (potentially driver mutations)

#### 🔵 ```min_PolyPhen2_score```
- The minimum PolyPhen2 score allowed to keep the variant. If the variant has a PolyPhen2 score > min_PolyPhen2_score, it will be kept.
- Options:
  - ```0.5``` ☑️ 
  - ```<your_value>```
> provided by ANNOVAR (using dbnsfp41a for example), advised to keep only the most impactful mutations (potentially driver mutations)
[advised threshold](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4480630/)

####  🔵 ```PolyPhen2_preds_to_keep```
- You can choose to keep only the mutations predicted as probably damaging (D) or possibly damaging (P) or benign or all.
- `Options:
  - ```D``` ☑️ 
  - ```P```
  - ```B```
  - ```all```
> provided by ANNOVAR (using dbnsfp41a for example), advised to keep only the most impactful mutations (potentially driver mutations)

---------------------------

## SUMMARISE (S)
Here are the parameters you can set to summarise the filtered variants (statistics, plots, tables).

### Parameters
#### 🔵 ```reference_genome```
- The reference genome is required to get the gene names from the gene ids of reference genome. You must specifiy the file path and format (pk and fasta are available).
- Options:
  - ```hg38.pk``` ☑️ 
  - ```<your_reference_genome>```
> Make sure the version of your reference genome fits the version of the one you used to align your BAM before getting your VCF files.

#### 🔵 ```S_enrichment```
- You can choose to perform a Gene Ontology Enrichment Anlysis on the mutated genes found in the sample:
  - **Panther** > the Gene Ontology Enrichment Analysis will be performed using Panther database
  - **ToppGene** > the Gene Ontology Enrichment Analysis will be performed using ToppGene database
  - **both** > the Gene Ontology Enrichment Analysis will be performed using both Panther and ToppGene databases
  - **none** > no Gene Ontology Enrichment Analysis will be performed
- Options:
  - ```Panther``` ☑️ 
  - ```ToppGene```
  - ```both```
  - ```none```
> ```none``` is recommended for large datasets, Panther and ToppGene can only run if the list of genes contains less than 1000 names
>
> ### Parameters
