# LOTUS Pipeline

## Global Arguments

- **Module(s):** Merge
- **Sequencing Method:** WES (Whole Exome Sequencing)
- **VCF Annotation Method:** ANNOVAR
- **Discard Weak Variants:** Yes
- **Previous Folder:** annovar_filtered_PON_germline

### Input

- **Dataset File Name:** dataset.xlsx
- **Comparison Column 1:** comp1
- **Comparison Column 2:** comp2

### Secondary

- **Program Version:** 1.0
- **Log File:** LOTUS.log
- **Run Verbosity:** False
- **Colored Execution:** True

## Filter Module

- **Output Filtered Variants:** filtered.vcf
- **Output Passed Variants:** passed.vcf
- **Working Method:** Direct
- **DP (Minimum Coverage Depth):** 20
- **MBQ (Minimum Base Quality):** 10
- **MQSBZ (Mapping Quality Strand Bias Threshold):** 0.5
- **AF (Minimum Allele Frequency):** 5
- **POPAF (Population Allele Frequency):** 0.00001
- **Unpaired Reads:** True

## Summarise Module

### Parameters

- **Reference Genome File:** reference_genome.pk
- **Enrichment:** True
- **Supplementary Information from Filtered File:** False

### Outputs

- **Statistics File Name:** stats
- **Indel Deletion Profile:** None
- **Indel Insertion Profile:** None
- **Indel Profile Name:** indel_profile
- **Indel Profile Format:** png
- **SNP Profile Name:** SNP_profile
- **SNP Profile Formats:** tsv, png
- **Mutation Types Statistics Table File Name:** mutation_types
- **Mutation Types Statistics File Format:** xlsx
- **Mutated Genes File Name:** MutatedGenes
- **Mutated Genes File Format:** xlsx
- **ToppGene File Name:** ToppGene
- **ToppGene File Formats:** xlsx
- **Panther File Name:** Panther
- **Panther File Formats:** xlsx

## Compare Module

### Parameters

- **GFF3 File:** Homo_sapiens.GRCh38.108.chr.gff3
- **Enrichment:** True

### Outputs

- **SNP Profile File Name:** SNP_profile
- **SNP Profile File Formats:** png
- **Indel Profile File Name:** indel_profile
- **Indel Profile File Formats:** png
- **Compared Mutated Genes File Name:** compared
- **Compared Mutated Genes File Format:** xlsx
- **ToppGene File Name:** ToppGene
- **ToppGene File Formats:** xlsx
- **Panther File Name:** Panther
- **Panther File Formats:** xlsx

## Merge Module

### Inputs

- **Cytoband File Name:** hg38_cytoband.tsv
- **Pairs IDs Column Name:** patient id

### Parameters

- **Chromosome Step:** 500000
- **Perform Enrichment Analysis:** True
- **Min Degree of Mutated Genes:** 1
- **Max Degree of Mutated Genes:** 0
- **Min Subset Size:** 1
- **Max Subset Size:** 0

### Outputs

- **Chromosome Plot File Name:** chromosomes
- **Chromosome Plot File Formats:** png
- **Mutated Genes File Name:** genes_union
- **Mutated Genes File Formats:** xlsx
- **ToppGene File Name:** Panther
- **ToppGene File Formats:** xlsx
- **Panther File Name:** ToppGene
- **Panther File Formats:** xlsx
- **Upset Plot File Name:** upset
- **Upset Plot File Formats:** png
- **Upset Plot Threshold:** 3
