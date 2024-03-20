# Parameters

This is an exhaustive description of all possible parameters of [SLOOP](../).\
Each parameter must be filled in the configuration file before running the pipeline.\

```<Name of parameter> = <required option>```\
module(s) = Filter, Summarise\
min_DP = 40

☑️ Default value\
🔵 Key parameter\
🟢 Secondary parameter

### Global arguments
#### LOTUS
##### 🔵 module(s)
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

##### 🔵 discard_weak_variants
- Choose if you want to keep the variants that don't pass the filter you chose:
    - **Yes** > the variants that don't pass the filter you chose are removed (strong variants)
    - **Summarise** > the variants that don't pass the filter you chose are kept (weak variants)
- Options: 
  - ```yes``` ☑️ 
  - ```no```

##### 🔵 output_folder_path
- Choose the path where you want to save the results of the analysis:
    - **new** > creates a new folder in "output" folder with date as a name (format : YYYY-MM-DD)
    - **<your path>** > creates a new folder located at the path specified
    > Note: if you already ran Filter and Summarise in a folder, you can use its path to complete it with Compare and Merge analysis
    (if you choose to run Filter or Summarise again, it will overwrite the previous results if you don't change the folderpath)
- Options: 
  - ```new``` ☑️ 
  - ```<you path>```

##### 🔵 keep_filtered_vcf_after_run
- After FILTER execution, filtered vcf files are created (indicate which variants did pass the filter), but
  you can save a lot of space by removing them when no more required:
    - **yes** > the filtered vcf files are kept in the "samples" subfolder (in your output folder)
    - **no** > the filtered vcf files are removed after the execution of Compare module
- Options: 
  - ```no``` ☑️ 
  - ```yes```
    > Caution! Summarise and Compare modules need the filtered vcf files to be correctly run!
    (if you removed the filtered vcf files by error, and you run Summarise or Compare, it will re-run Filter module)

#### Inputs
##### 🔵 vcf_folder_path
- We need to know the path of the folder containing all of your annotated vcf, as LOTUS can find them.
- Options: 
  - ```input/vcf/``` ☑️ 
  - ```<your_path>```

##### 🔵 dataset_path
- We need to know which files are paired within your dataset, so you need to fill a table with the filenames of each pair.
- Options: 
  - ```input/dataset.xlsx``` ☑️ 
  - ```<your_path>```

##### 🔵 time1_column_name
- In your table, all filenames from time1 must be in the same column, so we need to know its name.
- Options: 
  - ```time1``` ☑️ 
  - ```<your_column_name_for_time1>```
  
##### 🔵 time2_column_name
- In your table, all filenames from time2 must be in the same column, so we need to know its name.
- Options: 
  - ```time2``` ☑️ 
  - ```<your_column_name_for_time2>```

##### 🟢 pair_names_column_name
- In your table, you can add an optional column of pair names, that LOTUS will display in analysis results (plots, etc):
  - **<your_column_name_for_pairs_names>** > you can put your patients ids there, the outputs will be clearer (short names)
  - **none** > if you don't have a column with pair names, we will use ```file_time1_name___file_time2_name``` as pair id (very long names)
- Options: 
  - ```patients``` ☑️ 
  - ```<your_column_name_for_pairs_names>```
  - ```None```

#### Execution
##### 🟢 log
- For people running LOTUS on linux, you can specify a log file to save the execution logs:
  - **<your_log_file_name>** > the log file will be created in the output folder with this name (include .log)
  - **<none>** > no log file will be created

- Options: 
  - ```LOTUS.log``` ☑️ 
  - ```<your_log_file_name>```
  - ```none```

##### 🟢 verbose_prints
- During LOTUS running, you can choose to display synthetic or detailed information in the console:
    - **yes** > detailed information will be displayed in the console
    - **no** > synthetic information will be displayed in the console
- Options: 
  - ```no``` ☑️ 
  - ```yes```

##### 🟢 colored_execution
- During LOTUS running, the information appearing in the console can be colored or not:
    - **yes** > the information will be colored in the console (recommanded for recent consoles such as VSCode, PyCharm, etc)
    - **no** > the information will not be colored in the console (recommanded for old consoles such as PyScripter, etc)
- Options: 
  - ```yes``` ☑️ 
  - ```no```