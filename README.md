<p align="center">
  <img src="https://i.postimg.cc/C1S6DQwH/work-in-progress.png" alt="Image" style="width:15%; height:auto;">
</p>

<p align="center"> 
    <a href="#contributors" alt="Contributors">
        <img src="https://img.shields.io/badge/contributors-3-lightblue" /></a>
    <a href="#backers" alt="Backers">
        <img src="https://img.shields.io/badge/backers-2-lightblue" /></a>
    <a href="#commits" alt="Commits">
        <img src="https://img.shields.io/badge/commits-87-lightblue" /></a>
    <a href="#coverage" alt="Coverage">
        <img src="https://img.shields.io/badge/coverage-85%25-lightgreen" /></a>
    <a href="#codacy" alt="Codacy">
        <img src="https://img.shields.io/badge/codacy-B-lightgreen" /></a>
    <a href="#languages" alt="Languages">
        <img src="https://img.shields.io/badge/language-Python 3-brightgreen" /></a>
    <a href="#version" alt="Version">
        <img src="https://img.shields.io/badge/version-2.0-brightgreen" /></a>
    <a href="#packages" alt="Packages">
        <img src="https://img.shields.io/badge/packages-conda, venv-brightgreen" /></a>
    <a href="#parameters" alt="Parameters">
        <img src="https://img.shields.io/badge/parameters-100-lightyellow" /></a>
    <a href="#example" alt="Example">
        <img src="https://img.shields.io/badge/test dataset-available-lightyellow" /></a>
    <a href="#validation" alt="Validation">
        <img src="https://img.shields.io/badge/validation-TNBC, glioblastoma-yellow" /></a>
    <a href="#institute" alt="Institute">
        <img src="https://img.shields.io/badge/institute-Institut de Cancérologie de l'Ouest%20-orange" /></a>
    <a href="#country" alt="Country">
        <img src="https://img.shields.io/badge/made in-🇫🇷France-black" /></a>
</p>

### $${\color{lightblue}LOncoG: \space a \space software \space for \space Longitudinal \space OncoGenomics \space analysis}$$
This software plots, compare and merge information from all exomes of a cohort of cancer patients before and after treatment. It also includes a customizable filter to help you removing remaining germline and/or non driver mutations. The Filter, Summarise, Compare and Merge modules can be run separately or all together. If your study is not longitudinal, you can just run the first two modules to get a graphical and statistical sumup of the most impactant variants from your WES data. The software is designed to be user-friendly and to be used by bioinformaticians, biologists and even clinicians research teams.

![LOncoG simple workflow](tutorials/pictures/LOncoG_workflow.png)

## Versions
☑️ ```1.0```: the first version (06/2023) is available on [GitHub](https://github.com/gsiekaniec/LOTUS). \
✅ ```2.0```: the new version (04/2024) is a pre-release and is the current version. \
☑️ ```3.0```: the final version (07/2024) will be the version associated with the futur potential paper. \
3.0 version will input more annotations formats, perform more API requests and a wider choice of plot shapes. 
LOncoG will perform clustering on patients similarities and will run with a multiprocessing option.

## Warnings
- LOncoG is still in development and should not be used for clinical purposes yet. 
- The software is not yet optimized for large cohorts of patients.
- Protein impacts predictions only works with ANNOVAR dbsnfp41a database for the moment.
- Variant allelic frequency in population only works with ANNOVAR gnomad40 database for the moment.

## Project Organization
The project is organized as follows:
```Project/
├── python_scripts/
│   ├── reusable_functions/         -> Functions used in the main scripts to parse vcf, etc.
│   ├── api_requests/               -> Functions to request information from external databases.
│   └── modules/                    -> Main modules scripts to run the software.
├── input/
│   ├── resources/                  -> Resources used in the software (reference genome, etc).
│   └── vcf/                        -> VCF default input folder.
├── environment/
│   ├── conda/                      -> Conda environment files (yml, txt).
│   └── venv/                       -> Virtual environment file (txt).
├── tutorials/
│   ├── pictures/                   -> Pictures used in the tutorials.
│   └── examples/                   -> README for parameters, input examples.
├── toy_dataset/
│   ├── toy_vcf/                    -> VCF toy dataset.
│   └── toy_output/                 -> Output of the toy dataset.
├── output/
│   └── hour_month_year/            -> Auto-generated folder with the date and time of the run
|       ├── samples/                -> Filtered vcf and Summarise plots, one subfolder per sample/exome.
|       ├── comparisons/            -> Plots from Comparison module, one subfolder per patient (pair).
|       └── merge/                  -> Plots from Merge module, no subfolder.        
├── README.md                       -> This file, the main README.
└── logs/                           -> Logs are created here.
```
Being at ease with the project organization is important to get a good understanding of the software functioning.
*Please make sure to read the [PARAMETERS.md](tutorials/PARAMETERS.md) file to understand how to choose the parameters for the software.*

## Installation
### Prerequisites
Python 3.9 and Conda need to be installed on your machine. You can also use pip instead of Conda.

### Conda (recommended)
Clone the repository and navigate to the project directory.
Create a Conda environment using the [loncog_conda.yml](environment/conda/loncog_conda.yml) file with the following command: 
```bash
conda env create -f environment/conda/loncog_conda.yml --name loncog
```
Then, activate the environment with: 
```bash
conda activate loncog
```

You are also free to use the [loncog_conda.txt](environment/conda/loncog_conda.txt) file if you prefer to install the packages manually with: 
```bash
conda env create --file environment/conda/loncog_conda.txt --name loncog
```
Then, activate the environment with: 
```bash
conda activate loncog
```

### Virtual environment
Clone the repository and navigate to the project directory.
Create a virtual environment using the [loncog_venv.txt](environment/venv/loncog_venv.txt) file with the following commands: 
```bash
python3 -m venv venv 
source venv/bin/activate  
pip install -r environment/venv/loncog_venv.txt
```
## Compatible annotators
LOncoG is based on vcf annotations from input files, here is a list of compatible annotators for v2.0 (more to come):
| Colonne 1 | Colonne 2 | Colonne 3 |
| --------- | --------- | --------- |
| Donnée 1  | Donnée 2  | Donnée 3  |
| Donnée 4  | Donnée 5  | Donnée 6  |

## Usage
This section describes the quick steps you need to follow to run the software.\

### Parameters
Choose the parameters for the software in the [config.txt](config.txt) file.\
An exhaustive description with examples and advices is available in [PARAMETERS.md](tutorials/PARAMETERS.md) file.\

### Run the software