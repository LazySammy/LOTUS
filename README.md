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

## $${\color{lightblue}LOncoG: \space a \space software \space for \space Longitudinal \space OncoGenomics \space analysis \n version \space 2.0}$$
This software plots, compare and merge information from all exomes of a cohort of cancer patients before and after treatment. It also includes a customizable filter to help you removing remaining germline and/or non driver mutations. The Filter, Summarise, Compare and Merge modules can be run separately or all together. If your study is not longitudinal, you can just run the first two modules to get a graphical and statistical sumup of the most impactant variants from your WES data. The software is designed to be user-friendly and to be used by bioinformaticians, biologists and even clinicians research teams.

## Warnings
- LOncoG is still in development and should not be used for clinical purposes yet. 
- The software is not yet optimized for large cohorts of patients.
- Protein impacts predictions only works with ANNOVAR dbsnfp41a database for the moment.
- Variant allelic frequency in population only works with ANNOVAR gnomad40 database for the moment.

## Project Organization
The project is organized as follows:
```Project/
├── python_scripts/
│   ├── reusable_functions/         -> Functions used in the main scripts to parse vcf, retrive information from df, etc.
│   ├── api_requests/               -> Functions to request information from external databases (Gene Ontology, etc). COSMIC, OncoKB, ClinVar coming soon.
│   └── modules/                    -> Main scripts to run the software (Filter, Summarise, Compare, Merge).
├── input/
│   ├── resources/                  -> Resources used in the software (databases, reference genome, etc).
│   └── vcf/                        -> VCF default input folder.
├── tutorials/
│   ├── pictures/                   -> Pictures used in the tutorials.
│   └── examples/                   -> README for parameters, input examples.
├── output/                         -> Default output folder is create here.
├── README.md                       -> This file, the main README.
└── logs/                           -> Logs are created here.
```

### Heading 1
#### Heading 2
##### Heading 3
###### Heading 4
####### Heading 5
######## Heading 6

## Emphasis

You can add emphasis to your text using asterisks or underscores.

*Italic text*
_Italic text_

**Bold text**
__Bold text__

## Lists

You can create ordered and unordered lists.

### Unordered List

- Item 1
- Item 2
- Item 3

### Ordered List

1. First item
2. Second item
3. Third item

## [Docs Contribution Guide](https://www.codecademy.com/pages/contribute-docs)

You can add links to your text.

[GitHub](https://github.com)

## Images

You can also add images to your README.

![Image](https://example.com/image.jpg)

## Code Blocks

You can include code blocks using triple backticks

## Parameters
[Exhaustive description](tutorials/PARAMETERS.md)

## Table
First Header  | Second Header
------------- | -------------
Content Cell  | Content Cell
Content Cell  |  \| 
