# North Dakota Bumble Bee Landscape Genetics

---

## Overview

This repository contains all code used for my graduate research project, master’s thesis, and associated manuscript currently in preparation for publication. The analyses focus on ecological and genetic data processing, statistical modeling, and visualization workflows developed to support the results presented in the thesis and forthcoming paper.

The goal of this repository is to ensure **reproducibility, transparency, and long-term accessibility** of all analytical methods used in this research.


---

## Repository Structure
├── NeEstimator/              # Effective population size (Ne) analyses and outputs

├── R-code/                   # R scripts for data cleaning, statistical analyses, and figures

├── UMGC_scripts/             # Perl scripts provided by UMGC for bioinformatics

├── bioinformatics_scripts/   # Genomic data processing (VCF, PLINK, ADMIXTURE, etc.)

├── data/                     # Metadata associated with the project

├── unix_code/                # Standalone Unix/bash commands and text-processing workflows

└── README.md

---

## Methods Summary

This repository includes code for:

* Data cleaning and formatting
* Genetic data processing (e.g., VCF handling, filtering, summary statistics)
* Statistical analyses and modeling
* Figure and table generation for thesis and manuscript

Detailed methodological descriptions are provided in the thesis and manuscript; scripts are commented where appropriate to clarify analytical steps.

---

## Software & Dependencies

Primary languages and tools used:

* **R** (tidyverse, hierfstat, and other domain-specific packages)
* **Command-line tools** (e.g., bcftools, PLINK, BWA where applicable)
* **Excel / CSV workflows** for metadata management

Specific package versions may be documented within individual scripts.

---

## Reproducibility Notes

* Scripts are designed to be run in sequence where applicable.
* File paths may need to be updated to match local directory structures.
* Some analyses require large datasets not included in this repository.

If you have questions about reproducing specific results, feel free to contact me.

---

## Author

**Rhiannon Hall**
Graduate Researcher
rhiannonchall@gmail.com


