# protein-annotation-bias

# Protein Annotation Bias Analysis

## Introduction
This project analyzes biases in human proteome experimental GO annotations and their correlation with research interest over time (2013-2022).

## Dependencies
- Python 3.7+
- R 4.0+
- NCBI EDirect for PubMed data retrieval
- Required Python packages: `Bio`, `obonet`, `pandas`, `subprocess`, `datetime`
- Required R packages: `tidyverse`,`reshape2`, `ineq`, `ggplot2`, `ggpubr`, `cowplot`

## Data Availability
The data associated with this project is available on Figshare:
- DOI: [10.6084/m9.figshare.27252282](https://doi.org/10.6084/m9.figshare.27252282)
- Download and place the contents in the `0_data/ext_data` directory

## Usage

### 0. Data Generation

#### Knowledge Metrics
Scripts in this group process Gene Ontology Annotation (GOA) files to extract knowledge metrics: term count, unique term count, information content.
* `loop_years.py` - Iterates through GOA files placed in `0_data/ext_data` folder to process each GOA file
* `extract_GOA.py` - Processes GAF files to calculate GO-based metrics. Information content metric is calculated using pivot year of 2022 (GO structure and GO term frequencies in 2022)
* `ia.py` - Helper script for Information Accretion calculations
* `process_GOA_years.py` - Combines temporal GO-based metrics for each GO aspect and experiment type (articles describing high-throughput or low-throughput experiments)

#### Interest Metrics
* `interest_data.R` - Processes protein count data from Knowledge Management Center, each protein is identified with STRING ID
* `process_interest_data.R` - Preprocesses data by calculating cumulative full publication equivalents (FPEs), filtering clinical assay proteins (retrieved from MarkerDB), resolving duplicated protein aliases using GO-based metrics

Knowledge and interest metrics are stored in folder `1_metrics` for further analysis.

### 1. Analyzing Metrics
Scripts in this folder are used to analyze and visualize the calculated knowledge (GO-based) and interest (number of FPEs) metrics of proteins.
* `1.1_all_metrics.R`: 
    * Compare GO-based metrics of proteins that receive annotations from high-throughput articles and low-throughput articles
    * Examine temporal changes in protein information content and number of FPEs, with an emphasis on most informative proteins
* `1.2_gini.R`: Gini coefficient calculations using knowledge and interest metrics, comparing proteins first annotated before 2013 with all proteins

### 2. Correlation & Disease Association
Scripts in this folder are used to examine correlations between knowledge and interest metrics.
* `2.1_knowledge_vs_interest.R `: Knowledge-interest correlations over time with Spearman correlation tests
* `2.2_disease_association.R`: Knowledge-interest correlations among proteins associated with diseases (retrieved from DISEASES database)

### 3. Curation Delay Analysis
Scripts for analyzing the time gap between publication and annotation:
* `3.1_extract_pubyear.py`: EDirect-based scripts for retrieving PubMed publication years of the publications recorded in GOA files
* `3.2_delay_analysis.R`: Analysis and visualization of temporal differences between publication and annotation dates

## Installation

1. Clone this repository
2. Install dependencies
3. Download data from Figshare and place in `0_data/ext_data`
4. For EDirect data retrieval:
```bash
sh -c "$(wget -q https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh -O -)"
export PATH=${HOME}/edirect:${PATH}
```
5. Create `4_figures` directory to store generated plots. Run the scripts based on the folder order (from 0 to 3).
