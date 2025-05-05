# Population Disparities in Genetic Research: Code & Data

This repository contains all code and data used in the study:  
**"The Global Focus on European-Associated Genes Across Biomedical Research"**  
Authors: Christopher Ebuka Ojukwu, Daniel Acuna

---

## Overview

This project investigates how genetic research attention (measured via publications, citations, funding, clinical studies, and pharmacogenomic annotations) varies across populations for disease-associated genes identified in GWAS studies.

The analysis focuses on five diseases: Type 2 Diabetes, Breast Cancer, Prostate Cancer, Colorectal Cancer, and Schizophrenia.

---

## Dependencies

- `argparse`
- `requests`
- `biopython`
- `pandas`

---

## Data Sources

All data used are publicly available:

| Source | Description |
|--------|-------------|
| [GWAS Catalog (v1.0)](https://www.ebi.ac.uk/gwas/docs/file-downloads) | Disease-gene associations |
| [OpenAlex](https://api.openalex.org/works) | Publication + citation counts |
| [EuropePMC](https://europepmc.org/RestfulWebService) | Secondary publication metrics |
| [PubMed](https://pubmed.ncbi.nlm.nih.gov) | NIH-linked data |
| [PharmGKB](https://www.pharmgkb.org/downloads) | Pharmacogenomics annotations |
| [ClinicalTrials.gov](https://clinicaltrials.gov/data-api/api) | Clinical study counts |

---

## Usage

All analyses begin with a filtered GWAS file that lists disease–gene–population triplets. Example steps:

### 1. **Filter by Disease**
python disease_separator.py gwas_catalog_v1.0.tsv filtered_output.tsv --disease "Type 2 Diabetes"
### 2. **Get publication counts**
python publication_count.py genes.csv publication_counts.csv --disease "Type 2 Diabetes"
### 3. **Get citation metrics**
python citation.py genes.csv citations_results.csv --disease "Type 2 Diabetes"
### 4. **Fetch funding information**
python fetch_funding_info.py genes.csv funding_results.csv --disease "Type 2 Diabetes"
### 5. **Get country and institution data**
python countries_institutions.py genes.csv country_institution_results.csv --disease "Type 2 Diabetes"
### 6. **Retrieve clinical studies**
python fetch_clinical_studies.py genes.csv clinical_study_results.csv --disease "Type 2 Diabetes"
### 7. **Run pharmacogenomics analysis**
Download the following PharmGKB files first:
clinicalVariants.tsv
clinical_annotations.tsv
relationships.tsv
drugLabels.byGene.tsv
python pharmacogenomics.py genes.csv results/ pharmgkb_files/clinicalVariants.tsv \
pharmgkb_files/clinical_annotations.tsv pharmgkb_files/relationships.tsv \
pharmgkb_files/drugLabels.byGene.tsv
