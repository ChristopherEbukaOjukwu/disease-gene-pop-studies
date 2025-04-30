# Overview
This repository houses the scripts to ascertain how diseases and genes are studied in different populations.

## Dependencies
Ensure you have Python installed, then install the required dependencies: 

- pip install argparse
- pip install requests 
- biopython 
- pandas 
- argparse

## Data
The data is publicly available:

1. NHGRI-EBI GWAS Catalog (v1.0): https://www.ebi.ac.uk/gwas/docs/file-downloads.
2. OpenAlex API: https://api.openalex.org/works.
3. EuropePMC: https://europepmc.org/RestfulWebService.
4. PubMed: https://pubmed.ncbi.nlm.nih.gov.
5. Pharmacogenomics: https://www.pharmgkb.org/downloads.
6. Clinical studies: https://clinicaltrials.gov/data-api/api

## Usage
1. Separate the disease of interest from the GWAS file: python disease_separator.py gwas_catalog_v1.0.tsv filtered_output.tsv --disease "Disease"
2. Get publication count: python publication_count.py genes.csv publication_counts.csv --disease "Disease"
3. Funding: python fetch_funding_info.py genes.csv funding_results.csv --disease "Disease"
4. For pharmacogenomics, download CinicalVariants, ClinicalAnnotations, Relationships, and DrugLabels, then run: python pharmacogenomics.py genes.csv results/ pharmgkb_files/clinicalVariants.tsv pharmgkb_files/clinical_annotations.tsv pharmgkb_files/relationships.tsv pharmgkb_files/drugLabels.byGene.tsv
5. Citation: python citation.py genes.csv citations_results.csv --disease "Disease"
6. Country and institution: python countries_institutions.py genes.csv country_institution_results.csv --disease "Disease"
7. Clinical studies: python fetch_clinical_studies.py genes.csv clinical_study_results.csv --disease "Disease"
