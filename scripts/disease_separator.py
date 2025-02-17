# This file separates the diseases of interest from the GWAS file.

import csv
import argparse
import os

def filter_gwas_data(input_file, output_file, disease_name):
    """
    Filters GWAS data for a specified disease and saves the results.

    Parameters:
    - input_file (str): Path to the input GWAS Catalog TSV file.
    - output_file (str): Path to save the filtered output TSV file.
    - disease_name (str): Disease name to filter for.
    """
    # Define the columns of interest
    columns_of_interest = [
        "PUBMEDID",
        "REPORTED GENE(S)",
        "DISEASE/TRAIT",
        "INITIAL SAMPLE SIZE",
        "REPLICATION SAMPLE SIZE",
        "STRONGEST SNP-RISK ALLELE",
        "RISK ALLELE FREQUENCY",
        "P-VALUE",
        "OR or BETA",
        "95% CI (TEXT)"
    ]

    # Open the input TSV file and the output file
    with open(input_file, newline='', encoding='utf-8') as infile, open(output_file, 'w', newline='', encoding='utf-8') as outfile:
        reader = csv.DictReader(infile, delimiter='\t')
        writer = csv.DictWriter(outfile, fieldnames=columns_of_interest, delimiter='\t')

        # Write the header to the output file
        writer.writeheader()

        # Iterate through each row and write only the columns of interest for the specified disease
        for row in reader:
            if disease_name.lower() in row["DISEASE/TRAIT"].lower():  # Case-insensitive filtering
                writer.writerow({col: row[col] for col in columns_of_interest})

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter GWAS data for a specified disease.")
    parser.add_argument("input_file", help="Path to the input GWAS Catalog TSV file.")
    parser.add_argument("output_file", help="Path to save the filtered output TSV file.")
    parser.add_argument("--disease", required=True, help="Disease name to filter for (e.g., 'Type 2 Diabetes').")

    args = parser.parse_args()

    filter_gwas_data(args.input_file, args.output_file, args.disease)
