import pandas as pd
import argparse
import os

def search_genes_in_file(file_path, genes, column_name, delimiter='\t'):
    """Search for genes in a large TSV file."""
    found_genes = []
    try:
        for chunk in pd.read_csv(file_path, delimiter=delimiter, chunksize=100000):
            if column_name in chunk.columns:
                matching_rows = chunk[chunk[column_name].isin(genes)]
                found_genes.append(matching_rows)
            else:
                print(f"Column '{column_name}' not found in {file_path}")
                return pd.DataFrame()
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        return pd.DataFrame()

    return pd.concat(found_genes) if found_genes else pd.DataFrame()

def main(gene_file, output_dir, clinical_variants, clinical_annotations, relationships, drug_labels):
    """Main function to search for genes in clinical data files."""
    os.makedirs(output_dir, exist_ok=True)  # Ensure output directory exists

    # Load gene list
    gene_data = pd.read_csv(gene_file)
    if "Gene" not in gene_data.columns:
        print("Error: 'Gene' column not found in input file.")
        return

    genes = gene_data["Gene"].unique()
    print(f"Searching for {len(genes)} unique genes...")

    # Search for genes in each file
    results = {
        "clinical_variants": search_genes_in_file(clinical_variants, genes, "gene"),
        "clinical_annotations": search_genes_in_file(clinical_annotations, genes, "Gene"),
        "relationships": search_genes_in_file(relationships, genes, "Entity1_name"),
        "drug_labels": search_genes_in_file(drug_labels, genes, "Gene Symbol"),
    }

    # Save results
    for name, df in results.items():
        if not df.empty:
            output_path = os.path.join(output_dir, f"{name}.csv")
            df.to_csv(output_path, index=False)
            print(f"{len(df)} rows saved to {output_path}")
        else:
            print(f"No matches found in {name}. Skipping file.")

    print("\n Search completed.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract gene-related clinical data from PharmGKB files.")
    parser.add_argument("gene_file", help="Path to CSV file containing gene names.")
    parser.add_argument("output_dir", help="Directory to save output files.")
    parser.add_argument("clinical_variants", help="Path to clinicalVariants.tsv file.")
    parser.add_argument("clinical_annotations", help="Path to clinical_annotations.tsv file.")
    parser.add_argument("relationships", help="Path to relationships.tsv file.")
    parser.add_argument("drug_labels", help="Path to drugLabels.byGene.tsv file.")

    args = parser.parse_args()
    main(args.gene_file, args.output_dir, args.clinical_variants, args.clinical_annotations, args.relationships, args.drug_labels)
