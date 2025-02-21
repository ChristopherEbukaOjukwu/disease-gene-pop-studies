import requests
import pandas as pd
import argparse
from Bio import Entrez

# Set email (required by NCBI Entrez API)
Entrez.email = "your_email@example.com"

def fetch_pubmed_data(pubmed_id):
    """Fetches PubMed metadata for a given PubMed ID."""
    try:
        handle = Entrez.efetch(db="pubmed", id=str(pubmed_id), rettype="xml")
        records = Entrez.read(handle)
        return records
    except Exception as e:
        print(f"Error fetching PubMed metadata for {pubmed_id}: {e}")
        return None

def search_pubmed_for_gene_disease(gene_name, disease):
    """Searches PubMed for publications mentioning a gene and disease."""
    try:
        term = f"{gene_name} AND {disease}"
        handle = Entrez.esearch(db="pubmed", term=term, retmax=1000)
        record = Entrez.read(handle)
        return int(record.get("Count", 0))
    except Exception as e:
        print(f"Error searching PubMed for {gene_name}: {e}")
        return 0

def search_europe_pmc_for_gene_disease(gene_name, disease):
    """Searches Europe PMC for publications mentioning a gene and disease."""
    try:
        query = f"{gene_name} AND {disease}"
        url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/search?query={query}&resulttype=lite&format=json"
        response = requests.get(url)
        if response.status_code == 200:
            data = response.json()
            return data.get('hitCount', 0)
        else:
            return 0
    except Exception as e:
        print(f"Error searching Europe PMC for {gene_name}: {e}")
        return 0

def search_openalex_for_gene_disease(gene_name, disease):
    """Searches OpenAlex for publications mentioning a gene and disease."""
    query = f"{gene_name} AND {disease}"
    url = f"https://api.openalex.org/works?filter=fulltext.search:{query}&per-page=200"
    total_count = 0
    cursor = "*"

    try:
        while cursor:
            response = requests.get(url + f"&cursor={cursor}")
            if response.status_code == 200:
                data = response.json()
                total_count += len(data.get("results", []))
                cursor = data["meta"].get("next_cursor")  # Get next page cursor
                if not cursor:  # Stop if there's no next page
                    break
            else:
                print(f"Error {response.status_code} for {gene_name}")
                break
    except Exception as e:
        print(f"Error searching OpenAlex for {gene_name}: {e}")

    return total_count

def main(input_file, output_file, disease):
    """Main function to fetch publication counts for genes."""
    df = pd.read_csv(input_file)
    gene_pubmed_count = {}
    gene_europe_pmc_count = {}
    gene_openalex_count = {}

    for gene in df['Gene'].unique():
        print(f"Processing {gene}...")
        gene_pubmed_count[gene] = search_pubmed_for_gene_disease(gene, disease)
        gene_europe_pmc_count[gene] = search_europe_pmc_for_gene_disease(gene, disease)
        gene_openalex_count[gene] = search_openalex_for_gene_disease(gene, disease)

        print(f"{gene}: PubMed={gene_pubmed_count[gene]}, EuropePMC={gene_europe_pmc_count[gene]}, OpenAlex={gene_openalex_count[gene]}")

    # Create and save results
    final_gene_count_df = pd.DataFrame({
        'Gene': list(gene_pubmed_count.keys()),
        'PubMed Count': list(gene_pubmed_count.values()),
        'Europe PMC Count': list(gene_europe_pmc_count.values()),
        'OpenAlex Count': list(gene_openalex_count.values()),
    })

    final_gene_count_df.to_csv(output_file, index=False)
    print(f"\n The results have been saved to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fetch publication counts for genes from PubMed, Europe PMC, and OpenAlex.")
    parser.add_argument("input_file", help="Path to the input CSV file containing gene names.")
    parser.add_argument("output_file", help="Path to save the output CSV file.")
    parser.add_argument("--disease", required=True, help="Disease name to filter for (e.g., 'Type 2 Diabetes').")

    args = parser.parse_args()
    main(args.input_file, args.output_file, args.disease)
