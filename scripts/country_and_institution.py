import pandas as pd
import argparse
import time
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
        print(f"Error fetching PubMed data for ID {pubmed_id}: {e}")
        return None

def search_pubmed_for_gene_disease(gene_name, disease):
    """Searches PubMed for publications mentioning a gene and disease."""
    try:
        term = f"{gene_name} AND {disease}"
        print(f"Searching PubMed for term: {term}")

        all_pubmed_ids = []
        retstart = 0
        batch_size = 1000

        while True:
            handle = Entrez.esearch(db="pubmed", term=term, retmax=batch_size, retstart=retstart, rettype="xml")
            record = Entrez.read(handle)
            pubmed_ids = record.get("IdList", [])
            all_pubmed_ids.extend(pubmed_ids)

            if len(pubmed_ids) < batch_size:
                break  # Stop when all records are retrieved

            retstart += batch_size

        print(f"Found {len(all_pubmed_ids)} articles for {gene_name} and {disease}")
        return all_pubmed_ids
    except Exception as e:
        print(f"Error searching PubMed for {gene_name} and {disease}: {e}")
        return []

def extract_institution_country(pubmed_metadata):
    """Extracts institutions and countries from PubMed metadata."""
    institutions = set()
    countries = set()

    try:
        for article in pubmed_metadata.get("PubmedArticle", []):
            author_list = article["MedlineCitation"]["Article"].get("AuthorList", [])
            for author in author_list:
                if "AffiliationInfo" in author:
                    for affiliation in author["AffiliationInfo"]:
                        if "Affiliation" in affiliation:
                            institution = affiliation["Affiliation"]
                            institutions.add(institution)
                            country = institution.split(",")[-1].strip()  # Extract last part as country
                            countries.add(country)
    except Exception as e:
        print(f"Error parsing metadata: {e}")

    return list(institutions), list(countries)

def main(input_file, output_file, disease):
    """Main function to fetch PubMed IDs, institutions, and countries for genes."""
    df = pd.read_csv(input_file)

    if "Gene" not in df.columns:
        print("Error: 'Gene' column missing in input file.")
        return

    df["PubMed IDs"] = ""
    df["Institutions"] = ""
    df["Countries"] = ""

    for index, row in df.iterrows():
        gene_name = row["Gene"]
        print(f"\nProcessing gene: {gene_name}")

        # Search PubMed
        pubmed_ids = search_pubmed_for_gene_disease(gene_name, disease)
        df.at[index, "PubMed IDs"] = "; ".join(pubmed_ids)

        all_institutions = set()
        all_countries = set()

        # Fetch metadata for each PubMed ID
        for pubmed_id in pubmed_ids:
            metadata = fetch_pubmed_data(pubmed_id)
            if metadata:
                institutions, countries = extract_institution_country(metadata)
                all_institutions.update(institutions)
                all_countries.update(countries)
            time.sleep(0.5)  # Avoid API rate limits

        # Store results
        df.at[index, "Institutions"] = "; ".join(sorted(all_institutions))
        df.at[index, "Countries"] = "; ".join(sorted(all_countries))

        print(f"Finished processing gene: {gene_name}")

    df.to_csv(output_file, index=False)
    print(f"\nResults saved to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fetch institutions and countries from PubMed metadata.")
    parser.add_argument("input_file", help="Path to the input CSV file containing gene names.")
    parser.add_argument("output_file", help="Path to save the output CSV file.")
    parser.add_argument("--disease", required=True, help="Disease name to search for (e.g., 'Type 2 Diabetes').")

    args = parser.parse_args()
    main(args.input_file, args.output_file, args.disease)
