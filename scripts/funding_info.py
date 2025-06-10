import time
import pandas as pd
import argparse
from Bio import Entrez

# Set email 
Entrez.email = "your_email@example.com"

def search_pubmed_for_gene_disease(gene_name, disease):
    """Searches PubMed for publications mentioning a gene and disease."""
    try:
        term = f"{gene_name} AND {disease}"
        print(f"🔍 Searching PubMed for: {term}")

        all_pubmed_ids = []
        retstart = 0
        batch_size = 1000

        while True:
            handle = Entrez.esearch(db="pubmed", term=term, retmax=batch_size, retstart=retstart, rettype="xml")
            record = Entrez.read(handle)
            pubmed_ids = record.get("IdList", [])
            all_pubmed_ids.extend(pubmed_ids)

            if len(pubmed_ids) < batch_size:
                break 

            retstart += batch_size  

        print(f"Found {len(all_pubmed_ids)} articles for {gene_name} ({disease})")
        return all_pubmed_ids
    except Exception as e:
        print(f"Error searching PubMed for {gene_name}: {e}")
        return []

def fetch_pubmed_funding(pubmed_id):
    """Fetches funding agency information from a PubMed article."""
    try:
        handle = Entrez.efetch(db="pubmed", id=str(pubmed_id), rettype="xml")
        records = Entrez.read(handle)
        funding_info = []

        for article in records.get("PubmedArticle", []):
            grant_list = article["MedlineCitation"]["Article"].get("GrantList", [])
            for grant in grant_list:
                funding_info.append(grant.get("Agency", "Unknown Agency"))

        return "; ".join(funding_info) if funding_info else "No funding info"
    except Exception as e:
        print(f"Error fetching funding for PubMed ID {pubmed_id}: {e}")
        return "Error"

def main(input_file, output_file, disease):
    """Main function to fetch PubMed IDs and funding data for genes."""
    df = pd.read_csv(input_file)
    
    # Add new columns
    df["PubMed IDs"] = ""
    df["Funding Information"] = ""

    for index, row in df.iterrows():
        gene_name = row["Gene"]
        print(f"\n Processing gene: {gene_name}")

        # Search PubMed
        pubmed_ids = search_pubmed_for_gene_disease(gene_name, disease)
        df.at[index, "PubMed IDs"] = "; ".join(pubmed_ids)

        all_funding_info = []
        
        # Fetch funding details
        for pubmed_id in pubmed_ids:
            funding_info = fetch_pubmed_funding(pubmed_id)
            all_funding_info.append(funding_info)
            time.sleep(0.5) 

        df.at[index, "Funding Information"] = "; ".join(all_funding_info)

        print(f" Finished processing {gene_name}: {len(pubmed_ids)} articles found.")

    # Save results
    df.to_csv(output_file, index=False)
    print(f"\n Results saved to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fetch PubMed IDs and funding information for genes.")
    parser.add_argument("input_file", help="Path to the input CSV file containing gene names.")
    parser.add_argument("output_file", help="Path to save the output CSV file.")
    parser.add_argument("--disease", required=True, help="Disease name to search for (e.g., 'Type 2 Diabetes').")

    args = parser.parse_args()
    main(args.input_file, args.output_file, args.disease)
