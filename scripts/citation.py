import requests
import pandas as pd
import argparse
import time
from Bio import Entrez
from requests.exceptions import RequestException

# Set email for NCBI Entrez API (Required)
Entrez.email = "your_email@example.com"

def search_europe_pmc_citations(gene_name, disease):
    """Fetches total citation count from Europe PMC for a given gene-disease pair."""
    query = f'"{gene_name}" AND "{disease}"'
    page_size = 1000
    cursor_mark = '*'
    total_citations = 0

    while True:
        url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/search?query={query}&resultType=core&cursorMark={cursor_mark}&pageSize={page_size}&format=json"
        try:
            response = requests.get(url)
            response.raise_for_status()
            data = response.json()
            results = data.get('resultList', {}).get('result', [])

            if not results:
                break  # No more results

            total_citations += sum(article.get('citedByCount', 0) for article in results)

            next_cursor_mark = data.get('nextCursorMark')
            if cursor_mark == next_cursor_mark:
                break  # No new cursor, exit loop
            cursor_mark = next_cursor_mark
        except RequestException as e:
            print(f" Error fetching Europe PMC citations for {gene_name}: {e}")
            return 0

    return total_citations

def search_openalex_citations(gene_name, disease, max_retries=3, timeout=10):
    """Fetches total citation count from OpenAlex for a given gene-disease pair."""
    query = f"{gene_name} AND {disease}"
    url = f"https://api.openalex.org/works?filter=fulltext.search:{query}&per-page=200"
    total_citations = 0
    cursor = "*"

    while cursor:
        retries = 0
        while retries < max_retries:
            try:
                response = requests.get(url + f"&cursor={cursor}", timeout=timeout)
                response.raise_for_status()
                data = response.json()
                results = data.get("results", [])

                if results:
                    total_citations += sum(work.get("cited_by_count", 0) for work in results)

                cursor = data.get("meta", {}).get("next_cursor")
                if not cursor:
                    return total_citations
                break  # Exit retry loop on success

            except RequestException as e:
                retries += 1
                print(f" Retry {retries}/{max_retries} for {gene_name}: {e}")
                time.sleep(2)  # Wait before retrying

    return total_citations

def main(input_file, output_file, disease):
    """Main function to fetch citation counts from Europe PMC and OpenAlex."""
    df = pd.read_csv(input_file)

    if "Gene" not in df.columns or "Population" not in df.columns:
        print(" Error: 'Gene' or 'Population' column missing in input file.")
        return

    gene_europe_pmc_citation_count = {}
    gene_openalex_citation_count = {}

    for index, row in df.iterrows():
        gene = row["Gene"]
        population = row["Population"]

        print(f"\n🔍 Processing {gene} (Population: {population})...")
        europe_pmc_citation_count = search_europe_pmc_citations(gene, disease)
        openalex_citation_count = search_openalex_citations(gene, disease)

        gene_europe_pmc_citation_count[gene] = europe_pmc_citation_count
        gene_openalex_citation_count[gene] = openalex_citation_count

        print(f"{gene} done. Europe PMC: {europe_pmc_citation_count}, OpenAlex: {openalex_citation_count}")

    # Create and save results
    final_gene_citation_count_df = pd.DataFrame({
        "Gene": df["Gene"],
        "Population": df["Population"],
        "Europe PMC Citations": [gene_europe_pmc_citation_count.get(gene, 0) for gene in df["Gene"]],
        "OpenAlex Citations": [gene_openalex_citation_count.get(gene, 0) for gene in df["Gene"]],
    })

    final_gene_citation_count_df.to_csv(output_file, index=False)
    print(f"\n The results have been saved to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fetch citation counts for genes from Europe PMC and OpenAlex.")
    parser.add_argument("input_file", help="Path to the input CSV file containing gene names and populations.")
    parser.add_argument("output_file", help="Path to save the output CSV file.")
    parser.add_argument("--disease", required=True, help="Disease name to search for (e.g., 'Type 2 Diabetes').")

    args = parser.parse_args()
    main(args.input_file, args.output_file, args.disease)
