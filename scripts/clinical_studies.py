import pandas as pd
import requests
import argparse
import csv

def get_clinical_studies_count_and_countries(gene, condition_filters):
    """Fetches the number of clinical studies and involved countries for a given gene."""
    url = "https://clinicaltrials.gov/api/v2/studies"
    page_size = 100
    page_token = None
    total_studies = 0
    filtered_studies = 0
    country_data = []

    while True:
        params = {
            "query.term": gene,
            "format": "json",
            "pageSize": page_size,
            "pageToken": page_token
        }

        response = requests.get(url, params=params)
        if response.status_code == 200:
            data = response.json()
            studies = data.get("studies", [])
            total_studies += len(studies)

            for study in studies:
                conditions = study.get("protocolSection", {}).get("conditionsModule", {}).get("conditions", [])
                meshes = study.get("derivedSection", {}).get("conditionBrowseModule", {}).get("meshes", [])
                mesh_terms = [mesh["term"] for mesh in meshes]
                structured_terms = set(conditions + mesh_terms)

                # Match condition using structured fields
                if any(condition.lower() in (term.lower() for term in structured_terms) for condition in condition_filters):
                    filtered_studies += 1
                    locations = study.get("protocolSection", {}).get("contactsLocationsModule", {}).get("locations", [])
                    country_data.extend(location.get("country", "Unknown") for location in locations if location.get("country"))

            page_token = data.get("nextPageToken")
            if not page_token:
                break
        else:
            print(f"Error: {response.status_code} for gene {gene}")
            print("Details:", response.text)
            return 0, []

    return filtered_studies, sorted(set(country_data))

def main(input_file, output_file, disease):
    """Main function to fetch clinical study counts and country data for genes."""
    df = pd.read_csv(input_file)

    if "Gene" not in df.columns or "Population" not in df.columns:
        print("Error: 'Gene' or 'Population' column missing in input file.")
        return

    condition_filters = [disease.lower(), disease.capitalize(), disease.upper()]

    study_count_data = []
    for index, row in df.iterrows():
        gene = row["Gene"]
        population = row["Population"]
        print(f"Searching clinical studies for {gene}...")

        study_count, countries = get_clinical_studies_count_and_countries(gene, condition_filters)
        study_count_data.append([gene, population, study_count, ", ".join(countries)])

    # Save results
    with open(output_file, mode="w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(["Gene", "Population", "Number of Clinical Studies", "Countries"])
        writer.writerows(study_count_data)

    print(f"Results saved to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fetch clinical study counts and associated countries from ClinicalTrials.gov.")
    parser.add_argument("input_file", help="Path to the input CSV file containing gene names and populations.")
    parser.add_argument("output_file", help="Path to save the output CSV file.")
    parser.add_argument("--disease", required=True, help="Disease name to search for (e.g., 'Type 2 Diabetes').")

    args = parser.parse_args()
    main(args.input_file, args.output_file, args.disease)
