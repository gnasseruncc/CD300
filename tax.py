import requests, sys

"""
How to use this script. run python3 tax.py, it will prompt you repeatedly until you type exit. Each prompt asks for species
name in the format of "Homo sapiens" make sure the first part of the name is capitol, it will return a link where you find fill taxonomy
including order and sub order in the bottom of the right hand box on the taxonomy website
"""

def get_species_taxonomy(species_name):
    # Query Ensembl API for species taxonomy info
    url = f"https://rest.ensembl.org/taxonomy/name/{species_name}?content-type=application/json"
    headers = {"Content-Type": "application/json"}
    
    response = requests.get(url, headers=headers)
    
    if not response.ok:
        print(f"Failed to retrieve data: {response.status_code}")
        return None

    # Parse response as a list of dictionaries
    species_data = response.json()
    
    # Find the matching species in the response
    for species in species_data:
        if species.get("scientific_name", "").lower() == species_name.lower():
            # Return only the NCBI taxonomy ID
            return species.get("id")
    
    return None

def get_ncbi_lineage(taxonomy_id):
    # Query NCBI Taxonomy API for lineage using the Taxonomy ID
    url = f"https://www.ncbi.nlm.nih.gov/datasets/taxonomy/{taxonomy_id}/"
    print(url)

def main():
    while True:
        user_input = input("Enter a string (type 'exit' to quit): ")
        if user_input.lower() == "exit":
            print("Exiting...")
            break
        else:
            taxid = get_species_taxonomy(user_input)
            get_ncbi_lineage(taxid)

if __name__ == "__main__":
    main()
