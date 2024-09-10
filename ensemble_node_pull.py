import requests
import sys
import re

server = "https://rest.ensembl.org"
ext = "/info/genomes/taxonomy/Euarchontoglires?"

r = requests.get(server + ext, headers={"Content-Type": "application/json"})

if not r.ok:
    r.raise_for_status()
    sys.exit()

decoded = r.json()
output = repr(decoded)

print(output)
print(type(output))

# Step 1: Remove the outer brackets
inner_string = output.strip('[]')

# Step 2: Split the string into individual dictionary-like substrings
dict_strings = re.split(r'\s*,\s*(?=\{)', inner_string)

# Initialize a list to hold the extracted data
extracted_data = []

# Step 3: Extract fields from each dictionary-like string
for dict_string in dict_strings:
    # Extract fields using regex
    scientific_name = re.search(r"'scientific_name':\s*'([^']*)'", dict_string)
    display_name = re.search(r"'display_name':\s*'([^']*)'", dict_string)
    assembly_accession = re.search(r"'assembly_accession':\s*'([^']*)'", dict_string)
    base_count = re.search(r"'base_count':\s*'({[0-9]+:[0-9]+})'", dict_string)
    
    # Construct a dictionary with the extracted fields
    data = {
        'scientific_name': scientific_name.group(1) if scientific_name else None,
        'display_name': display_name.group(1) if display_name else None,
        'assembly_accession': assembly_accession.group(1) if assembly_accession else None,
        'base_count': base_count.group(1) if base_count else None
    }
    
    extracted_data.append(data)

# Print the results
for data in extracted_data:
    print(data)
