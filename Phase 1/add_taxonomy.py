#!/bin/bash

python3 << EOF
import csv
import requests
import pandas as pd
import time

def get_taxonomy(genus, species):
    url = f"https://www.itis.gov/ITISWebService/jsonservice/searchByScientificName?srchKey={genus}%20{species}"
    
    try:
        response = requests.get(url)
        response.raise_for_status()  # Raises an error if the request failed
        data = response.json()
        
        # Ensure the 'scientificNames' key exists and is not empty
        if data and 'scientificNames' in data and len(data['scientificNames']) > 0:
            tsn = data['scientificNames'][0]['tsn']
            hierarchy_url = f"https://www.itis.gov/ITISWebService/jsonservice/getFullHierarchyFromTSN?tsn={tsn}"
            hierarchy_response = requests.get(hierarchy_url)
            hierarchy_response.raise_for_status()  # Raises an error if the request failed
            hierarchy_data = hierarchy_response.json()
            
            # Extract family and order from hierarchy list if available
            family = next((item['taxonName'] for item in hierarchy_data['hierarchyList'] if item['rankName'] == 'Family'), None)
            order = next((item['taxonName'] for item in hierarchy_data['hierarchyList'] if item['rankName'] == 'Order'), None)
            
            return family, order
        else:
            print(f"No results found for {genus} {species}")
    
    except (requests.RequestException, KeyError, TypeError) as e:
        print(f"Error retrieving taxonomy for {genus} {species}: {e}")
    
    return None, None

# Read your specific CSV file
df = pd.read_csv('all_vertebrate_cd300_hits.csv')

# Add new columns for family and order
df['Family'] = ''
df['Order'] = ''

# Process each row with rate limiting (to avoid API overload)
for index, row in df.iterrows():
    family, order = get_taxonomy(row['Genus'], row['Species'])
    df.at[index, 'Family'] = family
    df.at[index, 'Order'] = order
    
    # Add a sleep to prevent overwhelming the server
    time.sleep(0.5)

# Save the updated DataFrame to a new CSV file
df.to_csv('updated_all_vertebrate_cd300_hits.csv', index=False)
EOF

echo "Script completed!"
