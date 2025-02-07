#!/bin/bash

# Make sure to provide the following arguments when running this script:
# $1 - BLAST output file
# $2 - BLAST database
# $3 - Query CD300 proteins FASTA file 
# $4 - Output CSV file

if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <blast_output> <blast_database> <query_cd300_fasta> <output_csv>"
    exit 1
fi

# Assigning command-line arguments to variables
blast_output="$1"
fasta_db="$2"
query_cd300_fasta="$3"
output_csv="$4"

python3 << EOF
import os
import csv
from Bio import SeqIO
import re
import requests
import pandas as pd
import time

# Step 1: Parse BLAST output
def parse_blast(blast_file):
    all_hits = []
    
    with open(blast_file, 'r') as f:
        for line in f:
            cols = line.strip().split("\t")
            query_id = cols[0]   # CD300 protein ID (e.g., CD300A, CD300B)
            subject_id = cols[1]  # Full subject ID
            
            # Extract genus and species from subject_id (assuming format: Genus_species_X)
            subject_info = subject_id.split("_")
            genus = subject_info[0]
            species = subject_info[1]
            
            # Store hit information
            all_hits.append({
                'query_id': query_id,
                'genus': genus,
                'species': species,
                'subject_id': subject_id,
                'percent_identity': cols[2],
                'alignment_length': cols[3],
                'evalue': cols[10],
                'bit_score': cols[11],
            })
    
    return all_hits

# Step 2: Get protein sequences from the FASTA file
def get_sequences(fasta_file):
    sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        header_id = record.id  # Full header ID
        sequences[header_id] = str(record.seq)
    return sequences

# Step 3: Get taxonomy information from ITIS
def get_taxonomy(genus, species):
    url = f"https://www.itis.gov/ITISWebService/jsonservice/searchByScientificName?srchKey={genus}%20{species}"
    try:
        response = requests.get(url)
        response.raise_for_status()
        data = response.json()
        if data and 'scientificNames' in data and len(data['scientificNames']) > 0:
            tsn = data['scientificNames'][0]['tsn']
            hierarchy_url = f"https://www.itis.gov/ITISWebService/jsonservice/getFullHierarchyFromTSN?tsn={tsn}"
            hierarchy_response = requests.get(hierarchy_url)
            hierarchy_response.raise_for_status()
            hierarchy_data = hierarchy_response.json()
            family = next((item['taxonName'] for item in hierarchy_data['hierarchyList'] if item['rankName'] == 'Family'), None)
            order = next((item['taxonName'] for item in hierarchy_data['hierarchyList'] if item['rankName'] == 'Order'), None)
            return family, order
    except (requests.RequestException, KeyError, TypeError) as e:
        print(f"Error retrieving taxonomy for {genus} {species}: {e}")
    return None, None

# Step 4: Parse the query CD300 FASTA file to get relevant query information
def parse_query_cd300_fasta(query_cd300_fasta):
    fasta_dict = {}
    for record in SeqIO.parse(query_cd300_fasta, "fasta"):
        seq_id = record.id
        description = record.description
        gene_name = description.split(" ", 2)[1]
        isoform_match = re.search(r'\[isoform=.*?\]', description)
        isoform_info = isoform_match.group(0) if isoform_match else ""
        relevant_info = f"{gene_name} {isoform_info}".strip()
        fasta_dict[seq_id] = relevant_info
    return fasta_dict

# Step 5: Main function to process BLAST results, taxonomy, and query CD300 information
def process_blast_results(blast_file, fasta_file, query_cd300_fasta, output_csv):
    all_hits = parse_blast(blast_file)
    sequences = get_sequences(fasta_file)
    query_cd300_info = parse_query_cd300_fasta(query_cd300_fasta)

    # Create dataframe from BLAST hits
    df = pd.DataFrame(all_hits)

    # Add protein sequences to dataframe
    df['Protein Sequence'] = df['subject_id'].apply(lambda x: sequences.get(x, "Sequence not found"))

    # Add taxonomy info (Family and Order)
    df['Family'], df['Order'] = zip(*df.apply(lambda row: get_taxonomy(row['genus'], row['species']), axis=1))

    # Add query CD300 information
    df['Query CD300'] = df['query_id'].apply(lambda x: query_cd300_info.get(x, ""))

    # Save the final dataframe to CSV
    df.to_csv(output_csv, index=False)
    print(f"CSV file '{output_csv}' generated successfully.")

# Execute the main function
process_blast_results("$blast_output", "$fasta_db", "$query_cd300_fasta", "$output_csv")

EOF
