import os
import csv
from Bio import SeqIO

# Input files
blast_output = "cd300_hits.out"
fasta_db = "vertebrate_proteins_cleaned.faa"  # Your protein sequence file
output_csv = "all_vertebrate_cd300_hits.csv"

# Function to parse BLAST output and store all hits for each CD300 type per species
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

# Function to retrieve protein sequences from the FASTA file
def get_sequences(fasta_file):
    sequences = {}
    
    # Parse the FASTA file and store sequences in a dictionary with header as key
    for record in SeqIO.parse(fasta_file, "fasta"):
        header_id = record.id  # This gets the full header ID (without the '>')
        sequences[header_id] = str(record.seq)
    
    return sequences

# Main function to parse BLAST results, get sequences, and write to CSV
def blast_to_csv(blast_file, fasta_file, output_file):
    # Parse BLAST output and extract all hits
    all_hits = parse_blast(blast_file)
    
    # Retrieve sequences from FASTA database
    sequences = get_sequences(fasta_file)
    
    # Open CSV file for writing
    with open(output_file, 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        
        # Write CSV header
        csvwriter.writerow(['Query ID', 'Genus', 'Species', 'Subject ID', 
                            'Percent Identity', 'Alignment Length', 
                            'E-value', 'Bit Score', 'Protein Sequence'])
        
        # Write all hits with sequences to the CSV file
        for hit in all_hits:
            subject_id = hit['subject_id']
            
            # Match the subject ID to the header in the FASTA file
            # Assuming headers in the FASTA file use the same format as subject IDs
            # For example: "Homo_sapiens_X" in BLAST output should match ">Homo_sapiens_X" in the FASTA file
            sequence = sequences.get(subject_id, "Sequence not found")
            
            # Write row to CSV
            csvwriter.writerow([hit['query_id'], hit['genus'], hit['species'], hit['subject_id'],
                                hit['percent_identity'], hit['alignment_length'], hit['evalue'], 
                                hit['bit_score'], sequence])

# Run the main function to generate the CSV
blast_to_csv(blast_output, fasta_db, output_csv)

print(f"CSV file '{output_csv}' generated successfully.")
