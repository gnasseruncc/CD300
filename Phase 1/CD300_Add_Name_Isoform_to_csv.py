import pandas as pd
from Bio import SeqIO
import re

# Load your CSV file into a pandas DataFrame
csv_file = "updated_all_vertebrate_cd300_hits.csv"
df = pd.read_csv(csv_file)

# Create a dictionary to store the relevant parts from the fasta file
fasta_file = "human_CD300_proteins.fasta"
fasta_dict = {}

# Parse the fasta file and extract relevant data
for record in SeqIO.parse(fasta_file, "fasta"):
    # Extract the ID and description
    seq_id = record.id
    description = record.description
    
    # Extract the gene name (second token after the sequence ID)
    gene_name = description.split(" ", 2)[1]
    
    # Extract the isoform information if present
    isoform_match = re.search(r'\[isoform=.*?\]', description)
    isoform_info = isoform_match.group(0) if isoform_match else ""
    
    # Combine the gene name and isoform info
    relevant_info = f"{gene_name} {isoform_info}".strip()  # This keeps both if present
    
    # Store it in the dictionary using the ID as the key
    fasta_dict[seq_id] = relevant_info

# Now match the Query ID in the DataFrame with the fasta IDs
def get_cd300_info(query_id):
    return fasta_dict.get(query_id, "")

# Apply the function to create a new column 'Query CD300'
df['Query CD300'] = df['Query ID'].apply(get_cd300_info)

# Save the updated DataFrame to a new CSV
output_csv_file = 'CD300_all_vertebrates_NCBI_BLAST_Final'  # replace with your desired output file name
df.to_csv(output_csv_file, index=False)

print(f"Updated CSV saved as {output_csv_file}")