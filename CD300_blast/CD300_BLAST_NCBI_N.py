from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO
import os
import logging
import time

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Define input and output paths
fasta_file = "./query_sequences.fasta"
output_dir = "./blastn_results"
os.makedirs(output_dir, exist_ok=True)

# Load nucleotide sequences
nucleotide_sequences = list(SeqIO.parse(fasta_file, "fasta"))
logging.info(f"Loaded {len(nucleotide_sequences)} nucleotide sequences.")

# Set BLAST parameters
e_value_thresh = 1e-3 # Change to -10
pause_time = 5  # Pause time between BLAST requests

blast_results = []

# Function to add BLAST results to the list
def collect_blast_results(record_id, blast_records_list, e_value_thresh):
    for blast_record in blast_records_list:
        if blast_record.alignments:
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    # Define criteria for a "good match"
                    is_good_match = (hsp.expect <= e_value_thresh and 
                                     (hsp.identities / alignment.length * 100) > 80 and
                                     alignment.length >= 50)  # You can adjust the 50 bp minimum alignment length as needed
                    
                    match_status = "Good Match" if is_good_match else "Poor Match"
                    
                    # Append the result as a dictionary
                    blast_results.append({
                        'Query ID': record_id,
                        'Subject Title': alignment.title,
                        'Alignment Length': alignment.length,
                        'E-value': hsp.expect,
                        'Score': hsp.score,
                        'Percentage Identity': (hsp.identities / alignment.length * 100),
                        'Match Status': match_status
                    })

# Perform BLASTN search and process results
for i, record in enumerate(nucleotide_sequences, start=1):
    logging.info(f"BLASTN search for sequence {i}: {record.id}")
    try:
        result_handle = NCBIWWW.qblast("blastn", "refseq_rna", record.seq, 
                                       entrez_query="Homo sapiens[organism]", expect=e_value_thresh)
        blast_records_list = list(NCBIXML.parse(result_handle))
        
        if blast_records_list:
            collect_blast_results(record.id, blast_records_list, e_value_thresh)
            logging.info(f"Collected results for sequence {record.id}")
        else:
            logging.warning(f"No BLAST results for sequence {record.id}")
        
        time.sleep(pause_time)  # Pause between requests to avoid rate limits

    except Exception as e:
        logging.error(f"Error processing sequence {record.id}: {e}")

logging.info("BLASTN searches completed.")

# Create a pandas DataFrame from the blast_results list
df = pd.DataFrame(blast_results)

# Save the DataFrame to an Excel file
output_excel_file = os.path.join(output_dir, "blastn_results.xlsx")
df.to_excel(output_excel_file, index=False)

logging.info(f"Results saved to {output_excel_file}")