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

# Function to write BLAST results to file
def write_blast_results(output_file, record_id, blast_records_list, e_value_thresh):
    with open(output_file, 'w') as f:
        f.write(f"BLASTN results for {record_id}:\nE-value threshold: {e_value_thresh}\n\n")
        for blast_record in blast_records_list:
            if blast_record.alignments:
                for alignment in blast_record.alignments:
                    for hsp in alignment.hsps:
                        if hsp.expect <= e_value_thresh:
                            f.write(f"Sequence: {alignment.title}\nLength: {alignment.length}\n"
                                    f"E-value: {hsp.expect}\nScore: {hsp.score}\n"
                                    f"Identities: {hsp.identities}/{alignment.length} "
                                    f"({hsp.identities / alignment.length * 100:.2f}%)\n"
                                    f"Query: {hsp.query[:50]}...\nMatch: {hsp.match[:50]}...\n"
                                    f"Sbjct: {hsp.sbjct[:50]}...\n\n")
            else:
                f.write("No significant alignments found.\n")

# Perform BLASTN search and process results
for i, record in enumerate(nucleotide_sequences, start=1):
    logging.info(f"BLASTN search for sequence {i}: {record.id}")
    try:
        result_handle = NCBIWWW.qblast("blastn", "refseq_rna", record.seq, 
                                       entrez_query="Homo sapiens[organism]", expect=e_value_thresh)
        blast_records_list = list(NCBIXML.parse(result_handle))
        
        if blast_records_list:
            output_file = os.path.join(output_dir, f"CD300H_human_nucleotide_blast_results_{i}.txt")
            write_blast_results(output_file, record.id, blast_records_list, e_value_thresh)
            logging.info(f"Results saved to {output_file}")
        else:
            logging.warning(f"No BLAST results for sequence {record.id}")
        
        time.sleep(pause_time)  # Pause between requests to avoid rate limits

    except Exception as e:
        logging.error(f"Error processing sequence {record.id}: {e}")

logging.info("BLASTN searches completed. Results saved in 'blastn_results' directory.")
