# Uses NCBI API to grab human sequences from refseq_protein on NCBI and blast it against cd300 protein file

from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO
import os
import logging
import time

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Load protein sequences from FASTA file
fasta_file = "C:/Users/ibirc/OneDrive/Documents/Python Scripts/cd300h.faa" #load with your protein sequences
protein_sequences = list(SeqIO.parse(fasta_file, "fasta"))
logging.info(f"Loaded {len(protein_sequences)} protein sequences from {fasta_file}")

# Create a directory for output files if it doesn't exist
output_dir = "C:/Users/ibirc/OneDrive/Documents/Python Scripts/blastp_results"
os.makedirs(output_dir, exist_ok=True)

# Set E-value threshold
e_value_thresh = 1e-3

# Perform BLAST search for each protein sequence and write results to file
for i, record in enumerate(protein_sequences):
    logging.info(f"Performing BLASTP search for protein {i+1} ({record.id})")
    try:
        # Check if the result handle has any content
        # Use BLASTP to search protein sequences against the human protein database
        result_handle = NCBIWWW.qblast("blastp", "refseq_protein", record.seq, 
                                         entrez_query="Homo sapiens[organism]", #load with your organism of choice
                                         expect=e_value_thresh)

        # Parse the XML results directly from the result handle
        blast_records = NCBIXML.parse(result_handle)

        # Check if any records were returned
        blast_records_list = list(blast_records)
        if not blast_records_list:
            logging.error(f"No results returned from BLAST search for protein {i+1} ({record.id}).")
            continue  # Skip to the next protein

        # Write results to text file
        output_file = os.path.join(output_dir, f"CD300H_human_protein_blast_results_{i+1}.txt")
        with open(output_file, 'w') as f:
            f.write(f"BLASTP results for protein {i+1} ({record.id}) against human protein sequences:\n")
            f.write(f"E-value threshold: {e_value_thresh}\n\n")
            
            # Check if there are any alignments
            has_alignments = False
            alignment_count = 0
            for blast_record in blast_records_list:
                if len(blast_record.alignments) > 0:
                    has_alignments = True
                    for alignment in blast_record.alignments:
                        for hsp in alignment.hsps:
                            if hsp.expect <= e_value_thresh:
                                alignment_count += 1
                                f.write(f"Sequence: {alignment.title}\n")
                                f.write(f"Length: {alignment.length}\n")
                                f.write(f"E-value: {hsp.expect}\n")
                                f.write(f"Score: {hsp.score}\n")
                                f.write(f"Identities: {hsp.identities}/{alignment.length} "
                                        f"({hsp.identities/alignment.length*100:.2f}%)\n")
                                f.write(f"Query: {hsp.query[:50]}...\n")
                                f.write(f"Match: {hsp.match[:50]}...\n")
                                f.write(f"Sbjct: {hsp.sbjct[:50]}...\n")
                                f.write("\n")
            logging.info(f"Found {alignment_count} significant alignments for protein {i+1} ({record.id}).")

            # If no alignments were found
            if not has_alignments:
                f.write("No significant alignments found.\n")

        
        logging.info(f"Results written to {output_file}")
        
        # Pause between requests to avoid hitting NCBI rate limits
        time.sleep(5)
        
    except Exception as e:
        logging.error(f"Error processing protein {i+1}: {e}")

logging.info("BLASTP searches completed. Results have been saved in the 'blastp_results' directory.")
