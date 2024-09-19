import requests, re, subprocess, os

print("Current Working Directory:", os.getcwd())

def fetch_sequence(ensembl_accession):
    # Ensure the accession is formatted correctly for the Ensembl API
    url = f"https://rest.ensembl.org/sequence/id/{ensembl_accession}?type=seq&content-type=text/x-fasta"
    response = requests.get(url)

    if response.ok:
        return response.text
    else:
        print(f"Error fetching data for {ensembl_accession}: {response.status_code} - {response.text}")
        return None

def run_blast(query_sequence, db_name):
    # Create a temporary file for the query sequence
    with open("temp_query.fasta", "w") as temp_file:
        temp_file.write(query_sequence)

    # Run BLAST against the local database
    blast_command = [
        "blastn",
        "-query", "temp_query.fasta",
        "-db", db_name,
        "-out", "blast_results.txt",
        "-outfmt", "6",  # Tabular format for easier parsing
    ]
    
    # Execute the command
    subprocess.run(blast_command)

def fetch_query_sequences(query_file):
    with open(query_file, "r") as f:
        sequences = f.read().split('>')[1:]  # Split by FASTA header
    return [">" + seq.strip() for seq in sequences]

def parse_fasta_header(header):
    # Extract the accession number
    match = re.match(r">(NC_[0-9.]+)", header)
    if not match:
        return None
    
    accession = match.group(1)
    details = {}
    
    # Remove the accession part from the header for further processing
    clean_header = header[len(match.group(0)):].strip()
    
    # Extract various attributes using regex
    gene_match = re.search(r'\[gene=([^]]+)\]', clean_header)
    db_xref_match = re.search(r'\[db_xref=([^]]+)\]', clean_header)
    protein_match = re.search(r'\[protein=([^]]+)\]', clean_header)
    protein_id_match = re.search(r'\[protein_id=([^]]+)\]', clean_header)
    location_match = re.search(r'\[location=([^]]+)\]', clean_header)
    gbkey_match = re.search(r'\[gbkey=([^]]+)\]', clean_header)
    
    # Populate the details dictionary
    if gene_match:
        details['gene'] = gene_match.group(1)
    if db_xref_match:
        details['db_xref'] = db_xref_match.group(1).split(',')
    if protein_match:
        details['protein'] = protein_match.group(1)
    if protein_id_match:
        details['protein_id'] = protein_id_match.group(1)
    if location_match:
        details['location'] = location_match.group(1)
    if gbkey_match:
        details['gbkey'] = gbkey_match.group(1)

    return accession, details

def process_fasta_headers(fasta_file):
    fasta_dict = {}
    
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                accession, details = parse_fasta_header(line)
                if accession:
                    fasta_dict[accession] = details
    
    return fasta_dict

def main():
    pr_accession_numbers = [
    "GCA_003255815.1", "GCA_000951035.1", "GCA_001698545.1", 
    "GCA_000235385.1", "GCA_000258655.2", "GCA_000181295.3", 
    "GCA_000001515.5", "GCA_000956105.1", "GCA_011100615.1", 
    "GCA_000951045.1", "GCA_000146795.3", "GCA_000769185.1", 
    "GCA_000151905.3", "GCA_003258685.1", "GCA_000001405.29", 
    "GCA_000952055.2", "GCA_003339765.3", "GCA_000165445.3", 
    "GCA_008728515.1", "GCA_001604975.1", "GCA_000956065.1", 
    "GCA_000955945.1", "GCA_002880775.3", "GCA_000164805.2", 
    "GCA_002776525.2", "GCA_000409795.2", "GCA_011100555.1"
    ]
    query_file = "query_sequences.fasta"  # Input your query file here

    # Fetch query sequences from the FASTA file
    query_sequences = fetch_query_sequences(query_file)
    fasta_data = process_fasta_headers(query_file)

    # Print the organized dictionary
    for accession, details in fasta_data.items():
        print(f"{accession}: {details}")

    for accession in pr_accession_numbers:
        print(f"Processing genome: {accession}")

        # Fetch the genome sequence
        genome_sequence = fetch_sequence(accession)
        if not genome_sequence:
            continue  # Skip to the next accession if fetching failed

        # Create a temporary FASTA file for the genome
        with open("temp_genome.fasta", "w") as genome_file:
            genome_file.write(genome_sequence)

        # Create a BLAST database from the temporary genome file
        subprocess.run(["makeblastdb", "-in", "temp_genome.fasta", "-dbtype", "nucl", "-out", "temp_genome_db"])

        # Iterate through each query sequence and perform BLAST
        for query_sequence in query_sequences:
            run_blast(query_sequence, "temp_genome_db")
        
        # Clean up temporary files
        os.remove("temp_genome.fasta")
        os.remove("temp_genome_db.nhr")
        os.remove("temp_genome_db.nin")
        os.remove("temp_genome_db.nsq")
        os.remove("temp_query.fasta")  # Clean up the query temp file

    print("BLAST searches completed. Results saved in blast_results.txt.")

if __name__ == "__main__":
    main()
