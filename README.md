# Vertebrate Protein BLAST Pipeline

## Overview
This pipeline automates the process of downloading, filtering, and analyzing vertebrate protein datasets to identify CD300 protein homologs across species using BLAST. Steps 1-4 were adapted from the [SIRP-Seeker-Pypeline](https://github.com/Rittika1/SIRP-Seeker-Pypeline).

## Pipeline Steps

### 1. Download Vertebrate Protein Datasets
```./datasets download genome taxon vertebrates --annotated --include protein --filename vertebrates_proteins.zip```

### 2. Ensure Only One Assembly Per Species
```./remove-multiple-assemblies-persp.py```

### 3. Unzip and Concatenate Protein Files
```unzip vertebrates_proteins.zip```

```cat ncbi_dataset/data/*/protein.faa >> vertebrate_proteins.faa```

### 4. Filter Selected Species
```python3 SIRP-Seeker-Pypeline/filteringproteinfiles.py vertebrate_proteins.faa vertebrate_proteins_cleaned.faa```

### 5. Remove Specific Species from the Protein File
```sbatch --mem=40G remove_species_from_list.sh vertebrate_proteins_cleaned.faa species_list.txt vertebrate_proteins_cleaned_filtered.faa```

### 6. Create BLAST Database
```makeblastdb -in vertebrate_proteins_cleaned_filtered.faa -dbtype prot```

### 7. Run BLAST with Human CD300 Proteins
Manually collect human CD300 proteins into `human_CD300_proteins.fasta` and run:

```blastp -query human_CD300_proteins.fasta -db vertebrate_proteins_cleaned.faa -out all_vertebrate_cd300_hits.out -evalue 1e-25 -outfmt 6```

### 8-10. Process BLAST Results
Steps 8, 9, and 10 have been combined into a single script:

```sbatch --mem=60G Process_BLASTP.sh [BLAST output file] [BLASTP Database] [query CD300 Fasta] [output csv file name]```

This script can potentially run overnight.

#### Individual Scripts Breakdown:
- **8. Convert BLAST Output to CSV:**
  ```python3 blast_to_csv.py```
- **9. Add Order and Family Information:**
  ```sbatch --mem=60G add_taxonomy.sh ```
  Uses the ITIS (Integrated Taxonomic Information System) API to add taxonomic details (may not work for all species).
- **10. Add CD300 Names and Isoforms to CSV:**
  ```sbatch --mem=60G CD300_add_name_isoform_to_csv.sh```

## Notes
- Ensure that all dependencies, including Python and BLAST+, are installed before running the pipeline.
- Adjust memory allocation (`--mem=XXG`) as needed for large datasets.

