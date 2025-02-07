# CD300 Analyzation Pipeline
The goal of this project is to analyze ortholog/paralog functional evolution of the CD300 immune gene family across all major vertebrate lineages.


## Phase 1: Generate CD300 Verebrate Dataset
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
```python3 filteringproteinfiles.py vertebrate_proteins.faa vertebrate_proteins_cleaned.faa```

### 5. Remove Specific Species from the Protein File
The species we removed from this dataset were already searched for their CD300 sequences manually on Ensembl, so there is no reason to search  for them again and can be excluded from this specific pipeline. 

```sbatch --mem=40G remove_species_from_list.sh vertebrate_proteins_cleaned.faa species_list.txt vertebrate_proteins_cleaned_filtered.faa```

### 6. Create BLAST Database
```makeblastdb -in vertebrate_proteins_cleaned_filtered.faa -dbtype prot```

### 7. Run BLAST with Human CD300 Proteins
Manually collected human CD300 proteins into `human_CD300_proteins.fasta` and run:

```blastp -query human_CD300_proteins.fasta -db vertebrate_proteins_cleaned.faa -out all_vertebrate_cd300_hits.out -evalue 1e-25 -outfmt 6```

You can manually collect other species' CD300 sequences and BLAST them when you want to search for other clades. 

### 8-10. Process BLAST Results into Properly Formatted CSV
Steps 8, 9, and 10 have been combined into a single script:

```sbatch --mem=60G Process_BLASTP.sh [BLAST output file] [BLASTP Database] [query CD300 Fasta] [output csv file name]```

This script can potentially run overnight.

#### Breakdown of individual scripts:
- **8. Convert BLAST Output to CSV:**
  ```python3 blast_to_csv.py```
- **9. Add Order and Family Information:**
  ```sbatch --mem=60G add_taxonomy.py ```
  Uses the ITIS (Integrated Taxonomic Information System) API to add taxonomic details (may not work for all species).
- **10. Add CD300 Names and Isoforms to CSV:**
  ```sbatch --mem=60G CD300_add_name_isoform_to_csv.py```


  ## Phase 2: Check CD300 Dataset
  Now that you have your CSV file with hopefully a lot of putative CD300 sequences, you now need to actually check to see if they are in fact CD300s. The easiest way to do this is by creating phylogenetic trees and BLASTing your results back to humans (or to another species with solid analysis of CD300 gene structures - like chickens or dogs).


## Notes
- Ensure that all dependencies, including Python and BLAST+, are installed before running the pipeline.
- Adjust memory allocation (`--mem=XXG`) as needed for large datasets.

