#!/bin/bash

# Check if the correct number of arguments are provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <fasta_file> <species_list.txt> <output_file>"
    exit 1
fi

fasta_file="$1"
species_list="$2"
output_file="$3"

# Create a temporary file to hold the headers we want to exclude
temp_file=$(mktemp)

# Read species_list.txt and create a pattern to exclude from the FASTA file
grep_pattern=$(awk '{printf("%s|", $0)}' "$species_list" | sed 's/|$//')

# Find headers in the fasta file that match any species from the list and store in temp file
grep -E "^>.*(${grep_pattern})" "$fasta_file" > "$temp_file"

# Use awk to extract all sequences NOT matching the excluded headers
awk 'BEGIN {RS=">"; ORS=""} NR==FNR {headers[$1]; next} !($1 in headers) {print ">"$0}' "$temp_file" "$fasta_file" > "$output_file"

# Remove the temporary file
rm "$temp_file"

echo "Sequences excluding species in $species_list have been saved to $output_file"
