# Import necessary modules
import os
import sys
import re  # Import the regular expressions module

# Function to convert a FASTA file to a dictionary
def convertfastatodict(filename):
    # Open the input FASTA file for reading
    fi = open(filename, 'r')
    
    # Initialize variables to store the header and sequence
    key = ''
    value = ''

    # Iterate through each line in the file
    for line in fi:
        # Remove leading and trailing whitespaces from the line
        line = line.strip()

        # Check if the line starts with ">"
        if line.startswith(">"):  # Indicates a header line
            # Replace spaces with underscores in the header
            line = re.sub(" ", "_", line)

            # If key is not empty, update the dictionary with the previous header and sequence
            if key != '':
                fastadict[key] = value

            # Find the species name enclosed in square brackets at the end of the line
            matches = re.findall(r'\[[^\]]*\]$', line)
            
            # If a match is found, extract the species name
            if matches:
                species_name = re.sub(' ', '_', matches[0][1:-1])
                # Extract descriptors from the line and construct the new key
                sequence_descriptors = str(line[1:-len(matches[0])]).split("_")
                key = species_name + "_" + "_".join(sequence_descriptors[0:-1])
            else:
                # If no match is found, use the entire line as the key (excluding the ">")
                key = line[1:]

            # Reset the value for the new header
            value = ''
            continue
        else:
            # If the line does not start with ">", append it to the current sequence
            value += line.strip()  # Append to dictionary

    # Close the input file
    fi.close()

# Function to print a dictionary of FASTA sequences to a new FASTA file
def PrintDictToFastaFile(dictionary, output):
    # Set the maximum length for each line in the output file
    n = 70
    
    # Open the output file for writing using a 'with' statement
    with open(output, 'w') as fo:
        # Iterate through key-value pairs in the dictionary
        for key, value in dictionary.items():
            # Strip leading and trailing whitespaces from the key
            key = key.strip()
            
            # Strip leading whitespaces from the value
            value = value.lstrip()
            
            # Split the value into chunks of length 'n'
            chunks = [value[i:i + n] for i in range(0, len(value), n)]
            
            # Write the modified header and sequence to the output file
            print('>' + key + '\n' + '\n'.join(chunks), file=fo)

# Check if the script is being run as the main program
if __name__ == '__main__':
    # Initialize an empty dictionary to store FASTA sequences
    fastadict = {}
    
    # Retrieve the input and output filenames from command line arguments
    fastafile = sys.argv[1]
    outputfile = sys.argv[2]

    # Call the function to convert the input FASTA file to a dictionary
    convertfastatodict(fastafile)

    # Call the function to print the dictionary to a new FASTA file
    PrintDictToFastaFile(fastadict, outputfile)


##--how to run
##--python3 filteringproteins.py proteins.faa filtered_proteins.faa
##--python3 scripts/filteringproteinfiles.py dataset/vertebrates_protein.faa dataset/vertebrates_protein_filtered.faa
##--python3 scripts/filteringproteinfiles.py dataset/amphibians_protein.faa dataset/amphibians_protein_filtered.faa
##-- python3 scripts/filteringproteinfiles.py dataset/mammals_protein.faa dataset/mammals_protein_filtered.faa
