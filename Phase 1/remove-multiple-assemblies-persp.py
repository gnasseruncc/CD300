'''
This code will comb through folders and subfolders, and get a list of the GCA acession code and the species name. 

Structure of the folder:
ncbi_dataset -> data -> multiple folders named as GC...
inside every GCA... folder is a file named as protein.faa
read the sequence header of the protein.faa, esp the first line. and see the pattern. 
>XP_034952902.1 protein Hook homolog 3 isoform X1 [Zootoca vivipara]
in this line 'Zootoca vivipara' is the species name. 
so the output should be folddername, species name

'''
import os
import re

def extract_species_name(file_path):
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('>'):
                matches = re.findall(r'\[([^\]]+)\]', line)
                if matches:
                    return matches[-1]
                break
    return None

def concatenate_proteins(folders, megafile_path, desired_species):
    with open(megafile_path, 'w') as megafile:
        for folder in folders:
            protein_file_path = os.path.join(folder, 'protein.faa')
            if os.path.exists(protein_file_path):
                with open(protein_file_path, 'r') as protein_file:
                    for line in protein_file:
                        if line.startswith('>'):
                            species_name = re.findall(r'\[([^\]]+)\]', line)[-1]
                            if species_name not in desired_species:
                                continue
                        megafile.write(line)

def read_desired_species(file_path):
    with open(file_path, 'r') as file:
        return [line.strip().replace('_', ' ') for line in file]

def main():
    base_path = '/purple/SIRP_newdata/ncbi_dataset/data'
    output = []
    species_to_folder = {}
    desired_species_file = 'species_list.txt'
    desired_species = read_desired_species(desired_species_file)

    for root, dirs, files in os.walk(base_path):
        for file in files:
            if file == 'protein.faa':
                folder_name = os.path.basename(root)
                species_name = extract_species_name(os.path.join(root, file))
                if species_name and species_name in desired_species and species_name not in species_to_folder:
                    species_to_folder[species_name] = folder_name
                    output.append((folder_name, species_name))

    # Print the results
    for folder, species in output:
        print(f'{folder}, {species}')

    # Concatenate protein.faa files into a megafile
    megafile_path = 'vertebrate_proteins.faa'
    concatenate_proteins([os.path.join(base_path, folder) for folder, _ in output], megafile_path, desired_species)
    print(f'All protein.faa files for the desired species have been concatenated into {megafile_path}')

if __name__ == '__main__':
    main()


##--usage
##-- python3 SIRP-Seeker-Pypeline/remove-multiple-assemblies-persp.py
#  root_dir = '/purple/SIRP_newdata/ncbi_dataset/data'
