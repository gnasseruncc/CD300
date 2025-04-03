import pandas as pd
import os
import re

# File path
data_file = 'insert/your/csv/file/here'

# Function to safely load CSV file
def load_csv(file_path):
    if not os.path.exists(file_path):
        print(f"Error: File not found - {file_path}")
        return None
    return pd.read_csv(file_path)

# Load the CSV file
df = load_csv(data_file)

# Check if DataFrame was loaded successfully
if df is None:
    print("Error: Unable to proceed due to missing file")
else:
    # Extract accession number function
    def extract_accession_number(accession):
        match = re.search(r'([XN]P_\d+\.\d+)', accession)
        return match.group(1) if match else accession

    # Apply extraction to the DataFrame
    df['Extracted_Accession'] = df['subject_id'].apply(extract_accession_number)

    # Create a dictionary mapping extracted accessions to BioProject_ID
    bioproject_dict = dict(zip(df['Extracted_Accession'], df['BioProject_ID']))

    # Fill in missing BioProject_ID
    df['BioProject_ID'] = df.apply(
        lambda row: bioproject_dict.get(row['Extracted_Accession'], row['BioProject_ID']) 
        if pd.isna(row['BioProject_ID']) else row['BioProject_ID'], 
        axis=1
    )

    # Remove the 'Extracted_Accession' column
    df = df.drop('Extracted_Accession', axis=1)

    # Save the updated DataFrame to a new CSV file
    output_file = '/insert/output/csv'
    df.to_csv(output_file, index=False)
    print(f"Updated file saved as: {output_file}")
