import pandas as pd
import os

# Directory where protein prediction files are stored
residue_data_dir = snakemake.input['p2rank_dir']
pocket_data_dir = snakemake.input['fpocket_dir']

# Create an empty DataFrame to store merged data from all proteins
final_merged_data = pd.DataFrame()

# Loop through all files in the folder
for file in os.listdir(residue_data_dir):
    if file.endswith('_residues.csv'):  # Process only residue files
        # Extract the base name of the file (e.g., A0A023W421.pdb_residues.csv -> A0A023W421.pdb)
        base_name = file.replace('_residues.csv', '')

        # Define paths for residue-level and pocket-level files
        residue_file = os.path.join(residue_data_dir, file)
        pocket_file = os.path.join(pocket_data_dir, f"{base_name}_predictions.csv")

        # Check if the matching predictions file exists
        if os.path.exists(pocket_file):
            # Load the residue-level predictions from P2Rank
            residue_data = pd.read_csv(residue_file)

            # Load the pocket-level predictions from fpocket rescored with PRANK
            pocket_data = pd.read_csv(pocket_file)

            # Remove leading and trailing spaces from column names
            residue_data.columns = residue_data.columns.str.strip()
            pocket_data.columns = pocket_data.columns.str.strip()

            # Split the 'residue_ids' column into a list of residue identifiers
            pocket_data['residue_ids'] = pocket_data['residue_ids'].apply(lambda x: x.split())

            # Create a dictionary to map residue_ids to pocket names
            residue_to_pocket = {}
            for _, row in pocket_data.iterrows():
                for residue_id in row['residue_ids']:
                    residue_to_pocket[residue_id] = row['name']  # Store the pocket name


            # Define a function to get the pocket name for each residue
            def get_pocket_name(residue_label, chain):
                residue_id = f"{chain}_{residue_label}"  # Format residue_id as 'chain_residue_label'
                return residue_to_pocket.get(residue_id, None)  # Return None if no match is found


            # Add a new 'fpocket.prediction' column in the residue_data
            residue_data['fpocket.prediction'] = residue_data.apply(lambda row: get_pocket_name(row['residue_label'], row['chain']),
                                                        axis=1)

            # Merge the residue-level data with the pocket-level data on the 'pocket' column
            merged_data = pd.merge(residue_data, pocket_data, how='left', left_on='fpocket.prediction', right_on='name')

            # Drop the 'name' column from pocket_data as it's redundant after the merge
            merged_data = merged_data.drop(columns=['name'])

            # Define a dictionary mapping old column names to new column names
            column_rename_map = {
                'residue_label': 'residue_number',
                'score_x': 'residue_score',
                'probability_x': 'residue_probability',
                'score_y': 'pocket_score',
                'probability_y': 'pocket_probability',
            }

            # Rename columns in the DataFrame
            merged_data = merged_data.rename(columns=column_rename_map)

            # Add the UniProtKB ID as a new column (extracted from the base name of the protein)
            merged_data['UniProtKB'] = base_name.split('.')[0]  # Example: Extract UniProtKB ID from the filename

            # Append the merged data for this protein structure to the final DataFrame
            final_merged_data = pd.concat([final_merged_data, merged_data], ignore_index=True)

# Save the final merged data to a new CSV file
final_merged_data.to_csv(snakemake.output['merged_data'], index=False)


