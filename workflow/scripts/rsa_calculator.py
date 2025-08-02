import os
import subprocess
from Bio.PDB import make_dssp_dict

# Function to calculate RSA (normalized ASA)
def calculate_rsa(asa, residue_type):
    """Calculate RSA as ASA divided by the max ASA for a given residue type."""
    max_asa = {
        'A': 121.0, 'R': 265.0, 'N': 187.0, 'D': 187.0, 'C': 148.0, 'E': 214.0,
        'Q': 214.0, 'G': 97.0, 'H': 216.0, 'I': 195.0, 'L': 191.0, 'K': 230.0,
        'M': 203.0, 'F': 228.0, 'P': 154.0, 'S': 143.0, 'T': 163.0, 'W': 264.0,
        'Y': 255.0, 'V': 165.0
    }
    return asa / max_asa[residue_type] if residue_type in max_asa else None

# Function to extract RSA and secondary structure features from a DSSP file
def process_dssp_file(dssp_file, uniprot_id):
    """Extract RSA and secondary structure from a DSSP file."""
    dssp_data, _ = make_dssp_dict(dssp_file)
    residue_data = []

    for key, dssp_info in dssp_data.items():
        chain_id, res_id = key
        residue_id = res_id[1]
        residue_type = dssp_info[0]
        secondary_structure = dssp_info[1]
        asa = dssp_info[2]
        phi = dssp_info[3]
        psi = dssp_info[4]
        rsa = calculate_rsa(asa, residue_type)

        residue_data.append({
            'uniprot_id': uniprot_id,
            'chain_id': chain_id,
            'residue_id': residue_id,
            'residue_type': residue_type,
            'secondary_structure': secondary_structure,
            'asa': asa,
            'rsa': rsa,
            'phi': phi,
            'psi': psi
        })

    return residue_data

# Function to process multiple DSSP files
def process_all_dssp_files(dssp_dir, output_file):
    """Process DSSP files and save RSA data to a CSV file."""
    all_residue_data = []
    dssp_files = [f for f in os.listdir(dssp_dir) if f.endswith('.dssp')]

    for dssp_file in dssp_files:
        uniprot_id = os.path.splitext(dssp_file)[0]
        dssp_file_path = os.path.join(dssp_dir, dssp_file)
        residue_data = process_dssp_file(dssp_file_path, uniprot_id)
        all_residue_data.extend(residue_data)

    save_residue_data_to_csv(all_residue_data, output_file)

# Function to save residue data to a CSV file
def save_residue_data_to_csv(data, output_file):
    """Save residue data to a CSV file."""
    import csv
    keys = data[0].keys() if data else []
    with open(output_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        writer.writerows(data)

# Main function
if __name__ == "__main__":
    # Paths
    dssp_output_dir = snakemake.input["dssp_dir"]
    rsa_output_file = snakemake.output["rsa"]

    # Process DSSP files and save RSA data
    process_all_dssp_files(dssp_output_dir, rsa_output_file)
    print(f"RSA data saved to {rsa_output_file}")
