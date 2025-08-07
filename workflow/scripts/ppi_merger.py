import os
import csv
import logging
from collections import defaultdict

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("logs/pesto_conversion.log"),
        logging.StreamHandler()
    ]
)


def parse_pesto_score_pdb(file_path):
    """Extract per-residue pesto scores from PDB file."""
    residue_scores = defaultdict(float)

    try:
        with open(file_path, 'r') as f:
            for line in f:
                if line.startswith("ATOM"):
                    residue_index = int(line[22:26].strip())  # Residue index
                    chain_id = line[21]  # Chain ID
                    b_factor = line[60:66].strip()  # B-factor column

                    try:
                        pesto_score = float(b_factor)
                    except ValueError:
                        continue  # Skip invalid values


                    # Store the first pesto score per residue (since they are the same for all atoms)
                    if (chain_id, residue_index) not in residue_scores:
                        residue_scores[(chain_id, residue_index)] = pesto_score

            return residue_scores
    except Exception as e:
        logging.error(f"Error parsing file {file_path}: {e}")
        return {}


def process_pdb_files(input_folder, output_csv):
    """Process all PDB files in the input folder and export data to CSV."""
    data = []

    for filename in os.listdir(input_folder):
        if filename.endswith("i0.pdb"):
            file_path = os.path.join(input_folder, filename)
            uniprot_id = filename.split("_")[0]  # Extract UniProtKB ID

            residue_data = parse_pesto_score_pdb(file_path)
            for (chain, residue), score in residue_data.items():
                data.append([uniprot_id, chain, residue, score])

    # Save results to CSV
    with open(output_csv, mode='w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["UniProtKB_ID", "Chain", "Residue", "Confidence_Score"])
        writer.writerows(data)

    logging.info(f"Results saved to {output_csv}")


def main():
    input_folder = snakemake.input["pesto_dir"]  # Input folder containing PDB files
    output_csv = snakemake.output["merged_data"]  # Output CSV file
    os.makedirs(os.path.dirname(output_csv), exist_ok=True)

    process_pdb_files(input_folder, output_csv)


if __name__ == "__main__":
    main()
