import os
import requests
import logging
import json
import matplotlib.pyplot as plt
import re
from concurrent.futures import ThreadPoolExecutor, as_completed
from collections import defaultdict, Counter
from threading import Lock
from Bio.PDB.MMCIFParser import MMCIFParser


# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("alphafold_script.log"),
        logging.StreamHandler()
    ]
)

# Constants for AlphaFold API
ALPHAFOLD_SUM_API_URL = "https://www.alphafold.ebi.ac.uk/api/uniprot/summary/"
ALPHAFOLD_API_URL = "https://www.alphafold.ebi.ac.uk/api/prediction/"

# Thread-safe data structures
avg_confidence_scores = []
filtered_uniprot_ids = []
organisms = []  # To store organism names for filtered structures
lock = Lock()  # Lock to synchronize access to shared resources


# Function to fetch metadata and filter based on pLDDT score
def fetch_metadata_and_filter(uniprot_id, threshold):
    """Fetch metadata from AlphaFold API and filter based on pLDDT score."""
    url = f"{ALPHAFOLD_SUM_API_URL}{uniprot_id}.json"
    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status()  # Raise exception for HTTP errors (non-2xx responses)
        metadata = response.json()

        # Log the metadata for debugging if needed
        logging.debug(f"Full metadata for {uniprot_id}: {json.dumps(metadata, indent=2)}")

        # Check for structures and filter based on pLDDT
        if 'structures' in metadata:
            structures = metadata['structures']
            for structure in structures:
                summary = structure.get('summary', {})
                avg_confidence_score = summary.get('confidence_avg_local_score')

                if avg_confidence_score is not None:
                    logging.info(f"avg_pLDDT for {uniprot_id}: {avg_confidence_score}")

                    # Filter by pLDDT threshold
                    if avg_confidence_score >= threshold:
                        logging.info(f"Structure for {uniprot_id} passes pLDDT filter: {avg_confidence_score:.2f}")
                        return summary  # Return the structure's metadata
                    else:
                        logging.info(f"Structure for {uniprot_id} does not pass pLDDT filter: {avg_confidence_score}")
        else:
            logging.warning(f"No structures found for {uniprot_id}.")
    except requests.RequestException as e:
        logging.error(f"Request failed for {uniprot_id}: {e}")
    return None


# Function to download the structure if it passes the pLDDT filter
def download_structure_pdb(uniprot_id, output_dir):
    """Download the structure file using the pdburl from the AlphaFold API response."""
    url = f"{ALPHAFOLD_API_URL}{uniprot_id}"  # API endpoint for the UniProt ID
    try:
        # Fetch metadata from AlphaFold API
        response = requests.get(url, timeout=10)
        response.raise_for_status()  # Raise exception for HTTP errors (non-2xx responses)

        data = response.json()
        pdburl = data[0].get("pdbUrl")  # Assuming "pdb_url" is part of the API response
        organism = data[0].get("organismScientificName")

        if not pdburl:
            logging.warning(f"No pdburl found for {uniprot_id}.")
            return None

        # Download PDB file
        pdb_response = requests.get(pdburl, timeout=10)
        pdb_response.raise_for_status()  # Raise exception for HTTP errors
        pdb_file_path = os.path.join(output_dir, f"{uniprot_id}.pdb")

        with open(pdb_file_path, 'wb') as pdb_file:
            pdb_file.write(pdb_response.content)
        logging.info(f"Downloaded structure for {uniprot_id} to {pdb_file_path}")
        return pdb_file_path, organism

    except requests.RequestException as e:
        logging.error(f"Error downloading structure for {uniprot_id}: {e}")
    return None


def download_structure_mmcif(uniprot_id, output_dir):
    """Download the structure file using the mmCIF URL from the AlphaFold API response."""
    url = f"{ALPHAFOLD_API_URL}{uniprot_id}"  # API endpoint for the UniProt ID
    try:
        # Fetch metadata from AlphaFold API
        response = requests.get(url, timeout=10)
        response.raise_for_status()  # Raise exception for HTTP errors

        data = response.json()
        cif_url = data[0].get("cifUrl")  # Assuming "cifUrl" is part of the API response
        organism = data[0].get("organismScientificName")

        if not cif_url:
            logging.warning(f"No mmCIF URL found for {uniprot_id}.")
            return None

        # Download mmCIF file
        cif_response = requests.get(cif_url, timeout=10)
        cif_response.raise_for_status()  # Raise exception for HTTP errors
        cif_file_path = os.path.join(output_dir, f"{uniprot_id}.mmcif")

        with open(cif_file_path, 'wb') as cif_file:
            cif_file.write(cif_response.content)
        logging.info(f"Downloaded structure for {uniprot_id} to {cif_file_path}")
        return cif_file_path, organism

    except requests.RequestException as e:
        logging.error(f"Error downloading structure for {uniprot_id}: {e}")
    return None


# Function to parse confidence scores from PDB file
def parse_confidence_score_pdb(file_path):
    """Extract per-residue pLDDT scores from PDB file."""
    res_confidence_scores = []
    residue_scores = defaultdict(list)  # To store atom confidence scores for each residue
    try:
        with open(file_path, 'r') as f:
            for line in f:
                if line.startswith("ATOM"):
                    residue_index = int(line[22:26].strip())  # Residue index (adjust based on PDB format)
                    chain_id = line[21]  # Chain ID
                    b_factor = line[60:66].strip()  # B-factor column, which represents confidence score

                    try:
                        atom_confidence_score = float(b_factor)
                    except ValueError:
                        continue  # Skip lines where B-factor isn't a valid number

                    # Group by residue (residue_index, chain_id)
                    residue_scores[(residue_index, chain_id)].append(atom_confidence_score)

        # Extract pLDDT scores for each residue (using first atom in residue)
        for residue, scores in residue_scores.items():
            res_confidence_scores.append(scores[0])  # Take the confidence score of the first atom

        # Compute average pLDDT score
        avg_confidence_score = sum(res_confidence_scores) / len(res_confidence_scores) if res_confidence_scores else 0
        return res_confidence_scores, avg_confidence_score

    except Exception as e:
        logging.error(f"Error parsing file {file_path}: {e}")
        return [], 0


def parse_confidence_score_mmcif(file_path):
    """Extract per-residue pLDDT scores from mmCIF file."""
    res_confidence_scores = []
    try:
        # Use Biopython's MMCIFParser to parse the file
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure("structure", file_path)

        # Loop through atoms and extract pLDDT scores from the B-factor
        for model in structure:
            for chain in model:
                for residue in chain:
                    if not residue.id[0].strip():  # Ignore heteroatoms (e.g., water, ligands)
                        atom_scores = [
                            atom.bfactor for atom in residue if atom.bfactor is not None
                        ]
                        if atom_scores:
                            # Use the first atom's B-factor as the residue's pLDDT score
                            res_confidence_scores.append(atom_scores[0])

        # Compute average pLDDT score
        avg_confidence_score = (
            sum(res_confidence_scores) / len(res_confidence_scores)
            if res_confidence_scores
            else 0
        )
        return res_confidence_scores, avg_confidence_score

    except Exception as e:
        logging.error(f"Error parsing file {file_path}: {e}")
        return [], 0


# Function to create a plot for pLDDT distribution
def plot_confidence_score_distribution(avg_con_scores, output_dir):
    """Plot the distribution of average pLDDT scores."""
    plt.figure(figsize=(10, 6))
    plt.hist(avg_con_scores, bins=20, color='skyblue', edgecolor='black')
    plt.title("Distribution of Average pLDDT Scores", fontsize=16)
    plt.xlabel("Average pLDDT Score", fontsize=14)
    plt.ylabel("Frequency", fontsize=14)
    plt.grid(True)
    plot_file = os.path.join(output_dir, "average_pLDDT_distribution.png")
    plt.savefig(plot_file)
    plt.close()
    logging.info(f"pLDDT distribution plot saved to: {plot_file}")


def shorten_organism_name(organism_name):
    """Shorten organism name to 'Genus species' format."""
    # Regular expression to match the genus and species, and optionally remove strain info
    match = re.match(r"([A-Za-z]+)\s([a-z]+)(?:\s\(strain.*\))?", organism_name)
    if match:
        genus = match.group(1)[0]  # First letter of the genus
        species = match.group(2)  # Full species name
        return f"{genus}. {species}"  # Shortened format like "B. subtilis"
    return organism_name  # If the format doesn't match, return the original name


# Function to plot organism distribution with shortened names
def plot_organism_distribution(organism_list, output_dir):
    """Plot the distribution of organisms for filtered structures with shortened names."""
    # Shorten all organism names
    shortened_organisms = [shorten_organism_name(organism) for organism in organism_list]

    # Count the occurrences of each organism
    organism_counts = Counter(shortened_organisms)

    # Create a bar plot
    plt.figure(figsize=(12, 8))  # Increase the figure size to allow more space
    plt.bar(organism_counts.keys(), organism_counts.values(), color='lightgreen', edgecolor='black')
    plt.title("Distribution of Organisms for Filtered Structures", fontsize=16)
    plt.xlabel("Organism", fontsize=14)
    plt.ylabel("Number of Structures", fontsize=14)

    # Rotate x-axis labels for better readability
    plt.xticks(rotation=90, fontsize=12)

    # Adjust layout to prevent label clipping and add extra padding
    plt.tight_layout(pad=2.0)  # Increase padding for more space around the plot

    # Save plot
    plot_file = os.path.join(output_dir, "organism_distribution.png")
    plt.savefig(plot_file)
    plt.close()
    logging.info(f"Organism distribution plot saved to: {plot_file}")



# Main workflow
def main():
    # Input and output configuration
    input_file = snakemake.input.ids  # Path to UniProt IDs
    output_dir = snakemake.output.structures_dir  # Directory to save data
    report_dir = snakemake.output.report_dir
    threshold = snakemake.params.threshold  # pLDDT threshold
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(report_dir, exist_ok=True)

    # Read UniProt IDs
    with open(input_file, 'r') as f:
        uniprot_ids = [line.strip() for line in f]

    logging.info(f"Loaded {len(uniprot_ids)} UniProt IDs from {input_file}.")

    # Output files
    scores_file = os.path.join(report_dir, "pLDDT_scores.tsv")
    filtered_file = os.path.join(report_dir, "afdb_filtered_structures.txt")

    # Open files for writing
    with open(scores_file, "w") as sf, open(filtered_file, "w") as ff:
        sf.write("UniProt_ID\tResidue_Index\tpLDDT\tAverage_pLDDT\n")  # Header

        # Concurrent processing
        def process_uniprot_id(uniprot_id):
            metadata = fetch_metadata_and_filter(uniprot_id, threshold)
            if metadata:
                logging.info(f"Processing structure for {uniprot_id}.")
                pdb_file, organism = download_structure_pdb(uniprot_id, output_dir)
                if pdb_file:
                    # Parse confidence scores
                    res_confidence_scores, avg_confidence_score = parse_confidence_score_pdb(pdb_file)

                    with lock:  # Ensure thread safety when modifying shared resources
                        avg_confidence_scores.append(avg_confidence_score)

                    # Write per-residue scores to CSV
                    with lock:
                        for idx, score in enumerate(res_confidence_scores, start=1):
                            sf.write(f"{uniprot_id}\t{idx}\t{score:.2f}\t{avg_confidence_score:.2f}\n")

                    # Filter and save UniProt ID if avg pLDDT passes threshold
                    if avg_confidence_score >= threshold:
                        with lock:
                            filtered_uniprot_ids.append(uniprot_id)
                            ff.write(f"{uniprot_id}\n")
                            if organism:
                                organisms.append(organism)  # Add organism to the list

        # Use ThreadPoolExecutor for concurrent processing
        with ThreadPoolExecutor(max_workers=10) as executor:
            futures = {executor.submit(process_uniprot_id, uniprot_id): uniprot_id for uniprot_id in uniprot_ids}

            # Collect results
            for future in as_completed(futures):
                try:
                    future.result()  # Trigger exception if the task failed
                except Exception as e:
                    logging.error(f"Error processing UniProt ID {futures[future]}: {e}")

    # Generate distribution plot
    if avg_confidence_scores:
        plot_confidence_score_distribution(avg_confidence_scores, report_dir)

    # Plot the organism distribution
    if organisms:
        plot_organism_distribution(organisms, report_dir)

    # Summary
    logging.info(f"Total UniProt IDs: {len(uniprot_ids)}")
    logging.info(f"Filtered structures with average pLDDT >= {threshold}: {len(filtered_uniprot_ids)}")
    logging.info(f"Results saved to {output_dir}.")


if __name__ == "__main__":
    main()

