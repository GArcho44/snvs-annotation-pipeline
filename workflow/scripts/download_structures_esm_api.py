import requests
import os

# === Configuration ===
input_file = snakemake.input.ids  # Your text file with MGYP IDs
output_folder = snakemake.output.structures_dir
plddt_output_file = snakemake.output.plddt_scores  # New file for pLDDT scores
plddt_threshold = snakemake.params.plddt_threshold  # Keep this as a float (0–1 scale)

os.makedirs(output_folder, exist_ok=True)

# Load IDs
with open(input_file, "r") as f:
    mgyp_ids = [line.strip() for line in f if line.strip()]

# Create / overwrite pLDDT scores file
with open(plddt_output_file, "w") as pf:
    pf.write("MGYP_ID\tResidue_Index\tpLDDT\tAverage_pLDDT\n")  # Write header

# Function to compute average pLDDT
def get_avg_plddt(plddt_list):
    return sum(plddt_list) / len(plddt_list) if plddt_list else 0

# Initialize counters
total = len(mgyp_ids)
available = 0
downloaded = 0
skipped = 0
failed = 0

# Add basic headers to prevent potential 403s
headers = {"User-Agent": "Mozilla/5.0"}

# Main loop
for mgyp_id in mgyp_ids:
    try:
        conf_url = f"https://api.esmatlas.com/fetchConfidencePrediction/{mgyp_id}"
        conf_response = requests.get(conf_url, headers=headers)
        if conf_response.status_code in (403, 404):
            failed += 1
            print(f"Not available: {mgyp_id}")
            continue

        conf_response.raise_for_status()
        conf_data = conf_response.json()

        avg_plddt = get_avg_plddt(conf_data.get("plddt", []))
        available += 1

        if avg_plddt > plddt_threshold:  # Keep if pLDDT > 0.80
            struct_url = f"https://api.esmatlas.com/fetchPredictedStructure/{mgyp_id}"
            struct_response = requests.get(struct_url, headers=headers)
            struct_response.raise_for_status()

            # Save PDB file
            pdb_path = os.path.join(output_folder, f"{mgyp_id}.pdb")
            with open(pdb_path, "w") as f:
                f.write(struct_response.text)

            # Save pLDDT per residue
            with open(plddt_output_file, "a") as pf:
                for idx, plddt_value in enumerate(conf_data.get("plddt", []), start=1):
                    pf.write(f"{mgyp_id}\t{idx}\t{plddt_value:.2f}\t{avg_plddt:.2f}\n")

            print(f"Downloaded: {mgyp_id} (Avg pLDDT: {avg_plddt:.2f})")
            downloaded += 1
        else:
            print(f"Skipped: {mgyp_id} (Low pLDDT: {avg_plddt:.2f})")
            skipped += 1

    except Exception as e:
        failed += 1
        print(f"Error for {mgyp_id}: {e}")

# === Summary ===
print("\n=== Summary ===")
print(f"Total IDs processed     : {total}")
print(f"Structures available    : {available} ({(available/total)*100:.2f}%)")
print(f"High-quality downloaded : {downloaded} ({(downloaded/total)*100:.2f}%)")
print(f"Skipped (low pLDDT)      : {skipped} ({(skipped/total)*100:.2f}%)")
print(f"Failed / Not available  : {failed} ({(failed/total)*100:.2f}%)")

