import os
import glob
import pandas as pd
import numpy as np
from Bio.PDB import PDBParser, is_aa

def parse_fpocket_centers(file_path):
    headers = ["name", "rank", "score", "probability", "sas_points", "surf_atoms",
               "center_x", "center_y", "center_z", "residue_ids", "surf_atom_ids"]
    df = pd.read_csv(file_path, names=headers, skiprows=1)

    # keep confident pockets only
    df = df[df["probability"].astype(float) > 0.5].reset_index(drop=True)

    centers = df[["center_x", "center_y", "center_z"]].astype(float).values
    return centers, df

def get_sidechain_center(residue):
    sidechain_atoms = [atom for atom in residue if atom.get_name() not in ["N", "CA", "C", "O"]]
    if not sidechain_atoms:
        return residue['CA'].get_coord()
    coords = np.array([atom.get_coord() for atom in sidechain_atoms])
    return coords.mean(axis=0)

def extract_structure_coords(pdb_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("model", pdb_file)
    model = structure[0]
    coords = []
    res_ids = []
    for chain in model:
        for residue in chain:
            if not is_aa(residue):
                continue
            coords.append(get_sidechain_center(residue))
            res_ids.append((chain.id, residue.id[1]))  # (chain_id, residue_number)
    return np.array(coords), res_ids

def calculate_dtl(protein_id, pdb_path, fpocket_path, output_path):
    all_coords, all_residues = extract_structure_coords(pdb_path)

    # Use only Fpocket pocket centers with probability > 0.5
    binding_coords, pocket_df = parse_fpocket_centers(fpocket_path)
    if len(binding_coords) == 0:
        binding_coords = np.empty((0, 3))
    else:
        binding_coords = np.vstack(binding_coords)

    dtl_list = []
    raw_distances = []
    for res_id, coord in zip(all_residues, all_coords):
        if len(binding_coords) > 0:
            distances = np.linalg.norm(binding_coords - coord, axis=1)
            dtl = np.min(distances)
            idx = int(np.argmin(distances))  # nearest pocket index
            raw_distances.append((res_id, dtl))
        else:
            continue

        max_dist = max(d for _, d in raw_distances)

        norm_dtl = dtl / max_dist if max_dist > 0 else 0.0

        pocket_row = pocket_df.iloc[idx]  # pocket metadata

        dtl_list.append((
            protein_id,
            res_id[0],
            res_id[1],
            dtl,
            norm_dtl,
            *pocket_row.tolist()  # keep pocket metadata
        ))

    if dtl_list:
        base_cols = ["protein_id", "chain_id", "residue_id", "DTL", "normalized_DTL"]
        df = pd.DataFrame(dtl_list, columns=base_cols + pocket_df.columns.tolist())
        return df
    else:
        return None

def batch_dtl(structure_dir, fpocket_dir, output_dir, combined_output_file):
    pdb_files = glob.glob(os.path.join(structure_dir, "*.pdb"))
    print(f"Found {len(pdb_files)} structures.")

    all_dtl_data = []

    for pdb_file in pdb_files:
        protein_id = os.path.splitext(os.path.basename(pdb_file))[0]
        fpocket_file = os.path.join(fpocket_dir, f"{protein_id}.pdb_predictions.csv")
        per_protein_output_file = os.path.join(output_dir, f"{protein_id}.csv")

        if not os.path.exists(fpocket_file):
            print(f"[!] Skipping {protein_id}: Fpocket file missing.")
            continue

        try:
            df = calculate_dtl(protein_id, pdb_file, fpocket_file, per_protein_output_file)
            if df is not None:
                os.makedirs(output_dir, exist_ok=True)
                df.to_csv(per_protein_output_file, index=False)
                print(f"[✓] {protein_id} → {per_protein_output_file}")
                all_dtl_data.append(df)
            else:
                print(f"[!] {protein_id}: No DTL values calculated.")
        except Exception as e:
            print(f"[!] Error processing {protein_id}: {e}")

    if all_dtl_data:
        combined_df = pd.concat(all_dtl_data, ignore_index=True)
        os.makedirs(os.path.dirname(combined_output_file), exist_ok=True)
        combined_df.to_csv(combined_output_file, index=False)
        print(f"[✓] Combined DTL saved to: {combined_output_file}")

# === Run it ===
if __name__ == "__main__":
    structure_dir = snakemake.input["structures_dir"]
    fpocket_dir = snakemake.input["fpocket_dir"]
    output_dir = snakemake.output["dtl_dir"]
    combined_output_file = snakemake.output["dtl_combined"]

    batch_dtl(structure_dir, fpocket_dir, output_dir, combined_output_file)

