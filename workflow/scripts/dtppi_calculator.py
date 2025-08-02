import os
import glob
import pandas as pd
import numpy as np
from Bio.PDB import PDBParser, is_aa
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import fcluster, linkage


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
            res_ids.append((chain.id, residue.id[1]))
    return np.array(coords), res_ids


def group_interface_residues(pdb_file, pred_df, protein_id):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("model", pdb_file)
    model = structure[0]

    coords = []
    confidences = []
    res_ids = []

    pred_subset = pred_df[(pred_df["protein_id"] == protein_id) & (pred_df["confidence"] > 0.5)]
    pred_lookup = {
        (row["chain_id"], int(row["residue_id"])): row["confidence"]
        for _, row in pred_subset.iterrows()
    }

    for chain in model:
        for residue in chain:
            if not is_aa(residue) or "CA" not in residue:
                continue

            res_id = (chain.id, residue.id[1])
            if res_id not in pred_lookup:
                continue

            pLDDT = residue["CA"].get_bfactor()
            if pLDDT <= 70:
                continue

            coords.append(residue["CA"].get_coord())
            confidences.append(pred_lookup[res_id])
            res_ids.append(res_id)

    if len(coords) < 2:
        # Nothing to cluster if fewer than 2 residues
        return set(), []

    dist_matrix = pdist(coords)
    linkage_matrix = linkage(dist_matrix, method="single")
    cluster_labels = fcluster(linkage_matrix, t=10, criterion="distance")

    cluster_map = {}
    for label, res_id, conf in zip(cluster_labels, res_ids, confidences):
        cluster_map.setdefault(label, []).append((res_id, conf))

    high_quality_residues = set()
    
    cluster_summaries = []
    
    for cluster_id, cluster in cluster_map.items():
        avg_conf = np.mean([c[1] for c in cluster])
        if avg_conf > 0.8:
            res_ids = [c[0] for c in cluster]
            avg_plddt = np.mean([
            structure[0][res_id[0]][(" ", res_id[1], " ")]['CA'].get_bfactor()
            for res_id in res_ids
            ])
            summary = {
            "protein_id": protein_id,
            "cluster_id": cluster_id,
            "num_residues": len(res_ids),
            "avg_confidence": avg_conf,
            "avg_pLDDT": avg_plddt,
            "residues": ";".join([f"{chain}:{res}" for chain, res in res_ids])
            }
            cluster_summaries.append(summary)
            high_quality_residues.update([c[0] for c in cluster])

    return high_quality_residues, cluster_summaries


def calculate_interface_distance(protein_id, pdb_path, interface_residues):
    all_coords, all_residues = extract_structure_coords(pdb_path)

    interface_coords = [coord for res, coord in zip(all_residues, all_coords) if res in interface_residues]
    interface_coords = np.array(interface_coords)

    dist_list = []
    raw_distances = []
    for res_id, coord in zip(all_residues, all_coords):
        if res_id in interface_residues:
            dist = 0.0
        elif len(interface_coords) > 0:
            dist = np.min(np.linalg.norm(interface_coords - coord, axis=1))
            raw_distances.append((res_id, dist))
        else:
            continue

        max_dist = max(d for _, d in raw_distances)

        norm_dist = dist / max_dist if max_dist > 0 else 0.0
        dist_list.append((protein_id, res_id[0], res_id[1], dist, norm_dist))

    if dist_list:
        df = pd.DataFrame(dist_list, columns=["protein_id", "chain_id", "residue_id", "DTPPI", "normalized_DTPPI"])
        return df
    else:
        return None


def batch_interface_distance(structure_dir, pesto_csv_path, output_dir, combined_output_file):
    pdb_files = glob.glob(os.path.join(structure_dir, "*.pdb"))
    print(f"Found {len(pdb_files)} structures.")

    all_data = []
    all_cluster_summaries = []
    pred_df = pd.read_csv(pesto_csv_path)
    pred_df.columns = ["protein_id", "chain_id", "residue_id", "confidence"]
    pred_df["confidence"] = pred_df["confidence"].astype(float)

    for pdb_file in pdb_files:
        protein_id = os.path.splitext(os.path.basename(pdb_file))[0]
        output_file = os.path.join(output_dir, f"{protein_id}.csv")

        interface_residues, cluster_summaries = group_interface_residues(pdb_file, pred_df, protein_id)

        if not interface_residues:
            print(f"[!] Skipping {protein_id}: No high-quality interfaces found.")
            continue

        try:
            df = calculate_interface_distance(protein_id, pdb_file, interface_residues)
            if df is not None:
                os.makedirs(output_dir, exist_ok=True)
                df.to_csv(output_file, index=False)
                print(f"[✓] {protein_id} → {output_file}")
                all_data.append(df)
                all_cluster_summaries.extend(cluster_summaries)
            else:
                print(f"[!] {protein_id}: No distances calculated.")
        except Exception as e:
            print(f"[!] Error processing {protein_id}: {e}")

    if all_data:
        combined_df = pd.concat(all_data, ignore_index=True)
        os.makedirs(os.path.dirname(combined_output_file), exist_ok=True)
        combined_df.to_csv(combined_output_file, index=False)
        print(f"[✓] Combined interface distances saved to: {combined_output_file}")
        
    if all_cluster_summaries:
        cluster_summary_df = pd.DataFrame(all_cluster_summaries)
        cluster_summary_file = os.path.join(os.path.dirname(combined_output_file), "interface_clusters_summary.csv")
        cluster_summary_df.to_csv(cluster_summary_file, index=False)
        print(f"[✓] Cluster summaries saved to: {cluster_summary_file}")
    
    


# === Run it ===
if __name__ == "__main__":
    structure_dir = snakemake.input["structures_dir"]
    pesto_csv_path = snakemake.input["pesto_pred"]  # only includes: protein_id, chain_id, residue_id, confidence
    output_dir = snakemake.output["sum_dir"]
    combined_output_file = snakemake.output["dtppi_combined"]

    batch_interface_distance(structure_dir, pesto_csv_path, output_dir, combined_output_file)

