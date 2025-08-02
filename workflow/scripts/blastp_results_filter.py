import matplotlib.pyplot as plt
import pandas as pd
import os

# Inputs
afdb_blast = snakemake.input.afdb_blast
esm_blast = snakemake.input.esm_blast
gene_annotation = snakemake.input.gene_annotation

# Outputs
afdb_ids_out = snakemake.output.afdb_ids
esm_ids_out = snakemake.output.esm_ids
combined_out = snakemake.output.combined

# --- Step 1: Load filtered gene annotation ---
genes_df = pd.read_csv(gene_annotation, sep="\t")
afdb_genes_uniref100 = (
    genes_df
    .dropna(subset=["UniRef100", "AlphaFoldDB"])
    .loc[~genes_df["UniRef100"].str.startswith("UniRef100_UPI", na=False)]
    [["ID", "AlphaFoldDB"]]
    .rename(columns={"ID": "query_id", "AlphaFoldDB": "subject_id"})
)
afdb_genes_uniref100["database"] = "AFDB"
afdb_genes_uniref100["method"] = "AFDB-UNIREF100"
afdb_genes_uniref100["s_start"] = 0
afdb_genes_uniref100["q_start"] = 0
afdb_genes_uniref100["subject_id"] = afdb_genes_uniref100["subject_id"].str.replace(";", "", regex=False)

# --- Step 2: Load ESM blast results ---
esm_cols = ["query_id", "subject_id", "perc_identity", "alignment_length", 
            "mismatches", "gap_openings", "q_start", "q_end", "s_start", 
            "s_end", "evalue", "bit_score"]

esm_df = pd.read_csv(esm_blast, sep="\t", names=esm_cols)
esm_df = esm_df[
    (esm_df["q_start"] == 1) &
    (esm_df["q_end"] == esm_df["s_end"]) &
    (esm_df["s_start"] == 1) &
    (esm_df["gap_openings"] == 0) &
    (esm_df["perc_identity"] >= 98)
]
esm_df["subject_id"] = esm_df["subject_id"].str.replace(r"^tr\|([^|]+)\|.*$", r"\1", regex=True)
esm_df["database"] = "ESMATLAS"
esm_df["method"] = "ESMATLAS-BLASTP"

# --- Step 3: Load AFDB blast results ---
afdb_df = pd.read_csv(afdb_blast, sep="\t", names=esm_cols)
afdb_df["subject_id"] = afdb_df["subject_id"].str.replace(r".*AF-([A-Z0-9]+)-F1", r"\1", regex=True)
afdb_df = afdb_df[
    (afdb_df["q_start"] == 1) &
    (afdb_df["q_end"] == afdb_df["s_end"]) &
    (afdb_df["s_start"] == 1) &
    (afdb_df["gap_openings"] == 0) &
    (afdb_df["perc_identity"] >= 98)
]
afdb_df["database"] = "AFDB"
afdb_df["method"] = "AFDB-BLAST"

# --- Step 4: Combine ESM + AFDB blast ---
combined_df = pd.concat([afdb_df, esm_df], ignore_index=True)

# Full-length match flag
combined_df["priority_match"] = (
    (combined_df["q_start"] == 1) &
    (combined_df["q_end"] == combined_df["s_end"]) &
    (combined_df["s_start"] == 1)
)

# Tie-break preference: AFDB over ESM
combined_df["db_priority"] = combined_df["database"].map({"AFDB": 1, "ESMATLAS": 0})

# Sort with explicit AFDB preference in ties
combined_df = (
    combined_df.sort_values(
        ["query_id", "priority_match", "alignment_length", "perc_identity", "db_priority"],
        ascending=[True, False, False, False, False]
    )
    .groupby("query_id", as_index=False)
    .first()
)

# --- Step 5: Remove BLAST matches for genes with AFDB-UniRef100 ---
blast_without_afdb_uniref = combined_df[
    ~combined_df["query_id"].isin(afdb_genes_uniref100["query_id"])
]

# --- Step 6: Add AFDB-UniRef100 mappings (replacement step) ---
combined_final = pd.concat([blast_without_afdb_uniref, afdb_genes_uniref100], ignore_index=True)

# --- Step 7: Save outputs ---
os.makedirs(os.path.dirname(combined_out), exist_ok=True)
combined_final.to_csv(combined_out, sep="\t", index=False)

esm_ids = combined_final.loc[combined_final["database"] == "ESMATLAS", "subject_id"]
afdb_ids = combined_final.loc[combined_final["database"] == "AFDB", "subject_id"]

os.makedirs(os.path.dirname(esm_ids_out), exist_ok=True)
esm_ids.to_csv(esm_ids_out, index=False, header=False)
afdb_ids.to_csv(afdb_ids_out, index=False, header=False)

# --- Step 8: Plot distribution of methods ---
method_counts = combined_final["method"].value_counts()

plt.figure(figsize=(8, 5))
method_counts.plot(kind="bar")
plt.xlabel("Mapping Method / Database")
plt.ylabel("Number of Genes")
plt.title("Distribution of Structure Mapping Methods")
plt.xticks(rotation=45, ha="right")
plt.tight_layout()

# Save plot (species-specific)
plot_path = snakemake.output.method_plot
os.makedirs(os.path.dirname(plot_path), exist_ok=True)
plt.savefig(plot_path, dpi=300)
plt.close()
