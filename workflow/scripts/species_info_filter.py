import os
import csv
from Bio import SeqIO

# Inputs from Snakemake
input_snv_file = snakemake.input.snvs_tsv
input_gene_file = snakemake.input.genes_tsv
input_fasta_file = snakemake.input.protein_fasta
species_code = str(snakemake.config["species"])

output_snvs_list = snakemake.output.snvs_list
output_genes_list = snakemake.output.genes_list
output_snv_annot = snakemake.output.snv_annotated
output_gene_annot = snakemake.output.gene_annotated
output_fasta = snakemake.output.protein_fasta_filtered

# Variant types to keep
valid_variants = {"synonymous_variant", "missense_variant"}

snv_set = set()
gene_set = set()

# --- STEP 1: Filter SNV annotation ---
os.makedirs(os.path.dirname(output_snv_annot), exist_ok=True)

with open(input_snv_file, "r", encoding="utf-8") as infile, open(output_snv_annot, "w", encoding="utf-8") as outfile:
    header = infile.readline()
    outfile.write(header)

    for line in infile:
        parts = line.rstrip("\n").split("\t")
        if len(parts) > 10:
            if parts[0] == species_code and parts[10] in valid_variants:
                snv_set.add(parts[1])
                gene_set.add(parts[7])
                outfile.write(line)

# Save unique SNV & Gene IDs
os.makedirs(os.path.dirname(output_snvs_list), exist_ok=True)
with open(output_snvs_list, "w") as f:
    f.write("\n".join(sorted(snv_set)))
with open(output_genes_list, "w") as f:
    f.write("\n".join(sorted(gene_set)))


# --- STEP 2: Filter Gene annotation ---
found_genes = set()

os.makedirs(os.path.dirname(output_gene_annot), exist_ok=True)
with open(input_gene_file, newline='', encoding='utf-8') as infile, open(output_gene_annot, 'w', newline='', encoding='utf-8') as outfile:
    reader = csv.DictReader(infile, delimiter='\t')
    writer = csv.DictWriter(outfile, fieldnames=reader.fieldnames, delimiter='\t')
    writer.writeheader()

    for row in reader:
        if row['ID'] in gene_set and row['UniRef100'].strip():
            writer.writerow(row)
            found_genes.add(row['ID'])

# --- STEP 3: Report missing genes ---
missing_genes = gene_set - found_genes
if missing_genes:
    print(f"[INFO] {len(missing_genes)} gene IDs from SNVs not found in gene annotation file.")
else:
    print("[INFO] All gene IDs from SNVs were found in the gene annotation file.")


# --- STEP 4: Filter FASTA sequences ---
os.makedirs(os.path.dirname(output_fasta), exist_ok=True)
filtered_sequences = []

for record in SeqIO.parse(input_fasta_file, "fasta"):
    # Extract the gene ID from header (first two parts before underscore)
    gene_id = '_'.join(record.id.split('_')[:2])
    
    if gene_id in gene_set:
        record.id = gene_id  # Set header ID
        record.description = ''  # Description can be left blank or customized
        filtered_sequences.append(record)

SeqIO.write(filtered_sequences, output_fasta, "fasta")

print(f"[INFO] Filtered FASTA contains {len(filtered_sequences)} sequences.")

