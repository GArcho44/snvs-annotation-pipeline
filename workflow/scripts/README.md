** Overview **

This folder contains all auxiliary scripts used by the Snakemake workflow.
Each script performs a specific processing, annotation, filtering, or computation step within the pipeline.

## Table of Contents

- [species_info_filter.py](#species_info_filterpy)
- [blastp_results_filter.py](#blastp_results_filterpy)
- [download_structures_afdb_api.py](#download_structures_afdb_apipy)
- [download_structures_esm_api.py](#download_structures_esm_apipy)
- [dssp_generator.py](#dssp_generatorpy)
- [rsa_calculator.py](#rsa_calculatorpy)
- [ppi_merger.py](#ppi_mergerpy)
- [dtl_calculator.py](#dtl_calculatorpy)
- [dtppi_calculator.py](#dtppi_calculatorpy)
- [all_features_merger.py](#all_features_mergerpy)


## species_info_filter.py

**Description**

This script filters species-specific information required for downstream analysis including: annotated genes and SNPs, protein sequences fasta files and corresponding lists.

**Input**
-	Annotated snps (_.tsv_)
-	Annotated genes (_.tsv_)
-	Protein fasta database (_.fasta_)
-	Species code (_unique species identifier_)

**Output**
-	Snps IDs list (_snps_ids.txt_)
-	Gene IDs list (_genes_ids.txt_)
-	Filtered annotated SNPs (_snp_annotated.tsv_)
-	Filtered annotated genes (_genes_annotated.tsv_)
-	Filtered protein fasta (_protein_sequences.fasta_)

## blastp_results_filter.py

**Description**

This script filters species-specific BLASTP results of top-1 hits against AFDB and ESMAtlas protein sequence databases.

**Input**
-	AFDB blastp results (_blastp_afdb_top1.m8_)
-	ESMAtlas blastp results (_blastp_esm_top1.m8_)
-	Filtered annotated genes (_genes_annotated.tsv_)
  
**Output**
-	AFDB protein UniProtKB IDs list (_afdb_ids.txt_)
-	ESMatlas protein Mgnify IDs list (_esm_ids.txt_)
-	Combined results (_combined_mapping.tsv_)

## download_structures_afdb_api.py

**Description**

This script filters high-quality protein structures and uses APIs for species-specific download of protein structures from AFDB based on UniProtKB IDs.

**Input**
-	AFDB protein UniProtKB IDs list (_afdb_ids.txt_)
-	pLDDT threshold (_Default: 80_)

**Output**
-	HQ protein structures (_afdb/*.pdb_)
-	Brief report (_reports/species_code/afdb_report_)

## download_structures_esm_api.py

**Description**

This script filters high-quality protein structures and uses APIs for species-specific download of protein structures from ESMAtlas based on Mgnify IDs.

**Input**
-	ESMAtlas protein Mgnify IDs list (_afdb_ids.txt_)
-	pLDDT threshold (_Default: 80_)

**Output**
-	HQ protein structures (_esm/*.pdb_)
-	Brief report (_reports/species_code/esm_report_)

## dssp_generator.py

**Description**

This script generates dssp information for each HQ protein structure available.

**Input**
-	HQ protein structures (_structures/*.pdb_)

**Output**
-	DSSP per protein structure (_dssp_output/*.dssp_)


## rsa_calculator.py

**Description**

This script calculates Relative Solvent Accessibility (RSA) per residue per HQ protein structure.

**Input**
-	DSSP per protein structure (_dssp_output/*.dssp_)

**Output**
-	Species combined RSA calculations for all proteins (_rsa_output/species_code_dssp_rsa_combined.csv_)

## ppi_merger.py

**Description**

This script merges PeSTo predictions for all the proteins per species into a single.

**Input**
-	Protein structures with PeSTo predictions (_pesto_output/*.pdb_)

**Output**
-	Species combined PPI calculations for all proteins (_species_code_ppi_predictions.csv_)

## dtl_calculator.py

**Description**

This script calculates the Distance to Ligand (DTL) of each protein residue to the closest ligand-binding pocket detected with fpocketPRANK.

**Input**
-	Fpocket predictions output file (_*.pdb_predictions.csv_)
-	HQ protein structures (_structures/*.pdb_)

**Output**
-	Species combined DTL calculations for all proteins (_species_code_fpocket_dtl_combined.csv_)

## dtppi_calculator.py

**Description**

This script calculates the Distance to Protein-protein Interface (DTPPI) of each protein residue to the closest PPI predicted with PeSTo.

**Input**
-	Species combined PPI calculations for all proteins (_species_code_ppi_predictions.csv_)
-	HQ protein structures (_structures/*.pdb_)

**Output**
-	Species combined DTPPI calculations for all proteins (_species_code_pesto_dtppi_combined.csv_)

## all_features_merger.py

**Description**

This script merges all the assigned, and predicted properties in gene, SNP, affected residues and proteins levels.

**Input**
-	Filtered annotated SNPs (_snp_annotated.tsv_)
-	Filtered annotated genes (_genes_annotated.tsv_)
-	Combined BLASTP results (_combined_mapping.tsv_)
-	Species combined RSA calculations for all proteins (_species_code_dssp_rsa_combined.csv_)
-	Species combined DTL calculations for all proteins (_species_code_fpocket_dtl_combined.csv_)
-	Species combined DTPPI calculations for all proteins (_species_code_pesto_dtppi_combined.csv_)

**Output**
-	Species-specific predictions/calculations for all SNPs, genes, structures and residues available combined into a single file. This is the final output that includes all the functional annotations generated per species (_snp_structural_features_merged.tsv_)
