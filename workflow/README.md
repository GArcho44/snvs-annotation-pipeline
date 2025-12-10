# Snakefile Overview

The Snakefile (snakemake workflow) includes 19 rules in total (see list below).

## List of rules

- [Rule 01: Retrieve species-level information including SNPs, Genes annotations and protein sequences](#smk_rule01)
- [Rule 02: Retrieve UniProtKB IDs based on BLASTP top-1 hits of query sequences against AFDB](#smk_rule02)
- [Rule 03: Retrieve MGnify IDs based on BLASTP top-1 hits of query sequences against ESM atlas database](#smk_rule03)
- [Rule 04: Process BLASTP results from AFDB and ESM atlas databases](#smk_rule04)
- [Rule 05: Download HQ structures from AFDB with UniProtKB IDs using APIs](#smk_rule05)
- [Rule 06: Download HQ structures from ESM atlas with MGnify IDs using APIs](#smk_rule06)
- [Rule 07: Move AFDB & ESM structure files into a single directory for downstream steps](#smk_rule07)
- [Rule 08: Run DSSP per structure for secondary structure features assignment](#smk_rule08)
- [Rule 09: Calculate RSA values per residue per structure and save them into a single file](#smk_rule09)
- [Rule 10: Create list of relative paths of structures inside P2Rank tool directory](#smk_rule10)
- [Rule 11: Run P2Rank to predict pockets](#smk_rule11)
- [Rule 12: Fpocket (PRANK re-scored) to detect pockets](#smk_rule12)
- [Rule 13: Merge the output from the two pocket prediction/detection tools](#smk_rule13)
- [Rule 14: Predict protein-protein interaction interfaces with PeSTo](#smk_rule14)
- [Rule 15: Merge PPI scores for all proteins](#smk_rule15)
- [Rule 16: Calculate Distance to Ligand (DTL) based on Fpocket (PRANK re-scored) output](#smk_rule16)
- [Rule 17: Calculate Distance to Protein-protein Interface (DTPPI) based on PeSTo output](#smk_rule17)
- [Rule 18: Merge all the information](#smk_rule18)
- [Rule 19: Visualize information](#smk_rule19)

## Rule 01: Retrieve species-level information including SNPs, Genes annotations and protein sequences {#smk_rule01}

**Input**
- **snps_tsv**: the SNP annotation file generated with SnpEff
- **genes_tsv**: the Genes annotation file generated with bakta
- **protein_fasta**: the protein sequence fasta file generated with transeq for all the gene sequences

**Output**
- **snps_list**: list of SNPs unique IDs
- **genes_list**: list of Genes unique IDs
- **snps_annotated**: species filtered annotated snps
- **genes_annotated**: species filtered annotated genes
- **protein_fasta_filtered**: species filtered protein sequences

## Rule 02: Retrieve UniProtKB IDs based on BLASTP top-1 hits of query sequences against AFDB {#smk_rule02}

**Input**
- **db**: the AFDB Diamond database to be used as reference
- **query**: fasta file of the query protein sequences to be used for the BLASTP search

**Output**
- **blastp_results**: Diamond BLASTP results of top-1 hit per protein sequence for AFDB

## Rule 03: Retrieve MGnify IDs based on BLASTP top-1 hits of query sequences against ESM atlas database {#smk_rule03}

**Input**
- **db**: the AFDB Diamond database to be used as reference
- **query**: fasta file of the query protein sequences to be used for the BLASTP search

**Output**
- **blastp_results**: Diamond BLASTP results of top-1 hit per protein sequence for AFDB

## Rule 04: Process BLASTP results from AFDB and ESM atlas databases {#smk_rule04}

**Input**
- **afdb_blast**: Diamond BLASTP results of top-1 hit per protein sequence for AFDB
- **esm_blast**: Diamond BLASTP results of top-1 hit per protein sequence for ESMatlas database
- **gene_annotation**: species filtered annotated genes

**Output**
- **afdb_ids**: UniProtKB IDs of the protein sequences
- **esm_ids**: Mgnify IDs of the protein sequences
- **combined_results**: Combined BLASTP results with best hit per query sequence based on pre-defined criteria (sequence identity, full length alignment and AFDB priority)
- **method_plot**: summary plot of the distribution per method and database

## Rule 05: Retrieve Download HQ structures from AFDB with UniProtKB IDs using APIs {#smk_rule05}
**Input**
- **ids**: UniProtKB IDs of the protein sequences

**Output**
- **structures_dir**: directory that contains all the downloaded HQ protein models from AFDB
- **report_dir**: directory that contains pLDDT confidence scores information per downloaded protein & per protein residue

## Rule 06: Download HQ structures from ESM atlas with MGnify IDs using APIs {#smk_rule06}

**Input**
- **ids**: UniProtKB IDs of the protein sequences

**Output**
- **structures_dir**: directory that contains all the downloaded HQ protein models from AFDB
- **report_dir**: directory that contains pLDDT confidence scores information per downloaded protein & per protein residue

## Rule 07: Move AFDB & ESM structure files into a single directory for downstream steps {#smk_rule07}

**Input**
- **afdb_dir**: directory that contains all the downloaded HQ protein models from AFDB
- **esm_dir**: directory that contains all the downloaded HQ protein models from ESMatlas

**Output**
- **merged_dir**: directory that contains all the downloaded HQ protein models (AFDB + ESMatlas)

## Rule 08: Run DSSP per structure for secondary structure features assignment {#smk_rule08}

**Input**
- **structures_dir**: directory that contains all the downloaded HQ protein models (AFDB + ESMatlas)

**Output**
- **dsspi_dir**: directory with the dssp secondary structure features assignments per protein structure

## Rule 09: Calculate RSA values per residue per structure and save them into a single file {#smk_rule09}

**Input**
- **dsspi_dir**: directory with the dssp secondary structure features assignments per protein structure

**Output**
- **rsa**: single file containing secondary elements assignment and calculated RSA per residue per structure

## Rule 10: Create list of relative paths of structures inside P2Rank tool directory {#smk_rule10}

**Input**
- **structures_dir**: directory that contains all the downloaded HQ protein models (AFDB + ESMatlas)

**Output**
- **structures_list**: A list of relative paths for the structures inside the P2Rank directory

## Rule 11: Run P2Rank to predict pockets {#smk_rule11}

**Input**
- **structures_list**: A list of relative paths for the structures inside the P2Rank directory

**Output**
- **predictions_dir**: directory that contains P2Rank predictions output

## Rule 12: Fpocket (PRANK re-scored) to detect pockets {#smk_rule12}

**Input**
- **structures_list**: A list of relative paths for the structures inside the P2Rank directory

**Output**
- **predictions_dir**: directory that contains Fpocket (PRANK re-scored) predictions output

## Rule 13: Merge the output from the two pocket prediction/detection tools {#smk_rule13}

**Input**
- **p2rank_dir**: directory that contains P2Rank predictions output 
- **fpocket_dir**: directory that contains Fpocket (PRANK re-scored) predictions output

**Output**
- **merged_data**: file that contains P2Rank and Fpocket (PRANK re-scored) predictions merged

## Rule 14: Predict protein-protein interaction interfaces with PeSTo {#smk_rule14}

**Input**
- **structures_dir**: directory that contains all the downloaded HQ protein models (AFDB + ESMatlas)

**Output**
- **predictions_dir**: directory that contains PeSTo predictions output

## Rule 15: Merge PPI scores for all proteins {#smk_rule15}
**Input**
- **pesto_dir**: directory that contains PeSTo predictions output

**Output**
- **merged_data**: PPI predictions for all proteins in a single file to be used for DTPPI calculations

## Rule 16: Calculate Distance to Ligand (DTL) based on Fpocket (PRANK re-scored) output {#smk_rule16}

**Input**
- **structures_dir**: directory that contains all the downloaded HQ protein models (AFDB + ESMatlas)
- **fpocket_dir**: directory that contains Fpocket (PRANK re-scored) predictions output

**Output**
- **dtl_dir**: directory of per protein DTL calculations
- **dtl_combined**: file of combined DTL calculations per residue per protein structure for all proteins with detected pockets

## Rule 17: Calculate Distance to Protein-protein Interface (DTPPI) based on PeSTo output {#smk_rule17}

**Input**
- **structures_dir**: directory that contains all the downloaded HQ protein models (AFDB + ESMatlas)
- **pesto_pred**: PeSTo PPI predictions for all proteins in a single file 

**Output**
- **dtppi_dir**: directory of per protein DTPPI calculations
- **dtppi_combined**: file of combined DTPPI calculations per residue per protein structure for all proteins with detected interfaces

## Rule 18: Merge all the information {#smk_rule18}

**Input**
- **snps_annotated**: species filtered annotated snps
- **genes_annotated**: species filtered annotated genes
- **mapping**: Combined BLASTP results with best hit per query sequence based on pre-defined criteria (sequence identity, full length alignment and AFDB priority)
- **rsa**: single file containing secondary elements assignment and calculated RSA per residue per structure
- **dtl**: file of combined DTL calculations per residue per protein structure for all proteins with detected pockets
- **dtppi**: file of combined DTPPI calculations per residue per protein structure for all proteins with detected interfaces

**Output**
- **dtppi_dir**: directory of per protein DTPPI calculations
- **combined**: Single file that contains all generated information per SNP

## Rule 19: Visualize information {#smk_rule19}
**Input**
- **merged**: Single file that contains all generated information per SNP
- **gene_counts**: file containing information about the number of unique genes per species

**Output**
- **combined**: SNP/protein structure coverage summary plot
