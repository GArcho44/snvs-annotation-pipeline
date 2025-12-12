## Workflow overview

This workflow is a best-practice workflow for `<structure-aware annotation of microbial SNPs>`.
The workflow is built using [snakemake](https://snakemake.readthedocs.io/en/stable/) and consists of the following 3 main configurable modules:

- A. Sequence-based functional annotation of SNPs & Genes (later release)
- B. Retrieval of high-confidence structure models from AFDB and ESMatlas protein structure databases
- C. Measurement and prediction of structure-aware properties for protein residues affected from missense SNPs

## Running the workflow

### Input data

1) Species-specific identifier

2) A variant calling file (.vcf) that follows the required format (see example)

3) A gene annotation file (.gff or .tsv) that follows the required format (see example)

4) All protein sequences for the annotated genes in a single fasta file (.fasta)

5) Diamond database (.dmnd) of protein sequences in AlphaFold Database

6) Diamond database (.dmnd) of protein sequences in ESMatlas Database

*Note: later releases will provide more generalized input files

### Parameters

This table lists all parameters that can be adjusted/modified to run the workflow.

| parameter               | type | details                               | default                        |
| ----------------------- | ---- | ------------------------------------- | ------------------------------ |
| **species**             |      |                                       |                                |
| species_id              | str  | species-specific identifier           |                                |
| **snps_annotation**     |      |                                       |                                |
| path                    | str  | path to snps annotation file          |                                |
| **genes_annotation**    |      |                                       |                                |
| path                    | str  | path to genes annotation file         |                                |
| **protein_db**          |      |                                       |                                |
| path                    | str  | path to protein sequences .fasta file |                                |
| **afdb_diamond_db**     |      |                                                      |                                |
| path                    | str  | path to AlphaFold Diamond database file (.dmnd)      |                                |
| **esm_diamond_db**      |      |                                                      |                                |
| path                    | str  | path to ESMatlas Diamond database file (.dmnd)       |                                |

### Execute the workflow
snakemake --software-deployment-method conda --cores 1 all --config species=101157
