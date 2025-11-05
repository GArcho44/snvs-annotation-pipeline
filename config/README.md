## Workflow overview

This workflow is a best-practice workflow for `<structure-aware annotation of microbial SNPs>`.
The workflow is built using [snakemake](https://snakemake.readthedocs.io/en/stable/) and consists of the following 3 main configurable modules:

A. Sequence-based functional annotation of SNPs & Genes 
B. retrieval of high-confidence structure models from AFDB and ESMatlas protein structure databases
C. Measurement and prediction of structure-aware properties for protein residues affected from missense SNPs

## Running the workflow

### Input data

A) A variant calling file (.vcf) that follows the required format (see Experimental Design)

B) Species reference genome(s)

### Parameters

This table lists all parameters that can be adjusted/modified to run the workflow.

| parameter          | type | details                               | default                        |
| ------------------ | ---- | ------------------------------------- | ------------------------------ |
| **species**        |      |                                       |                                |
| species_id               | str  | species specific identifier, mandatory|                                |
| **ref_genome**     |      |                                 |                                |
| path           | str  | path to reference genome file             |                                |
| **conf_score** |      |                                       |                                |
| confidence_score        | num  | pLDDT confidence score filter          | 80                            |
