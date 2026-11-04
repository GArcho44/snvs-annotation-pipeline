# proStAMP

**pro**tein **st**ructure-based **a**nnotation of **m**icrobial single nucleotide **p**olymorphisms

[![Snakemake](https://img.shields.io/badge/Snakemake-%E2%89%A58.0.0-brightgreen.svg)](https://snakemake.github.io/)
![Development status](https://img.shields.io/badge/status-under%20development-2A9585)
[![Run with Conda](https://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/)
[![GitHub Actions](https://github.com/GArcho44/snvs-annotation-pipeline/workflows/Tests/badge.svg?branch=main)](https://github.com/GArcho44/snvs-annotation-pipeline/actions)

**proStAMP** is a Snakemake-based workflow for annotating microbial single nucleotide polymorphisms (SNPs) in protein-coding genes and characterizing their positions in predicted protein structures. Starting from a reference genome and an SNP VCF, it combines gene and variant annotation, protein-to-structure mapping, and residue-level structural analyses into a species-specific annotation table.

## Workflow overview

![proStAMP workflow overview](docs/images/proStAMP-v0.0.3.svg)

The workflow comprises these steps:

1. **Genome and variant annotation:** [Bakta](https://github.com/oschwengers/bakta) annotates the reference genome, while [SnpEff](https://github.com/pcingola/SnpEff) assigns SNV consequences. Synonymous and missense SNVs, gene annotations, and corresponding protein sequences are selected.
2. **Protein-to-structure mapping:** [DIAMOND](https://github.com/bbuchfink/diamond) searches protein sequences against local databases representing [AlphaFold DB](https://alphafold.ebi.ac.uk/) and [ESMAtlas](https://esmatlas.com/). Stringent alignment criteria are applied to select suitable structural references.
3. **Structure retrieval:** Selected models are downloaded through database APIs or extracted from local [Foldcomp](https://github.com/steineggerlab/foldcomp) databases, according to the configuration.
4. **Structural feature annotation:** [DSSP](https://github.com/PDB-REDO/dssp) provides solvent accessibility information for RSA calculation; [Fpocket](https://github.com/Discngine/fpocket) with [PRANK/P2Rank](https://github.com/rdk/p2rank) identifies and rescores candidate ligand-binding pockets; and [PeSTo](https://github.com/LBM-EPFL/PeSTo) predicts protein interaction interfaces.
5. **Feature integration:** Gene and SNV annotations, protein mappings, and structural descriptors are combined in a final per-species TSV file.

### Structural annotations

- **Relative solvent accessibility (RSA):** Degree of residue exposure to solvent.
- **Ligand-binding pockets (LBP):** Candidate pocket locations identified by Fpocket and rescored by PRANK.
- **Distance to ligand-binding pockets (DTLBP):** Minimum 3D distance from a residue to a selected predicted pocket center, in Å.
- **Predicted protein–protein interaction interfaces:** Candidate interface regions inferred by PeSTo.
- **Distance to protein–protein interaction interfaces (DTPPI):** 3D distance to selected predicted interface residues, in Å; residues assigned to an interface have DTPPI = 0.

**Main output:** `data/merged_features/{species}/snv_structural_features_merged.tsv`

## Software and resources

proStAMP integrates the following tools and reference resources. Please consult their repositories for documentation and appropriate software citations.

| Tool or resource | Role in proStAMP | GitHub / resource |
| --- | --- | --- |
| **Snakemake** | Workflow orchestration and reproducibility | [snakemake/snakemake](https://github.com/snakemake/snakemake) |
| **Bakta** | Bacterial genome annotation | [oschwengers/bakta](https://github.com/oschwengers/bakta) |
| **SnpEff** | Variant consequence annotation | [pcingola/SnpEff](https://github.com/pcingola/SnpEff) |
| **DIAMOND** | Protein sequence alignment against structural reference databases | [bbuchfink/diamond](https://github.com/bbuchfink/diamond) |
| **AlphaFold DB** | Predicted structures of UniProt-associated proteins | [Database](https://alphafold.ebi.ac.uk/) · [AlphaFold/AFDB resources](https://github.com/google-deepmind/alphafold/tree/main/afdb) |
| **ESMAtlas / ESMFold** | Predicted structures of metagenomic proteins | [Atlas](https://esmatlas.com/) · [facebookresearch/esm](https://github.com/facebookresearch/esm) |
| **Foldcomp** *(optional retrieval backend)* | Selective extraction from locally stored compressed structure databases | [steineggerlab/foldcomp](https://github.com/steineggerlab/foldcomp) |
| **DSSP** | Secondary structure and solvent accessibility assignment | [PDB-REDO/dssp](https://github.com/PDB-REDO/dssp) |
| **Fpocket** | Candidate ligand-binding pocket detection | [Discngine/fpocket](https://github.com/Discngine/fpocket) |
| **PRANK / P2Rank** | Rescoring and ranking Fpocket predictions | [rdk/p2rank](https://github.com/rdk/p2rank) |
| **PeSTo** | Predicted protein–protein interaction interfaces | [LBM-EPFL/PeSTo](https://github.com/LBM-EPFL/PeSTo) |
| **Python scientific libraries** | Sequence/structure parsing, statistics, and data integration | [Biopython](https://github.com/biopython/biopython) · [pandas](https://github.com/pandas-dev/pandas) · [NumPy](https://github.com/numpy/numpy) · [SciPy](https://github.com/scipy/scipy) |

## Usage

The workflow processes one species/genome at a time, using the identifier supplied in the configuration. The main inputs are a **reference genome FASTA** and an **SNV VCF** (`.vcf` or `.vcf.gz`). Local reference databases and the software resources listed below must also be configured.

Edit `config/config.yml` with your species identifier, input files, and database paths. The following illustrates the principal settings (paths are placeholders):

```yaml
species: "100344"
genome_fasta: "/path/to/reference_genome.fasta"
snvs_vcf: "/path/to/variants.vcf.gz"
bakta_db: "/path/to/bakta_db"
afdb_diamond_db: "/path/to/afdb.dmnd"
esm_diamond_db: "/path/to/esmatlas.dmnd"

perc_identity_threshold: 98
plddt_threshold: 80
structure_retrieval: "api"  # "api" or "foldcomp"

# Required only when structure_retrieval: "foldcomp"
foldcomp_afdb_db: "/path/to/afdb_foldcomp_prefix"
foldcomp_esmatlas_db: "/path/to/esmatlas_foldcomp_prefix"
```

**Reference setup:** Provide a species-specific `snpEff.config` in `data/snpeff_preparation/{species}/` with the correct genome identifier and bacterial genetic code, as described in the pipeline configuration guide. The Bakta reference database, DIAMOND databases, and local PeSTo and P2Rank/PRANK installations must also be available. API retrieval requires network access; Foldcomp retrieval requires locally installed Foldcomp databases.

For detailed configuration instructions, see [`config/README.md`](config/README.md).

### Choosing a structure retrieval method

Change **one parameter** in `config/config.yml`:

```yaml
structure_retrieval: "api"
```

or

```yaml
structure_retrieval: "foldcomp"
```

Both methods use the same selected protein references and produce the structure directories expected by downstream analyses. When switching methods for an already processed species, rerun the retrieval steps so previous outputs are not inadvertently reused.

## Deployment options

Install [Snakemake](https://snakemake.readthedocs.io/) (version 8 or later) and a Conda-compatible package manager, then clone the repository:

```bash
git clone https://github.com/GArcho44/snvs-annotation-pipeline.git
cd snvs-annotation-pipeline
```

Configure the required inputs and databases in `config/config.yml`. Run the workflow from its `workflow/` directory, where the species-specific output paths are defined:

```bash
cd workflow

# Preview the workflow for species 100344
snakemake -n -p \
    data/merged_features/100344/snv_structural_features_merged.tsv

# Run with Conda-managed environments
snakemake --cores 4 --sdm conda \
    data/merged_features/100344/snv_structural_features_merged.tsv
```

The example target and species identifier must agree with `config/config.yml`. HPC installations may additionally require site-specific module configuration (for example, Java for PRANK) and sufficient local storage for the structural reference databases.

## Authors

- **Archontis Goumagias** — University of Groningen, University Medical Center Groningen, Department of Genetics/Pediatrics, Groningen, The Netherlands. [ORCID: 0009-0007-1430-6210](https://orcid.org/0009-0007-1430-6210)

## References

If you use proStAMP, please cite this repository (or its DOI, if available) and the relevant original publications for the tools used in your analysis.
