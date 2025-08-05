# microbial-snvs-annotation-pipeline

[![Snakemake](https://img.shields.io/badge/snakemake-≥8.0.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/GArcho44/snvs-annotation-pipeline/workflows/Tests/badge.svg?branch=main)](https://github.com/GArcho44/snvs-annotation-pipeline/actions?query=branch%3Amain+workflow%3ATests)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)

This Snakemake-based bioinformatics pipeline functionally annotates microbial single-nucleotide variants (SNVs) in protein-coding genes by leveraging protein-structure information.

It integrates standard metagenomic outputs (SNVs + gene annotations) with high-quality structures from AlphaFold DB (AFDB) and ESMFold, then derives the following structure-aware features for each variant:

- **Relative solvent accessibility (RSA)**
- **Ligand-binding sites (LBS)**
- **Distance to ligand (DTL)**
- **Protein–protein interaction interfaces (PPIIs)**
- **Distance to protein–protein interfaces (DTPPI)**

These features are combined to infer the potential functional impact of each SNV.

*An overview of the pipeline is shown below.*


- [Snakemake workflow: `microbial-snvs-annotation-pipeline`](#snakemake-workflow-name)
  - [Usage](#usage)
  - [Deployment options](#deployment-options)
  - [Authors](#authors)
  - [References](#references)
  - [TODO](#todo)

## Usage

Detailed information about input data and workflow configuration can also be found in the [`config/README.md`](config/README.md).

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this repository or its DOI.

## Deployment options

To run the workflow from command line, change the working directory.

```bash
cd path/to/snakemake-workflow-name
```

Adjust options in the default config file `config/config.yml`.
Before running the complete workflow, you can perform a dry run using:

```bash
snakemake --dry-run
```

## Authors

- Archontis Goumagias
  - University of Groningen, University Medical Center Groningen, Department of Genetics/Pediatrics, Groningen, The Netherlands
  - ORCID profile

## References

> Köster, J., Mölder, F., Jablonski, K. P., Letcher, B., Hall, M. B., Tomkins-Tinch, C. H., Sochat, V., Forster, J., Lee, S., Twardziok, S. O., Kanitz, A., Wilm, A., Holtgrewe, M., Rahmann, S., & Nahnsen, S. _Sustainable data analysis with Snakemake_. F1000Research, 10:33, 10, 33, **2021**. https://doi.org/10.12688/f1000research.29032.2.
