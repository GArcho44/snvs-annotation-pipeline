# Environments Overview

This folder contains **Conda environments** used by the workflow. Each environment is used in specified rules (Table 1), keeping dependencies minimal and builds reproducible.

- Use `conda env create -f <env>.yaml` to create locally.
- In the pipeline, Snakemake activates environments per rule via  
  `conda: "envs/<env>.yaml"`.

## Table 1 - Environments and usage

| Environment name   | Rules IDs activated                         | Key packages |
|--------------------|---------------------------------------------|--------------|
| `python_env.yaml`  | 1, 4, 5, 6, 9, 13, 15, 16, 17, 18, 19          | `python=3.9`, `biopython`, `requests`, `matplotlib`, `threadpoolctl`, `pandas`, `numpy` |
| `diamond.yaml`     | 2, 3                                         | `diamond=2.0.15` |
| `dssp_env.yaml`    | 8                                           | `dssp=3.0.0`, `libboost=1.73.0` |
| `fpocket_env.yaml` | 12                                          | `fpocket=4.2.1` |
| `pesto_env.yaml`   | 14                                          | `python=3.9`, `pytorch`, `torchvision`, `torchaudio`, `numpy=1.22.4`, `scipy`, `scikit-learn`, `pandas`, `biopython`, `mdtraj`, `h5py`, `gemmi`, `tqdm`, `matplotlib`, *(pip)* `tensorboard` |
