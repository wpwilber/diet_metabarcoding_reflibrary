This is the project directory for preparing a reference library to be used in the diet metabarcoding 24480 pipeline. This library construction is under early development.

Reproducing this analysis requires Conda and a Unix environment.

This project is designed for HPC deployment on Univa Grid Engine at Notre Dame CRC. It can be configured to run on other HPC clusters by editing the environments and profiles as necessary. 

2026-04-23: Currently I am having a hard time consistently resolving a qiime2 environment with conda from yaml on Notre Dame's HPC. Temporarily I have hardcoded a path to an existing working qiime2 environment stored locally. This environment is not contained in this repo; executing this snakefile therefore requires a locally installed qiime2 environment and editing paths in Snakefile rules that call on qiime2. 

How to execute this project: 

1: Install and activate the contained snakemake environment

```
conda env create -f envs/snakemake-env.yaml
conda activate snakemake-env
```

2: Revise parameters at the top of the snakefile as necessary.

3: Execute the snakefile

```
snakemake --profile profiles/hpc --jobs 1
```

