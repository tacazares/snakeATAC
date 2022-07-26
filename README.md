# snakeATAC

A snakemake workflow for ATAC-seq data processing

This pipeline was created from code developed by [Crazy Hot Tommy!](https://github.com/crazyhottommy?tab=repositories)

## Installation

This pipeline uses Anaconda and Snakemake. Follow the [Snakemake install instructions](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) for the best experience. Below is a brief overview of how to install Snakemake.

### Create environment

Create a `conda` environment and download `mamba`:

```bash
conda install -n snakemake -c conda-forge -bioconda mamba snakemake
```

Activate the `snakemake` environment:

```bash
conda activate snakemake
```

### Clone the snakeATAC repository

In your favorite directory clone the snakeATAC repo:

```bash
git clone https://github.com/tacazares/snakeATAC.git
```

## Workflow Overview

Snakemake pipelines promote experimental reproducibility. For each project that you have, you should have a seperate [config.yaml](docs/config_yaml.md), [tab-delimited sample meta file](docs/meta_file.md), and a unique output directory.

This workflow assumes that your run parameters are stored in the [config.yaml](docs/config_yaml.md) file and meta data for the experiments are found in a [tab-delimited sample meta file](docs/meta_file.md).

You will need to modify the [config.yaml](docs/config_yaml.md) and create a [tab-delimited sample meta file](docs/meta_file.md) before running the pipeline.

A detailed overview of the steps in the ATAC-seq data processing are found [here](docs/ATAC_processing.md).

## Process GM12878 ATAC-seq data from [Corces (2017)](https://www.nature.com/articles/nmeth.4396) and [Buenrostro (2013)](https://www.nature.com/articles/nmeth.2688)

This is a brief example of how to process ATAC-seq data from public sources. The SRA accession list for public GM12878 data is available at [./snakeATAC/inputs/GM12878_sample.tsv](./inputs/GM12878_sample.tsv).

If you are running this pipeline for your first time, you will need to install all the `conda` environments used and perform a dry-run to make sure that everything was installed right.

1) Adjust the [config.yaml](docs/config_yaml.md) and the [tab-delimited sample meta file](docs/meta_file.md) for your specific experiment.

2) Change to the working directory for snakeATAC:

```bash
cd ./snakeATAC/
```

3) Next, use the `--conda-create-envs-only` flag to create the environments.

```bash
snakemake --cores 14 --use-conda --conda-frontend mamba --conda-create-envs-only --configfile config.yaml
```

4) Test the workflow and scripts are correctly stet up by performing a dry-run with the `--dry-run` flag.

```bash
snakemake --cores 14 --use-conda --conda-frontend mamba  --configfile config.yaml --dry-run
```

5) Then you can run the full run using your favorite HPC system. I use LSF and below is an example script:

```bash
snakemake --cores 14 --use-conda --conda-frontend mamba  --configfile config.yaml
```
