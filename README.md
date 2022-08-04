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

## Process GM12878 ATAC-seq data from [Corces (2017)](https://www.nature.com/articles/nmeth.4396)

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

4) Test the workflow and scripts are correctly set up by performing a dry-run with the `--dry-run` flag.

```bash
snakemake --cores 14 --use-conda --conda-frontend mamba  --configfile config.yaml --dry-run
```

5) Then you can run the full run using your favorite HPC system. I use LSF and below is an example script:

```bash
snakemake --cores 14 --use-conda --conda-frontend mamba  --configfile config.yaml
```

## Use Snakemake to submit jobs through SLURM

If you want to use Snakemake to submit jobs to slurm, you will need to follow the instruction described by [jdblischak/smk-simple-slur repo](https://github.com/jdblischak/smk-simple-slurm). The directory and scripts are included in this repository, but you will need to adjust the `account` information. You can also adjust any defaults that you wish to use with your job submissions.

Example `.bat` file to drive the snakeATAC workflow

```bash
#!/bin/bash
#SBATCH -D ./outputs
#SBATCH -J dmnd_snake 
#SBATCH -t 96:00:00
#SBATCH --ntasks=8
#SBATCH --mem=16gb
#SBATCH --account={YOUR_ACCOUNT}
#SBATCH --output ./outputs/snakeatac-%j.out
#SBATCH --error ./outputs/snakeatac-%j.err

# Load modules
module load python/3.7-2019.10

# Load the snakemake/mamba env
source activate mamba

# go to a particular directory
cd ./snakeATAC

# make things fail on errors
set -o nounset
set -o errexit
set -x

### run your commands here!
# Develop from the below links
# https://bluegenes.github.io/snakemake-via-slurm/
# https://github.com/jdblischak/smk-simple-slurm

snakemake -s /snakeATAC/Snakefile \
--use-conda \
--conda-frontend mamba \
--configfile /snakeATAC/inputs/config.yaml \
--profile simple/
```