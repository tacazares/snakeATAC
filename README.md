# snakeATAC

Yet another snakemake workflow for ATAC-seq data processing. This pipeline was created from code developed by:

* [Crazy Hot Tommy!](https://github.com/crazyhottommy?tab=repositories)'s many instructional guides
* [TOBIAS](https://github.molgen.mpg.de/loosolab/TOBIAS_snakemake) ATAC-seq footprinting Snakemake workflow

For SLURM setup we reference:

* [jdblischak/smk-simple-slur repo](https://github.com/jdblischak/smk-simple-slurm) for simple submitting snakemake on SLURM
* [Tessa Pierce](https://bluegenes.github.io/snakemake-via-slurm/) blog for example templates

## Workflow Overview

Snakemake pipelines promote experimental reproducibility. For this project, you should have the following inputs customized for your analysis:

1) A [config.yaml](docs/config_yaml.md) that describes the run parameters and location of reference data.
2) A [tab-delimited sample meta file](docs/meta_file.md) file that describes the experiments to download from SRA and how to group them.
3) A unique output directory.

A detailed overview of the steps in the ATAC-seq data processing are found [on the maxATAC wiki site]([docs/ATAC_processing.md](https://github.com/MiraldiLab/maxATAC/wiki/ATAC-seq-Data-Processing)).

### Data Processing

* Trim Galore! for `.fastq` QC and adapter trimming.
* Bowtie2 for read alignment

### Data QC

### TFBS Analysis

## Installation

This pipeline uses Anaconda and Snakemake. Follow the [Snakemake install instructions](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) for the best experience. Below is a brief overview of how to install Snakemake.

### Create environment

Create a `conda` environment and download `mamba`:

```bash
conda create -n snakeatac -c conda-forge -c bioconda mamba snakemake
```

Activate the `snakeatac` environment:

```bash
conda activate snakeatac
```

### Clone the snakeATAC repository

In your favorite directory clone the snakeATAC repo:

```bash
git clone https://github.com/tacazares/snakeATAC.git
```

### Set up run-specific parameters

If you are running this pipeline for your first time, you will need to install all the `conda` environments used and perform a dry-run to make sure that everything was installed right.

1) Adjust the [config.yaml](docs/config_yaml.md) and the [tab-delimited sample meta file](docs/meta_file.md) for your specific experiment.

2) Change to the working directory for snakeATAC. By default, Snakemake will look for a file called `Snakefile` with the rules and run information. You can use a custome `Snakefile` with the `-s` flag followed by the path to the file.

   ```bash
   cd ./snakeATAC/
   ```

3) Next, use the `--conda-create-envs-only` flag to create the environments.

   ```bash
   snakemake --cores 14 --use-conda --conda-frontend mamba --conda-create-envs-only --configfile ./inputs/config.yaml
   ```

4) Test the workflow and scripts are correctly set up by performing a dry-run with the `--dry-run` flag.

   ```bash
   snakemake --cores 14 --use-conda --conda-frontend mamba  --configfile ./inputs/config.yaml --dry-run
   ```

### Test snakeATAC

The [./snakeATAC/inputs/GM12878_sample.tsv](./inputs/GM12878_sample.tsv) contains information for a test run to process GM12878 OMNI ATAC-seq data.

After install, you can run the full run using your favorite HPC system.

```bash
snakemake --cores 14 --use-conda --conda-frontend mamba  --configfile ./inputs/config.yaml
```

## Use Snakemake to submit jobs through SLURM

If you want to use Snakemake to submit jobs to slurm, you will need to follow the instruction described by [jdblischak/smk-simple-slur repo](https://github.com/jdblischak/smk-simple-slurm). The directory and scripts are included in this repository, but you will need to adjust the `account` information. You can also adjust any defaults that you wish to use with your job submissions. NOTE: You will need to use `chmod +x status-sacct.sh` to make the script executable.

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
