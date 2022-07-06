# snakeATAC

A snakemake workflow for ATAC-seq data processing

This pipeline was created from code developed by [Crazy Hot Tommy!](https://github.com/crazyhottommy?tab=repositories)

## Installation

This pipeline uses Anaconda.

### Create environment

Create a `conda` environment with:

```bash
conda create -n snakemake python=3.9
```

Activate the environment:

```bash
conda activate snakemake
```

### Install snakemake

Pip install snakemake:

```bash
pip install snakemake
```

### Clone the snakeATAC repository

In your favorite repo directory clone the snakeATAC repo:

```bash
git clone https://github.com/tacazares/snakeATAC.git
```

## Workflow Overview

Snakemake pipelines promote experimental reproducibility. For each project that you have, you should have a seperate [config.yaml](docs/config_yaml.md), [tab-delimited sample meta file](docs/meta_file.md), and a unique output directory.

This workflow assumes that your run parameters are stored in the [config.yaml](docs/config_yaml.md) file and meta data for the experiments are found in a [tab-delimited sample meta file](docs/meta_file.md).

You will need to modify the [config.yaml](docs/config_yaml.md) and create a [tab-delimited sample meta file](docs/meta_file.md) before running the pipeline.

A detailed overview of the steps in the ATAC-seq data processing are found [here](docs/ATAC_processing.md).

## Run snakeATAC

1) Adjust the [config.yaml](docs/config_yaml.md) and the [tab-delimited sample meta file](docs/meta_file.md) for your specific experiment.
2) Change to the directory containing the `Snakefile` and execute the following command where `{threads}` is the # of threads available:

```bash
snakemake --cores {threads} --use-conda --conda-frontend conda
```

## Process GM12878 ATAC-seq data from [Corces (2017)](https://www.nature.com/articles/nmeth.4396) and [Buenrostro (2013)](https://www.nature.com/articles/nmeth.2688)

This is a brief example of how to process ATAC-seq data from public sources. 

First, create a directory to store your fastq files. The SRA accession list for public GM12878 data is available at [./snakeATAC/inputs/GM12878_sample.tsv](./inputs/GM12878_sample.tsv).

### Download fastq files

First, download the fastq files from sra using `fasterq-dump`. [This](https://rnnh.github.io/bioinfo-notebook/docs/fasterq-dump.html) is a good overview for new users.

Below is example code of how to download data from SRA using `fasterq-dump` and `pigz` for file compression.

```bash
# Make the directory where the fastq files will be stored
mkdir -p ./data/fastq

# Change into directory
cd ./data/fastq

# Loop through sample names and downlad
for SRA in $(tail -n +2 ./inputs/GM12878_sample.tsv | cut -f1);
do
# Fastq dump
fasterq-dump -e 6 -p ${SRA}

# Compress fastq file
pigz  ${SRA}*.fastq
done
```

### Run snakemake

First, if you are running this pipeline for your first time, you will need to probably do a dry-run to make sure that everything was installed right. In order to do the dry-run, you will need to have your [`config.yaml`](docs/config_yaml.md) file and [`sample.tsv`](docs/meta_file.md) file correctly set up.

1) Change to the working directory for snakeATAC:

   ```bash
   cd ./snakeATAC/
   ```

2) Then you can run the full run using your favorite HPC system. I use LSF and below is an example script:

   ```bash
   snakemake --cores 24 --use-conda --conda-frontend conda
   ```
