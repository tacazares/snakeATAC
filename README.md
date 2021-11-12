# snakeATAC

A snakemake workflow for ATAC-seq data processing

This pipeline was developed from code used in https://github.com/tacazares/pyflow-ATACseq

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

In your favorite repo directory, clone `[snakeATAC](https://github.com/tacazares/snakeATAC.git)`. You will modify the `[config.yaml](docs/config_yaml.md)` and create a `[tab-delimited sample meta file](docs/meta_file.md)`.

## Workflow Overview

Snakemake pipelines promote experimental reproducibility. For each project that you have, you should have a seperate `[config.yaml](docs/config_yaml.md)` file and a unique input/output directory.

This workflow assumes that your run parameters are stored in the `[config.yaml](docs/config_yaml.md)` file and meta data for the experiments are found in a `[tab-delimited sample meta file](docs/meta_file.md)`.


## Run snakeATAC

First you must adjust the [config.yaml](docs/config_yaml.md)` and the `[tab-delimited sample meta file](docs/meta_file.md)` for your specific experiment. Then, you must change to the directory containing the `Snakefile` and execute the following command where `{threads}` is the # of threads available:

```bash
snakemake --cores {threads} --use-conda --conda-frontend conda
```
