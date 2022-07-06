# config.yaml

The `config.yaml` file defines the parameters for the run. The input `config.yaml` file must be in the inputs directory relative to the `Snakefile`. If you clone this directory from github, you will only need to modify the paths in the `config.yaml` file found in the `./inputs` directory.

## Required Fields

The `config.yaml` must contain the following fields:

| variable         | description                                            |
|------------------|--------------------------------------------------------|
| `fastq_dir`      | Path to directory containing fastq files               |
| `output_dir`     | Path to the output directory                           |
| `idx_bt2`        | Path to the bowtie2 index                              |
| `SAMPLES_TSV`    | Path to the sample meta file in `.tsv` format          |
| `chrom_sizes`    | Path to the chrom sizes `.txt` file                    |
| `blacklist`      | Path to the blacklist file `.bed` file                 |
| `slop`           | Slop size to use in bp                                 |
| `million_factor` | Millions factor to use, i.e. 1000000, 20000000, etc... |
| `keepChr`        | List of chromosomes to keep in the analysis            |
| `species`        | Species of the experiments                             |

Example:

```bash
# FASTQ directory
fastq_dir: ./inputs/fastq
output_dir: ./outputs

# Bowtie2 index
idx_bt2: ./hg38/bowtie2_index/hg38

# Path to a JSON file with samples and their corresponding FASTQ files.
SAMPLES_TSV: './ATAC_test_sample.tsv'

# Path to chromosome sizes file
chrom_sizes: "./genome_inf/hg38.22XY.chrom.sizes"

# Path to blacklisted regions file
blacklist: "./genome_inf/hg38_maxatac_blacklist_merged.bed"

# Slop size to use around the Tn5 insertion sites
slop: 20

# Species
species: hs

# What is the millions factor that you want to use to scale the read counts
million_factor: 1000000

# List of chromosomes to limit the study to
keepChr: 'chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22'
```
