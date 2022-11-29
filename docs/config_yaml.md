# Configuration `.yaml` files

The [`config.yaml`](../inputs/config.yaml) and the cluster [`config.yaml`](../simple/config.yaml) files define the parameters for the run. The input `config.yaml` file must be in the [`inputs`](../inputs) directory relative to the `Snakefile` and is the only required `.yaml` file.  The `config.yaml` file can be specified from the command line using `--configfile` instead of making sure it is placed in the correct directory. The cluster `config.yaml` contains information specific to using Snakemake to submit jobs to the SLURM scheduler and can be found in the [`simple`](../simple/) directory.

## `config.yaml`

The `config.yaml` must contain the following fields:

| Variable         | Description                                            |
|------------------|--------------------------------------------------------|
| `output_dir`     | Path to the output directory                           |
| `idx_bt2`        | Path to the bowtie2 index                              |
| `SAMPLES_TSV`    | Path to the sample meta file in `.tsv` format          |
| `chrom_sizes`    | Path to the chrom sizes `.txt` file                    |
| `blacklist`      | Path to the blacklist file `.bed` file                 |
| `slop`           | Slop size to use in bp                                 |
| `million_factor` | Millions factor to use, i.e. 1000000, 20000000, etc... |
| `keepChr`        | List of chromosomes to keep in the analysis            |
| `species`        | Species of the experiments                             |
| `maxatac_tfs`    | List of TFs available from maxATAC to make predictions |
| `flags`          | Flags used for limiting the workflow steps performed   |

Example:

```yaml
# Output directory
output_dir: /fs/project/PES0738/snakemake/runs/snakeatac_test_tobias

# Bowtie2 index
idx_bt2: /fs/project/PES0738/maxATAC_inputs/genome_inf/hg38/bowtie2_index/hg38

# STAR index: only required if TOBIAS analysis is performed
idx_star: /fs/project/PES0738/maxATAC_inputs/genome_inf/hg38/STAR_hg38_index

# Path to a TSV file with samples and their corresponding FASTQ file paths.
SAMPLES_TSV: '/fs/project/PES0738/snakemake/snakeATAC/inputs/sample.tsv'

# Path to chromosome sizes file
chrom_sizes: "/fs/project/PES0738/maxATAC_inputs/genome_inf/hg38.22XY.chrom.sizes"

# Path to blacklisted regions file
blacklist: "/fs/project/PES0738/maxATAC_inputs/genome_inf/hg38_maxatac_blacklist_V2.bed"

# Slop size to use around the Tn5 insertion sites
slop: 20

# What is the millions factor that you want to use to scale the read counts
million_factor: 20000000

# MACS2 species
species: "hs"

# List of chromosomes to limit the study to
keepChr: 'chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY'

maxatac_tfs:
  TCF7
  TCF12
  LEF1

# Analysis
flags:
  fastqc: True
  maxatac: True
  stranded: False
  tobias: False
  atacseqqc: False
```

## Cluster `config.yaml` for use with SLURM

The cluster [`config.yaml`](../simple/config.yaml) is used to define default cluster configurations if you want to use Snakemake to submit jobs for you. If the rules do not have a resource limit described, they will default to these values when the job is submitted. This information is expanded on in the [jdblischak/smk-simple-slur repo](https://github.com/jdblischak/smk-simple-slurm).

Example:

```yaml
cluster:
  mkdir -p logs/{rule} &&
  sbatch
    --cpus-per-task={threads}
    --mem={resources.mem_mb}
    --job-name=smk-{rule}-{wildcards}
    --account={YOUR_ACCOUNT}
    --output=logs/{rule}/{rule}-{wildcards}-%j.out
    --time={resources.time}
    --parsable
default-resources:
  - time="12:00:00"
  - mem_mb=32000
restart-times: 3
max-jobs-per-second: 1
max-status-checks-per-second: 1
local-cores: 1
latency-wait: 60
jobs: 15
keep-going: True
rerun-incomplete: True
printshellcmds: True
scheduler: greedy
use-conda: True
cluster-status: status-sacct.sh
```
