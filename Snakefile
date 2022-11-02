"""
~~~~~~~~~~ snakeATAC ~~~~~~~~~~~
Description:
A snakemake workflow for processing ATAC-seq data. This workflow was
developed using code from TOBIAS, HINT-ATAC, maxATAC, and many more.

See detailed introduction: https://github.com/tacazares/snakeATAC

Inputs:
config.yaml: run parameters
sample.tsv: experimental meta data wil columns (condition, srx, srr)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Environment Setup
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import os
from scripts.meta import MetaTable

# Default path to sample config file. This can be overwritten by the snakemake CLI
configfile: "./inputs/config.yaml"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load Meta
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Initialize the meta table. This will be used later to help build filenames
meta = MetaTable(meta_path=config["SAMPLES_TSV"], output_dir=config["output_dir"])

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Build Output Filenames
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

outputs = []

ALL_TRIMMED_FASTQ = expand(os.path.join(config["output_dir"], "{condition}", "replicate_data/{sample}/trimmed_fastq/{sample}_{read}_val_{read}.fq.gz"), zip, condition = meta.name_df.condition, sample = meta.name_df.srx, read = meta.name_df.read)

outputs.append(ALL_TRIMMED_FASTQ)

# If the fastqc argument is true generate the files. This will currently perform pre and post trimming QC. 
if config["flags"]["fastqc"]:
    # Build all fastqc names for pre and post trimmed files
    ALL_PRE_TRIM_FASTQC  = expand(os.path.join(config["output_dir"], "{condition}", "replicate_data/{sample}/qc/fastqc/pre_trim/{sample}_{read}_fastqc.zip"), zip, condition = meta.name_df.condition, sample = meta.name_df.srx, read = meta.name_df.read)
    outputs.append(ALL_PRE_TRIM_FASTQC)

    ALL_POST_TRIM_FASTQC  = expand(os.path.join(config["output_dir"], "{condition}", "replicate_data/{sample}/qc/fastqc/post_trim/{sample}_{read}_fastqc.zip"), zip, condition = meta.name_df.condition, sample = meta.name_df.srx, read = meta.name_df.read)
    outputs.append(ALL_POST_TRIM_FASTQC)

## Build BAM, BIGWIG, and peak names
if config["flags"]["maxatac"]:
    #ALL_SAM = expand(os.path.join(config["output_dir"], "{condition}/replicate_data/{sample}/alignments/bowtie2/{sample}.sam"), zip, condition = meta.name_df.condition, sample = meta.name_df.srx)
    #ALL_SORTED_BAM = expand(os.path.join(config["output_dir"], "{condition}/replicate_data/{sample}/alignments/bowtie2/{sample}.namesorted.bam"), zip, condition = meta.name_df.condition, sample = meta.name_df.srx)
    ALL_BAM = expand(os.path.join(config["output_dir"], "{condition}/replicate_data/{sample}/alignments/bowtie2/{sample}.bam"), zip, condition = meta.name_df.condition, sample = meta.name_df.srx)
    #ALL_CUT_SITES = expand(os.path.join(config["output_dir"], "{condition}/replicate_data/{sample}/cut_sites/{sample}_Tn5_slop" + str(config["slop"]) + "_blacklisted.bed.gz"), zip, condition = meta.name_df.condition, sample = meta.name_df.srx)
    ALL_BIGWIG = expand(os.path.join(config["output_dir"], "{condition}/replicate_data/{sample}/rpm_signal/{sample}_Tn5_slop" + str(config["slop"]) + "_blacklisted.bw"), zip, condition = meta.name_df.condition, sample = meta.name_df.srx)
    ALL_MAXATAC_PEAKS = expand(os.path.join(config["output_dir"], "{condition}/replicate_data/{sample}/peaks/maxatac/{sample}_ext40_q05_peaks.narrowPeak"), zip, condition = meta.name_df.condition, sample = meta.name_df.srx)
    ALL_SAM_FLAGSTAT = expand(os.path.join(config["output_dir"], "{condition}/replicate_data/{sample}/qc/maxatac/flagstat/{sample}.sam.txt"), zip, condition = meta.name_df.condition, sample = meta.name_df.srx)
    ALL_BAM_FLAGSTAT = expand(os.path.join(config["output_dir"], "{condition}/replicate_data/{sample}/qc/maxatac/flagstat/{sample}.final_bam.txt"), zip, condition = meta.name_df.condition, sample = meta.name_df.srx)

    ALL_FRIPS = expand(os.path.join(config["output_dir"], "{condition}/replicate_data/{sample}/qc/maxatac/frip/{sample}_ext40_q05_FRIP.txt"), zip, condition = meta.name_df.condition, sample = meta.name_df.srx)
    ALL_FRAG_DIST = expand(os.path.join(config["output_dir"], "{sample}/qc/maxatac/fragment_dist/{sample}.tsv"), zip, condition = meta.name_df.condition, sample = meta.name_df.srx)
    ALL_CHROM_READ_COUNTS = expand(os.path.join(config["output_dir"], "{sample}/qc/maxatac/chrom_counts/{sample}_chrom_read_counts.txt"), zip, condition = meta.name_df.condition, sample = meta.name_df.srx)

    outputs.append(ALL_BAM + ALL_BIGWIG + ALL_MAXATAC_PEAKS + ALL_SAM_FLAGSTAT + ALL_BAM_FLAGSTAT + ALL_FRIPS + ALL_FRAG_DIST + ALL_CHROM_READ_COUNTS)

# Run ATAC-seq qc
if config["flags"]["atacseqqc"]:
    ALL_ATACSEQQC = expand(os.path.join(config["output_dir"], "{sample}/logs/ATACseqQC/atacseqqc.txt"), sample=ALL_SAMPLES)
    outputs.append(ALL_ATACSEQQC)

# Generate stranded signal outputs
if config["flags"]["stranded"]:
    ALL_STRAND = expand(os.path.join(config["output_dir"], "{condition}/replicate_data/{sample}/cut_sites/{sample}_Tn5_slop" + str(config["slop"]) + "_{strand}Strand_blacklisted.bed.gz"), sample=meta.sample_list, strand = ["pos", "neg"])
    outputs.append(ALL_STRAND)

shift_dict = {"40": "0", "200": "80"}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Import functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

include: os.path.join("./snakefiles", "process_fastq.sf")
include: os.path.join("./snakefiles", "maxatac.sf")
include: os.path.join("./snakefiles", "tobias.sf")
#include: os.path.join("./snakefiles", "replicate_analysis.sf")
#include: os.path.join("./snakefiles", "atacseqqc.sf")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Start Run
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rule all:
    input:  outputs
