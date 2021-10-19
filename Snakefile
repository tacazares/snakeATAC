configfile: "./inputs/config.yaml"

# localrules will let the rule run locally rather than submitting to cluster
# computing nodes, this is for very small jobs
localrules: all

import pandas as pd
import os

samples_df = pd.read_table(config['SAMPLES_TSV']).set_index("Sample_ID", drop=False)

ALL_SAMPLES = samples_df["Sample_ID"].unique().tolist()

ALL_FASTQC  = expand(os.path.join(config["output_dir"], "{sample}/qc/fastqc/{types}/{sample}_{read}_fastqc.zip"), sample = ALL_SAMPLES, read = ["1", "2"], types=["pre_trim", "post_trim"])
ALL_TRIMMED_FASTQ = expand(os.path.join(config["output_dir"], "{sample}/trimmed_fastq/{sample}_{read}_val_{read}.fq.gz"), sample = ALL_SAMPLES, read = ["1", "2"])
ALL_BAM = expand(os.path.join(config["output_dir"], "{sample}/alignments/{sample}.bam"), sample = ALL_SAMPLES)
ALL_FLAGSTAT = expand(os.path.join(config["output_dir"], "{sample}/flagstat/{sample}.{type}.flagstat"), sample = ALL_SAMPLES, type= ["sam", "final_bam"])
ALL_BIGWIG = expand(os.path.join(config["output_dir"], "{sample}/tags/{sample}_Tn5_slop" + str(config["slop"]) + "_blacklisted.bw"), sample= ALL_SAMPLES)
ALL_PEAKS = expand(os.path.join(config["output_dir"], "{sample}/peaks/{sample}_ext{ext}_{qval}_peaks.narrowPeak"), sample=ALL_SAMPLES, ext = ["200", "40"], qval=["q05", "q01"])

rule all:
    input:  ALL_BIGWIG  + ALL_FLAGSTAT + ALL_FASTQC + ALL_PEAKS

# Use fastqc to get pre-trim stats on fastq
rule fastqc_pre_trim:
    input: os.path.join(config["fastq_dir"], "{sample}_{read}.fastq.gz")

    output:
        html=os.path.join(config["output_dir"], "{sample}/qc/fastqc/pre_trim/{sample}_{read}.html"),
        zip=os.path.join(config["output_dir"], "{sample}/qc/fastqc/pre_trim/{sample}_{read}_fastqc.zip") # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename

    params: "--quiet"

    log:
        os.path.join(config["output_dir"], "{sample}/logs/fastqc/{sample}_{read}_pre_trim.log")

    threads: 4

    wrapper:
        "0.77.0/bio/fastqc"

# Trim the fastq reads if there is adapter contamination
rule trim_galore_pe:
    input:  
        fq1=os.path.join(config["fastq_dir"], "{sample}_1.fastq.gz"),
        fq2=os.path.join(config["fastq_dir"], "{sample}_2.fastq.gz")

    output:
        os.path.join(config["output_dir"], "{sample}/trimmed_fastq/{sample}_1_val_1.fq.gz"),
        os.path.join(config["output_dir"], "{sample}/trimmed_fastq/{sample}_2_val_2.fq.gz")

    threads: 16

    params: TRIM_DIR = os.path.join(config["output_dir"], "{sample}/trimmed_fastq")

    message: "Trimming adaptors for {input} using trim_galore"

    log: os.path.join(config["output_dir"], "{sample}/logs/trimmed_fastq/{sample}.log")

    conda: "./envs/trim_galore.yaml"

    benchmark: os.path.join(config["output_dir"], "{sample}/logs/benchmark/{sample}_trim_adapter.benchmark")

    shell:
        """
        trim_galore -q 30 -paired -j 4 -o {params.TRIM_DIR} {input.fq1} {input.fq2} 2> {log}
        """

# Get reads statistics post adapter trimming
rule fastqc_post_trim:
    input: os.path.join(config["output_dir"], "{sample}/trimmed_fastq/{sample}_{read}_val_{read}.fq.gz")

    output:
        html=os.path.join(config["output_dir"], "{sample}/qc/fastqc/post_trim/{sample}_{read}.html"),
        zip=os.path.join(config["output_dir"], "{sample}/qc/fastqc/post_trim/{sample}_{read}_fastqc.zip") # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename

    params: "--quiet"

    log:
        os.path.join(config["output_dir"], "{sample}/logs/fastqc/{sample}_{read}_post_trim.log")

    threads: 4

    wrapper:
        "0.77.0/bio/fastqc"

# Align reads with bowtie2
rule bowtie2_align_PE:
    input:
        os.path.join(config["output_dir"], "{sample}/trimmed_fastq/{sample}_1_val_1.fq.gz"), 
        os.path.join(config["output_dir"], "{sample}/trimmed_fastq/{sample}_2_val_2.fq.gz")

    output: temp(os.path.join(config["output_dir"], "{sample}/alignments/{sample}.sam"))

    threads: 8

    params: bowtie = "--very-sensitive --maxins 2000"

    message: "Align reads with bowtie2 {input}: {threads} threads"

    conda: "./envs/bowtie2.yaml"

    benchmark: os.path.join(config["output_dir"], "{sample}/logs/benchmark/{sample}_bowtie2.benchmark")

    log: os.path.join(config["output_dir"], "{sample}/logs/bowtie2/{sample}.log")

    shell:
        """
        bowtie2 {params.bowtie} -p {threads} -x {config[idx_bt2]} -1 {input[0]} -2 {input[1]} -S {output} 2> {log}
        """

# Remove low quality from the file
rule quality_filter_namesort_sam2bam_pe:
    input:  os.path.join(config["output_dir"], "{sample}/alignments/{sample}.sam")

    output: temp(os.path.join(config["output_dir"], "{sample}/alignments/{sample}.namesorted.bam"))

    log:    os.path.join(config["output_dir"], "{sample}/logs/filter_namesort/{sample}.namesorted_bam")

    conda: "./envs/samtools.yaml"

    threads: 8

    benchmark: os.path.join(config["output_dir"], "{sample}/logs/benchmark/{sample}_filter_namesort.benchmark")

    message: "Quality filtering and name sorting BAM {input}: {threads} threads"

    shell:
        """
        samtools view -@ {threads} -b -u -q 30 {input} | \
        samtools sort -@ {threads} -n -o {output} -
        """

# Namesort bam and remove PCR duplicates
rule fixmates_pre_dedup:
    input: os.path.join(config["output_dir"], "{sample}/alignments/{sample}.namesorted.bam")

    output: fixmate_bam=temp(os.path.join(config["output_dir"], "{sample}/alignments/{sample}.fixmate.bam")),
            fixmate_bai=temp(os.path.join(config["output_dir"], "{sample}/alignments/{sample}.fixmate.bam.bai"))

    log:    os.path.join(config["output_dir"], "{sample}/logs/deduplicate/{sample}.fixmate_bam")

    threads: 8

    conda: "./envs/samtools.yaml"

    message: "Fixmates pre PCR deduplication from {input}: {threads} threads"

    shell:
        """
        samtools fixmate -@ {threads} -r -m {input} - 2> {log} | \
        samtools sort -@ {threads} -o {output.fixmate_bam} - 2> {log}
        
        samtools index -@ {threads} {output.fixmate_bam} 2> {log}

        """

# Remove PCR duplicates from BAM
rule remove_PCR_duplicates:
    input: os.path.join(config["output_dir"], "{sample}/alignments/{sample}.fixmate.bam")

    output: dedup_bam=temp(os.path.join(config["output_dir"], "{sample}/alignments/{sample}.deduplicated.bam")),
            dedup_bai=temp(os.path.join(config["output_dir"], "{sample}/alignments/{sample}.deduplicated.bam.bai")),
            final_bam=os.path.join(config["output_dir"], "{sample}/alignments/{sample}.bam"),
            final_bai=os.path.join(config["output_dir"], "{sample}/alignments/{sample}.bam.bai")

    params: config["keepChr"]

    log:    os.path.join(config["output_dir"], "{sample}/logs/deduplicate/{sample}.deduplicate_bam")

    threads: 8

    conda: "./envs/samtools.yaml"

    message: "Removing PCR duplicates and unwanted chromosomes from {input}: {threads} threads"

    shell:
        """
        samtools markdup -@ {threads} -r -s {input} - 2> {log} | \
        samtools sort -@ {threads} -o {output.dedup_bam} - 2> {log}

        samtools index -@ {threads} {output.dedup_bam} 2> {log}
        
        # remove unmapped reads with -F 4
        samtools view -@ {threads} -b -f 3 -F 4 {output.dedup_bam} {params} | \
        samtools sort -@ {threads} -o {output.final_bam} - 2> {log}
        
        samtools index -@ {threads} {output.final_bam} 2> {log}

        """

# Convert alignments to bed file and bedpe file
rule generate_Tn5_sites_bed:
    input: os.path.join(config["output_dir"], "{sample}/alignments/{sample}.bam")

    output: os.path.join(config["output_dir"], "{sample}/tags/{sample}_Tn5_slop" + str(config["slop"]) + "_blacklisted.bed.gz")

    log:    os.path.join(config["output_dir"], "{sample}/logs/bamtobed/{sample}.tagalign")

    threads: 4

    params: chrom_sizes=config["chrom_sizes"],
            slop=config["slop"],
            blacklist=config["blacklist"]

    conda: "./envs/genome_tools.yaml"

    message: "Generate Tn5 sites for {input}: {threads} threads"

    shell:
        """        
        sh ./scripts/shift_reads.sh {input} {params.chrom_sizes} {params.slop} {params.blacklist} {output}
        """

# check number of reads mapped by samtools flagstat
rule flagstat_sam:
    input:  os.path.join(config["output_dir"], "{sample}/alignments/{sample}.sam")

    output: os.path.join(config["output_dir"], "{sample}/flagstat/{sample}.sam.flagstat")

    log:    os.path.join(config["output_dir"], "{sample}/logs/flagstat/{sample}_sorted.log")

    threads: 4

    params: jobname = "{sample}"

    conda: "./envs/samtools.yaml"

    message: "Get flagstats of sam {input}: {threads} threads"

    shell:
        """
        samtools flagstat {input} > {output} 2> {log}
        """

# check number of reads mapped by samtools flagstat
rule flagstat_final_bam:
    input:  os.path.join(config["output_dir"], "{sample}/alignments/{sample}.bam")

    output: os.path.join(config["output_dir"], "{sample}/flagstat/{sample}.final_bam.flagstat")

    log:    os.path.join(config["output_dir"], "{sample}/logs/flagstat/{sample}_final_bam.log")

    threads: 4

    params: jobname = "{sample}"

    conda: "./envs/samtools.yaml"

    message: "Get flagstats of final bam {input}: {threads} threads"

    shell:
        """
        samtools flagstat {input} > {output} 2> {log}
        """

# Generate bigwig coverage track of Tn5 signal scaled to the sequencing depth
rule get_Tn5_coverage:
    input: bam=os.path.join(config["output_dir"], "{sample}/alignments/{sample}.bam"),
           bam_index=os.path.join(config["output_dir"], "{sample}/alignments/{sample}.bam.bai"),
           tn5_bed=os.path.join(config["output_dir"], "{sample}/tags/{sample}_Tn5_slop" + str(config["slop"]) + "_blacklisted.bed.gz"),

    params: millions_factor=config["million_factor"],
            chrom_sizes=config["chrom_sizes"]

    output: tn5_bedgraph=temp(os.path.join(config["output_dir"], "{sample}/tags/{sample}_Tn5_slop" + str(config["slop"]) + "_blacklisted.bedraph")),
            tn5_bigwig=os.path.join(config["output_dir"], "{sample}/tags/{sample}_Tn5_slop" + str(config["slop"]) + "_blacklisted.bw")

    log: os.path.join(config["output_dir"], "{sample}/logs/bigwig/{sample}.bigwig")

    threads: 4

    conda: "./envs/genome_tools.yaml"

    message: "Convert {input} to bigwig: {threads} threads"

    shell:
        """
        sh ./scripts/generate_bigwig.sh {input.bam} {params.millions_factor} {input.tn5_bed} {params.chrom_sizes} {output.tn5_bedgraph} {output.tn5_bigwig}
        """

# Call Peaks with MACS2
rule macs2_call_peaks:
    input: os.path.join(config["output_dir"], "{sample}/tags/{sample}_Tn5_slop" + str(config["slop"]) + "_blacklisted.bed.gz")

    output: os.path.join(config["output_dir"], "{sample}/peaks/{sample}_ext200_q05_peaks.narrowPeak"),
            os.path.join(config["output_dir"], "{sample}/peaks/{sample}_ext40_q05_peaks.narrowPeak"),
            os.path.join(config["output_dir"], "{sample}/peaks/{sample}_ext200_q01_peaks.narrowPeak"),
            os.path.join(config["output_dir"], "{sample}/peaks/{sample}_ext40_q01_peaks.narrowPeak")

    log:    os.path.join(config["output_dir"], "{sample}/logs/macs2/{sample}.macs2")
    params: PEAK_DIR = os.path.join(config["output_dir"], "{sample}/peaks"),
            NAME_ext200_q05 = "{sample}_ext200_q05",
            NAME_ext40_q05 = "{sample}_ext40_q05",
            NAME_ext200_q01 = "{sample}_ext200_q01",
            NAME_ext40_q01 = "{sample}_ext40_q01"

    threads: 4
    conda: "./envs/macs2.yaml"
    message: "call peaks {input}: {threads} threads"
    shell:
        """
        # ext=200 q=.05
        macs2 callpeak -t {input} --name {params.NAME_ext200_q05} -g hs --outdir {params.PEAK_DIR} --nomodel --shift -80 --extsize 200 --keep-dup=all

        # ext=40 q=.05
        macs2 callpeak -t {input} --name {params.NAME_ext40_q05} -g hs --outdir {params.PEAK_DIR} --nomodel --extsize 40 --keep-dup=all
        
        # ext=40 q=.01
        macs2 callpeak -t {input} --name {params.NAME_ext40_q01} -g hs --outdir {params.PEAK_DIR} --nomodel --extsize 40 --keep-dup=all -q 0.01

        # ext=200 q=.01
        macs2 callpeak -t {input} --name {params.NAME_ext200_q01} -g hs --outdir {params.PEAK_DIR} --nomodel --shift -80 --extsize 200 --keep-dup=all -q 0.01
        """