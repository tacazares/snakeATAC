shell.prefix("set -eo pipefail; echo BEGIN at $(date); ")
shell.suffix("; exitstat=$?; echo END at $(date); echo exit status was $exitstat; exit $exitstat")

configfile: "config.yaml"

# localrules will let the rule run locally rather than submitting to cluster
# computing nodes, this is for very small jobs
localrules: all

samples = pd.read_table(config['SAMPLES_TSV']).set_index("Sample_ID", drop=False)

ALL_SAMPLES = sample["Sample_ID"].unique().tolist()

ALL_FASTQC  = expand(os.join.path(config["output_dir"], "{sample}/fastqc/{sample}_{read}_fastqc.zip", sample = ALL_SAMPLES))
ALL_TRIMMED_FASTQ = expand(os.join.path(config["output_dir"], "{sample}/trimmed_fastq/{sample}_{read}_val_{read}.fq.gz", sample = ALL_SAMPLES, read = ["1", "2"]))
ALL_BAM = expand(os.join.path(config["output_dir"], "{sample}/alignments/{sample}.sorted.bam", sample = ALL_SAMPLES))
ALL_FLAGSTAT = expand(os.join.path(config["output_dir"], "{sample}/alignments/{sample}.sorted.bam.flagstat", sample = ALL_SAMPLES))


rule all:
	input: ALL_TRIMMED_FASTQ + ALL_BAM  + ALL_FLAGSTAT 

rule fastqc_pre_trim:
    input: os.join.path(config["fastq_dir"], "{sample}_{read}.fastq.gz")
    output:
        html=os.path.join(config["output_dir"], "{sample}/qc/fastqc/pre_trim/{sample}_{read}.html"),
        zip=os.path.join(config["output_dir"], "{sample}/qc/fastqc/pre_trim/{sample}_{read}_fastqc.zip") # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params: "--quiet"
    log:
        os.path.join(config["output_dir"], "{sample}/logs/fastqc/{sample}_{read}_pre_trim.log")
    threads: 4
    wrapper:
        "0.77.0/bio/fastqc"

rule trim_galore_pe:
    input:  
        os.join.path(config["fastq_dir"], "{sample}_{read}.fastq.gz"), 
        os.join.path(config["fastq_dir"], "{sample}_{read}.fastq.gz")

    output:
        temp(os.path.join(config["output_dir"], "{sample}/trimmed_fastq/{sample}_1_val_1.fq.gz")),
        temp(os.path.join(config["output_dir"], "{sample}/trimmed_fastq/{sample}_2_val_2.fq.gz"))

    threads: 6
    params: TRIM_DIR = os.path.join(config["output_dir"], "{sample}/trimmed_fastq")
    message: "trimming adaptors for {input} using trim_galore"
    log: os.path.join(config["output_dir"], "{sample}/logs/trimmed_fastq/{sample}.log")
    conda: "./envs/trim_galore.yaml"
    benchmark: os.path.join(config["output_dir"], "{sample}/logs/benchmark/{sample}_trim_adapter.benchmark")
    shell:
        """
        trim_galore -j {threads} -o {params.TRIM_DIR} -paired {input[0]} {input[1]} 2> {log}
        """

rule bowtie2_align_PE:
    input:
        os.path.join(config["output_dir"], "{sample}/trimmed_fastq/{sample}_1_val_1.fq.gz"), 
        os.path.join(config["output_dir"], "{sample}/trimmed_fastq/{sample}_2_val_2.fq.gz")
    
    output: temp(os.path.join(config["output_dir"], "{sample}/alignments/{sample}.sam"))
    
    threads: 8
    params: bowtie = "--very-sensitive --maxins 2000"
    message: "bowtie2_align_PE {input}: {threads} threads"
    conda: "./envs/bowtie2.yaml"
    benchmark: os.path.join(config["output_dir"], "{sample}/logs/benchmark/{sample}_bowtie2.benchmark")
    log: os.path.join(config["output_dir"], "{sample}/logs/bowtie2/{sample}.log")
    shell:
        """
        bowtie2 {params.bowtie} -p {threads} -x {config[idx_bt2]} -1 {input[0]} -2 {input[1]} -S {output} 2> {log}
        """

rule quality_filter_namesort_sam2bam_pe:
    input:  os.path.join(config["output_dir"], "{sample}/alignments/{sample}.sam")

    output: os.path.join(config["output_dir"], "{sample}/alignments/{sample}.sorted.bam")

    log:    os.path.join(config["output_dir"], "{sample}/logs/filter_namesort/{sample}.sorted_bam")
    conda: "./envs/samtools.yaml"
    threads: 8
    benchmark: os.path.join(config["output_dir"], "{sample}/logs/benchmark/{sample}_filter_namesort.benchmark")
    message: "quality_filter_sam2bam {input}: {threads} threads"
    shell:
        """
        samtools view -@ {threads} -b -u -q 30 {input} | \
        samtools sort -@ {threads} -o {output} -
        samtools index -@ {threads} {output}
        """

# check number of reads mapped by samtools flagstat
rule flagstat_bam:
    input:  os.path.join(config["output_dir"], "{sample}/alignments/{sample}.sorted.bam")
    output: os.path.join(config["output_dir"], "{sample}/alignments/{sample}.sorted.bam.flagstats")
    log:    os.path.join(config["output_dir"], "{sample}/logs/flagstat/{sample}_{read}_pre_trim.log")
    threads: 1
    params: jobname = "{sample}"
    message: "flagstat_bam {input}: {threads} threads"
    shell:
        """
        samtools flagstat {input} > {output} 2> {log}
        """
