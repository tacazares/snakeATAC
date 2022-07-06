# Meta file

The meta file is a `.tsv` file that has 1 fastq filename prefix per row.

The code will use the filename prefix to find the fastq files that are to be processed. This pipeline assumes that your filenames end with `{filename}_1.fastq.gz` and `{filename}_1.fastq.gz`.

```tsv
Sample_ID
SRX298000
SRX298001
```

These files should be located in the directory that is being pointed to by the `fastq_dir` variable in the `config.yaml` file: 

```yaml
# FASTQ directory
fastq_dir: ./inputs/fastq
```