# Meta file

The meta file is a `.tsv` file that has 1 fastq filename prefix per row:

```bash
Sample_ID
SRX298000
SRX298001
```

The code will use these names to find the fastq files that are to be processed. This pipeline assumes that your filenames end with `{filename}_1.fastq.gz` and `{filename}_1.fastq.gz`.
