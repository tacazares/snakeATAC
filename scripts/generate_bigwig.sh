#!/bin/bash

# $1: Input BAM
# $2: Millions factor
# $3: Tn5 bed
# $4: Chromosome sizes file
# $5: Output bedgraph
# $6: Output bigwig

mapped_reads=$(samtools view -c -F 260 ${1})
reads_factor=$(bc -l <<< "1/${mapped_reads}")
rpm_factor=$(bc -l <<< "${reads_factor} * ${2}")

bedtools genomecov -i ${3} -g ${4} -bg -scale ${rpm_factor} | LC_COLLATE=C sort -k1,1 -k2,2n > ${5}

bedGraphToBigWig ${5} ${4} ${6}
