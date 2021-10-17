#!/bin/bash

# $1: Input BAM
# $2: Chromosome sizes text file
# $3: Slop size to use
# $4: Path to blacklist
# $5: Output filename

bedtools bamtobed -i ${1} | \
awk 'BEGIN {OFS = "\t"} ; {if ($6 == "+") print $1, $2 + 4, $2 + 5, $4, $5, $6; else print $1, $3 - 5, $3 - 4, $4, $5, $6}' | \
bedtools slop  -i - -g ${2} -b ${3} | \
bedtools intersect -a - -b ${4} -v | \
pigz > ${5}