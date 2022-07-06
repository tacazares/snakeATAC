#!/bin/bash

# ${1} Tn5 sites
# ${2} Peaks

# Use wc to count the number of rows
read_counts=$(gunzip -c ${1} | wc -l)

# Sort and merge input peaks. Then intersect with list of tags
# Bedtools docs for -u: Write original A entry (tag) once if any overlaps found in B (peaks). In other words, just report the fact at least one overlap was found in B
reads_in_peaks=$(bedtools sort -i ${2} | bedtools merge -i - | bedtools intersect -u -a ${1} -b - | wc -l)

echo ${reads_in_peaks} > ${3}
echo ${read_counts} >> ${3}

echo $(bc -l <<< "${reads_in_peaks}/${read_counts}") >> ${3}