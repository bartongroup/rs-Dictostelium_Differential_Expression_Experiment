#!/bin/bash

IFS=''
while read -r line
do
var1=$(echo $line | awk '{print $1}')

eval 'qsub -q c6145.q -cwd -V -b y samtools sort -n' "$var1/accepted_hits.bam" $var1'_sort'
eval 'qsub -q c6145.q -cwd -V -b y samtools sort ' "$var1/accepted_hits.bam" $var1'-sortIGV'
done <targets
