#!/bin/bash

IFS=''
while read -r line
do
var1=$(echo $line | awk '{print $1}')

eval 'qsub -q c6145.q -cwd -V -b y -pe smp 10 samtools view -o' $var1'_sort.sam' $var1'_sort.bam'
eval 'qsub -q c6145.q -cwd -V -b y -pe smp 10 samtools index' $var1'-sortIGV.bam'
done <targets
