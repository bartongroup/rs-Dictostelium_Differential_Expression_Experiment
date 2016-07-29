#!/bin/bash

IFS=''
while read -r line
do
var1=$(echo $line | awk '{print $1}')
var2=$(echo $line | awk '{print $2}')
var3=$(echo $line | awk '{print $3}')

eval 'qsub -q c6145.q -cwd -V -b y -pe smp 10 tophat -G DictyEnsemble.gtf -p 10 -I 146 -o' $var1 'DictyEnsemble.fas' $var2 $var3
done <targets
