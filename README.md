# Dictostelium_Differential_Expression_Experiment
Identify Differential Expressed Genes in DgcA null mutant and Knock-out samples in Dictyostelium Discoideum. The same workflow has been used for the differential expression analysis of Polysphondylium pallidum.

## Quality Control

  R CMD BATCH QualityControl.R &

## Alignmnet using Tophat

  ./try.sh
  ./try1.sh
  ./try2.sh
  
## Alignmnet Visualization

   R CMD BATCH InputIGV-HtseqCount.R

## Read Count

  R CMD BATCH htseq-count.R

## Differential Gene Expression

  

## Gene-Set Enrichment Analysis
