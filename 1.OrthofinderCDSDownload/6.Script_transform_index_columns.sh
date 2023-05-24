#!/bin/bash
awk '{print $1,$2,$3,$4}' *_df.txt > genomes_index_geneID.txt
awk -v OFS="\t" '$1=$1' genomes_index_geneID.txt > genomes_index_geneID_tab.txt
rm genomes_index_geneID.txt

