# genome-scripts
Some hodge podge scripts for dealing with genome stuff

#Transcriptomes
This folder contains instructions on the KU BI transcriptome pipeline

#Other stuff
linebylineblast.sh (and the associated R-scripts: onelining_tempseq.R and linebyline.R) blast an assembly against itself to look for regions of an assembly that made have been subject to a "false tandem duplication" (e.g. instead of alleles being assembled together, they have been assembled end to end).

summarizing_blast_hits.R summarizes the results of linebylineblast.sh

length_dist.R can be used to get the lengths of the contigs/scaffolds in the assembly and to generate histogram summaries of this data, as well as the percentage of the assembly located in each contig/scaffold length bin.
