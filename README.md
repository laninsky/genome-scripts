# genome-scripts
Some hodge podge scripts for dealing with genome stuff

#Transcriptomes
This folder contains instructions on the KU BI transcriptome pipeline

#Other stuff

filter_for_length_fastx.R filters scaffolds/contigs below a certain length, and also reports total assembly size/number of contigs, N10-N90/L10-L90. It requires a fastx file as input, which you can generate using fastx-tool's fasta_formatter with the -t flag. Paste the R-script into R and then call it by filter_for_length_fastx(your_fastx_file,contig_length).

filter_contig_by_length.R does the same thing but using a fasta file as input. It is MUCH slower than filter_for_length_fastx.R, so is not recommended. 

onelining.R removes the line breaks within the sequence associated with a single contig/scaffold (i.e. gets rid of the "hard wrap" that a lot of the default programs place in fasta files). The resulting fasta file should be twice as long as the number of contigs (one line for the header row, one line for the sequence).

linebylineblast.sh (and the associated R-scripts: onelining_tempseq.R and linebyline.R) blast an assembly against itself to look for regions of an assembly that made have been subject to a "false tandem duplication" (e.g. instead of alleles being assembled together, they have been assembled end to end).

summarizing_blast_hits.R summarizes the results of linebylineblast.sh

length_dist.R can be used to get the lengths of the contigs/scaffolds in the assembly and to generate histogram summaries of this data, as well as the percentage of the assembly located in each contig/scaffold length bin.

vectorcontam.sh (and the associated R-script: scrub_genome.R): searches your genome for contamination, and scrubs this from the genome, outputting contigs/scaffolds greater than 100 bp in size.

restart_vector_contam.sh (and the associated R-script: restart_scrub_genome.R): how to restart vectorcontam.sh if the job crashes part-way through.
