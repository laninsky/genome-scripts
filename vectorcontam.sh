PBS script for running VecScreen-like search of UniVec database (univec.fasta) against hirise assembly (trunk_anole_19Jun2016.xkeD9.fasta)

```
#PBS -N univecblast
#PBS -l nodes=1:ppn=1:avx,mem=2000m,walltime=24:00:00
#PBS -M a499a400@ku.edu
#PBS -r n
#PBS -m n
#PBS -j oe
#PBS -o /scratch/a499a400/anolis/dovetail/blasterror.log
#PBS -d /scratch/a499a400/anolis/dovetail

blastn -task blastn -reward 1 -penalty -5 -gapopen 3 -gapextend 3 -dust yes -soft_masking true -evalue 700 -searchsp 1750000000000 -db trunk_anole_19Jun2016_xkeD9.fasta -query univec.fasta -outfmt 6 >> univecblast.txt
```

Interactive code for scrubbing portions of scaffolds/contigs identified as vector contamination, using an e-value threshold of 0.001.
```
#Define the name of your final scaffolded genome assembly
genome=trunk_anole_19Jun2016_xkeD9.fasta
echo $genome > temp

#Rscript for removing portions of scaffolds 
Rscript scrub_genome.R
