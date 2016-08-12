#Step 1: run the onelining.R script on your original genome fasta file
 
#Step 2: make the blast database. Change the file name following -in to whatever your original fasta file was called.
~/bin/ncbi-blast-2.2.30+/bin/makeblastdb -in a.lines.fasta -dbtype nucl

#Step 3: line by line blast. Make sure that linebyline.R is in the directory with your fasta file and tempout file.
noseqs=`wc -l tempout | awk '{print $1}'`
echo "seqid" "pident" "match_length" "mismatch" "gapopen" "qstart" "qend" "sstart" "send" "evalue" "bitscore" "scaffold_length" > blastoutput.txt
for i in `seq 1 2 $noseqs`;
do tail -n+$i tempout | head -n1 > tempseq;
j=`expr $i + 1`;
tail -n+$j tempout | head -n1 >> tempseq;
tail -n+$j tempout | head -n1 | awk '{print length}' > length.txt;
echo "qseqid""sqeqid" "pident" "match_length" "mismatch" "gapopen" "qstart" "qend" "sstart" "send" "evalue" "bitscore" > rawblast.txt;
~/bin/ncbi-blast-2.2.30+/bin/blastn -db a.lines.fasta -query tempseq -perc_identity 75 -outfmt 6 >> rawblast.txt;
Rscript linebyline.R;
done;
