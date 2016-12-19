#What happens if my job dies?
#Find the last sequence to be written out to the scrubbed_genome.fasta, and make a note of its name. Just in case this sequence is incomplete, we are going to remove it from the scrubbed_genome.fasta file and restart the code at this point. If just the sequence name is written out:

tail -n -1 scrubbed_genome.fasta > test
less test #To make sure everything looks good before you over-write scrubbed_genome.fasta
mv test scrubbed_genome.fasta

#If the sequene name AND some associated sequence has been written out:
tail -n -2 scrubbed_genome.fasta > test
less test #To make sure everything looks good before you over-write scrubbed_genome.fasta
mv test scrubbed_genome.fasta

#If this last sequence is also present in the name_lengths_Ns.txt file:
tail -n -1 name_lengths_Ns.txt > test
less test #To make sure everything looks good before you over-write name_lengths_Ns.txt
mv test name_lengths_Ns.txt

#We then need to modify our input genome file. Search for the name of the last sequence (that we have now deleted from the scrubbed_genome.fasta file), and find the line number of this in the original genome assembly file e.g. 49976687. In the code below, replace trunk_anole_19Jun2016_xkeD9.fasta with the name of your assembly, and the number on lastlines with the line number of the last sequence in your original genome assembly.
totalline=`wc -l trunk_anole_19Jun2016_xkeD9.fasta | awk '{print $1}'`
lastlines=$(($totalline-49976687+1))

tail -n $lastlines trunk_anole_19Jun2016_xkeD9.fasta > test
less test #To make sure everything looks good before you move to mod.fasta
mv test mod.fasta

#Define the name of your modified final scaffolded genome assembly
genome=mod.fasta
echo $genome > temp

#Making a backup just in case anything goes wrong
cp scrubbed_genome.fasta backup_scrubbed_genome.fasta
cp name_lengths_Ns.txt backup_name_lengths_Ns.txt

#Rscript for removing portions of scaffolds 
Rscript restart_scrub_genome.R
