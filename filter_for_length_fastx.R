#file_to_filter <- "robin_genome_28May2018.fastx"
# temp <- temp[600000:600010,]

filter_for_length_fastx <- function(file_to_filter,contig_length) {

if(missing(file_to_filter) | missing(contig_length)) {
  print(paste("This script needs you to define the location and name of your fastx file and the minimum contig length you want to retain in the outfile: length_filtered_genome.fasta"))
  print(paste("Example of calling filter_contig_by_length:"))
  cat('filter_for_length_fastx("robin_genome_28May2018.fastx",2000)\n\n')
  stop("Please fill in the missing info in your function call")
}

library(dplyr)
library(readr)
library(stringr)
  
#Reading in the fastx file  
temp <- read_delim(file_to_filter,col_names=FALSE,delim="\t")

#Creating a column with the sequence length and number of Ns for each contig
temp <- temp %>% mutate(seq_length = str_length(X2))
temp <- temp %>% mutate(no_of_Ns = str_length(gsub("N","",X2,fixed=TRUE)))

  
  
  
  
temp <- mutate(temp,nchar(temp[,2]))
temp <- mutate(temp, str_length(gsub("N","",X2,fixed=TRUE)))  

names(temp) <- c("scaffold_name","sequence","length","no_of_Ns")  

#Writing out the length of each scaffold/contig
write.table(temp[,c(1,3)],  
  
#Setting the variable to record sequence (we are going to write out all the sequence on one line for each scaffold, rather than having linebreaks in teh sequence)
sequencerec <- NULL

con <- file(file_to_filter)
open(con)

while ( TRUE ) {
    line = readLines(con, n = 1)
    #If the line is longer than length 0 (i.e. empty carriage returns at the end of the file)
    if ( length(line) == 0 ) {
      break
    }
    #If there is a ">" in the line, and we've previously recorded sequence (i.e. we won't do this for the first line in the file)
    if (grepl(">",line,fixed=TRUE) && (!(is.null(sequencerec)))) {
      #Getting the length of the sequence
      sequencelength <- nchar(sequencerec)
      #If the sequence is longer than our defined length, finding the non-gapped length, and writing out the scaffold name, length and non-gapped length to "name_lengths_Ns.txt". Writing out the scaffold name and sequence to "scrubbed_genome.fasta"
      if(sequencelength>=contig_length) {
        sequencelengthN <- nchar(gsub("N","",sequencerec,fixed=TRUE))
        lengthdist <- t(as.matrix(c(seqname,sequencelength,sequencelengthN)))
        write.table(lengthdist,"name_lengths_Ns_filtered.txt",append=TRUE,quote=FALSE,row.names=FALSE,col.names=FALSE)
        write.table(seqname,"length_filtered_genome.fasta",append=TRUE,quote=FALSE,row.names=FALSE,col.names=FALSE)
        write.table(sequencerec,"length_filtered_genome.fasta",append=TRUE,quote=FALSE,row.names=FALSE,col.names=FALSE)
      }
      #resetting the sequencerec value and changing the sequencename to the line we just read in
      sequencerec <- NULL
      seqname <- line
    } else {
        #or else, if we didn't grep ">", appending sequence to the previous sequence
        if (!(grepl(">",line,fixed=TRUE))) {
          sequencerec <- paste(sequencerec,line,sep="")
        } else {
        #This condition (sequence is null, we've grepped ">"), should only be true for the first line  
          seqname <- line
        }
    }
}        

# The following code writes out the last sequence in the file if it is over 100bp
sequencelength <- nchar(sequencerec)
  if(sequencelength>=contig_length) {
    sequencelengthN <- nchar(gsub("N","",sequencerec,fixed=TRUE))
    lengthdist <- t(as.matrix(c(seqname,sequencelength,sequencelengthN)))
        write.table(lengthdist,"name_lengths_Ns_filtered.txt",append=TRUE,quote=FALSE,row.names=FALSE,col.names=FALSE)
        write.table(seqname,"length_filtered_genome.fasta",append=TRUE,quote=FALSE,row.names=FALSE,col.names=FALSE)
        write.table(sequencerec,"length_filtered_genome.fasta",append=TRUE,quote=FALSE,row.names=FALSE,col.names=FALSE)
  }

#Closing the connection to the file
close(con)
}
