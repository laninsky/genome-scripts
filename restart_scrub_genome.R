#Setting the variable to record sequence (we are going to write out all the sequence on one line for each scaffold, rather than having linebreaks in teh sequence)
sequencerec <- NULL

#Getting your genome file from the temp file set up in the bash script (echo $genome > temp)
filename <- as.matrix(read.table("temp"))

#Creating a connection to the file so we can read in line at a time (to reduce insane memory requirements)
con <- file(filename[1,1])
open(con)

#while there are lines to read...
while ( TRUE ) {
    line = readLines(con, n = 1)
    #If the line is longer than length 0 (i.e. empty carriage returns at the end of the file)
    if ( length(line) == 0 ) {
      break
    }
    #If there is a ">" in the line, and we've previously recorded sequence (i.e. we won't do this for the first line in the file)
    if (grepl(">",line,fixed=TRUE) && (!(is.null(sequencerec)))) {
        #If the previous contig needs to be scrubbed for contamination
        if (gsub(">","",seqname,fixed=TRUE) %in% range_to_delete[,1]) {
        
          #Splitting the DNA string into a vector so we can extract the contaminated positions
          sequencerec <- unlist(strsplit(sequencerec,split=""))
        
          #Obtaining the contaminated positions in an array from the range_to_delete matrix  
          toremove <- unlist(strsplit(range_to_delete[(which(range_to_delete[,1]==gsub(">","",seqname,fixed=TRUE))),2],split=","))
          toremovemod <- NULL  
          for (i in toremove) {
            temp <- as.numeric(unlist(strsplit(i,split=":")))
            toremovemod <- c(toremovemod,temp[1]:temp[2])
          }
        
        #Removing those positions and collapsing the sequence back into a single string
          sequencerec <- sequencerec[-toremovemod]
          sequencerec <- paste(sequencerec,collapse="")
        }
      #Getting the length of the sequence
      sequencelength <- nchar(sequencerec)
      
      #If the sequence is longer than 100, finding the non-gapped length, and writing out the scaffold name, length and non-gapped length to "name_lengths_Ns.txt". Writing out the scaffold name and sequence to "scrubbed_genome.fasta"
      if(sequencelength>=100) {
        sequencelengthN <- nchar(gsub("N","",sequencerec,fixed=TRUE))
        lengthdist <- t(as.matrix(c(seqname,sequencelength,sequencelengthN)))
        write.table(lengthdist,"name_lengths_Ns.txt",append=TRUE,quote=FALSE,row.names=FALSE,col.names=FALSE)
        write.table(seqname,"scrubbed_genome.fasta",append=TRUE,quote=FALSE,row.names=FALSE,col.names=FALSE)
        write.table(sequencerec,"scrubbed_genome.fasta",append=TRUE,quote=FALSE,row.names=FALSE,col.names=FALSE)
      }
      #resetting the sequencerec value ang changing the sequencename to the line we just read in
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
if (gsub(">","",seqname,fixed=TRUE) %in% range_to_delete[,1]) {
  sequencerec <- unlist(strsplit(sequencerec,split=""))
        
  toremove <- unlist(strsplit(range_to_delete[(which(range_to_delete[,1]==gsub(">","",seqname,fixed=TRUE))),2],split=","))
  toremovemod <- NULL  
  for (i in toremove) {
    temp <- as.numeric(unlist(strsplit(i,split=":")))
    toremovemod <- c(toremovemod,temp[1]:temp[2])
  }
        
  sequencerec <- sequencerec[-toremovemod]
  sequencerec <- paste(sequencerec,collapse="")
}

sequencelength <- nchar(sequencerec)
  if(sequencelength>=100) {
    sequencelengthN <- nchar(gsub("N","",sequencerec,fixed=TRUE))
    lengthdist <- t(as.matrix(c(seqname,sequencelength,sequencelengthN)))
    write.table(lengthdist,"name_lengths_Ns.txt",append=TRUE,quote=FALSE,row.names=FALSE,col.names=FALSE)
    write.table(seqname,"scrubbed_genome.fasta",append=TRUE,quote=FALSE,row.names=FALSE,col.names=FALSE)
    write.table(sequencerec,"scrubbed_genome.fasta",append=TRUE,quote=FALSE,row.names=FALSE,col.names=FALSE)
  }

write.table(range_to_delete,"ranges_deleted_from_scaffolds.txt",append=TRUE,quote=FALSE,row.names=FALSE,col.names=FALSE)

#Closing the connection to the file
close(con)

