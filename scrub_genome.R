#Reading in the file with our vector contamination matches, and taking just the matches at an evalue of 0.001 or less
contam <- as.matrix(read.table("univecblast.txt"))
contam <- contam[(which(as.numeric(contam[,11])<=0.001)),]

#Writing out tables on the frequency that each vector was found (freq_of_contam_matches) and the specific matches between scaffolds and contamination (identity_of_contam_matches)
table_contams <- as.matrix(table(contam[,1]))
write.table(table_contams,"freq_of_contam_matches.txt")
write.table(contam,"identity_of_contam_matches.txt")

#Getting a list of the scaffold names with identified contamination and setting up a matrix to record the positions in the scaffold we need to delete
scaffoldnames <- names(table(contam[,2]))
range_to_delete <- matrix(nrow=(length(scaffoldnames)),ncol=2)
range_to_delete[,1] <- scaffoldnames

#A loop to record the ranges we need to delete for each scaffold
for (i in scaffoldnames) {
  to_delete <- NULL

  #Positions 9 and 10 have the start and end of the match to the subject
  subscaff <- contam[(which(contam[,2]==i)),9:10]

  #Where only one stretch of contamination is found, this if loop is used
  if (!(is.matrix(subscaff))) {  
        #For matches to the reverse strand, this swaps the beginning and end coordinates
        if(as.numeric(subscaff[1])>as.numeric(subscaff[2])) {
          tempsubscaff <- subscaff[1]
          subscaff[1] <- subscaff[2]
          subscaff[2] <- tempsubscaff
        }
    #Recording beginning and end coordinates in a temporary vector
    to_delete <- paste(as.numeric(subscaff[1]),":",as.numeric(subscaff[2]),sep="")
  
  #If there is more than one stretch of contamination, this loop is used
  } else {
    #Converting the matrix to numerical values
    subscaff <- cbind(as.numeric(subscaff[,1]),as.numeric(subscaff[,2]))
    
    #Reversing the order of the beginning and end coordinates where matches to reverse strand occured
    for (j in 1:(dim(subscaff)[1])) {
      if(subscaff[j,1]>subscaff[j,2]) {
        tempsubscaff <- subscaff[j,1]
        subscaff[j,1] <- subscaff[j,2]
        subscaff[j,2] <- tempsubscaff
      }
    }
    #Ordering the matrix of matches from beginning of scaffold to end of scaffold, sorted first on beginning and then on end coordinates
    subscaff <- subscaff[order(subscaff[,2]),]
    subscaff <- subscaff[order(subscaff[,1]),]
    
    #Taking the first coordinates
    tempscaff <- as.matrix(subscaff[1,],nrow=2)
    
    #Looping through the remaining ranges. Where these ranges overlap, or are within 200 bp of each other, extending the sequence to move out over the overlapping segments
    for (j in 2:(dim(subscaff)[1])) {
      if((subscaff[j,1]>=tempscaff[1,(dim(tempscaff)[2])]) && (subscaff[j,1]<tempscaff[2,(dim(tempscaff)[2])]) && (subscaff[j,2]>tempscaff[1,(dim(tempscaff)[2])]) && (subscaff[j,2]<=tempscaff[2,(dim(tempscaff)[2])])) {
        next
      } else {
        if((subscaff[j,1]>=tempscaff[1,(dim(tempscaff)[2])]) && (subscaff[j,1]<=(tempscaff[2,(dim(tempscaff)[2])]+200))) {
          tempscaff[2,(dim(tempscaff)[2])] <- subscaff[j,2]
          next
        }
        temp <- rbind(subscaff[j,1],subscaff[j,2])
        tempscaff <- cbind(tempscaff,temp)
        print("binding!")
      }
    }
  
   to_delete <- paste(tempscaff[1,1],":",tempscaff[2,1],sep="")
   
   #If, after combining overlapping (or nearby) ranges, there are still more than one stretches of separate contamination identified, mushing these into a single string to write out
   if(dim(tempscaff)[2]>1) { 
      for (j in 2:(dim(tempscaff)[2])) {
        to_delete <- paste(to_delete,paste(tempscaff[1,j],":",tempscaff[2,j],sep=""),sep=",")
      } 
   }
      
      }
#Adding the stretches to delete to the appropriate row  
range_to_delete[which(scaffoldnames==i),2] <- to_delete

} 

#Removing variables not needed for the second part of this script
rm(contam)
rm(table_contams)
rm(temp)
rm(to_delete)
rm(scaffoldnames)
rm(tempscaff)
rm(subscaff)
rm(tempsubscaff)

#Setting up a table where we will record scaffolds, lengths and lengths excluding gapped positions (non-N-length)
lengthdist <- matrix(c("scaffoldname","length","non-N-length"),ncol=3)
write.table(lengthdist,"name_lengths_Ns.txt",append=FALSE,quote=FALSE,row.names=FALSE,col.names=FALSE)

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

