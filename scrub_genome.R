contam <- as.matrix(read.table("univecblast.txt"))
contam <- contam[(which(as.numeric(contam[,11])<=0.001)),]

table_contams <- as.matrix(table(contam[,1]))
write.table(table_contams,"freq_of_contam_matches.txt")
write.table(contam,"identity_of_contam_matches.txt")

scaffoldnames <- names(table(contam[,2]))

range_to_delete <- matrix(nrow=(length(scaffoldnames)),ncol=2)
range_to_delete[,1] <- scaffoldnames

for (i in scaffoldnames) {
  to_delete <- NULL
  subscaff <- contam[(which(contam[,2]==i)),9:10]

  if (!(is.matrix(subscaff))) {  
        if(as.numeric(subscaff[1])>as.numeric(subscaff[2])) {
          tempsubscaff <- subscaff[1]
          subscaff[1] <- subscaff[2]
          subscaff[2] <- tempsubscaff
        }
    to_delete <- paste(as.numeric(subscaff[1]),":",as.numeric(subscaff[2]),sep="")
  
  } else {
   
    subscaff <- cbind(as.numeric(subscaff[,1]),as.numeric(subscaff[,2]))
    
    for (j in 1:(dim(subscaff)[1])) {
      if(subscaff[j,1]>subscaff[j,2]) {
        tempsubscaff <- subscaff[j,1]
        subscaff[j,1] <- subscaff[j,2]
        subscaff[j,2] <- tempsubscaff
      }
    }

    subscaff <- subscaff[order(subscaff[,2]),]
    subscaff <- subscaff[order(subscaff[,1]),]
    
    tempscaff <- as.matrix(subscaff[1,],nrow=2)
    
    for (j in 2:(dim(subscaff)[1])) {
      if((subscaff[j,1]>=tempscaff[1,(dim(tempscaff)[2])]) && (subscaff[j,1]<tempscaff[2,(dim(tempscaff)[2])]) && (subscaff[j,2]>tempscaff[1,(dim(tempscaff)[2])]) && (subscaff[j,2]<=tempscaff[2,(dim(tempscaff)[2])])) {
        break
      } else {
        if((subscaff[j,1]>=tempscaff[1,(dim(tempscaff)[2])]) && (subscaff[j,1]<=(tempscaff[2,(dim(tempscaff)[2])]+200))) {
          tempscaff[2,(dim(tempscaff)[2])] <- subscaff[j,2]
          break
        }
        temp <- rbind(subscaff[j,1],subscaff[j,2])
        tempscaff <- cbind(tempscaff,temp)
      }
    }
  
   to_delete <- paste(tempscaff[1,1],":",tempscaff[2,1],sep="")
   
   if(dim(tempscaff)[2]>1) { 
      for (j in 2:(dim(tempscaff)[2])) {
        to_delete <- paste(to_delete,paste(tempscaff[1,j],":",tempscaff[2,j],sep=""),sep=",")
      } 
   }
      
      }
range_to_delete[which(scaffoldnames==i),2] <- to_delete

} 

rm(contam)
rm(table_contams)
rm(temp)
rm(to_delete)
rm(scaffoldnames)
rm(tempscaff)
rm(subscaff)
rm(tempsubscaff)

lengthdist <- matrix(c("scaffoldname","length","non-N-length"),ncol=3)
write.table(lengthdist,"name_lengths_Ns.txt",append=FALSE,quote=FALSE,row.names=FALSE,col.names=FALSE)

sequencerec <- NULL

filename <- as.matrix(read.table("temp"))
con <- file(filename[1,1])
open(con)

  while ( TRUE ) {
    line = readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    if (grepl(">",line,fixed=TRUE) && (!(is.null(sequencerec)))) {
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
        lengthdist <- c(seqname,sequencelength,sequencelengthN)
        write.table(lengthdist,"name_lengths_Ns.txt",append=TRUE,quote=FALSE,row.names=FALSE,col.names=FALSE)
          
        

        
        
        
    
        
    seqname <- line    
    }



  
close(con)


https://stat.ethz.ch/pipermail/r-help/2003-April/033097.html

