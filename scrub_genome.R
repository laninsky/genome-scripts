contam <- as.matrix(read.table("univecblast.txt"))
contam <- contam[(which(as.numeric(contam[,11])<=0.001)),]

table_contams <- as.matrix(table(contam[,1]))
write.table(table_contams,"freq_of_contam_matches.txt")

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


