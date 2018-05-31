#file_to_filter <- "robin_genome_28May2018.fastx"
# temp <- temp[600000:600010,]
# contig_length <- 201

filter_for_length_fastx <- function(file_to_filter,contig_length) {

if(missing(file_to_filter) | missing(contig_length)) {
  print(paste("This script needs you to define the location and name of your fastx file and the minimum contig length you want to retain in the outfile: length_filtered_genome.fasta"))
  print(paste("Example of calling filter_contig_by_length:"))
  cat('filter_for_length_fastx("robin_genome_28May2018.fastx",201)\n\n')
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
names(temp) <- c("scaffold_name","sequence","length","length_without_Ns")
temp <- temp %>% arrange(desc(length))  
  
#Writing out the length of each scaffold/contig and number of Ns
write.table(temp[,c(1,3,4)],"name_lengths_Ns_filtered.txt",quote=FALSE,row.names=FALSE)

#Writing out the N50,L50 and total assembly size etc for the input assembly
assemblylength <- sum(temp[,3])
print(paste("The input assembly is ",assemblylength," bp in length",sep=""))
print(paste("The input assembly is ",sum(temp[,4])," bp in length when excluding Ns",sep=""))  
print(paste("The assembly consists of ",dim(temp)[1]," total scaffolds/contigs",sep=""))

sumseqlength <- 0
N80 <- "DO"
N70 <- "DO"
N60 <- "DO"
N50 <- "DO"
N40 <- "DO"
N30 <- "DO"
N20 <- "DO"
N10 <- "DO"   
for (j in 1:(dim(temp)[1])) {
  sumseqlength <- sumseqlength + temp[j,3]
  if(sumseqlength>(0.9*assemblylength)) {
    print(paste("The assembly N90 is ",temp[j,3]," bp",sep=""))
    print(paste("The assembly L90 is ",j," scaffolds/contigs",sep=""))
    break
  }  
  if(sumseqlength>(0.8*assemblylength) & N80=="DO") {
    print(paste("The assembly N80 is ",temp[j,3]," bp",sep=""))
    print(paste("The assembly L80 is ",j," scaffolds/contigs",sep=""))
    N80 <- "DONE"
  }
  if(sumseqlength>(0.7*assemblylength) & N70=="DO") {
    print(paste("The assembly N70 is ",temp[j,3]," bp",sep=""))
    print(paste("The assembly L70 is ",j," scaffolds/contigs",sep=""))
    N70 <- "DONE"
  }
  if(sumseqlength>(0.6*assemblylength) & N60=="DO") {
    print(paste("The assembly N60 is ",temp[j,3]," bp",sep=""))
    print(paste("The assembly L60 is ",j," scaffolds/contigs",sep=""))
    N60 <- "DONE"
  }  
  if(sumseqlength>(0.5*assemblylength) & N50=="DO") {
    print(paste("The assembly N50 is ",temp[j,3]," bp",sep=""))
    print(paste("The assembly L50 is ",j," scaffolds/contigs",sep=""))
    N50 <- "DONE"
  }
  if(sumseqlength>(0.4*assemblylength) & N40=="DO") {
    print(paste("The assembly N40 is ",temp[j,3]," bp",sep=""))
    print(paste("The assembly L40 is ",j," scaffolds/contigs",sep=""))
    N40 <- "DONE"
  }
  if(sumseqlength>(0.3*assemblylength) & N30=="DO") {
    print(paste("The assembly N30 is ",temp[j,3]," bp",sep=""))
    print(paste("The assembly L30 is ",j," scaffolds/contigs",sep=""))
    N30 <- "DONE"
  }
  if(sumseqlength>(0.2*assemblylength) & N20=="DO") {
    print(paste("The assembly N20 is ",temp[j,3]," bp",sep=""))
    print(paste("The assembly L20 is ",j," scaffolds/contigs",sep=""))
    N20 <- "DONE"
  }
  if(sumseqlength>(0.1*assemblylength) & N10=="DO") {
    print(paste("The assembly N10 is ",temp[j,3]," bp",sep=""))
    print(paste("The assembly L10 is ",j," scaffolds/contigs",sep=""))
    N10 <- "DONE"
  }
}  
    
#Filtering for scaffold/contigs above the minimum size
temp <- temp %>% filter(.,length>=contig_length)   
  
#Writing out the N50,L50 and total assembly size etc for the output assembly
print(paste("After filtering out contigs shorter than ",contig_length," bp, the assembly had the following statistics:",sep=""))
assemblylength <- sum(temp[,3])
print(paste("The input assembly is ",assemblylength," bp in length",sep=""))
print(paste("The input assembly is ",sum(temp[,4])," bp in length when excluding Ns",sep=""))  
print(paste("The assembly consists of ",dim(temp)[1]," total scaffolds/contigs",sep=""))

sumseqlength <- 0
N80 <- "DO"
N70 <- "DO"
N60 <- "DO"
N50 <- "DO"
N40 <- "DO"
N30 <- "DO"
N20 <- "DO"
N10 <- "DO"   
for (j in 1:(dim(temp)[1])) {
  sumseqlength <- sumseqlength + temp[j,3]
  if(sumseqlength>(0.9*assemblylength)) {
    print(paste("The assembly N90 is ",temp[j,3]," bp",sep=""))
    print(paste("The assembly L90 is ",j," scaffolds/contigs",sep=""))
    break
  }  
  if(sumseqlength>(0.8*assemblylength) & N80=="DO") {
    print(paste("The assembly N80 is ",temp[j,3]," bp",sep=""))
    print(paste("The assembly L80 is ",j," scaffolds/contigs",sep=""))
    N80 <- "DONE"
  }
  if(sumseqlength>(0.7*assemblylength) & N70=="DO") {
    print(paste("The assembly N70 is ",temp[j,3]," bp",sep=""))
    print(paste("The assembly L70 is ",j," scaffolds/contigs",sep=""))
    N70 <- "DONE"
  }
  if(sumseqlength>(0.6*assemblylength) & N60=="DO") {
    print(paste("The assembly N60 is ",temp[j,3]," bp",sep=""))
    print(paste("The assembly L60 is ",j," scaffolds/contigs",sep=""))
    N60 <- "DONE"
  }  
  if(sumseqlength>(0.5*assemblylength) & N50=="DO") {
    print(paste("The assembly N50 is ",temp[j,3]," bp",sep=""))
    print(paste("The assembly L50 is ",j," scaffolds/contigs",sep=""))
    N50 <- "DONE"
  }
  if(sumseqlength>(0.4*assemblylength) & N40=="DO") {
    print(paste("The assembly N40 is ",temp[j,3]," bp",sep=""))
    print(paste("The assembly L40 is ",j," scaffolds/contigs",sep=""))
    N40 <- "DONE"
  }
  if(sumseqlength>(0.3*assemblylength) & N30=="DO") {
    print(paste("The assembly N30 is ",temp[j,3]," bp",sep=""))
    print(paste("The assembly L30 is ",j," scaffolds/contigs",sep=""))
    N30 <- "DONE"
  }
  if(sumseqlength>(0.2*assemblylength) & N20=="DO") {
    print(paste("The assembly N20 is ",temp[j,3]," bp",sep=""))
    print(paste("The assembly L20 is ",j," scaffolds/contigs",sep=""))
    N20 <- "DONE"
  }
  if(sumseqlength>(0.1*assemblylength) & N10=="DO") {
    print(paste("The assembly N10 is ",temp[j,3]," bp",sep=""))
    print(paste("The assembly L10 is ",j," scaffolds/contigs",sep=""))
    N10 <- "DONE"
  }
}  

#Populating output matrix  
output <- matrix(NA,ncol=1,nrow=(2*(dim(temp)[1])))
output[(seq(1,(dim(output)[1]),2)),1] <- paste(">",as.matrix(temp[(1:dim(temp)[1]),1]),sep="")
output[(seq(2,(dim(output)[1]),2)),1] <- as.matrix(temp[(1:dim(temp)[1]),2])  
  

#Writing this out
write.table(output,paste(file_to_filter,".",contig_length,"bp_filtered.fa",sep=""),quote=FALSE,col.names=FALSE,row.names=FALSE)
  
