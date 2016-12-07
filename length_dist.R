#After previously converting genome files to not have line breaks
temp <- as.matrix(read.table("tempout"))
backuptemp <- temp #just in case I muffed the original because importing this file takes forever
linestoremove <- which(temp[,1]=="circular")
temp <- as.matrix(temp[-(linestoremove),1],ncol=1)
lengths <- nchar(temp[(seq(2,(dim(temp)[1]),2)),1])

# Discovar outputs "circular" as an additional header for some of the contigs, and potentially there are other
names <- (temp[(seq(1,(dim(temp)[1]),2)),1])
test <- grep(">",names[1:length(names)])
if (length(test)==length(names)) {
print("All good. Each title row has a '>'") } else {
print("Somewhere in your file there are sequences which have multiple header rows per sequence")
print("Find these and remove them before proceeding")
}

newmatrix <- as.matrix(cbind(names,lengths),ncol=1)
write.table(newmatrix,"name_lengths.txt",quote=FALSE,row.names=FALSE)

#To plot the distribution of contig lengths, getting a log scale appropriate to our data
upperno <- round(max(lengths),digits=-((nchar(max(lengths))-1)))

histbreaks <- NULL
incr <- 10
x <- 0
while (x < upperno) {
x <- x+incr
histbreaks <- c(histbreaks,x)

if (length(histbreaks)>9) {
if (((histbreaks[length(histbreaks)])/(histbreaks[(length(histbreaks)-9)])==10)) {
while (x < upperno) {
temphistbreaks <- histbreaks[(length(histbreaks)-8):(length(histbreaks))]*incr
histbreaks <- c(histbreaks,temphistbreaks)
x <- histbreaks[length(histbreaks)]
}
}
}
}

histbreaks <- c(0,histbreaks)

#Calculating the number of contigs/scaffolds in each size category
fragdist <- NULL
for (i in 2:length(histbreaks)) {
bin_no <- length(which(lengths<=histbreaks[i]))-length(which(lengths<=histbreaks[i-1]))
fragdist <- c(fragdist,bin_no)
}

fragdist <- c("NA",fragdist)

#Writes out a table with left column = the frag size cut offs, and right hand = the number of fragments less than or equal to the bin size, but greater than the previous bin size
newmatrix <- as.matrix(cbind(histbreaks,fragdist),ncol=1)
write.table(newmatrix,"hist_lengths.txt",quote=FALSE,row.names=FALSE)


#Slightly different calcualtion: calculating the proportion of bases in the total assembly found in each size category
totallength <- sum(as.numeric(lengths))

fragdist <- NULL
for (i in 2:length(histbreaks)) {
bin_no <- (sum(as.numeric(lengths[(which(lengths<=histbreaks[i]))]))-sum(as.numeric(lengths[(which(lengths<=histbreaks[i-1]))])))/totallength
fragdist <- c(fragdist,bin_no)
}

fragdist <- c("NA",fragdist)

#Writes out a table with left column = the frag size cut offs, and right hand = the % of the assembly contained  in fragments less than or equal to the bin size, but greater than the previous bin size
newmatrix <- as.matrix(cbind(histbreaks,fragdist),ncol=1)
write.table(newmatrix,"hist_perc_assembly.txt",quote=FALSE,row.names=FALSE)
