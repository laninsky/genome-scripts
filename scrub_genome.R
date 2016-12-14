contam <- as.matrix(read.table("univecblast.txt"))
contam <- contam[(which(as.numeric(contam[,11])<=0.001)),]

table_contams <- as.matrix(table(contam[,1]))
write.table(table_contams,"freq_of_contam_matches.txt")

scaffoldnames <- names(table(contam[,2]))

range_to_delete <- matrix(nrow=(length(scaffoldnames)),ncol=2)
range_to_delete[,1] <- scaffoldnames

for (i in scaffoldnames) {
subscaff <- contam[(which(contam[,2]==i)),9:10]

for (j in 1:(dim(subscaff)[1])) {
if(as.numeric(subscaff[j,1])>as.numeric(subscaff[j,2])) {
tempsubscaff <- subscaff[j,1]
subscaff[j,1] <- subscaff[j,2]
subscaff[j,2] <- tempsubscaff
}
}


for (j in 1:(dim(subscaff)[1])) {




range_to_delete[which(scaffoldnames==i),2] <- to_delete
