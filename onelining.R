#Unlike other versions of this script, writes out as it goes so it doesn't end up being such a memory hog. Change input name to whatever your file is called.

intable <- read.table("inputname",header=FALSE,stringsAsFactors=FALSE,sep="\t")
# e.g. intable <- read.table("a.lines.fasta",header=FALSE,stringsAsFactors=FALSE,sep="\t")
# e.g. intable <- read.table("a.lines.fasta",header=FALSE,stringsAsFactors=FALSE,sep="\t")

rows <- dim(intable)[1]

to_write <- intable[1,1]
write.table(to_write, "tempout",quote=FALSE, col.names=FALSE,row.names=FALSE)

sequencepaste <- NULL

for (j in 2:rows) {
if ((length(grep(">",intable[j,1])))>0) {
to_write <- toupper(sequencepaste)
to_write <- rbind(to_write,intable[j,1])
write.table(to_write, "tempout",quote=FALSE, col.names=FALSE,row.names=FALSE,append=TRUE)
sequencepaste <- NULL
} else {
sequencepaste <- paste(sequencepaste,intable[j,1],sep="")
}
}

write.table(toupper(sequencepaste), "tempout",quote=FALSE, col.names=FALSE,row.names=FALSE,append=TRUE)

q()
