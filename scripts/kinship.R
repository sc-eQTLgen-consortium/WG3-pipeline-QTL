args = commandArgs(trailingOnly=TRUE)

n = read.delim(args[1],header=F)

kinIds = read.delim(args[2])

kin=kin*2

rownames(kin)=colnames(kin)=kinIds[,2]

write.table(kin,args[3],sep="\t",quote=F,col.names=NA)

##End R
