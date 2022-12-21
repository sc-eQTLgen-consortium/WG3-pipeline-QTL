#Rscript to process the output of the run.

args = commandArgs(trailingOnly=TRUE)

mafFilter <- 0.01
machR2Filter <- 0.6
callFilter <- 0.99

##Read information
variants = data.table::fread(args[1])

##Filter on MAF
variants = variants[which(variants$MAF>=mafFilter),]

##Filter on imputation score
variants = variants[which(variants$MACH_R2>=machR2Filter),]

##Filter on call rate
variants = variants[which(variants$CALL>=callFilter),]

#To upload to central site.
write.table(variants, args[2],quote=F,sep="\t",row.names=F)

for(chr in 1:22){
	#To filter in genotype harmonizer
	write.table(variants$ID[which(variants$CHR==chr)], paste0(args[3],chr,".vars"),quote=F,sep="\t",row.names=F)

}

##End R

