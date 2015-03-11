# setwd("H:\\RNA-seq")

library("DESeq")
datatype <- "file_path"	
Nrep1 <- control_num						# number of replicates for control
Nrep2 <- treatment_num						# number of replicates for treatment

# read in data information #	
datafile <- paste(datatype, ".txt", sep="")
countTable <- read.table(datafile, header=TRUE, row.names=1)
condition = c(rep("control", times=Nrep1), rep("treatment", times=Nrep2))

cds = newCountDataSet(countTable,condition)

# perfor normalization #
cds = estimateSizeFactors(cds)
sizeFactors(cds)
head(counts(cds,normalized=TRUE))

# estimate variance #
if(Nrep1==1 || Nrep2==1) cds=estimateDispersions(cds,method='blind',sharingMode="fit-only") else cds=estimateDispersions(cds)
plotDispEsts(cds)

# call differential expression
res = nbinomTest(cds,"control","treatment")
plotMA(res)

# output differential expression result
outputfile <- paste(datatype, "_DEseq_DE_output.txt", sep="")
write.table(res, outputfile, sep = '\t')