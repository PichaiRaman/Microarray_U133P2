###########################################
#Script to process U133 P2 microarray Data
#Pichai Raman
#March 7th 2016
#Just call this function Rscript process_U133P2.R CEL FILE DIRECTORY
#
#First Argument - Cel File directory
#Next Argument - Annotation File
#Last Argumnet - prefix of name
#
###########################################

#This is to pass arguments to R
options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)


#Call libraries
library(affy);

curDir <- getwd(); 
celDir = args[1];
setwd(celDir);

#Get annotation file
annot <- read.delim(args[2]);
rownames(annot) <- annot[,1];

#Process
tmpData <- ReadAffy();
tmpData <- mas5(tmpData, sc=150);
dataExp <- data.frame(exprs(tmpData));

#Update with gene names instead of probes 
dataExp[,"MAX"] <- apply(dataExp, FUN=max, MARGIN=1);
dataExp <- dataExp[order(-dataExp[,"MAX"]),]
annot <- annot[rownames(dataExp),];
dataExp[,"Gene"] <- as.character(annot[,"Gene.Symbol"]);
dataExp <- dataExp[!duplicated(dataExp[,"Gene"]),]
dataExp <- dataExp[!grepl("\\//", dataExp[,"Gene"]),];
rownames(dataExp) <- dataExp[,"Gene"];


zLog <- function(data)
{
tmp <- log2(data);
tmp <- (tmp-mean(tmp))/sd(tmp);
}

zDataExp <- zLog(data);
setwd(curDir);
expFileName <- paste(args[3], "_Exprs.txt", sep="");
zExpFileName <- paste(args[3], "_Z_Scores.txt", sep="");

write.table(dataExp, expFileName, sep="\t", row.names=T);
write.table(dataExp, zExpFileName, sep="\t", row.names=T);







