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
numCol <- ncol(dataExp);

#Update with gene names instead of probes 
dataExp[,"MAX"] <- apply(dataExp, FUN=max, MARGIN=1);
dataExp <- dataExp[order(-dataExp[,"MAX"]),]
annot <- annot[rownames(dataExp),];
dataExp[,"Hugo_Gene_Symbol"] <- as.character(annot[,"Gene.Symbol"]);
dataExp[,"Entrez_Gene_Id"] <- as.character(annot[,"Entrez.Gene"]);

dataExp <- dataExp[!duplicated(dataExp[,"Hugo_Gene_Symbol"]),]
dataExp <- dataExp[!grepl("\\//", dataExp[,"Hugo_Gene_Symbol"]),];
rownames(dataExp) <- dataExp[,"Hugo_Gene_Symbol"];
dataAnnot <- dataExp[,c((numCol+2):ncol(dataExp))];
dataExp <- dataExp[,c(1:numCol)];


zLog <- function(x)
{
tmp <- log2(x);
tmp <- (tmp-mean(tmp))/sd(tmp);
}

zDataExp <- data.frame(t(apply(dataExp, FUN=zLog, MARGIN=1)));
dataExp <- data.frame(dataAnnot, dataExp);
zDataExp <- data.frame(dataAnnot, zDataExp);

colnames(dataExp) <- gsub(".CEL", "", colnames(dataExp));
colnames(zDataExp) <- gsub(".CEL", "", colnames(zDataExp));


setwd(curDir);
expFileName <- paste(args[3], "_Exprs.txt", sep="");
zExpFileName <- paste(args[3], "_Z_Scores.txt", sep="");

write.table(dataExp, expFileName, sep="\t", row.names=T);
write.table(zDataExp, zExpFileName, sep="\t", row.names=T);







