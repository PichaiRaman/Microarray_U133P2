###########################################
#Script to process U133 P2 microarray Data
#Pichai Raman
#March 7th 2016
#Just call this function Rscript process_U133P2.R CEL FILE DIRECTORY
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

#Process
tmpData <- ReadAffy();
tmpData <- mas5(tmpData, sc=150);
dataExp <- exprs(tmpData);

zLog <- function(data)
{
tmp <- log2(data);
tmp <- (tmp-mean(tmp))/sd(tmp);
}

zDataExp <- zLog(data);
setwd(curDir);
expFileName <- paste(args[1], "_Exprs.txt", sep="");
zExpFileName <- paste(args[1], "_Z_Scores.txt", sep="");

write.table(dataExp, expFileName, sep="\t", row.names=T);
write.table(dataExp, zExpFileName, sep="\t", row.names=T);







