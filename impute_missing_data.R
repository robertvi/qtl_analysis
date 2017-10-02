#!/usr/bin/Rscript

#
# impute all missing data using k nearest neighbour
#

args <- commandArgs(trailingOnly = TRUE)

if(length(args) != 7)
{
    cat("usage: impute_missing_data.R <input_file> <output_file>\n\n")
    cat("input_file/output_file: plink --recode A type .raw file, space separated columns\n")
    cat("headings: FID IID PAT MAT SEX PHENOTYPE <marker1> [<marker2>...]\n")
    quit()
}

inpfile   <- args[1]
outfile   <- args[2]

library(DMwR)

#load from plink --recode A .raw file
#columns must be: FID IID PAT MAT SEX PHENOTYPE marker1 marker2...
data <- read.csv(inpfile,header=T,sep='')

n_markers <- ncol(data) - 6
n_samples <- nrow(data)

cat("loaded markers:",n_markers," samples:",n_samples,"\n")

#extract just the genotypes as a numerical matrix
mat_x <- data.matrix(data[,7:ncol(data)])

#impute missing genotype values
df <- as.data.frame(t(mat_x))
df_imp <- knnImputation(df,k=1)

#reinsert into original dataframe
data[,7:ncol(data)] <- t(df_imp)

#save just the estimated effects information to csv file
write.table(data,file=outfile,quote=F,row.names=F,sep=' ')
