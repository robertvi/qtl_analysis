#!/usr/bin/Rscript

#
# add a second numerical encoding for each marker set to 1 if the genotype is a het
# otherwise zero
#

args <- commandArgs(trailingOnly = TRUE)

if(length(args) != 2)
{
    cat("usage: add_dominance_info.R <input_file> <output_file>\n\n")
    cat("input_file/output_file: plink --recode A type .raw file, space separated columns\n")
    cat("headings: FID IID PAT MAT SEX PHENOTYPE <marker1> [<marker2>...]\n")
    quit()
}

inpfile   <- args[1]
outfile   <- args[2]

#load from plink --recode A .raw file
#columns must be: FID IID PAT MAT SEX PHENOTYPE marker1 marker2...
data <- read.csv(inpfile,header=T,sep='')

n_markers <- ncol(data) - 6
n_samples <- nrow(data)

cat("loaded markers:",n_markers," samples:",n_samples,"\n")

marker_names <- colnames(data)[7:ncol(data)]

for ( name in marker_names )
{
    new_name <- paste0(unlist(strsplit(name,"_"))[1],"_het")
    data[new_name] <- 0
    wh <- data[name]==1
    data[wh,new_name] <- 1
}

#save imputed data
write.table(data,file=outfile,quote=F,row.names=F,sep=' ')
