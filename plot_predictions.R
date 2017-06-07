#!/usr/bin/Rscript

#
# plot multiple predictions versus true phenotype values
#

library(ggplot2)
#source("~/rjv_mnt/cluster/git_repos/qtl_analysis/qtl_funcs.R")

args <- commandArgs(trailingOnly = TRUE)

if(length(args) < 5)
{
    cat("usage: plot_predictions.R <true_phenotypes_file> <output_image_file> [<predictions_file> <phenotype_label> <prediction_label>]... \n\n")
    cat("true_phenotypes_file: csv with header: sample, <phenotype1_name> [<phenotype2_name>...]\n")
    cat("output_image_file: PNG/JPG file to output plot to\n")
    cat("predictions_file: two column csv with header: sample, prediction\n")
    cat("phenotype_label: which phenotype to compare against - must match a column label from true_phenotypes_file\n")
    cat("prediction_label: plot label for this prediction method\n")

    quit()
}

phenofile <- args[1]
outfile <- args[2]

#phenofile <- "phenotypes.csv"
#outfile <- "test_pred_plot.png"

#load genotypes from csv file
## column headers required: marker,type,<sample1name>,<sample2name>...
##marker=marker name string
##type=any string
##sample names are strings
##then one marker per line, genotype codes as integers
phenotypes <- read.csv(phenofile,header=T)

#load qtl data
i <- 3
df <- NULL
while(i < length(args))
{
    #get next three options from command line
    predfile <- args[i]
    phenolabel <- args[i+1]
    predlabel <- args[i+2]
    i <- i + 3

    cat(predfile,"\n")

    #load predictions
    newpred <- read.csv(predfile,header=T)

    newdf <- data.frame(sample=phenotypes$sample,true_val=phenotypes[phenolabel])

    newdf <- merge(newdf,newpred,by.x="sample",by.y="sample")
    colnames(newdf) <- c("sample","true_val","pred_val")
    newdf <- newdf[!is.na(newdf$true_val),]

    newdf$phenotype = phenolabel
    newdf$method = predlabel

    df <- rbind(df,newdf)
}

gg <- ggplot(df,aes(x=true_val,y=pred_val)) + geom_point(aes(colour=method),alpha=0.5) + facet_wrap(~phenotype,scales="free")
ggsave(file=outfile,gg,width=8,height=6,dpi=300)
