#!/usr/bin/Rscript

#
# calculate Kruskal Wallis statistics for all markers
# output to a csv file
#

args <- commandArgs(trailingOnly = TRUE)

if(length(args) != 4)
{
    cat("usage: calc_kruskalwallis.R <genotype_file> <phenotype_file> <phenotype> <output_file>\n\n")
    cat("genotype file csv column headings: marker,type,<sample1name>,<sample2name>\n")
    cat("phenotype file csv column headings: sample,<phenotype1name>[,<phenotype2name>...]\n")
    cat("phenotype: which phenotype columns to process (must match a column heading in the phenotype file)\n")
    cat("output_file: name for output file\n")
    quit()
}


#source("~/git_repos/qtl_analysis/qtl_funcs.R")
source("~/rjv_mnt/cluster/git_repos/qtl_analysis/qtl_funcs.R")

genofile <- args[1]
phenofile <- args[2]
pheno_col <- args[3]
outfile <- args[4]

#load genotypes from csv file
## column headers required: marker,type,<sample1name>,<sample2name>...
##marker=marker name string
##type=any string
##sample names are strings, must be unique and match exactly those in the phenotype file
##then one marker per line, genotype codes as integers
genotypes <- read.csv(genofile,header=T)

#load phenotypes
##column headers required: sample,<phenotype1name>[,<phenotype2name>...]
phenotypes <- read.csv(phenofile,header=T)

#match up genotype with phenotype by progeny name
#exclude progeny with missing phenotype values
res <- matchup_and_extract(genotypes,phenotypes,pheno_col)

#for each marker
df <- NULL
for( marker in colnames(res$comb)[2:ncol(res$comb)] )
{
    kw <- kruskal.test(res$vec_y,as.factor(res$mat_x[,marker]))
    bycall <- aggregate(res$vec_y,list(res$mat_x[,marker]),mean)
    colnames(bycall) <- c("call","mean")
    bycall <- bycall[order(bycall$call),,drop=F]

    newrow <- data.frame(marker,kw$parameter,kw$p.value,kw$statistic,paste(bycall$call,collapse="/"),paste(bycall$mean,collapse="/"))
    colnames(newrow) <- c("marker","dof","pvalue","kvalue","calls","means")

    df <- rbind(df,newrow)
}

write.csv(df,outfile,quote=FALSE,row.names=FALSE)
