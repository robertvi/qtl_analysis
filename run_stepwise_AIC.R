#!/usr/bin/Rscript

#
# estimate SNP effects using stepwise AIC
# uses already calculated Kruskal Wallis values
#

args <- commandArgs(trailingOnly = TRUE)

if(length(args) != 5)
{
    cat("usage: run_stepwise_AIC.R <genotype_file> <phenotype_file> <phenotype> <kw_file> <output_basename>\n\n")
    cat("genotype_file csv, numerical genotype codes, column headings: marker,type,<sample1name>,<sample2name>\n")
    cat("phenotype_file csv, numerical phenotype values, column headings: sample,<phenotype1name>[,<phenotype2name>...]\n")
    cat("phenotype: which phenotype column to process (must match a column heading in the phenotype file)\n")
    cat("kw_file: file containing output from calc_kruskalwallis, columns are: marker,dof,pvalue,kvalue,calls,means\n")
    cat("outbase: base name for output files\n")
    quit()
}

genofile  <- args[1]
phenofile <- args[2]
pheno_col <- args[3]
kwfile    <- args[4]
outbase   <- args[5]

#example values
if(FALSE)
{
    genofile  <- "uniq_genotypes_numerical.csv"
    phenofile <- "phenotypes.csv"
    pheno_col <- "score"
    kwfile    <- "uniq_kw.csv"
    outbase   <- "test_AIC.csv"
    maxp <- 0.005
}

#source("~/git_repos/qtl_analysis/qtl_funcs.R")
source("~/rjv_mnt/cluster/git_repos/qtl_analysis/qtl_funcs.R")

#load genotypes from csv file
## column headers required: marker,type,<sample1name>,<sample2name>...
##marker=marker name string
##type=any string
##sample names are strings, must be unique and match exactly those in the phenotype file
##then one marker per line, genotype codes as integers
genotypes <- read.csv(genofile,header=T)
genotypes$marker <- gsub("-","_",genotypes$marker)

#load phenotypes
##column headers required: sample,<phenotype1name>[,<phenotype2name>...]
phenotypes <- read.csv(phenofile,header=T)

#match up genotype to phenotype using sample names, retain only the phenotype column of interest
res <- matchup_and_extract(genotypes,phenotypes,pheno_col)

#load KruskalWallis values
kwvals <- read.csv(kwfile,header=T)
kwvals$marker <- gsub("-","_",kwvals$marker)

#select significant markers by pvalue
sig_markers <- kwvals$marker[kwvals$pvalue<=maxp]

marker_names <- paste0(sig_markers,collapse='+')
fullformula <- as.formula(paste0("score~",marker_names))

qtl_model <- lm(formula=fullformula,data=res$comb)

step_model <- step(qtl_model,
                   scope=list(upper=fullformula,lower=as.formula(score~1)),
                   direction="both",
                   steps=9999)
