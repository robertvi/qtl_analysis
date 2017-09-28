#!/usr/bin/Rscript

#
# estimate SNP effects using stepwise AIC
# uses already calculated Kruskal Wallis values
#

args <- commandArgs(trailingOnly = TRUE)

if(length(args) != 6)
{
    cat("usage: run_stepwise_AIC.R <genotype_file> <phenotype_file> <phenotype> <kw_file> <output_basename>\n\n")
    cat("genotype_file csv, numerical genotype codes, column headings: marker,type,<sample1name>,<sample2name>\n")
    cat("phenotype_file csv, numerical phenotype values, column headings: sample,<phenotype1name>[,<phenotype2name>...]\n")
    cat("phenotype: which phenotype column to process (must match a column heading in the phenotype file)\n")
    cat("kw_file: file containing output from calc_kruskalwallis, columns are: marker,dof,pvalue,kvalue,calls,means\n")
    cat("pvalue: Kruskal Wallis pvalue threshold for consideration of marker in the model\n")
    cat("outbase: base name for output files\n")
    quit()
}

genofile  <- args[1]
phenofile <- args[2]
pheno_col <- args[3]
kwfile    <- args[4]
pvalue    <- as.double(args[5])
outbase   <- args[6]

#example values
#genofile  <- "uniq_genotypes_numerical.csv"
#phenofile <- "phenotypes.csv"
#pheno_col <- "score"
#kwfile    <- "uniq_kw.csv"
#outbase   <- "test_AIC.csv"
#pvalue    <- 0.001

#source("~/git_repos/qtl_analysis/qtl_funcs.R")
source("~/rjv_mnt/cluster/git_repos/qtl_analysis/qtl_funcs.R")

#load genotypes from csv file
## column headers required: marker,type,<sample1name>,<sample2name>...
##marker=marker name string
##type=any string
##sample names are strings, must be unique and match exactly those in the phenotype file
##then one marker per line, genotype codes as integers
genotypes <- read.csv(genofile,header=T)
genotypes$marker <- gsub("-","_",genotypes$marker) #make sure marker names are valid variable names

#load phenotypes
##column headers required: sample,<phenotype1name>[,<phenotype2name>...]
phenotypes <- read.csv(phenofile,header=T)

#match up genotype to phenotype using sample names, retain only the phenotype column of interest
res <- matchup_and_extract(genotypes,phenotypes,pheno_col)

#load KruskalWallis values
kwvals <- read.csv(kwfile,header=T)
kwvals$marker <- gsub("-","_",kwvals$marker)

#select significant markers by kw pvalue
sig_markers <- kwvals$marker[ kwvals$pvalue<=pvalue ]

#run step once to get list of significant markers using AIC criterion
marker_names <- paste0(sig_markers,collapse='+')
fullformula <- as.formula(paste0(pheno_col,"~",marker_names))
qtl_model <- lm(formula=fullformula,data=res$comb)
step_model <- step(qtl_model,
                   scope=list(upper=fullformula,lower=as.formula(score~1)),
                   direction="both",
                   steps=9999)

#iteratively prune markers based on significance
while(TRUE)
{
    #extract the chosen markers and their pvalues
    df <- as.data.frame(coef(summary(step_model))[-1,])
    df$marker <- rownames(df)
    df$pvalue <- df[,4]
    df <- df[,c("marker","pvalue")]
    df <- df[order(-df$pvalue),]
    rownames(df) <- 1:nrow(df)
    sig_markers <- df$marker

    print(paste0("largest p value=",df$pvalue[1]))
    print(paste0("number of markers=",length(df$marker)))

    #stop if least significant marker has pvalue <= 0.05
    if(df$pvalue[1] <= 0.05) break

    #exclude least significant marker
    sig_markers <- sig_markers[sig_markers!=df$marker[1]]

    #run step again
    marker_names <- paste0(sig_markers,collapse='+')
    fullformula <- as.formula(paste0(pheno_col,"~",marker_names))
    qtl_model <- lm(formula=fullformula,data=res$comb)
    step_model <- step(qtl_model,
                       scope=list(upper=fullformula,lower=as.formula(score~1)),
                       direction="both",
                       steps=9999)
}

mod <- as.data.frame(step_model$coefficients)
mod$marker <- gsub("_","-",rownames(mod))        #return marker names to original form
colnames(mod) <- c("value","marker")
mod <- mod[,c("marker","value")]
rownames(mod) <- 1:nrow(mod)

#save all data to file
routput = gsub(" ","_",paste0(outbase,"_",pheno_col,"_AICstep.Rdata"))
save.image(file=routput)

#save just the estimated effects information to csv file
csvoutput = gsub(" ","_",paste0(outbase,"_",pheno_col,"_AICstep.csv"))
write.csv(mod,file=csvoutput,quote=F,row.names=F)
