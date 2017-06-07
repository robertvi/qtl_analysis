#!/usr/bin/Rscript

#
# estimate SNP effects using empirical bayesian elastic net / lasso
# EBlasso-NE assigns more probability mass to the tails, resulting in more nonzero estimates having large p-values
# EBlasso-NEG results in higher probability mass around zero, thus more sparse results in the final outcome
# EBEN encourages a grouping effect such that highly correlated variables can be selected as a group
# model is one of lassoNEG lasso elasticnet
#

args <- commandArgs(trailingOnly = TRUE)

if(length(args) != 6)
{
    cat("usage: run_ebglmnet.R <genotype_file> <phenotype_file> <phenotype> <model> <nfold> <output_basename>\n\n")
    cat("genotype file csv column headings: marker,type,<sample1name>,<sample2name>\n")
    cat("phenotype file csv column headings: sample,<phenotype1name>[,<phenotype2name>...]\n")
    cat("phenotype: which phenotype columns to process (must match a column heading in the phenotype file)\n")
    cat("model: one of lassoNEG, lasso, elasticnet\n")
    cat("nfold: nfold cross validation, eg 10 means split into ten parts, nine for fitting and one for testing\n")
    cat("outbase: base name for output files\n")
    quit()
}

genofile  <- args[1]
phenofile <- args[2]
pheno_col <- args[3]
model     <- args[4]
nfold     <- as.integer(args[5])
outbase   <- args[6]

#model    <- "lasso"
#crotfile <- "crot_scores_alt.csv"
#genofile <- "emxfe_crot_vescax4/emxfe_genotypes_cult.csv"
#nfold <- 10  #how many parts to split the data into for cross validation
#outfile  <- "emxfe_genotypes_cult_lasso.Rdata"

library(EBglmnet) #bayesian version of glmnet
#source("~/git_repos/qtl_analysis/qtl_funcs.R")
source("~/rjv_mnt/cluster/git_repos/qtl_analysis/qtl_funcs.R")

if(model == "elasticnet") model <- "elastic net";

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

res <- matchup_and_extract(genotypes,phenotypes,pheno_col)

#fit model
cvfit <- cv.EBglmnet(res$mat_x,res$vec_y,prior=model,nfold=nfold)

#save all data to file
routput = gsub(" ","_",paste0(outbase,"_",pheno_col,"_ebglmnet",model,".Rdata"))
save.image(file=routput)

#save just the estimated effects information to csv file
fit <- as.data.frame(cvfit$fit)
fit <- fit[order(-fit["t-value"]),]
fit$marker <- colnames(res$mat_x)[fit$locus1]
fit <- fit[c("marker","p-value","t-value","beta","posterior variance")]
colnames(fit) <- c("marker","pvalue","tvalue","beta","posterior variance") #prevent - being changed to . when reloading
csvoutput = gsub(" ","_",paste0(outbase,"_",pheno_col,"_ebglmnet",model,".csv"))
write.csv(fit,quote=F,file=csvoutput)

#save just the info needed to predict phenotype values
coeffs <- data.frame(as.character(fit$marker),fit$beta)
colnames(coeffs) <- c("marker","value")
coeffs$marker <- as.character(coeffs$marker)
coeffs <- rbind(c("(Intercept)",cvfit$Intercept),coeffs)
modeloutput = gsub(".csv","_model.csv",csvoutput)
write.csv(coeffs,file=modeloutput,quote=F,row.names=F)
