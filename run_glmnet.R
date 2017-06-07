#!/usr/bin/Rscript

#
# estimate SNP effects using elastic net / lasso
#

args <- commandArgs(trailingOnly = TRUE)

if(length(args) != 7)
{
    cat("usage: run_glmnet.R <genotype_file> <phenotype_file> <phenotype> <alpha> <nfold> <output_basename>\n\n")
    cat("genotype file csv column headings: marker,type,<sample1name>,<sample2name>\n")
    cat("phenotype file csv column headings: sample,<phenotype1name>[,<phenotype2name>...]\n")
    cat("phenotype: which phenotype columns to process (must match a column heading in the phenotype file)\n")
    cat("alpha: 0=>ridge regression, 1=>lasso 0<alpha<1=>elastic net\n")
    cat("nfold: nfold cross validation, eg 10 means split into ten parts, nine for fitting and one for testing\n")
    cat("reps: number of replicates of cross validation used to pick best lambda value\n")
    cat("outbase: base name for output files\n")
    quit()
}

genofile  <- args[1]
phenofile <- args[2]
pheno_col <- args[3]
alpha     <- as.double(args[4])
nfold     <- as.integer(args[5])
reps      <- as.integer(args[6])
outbase   <- args[7]


library(glmnet) #bayesian version of glmnet
#source("~/git_repos/qtl_analysis/qtl_funcs.R")
source("~/rjv_mnt/cluster/git_repos/qtl_analysis/qtl_funcs.R")

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
lambda_1se <- c()
for(i in 1:reps)
{
    cvfit <- cv.glmnet(res$mat_x,res$vec_y,alpha=alpha,nfold=nfold,family="gaussian")
    lambda_1se[i] <- cvfit$lambda.1se
}
best_lambda <- mean(lambda_1se)

#fit the model for real
fit <- glmnet(res$mat_x,res$vec_y,lambda=best_lambda)

#pred <- predict(fit, res$mat_x, type="link")[,1]
coeffs <- as.data.frame(as.matrix(coef(fit)))
coeffs <- subset(coeffs,s0 != 0.0)
coeffs <- coeffs[order(-coeffs$s0), , drop = FALSE]
colnames(coeffs) <- "value"

new_row <- data.frame(best_lambda)
colnames(new_row) <- "value"
rownames(new_row) <- c("lambda")

coeffs <- rbind(new_row,coeffs)
coeffs$marker <- rownames(coeffs)
coeffs <- coeffs[c("marker","value")]

#save all data to file
routput = gsub(" ","_",paste0(outbase,"_",pheno_col,"_glmnet",alpha,".Rdata"))
save.image(file=routput)

#save just the estimated effects information to csv file
csvoutput = gsub(" ","_",paste0(outbase,"_",pheno_col,"_glmnet",alpha,".csv"))
write.csv(coeffs,file=csvoutput,quote=F,row.names=F)
