#!/usr/bin/Rscript

#
# estimate SNP effects using elastic net
# use cross validation to select lambda and alpha
# load data from plink's "--recode A" .raw output format
# all missing genotypes must have been imputed before running this script
# (see impute_missing_data.R)
#

args <- commandArgs(trailingOnly = TRUE)

if(length(args) != 5)
{
    cat("usage: run_glmnet_cv_alpha_lambda_Araw.R <input_file> <nfold> <reps> <alpha_inc> <output_basename>\n\n")
    cat("genotype_file csv column headings: marker,type,<sample1name>,<sample2name>\n")
    cat("nfold: nfold cross validation, eg 10 means split into ten parts, nine for fitting and one for testing\n")
    cat("reps: number of replicates of cross validation used to pick best lambda value\n")
    cat("alpha_inc: how to space alpha values apart, eg 0.1 -> 0.1,0.2,0.3...1.0\n")
    cat("outbase: base name for output files\n")
    quit()
}

inpfile   <- args[1]
nfold     <- as.integer(args[2])
reps      <- as.integer(args[3])
alpha_inc <- as.double(args[4])
outbase   <- args[5]

library(glmnet)
#source("~/git_repos/qtl_analysis/qtl_funcs.R")

#load from plink --recode A .raw file
#columns must be: FID IID PAT MAT SEX PHENOTYPE marker1 marker2...
data <- read.csv(inpfile,header=T,sep='')

n_markers <- ncol(data) - 6
n_samples <- nrow(data)

cat("loaded markers:",n_markers," samples:",n_samples,"\n")

#drop progeny with missing phenotypes
data <- data[!is.na(data$PHENOTYPE),]

cat(n_samples-nrow(data),"samples dropped due to missing phenotype\n")

n_samples <- nrow(data)

marker_names <- colnames(data)[7:ncol(data)]

#extract just the genotypes as a numerical matrix
mat_x <- data.matrix(data[,7:ncol(data)])

#extract phenotype as a vector
vec_y <- as.numeric(data$PHENOTYPE)

#assign sample to cv folds in advance to allow the same folds to be used for different alphas
foldid <- NULL
for(i in 1:reps)
{
    foldid <- rbind(foldid,sample(1:nfold,n_samples,replace=T))
}

#try a range of alpha values
alpha_list <- seq(alpha_inc,1.0,alpha_inc)

df <- NULL

for(alpha in alpha_list)
{
    #multiple repeats to average out affects of different cross validation folds
    for(i in 1:reps)
    {
        cvfit <- cv.glmnet(mat_x,vec_y,alpha=alpha,foldid=foldid[i,],family="gaussian")

        #collect results into a single dataframe
        for(j in 1:length(cvfit$lambda))
        {
            new_row <- data.frame(j,alpha,cvfit$lambda[j],cvfit$cvm[j],"cvm")
            colnames(new_row) <- c("rep","alpha","lambda","value","type")
            df <- rbind(df,new_row)

            new_row <- data.frame(j,alpha,cvfit$lambda[j],cvfit$cvsd[j],"cvsd")
            colnames(new_row) <- c("rep","alpha","lambda","value","type")
            df <- rbind(df,new_row)
        }
    }
}

df2 <- NULL
for(alpha in alpha_list)
{
    #average cvm (cross validation MSE) across replicates
    mean_by_lambda <- aggregate(. ~ lambda, df[df$alpha == alpha & df$type == "cvm",], mean)
    mean_by_lambda <- mean_by_lambda[order(mean_by_lambda$lambda),,drop=F]

    #average cvsd (standard deviations of cross validation MSE) across replicates
    sd_by_lambda <- aggregate(. ~ lambda, df[df$alpha == alpha & df$type == "cvsd",], mean)
    sd_by_lambda <- sd_by_lambda[order(sd_by_lambda$lambda),,drop=F]

    #find minimum cv error for this value of alpha
    i <- which.min(mean_by_lambda$value)
    cvm_mse <- mean_by_lambda$value[i]

    #calculate minimum cv error plus one sd
    max_cvm <- mean_by_lambda$value[i] + sd_by_lambda$value[i]

    #find largest lambda where cvm is less than max_cvm
    candidate_lambda <- mean_by_lambda[which(mean_by_lambda$value < max_cvm),]

    cvm_1se <- mean_by_lambda$value[nrow(candidate_lambda)]
    lambda_1se <- candidate_lambda$lambda[nrow(candidate_lambda)]
    cvsd_1se <- sd_by_lambda$value[nrow(candidate_lambda)]

    #record alpha,lambda,cvm,cvsd
    new_row <- data.frame(alpha=alpha,lambda=lambda_1se,mse=cvm_1se,sd=cvsd_1se)

    df2 <- rbind(df2,new_row)
}

#find best alpha and lambda
i <- which.min(df2$mse)
best_alpha <- df2$alpha[i]
best_lambda <- df2$lambda[i]

#fit the model to all the data
fit <- glmnet(mat_x,vec_y,lambda=best_lambda,alpha=best_alpha)

coeffs <- as.data.frame(as.matrix(coef(fit)))
coeffs <- subset(coeffs,s0 != 0.0)
coeffs <- coeffs[order(-abs(coeffs$s0)), , drop = FALSE]
colnames(coeffs) <- "value"

new_row <- data.frame(best_lambda)
colnames(new_row) <- "value"
rownames(new_row) <- c("lambda")
coeffs <- rbind(new_row,coeffs)

new_row <- data.frame(best_alpha)
colnames(new_row) <- "value"
rownames(new_row) <- c("alpha")
coeffs <- rbind(new_row,coeffs)

coeffs$marker <- rownames(coeffs)
row.names(coeffs) <- NULL
coeffs <- coeffs[c("marker","value")]

#save all data to file
routput = gsub(" ","_",paste0(outbase,"_glmnet_cvalpha.Rdata"))
save.image(file=routput)

#save just the estimated effects information to csv file
csvoutput = gsub(" ","_",paste0(outbase,"_glmnet_cvalpha.csv"))
write.csv(coeffs,file=csvoutput,quote=F,row.names=F)
