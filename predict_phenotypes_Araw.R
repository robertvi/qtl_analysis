#!/usr/bin/Rscript

#
# load sample genotypes and a linear model
# prediction phenotype values for all samples
# this version loads data from the plink --recodeA .raw format

args <- commandArgs(trailingOnly = TRUE)

if(length(args) != 3)
{
    cat("usage: predict_phenotypes_Araw.R <genotype_file> <model_file> <output_file>\n\n")
    cat("genotype_file csv column headings: marker,type,<sample1name>,<sample2name>\n")
    cat("model_file csv column headings: marker,value\n")
    cat("output_file: output file to receive predicted values\n")
    quit()
}

genofile  <- args[1]
modelfile <- args[2]
outfile   <- args[3]

#load genotypes from csv file
## column headers required: marker,type,<sample1name>,<sample2name>...
##marker=marker name string
##type=any string
##sample names are strings
##then one marker per line, genotype codes as integers
#genotypes <- read.csv(genofile,header=T)

#load from plink --recode A .raw file
#columns must be: FID IID PAT MAT SEX PHENOTYPE marker1 marker2...
data <- read.csv(genofile,header=T,sep='')
n_markers <- ncol(data) - 6
n_samples <- nrow(data)
cat("loaded markers:",n_markers," samples:",n_samples,"\n")

samples <- data$IID
genotypes <- data.frame(t(data[,7:ncol(data)]))
colnames(genotypes) <- samples

#load model parameters
##column headers required: marker,value
##the special marker name (Intercept) indicates the start of the model parameters
##items above are ignored (may be used for meta data like lambda value)
##items below are marker name strings and must match those in the genotype file exactly
model <- read.csv(modelfile,header=T)
rownames(model) <- 1:nrow(model)

#extract baseline and marker effects
i <- which(model$marker == "(Intercept)")
baseline <- model$value[i]
model <- model[-c(1:i),]

genotypes$marker <- rownames(genotypes)

#match up effects with marker values
comb <- merge(model,genotypes,by.x="marker",by.y="marker")

#multiply genotype by effect
for(col in colnames(comb)[3:ncol(comb)])
{
    comb[col] <- comb$value * comb[col]
}

#sum and add baseline
pred <- as.data.frame(colSums(comb[3:ncol(comb)])+baseline)
colnames(pred) <- c("prediction")
pred$sample <- rownames(pred)
pred <- pred[,c("sample","prediction")]
rownames(pred) <- 1:nrow(pred)

#save predictions
write.csv(pred,file=outfile,quote=F,row.names=F)
