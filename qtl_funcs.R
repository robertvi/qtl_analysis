#match up phenotypes to genotypes via progeny names
#extract a single phenotype out of one or more columns of phenotype data
#drop progeny which are missing the selected phenotype
#genotypes are assumed to all be present or already imputed

#genotype dataframe fields required: marker,type,<sample1name>,<sample2name>...
#marker=marker name string
#type=any string
#sample names are strings, must be unique and match exactly those in the phenotype file
#then one marker per line, genotype codes as integers


#phenotype dataframe fields required: sample,<phenotype1name>[,<phenotype2name>...]
#pcolumn string must match one of the phenotype column names

matchup_and_extract <- function(genotypes,phenotypes,pcolumn)
{
    #match up phenotype to genotype by progeny name
    scores <- t(phenotypes[pcolumn])
    colnames(scores) <- phenotypes$sample
    scores <- cbind(NA,NA,scores)
    colnames(scores)[1:2] <- c("marker","type")
    comb <- rbind(scores,genotypes)

    #separate marker type info from phenotype and genotype data
    marker_type <- comb[2:nrow(comb),1:2]
    comb <- as.data.frame(t(comb[,3:ncol(comb)]))
    colnames(comb)[2:ncol(comb)] <- marker_type$marker

    #drop progeny with missing phenotypes
    comb <- comb[!is.na(comb[pcolumn]),]

    #extract just the genotypes as a numerical matrix
    mat_x <- data.matrix(comb[,2:ncol(comb)])

    #extract phenotype as a vector
    vec_y <- as.numeric(comb[,1])

    return (list(mat_x=mat_x, vec_y=vec_y, comb=comb, marker_type=marker_type))
}
