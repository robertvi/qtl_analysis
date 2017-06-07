#!/usr/bin/Rscript

#
# plot pvalues and/or effect sizes versus map position
#

library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)

if(length(args) < 7)
{
    cat("usage: plot_qtl_vs_mapposn.R <map_file> <output_file> [<data_file> <marker_column> <data_column> <method_label> <data_label>]... \n\n")
    cat("map_file: no header csv, three columns required: markerid, chromosome, position\n")
    cat("output_file: PNG/JPG file to output plot to\n")
    cat("data_file: csv with header\n")
    cat("marker_column: header of column containing marker ids (must match those in the map file)\n")
    cat("data_column: header of column containing the data to be plotted\n")
    cat("method_label: method used to create the data eg kruskalwallis or lasso, used to colour the plots\n")
    cat("data_label: type of data, eg pvalue, effect, used to group plots into facets\n")

    quit()
}

mapfile <- args[1]
outfile <- args[2]

cat(mapfile,"\n")
cat(outfile,"\n")

#load map
dfmap <- read.csv(mapfile,header=F)
colnames(dfmap) <- c("marker","chrom","posn")

#load qtl data
i <- 3
data <- NULL
while(i < length(args))
{
    #get next five options from command line
    datafile <- args[i]
    markercol <- args[i+1]
    datacol <- args[i+2]
    methodlab <- args[i+3]
    datalab <- args[i+4]
    i <- i + 5

    cat(datafile,"\n")

    #load and retain just the two required columns
    newdata <- read.csv(datafile,header=T)[c(markercol,datacol)]
    colnames(newdata) <- c("marker","data_value_column")
    newdata$method_label <- methodlab
    newdata$data_label <- datalab

    if(substr(datalab,1,4) == "pval") newdata$data_value_column <- -log10(newdata$data_value_column)

    #merge in map positions
    merged <- merge(newdata,dfmap,by.x="marker",by.y="marker")

    data <- rbind(data,merged)
}

ggplot(data,aes(x=posn,y=data_value_column,colour=method_label)) +
    geom_point(size=0.2,alpha=0.45) +
    facet_grid(data_label ~ chrom,scales="free_y") +
    geom_hline(yintercept = 0, colour="grey")
    #xlab("Position (bp)")
    #ylab("-log10(p)")
    #geom_line() +

ggsave(outfile,width=30,height=6)
