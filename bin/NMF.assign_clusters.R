#!/usr/bin/env Rscript
#
# Wrriten by: Binbin Shi <ltbyshi@gmail.com>
# Copyright 2017 Binbin Shi
#
# This R script is used to select cluster components from the coefficient matrix after NMF
#

suppressPackageStartupMessages(library("argparse"))
parser <- ArgumentParser()
parser$add_argument("-i", "--input", required=TRUE, help="coefficient matrix")
parser$add_argument("-o", "--output", required=TRUE, help="output file for cluster components")
parser$add_argument("-t", "--threshold", type="double", default=0.2,
    help="the threshold for coefficient matrix values to define cluster components. [default = %(default)s]")
args <- parser$parse_args()

coef <- as.matrix(read.table(args$input, header=TRUE))
coef <- t(t(coef)/colSums(coef))

components = character(nrow(coef))
for(i in 1:nrow(coef)){
    components[i] <- paste(colnames(coef)[coef[i,] > args$threshold], collapse=',')
}
df <- data.frame(seq(nrow(coef)), components)
write.table(df, file=args$output, sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)
