#!/usr/bin/env Rscript
#
# Wrriten by: Yang Li <ly-13@mails.tsinghua.edu.cn> and
#             Binbin Shi <ltbyshi@gmail.com>
# Copyright 2017 Yang Li
#
# This R script is used to estimate the rank (R) for NMF algorithm
#

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, help="input matrix")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
parser$add_argument("-s", "--start", type="integer", default=1,
    help="Number of start rank [default %(default)s]",
    metavar="RANK")
parser$add_argument("-e", "--end", type="integer", required=TRUE,
    help="Number of end rank",
    metavar="RANK")
parser$add_argument('-m', "--method", type="character", default="KL",
    metavar="STRING",
    help="the NMF algorithm to use. Should be one of KL,euclidean,KL_ortho,euclidean_ortho. [default = %(default)s]")
parser$add_argument("-a", "--alpha", type="double", default=10.0,
    metavar="NUMBER",
    help="regularization factor for orthogonality of the coefficient matrix [default = %(default)s]")
parser$add_argument("--seed", type="integer",
    help="Seed for the random number generator",
    metavar="NUMBER")
parser$add_argument("-n", "--runs", type="integer", default=30,
    help="Number of runs to perform [default = %(default)s]",
    metavar="NUMBER")
parser$add_argument("-p", "--processors", type="integer", default=1,
    help="Number of processors to use. This option is useful on multicore *nix or Mac machine only, when performing multiple runs (nrun > 1) [default %(default)s]",
    metavar="NUMBER")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

# perform NMF
# loading matrix
data <- as.matrix(read.table(args$input,head=T))
data <- data + 1e-100

library(NMF)
# Euclidean distance with regularization for orthognnality
objective.euclidean_ortho <- function(x, y, alpha=1.0, ...)
{
    w <- .basis(x)
    h <- .coef(x)
    rss(y, fitted(x)) + alpha*rss(h %*% t(h), diag(nrow=nrow(h)))
}
nmf_update.euclidean_ortho <- function(i, v, x, alpha=1.0, eps=1e-100, ...)
{
    # retrieve each factor
    w <- .basis(x)
    h <- .coef(x)
    # update W
    w.update <- (v %*% t(h))/(w %*% h %*% t(h))
    w <- w * sqrt(w.update)
    # update H
    h.update <- ((t(w) %*% v) + (2.0*alpha*h))/((t(w) %*% w %*% h) + 2.0*alpha*(h %*% t(h) %*% h))
    h <- h * sqrt(h.update)
    # avoid numerical underflow
    w <- pmax(w, eps)
    h <- pmax(h, eps)
    #message('cost = ', objective.euclidean_ortho(x, v, alpha, ...))
    # save the model
    .basis(x) <- w
    .coef(x) <- h
    return(x)
}
setNMFMethod('euclidean_ortho',
             objective=objective.euclidean_ortho,
             Update=nmf_update.euclidean_ortho,
             Stop='connectivity')

# KL distance with regularization term for orthgonality
objective.KL_ortho <- function(x, y, alpha=1.0, ...)
{
    w <- .basis(x)
    h <- .coef(x)
    wh <- w %*% h
    cost <- y*log(y/(wh)) - y + wh + alpha*rss(h %*% t(h), diag(nrow=dim(h)[1]))
    return(cost)
}
nmf_update.KL_ortho <- function(i, v, x, alpha=1.0, eps=1e-100, ...)
{
    # retrieve each factor
    w <- .basis(x)
    h <- .coef(x)
    # update W
    x_wh <- v / (w %*% h)
    h_sum <- rowSums(h)
    w.update <- t(t(x_wh %*% t(h))/h_sum)
    w <- w * sqrt(w.update)
    # update H
    x_wh <- v / (w %*% h)
    w_sum <- colSums(w)
    h.update <- (t(w) %*% x_wh + 4.0*alpha*h)/(w_sum + 4.0*alpha*(h %*% t(h) %*% h))
    h <- h * sqrt(h.update)
    # avoid numerical underflow
    w <- pmax(w, 1e-100)
    h <- pmax(h, 1e-100)
    .basis(x) <- w
    .coef(x) <- h
    return(x)
}
setNMFMethod('KL_ortho', objective=objective.KL_ortho,
             Update=nmf_update.KL_ortho,
             Stop='connectivity')

ranks <- args$start:args$end
estim <- lapply(ranks, function(r){
    nmf_args <- list(data, r, method=args$method, nrun=args$runs, maxIter=2000,
        .opt=paste("vmp",args$processors,sep=''))
    if(!is.null(args$seed)) nmf_args$seed <- args$seed
    if(args$method == "KL") {
        nmf_args$method <- "brunet"
    }else if(args$method == "euclidean") {
        nmf_args$method <- "lee"
    }else if(args$method == "euclidean_ortho" || args$method == "KL_ortho"){
        nmf_args$alpha <- args$alpha
    }
    fit <- do.call(nmf, nmf_args)
    list(fit = fit, coph = cophcor(fit), disp = dispersion(fit), rsdl = residuals(fit), consensus = consensus(fit))
})

names(estim) <- paste('rank', ranks)

# summary plots
pdf(paste(args$output,".pdf",sep=''))
plot(ranks, sapply(estim, '[[', 'coph'),
    main = "-Cophenetic coefficient-",
    xlab="Factorization rank",
    ylab="Quality measure: cophenetic",
    type = 'b', col='Blue', lwd=2)
plot(ranks, sapply(estim, '[[', 'disp'),
    main = "-Dispersion coefficient-",
    xlab="Factorization rank",
    ylab="Quality measure: dispersion",
    type = 'b', col='Blue', lwd=2)
plot(ranks, sapply(estim, '[[', 'rsdl'),
    main = "-Residuals-",
    xlab="Factorization rank",
    ylab="Quality measure: residuals",
    type = 'b', col='Blue', lwd=2)

## consensus matrix plots 
consensus <- sapply(estim, '[[', 'consensus', simplify = FALSE)
for (i in consensus){
    consensusmap(i, info = T)
}
dev.off()

cophcor <- as.numeric(sapply(estim, '[[', 'coph'))
dispersion <- as.numeric(sapply(estim, '[[', 'disp'))
residuals <- as.numeric(sapply(estim, '[[', 'rsdl'))
outmx <- rbind(cophcor,dispersion,residuals)
colnames(outmx) <- paste("Rank",ranks,sep='')
write.table(outmx, file=paste(args$output,".txt", sep=''),
    sep='\t', quote=F, row.names=T, col.names=T)

