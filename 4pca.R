#!/usr/bin/env Rscript

library(poolSeq)
F30 <- read.sync(file="~/RS/F30_all/result/F30_all_chr_mq20_bq20.sync",
               repl=rep(c(11,12,13,14), each=5), sub=rep(c(1,2,3,4,5),times=2), rising = FALSE)

F30_freq <- af(F30, repl=rep(c(11,12,13,14), each=5), sub=rep(c(1,2,3,4,5),times=2));
F30_freq <- cbind(F30_freq, splitLocusID(rownames(F30_freq)))
