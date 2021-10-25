#!/usr/bin/env Rscript

library(poolSeq)
F30_r11_r12 <- read.sync(file="~/RS/F30_all/result/F30_r11_r12_all_chr_mq20_bq20.sync",
               repl=rep(c(11,12), each=5), sub=rep(c(1,2,3,4,5),times=2), rising = FALSE)
F30_r11_r12_freq <- af(F30_r11_r12, repl=rep(c(11,12), each=5), sub=rep(c(1,2,3,4,5),times=2));
F30_r11_r12_freq <- cbind(F30_r11_r12_freq, splitLocusID(rownames(F30_r11_r12_freq)))
