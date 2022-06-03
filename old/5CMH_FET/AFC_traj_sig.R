#!/usr/bin/env Rscript

###find the allele frequency change in isgnificant sites

#### loading packages####
setwd("~/Dropbox (PopGen)/Yiwen/RS/data/")
require(ggplot2)
require(dplyr)
require(poolSeq)
require(parallel)
require(ACER)
library(viridis)
library(factoextra)
library(reshape2)
library(data.table)
require(gridExtra)
require(ggrepel)
require(permute)
library(ggpubr)
#get sig or not sites
res_r11 <- readRDS("./res_r11")
thres_cmh_genome <- 1e-7
thres_cmh_suggest <- 5e-7

sign_sites <- res_r11[which(res_r11$P<1e-7 & res_r11$P > 5e-8),]
top_sites <- res_r11[order(res_r11$P),][1:100,]
low_sites <- res_r11[which(res_r11$P>1e-7),]; low_sites <- low_sites[sort(sample.int(nrow(low_sites), 100)), ]

#get allele frequency matrix
Repl <- rep(c(1:5), times=2)
Gen <- rep(c(0,30),each=5)
Sync.data <- read.sync(file="./cold_F30_r11_all_chr_mq20_bq20_Neda_filtered_base_duplicated_nozerocov.sync", 
                       repl=Repl, gen=Gen, polarization = "rising")
Allele_f <- as.data.frame(af(Sync.data, repl=Repl, gen=Gen))
Coverage <- as.data.frame(coverage(Sync.data, repl=Repl, gen=Gen))
Alleles <- as.data.frame(alleles(Sync.data))
row.names(Alleles) <- row.names(Coverage)
Allele_f <- na.omit(Allele_f)
rm(Coverage, Alleles, Sync.data, Gen, Repl, sig_r11_s1, sig_r11_s2, sig_r11_s3, sig_r11_s4, sig_r11_s5)

sign_sites_af <- Allele_f[rownames(sign_sites),]
top_sites_af <- Allele_f[rownames(top_sites),]
low_sites_af <- Allele_f[rownames(low_sites),]

sign_sites_af$max.freq <- apply(sign_sites_af, 1, max); sign_sites_af$min.freq <- apply(sign_sites_af, 1, min)
top_sites_af$max.freq <- apply(top_sites_af, 1, max); top_sites_af$min.freq <- apply(top_sites_af, 1, min)
low_sites_af$max.freq <- apply(low_sites_af, 1, max); low_sites_af$min.freq <- apply(low_sites_af, 1, min)

sign_sites_af$mean.start.freq <- apply(sign_sites_af[,c(1,3,5,7,9)], 1, mean); sign_sites_af$mean.end.freq <- apply(sign_sites_af[,c(2,4,6,8,10)], 1, mean)
top_sites_af$mean.start.freq <- apply(top_sites_af[,c(1,3,5,7,9)], 1, mean); top_sites_af$mean.end.freq <- apply(top_sites_af[,c(2,4,6,8,10)], 1, mean)
low_sites_af$mean.start.freq <- apply(low_sites_af[,c(1,3,5,7,9)], 1, mean); low_sites_af$mean.end.freq <- apply(low_sites_af[,c(2,4,6,8,10)], 1, mean)


### plotting
sign_sites_af_start <- cbind(SNP = rownames(sign_sites_af), Freq = sign_sites_af[,13], State = 0)
sign_sites_af_end <- cbind(SNP = rownames(sign_sites_af), Freq = sign_sites_af[,14], State = 30)
sign_sites_af_start_end <- as.data.frame(rbind(sign_sites_af_start, sign_sites_af_end))
sign_sites_af_start_end$Freq <- as.numeric(sign_sites_af_start_end$Freq)

sign_plot <- ggplot(data = sign_sites_af_start_end, aes(x = State, y = Freq, group = SNP)) +
  geom_line( size = 0.1) +
  geom_point( size = 0.5) +
  theme_bw() +
  theme(legend.position="none") +
  xlab("Generation") + ylab("Allele frequency") +
  geom_hline(yintercept = 0.05, linetype="dotted", color = "yellow", size=0.5) + 
  geom_hline(yintercept = 0.10, linetype="dotted", color = "orange", size=0.5) + 
  geom_hline(yintercept = 0.15, linetype="dotted", color = "red", size=0.5) + 
  ggtitle("Allele frequency change from start to end in significant SNPs near threshold")
tiff("sign_SNPs_AFC.tiff",  width = 18, height = 10, res = 300, units = "cm")
print(sign_plot)
dev.off()

top_sites_af_start <- cbind(SNP = rownames(top_sites_af), Freq = top_sites_af[,13], State = 0)
top_sites_af_end <- cbind(SNP = rownames(top_sites_af), Freq = top_sites_af[,14], State = 30)
top_sites_af_start_end <- as.data.frame(rbind(top_sites_af_start, top_sites_af_end))
top_sites_af_start_end$Freq <- as.numeric(top_sites_af_start_end$Freq)

top_plot <- ggplot(data = top_sites_af_start_end, aes(x = State, y = Freq, group = SNP)) +
  geom_line( size = 0.1) +
  geom_point( size = 0.5) +
  theme_bw() +
  theme(legend.position="none") +
  xlab("Generation") + ylab("Allele frequency") +
  geom_hline(yintercept = 0.05, linetype="dotted", color = "yellow", size=0.5) + 
  geom_hline(yintercept = 0.10, linetype="dotted", color = "orange", size=0.5) + 
  geom_hline(yintercept = 0.15, linetype="dotted", color = "red", size=0.5) + 
  ggtitle("Allele frequency change from start to end in top 100 SNPs")
tiff("top_SNPs_AFC.tiff",  width = 18, height = 10, res = 300, units = "cm")
print(top_plot)
dev.off()


low_sites_af_start <- cbind(SNP = rownames(low_sites_af), Freq = low_sites_af[,13], State = 0)
low_sites_af_end <- cbind(SNP = rownames(low_sites_af), Freq = low_sites_af[,14], State = 30)
low_sites_af_start_end <- as.data.frame(rbind(low_sites_af_start, low_sites_af_end))
low_sites_af_start_end$Freq <- as.numeric(low_sites_af_start_end$Freq)

low_plot <- ggplot(data = low_sites_af_start_end, aes(x = State, y = Freq, group = SNP)) +
  geom_line( size = 0.1) +
  geom_point( size = 0.5) +
  theme_bw() +
  theme(legend.position="none") +
  xlab("Generation") + ylab("Allele frequency") +
  geom_hline(yintercept = 0.05, linetype="dotted", color = "yellow", size=0.5) + 
  geom_hline(yintercept = 0.10, linetype="dotted", color = "orange", size=0.5) + 
  geom_hline(yintercept = 0.15, linetype="dotted", color = "red", size=0.5) + 
  ggtitle("Allele frequency change from start to end in low 100 SNPs")
tiff("low_SNPs_AFC.tiff",  width = 18, height = 10, res = 300, units = "cm")
print(low_plot)
dev.off()


####plot dist of max af


#get sig or not sites
res_r11 <- readRDS("./res_r11")
thres_cmh_genome <- 1e-7
thres_cmh_suggest <- 5e-7

sign_sites <- res_r11[which(res_r11$P<1e-7 & res_r11$P > 5e-8),]
top_sites <- res_r11[order(res_r11$P),][1:100,]
low_sites <- res_r11[which(res_r11$P>1e-7),]; low_sites <- low_sites[sort(sample.int(nrow(low_sites), 100)), ]

#get allele frequency matrix
Repl <- rep(c(1:5), times=2)
Gen <- rep(c(0,30),each=5)
Sync.data <- read.sync(file="./cold_F30_r11_all_chr_mq20_bq20_Neda_filtered_base_duplicated_nozerocov.sync", 
                       repl=Repl, gen=Gen, polarization = "rising")
Allele_f <- as.data.frame(af(Sync.data, repl=Repl, gen=Gen))
Coverage <- as.data.frame(coverage(Sync.data, repl=Repl, gen=Gen))
Alleles <- as.data.frame(alleles(Sync.data))
row.names(Alleles) <- row.names(Coverage)
Allele_f <- na.omit(Allele_f)
rm(Coverage, Alleles, Sync.data, Gen, Repl, sig_r11_s1, sig_r11_s2, sig_r11_s3, sig_r11_s4, sig_r11_s5)

sign_sites_af <- Allele_f[rownames(sign_sites),]
top_sites_af <- Allele_f[rownames(top_sites),]
low_sites_af <- Allele_f[rownames(low_sites),]

sign_sites_af$max.freq <- apply(sign_sites_af, 1, max); sign_sites_af$min.freq <- apply(sign_sites_af, 1, min)
top_sites_af$max.freq <- apply(top_sites_af, 1, max); top_sites_af$min.freq <- apply(top_sites_af, 1, min)
low_sites_af$max.freq <- apply(low_sites_af, 1, max); low_sites_af$min.freq <- apply(low_sites_af, 1, min)


sign_max_af <- ggplot(sign_sites_af, aes(x=sign_sites_af$max.freq)) +
  geom_histogram(binwidth=0.002) +
  xlim(0, 0.30) + xlab("Maximum AF across samples of near threshold significant SNPs") +
  ylab("Count")
top_max_af <- ggplot(top_sites_af, aes(x=top_sites_af$max.freq)) + 
  geom_histogram(binwidth=0.002) +
  xlim(0, 0.30) + xlab("Maximum AF across samples of top significant SNPs") +
  ylab("Count")
low_max_af <- ggplot(low_sites_af, aes(x=low_sites_af$max.freq)) + 
  geom_histogram(binwidth=0.002) +
  xlim(0, 0.30) + xlab("Maximum AF across samples of non-significant SNPs") +
  ylab("Count")

tiff("Max_AF_three_types_SNPs.tiff",  width = 24, height = 16, res = 450, units = "cm")
ggarrange(sign_max_af, top_max_af, low_max_af, labels = c("A", "B", "C"), ncol = 2, nrow = 2)
dev.off()
