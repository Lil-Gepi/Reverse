#!/usr/bin/env Rscript
#### loading packages####
setwd("~/Dropbox (PopGen)/Yiwen/RS/data/")
require(ggplot2)
require(dplyr)

thres_fet <- 0.05
res_r11_s1 <- readRDS(file = "./res_r11_s1"); sig_r11_s1 <- res_r11_s1[which(res_r11_s1$P < thres_fet),4]; rm(res_r11_s1)
res_r11_s2 <- readRDS(file = "./res_r11_s2"); sig_r11_s2 <- res_r11_s2[which(res_r11_s2$P < thres_fet),4]; rm(res_r11_s2)
res_r11_s3 <- readRDS(file = "./res_r11_s3"); sig_r11_s3 <- res_r11_s3[which(res_r11_s3$P < thres_fet),4]; rm(res_r11_s3)
res_r11_s4 <- readRDS(file = "./res_r11_s4"); sig_r11_s4 <- res_r11_s4[which(res_r11_s4$P < thres_fet),4]; rm(res_r11_s4)
res_r11_s5 <- readRDS(file = "./res_r11_s5"); sig_r11_s5 <- res_r11_s5[which(res_r11_s5$P < thres_fet),4]; rm(res_r11_s5)

res_r11 <- readRDS("./res_r11")
thres_cmh_genome <- 1e-7
thres_cmh_suggest <- 5e-7


# Plotting----
####highlight s1----
result <- res_r11
don <- result %>% 
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  left_join(result, ., by=c("CHR"="CHR")) %>%
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot) %>%
# Add highlight and annotation information
  mutate( is_highlight=ifelse(SNP %in% sig_r11_s1, "yes", "no"))

axisdf = don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
pdf("./Plots/r11_plots/CMH_adapted_r11_F0-F30_fdr_adjusted_rising_highlighted_r11s1.pdf", width = 10, height = 6)
ggplot(don, aes(x=BPcum, y=-log10(P))) +
  ggtitle("Adapted CMH test in r11 contrasting F0 and F30, highlighting sig SNPs in FET of r11s1")+
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1) +
  scale_color_manual(values = rep(c("grey", "black"), 22 )) +
  scale_x_continuous(name = "Chromosome", label = c("X", "2L","2R","3L", "3R", "4"), breaks= axisdf$center ) +
  scale_y_continuous(name = "-log10(q)", expand = c(0, 0) ) +     
  geom_hline(yintercept=-log10(thres_cmh_genome), linetype="dashed", color = "red") +
  geom_hline(yintercept=-log10(thres_cmh_suggest), linetype="dashed", color = "blue") +
  geom_point(data=subset(don, is_highlight=="yes"), color="tomato", size=1) +
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
dev.off()

####highlight s2----
result <- res_r11
don <- result %>% 
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  left_join(result, ., by=c("CHR"="CHR")) %>%
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot) %>%
  # Add highlight and annotation information
  mutate( is_highlight=ifelse(SNP %in% sig_r11_s2, "yes", "no"))

axisdf = don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
pdf("./Plots/r11_plots/CMH_adapted_r11_F0-F30_fdr_adjusted_rising_highlighted_r11s2.pdf", width = 10, height = 6)
ggplot(don, aes(x=BPcum, y=-log10(P))) +
  ggtitle("Adapted CMH test in r11 contrasting F0 and F30, highlighting sig SNPs in FET of r11s2")+
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1) +
  scale_color_manual(values = rep(c("grey", "black"), 22 )) +
  scale_x_continuous(name = "Chromosome", label = c("X", "2L","2R","3L", "3R", "4"), breaks= axisdf$center ) +
  scale_y_continuous(name = "-log10(q)", expand = c(0, 0) ) +     
  geom_hline(yintercept=-log10(thres_cmh_genome), linetype="dashed", color = "red") +
  geom_hline(yintercept=-log10(thres_cmh_suggest), linetype="dashed", color = "blue") +
  geom_point(data=subset(don, is_highlight=="yes"), color="springgreen", size=1) +
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
dev.off()

####highlight s3----
result <- res_r11
don <- result %>% 
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  left_join(result, ., by=c("CHR"="CHR")) %>%
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot) %>%
  # Add highlight and annotation information
  mutate( is_highlight=ifelse(SNP %in% sig_r11_s3, "yes", "no"))

axisdf = don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
pdf("./Plots/r11_plots/CMH_adapted_r11_F0-F30_fdr_adjusted_rising_highlighted_r11s3.pdf", width = 10, height = 6)
ggplot(don, aes(x=BPcum, y=-log10(P))) +
  ggtitle("Adapted CMH test in r11 contrasting F0 and F30, highlighting sig SNPs in FET of r11s3")+
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1) +
  scale_color_manual(values = rep(c("grey", "black"), 22 )) +
  scale_x_continuous(name = "Chromosome", label = c("X", "2L","2R","3L", "3R", "4"), breaks= axisdf$center ) +
  scale_y_continuous(name = "-log10(q)", expand = c(0, 0) ) +     
  geom_hline(yintercept=-log10(thres_cmh_genome), linetype="dashed", color = "red") +
  geom_hline(yintercept=-log10(thres_cmh_suggest), linetype="dashed", color = "blue") +
  geom_point(data=subset(don, is_highlight=="yes"), color="tan3", size=1) +
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
dev.off()

####highlight s4----
result <- res_r11
don <- result %>% 
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  left_join(result, ., by=c("CHR"="CHR")) %>%
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot) %>%
  # Add highlight and annotation information
  mutate( is_highlight=ifelse(SNP %in% sig_r11_s4, "yes", "no"))

axisdf = don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
pdf("./Plots/r11_plots/CMH_adapted_r11_F0-F30_fdr_adjusted_rising_highlighted_r11s4.pdf", width = 10, height = 6)
ggplot(don, aes(x=BPcum, y=-log10(P))) +
  ggtitle("Adapted CMH test in r11 contrasting F0 and F30, highlighting sig SNPs in FET of r11s4")+
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1) +
  scale_color_manual(values = rep(c("grey", "black"), 22 )) +
  scale_x_continuous(name = "Chromosome", label = c("X", "2L","2R","3L", "3R", "4"), breaks= axisdf$center ) +
  scale_y_continuous(name = "-log10(q)", expand = c(0, 0) ) +     
  geom_hline(yintercept=-log10(thres_cmh_genome), linetype="dashed", color = "red") +
  geom_hline(yintercept=-log10(thres_cmh_suggest), linetype="dashed", color = "blue") +
  geom_point(data=subset(don, is_highlight=="yes"), color="plum", size=1) +
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
dev.off()

####highlight s5----
result <- res_r11
don <- result %>% 
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  left_join(result, ., by=c("CHR"="CHR")) %>%
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot) %>%
  # Add highlight and annotation information
  mutate( is_highlight=ifelse(SNP %in% sig_r11_s5, "yes", "no"))

axisdf = don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
pdf("./Plots/r11_plots/CMH_adapted_r11_F0-F30_fdr_adjusted_rising_highlighted_r11s5.pdf", width = 10, height = 6)
ggplot(don, aes(x=BPcum, y=-log10(P))) +
  ggtitle("Adapted CMH test in r11 contrasting F0 and F30, highlighting sig SNPs in FET of r11s5")+
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1) +
  scale_color_manual(values = rep(c("grey", "black"), 22 )) +
  scale_x_continuous(name = "Chromosome", label = c("X", "2L","2R","3L", "3R", "4"), breaks= axisdf$center ) +
  scale_y_continuous(name = "-log10(q)", expand = c(0, 0) ) +     
  geom_hline(yintercept=-log10(thres_cmh_genome), linetype="dashed", color = "red") +
  geom_hline(yintercept=-log10(thres_cmh_suggest), linetype="dashed", color = "blue") +
  geom_point(data=subset(don, is_highlight=="yes"), color="steelblue", size=1) +
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
dev.off()







