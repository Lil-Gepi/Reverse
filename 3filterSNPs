#!/bin/sh
# to filter the SNPs piled up in our mapped data based on Neda's SNP set.
# this awk command line is from Sheng-kai

awk 'BEGIN{OFS="\t"} NR==FNR {f1[$1$2] = $0; next} ($1$2 in f1) {print $0}' \
~/RS/F30_all/reference/Dsim_base_F60_bam.sorted.rmdup.filter.mpileup.Q20.polymorphic_indel_repeat_Ytranslocremoved_masked_MajorChr.sync_furthermasked.cmh \
~/RS/F30_all/result/F30_r11_r12_all_chr_mq20_bq20.sync > \
~/RS/F30_all/result/F30_r11_r12_all_chr_mq20_bq20_Neda_filtered.sync
