#!/bin/bash
samtools --version
bcftools --version
which gawk
which awk

for i in {10,20,30,40,50,60,80,100,120,140}; do
  for j in {11..14}; do
    for k in {1..5}; do
    mkdir "/Users/ychen/RS/inter_data/F${i}_r${j}sub${k}/"
    done
  done
done

for j in {11..20}; do
  mkdir "/Users/ychen/RS/inter_data/F0_r${j}/"
done
