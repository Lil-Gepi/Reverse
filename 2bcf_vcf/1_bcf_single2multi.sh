#!/bin/bash

in_folder="/Users/ychen/RS/inter_data"
out_bcf="/Users/ychen/RS/result/piled_bcf/RS_piled.bcf"

for bcf_files in $(ls ${in_folder}/*/*.bcf)
do
  echo $bcf_files
done | parallel -j 5 "time bcftools index -f --threads 4 {}"

bcftools merge --no-index --merge both --threads 16 -Ov -o - ${in_folder}/*/*.bcf |\
  /Users/ychen/RS/pipeline/rupert_pipe_v3/post-merging.awk |\
  bcftools view --thread 16 --no-version -Ob -o "$out_bcf"
