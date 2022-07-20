#!/bin/bash

in_folder="/Users/ychen/RS/inter_data"
out_bcf="/Users/ychen/RS/result/piled_bcf/RS_piled.bcf"

bcftools merge --force-samples --no-index --merge both --threads 14 -Ov -o - ${in_folder}/*/*.bcf |\
  /Users/ychen/RS/pipeline/rupert_pipe_v3/post-merging.awk |\
  bcftools view --thread 14 --no-version -Ob -o "$out_bcf"
