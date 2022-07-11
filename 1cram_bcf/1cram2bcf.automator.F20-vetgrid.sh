#!/bin/bash
for (( i = 11; i < 15; i++ )); do
  for (( j = 1; j < 6; j++ )); do
    sample=F20r${i}sub${j}
    echo $sample
  done
done | parallel -j 2 "time bash /Users/ychen/RS/pipeline/1cram_bcf/1cram2bcf.operator.sh {}"
