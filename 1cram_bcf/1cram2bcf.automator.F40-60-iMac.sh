#!/bin/bash
for (( i = 11; i < 15; i++ )); do
  for (( j = 1; j < 6; j++ )); do
      sample=F40r${i}sub${j}
      echo $sample
  done
done | parallel -j 5 "nohup bash ./1cram2bcf.operator.sh {}"

for (( i = 11; i < 15; i++ )); do
  for (( j = 1; j < 6; j++ )); do
      sample=F50r${i}sub${j}
      echo $sample
  done
done | parallel -j 5 "nohup bash ./1cram2bcf.operator.sh {}"

for (( i = 11; i < 15; i++ )); do
  for (( j = 1; j < 6; j++ )); do
      sample=F60r${i}sub${j}
      echo $sample
  done
done | parallel -j 5 "nohup bash ./1cram2bcf.operator.sh {}"
