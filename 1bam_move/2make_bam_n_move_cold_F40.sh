#!/bin/bash

# It could be better if we intake a txt containing .cram file name and coordinate biological meaning, so that he bam file name make sense


for cramfile in $(find /Volumes/Data/311/a/ -type f -maxdepth 1 -name "*.cram" -exec basename {} \;)
do
  cramfilename=${cramfile%.*}
  bamfilename=$(awk -v var="$cramfilename" '$1==var{print $2}' ./bam_name_convert.txt)
#  echo "$cramfilename"
#  echo "$bamfilename"
#  rm -r  ~/RS/F30_all/data/${bamfilename}/; mkdir ~/RS/F30_all/data/${bamfilename}/
  echo "samtools view -b /Volumes/Data/311/a/${cramfile} -T ~/RS/flight_test/data/reference/dsimM252v1.2+microbiome.fa -o ~/RS/F30_all/data/${bamfilename}/${bamfilename}.bam"
done
#samtools view -b /Volumes/Data/311/a/LB_311.TI02.cram -T ~/RS/flight_test/data/reference/dsimM252v1.2+microbiome.fa -o ~/RS/flight_test/data/F30/LB_311.TI02.bam
