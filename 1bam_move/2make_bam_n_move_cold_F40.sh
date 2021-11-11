


for cramfile in $(find /Volumes/Data/185/*/ -type f -maxdepth 1 -name "*.cram" -exec basename {} \;)
do
  cramfilename=${cramfile%.*}
  bamfilename=$(awk -v var="$cramfilename" '$1==var{print $2}' ~/RS/pipeline/1bam_move/cram2bam_name.txt)
#  echo "$cramfilename"
#  echo "$bamfilename"
  # if [ ! -d ~/RS/F30_all/data/${bamfilename}/ ]; then
  #   mkdir ~/RS/F30_all/data/${bamfilename}/
  #   samtools view -b -q 20 -F 0x400 -T ~/RS/flight_test/data/reference/dsimM252v1.2+microbiome.fa -o ~/RS/F30_all/data/${bamfilename}/${bamfilename}.bam /Volumes/Data/311/a/${cramfile}
  # fi
rm -r  ~/RS/data/${bamfilename}/; mkdir ~/RS/data/${bamfilename}/
samtools view -b -q 20 -F 0x400 -T ~/RS/reference/dsimM252v1.2+microbiome.fa \
-o ~/RS/data/${bamfilename}/${bamfilename}.bam /Volumes/Data/${pool}/a/${cramfile}
done | parallel -j 10
