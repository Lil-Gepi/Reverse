#!/bin/sh

# index the reference genome fasta
samtools faidx ~/RS/F30_all/reference/dsimM252v1.2+microbiome.fa
# do the job from bam to cram to sync for each subreplicate
# bamfile are like F30r11sub1 etc.
for chr in 2L 2R 3L 3R 4 X
do
  for bamfile in $(ls ~/RS/F30_all/data/)
  do
    cd ~/RS/F30_all/data/${bamfile}/
    # for samtools version 1.9. Index Bam file (create bai file) before split.
    # samtools index -b ~/RS/F30_all/data/${bamfile}/${bamfile}.bam
    # split BAM into chromosomes of interest and index them,
    # and we do the job on each chromosome and later concate them
    # can also extract information for mitochondrial DNA which name is mtDNA_65039429 something
    samtools view -b ~/RS/F30_all/data/${bamfile}/${bamfile}.bam ${chr} > \
    ~/RS/F30_all/data/${bamfile}/${bamfile}_${chr}.bam
    samtools index ~/RS/F30_all/data/${bamfile}/${bamfile}_${chr}.bam
  done
  samtools mpileup -B -q 20 -Q 0 -f ~/RS/F30_all/reference/dsimM252v1.2+microbiome.fa \
  $(ls -R /Users/ychen/RS/F30_all/data/*/*_${chr}.bam) | \
  java -jar ~/RS/F30_all/pipeline/popoolation2/mpileup2sync.jar --input /dev/stdin \
  --output ~/RS/F30_all/result/F30_r11_r12_${chr}.sync --fastq-type sanger \
  --min-qual 20 --threads 8

done


cat \
~/RS/F30_all/result/F30_r11_r12_2L.sync ~/RS/F30_all/result/F30_r11_r12_2R.sync \
~/RS/F30_all/result/F30_r11_r12_3L.sync ~/RS/F30_all/result/F30_r11_r12_3R.sync \
~/RS/F30_all/result/F30_r11_r12_4.sync ~/RS/F30_all/result/F30_r11_r12_X.sync > \
~/RS/F30_all/result/F30_r11_r12_all_chr_mq20_bq20.sync




# old codes :
# samtools view -b /Users/ychen/RS/flight_test/data/F30/LB_311.TI02.bam 2L > ../data/F30/Pool_311.F30r11sub1_2L.bam &&\
# samtools view -b /Users/ychen/RS/flight_test/data/F30/LB_311.TI02.bam 2R > ../data/F30/Pool_311.F30r11sub1_2R.bam &&\
# samtools view -b /Users/ychen/RS/flight_test/data/F30/LB_311.TI02.bam 3L > ../data/F30/Pool_311.F30r11sub1_3L.bam &&\
# samtools view -b /Users/ychen/RS/flight_test/data/F30/LB_311.TI02.bam 3R > ../data/F30/Pool_311.F30r11sub1_3R.bam &&\
# samtools view -b /Users/ychen/RS/flight_test/data/F30/LB_311.TI02.bam 4 > ../data/F30/Pool_311.F30r11sub1_4.bam &&\
# samtools view -b /Users/ychen/RS/flight_test/data/F30/LB_311.TI02.bam X > ../data/F30/Pool_311.F30r11sub1_X.bam &&\
#
# samtools index ../data/F30/Pool_311.F30r11sub1_2L.bam &&\
# samtools index ../data/F30/Pool_311.F30r11sub1_2R.bam &&\
# samtools index ../data/F30/Pool_311.F30r11sub1_3L.bam &&\
# samtools index ../data/F30/Pool_311.F30r11sub1_3R.bam &&\
# samtools index ../data/F30/Pool_311.F30r11sub1_4.bam &&\
# samtools index ../data/F30/Pool_311.F30r11sub1_X.bam &&\
#
# samtools mpileup -B -q 20 -Q 0 -f ~/RS/flight_test/data/reference/dsimM252v1.2+microbiome.fa ../data/F30/Pool_311.F30r11sub1_2L.bam > ../data/F30/Pool_311.F30r11sub1_2L.pileup &&\
# samtools mpileup -B -q 20 -Q 0 -f ~/RS/flight_test/data/reference/dsimM252v1.2+microbiome.fa ../data/F30/Pool_311.F30r11sub1_2R.bam > ../data/F30/Pool_311.F30r11sub1_2R.pileup &&\
# samtools mpileup -B -q 20 -Q 0 -f ~/RS/flight_test/data/reference/dsimM252v1.2+microbiome.fa ../data/F30/Pool_311.F30r11sub1_3L.bam > ../data/F30/Pool_311.F30r11sub1_3L.pileup &&\
# samtools mpileup -B -q 20 -Q 0 -f ~/RS/flight_test/data/reference/dsimM252v1.2+microbiome.fa ../data/F30/Pool_311.F30r11sub1_3R.bam > ../data/F30/Pool_311.F30r11sub1_3R.pileup &&\
# samtools mpileup -B -q 20 -Q 0 -f ~/RS/flight_test/data/reference/dsimM252v1.2+microbiome.fa ../data/F30/Pool_311.F30r11sub1_4.bam > ../data/F30/Pool_311.F30r11sub1_4.pileup &&\
# samtools mpileup -B -q 20 -Q 0 -f ~/RS/flight_test/data/reference/dsimM252v1.2+microbiome.fa ../data/F30/Pool_311.F30r11sub1_X.bam > ../data/F30/Pool_311.F30r11sub1_X.pileup &&\
#
# # concate all the files into one
# cd ../data/F30/
# cat Pool_311.F30r11sub1_2L.pileup Pool_311.F30r11sub1_2R.pileup Pool_311.F30r11sub1_3L.pileup Pool_311.F30r11sub1_3R.pileup Pool_311.F30r11sub1_4.pileup Pool_311.F30r11sub1_X.pileup > Pool_311.F30r11sub1_allChr.pileup
