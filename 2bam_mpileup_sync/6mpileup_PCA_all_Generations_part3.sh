#!/bin/sh

# do the job from bam to cram to sync for each subreplicate
# bamfile are like F30r11sub1 etc.
bamfile_id_list=(coldF40r11 coldF40r12 coldF40r13 coldF40r14 \
F10r11sub1 F10r11sub2 F10r12sub1 F10r12sub2 F10r13sub1 F10r13sub2 F10r14sub1 F10r14sub2 \
F20r11sub1 F20r11sub2 F20r12sub1 F20r12sub2 F20r13sub1 F20r13sub2 F20r14sub1 F20r14sub2 \
F30r11sub1 F30r11sub2 F30r11sub3 F30r11sub4 F30r11sub5 F30r12sub1 F30r12sub2 F30r12sub3 F30r12sub4 F30r12sub5 F30r13sub1 F30r13sub2 F30r13sub3 F30r13sub4 F30r13sub5 F30r14sub1 F30r14sub2 F30r14sub3 F30r14sub4 F30r14sub5 \
F40r11sub1 F40r11sub2 F40r12sub1 F40r12sub2 F40r13sub1 F40r13sub2 F40r14sub1 F40r14sub2 \
F50r11sub1 F50r11sub2 F50r12sub1 F50r12sub2 F50r13sub1 F50r13sub2 F50r14sub1 F50r14sub2 \
F60r11sub1 F60r11sub2 F60r12sub1 F60r12sub2 \
F80r11sub1 F80r11sub2 F80r12sub1 F80r12sub2 F80r13sub1 F80r13sub2 F80r14sub1 F80r14sub2)
for bamfile in ${bamfile_id_list}
do
  if test -f "~/RS/data/${bamfile}/${bamfile}.bai"
  then
    echo "${bamfile} already indexed."
  else
    echo "${bamfile} is now being indexed."
    samtools index -b ~/RS/data/${bamfile}/${bamfile}.bam
  fi
done

for bamfile in ${bamfile_id_list}
do
  if test -f "~/RS/data/${bamfile}/${bamfile}_2L.bam"
  then
    echo "${bamfile} already seperated by chromosome."
  else
    for chr in 2L 2R 3L 3R 4 X
    do
      echo "seperating ${bamfile} for chromosome ${chr}."
      samtools view -b ~/RS/data/${bamfile}/${bamfile}.bam ${chr} > \
      ~/RS/data/${bamfile}/${bamfile}_${chr}.bam
      samtools index ~/RS/data/${bamfile}/${bamfile}_${chr}.bam
    done
  fi
done

for chr in 2L 2R 3L 3R 4 X
do
  samtools mpileup -B -q 20 -Q 0 -f ~/RS/reference/dsimM252v1.2+microbiome.fa \
  $(ls -R /Users/ychen/RS/data/*/*_${chr}.bam) | \
  # need to figure out how to specifz the Pop and Gen instead of * all bam files.
  # maybe partition the gen? and R? so that I can put them into the code.
  java -jar ~/RS/pipeline/popoolation2/mpileup2sync.jar --input /dev/stdin \
  --output ~/RS/result/reverse_${chr}.sync --fastq-type sanger \
  --min-qual 20 --threads 12

  samtools mpileup -B -q 20 -Q 0 -f ~/RS/reference/dsimM252v1.2+microbiome.fa \
  ~/RS/data/coldF40r11/coldF40r11_${chr}.bam ~/RS/data/coldF40r11/coldF40r12_${chr}.bam \
  ~/RS/data/coldF40r11/coldF40r13_${chr}.bam ~/RS/data/coldF40r11/coldF40r14_${chr}.bam \
  ~/RS/data/F30r11sub1/F10r11sub1_${chr}.bam ~/RS/data/F30r11sub2/F10r11sub2_${chr}.bam \
  ~/RS/data/F30r12sub1/F10r12sub1_${chr}.bam ~/RS/data/F30r12sub2/F10r12sub2_${chr}.bam \

  ~/RS/data/F30r11sub1/F30r11sub1_${chr}.bam \
  ~/RS/data/F30r11sub2/F30r11sub2_${chr}.bam ~/RS/data/F30r11sub3/F30r11sub3_${chr}.bam \
  ~/RS/data/F30r11sub4/F30r11sub4_${chr}.bam ~/RS/data/F30r11sub5/F30r11sub5_${chr}.bam | \
  # need to figure out how to specifz the Pop and Gen instead of * all bam files.
  # maybe partition the gen? and R? so that I can put them into the code.
  java -jar ~/RS/pipeline/popoolation2/mpileup2sync.jar --input /dev/stdin \
  --output ~/RS/result/cold_F30_cmh/cold_F30_r11_${chr}.sync --fastq-type sanger \
  --min-qual 20 --threads 12

done


cat \
~/RS/result/cold_F30_cmh/cold_F30_r11_2L.sync ~/RS/result/cold_F30_cmh/cold_F30_r11_2R.sync \
~/RS/result/cold_F30_cmh/cold_F30_r11_3L.sync ~/RS/result/cold_F30_cmh/cold_F30_r11_3R.sync \
~/RS/result/cold_F30_cmh/cold_F30_r11_4.sync ~/RS/result/cold_F30_cmh/cold_F30_r11_X.sync > \
~/RS/result/cold_F30_cmh/cold_F30_r11_all_chr_mq20_bq20.sync
