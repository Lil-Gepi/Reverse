#!/bin/bash
sample=$1
id_cram=$(awk -v var="$sample" '$2==var{print $1}' ~/RS/pipeline/1cram_bcf/cram2bam_name.txt)
pool=$(echo ${id_cram:3:3})
input_cram="/Volumes/Data/${pool}/a/${id_cram}.cram"

working_dir="~/RS"
cd "${working_dir}/pipeline/"
out_bcf_dir="${working_dir}/inter_data/${sample}/"
out_bcf="${out_bcf_dir}/${sample}.bcf"
ref_fa="${working_dir}/reference/dsimM252v1.2+microbiome.fa"

rl=$(samtools view -u ${input_cram} |\
  head -n 1000000 | cut -f 10 | perl -ne 'chomp;print length($_) . "\n"' |\
  sort | uniq -c | sort -n -r | head -n1 |awk -F ' ' '{print $2}')

cov=$(printf "%.0f\n" $( samtools idxstats ${input_cram} | awk -v readlen=${rl} '
    {
        len += $2
        nreads += $3
    }
    END {
        print nreads * readlen / len
    }
')| bc -l)
max_cov=$((cov * 5))
chroms="X 2L 2R 3L 3R 4"
samtools view -u ${input_cram} $chroms |\
samtools collate -Ouf - |\
 samtools fixmate -ru - - |\
  samtools view -f 0x2 -h - |\
./rupert_pipe_v3/flag-short.awk -v READ_LENGTH=$rl |\
./rupert_pipe_v3/remapq.awk -v NORD=7.5 |\
samtools sort -@3 -m 4G -u -T /tmp/sortsam.$(basename ${out_bcf%.bcf}).$tag - |\
bcftools mpileup -f ${ref_fa} -d $max_cov --skip-all-unset 0x2 -q 10 -Q 20 -D -a DP,AD,QS,SCR -Ou --threads 3 - |\
bcftools annotate -Ov -x INFO/IMF,INFO/VDB,INFO/RPBZ,INFO/MQBZ,INFO/MQSBZ,INFO/SCBZ,INFO/BQBZ,INFO/FS,INFO/SGB,INFO/I16,INFO/QS,INFO/MQ0F,FORMAT/PL - |\
./rupert_pipe_v3/info2fmt.awk -v tags=DP |\
bcftools view --no-version -Ob -o "$out_bcf"

exit
