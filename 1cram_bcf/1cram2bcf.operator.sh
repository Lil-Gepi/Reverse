#!/bin/bash
working_dir="/Users/ychen/RS"
cd "${working_dir}/pipeline/"

sample=$1
id_cram=$(awk -v var="$sample" '$2==var{print $1}' "${working_dir}/pipeline/1cram_bcf/cram2bam_name.txt")
pool=$(echo ${id_cram:3:3})
input_cram="/Volumes/Data/${pool}/a/${id_cram}.cram"
samtools index ${input_cram}
out_bcf_dir="${working_dir}/inter_data/${sample}/"
out_bcf="${out_bcf_dir}/${sample}.bcf"
ref_fa="${working_dir}/reference/dsimM252v1.2+microbiome.fa"

echo "now dealing with ${sample}, the data is from pool ${pool}, its id is ${id_cram}"

rl=$(samtools view -T ${ref_fa} -h -F 0x400 ${input_cram} |\
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
echo "the read length of ${sample} is ${rl}, and average coverage of ${cov}"
max_cov=$((cov * 5))
chroms="X 2L 2R 3L 3R 4"
samtools view -T ${ref_fa} -h -F 0x400 ${input_cram} $chroms |\
samtools collate -Ouf - |\
samtools fixmate -ru - - |\
samtools view -f 0x2 -h - |\
./rupert_pipe_v3/flag-short.awk -v READ_LENGTH=$rl |\
./rupert_pipe_v3/remapq.awk -v NORD=7.5 |\
samtools sort -@3 -m 4G -T /tmp/sortsam.${sample} -u - |\
bcftools mpileup -f ${ref_fa} -d $max_cov --skip-all-unset 0x2 -q 10 -Q 20 -D -a DP,AD,QS,SCR -Ou --threads 3 - |\
bcftools annotate -Ov -x INFO/IMF,INFO/VDB,INFO/RPBZ,INFO/MQBZ,INFO/MQSBZ,INFO/SCBZ,INFO/BQBZ,INFO/FS,INFO/SGB,INFO/I16,INFO/QS,INFO/MQ0F,FORMAT/PL |\
./rupert_pipe_v3/info2fmt.awk -v tags=DP |\
bcftools view --no-version -Ob -o "$out_bcf"
echo "done with ${sample}"
exit
