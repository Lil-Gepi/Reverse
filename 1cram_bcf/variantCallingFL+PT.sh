#!/bin/zsh

ref=/Volumes/Data/Dagny/Dropbox\ \(PopGen\)/DagnyPhD/Projects/PopulationSpecificAdaptation/reference/ 
refLength=$(head -4 ${ref}dsimM252v1.2+microbiome.fa.fai | awk 'BEGIN{FS="\t"; sum=0} {sum+=$2} END{print sum}')
input=/Volumes/Data/Dagny/Dropbox\ \(PopGen\)/DagnyPhD/Projects/PopulationSpecificAdaptation/data/
echo $refLength
bcfOutput=/Volumes/Temp2/preBcfHybridData/FLPT/
output=/Volumes/Temp/mergedHybridData/

#start for the folder which have multiple sequence files - sequneced on multiple lanes.
for folder in 517 522 524 525 526 #523
do
	for file in ${input}${folder}/a/*.cram
	do

		#define a file name of the new merged file and the bcf file which will be the end product of the loop
                bcfFile=${bcfOutput}${$(basename ${file})%".cram"}.bcf

		#get the read length of the file. Because the reads are not trimmed the first read can simply we checcked and saved into a variable
               	#tempReadLength=$(samtools view -T ${ref}dsimM252v1.2+microbiome.fa -F 0x400 ${file} | head -1 | awk '{ print length($10)}')
                #tempReadLength=$(samtools view -T ${ref}dsimM252v1.2+microbiome.fa -F 0x400 ${file} | awk '{print $10}' | head -100000 | python3 mostFreqReadLength.py)	
		tempReadLength=$(samtools view -T ${ref}dsimM252v1.2+microbiome.fa -F 0x400 ${file} | awk 'match($6,"H"){next}{print}' | head -1 | awk '{print length($10)}')

		#Get an estimate of the coverage by: coverage = (read count * read length ) / total genome size. First, index the merged bamfile
		samtools index ${file}
                readCount=$(samtools view -T ${ref}dsimM252v1.2+microbiome.fa -c -F 0x400 ${file} 2L 2R 3L 3R)
                coverage=$((readCount * tempReadLength))
                coverage=$((coverage / refLength))
                coveragex5=$((coverage * 5)) #max. depth to consider, set to 5x expected autosomal coverage
		
		echo "Pileing up " ${file} " using read length: " $tempReadLength " and coverage: " $coveragex5 

		#convert to qname-sorted SAM that only contains primary pairs (needed by remapq)
		samtools view -u -T ${ref}dsimM252v1.2+microbiome.fa ${file} 2L 2R 3L 3R X 4 | samtools collate -Ouf - | samtools fixmate -ru - - | samtools view -F 0x400 -f 0x2 -h - |\
		
		#flag short fragments  (replace $rl with the actual read length)
		./flag-short.awk -v READ_LENGTH=$tempReadLength |\

		#recalculate mapping quality
		./remapq.awk -v NORD=7.5 |\

		#resort by coordinate (maybe use -m 4G to use more memory and -T /Volumes/Temp/sortsam to use a larger temp dir)
		samtools sort -@1 -m 4G -u -T ${output}sortsamt.$(basename ${file}) - |\

		#make single-sample pileup
		bcftools mpileup -f ${ref}dsimM252v1.2+microbiome.fa -d $coveragex5 --skip-all-unset 0x2 -q 10 -Q 20 -D -a DP,AD,QS,SCR -Ou --threads 3 - |\
		
		#remove INFO tags we do not need to reduce file size
		bcftools annotate -Ov -x INFO/IMF,INFO/VDB,INFO/RPBZ,INFO/MQBZ,INFO/MQSBZ,INFO/SCBZ,INFO/BQBZ,INFO/FS,INFO/SGB,INFO/I16,INFO/QS,INFO/MQ0F,FORMAT/PL |\
		
		#annotate raw single-sample pileup with more FORMAT tags (specifically, save INFO/DP to FORMAT/DP)
		./info2fmt.awk -v tags=DP | bcftools view --no-version -Ob > ${bcfFile}

	done
done

#now combine multiple single-sample pilups into a multi-sample pileup
#echo "combining bcf files: \n" $(ls ${bcfOutput}LB*.bcf)  #shows the files in the same order in which you get the columns in the mpileup

#make a list for the merged pileup
#ls ${bcfOutput}LB*.bcf > ${bcfOutput}bcfList.txt

#the mpilup file which will be saved
#mbcfFile=${bcfOutput}FLPTmpileup.bcf

#merge multiple raw pileups into multi-sample mpileup
#bcftools merge --file-list ${bcfOutput}bcfList.txt --force-samples --no-index --merge both --threads 8 -Ov -o - |\
  # 2.b annotate with FORMAT/AF (observed allele frequencies) and FORMAT/XF (expected allele frequencies
  # under a multinomial sampling model) tags; also add  FORMAT/SAC (sum of allele counts) tag for convenience
  # finally, change REF of positions with a reference count of 0 to ALT{1} and modify all affected tags accordingly
  # tag such positions with the INFO/RMOD flag
#	./post-merging.awk |\
  # directly pipe into the step 4.a or
  # 2.c save VCF to compressed format
#	bcftools view --no-version -Ob > ${mbcfFile}

#bcftools reheader -s <(cut -f2 ${bcfOutput}bcfListNames.txt) --output ${bcfOutput}reheader.${$(basename ${mbcfFile})} ${mbcfFile}

#Filter mpileup to keep only SNPs
#echo "filtering SNPs from mpileup"

#BED file of regions to exclude (must match the reference)
#norepeats_bed=${ref}dsimM252v1.2_norepeats.bed

#filter expression removin non-SNPs and positions of dubious
#coverage depth (adjust as needed; here, ~2200x coverage was
#expected across all samples in the mpileup)
#averagePileupDepth=$(bcftools query -f '%DP\n' ${mbcfFile} | awk '{a+=$1}END{print a/NR}')
#echo ${averagePileupDepth}

#flt_expr="TYPE = \"snp\" & INFO/DP > $((averagePileupDepth / 2)) & INFO/DP < $((averagePileupDepth * 2))"
#echo ${flt_expr}

#output path (bgzipped vcf is the most widely used format)
#vcf=${bcfOutput}FLPT.flt.vcf.gz

# 3.a filter SNPs close to INDELs or in repeat regions
# -g 5       ... filter out SNPs that are this close to INDELs
# -T bedfile ... only output positions described by the provided BED file
#bcftools filter -g 5 -T "$norepeats_bed" -Ou ${mbcfFile} |\
  # 3.b keep only SNPs and adjust coverage depth thresholds as appropriate
  # (here ~2200x was expected across all samples)
#  bcftools filter -i "$flt_expr" -Ou - | \
  # 3.c make positions single-allelic for easier manipulation
#  bcftools norm -m- -Oz - > ${vcf}


# Some examples of what to do with the mpileups

# Example: extract coverages & allele frequencies as needed for CMH test, e.g., using ACER
#bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAD]\n' ${vcf} | bgzip > ${bcfOutput}FLPT.flt.cov.txt.gz
#bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%XF]\n' ${vcf} | bgzip > ${bcfOutput}FLPT.flt.alf.txt.gz





