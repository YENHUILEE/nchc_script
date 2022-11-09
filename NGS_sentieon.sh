#!/bin/bash
#SBATCH -p ngs48G
#SBATCH -c 14
#SBATCH --mem=46g
#SBATCH -A MST109178
#SBATCH -J MG244
#SBATCH -e /work/u3003390/err 
#SBATCH -o /work/u3003390/out
#SBATCH --mail-user=k32650805@gmail.com
#SBATCH --mail-type=FAIL,END

# This code is used for NGS germline variant calling (SNV + Indels)

echo "processes directories"
# before start, you should copy the fastq to fastq_folder
# fold to be set
sample_id="MG244" #change as you needed
user_dir="/work/u3003390"
fastq_dir="/work/u3003390/FASTQ" #data fold
temp_dir="/work/u3003390/TEMP"
release_dir="/work/u3003390/RESULT"
set -euo pipefail

# Update with the fullpath location of your sample fastq
fastq_1="${fastq_dir}/${sample_id}*R1*.gz" #NGS1_20170305A.R1.fastq.gz
fastq_2="${fastq_dir}/${sample_id}*R2*.gz"  #If using Illumina paired data
sample="SM_"${sample_id}
group="GP_"${sample_id}
platform="ILLUMINA"


# Update with the location of the reference data files (hg 19)
# May use hg 38 if WGS/WES 
ref_dir="/staging/reserve/paylong_ntu/AI_SHARE/reference/GATK_bundle/2.8/hg19" 
fasta="${ref_dir}/ucsc.hg19.fasta"
dbsnp="${ref_dir}/dbsnp_138.hg19.vcf"
known_Mills_indels="${ref_dir}/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"
known_1000G_indels="${ref_dir}/1000G_phase1.indels.hg19.sites.vcf"

# Other settings
nt=40 #number of threads to use in computation

# Set working directory
Time=`date +%Y%m%d%H%M`
Date=`date +%Y%m%d`
logfile=${temp_dir}/${Time}_run.log
set -x
exec 3<&1 4<&2 #???
exec >$logfile 2>&1

#module loading#
module load biology/BWA/0.7.17
module load biology/Picard/2.27.4
ml load biology/GATK/4.2.3.0
sentieon = "/staging/reserve/paylong_ntu/AI_SHARE/software/Sentieon/sentieon-genomics-201808/bin/sentieon driver"
module load biology/ANNOVAR/2020-06-08
module load biology/Samtools/1.15.1

# ******************************************
# 0. Setup
# ******************************************
mkdir -p $workdir
DATE=`date +%Y%m%d%H%M%S`
logfile=$workdir/${DATE}_run.log
set -x
exec 3<&1 4<&2
exec >$logfile 2>&1

cd $workdir


# ******************************************
# 1. Mapping reads with BWA-MEM, sorting
# ******************************************

##fasta index file
#if [ ! -d "${ref_dir}/ucsc.hg19.fasta.bwt" ]; then
#$release_dir/bin/bwa index $fasta
#fi

(bwa mem -M -R "@RG\tID:$group\tSM:$sample\tPL:$platform" \
-t $nt -K 10000000 $fasta $fastq_1 $fastq_2 || echo -n 'error' ) | \
$sentieon util sort \
-r $fasta \
-o ${SampleName}.sorted.bam \
-t $nt --sam2bam -i-

### To output unmapped.bam
samtools view -b -f 4 ${SampleName}.sorted.bam > ${SampleName}.sorted.unmapped.bam

### To output in cram format
samtools view -C -T $fasta ${SampleName}.sorted.unmapped.bam > ${SampleName}.sorted.unmapped.cram
samtools index ${SampleName}.sorted.unmapped.cram
samtools view -C -T $fasta ${SampleName}.sorted.bam > ${SampleName}.sorted.cram
samtools index ${SampleName}.sorted.cram

# ******************************************
# 2. Metrics
# ******************************************
${sentieon} driver \
-r ${fasta} \
-i sorted.bam \
-t $nt \
--algo MeanQualityByCycle mq_metrics.txt \
--algo QualDistribution qd_metrics.txt \
--algo GCBias \
--summary gc_summary.txt gc_metrics.txt \
--algo AlignmentStat \
--adapter_seq '' aln_metrics.txt \
--algo InsertSizeMetricAlgo is_metrics.txt

${sentieon} plot metrics \
-o metrics-report.pdf gc=gc_metrics.txt qd=qd_metrics.txt mq=mq_metrics.txt isize=is_metrics.txt


# ******************************************
# 3. Remove Duplicate Reads
# ******************************************
${sentieon} driver  \
-t $nt -i ${SampleName}.sorted.bam \
--algo LocusCollector \
--fun score_info ${SampleName}.score.txt
${sentieon} driver  \
-t $nt \
-i ${SampleName}.sorted.bam \
--algo Dedup \
--rmdup \
--score_info ${SampleName}.score.txt \
--metrics ${SampleName}.dedup_metrics.txt ${SampleName}.deduped.bam


# ******************************************
# 4. Indel realigner
# ******************************************
${sentieon} driver \
-t $nt \
-r $fasta  \
-i ${SampleName}.deduped.bam \
--algo Realigner \
-k $known_Mills_indels \
-k $known_1000G_indels \
${SampleName}.realigned.bam


# ******************************************
# 5. Base recalibration
# ******************************************
${sentieon} driver\
-r $fasta \
-t $nt \
-i ${SampleName}.realigned.bam \
--algo QualCal \
-k $dbsnp -k $known_Mills_indels -k $known_1000G_indels \
${SampleName}.recal_data.table

# ${sentieon} driver \
# -r $fasta \
# -t $nt \
# -i realigned.bam \
# -q recal_data.table \
# --algo QualCal \
# -k $dbsnp -k $known_Mills_indels -k $known_1000G_indels \
# recal_data.table.post

# ${sentieon} driver \
# -t $nt \
# --algo QualCal \
# --plot --before recal_data.table \
# --after recal_data.table.post recal.csv

# ${sentieon} driver \
# plot bqsr \
# -o recal_plots.pdf recal.csv


# ******************************************
# 6a. UG Variant caller
# ******************************************
#$release_dir/bin/sentieon driver -r $fasta -t $nt -i realigned.bam -q recal_data.table --algo Genotyper -d $dbsnp --var_type BOTH --emit_conf=10 --call_conf=30 output-ug.vcf.gz


# ******************************************
# 6b. HC Variant caller
# ******************************************
${sentieon} driver \
-r $fasta \
-t $nt \
-i ${SampleName}.realigned.bam \
-q ${SampleName}.recal_data.table \
--algo Haplotyper \
-d $dbsnp \
--emit_conf=10 \
--call_conf=30 ${SampleName}.output-hc.vcf.gz

# gvcf
${sentieon} driver \
-r $fasta \
-t $nt \
-i ${SampleName}.realigned.bam \
-q ${SampleName}.recal_data.table \
--algo Haplotyper \
-d $dbsnp \
--emit_mode gvcf \
${SampleName}.output-hc.g.vcf.gz

# ******************************************
# 5b. ReadWriter to output recalibrated bam
# This stage is optional as variant callers
# can perform the recalibration on the fly
# using the before recalibration bam plus
# the recalibration table
# ******************************************
${sentieon} driver \
-r $fasta \
-t $nt \
-i ${SampleName}.realigned.bam \
-q ${SampleName}.recal_data.table \
--algo ReadWriter ${SampleName}.recaled.bam

### To get collect QC HsMetrics
module load biology/Picard/2.27.4
BED="/staging/reserve/paylong_ntu/AI_SHARE/GitHub/Germline_variant/A1_Panel/Deafness_bed/20210304.hg19.interval_list"
picard CollectHsMetrics \
INPUT=${SampleName}.recaled.bam \
OUTPUT=${SampleName}.hs_metrics.txt \
R=$fasta BAIT_INTERVALS=${BED} \
TARGET_INTERVALS=${BED}

### To output in cram format
samtools view -C -T $fasta ${SampleName}.recaled.bam > ${SampleName}.recaled.cram
samtools index ${SampleName}.recaled.cram


if [[ $? -eq 0 ]]; then
	ls ${SampleName}.sorted.unmapped.bam ${SampleName}.sorted.bam ${SampleName}.deduped.bam ${SampleName}.realigned.bam ${SampleName}.recaled.bam
	rm                                   ${SampleName}.sorted.bam ${SampleName}.deduped.bam ${SampleName}.realigned.bam 
fi

set +x
exec >&3 2>&4
exec 3<&- 4<&-

printf "#############################################################################\n"
printf "###                  Work completed: $(date +%Y-%m-%d:%H:%M:%S)                  ###\n"
printf "#############################################################################\n"
