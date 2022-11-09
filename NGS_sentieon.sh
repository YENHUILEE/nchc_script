#!/usr/bin/sh
#SBATCH -A MST109178        # Account name/project number
#SBATCH -J SAMPLE_NAME         # Job name
#SBATCH -p ngs48G           # Partition Name 等同PBS裡面的 -q Queue name
#SBATCH -c 14               # 使用的core數 請參考Queue資源設定 
#SBATCH --mem=46g           # 使用的記憶體量 請參考Queue資源設定
#SBATCH -o out.log          # Path to the standard output file 
#SBATCH -e err.log          # Path to the standard error ouput file
#SBATCH --mail-user=@gmail.com    # email
#SBATCH --mail-type=FAIL              # 指定送出email時機 可為NONE, BEGIN, END, FAIL, REQUEUE, ALL


JOBDIR="/staging/reserve/paylong_ntu/AI_SHARE/GitHub/Germline_variant/Sentieon/"

SampleName=SAMPLE_NAME

cd $JOBDIR

set -euo pipefail

printf "#############################################################################\n"
printf "###                  Work started:   $(date +%Y-%m-%d:%H:%M:%S)                  ###\n"
printf "#############################################################################\n"


# *******************************************
# Script to perform DNA seq variant calling
# using a single sample with fastq files
# named 1.fastq.gz and 2.fastq.gz
# *******************************************

# Update with the fullpath location of your sample fastq

fastq_folder="/staging/reserve/paylong_ntu/AI_SHARE/rawdata/Deafness/"
fastq_1="${fastq_folder}/${SampleName}/panel/${SampleName}*R1*.gz" #NGS1_20170305A.R1.fastq.gz
fastq_2="${fastq_folder}/${SampleName}/panel/${SampleName}*R2*.gz"  #If using Illumina paired data
sample="SM_"${SampleName}
group="GP_"${SampleName}
platform="ILLUMINA"


# Update with the location of the reference data files
ref_dir="/staging/reserve/paylong_ntu/AI_SHARE/reference/GATK_bundle/2.8/hg19"

fasta="${ref_dir}/ucsc.hg19.fasta"
dbsnp="${ref_dir}/dbsnp_138.hg19.vcf"
known_Mills_indels="${ref_dir}/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"
known_1000G_indels="${ref_dir}/1000G_phase1.indels.hg19.sites.vcf"

# Determine whether Variant Quality Score Recalibration will be run
# VQSR should only be run when there are sufficient variants called
#run_vqsr="yes"

# Update with the location of the resource files for VQSR
vqsr_Mill="${ref_dir}/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"
vqsr_1000G_omni="${ref_dir}/1000G_omni2.5.hg19.sites.vcf"
vqsr_hapmap="${ref_dir}/hapmap_3.3.hg19.sites.vcf"
vqsr_1000G_phase1="${ref_dir}/1000G_phase1.snps.high_confidence.hg19.sites.vcf"
vqsr_1000G_phase1_indel="${ref_dir}/1000G_phase1.indels.hg19.sites.vcf"
vqsr_dbsnp="${ref_dir}/dbsnp_138.hg19.vcf"

# Update with the location of the Sentieon software package and license file
export SENTIEON_LICENSE=140.110.16.119:8990
release_dir="/staging/reserve/paylong_ntu/AI_SHARE/software/Sentieon/sentieon-genomics-201808"

# Other settings
nt=40 #number of threads to use in computation
workdir=${JOBDIR}/${SampleName} #Determine where the output files will be stored

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
#The results of this call are dependent on the number of threads used. To have number of threads independent results, add chunk size option -K 10000000


#if [ ! -d "${ref_dir}/ucsc.hg19.fasta.bwt" ]; then
#$release_dir/bin/bwa index $fasta
#fi

($release_dir/bin/bwa mem -M -R "@RG\tID:$group\tSM:$sample\tPL:$platform" -t $nt -K 10000000 $fasta $fastq_1 $fastq_2 || echo -n 'error' ) | $release_dir/bin/sentieon util sort -r $fasta -o ${SampleName}.sorted.bam -t $nt --sam2bam -i-

### To output unmapped.bam
/opt/ohpc/Taiwania3/pkg/biology/SAMTOOLS/SAMTOOLS_v1.13/bin/samtools view -b -f 4 ${SampleName}.sorted.bam > ${SampleName}.sorted.unmapped.bam

### To output in cram format
/opt/ohpc/Taiwania3/pkg/biology/SAMTOOLS/SAMTOOLS_v1.13/bin/samtools view -C -T $fasta ${SampleName}.sorted.unmapped.bam > ${SampleName}.sorted.unmapped.cram
/opt/ohpc/Taiwania3/pkg/biology/SAMTOOLS/SAMTOOLS_v1.13/bin/samtools index ${SampleName}.sorted.unmapped.cram
/opt/ohpc/Taiwania3/pkg/biology/SAMTOOLS/SAMTOOLS_v1.13/bin/samtools view -C -T $fasta ${SampleName}.sorted.bam > ${SampleName}.sorted.cram
/opt/ohpc/Taiwania3/pkg/biology/SAMTOOLS/SAMTOOLS_v1.13/bin/samtools index ${SampleName}.sorted.cram

# ******************************************
# 2. Metrics
# ******************************************
#$release_dir/bin/sentieon driver -r $fasta -t $nt -i sorted.bam --algo MeanQualityByCycle mq_metrics.txt --algo QualDistribution qd_metrics.txt --algo GCBias --summary gc_summary.txt gc_metrics.txt --algo AlignmentStat --adapter_seq '' aln_metrics.txt --algo InsertSizeMetricAlgo is_metrics.txt
#$release_dir/bin/sentieon plot metrics -o metrics-report.pdf gc=gc_metrics.txt qd=qd_metrics.txt mq=mq_metrics.txt isize=is_metrics.txt


# ******************************************
# 3. Remove Duplicate Reads
# ******************************************
$release_dir/bin/sentieon driver  -t $nt -i ${SampleName}.sorted.bam --algo LocusCollector --fun score_info ${SampleName}.score.txt
$release_dir/bin/sentieon driver  -t $nt -i ${SampleName}.sorted.bam --algo Dedup --rmdup --score_info ${SampleName}.score.txt --metrics ${SampleName}.dedup_metrics.txt ${SampleName}.deduped.bam


# ******************************************
# 4. Indel realigner
# ******************************************
$release_dir/bin/sentieon driver -r $fasta  -t $nt -i ${SampleName}.deduped.bam --algo Realigner -k $known_Mills_indels -k $known_1000G_indels ${SampleName}.realigned.bam


# ******************************************
# 5. Base recalibration
# ******************************************
$release_dir/bin/sentieon driver -r $fasta -t $nt -i ${SampleName}.realigned.bam --algo QualCal -k $dbsnp -k $known_Mills_indels -k $known_1000G_indels ${SampleName}.recal_data.table
#$release_dir/bin/sentieon driver -r $fasta -t $nt -i realigned.bam -q recal_data.table --algo QualCal -k $dbsnp -k $known_Mills_indels -k $known_1000G_indels recal_data.table.post
#$release_dir/bin/sentieon driver -t $nt --algo QualCal --plot --before recal_data.table --after recal_data.table.post recal.csv
#$release_dir/bin/sentieon plot bqsr -o recal_plots.pdf recal.csv


# ******************************************
# 6a. UG Variant caller
# ******************************************
#$release_dir/bin/sentieon driver -r $fasta -t $nt -i realigned.bam -q recal_data.table --algo Genotyper -d $dbsnp --var_type BOTH --emit_conf=10 --call_conf=30 output-ug.vcf.gz


# ******************************************
# 6b. HC Variant caller
# ******************************************
$release_dir/bin/sentieon driver -r $fasta -t $nt -i ${SampleName}.realigned.bam -q ${SampleName}.recal_data.table --algo Haplotyper -d $dbsnp --emit_conf=10 --call_conf=30 ${SampleName}.output-hc.vcf.gz

# gvcf
$release_dir/bin/sentieon driver -r $fasta -t $nt -i ${SampleName}.realigned.bam -q ${SampleName}.recal_data.table --algo Haplotyper -d $dbsnp --emit_mode gvcf ${SampleName}.output-hc.g.vcf.gz

# ******************************************
# 5b. ReadWriter to output recalibrated bam
# This stage is optional as variant callers
# can perform the recalibration on the fly
# using the before recalibration bam plus
# the recalibration table
# ******************************************
$release_dir/bin/sentieon driver -r $fasta -t $nt -i ${SampleName}.realigned.bam -q ${SampleName}.recal_data.table --algo ReadWriter ${SampleName}.recaled.bam

### To get collect QC HsMetrics
PICARD=/opt/ohpc/Taiwania3/pkg/biology/Picard/picard_v2.26.0/picard.jar
BED=/staging/reserve/paylong_ntu/AI_SHARE/GitHub/Germline_variant/A1_Panel/Deafness_bed/20210304.hg19.interval_list
java -Xmx40g -jar ${PICARD} CollectHsMetrics INPUT=${SampleName}.recaled.bam OUTPUT=${SampleName}.hs_metrics.txt R=$fasta BAIT_INTERVALS=${BED} TARGET_INTERVALS=${BED}

### To output in cram format
/opt/ohpc/Taiwania3/pkg/biology/SAMTOOLS/SAMTOOLS_v1.13/bin/samtools view -C -T $fasta ${SampleName}.recaled.bam > ${SampleName}.recaled.cram
/opt/ohpc/Taiwania3/pkg/biology/SAMTOOLS/SAMTOOLS_v1.13/bin/samtools index ${SampleName}.recaled.cram


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
