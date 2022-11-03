#!/bin/bash
#!/usr/bin/sh
#SBATCH -A TRI1113046        # Account name/project number
#SBATCH -J SAMPLE_NAME         # Job name
#SBATCH -p ngs48G           # Partition Name 等同PBS裡面的 -q Queue name
#SBATCH -c 14               # 使用的core數 請參考Queue資源設定 
#SBATCH --mem=46g           # 使用的記憶體量 請參考Queue資源設定
#SBATCH -o out.log          # Path to the standard output file 
#SBATCH -e err.log          # Path to the standard error ouput file
#SBATCH --mail-user=k32650805@gmail.com    # email
#SBATCH --mail-type=FAIL              # 指定送出email時機 可為NONE, BEGIN, END, FAIL, REQUEUE, ALL


echo "processes directories"
# fold to be set
fastq_folder="/work/u3003390/FASTQ/" #data fold
JOBDIR="/work/u3003390/"
SampleName="MG244"
release_dir="/work/u3003390/RESUL/"

cd $JOBDIR
set -euo pipefail

# Update with the fullpath location of your sample fastq
fastq_1="${fastq_folder}/${SampleName}*R1*.gz" #NGS1_20170305A.R1.fastq.gz
fastq_2="${fastq_folder}/${SampleName}*R2*.gz"  #If using Illumina paired data
sample="SM_"${SampleName}
group="GP_"${SampleName}
platform="ILLUMINA"


# Update with the location of the reference data files (hg 19)
ref_dir="/staging/reserve/paylong_ntu/AI_SHARE/reference/GATK_bundle/2.8/hg19"
fasta="${ref_dir}/ucsc.hg19.fasta"
dbsnp="${ref_dir}/dbsnp_138.hg19.vcf"
known_Mills_indels="${ref_dir}/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"
known_1000G_indels="${ref_dir}/1000G_phase1.indels.hg19.sites.vcf"

# Determine whether Variant Quality Score Recalibration will be run
# VQSR should only be run when there are sufficient variants called
# run_vqsr="yes"
# Update with the location of the resource files for VQSR
vqsr_Mill="${ref_dir}/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"
vqsr_1000G_omni="${ref_dir}/1000G_omni2.5.hg19.sites.vcf"
vqsr_hapmap="${ref_dir}/hapmap_3.3.hg19.sites.vcf"
vqsr_1000G_phase1="${ref_dir}/1000G_phase1.snps.high_confidence.hg19.sites.vcf"
vqsr_1000G_phase1_indel="${ref_dir}/1000G_phase1.indels.hg19.sites.vcf"
vqsr_dbsnp="${ref_dir}/dbsnp_138.hg19.vcf"

# Update with the location of the Sentieon software package and license file
export SENTIEON_LICENSE=140.110.16.119:8990

# Other settings
nt=40 #number of threads to use in computation

# Set working directory
workdir=${JOBDIR}/${SampleName} #Determine where the output files will be stored
mkdir -p $workdir
DATE=`date +%Y%m%d%H%M%S`
logfile=$workdir/${DATE}_run.log
set -x
exec 3<&1 4<&2 #???
exec >$logfile 2>&1

cd $workdir


# -------------------
# STEP 1: QC - Run fastqc 
# -------------------

echo "STEP 1: QC - Run fastqc"

fastqc ${fastq_1} -o ${release_dir}/
fastqc ${fastq_2} -o ${release_dir}/

# No trimming required, quality looks okay.


# # --------------------------------------
# # STEP 2: Map to reference using BWA-MEM
# # --------------------------------------

# echo "STEP 2: Map to reference using BWA-MEM"

# # BWA index reference 
# bwa index ${ref}


# # BWA alignment
# bwa mem -t 4 -R "@RG\tID:SRR062634\tPL:ILLUMINA\tSM:SRR062634" ${ref} ${reads}/SRR062634_1.filt.fastq.gz ${reads}/SRR062634_2.filt.fastq.gz > ${aligned_reads}/SRR062634.paired.sam




# # -----------------------------------------
# # STEP 3: Mark Duplicates and Sort - GATK4
# # -----------------------------------------

# echo "STEP 3: Mark Duplicates and Sort - GATK4"

# gatk MarkDuplicatesSpark -I ${aligned_reads}/SRR062634.paired.sam -O ${aligned_reads}/SRR062634_sorted_dedup_reads.bam



# # ----------------------------------
# # STEP 4: Base quality recalibration
# # ----------------------------------


# echo "STEP 4: Base quality recalibration"

# # 1. build the model
# gatk BaseRecalibrator -I ${aligned_reads}/SRR062634_sorted_dedup_reads.bam -R ${ref} --known-sites ${known_sites} -O ${data}/recal_data.table


# # 2. Apply the model to adjust the base quality scores
# gatk ApplyBQSR -I ${aligned_reads}/SRR062634_sorted_dedup_reads.bam -R ${ref} --bqsr-recal-file {$data}/recal_data.table -O ${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam 



# # -----------------------------------------------
# # STEP 5: Collect Alignment & Insert Size Metrics
# # -----------------------------------------------


# echo "STEP 5: Collect Alignment & Insert Size Metrics"

# gatk CollectAlignmentSummaryMetrics R=${ref} I=${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam O=${aligned_reads}/alignment_metrics.txt
# gatk CollectInsertSizeMetrics INPUT=${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam OUTPUT=${aligned_reads}/insert_size_metrics.txt HISTOGRAM_FILE=${aligned_reads}/insert_size_histogram.pdf



# # ----------------------------------------------
# # STEP 6: Call Variants - gatk haplotype caller
# # ----------------------------------------------

# echo "STEP 6: Call Variants - gatk haplotype caller"

# gatk HaplotypeCaller -R ${ref} -I ${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam -O ${results}/raw_variants.vcf



# # extract SNPs & INDELS

# gatk SelectVariants -R ${ref} -V ${results}/raw_variants.vcf --select-type SNP -O ${results}/raw_snps.vcf
# gatk SelectVariants -R ${ref} -V ${results}/raw_variants.vcf --select-type INDEL -O ${results}/raw_indels.vcf

