#!/bin/bash
#SBATCH -p ngs96G
#SBATCH -c 28
#SBATCH --mem=92g
#SBATCH -A MST109178
#SBATCH -J MG244
#SBATCH -o /staging/reserve/paylong_ntu/Research/DF/DF_20220923_Research/MG244/MG244_20220923_bwamem_all_out.txt
#SBATCH -e /staging/reserve/paylong_ntu/Research/DF/DF_20220923_Research/MG244/MG244_20220923_bwamem_all_err.txt
#SBATCH --mail-user=sharonchiu0104@gmail.com
#SBATCH --mail-type=FAIL,END

# This code is used for NGS germline variant calling (SNV + Indels)
# 用 “&&” 分開兩個指令, 即第一道指令執行成功後, 才會執行第二道指令

echo "processes directories"
# before start, you should copy the fastq to fastq_folder
# fold to be set
sample_id="MG244" #change as you needed
user_dir="/work/u3003390/"
fastq_dir="/work/u3003390/FASTQ" #data fold
temp_dir="/work/u3003390/TEMP"
release_dir="/work/u3003390/RESULT"

cd $user_dir
set -euo pipefail

# Update with the fullpath location of your sample fastq
fastq_1="${fastq_dir}/${sample_id}*R1*.gz" #NGS1_20170305A.R1.fastq.gz
fastq_2="${fastq_dir}/${sample_id}*R2*.gz"  #If using Illumina paired data
sample="SM_"${SampleName}
group="GP_"${SampleName}
platform="ILLUMINA"


# Update with the location of the reference data files (hg 19)
# May use hg 38 if WGS/WES 
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

# Other settings
nt=40 #number of threads to use in computation

# Set working directory
Time=`date +%Y/%m/%d %H:%M`
Date=`date +%Y/%m/%d`
logfile=$temp_dir/${Time}_run.log
set -x
exec 3<&1 4<&2 #???
exec >$logfile 2>&1


#bwa mem (fastq -> sam)
module load biology/BWA/0.7.17
bwa mem \
-t 16 -R '@RG\tID:MG244_20220923_bwamem\tLB:MG244_20220923_bwamem\tSM:MG244_20220923_bwamem\tPL:ILLUMINA\' \
$fasta \
$fastq_1 \
$fastq_2 \
> $temp_dir/${sample_id}.$Date.bwamem.sam &&

#Sortsam (sam -> bwamem.bam)#
module load biology/Picard/2.27.4
picard SortSam \
INPUT=$temp_dir/${sample_id}.$Date.bwamem.sam \
OUTPUT=$temp_dir/${sample_id}.$Date.bwamem.bam \
SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true &&

#MarkDuplication (bwamem.bam->bwamem.marked.bam)#
picard MarkDuplicates \
INPUT=$temp_dir/${sample_id}.$Date.bwamem.bam \
OUTPUT=$temp_dir/${sample_id}.$Date.bwamem.marked.bam \
METRICS_FILE=$temp_dir/${sample_id}.$Date.bwamem_metrics #file to write duplication metrices\
VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true &&


# #RealignerTargetCreator (bwamem.marked.bam -> bwamem.marked.bam.intervals)
# java  -Xmx80g -jar /home/u1151339/software/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar \
# -T RealignerTargetCreator \
# -R /home/u1151339/reference/ucsc.hg19.NC_012920.fasta \
# -I /staging/reserve/paylong_ntu/Research/DF/DF_20220923_Research/MG244/MG244_20220923_bwamem.marked.bam \
# -o /staging/reserve/paylong_ntu/Research/DF/DF_20220923_Research/MG244/MG244_20220923_bwamem.marked.bam.intervals \
# -nt 16 &&

# #IndelRealigner (bwamem.marked.realigned.bam)
# java  -Xmx80g -jar /home/u1151339/software/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar \
# -T IndelRealigner -R /home/u1151339/reference/ucsc.hg19.NC_012920.fasta \
# -I /staging/reserve/paylong_ntu/Research/DF/DF_20220923_Research/MG244/MG244_20220923_bwamem.marked.bam \
# -targetIntervals /staging/reserve/paylong_ntu/Research/DF/DF_20220923_Research/MG244/MG244_20220923_bwamem.marked.bam.intervals \
# -o /staging/reserve/paylong_ntu/Research/DF/DF_20220923_Research/MG244/MG244_20220923_bwamem.marked.realigned.bam && 

#FixMateInformation (bwamem.marked.realigned.fixed.bam)
picard FixMateInformation \
INPUT=$temp_dir/${sample_id}.bwamem.marked.realigned.bam \
OUTPUT=$temp_dir/${sample_id}.bwamem.marked.realigned.fixed.bam \
SO=coordinate VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true && 

#BaseRecalibrator (bwamem.marked.realigned.fixed.bam.recal_data.grp)
ml load biology/GATK/4.2.3.0
gatk
-T BaseRecalibrator \
-R $fasta \
-knownSites $dbsnp \
-I $temp_dir/${sample_id}.bwamem.marked.realigned.fixed.bam \
-o $temp_dir/${sample_id}.bwamem.marked.realigned.fixed.bam.recal_data.table \
-rf BadCigar -nct 16 && 

# ApplyBQSR (to be edited)
gatk ApplyBQSR \
-I $temp_dir/${sample_id}.bwamem.marked.realigned.fixed.bam \
-R $fasta \
--bqsr-recal-file $temp_dir/${sample_id}.bwamem.marked.realigned.fixed.bam.recal_data.table \
-O $temp_dir/${sample_id}.bwamem.marked.realigned.fixed.recal.bam

#SortSam (bwamem.marked.realigned.fixed.recal.indexed.bam)
picard SortSam \
INPUT= $temp_dir/${sample_id}.bwamem.marked.realigned.fixed.recal.bam \
OUTPUT= $temp_dir/${sample_id}.bwamem.marked.realigned.fixed.recal.indexed.bam \
SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true &&

#change name (indexed.bai->indexed.bam.bai)?
cp $temp_dir/${sample_id}.bwamem.marked.realigned.fixed.recal.indexed.bai \
   $temp_dir/${sample_id}.bwamem.marked.realigned.fixed.recal.indexed.bam.bai

# #findSHv2_20170921
# perl /home/u1151339/software/BP/findSHv2_20170921.pl \
# /staging/reserve/paylong_ntu/Research/DF/DF_20220923_Research/MG244/MG244_20220923_bwamem.sam \
# /staging/reserve/paylong_ntu/Research/DF/DF_20220923_Research/MG244/MG244_20220923_bwamem_FindSH.txt

# #findBPv2_20170921
# perl /home/u1151339/software/BP/findBPv2_20170921.pl \
# /staging/reserve/paylong_ntu/Research/DF/DF_20220923_Research/MG244/MG244_20220923_bwamem_FindSH.txt \
# /staging/reserve/paylong_ntu/Research/DF/DF_20220923_Research/MG244/MG244_20220923_bwamem_findBP_original.txt \
# /staging/reserve/paylong_ntu/Research/DF/DF_20220923_Research/MG244/MG244_20220923_bwamem_findBP_merged.txt

# #bam-readcount (.bam -> bwamem_nucleotide_composition.txt)
# /home/u1151339/software/bam-readcount-master/bam-readcount \
# -f /home/u1151339/reference/ucsc.hg19.NC_012920.fasta \
# -l /staging/reserve/paylong_ntu/AI_SHARE/Bed/DF/DF_v3/DF_v3.bed \
# -w 0 /staging/reserve/paylong_ntu/Research/DF/DF_20220923_Research/MG244/MG244_20220923_bwamem.marked.realigned.fixed.recal.indexed.bam \
# > /staging/reserve/paylong_ntu/Research/DF/DF_20220923_Research/MG244/MG244_20220923_bwamem_nucleotide_composition.txt

# # copy bamReadCount_rearrangement to case file
# cp /staging/reserve/paylong_ntu/AI_SHARE/Bed/bamReadCount_rearrangement.py \
# /staging/reserve/paylong_ntu/Research/DF/DF_20220923_Research/MG244/
# # copy gene list (DF_v3_genelist.txt) to case file
# cp /staging/reserve/paylong_ntu/AI_SHARE/Bed/DF/DF_v3/DF_v3_genelist.txt \
# /staging/reserve/paylong_ntu/Research/DF/DF_20220923_Research/MG244/

# #copy total_average.py to case file
# cp /staging/reserve/paylong_ntu/AI_SHARE/Bed/total_average.py \
# /staging/reserve/paylong_ntu/Research/DF/DF_20220923_Research/MG244/
# #exec total_average.py
# /home/u1151339/software/Python/Python-3.6.1/python \
# /staging/reserve/paylong_ntu/Research/DF/DF_20220923_Research/MG244/total_average.py

# #bamReadCount_rearrangement.py
# /home/u1151339/software/Python/Python-3.6.1/python \
# /staging/reserve/paylong_ntu/Research/DF/DF_20220923_Research/MG244/bamReadCount_rearrangement.py

# #copy read_average_v2.py to case file
# cp /staging/reserve/paylong_ntu/AI_SHARE/Bed/read_average_v2.py \
# /staging/reserve/paylong_ntu/Research/DF/DF_20220923_Research/MG244/
# #exec read_average_v2.py
# /home/u1151339/software/Python/Python-3.6.1/python \
# /staging/reserve/paylong_ntu/Research/DF/DF_20220923_Research/MG244/read_average_v2.py

#HaplotypeCaller (bwamem.marked.realigned.fixed.recal.indexed.bam -> bwamem.haplotype.SnpIndel.g.vcf.gz)
gatk \
-T HaplotypeCaller  \
-l INFO \
#?
-R /home/u1151339/reference/ucsc.hg19.NC_012920.fasta \
-I $temp_dir/${sample_id}.bwamem.marked.realigned.fixed.recal.indexed.bam \
--emitRefConfidence GVCF \
# --variant_index_type LINEAR \
# --variant_index_parameter 128000  \
--dbsnp $dbsnp \
--max_alternate_alleles 30 \ 
#why 30?
-o $temp_dir/${sample_id}.bwamem.haplotype.SnpIndel.g.vcf.gz \
-nct 16 

#GenotypeGVCFs (g.vcf -> vcf)
gatk \
-T GenotypeGVCFs \
-R $fasta \
-V $temp_dir/${sample_id}.bwamem.haplotype.SnpIndel.g.vcf.gz \
-o $temp_dir/${sample_id}.bwamem.haplotype.SnpIndel.vcf.gz 

#VariantFiltration (vcf -> filtered.vcf)
gatk \
-T VariantFiltration \
-R $fasta \
--variant $temp_dir/${sample_id}.bwamem.haplotype.SnpIndel.vcf.gz\
-o $temp_dir/${sample_id}.bwamem.haplotype.SnpIndel.filtered.vcf.gz \
--clusterWindowSize 10 \
--filterExpression "DP < 5" --filterName "LowCoverage" --filterExpression "QUAL < 30.0" --filterName "VeryLowQual" --filterExpression "QUAL > 30.0 && QUAL < 50.0" --filterName "LowQual" --filterExpression "QD < 1.5" --filterName "LowQD" 

#exec table_annovar.pl 
module load biology/ANNOVAR/2020-06-08
table_annovar.pl \
$temp_dir/${sample_id}.bwamem.filtered.haplotype.SnpIndel.vcf.gz \
/staging/reserve/paylong_ntu/AI_SHARE/reference/annovar_2016Feb01/humandb/ \
-buildver hg19 \
-out $release_dir/${sample_id}.annotate \
-remove \
-protocol refGene,cytoBand,knownGene,ensGene,gnomad211_genome,avsnp150,TaiwanBiobank-official,gnomad211_exome,TWB_1497_joing_calling_AF,intervar_20180118,clinvar_20210501,cosmic_coding_GRCh37_v92,cosmic_noncoding_GRCh37_v92,icgc28,dbnsfp41a,cg69,kaviar_20150923,dbscsnv11,spidex,gwava,wgRna,targetScanS \
-operation gx,r,gx,gx,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f \
-arg '-splicing 10',,,,,,,,,,,,,,,,,,,,, -nastring . -vcfinput -polish --maxgenethread 20 --thread 20 

# Mutect2 (somatic variant calling)
# java  -Xmx80g -jar /home/u1151339/software/GATK_v4.1.8.0/gatk-package-4.1.8.0-local.jar Mutect2 -R /home/u1151339/reference/ucsc.hg19.NC_012920.fasta -I /staging/reserve/paylong_ntu/Research/DF/DF_20220923_Research/MG244/MG244_20220923_bwamem.marked.realigned.fixed.recal.indexed.bam -O /staging/reserve/paylong_ntu/Research/DF/DF_20220923_Research/MG244/MG244_20220923_bwamem.Mutect2.vcf.gz 
# java  -Xmx80g -jar /home/u1151339/software/GATK_v4.1.8.0/gatk-package-4.1.8.0-local.jar FilterMutectCalls -R /home/u1151339/reference/ucsc.hg19.NC_012920.fasta --variant /staging/reserve/paylong_ntu/Research/DF/DF_20220923_Research/MG244/MG244_20220923_bwamem.Mutect2.vcf.gz -O /staging/reserve/paylong_ntu/Research/DF/DF_20220923_Research/MG244/MG244_20220923_bwamem.filtered.Mutect2.vcf.gz 
# perl /work/opt/ohpc/Taiwania3/pkg/biology/ANNOVAR/annovar_20210819/table_annovar.pl /staging/reserve/paylong_ntu/Research/DF/DF_20220923_Research/MG244/MG244_20220923_bwamem.filtered.Mutect2.vcf.gz /staging/reserve/paylong_ntu/AI_SHARE/reference/annovar_2016Feb01/humandb/ -buildver hg19 -out /staging/reserve/paylong_ntu/Research/DF/DF_20220923_Research/MG244/MG244_20220923_bwamem_Mutect2 -remove -protocol refGene,cytoBand,knownGene,ensGene,gnomad211_genome,avsnp150,TaiwanBiobank-official,gnomad211_exome,TWB_1497_joing_calling_AF,intervar_20180118,clinvar_20210501,cosmic_coding_GRCh37_v92,cosmic_noncoding_GRCh37_v92,icgc28,dbnsfp41a,cg69,kaviar_20150923,dbscsnv11,spidex,gwava,wgRna,targetScanS -operation gx,r,gx,gx,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f -arg '-splicing 10',,,,,,,,,,,,,,,,,,,,, -nastring . -vcfinput -polish --maxgenethread 20 --thread 20 
