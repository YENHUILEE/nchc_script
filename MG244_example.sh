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

# 用 “&&” 分開兩個指令, 即第一道指令執行成功後, 才會執行第二道指令

#bwa mem (fastq -> sam)
/home/u1151339/software/bwa-0.7.12/bwa mem \
-t 16 -R '@RG\tID:MG244_20220923_bwamem\tLB:MG244_20220923_bwamem\tSM:MG244_20220923_bwamem\tPL:ILLUMINA\' \
/home/u1151339/reference/ucsc.hg19.NC_012920.fasta \
/staging/reserve/paylong_ntu/Research/DF/DF_20220923_Research/MG244_*R1*.fastq.gz \
/staging/reserve/paylong_ntu/Research/DF/DF_20220923_Research/MG244_*R2*.fastq.gz \
> /staging/reserve/paylong_ntu/Research/DF/DF_20220923_Research/MG244/MG244_20220923_bwamem.sam &&

#Sortsam (sam -> bwamem.bam)#
java  -Xmx80g -jar /home/u1151339/software/picard-tools-1.134/picard.jar SortSam \
INPUT=/staging/reserve/paylong_ntu/Research/DF/DF_20220923_Research/MG244/MG244_20220923_bwamem.sam \
OUTPUT=/staging/reserve/paylong_ntu/Research/DF/DF_20220923_Research/MG244/MG244_20220923_bwamem.bam \
SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true &&

#MarkDuplication (bwamem.bam->bwamem.marked.bam)#
java  -Xmx80g -jar /home/u1151339/software/picard-tools-1.134/picard.jar MarkDuplicates \
INPUT=/staging/reserve/paylong_ntu/Research/DF/DF_20220923_Research/MG244/MG244_20220923_bwamem.bam \
OUTPUT=/staging/reserve/paylong_ntu/Research/DF/DF_20220923_Research/MG244/MG244_20220923_bwamem.marked.bam \
METRICS_FILE=/staging/reserve/paylong_ntu/Research/DF/DF_20220923_Research/MG244/MG244_20220923_bwamem_metrics #file to write duplication metrices\
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
java  -Xmx80g -jar /home/u1151339/software/picard-tools-1.134/picard.jar FixMateInformation \
INPUT=/staging/reserve/paylong_ntu/Research/DF/DF_20220923_Research/MG244/MG244_20220923_bwamem.marked.realigned.bam \
OUTPUT=/staging/reserve/paylong_ntu/Research/DF/DF_20220923_Research/MG244/MG244_20220923_bwamem.marked.realigned.fixed.bam \
SO=coordinate VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true && 

#BaseRecalibrator (bwamem.marked.realigned.fixed.bam.recal_data.grp)
java  -Xmx80g -jar /home/u1151339/software/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar \
-T BaseRecalibrator -R /home/u1151339/reference/ucsc.hg19.NC_012920.fasta \
-knownSites /home/u1151339/reference/dbsnp_137.hg19.rmMT.vcf \
-I /staging/reserve/paylong_ntu/Research/DF/DF_20220923_Research/MG244/MG244_20220923_bwamem.marked.realigned.fixed.bam \
-o /staging/reserve/paylong_ntu/Research/DF/DF_20220923_Research/MG244/MG244_20220923_bwamem.marked.realigned.fixed.bam.recal_data.grp \
-rf BadCigar -nct 16 && 

# ApplyBQSR (to be edited)
gatk ApplyBQSR \
-I ${workdir}/${SampleName}_sorted_dedup_reads.bam \
-R $fasta --bqsr-recal-file ${workdir}/recal_data.table \
-O ${workdir}/${SampleName}_sorted_dedup_bqsr_reads.bam 

#SortSam (bwamem.marked.realigned.fixed.recal.indexed.bam)
java  -Xmx80g -jar /home/u1151339/software/picard-tools-1.134/picard.jar SortSam \
INPUT=/staging/reserve/paylong_ntu/Research/DF/DF_20220923_Research/MG244/MG244_20220923_bwamem.marked.realigned.fixed.recal.bam \
OUTPUT=/staging/reserve/paylong_ntu/Research/DF/DF_20220923_Research/MG244/MG244_20220923_bwamem.marked.realigned.fixed.recal.indexed.bam \
SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true &&

#change name (indexed.bai->indexed.bam.bai)?
cp /staging/reserve/paylong_ntu/Research/DF/DF_20220923_Research/MG244/MG244_20220923_bwamem.marked.realigned.fixed.recal.indexed.bai \
   /staging/reserve/paylong_ntu/Research/DF/DF_20220923_Research/MG244/MG244_20220923_bwamem.marked.realigned.fixed.recal.indexed.bam.bai

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
java  -Xmx80g -jar /home/u1151339/software/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar \
-T HaplotypeCaller  \
-l INFO \
#?
-R /home/u1151339/reference/ucsc.hg19.NC_012920.fasta \
-I /staging/reserve/paylong_ntu/Research/DF/DF_20220923_Research/MG244/MG244_20220923_bwamem.marked.realigned.fixed.recal.indexed.bam \
--emitRefConfidence GVCF \
# --variant_index_type LINEAR \
# --variant_index_parameter 128000  \
--dbsnp /home/u1151339/reference/dbsnp_137.hg19.rmMT.vcf \
--max_alternate_alleles 30 \ 
#why 30?
-o /staging/reserve/paylong_ntu/Research/DF/DF_20220923_Research/MG244/MG244_20220923_bwamem.haplotype.SnpIndel.g.vcf.gz \
-nct 16 

#GenotypeGVCFs (g.vcf -> vcf)
java  -Xmx80g -jar /home/u1151339/software/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar \
-T GenotypeGVCFs \
-R /home/u1151339/reference/ucsc.hg19.NC_012920.fasta \
-V /staging/reserve/paylong_ntu/Research/DF/DF_20220923_Research/MG244/MG244_20220923_bwamem.haplotype.SnpIndel.g.vcf.gz \
-o /staging/reserve/paylong_ntu/Research/DF/DF_20220923_Research/MG244/MG244_20220923_bwamem.haplotype.SnpIndel.vcf.gz 

#VariantFiltration (vcf -> filtered.vcf)
java  -Xmx80g -jar /home/u1151339/software/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar \
-T VariantFiltration \
-R /home/u1151339/reference/ucsc.hg19.NC_012920.fasta \
--variant /staging/reserve/paylong_ntu/Research/DF/DF_20220923_Research/MG244/MG244_20220923_bwamem.haplotype.SnpIndel.vcf.gz\
-o /staging/reserve/paylong_ntu/Research/DF/DF_20220923_Research/MG244/MG244_20220923_bwamem.filtered.haplotype.SnpIndel.vcf.gz \
--clusterWindowSize 10 \
--filterExpression "DP < 5" --filterName "LowCoverage" --filterExpression "QUAL < 30.0" --filterName "VeryLowQual" --filterExpression "QUAL > 30.0 && QUAL < 50.0" --filterName "LowQual" --filterExpression "QD < 1.5" --filterName "LowQD" 

#exec table_annovar.pl 
perl /work/opt/ohpc/Taiwania3/pkg/biology/ANNOVAR/annovar_20210819/table_annovar.pl \
/staging/reserve/paylong_ntu/Research/DF/DF_20220923_Research/MG244/MG244_20220923_bwamem.filtered.haplotype.SnpIndel.vcf.gz \
/staging/reserve/paylong_ntu/AI_SHARE/reference/annovar_2016Feb01/humandb/ \
-buildver hg19 \
-out /staging/reserve/paylong_ntu/Research/DF/DF_20220923_Research/MG244/MG244_20220923_bwamem_haplotype \
-remove \
-protocol refGene,cytoBand,knownGene,ensGene,gnomad211_genome,avsnp150,TaiwanBiobank-official,gnomad211_exome,TWB_1497_joing_calling_AF,intervar_20180118,clinvar_20210501,cosmic_coding_GRCh37_v92,cosmic_noncoding_GRCh37_v92,icgc28,dbnsfp41a,cg69,kaviar_20150923,dbscsnv11,spidex,gwava,wgRna,targetScanS \
-operation gx,r,gx,gx,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f \
-arg '-splicing 10',,,,,,,,,,,,,,,,,,,,, -nastring . -vcfinput -polish --maxgenethread 20 --thread 20 

# Mutect2 (somatic variant calling)
# java  -Xmx80g -jar /home/u1151339/software/GATK_v4.1.8.0/gatk-package-4.1.8.0-local.jar Mutect2 -R /home/u1151339/reference/ucsc.hg19.NC_012920.fasta -I /staging/reserve/paylong_ntu/Research/DF/DF_20220923_Research/MG244/MG244_20220923_bwamem.marked.realigned.fixed.recal.indexed.bam -O /staging/reserve/paylong_ntu/Research/DF/DF_20220923_Research/MG244/MG244_20220923_bwamem.Mutect2.vcf.gz 
# java  -Xmx80g -jar /home/u1151339/software/GATK_v4.1.8.0/gatk-package-4.1.8.0-local.jar FilterMutectCalls -R /home/u1151339/reference/ucsc.hg19.NC_012920.fasta --variant /staging/reserve/paylong_ntu/Research/DF/DF_20220923_Research/MG244/MG244_20220923_bwamem.Mutect2.vcf.gz -O /staging/reserve/paylong_ntu/Research/DF/DF_20220923_Research/MG244/MG244_20220923_bwamem.filtered.Mutect2.vcf.gz 
# perl /work/opt/ohpc/Taiwania3/pkg/biology/ANNOVAR/annovar_20210819/table_annovar.pl /staging/reserve/paylong_ntu/Research/DF/DF_20220923_Research/MG244/MG244_20220923_bwamem.filtered.Mutect2.vcf.gz /staging/reserve/paylong_ntu/AI_SHARE/reference/annovar_2016Feb01/humandb/ -buildver hg19 -out /staging/reserve/paylong_ntu/Research/DF/DF_20220923_Research/MG244/MG244_20220923_bwamem_Mutect2 -remove -protocol refGene,cytoBand,knownGene,ensGene,gnomad211_genome,avsnp150,TaiwanBiobank-official,gnomad211_exome,TWB_1497_joing_calling_AF,intervar_20180118,clinvar_20210501,cosmic_coding_GRCh37_v92,cosmic_noncoding_GRCh37_v92,icgc28,dbnsfp41a,cg69,kaviar_20150923,dbscsnv11,spidex,gwava,wgRna,targetScanS -operation gx,r,gx,gx,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f -arg '-splicing 10',,,,,,,,,,,,,,,,,,,,, -nastring . -vcfinput -polish --maxgenethread 20 --thread 20 

