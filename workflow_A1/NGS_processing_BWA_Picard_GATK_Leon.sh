#!/bin/bash
#SBATCH -p ngs96G
#SBATCH -c 28
#SBATCH --mem=92g
#SBATCH -A MST109178
#SBATCH -J A2974-PCR-NGS
#SBATCH -o Sbatch_BWAMEM_PICARD_GATK_ANNOVAR_Processing_Info.out
#SBATCH -e Sbatch_BWAMEM_PICARD_GATK_ANNOVAR_Processing_Info.err
#SBATCH --mail-user=leon9139@gmail.com
#SBATCH --mail-type=FAIL,END

#==========================Input Info (Taiwania)===================================#

Input_Deafness_PWD=/home/u8454250/NGS_analysis/A2796-PCR_S25_L001_NGS_20220704_test

#### while input file: FASTQ => BAM ####

Input_R1=A2796-PCR_S25_L001_R1_001.fastq.gz
Input_R2=A2796-PCR_S25_L001_R2_001.fastq.gz


### output file ###

Outputname=A2796-PCR_S25

mkdir /$Input_Deafness_PWD/$Outputname


### while input file: SAM/BAM => FASTQ => BAM ###

#Input_SAM=GJB2-50-Q30-3-end-read_S1.sam
#Input_BAM=GJB2-20-Q30-3-end-read_S1.bam

#Input_R1=$Outputname.R1_001.fastq
#Input_R2=$Outputname.R2_001.fastq


#===========================Other Fixed condition=================================#

MappingRef=/home/u8454250/NGS_reference_hg19.NC_012920/ucsc.hg19.NC_012920.fasta
AnnoVar_VCF_RefFILE=/home/u8454250/NGS_reference_hg19.NC_012920/dbsnp_137.hg19.rmMT.vcf
AnnoVar_database=/staging/reserve/paylong_ntu/AI_SHARE/reference/annovar_2016Feb01/humandb

BWAMEM_SOFT=/work/opt/ohpc/Taiwania3/pkg/biology/BWA/BWA_v0.7.17/bwa
PICARD_SOFT=/home/u8454250/NGS_software/picard-tools-1.134/picard.jar
GATK_SOFT=/work/opt/ohpc/Taiwania3/pkg/biology/GATK/gatk_v3.8.1.0/GenomeAnalysisTK.jar
GATK_v4_SOFT=/work/opt/ohpc/Taiwania3/pkg/biology/GATK/gatk_v4.2.0.0/gatk-package-4.2.0.0-local.jar
ANNOVAR_SOFT=/work/opt/ohpc/Taiwania3/pkg/biology/ANNOVAR/annovar_20210819/table_annovar.pl

#===============================Data Processing====================================#
######                                           ######
#############  SAM/BAM to FASTQ  ######################
######                                           ######

# Step0 # Transfering #~Picard~# Transfering initial SAM to FASTQ
#java  -Xmx16g -jar $PICARD_SOFT SamToFastq I=$Input_Deafness_PWD/$Input_SAM F=$Input_Deafness_PWD/$Outputname/$Input_R1 F2=$Input_Deafness_PWD/$Outputname/$Input_R2

# Step0 # Transfering #~Picard~# Transfering initial BAM to FASTQ
#java  -Xmx16g -jar $PICARD_SOFT SamToFastq I=$Input_Deafness_PWD/$Input_BAM F=$Input_Deafness_PWD/$Input_R1 F2=$Input_Deafness_PWD/$Input_R2

######                                           ######
############# FASTQ to mapped-BAM #####################
######                                           ######
	
# Step1 # Mapping #~BWAMEM~# Mapping fastq to sam

$BWAMEM_SOFT mem -t 16 -R '@RG\tID:NGS_sample\tLB:NGS_sample\tSM:NGS_sample\tPL:ILLUMINA\' $MappingRef $Input_Deafness_PWD/$Input_R1 $Input_Deafness_PWD/$Input_R2 > $Input_Deafness_PWD/$Outputname/$Outputname.bwamem.sam &&

	# StepM-2 # MergeSamFiles #~Picard~# merge SAM or BAM
	#java -Xmx16g -jar $PICARD_SOFT MergeSamFiles I=$Input_Deafness_PWD/$TempInput01.bwamem.sam I=$Input_Deafness_PWD/$TempInput02.bwamem.sam I=$Input_Deafness_PWD/$TempInput03.bwamem.sam I=$Input_Deafness_PWD/$TempInput04.bwamem.sam I=$Input_Deafness_PWD/$TempInput05.bwamem.sam O=$Input_Deafness_PWD/$Outputname.bwamem.sam

# Step2 # SortSAM #~Picard~# sam to bam
java -Xmx80g -jar $PICARD_SOFT SortSam INPUT=$Input_Deafness_PWD/$Outputname/$Outputname.bwamem.sam OUTPUT=$Input_Deafness_PWD/$Outputname/$Outputname.bwamem.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true &&

# Step3 # MarkDuplicate #~Picard~# bam to markduplicates
java -Xmx80g -jar $PICARD_SOFT MarkDuplicates INPUT=$Input_Deafness_PWD/$Outputname/$Outputname.bwamem.bam OUTPUT=$Input_Deafness_PWD/$Outputname/$Outputname.bwamem.marked.bam METRICS_FILE=$Input_Deafness_PWD/$Outputname/$Outputname.bwamem_metrics VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true &&

# Step4-1 # Realignment - Realigner Target Creator #~GATK~# Select the range of realignment
java -Xmx80g -jar $GATK_SOFT -T RealignerTargetCreator -R $MappingRef -I $Input_Deafness_PWD/$Outputname/$Outputname.bwamem.marked.bam -o $Input_Deafness_PWD/$Outputname/$Outputname.bwamem.marked.bam.intervals -nt 16 &&

# Step4-2 # Realignment - Indel Realigner #~GATK~# Indel realignment 
java -Xmx80g -jar $GATK_SOFT -T IndelRealigner -R $MappingRef -I $Input_Deafness_PWD/$Outputname/$Outputname.bwamem.marked.bam -targetIntervals $Input_Deafness_PWD/$Outputname/$Outputname.bwamem.marked.bam.intervals -o $Input_Deafness_PWD/$Outputname/$Outputname.bwamem.marked.realigned.bam &&

java  -Xmx80g -jar $PICARD_SOFT FixMateInformation INPUT=$Input_Deafness_PWD/$Outputname/$Outputname.bwamem.marked.realigned.bam OUTPUT=$Input_Deafness_PWD/$Outputname/$Outputname.bwamem.marked.realigned.fixed.bam SO=coordinate VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true &&

# Step5-1 # Base Recalibration - Base Recalibrator #~GATK~# Detect systematic errors in base quality scores
java  -Xmx80g -jar $GATK_SOFT -T BaseRecalibrator -R $MappingRef -knownSites $AnnoVar_VCF_RefFILE -I $Input_Deafness_PWD/$Outputname/$Outputname.bwamem.marked.realigned.fixed.bam -o $Input_Deafness_PWD/$Outputname/$Outputname.bwamem.marked.realigned.fixed.bam.recal_data.grp -rf BadCigar -nct 16 &&

# Step5-2 # Base Recalibration - Base quality Realigning #~GATK~# Realign recalibrate bam file
java  -Xmx80g -jar $GATK_SOFT -T PrintReads -R $MappingRef -I $Input_Deafness_PWD/$Outputname/$Outputname.bwamem.marked.realigned.fixed.bam -BQSR $Input_Deafness_PWD/$Outputname/$Outputname.bwamem.marked.realigned.fixed.bam.recal_data.grp -o $Input_Deafness_PWD/$Outputname/$Outputname.bwamem.marked.realigned.fixed.recal.bam -nct 16 &&

# Step6 # SortSam #~Picard~#
java  -Xmx80g -jar $PICARD_SOFT SortSam INPUT=$Input_Deafness_PWD/$Outputname/$Outputname.bwamem.marked.realigned.fixed.recal.bam OUTPUT=$Input_Deafness_PWD/$Outputname/$Outputname.bwamem.marked.realigned.fixed.recal.indexed.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true


######                                                          ######
######## mapped-BAM to vcf/txt (Normal Genotype) #####################
######                                                          ######

# Step7-1 # HaplotypeCaller #~GATK~# Variant calling by HaplotypeCaller

java  -Xmx80g -jar $GATK_SOFT -T HaplotypeCaller  -l INFO -R $MappingRef -I $Input_Deafness_PWD/$Outputname/$Outputname.bwamem.marked.realigned.fixed.recal.indexed.bam --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000  --dbsnp $AnnoVar_VCF_RefFILE --max_alternate_alleles 30 -o $Input_Deafness_PWD/$Outputname/$Outputname.bwamem.haplotype.SnpIndel.g.vcf.gz -nct 16 

# Step8 # GenotypeCaller #~GATK~# GenotpeGVCF

java  -Xmx80g -jar $GATK_SOFT -T GenotypeGVCFs -R $MappingRef -V $Input_Deafness_PWD/$Outputname/$Outputname.bwamem.haplotype.SnpIndel.g.vcf.gz -o $Input_Deafness_PWD/$Outputname/$Outputname.bwamem.haplotype.SnpIndel.vcf.gz 

# Step9-1 # Filtering #~GATK~# Filtering of Variants for Haplotype-Genotype

java  -Xmx80g -jar $GATK_SOFT -T VariantFiltration -R $MappingRef --variant $Input_Deafness_PWD/$Outputname/$Outputname.bwamem.haplotype.SnpIndel.vcf.gz  -o $Input_Deafness_PWD/$Outputname/$Outputname.bwamem.haplotype.SnpIndel.filtered.vcf.gz --clusterWindowSize 10 --filterExpression "DP < 5" --filterName "LowCoverage" --filterExpression "QUAL < 30.0" --filterName "VeryLowQual" --filterExpression "QUAL > 30.0 && QUAL < 50.0" --filterName "LowQual" --filterExpression "QD < 1.5" --filterName "LowQD" 

# Step10-1 # Annotation #~Perl~#annotation by Annovar for Haplotype-Genotype

perl $ANNOVAR_SOFT $Input_Deafness_PWD/$Outputname/$Outputname.bwamem.haplotype.SnpIndel.filtered.vcf.gz $AnnoVar_database -buildver hg19 -out $Input_Deafness_PWD/$Outputname/$Outputname.bwamem_haplotype -remove -protocol refGene,cytoBand,knownGene,ensGene,gnomad211_genome,avsnp150,TaiwanBiobank-official,gnomad211_exome,TWB_1497_joing_calling_AF,intervar_20180118,clinvar_20210501,cosmic_coding_GRCh37_v92,cosmic_noncoding_GRCh37_v92,icgc28,dbnsfp41a,cg69,kaviar_20150923,dbscsnv11,spidex,gwava,wgRna,targetScanS -operation gx,r,gx,gx,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f -arg '-splicing 10',,,,,,,,,,,,,,,,,,,,, -nastring . -vcfinput -polish --maxgenethread 20 --thread 20 

######                                                          ######
######## mapped-BAM to vcf/txt (Mutect2 Genotype) ####################
######                                                          ######

# Step7-2 # HaplotypeCaller #~GATK~# Variant calling by Mutect2

java  -Xmx80g -jar $GATK_v4_SOFT Mutect2 -R $MappingRef -I $Input_Deafness_PWD/$Outputname/$Outputname.bwamem.marked.realigned.fixed.recal.indexed.bam -O $Input_Deafness_PWD/$Outputname/$Outputname.bwamem.Mutect2.vcf.gz 


# Step9-2 # Filtering #~GATK~# Filtering of Variants for Mutect2

java  -Xmx80g -jar $GATK_v4_SOFT FilterMutectCalls -R $MappingRef --variant $Input_Deafness_PWD/$Outputname/$Outputname.bwamem.Mutect2.vcf.gz -O $Input_Deafness_PWD/$Outputname/$Outputname.bwamem.filtered.Mutect2.vcf.gz 


# Step10-2 # Annotation #~Perl~#annotation by Annovar for Mutect2

perl $ANNOVAR_SOFT $Input_Deafness_PWD/$Outputname/$Outputname.bwamem.filtered.Mutect2.vcf.gz $AnnoVar_database -buildver hg19 -out $Input_Deafness_PWD/$Outputname/$Outputname.bwamem_Mutect2 -remove -protocol refGene,cytoBand,knownGene,ensGene,gnomad211_genome,avsnp150,TaiwanBiobank-official,gnomad211_exome,TWB_1497_joing_calling_AF,intervar_20180118,clinvar_20210501,cosmic_coding_GRCh37_v92,cosmic_noncoding_GRCh37_v92,icgc28,dbnsfp41a,cg69,kaviar_20150923,dbscsnv11,spidex,gwava,wgRna,targetScanS -operation gx,r,gx,gx,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f -arg '-splicing 10',,,,,,,,,,,,,,,,,,,,, -nastring . -vcfinput -polish --maxgenethread 20 --thread 20 

echo "The job $Outputname  has been done at `date`"

mv $Input_Deafness_PWD/bwamem_all_out.txt $Input_Deafness_PWD/$Outputname
mv $Input_Deafness_PWD/bwamem_all_err.txt $Input_Deafness_PWD/$Outputname
mv $Input_Deafness_PWD/$Input_SAM $Input_Deafness_PWD/$Outputname
mv $Input_Deafness_PWD/$Input_R1 $Input_Deafness_PWD/$Outputname
mv $Input_Deafness_PWD/$Input_R2 $Input_Deafness_PWD/$Outputname