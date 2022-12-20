#!/usr/bin/sh
#SBATCH -A MST109178        # Account name/project number
#SBATCH -J SAMPLE_NAME         # Job name
#SBATCH -p ngs48G           # Partition Name 
#SBATCH -c 14               # Core numbers
#SBATCH --mem=46g           # Memory size
#SBATCH -o out.log          # Path to the standard output file 
#SBATCH -e err.log          # Path to the standard error ouput file
#SBATCH --mail-user=@gmail.com    # email
#SBATCH --mail-type=FAIL              # When to send an email = NONE, BEGIN, END, FAIL, REQUEUE, or ALL


### Read1 and Read2
ID=SAMPLE_ID
READ1=/project/GP1/u3710062/AI_SHARE/rawdata/Deafness/${ID}/panel/${ID}*R1_001.fastq.gz
READ2=/project/GP1/u3710062/AI_SHARE/rawdata/Deafness/${ID}/panel/${ID}*R2_001.fastq.gz
### Parameters
FOLDER=/project/GP1/u3710062/AI_SHARE/shared_scripts/A1_Panel
GATK=/project/GP1/u3710062/AI_SHARE/software/GATK/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar
PICARD=/pkg/biology/Picard/Picard_v2.18.11/picard.jar
DAY=`date +%Y%m%d`
cd ${FOLDER}
mkdir -p ${ID}
cd ${FOLDER}/${ID}
TIME=`date +%Y%m%d%H%M`
logfile=./${TIME}_${ID}_run.log
exec 3<&1 4<&2
exec >$logfile 2>&1
set -euo pipefail
set -x
################################################################ For hg38
Ver=hg38
HG38=/project/GP1/u3710062/AI_SHARE/reference/GATK_bundle/2.8/hg38/Homo_sapiens_assembly38.fasta
BED=/project/GP1/u3710062/AI_SHARE/shared_scripts/A1_Panel/20210304.hg38.interval_list
dbSNP=/project/GP1/u3710062/AI_SHARE/reference/GATK_bundle/2.8/hg38/dbsnp_146.hg38.vcf.gz
mkdir -p ${FOLDER}/${ID}/${Ver}
cd ${FOLDER}/${ID}/${Ver}
### BWA alignment
/pkg/biology/BWA/BWA_v0.7.17/bwa mem -t 16 -R '@RG\tID:'${ID}'_'${DAY}'_'${Ver}'_bwamem\tLB:'${ID}'_'${DAY}'_'${Ver}'_bwamem\tSM:'${ID}'_'${DAY}'_'${Ver}'_bwamem\tPL:ILLUMINA\' ${HG38} ${READ1} ${READ2} > ${ID}_${DAY}_${Ver}_bwamem.sam 

java  -Xmx40g -jar ${PICARD} SortSam INPUT=${ID}_${DAY}_${Ver}_bwamem.sam OUTPUT=${ID}_${DAY}_${Ver}_bwamem.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true 

### To output unmapped.bam
/pkg/biology/SAMtools/SAMtools_v1.10/bin/samtools view -b -f 4 ${ID}_${DAY}_${Ver}_bwamem.bam > ${ID}_${DAY}_${Ver}_bwamem.unmapped.bam

### GATK & Picard pre-process
java  -Xmx40g -jar ${PICARD} MarkDuplicates INPUT=${ID}_${DAY}_${Ver}_bwamem.bam OUTPUT=${ID}_${DAY}_${Ver}_bwamem.marked.bam METRICS_FILE=${ID}_${DAY}_${Ver}_bwamem_metrics VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true 

java  -Xmx40g -jar ${GATK} -T RealignerTargetCreator -R ${HG38} -I ${ID}_${DAY}_${Ver}_bwamem.marked.bam -o ${ID}_${DAY}_${Ver}_bwamem.marked.bam.intervals -nt 16 

java  -Xmx40g -jar ${GATK} -T IndelRealigner -R ${HG38} -I ${ID}_${DAY}_${Ver}_bwamem.marked.bam -targetIntervals ${ID}_${DAY}_${Ver}_bwamem.marked.bam.intervals -o ${ID}_${DAY}_${Ver}_bwamem.marked.realigned.bam 

java  -Xmx40g -jar ${PICARD} FixMateInformation INPUT=${ID}_${DAY}_${Ver}_bwamem.marked.realigned.bam OUTPUT=${ID}_${DAY}_${Ver}_bwamem.marked.realigned.fixed.bam SO=coordinate VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true 

java  -Xmx40g -jar ${GATK} -T BaseRecalibrator -R ${HG38} -knownSites ${dbSNP} -I ${ID}_${DAY}_${Ver}_bwamem.marked.realigned.fixed.bam -o ${ID}_${DAY}_${Ver}_bwamem.marked.realigned.fixed.bam.recal_data.grp -rf BadCigar -nct 16 

java  -Xmx40g -jar ${GATK} -T PrintReads -R ${HG38} -I ${ID}_${DAY}_${Ver}_bwamem.marked.realigned.fixed.bam -BQSR ${ID}_${DAY}_${Ver}_bwamem.marked.realigned.fixed.bam.recal_data.grp -o ${ID}_${DAY}_${Ver}_bwamem.marked.realigned.fixed.recal.bam -nct 16 

java  -Xmx40g -jar ${PICARD} SortSam INPUT=${ID}_${DAY}_${Ver}_bwamem.marked.realigned.fixed.recal.bam OUTPUT=${ID}_${DAY}_${Ver}_bwamem.marked.realigned.fixed.recal.indexed.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true 

cp ${ID}_${DAY}_${Ver}_bwamem.marked.realigned.fixed.recal.indexed.bai ${ID}_${DAY}_${Ver}_bwamem.marked.realigned.fixed.recal.indexed.bam.bai

### Custimized scripts

#perl /home/u4/u3923203/software/BP/findSHv2_20170921.pl ${ID}_${DAY}_${Ver}_bwamem.sam ${ID}_${DAY}_${Ver}_bwamem_FindSH.txt

#perl /home/u4/u3923203/software/BP/findBPv2_20170921.pl ${ID}_${DAY}_${Ver}_bwamem_FindSH.txt ${ID}_${DAY}_${Ver}_bwamem_findBP_original.txt ${ID}_${DAY}_${Ver}_bwamem_findBP_merged.txt

#/home/u4/u3923203/software/bam-readcount-master/bam-readcount -f ${HG38} -l ${BED} -w 0 ${ID}_${DAY}_${Ver}_bwamem.marked.realigned.fixed.recal.indexed.bam > ${ID}_${DAY}_${Ver}_bwamem_nucleotide_composition.txt

#cp /home/u4/u3923203/Bed/bamReadCount_rearrangement.py 

#/home/u4/u3923203/software/Python/Python-3.6.1/python bamReadCount_rearrangement.py

#cp /home/u4/u3923203/Bed/HRD/HRD_v3/HRD_v3_genelist.txt 

#cp /home/u4/u3923203/Bed/total_average.py 

#cp /home/u4/u3923203/Bed/read_average_v2.py 

#/home/u4/u3923203/software/Python/Python-3.6.1/python total_average.py

#/home/u4/u3923203/software/Python/Python-3.6.1/python read_average_v2.py


### To output in cram format
/pkg/biology/SAMtools/SAMtools_v1.10/bin/samtools view -C -T ${HG38} ${ID}_${DAY}_${Ver}_bwamem.marked.realigned.fixed.recal.indexed.bam  > ${ID}_${DAY}_${Ver}_bwamem.marked.realigned.fixed.recal.indexed.cram
/pkg/biology/SAMtools/SAMtools_v1.10/bin/samtools index ${ID}_${DAY}_${Ver}_bwamem.marked.realigned.fixed.recal.indexed.cram

### To get collect QC HsMetrics
java -Xmx40g -jar ${PICARD} CollectHsMetrics INPUT=${ID}_${DAY}_${Ver}_bwamem.marked.realigned.fixed.recal.indexed.bam OUTPUT=${ID}_${DAY}.hs_metrics.txt R=${HG38} BAIT_INTERVALS=${BED} TARGET_INTERVALS=${BED}

### Germline Variant calling by HaplotypeCaller
java  -Xmx40g -jar ${GATK} -T HaplotypeCaller  -l INFO -R ${HG38} -I ${ID}_${DAY}_${Ver}_bwamem.marked.realigned.fixed.recal.indexed.bam --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000  --dbsnp ${dbSNP} --max_alternate_alleles 30 -o ${ID}_${DAY}_${Ver}_bwamem.haplotype.SnpIndel.g.vcf.gz -nct 16 

java  -Xmx40g -jar ${GATK} -T GenotypeGVCFs -R ${HG38} -V ${ID}_${DAY}_${Ver}_bwamem.haplotype.SnpIndel.g.vcf.gz -o ${ID}_${DAY}_${Ver}_bwamem.haplotype.SnpIndel.vcf.gz 

java  -Xmx40g -jar ${GATK} -T VariantFiltration -R ${HG38} --variant ${ID}_${DAY}_${Ver}_bwamem.haplotype.SnpIndel.vcf.gz -o ${ID}_${DAY}_${Ver}_bwamem.filtered.haplotype.SnpIndel.vcf.gz --clusterWindowSize 10 --filterExpression "DP < 5" --filterName "LowCoverage" --filterExpression "QUAL < 30.0" --filterName "VeryLowQual" --filterExpression "QUAL > 30.0 && QUAL < 50.0" --filterName "LowQual" --filterExpression "QD < 1.5" --filterName "LowQD" 

### Low VAF variant calling by Mutect2
java  -Xmx40g -jar /pkg/biology/GATK/GATK_v4.1.8.0/gatk-package-4.1.8.0-local.jar Mutect2 -R ${HG38} -I ${ID}_${DAY}_${Ver}_bwamem.marked.realigned.fixed.recal.indexed.bam -O ${ID}_${DAY}_${Ver}_bwamem.Mutect2.vcf.gz 

java  -Xmx40g -jar /pkg/biology/GATK/GATK_v4.1.8.0/gatk-package-4.1.8.0-local.jar FilterMutectCalls -R ${HG38} --variant ${ID}_${DAY}_${Ver}_bwamem.Mutect2.vcf.gz -O ${ID}_${DAY}_${Ver}_bwamem.filtered.Mutect2.vcf.gz 

### Variant annotation by ANNOVAR

perl /pkg/biology/ANNOVAR/ANNOVAR_20191024/table_annovar.pl ${ID}_${DAY}_${Ver}_bwamem.filtered.haplotype.SnpIndel.vcf.gz /project/GP1/u3710062/AI_SHARE/reference/annovar_2016Feb01/humandb/ -buildver hg38 -out ${ID}_${DAY}_${Ver}_bwamem.filtered.haplotype -remove -protocol refGene,cytoBand,knownGene,ensGene,gnomad30_genome,avsnp150,gnomad211_exome,intervar_20180118,clinvar_20210123,icgc28,dbnsfp41a -operation gx,r,gx,gx,f,f,f,f,f,f,f -arg '-splicing 10',,,,,,,,,, -nastring . -vcfinput -polish --maxgenethread 20 --thread 20

perl /pkg/biology/ANNOVAR/ANNOVAR_20191024/table_annovar.pl ${ID}_${DAY}_${Ver}_bwamem.filtered.Mutect2.vcf.gz /project/GP1/u3710062/AI_SHARE/reference/annovar_2016Feb01/humandb/ -buildver hg38 -out ${ID}_${DAY}_${Ver}_bwamem_Mutect2 -remove -protocol refGene,cytoBand,knownGene,ensGene,gnomad30_genome,avsnp150,gnomad211_exome,intervar_20180118,clinvar_20210123,icgc28,dbnsfp41a -operation gx,r,gx,gx,f,f,f,f,f,f,f -arg '-splicing 10',,,,,,,,,, -nastring . -vcfinput -polish --maxgenethread 20 --thread 20

rm ${ID}_${DAY}_${Ver}_bwamem.sam ${ID}_${DAY}_${Ver}_bwamem.bam ${ID}_${DAY}_${Ver}_bwamem.marked.bam ${ID}_${DAY}_${Ver}_bwamem.marked.realigned.bam ${ID}_${DAY}_${Ver}_bwamem.marked.realigned.fixed.bam ${ID}_${DAY}_${Ver}_bwamem.marked.realigned.fixed.recal.bam ${ID}_${DAY}_${Ver}_bwamem.marked.realigned.fixed.recal.indexed.bam

################################################################ For hg19
Ver=hg19
A1_HG19=/project/GP1/u3710062/AI_SHARE/reference/hg19.NC_012920/ucsc.hg19.NC_012920.fasta
BED=/project/GP1/u3710062/AI_SHARE/shared_scripts/A1_Panel/20210304.A1_hg19.interval_list
dbSNP=/project/GP1/u3710062/AI_SHARE/reference/hg19.NC_012920/dbsnp_137.hg19.rmMT.vcf
mkdir -p ${FOLDER}/${ID}/${Ver}
cd ${FOLDER}/${ID}/${Ver}
### BWA alignment
/pkg/biology/BWA/BWA_v0.7.17/bwa mem -t 16 -R '@RG\tID:'${ID}'_'${DAY}'_'${Ver}'_bwamem\tLB:'${ID}'_'${DAY}'_'${Ver}'_bwamem\tSM:'${ID}'_'${DAY}'_'${Ver}'_bwamem\tPL:ILLUMINA\' ${A1_HG19} ${READ1} ${READ2} > ${ID}_${DAY}_${Ver}_bwamem.sam

java  -Xmx40g -jar ${PICARD} SortSam INPUT=${ID}_${DAY}_${Ver}_bwamem.sam OUTPUT=${ID}_${DAY}_${Ver}_bwamem.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true

### To output unmapped.bam
/pkg/biology/SAMtools/SAMtools_v1.10/bin/samtools view -b -f 4 ${ID}_${DAY}_${Ver}_bwamem.bam > ${ID}_${DAY}_${Ver}_bwamem.unmapped.bam

### GATK & Picard pre-process
java  -Xmx40g -jar ${PICARD} MarkDuplicates INPUT=${ID}_${DAY}_${Ver}_bwamem.bam OUTPUT=${ID}_${DAY}_${Ver}_bwamem.marked.bam METRICS_FILE=${ID}_${DAY}_${Ver}_bwamem_metrics VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true

java  -Xmx40g -jar ${GATK} -T RealignerTargetCreator -R ${A1_HG19} -I ${ID}_${DAY}_${Ver}_bwamem.marked.bam -o ${ID}_${DAY}_${Ver}_bwamem.marked.bam.intervals -nt 16

java  -Xmx40g -jar ${GATK} -T IndelRealigner -R ${A1_HG19} -I ${ID}_${DAY}_${Ver}_bwamem.marked.bam -targetIntervals ${ID}_${DAY}_${Ver}_bwamem.marked.bam.intervals -o ${ID}_${DAY}_${Ver}_bwamem.marked.realigned.bam

java  -Xmx40g -jar ${PICARD} FixMateInformation INPUT=${ID}_${DAY}_${Ver}_bwamem.marked.realigned.bam OUTPUT=${ID}_${DAY}_${Ver}_bwamem.marked.realigned.fixed.bam SO=coordinate VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true

java  -Xmx40g -jar ${GATK} -T BaseRecalibrator -R ${A1_HG19} -knownSites ${dbSNP} -I ${ID}_${DAY}_${Ver}_bwamem.marked.realigned.fixed.bam -o ${ID}_${DAY}_${Ver}_bwamem.marked.realigned.fixed.bam.recal_data.grp -rf BadCigar -nct 16

java  -Xmx40g -jar ${GATK} -T PrintReads -R ${A1_HG19} -I ${ID}_${DAY}_${Ver}_bwamem.marked.realigned.fixed.bam -BQSR ${ID}_${DAY}_${Ver}_bwamem.marked.realigned.fixed.bam.recal_data.grp -o ${ID}_${DAY}_${Ver}_bwamem.marked.realigned.fixed.recal.bam -nct 16

java  -Xmx40g -jar ${PICARD} SortSam INPUT=${ID}_${DAY}_${Ver}_bwamem.marked.realigned.fixed.recal.bam OUTPUT=${ID}_${DAY}_${Ver}_bwamem.marked.realigned.fixed.recal.indexed.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true

cp ${ID}_${DAY}_${Ver}_bwamem.marked.realigned.fixed.recal.indexed.bai ${ID}_${DAY}_${Ver}_bwamem.marked.realigned.fixed.recal.indexed.bam.bai

### Custimized scripts

#perl /home/u4/u3923203/software/BP/findSHv2_20170921.pl ${ID}_${DAY}_${Ver}_bwamem.sam ${ID}_${DAY}_${Ver}_bwamem_FindSH.txt

#perl /home/u4/u3923203/software/BP/findBPv2_20170921.pl ${ID}_${DAY}_${Ver}_bwamem_FindSH.txt ${ID}_${DAY}_${Ver}_bwamem_findBP_original.txt ${ID}_${DAY}_${Ver}_bwamem_findBP_merged.txt

#/home/u4/u3923203/software/bam-readcount-master/bam-readcount -f ${A1_HG19} -l ${BED} -w 0 ${ID}_${DAY}_${Ver}_bwamem.marked.realigned.fixed.recal.indexed.bam > ${ID}_${DAY}_${Ver}_bwamem_nucleotide_composition.txt

#cp /home/u4/u3923203/Bed/bamReadCount_rearrangement.py

#/home/u4/u3923203/software/Python/Python-3.6.1/python bamReadCount_rearrangement.py

#cp /home/u4/u3923203/Bed/HRD/HRD_v3/HRD_v3_genelist.txt

#cp /home/u4/u3923203/Bed/total_average.py

#cp /home/u4/u3923203/Bed/read_average_v2.py

#/home/u4/u3923203/software/Python/Python-3.6.1/python total_average.py

#/home/u4/u3923203/software/Python/Python-3.6.1/python read_average_v2.py

### To output in cram format
/pkg/biology/SAMtools/SAMtools_v1.10/bin/samtools view -C -T ${A1_HG19} ${ID}_${DAY}_${Ver}_bwamem.marked.realigned.fixed.recal.indexed.bam  > ${ID}_${DAY}_${Ver}_bwamem.marked.realigned.fixed.recal.indexed.cram
/pkg/biology/SAMtools/SAMtools_v1.10/bin/samtools index ${ID}_${DAY}_${Ver}_bwamem.marked.realigned.fixed.recal.indexed.cram

### To get collect QC HsMetrics
java -Xmx40g -jar ${PICARD} CollectHsMetrics INPUT=${ID}_${DAY}_${Ver}_bwamem.marked.realigned.fixed.recal.indexed.bam OUTPUT=${ID}_${DAY}.hs_metrics.txt R=${A1_HG19} BAIT_INTERVALS=${BED} TARGET_INTERVALS=${BED}

### Germline Variant calling by HaplotypeCaller
java -Xmx35g -jar ${GATK} -T HaplotypeCaller  -l INFO -R ${A1_HG19} -I ${ID}_${DAY}_${Ver}_bwamem.marked.realigned.fixed.recal.indexed.bam --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000  --dbsnp ${dbSNP} --max_alternate_alleles 30 -o ${ID}_${DAY}_${Ver}_bwamem.haplotype.SnpIndel.g.vcf.gz -nct 16

java -Xmx35g -jar ${GATK} -T GenotypeGVCFs -R ${A1_HG19} -V ${ID}_${DAY}_${Ver}_bwamem.haplotype.SnpIndel.g.vcf.gz -o ${ID}_${DAY}_${Ver}_bwamem.haplotype.SnpIndel.vcf.gz

java -Xmx35g -jar ${GATK} -T VariantFiltration -R ${A1_HG19} --variant ${ID}_${DAY}_${Ver}_bwamem.haplotype.SnpIndel.vcf.gz -o ${ID}_${DAY}_${Ver}_bwamem.filtered.haplotype.SnpIndel.vcf.gz --clusterWindowSize 10 --filterExpression "DP < 5" --filterName "LowCoverage" --filterExpression "QUAL < 30.0" --filterName "VeryLowQual" --filterExpression "QUAL > 30.0 &&  QUAL < 50.0" --filterName "LowQual" --filterExpression "QD < 1.5" --filterName "LowQD"

### Low VAF variant calling by Mutect2
java  -Xmx35g -jar /pkg/biology/GATK/GATK_v4.1.8.0/gatk-package-4.1.8.0-local.jar Mutect2 -R ${A1_HG19} -I ${ID}_${DAY}_${Ver}_bwamem.marked.realigned.fixed.recal.indexed.bam -O ${ID}_${DAY}_${Ver}_bwamem.Mutect2.vcf.gz

java  -Xmx35g -jar /pkg/biology/GATK/GATK_v4.1.8.0/gatk-package-4.1.8.0-local.jar FilterMutectCalls -R ${A1_HG19} --variant ${ID}_${DAY}_${Ver}_bwamem.Mutect2.vcf.gz -O ${ID}_${DAY}_${Ver}_bwamem.filtered.Mutect2.vcf.gz

### Variant annotation by ANNOVAR

perl /pkg/biology/ANNOVAR/ANNOVAR_20191024/table_annovar.pl ${ID}_${DAY}_${Ver}_bwamem.filtered.haplotype.SnpIndel.vcf.gz /project/GP1/u3710062/AI_SHARE/reference/annovar_2016Feb01/humandb/ -buildver hg19 -out ${ID}_${DAY}_${Ver}_bwamem_haplotype -remove -protocol refGene,cytoBand,knownGene,ensGene,gnomad211_genome,avsnp150,TaiwanBiobank-official,gnomad211_exome,TWB_1497_joing_calling_AF,intervar_20180118,clinvar_20210123,dbnsfp41a,dbscsnv11,gwava -operation gx,r,gx,gx,f,f,f,f,f,f,f,f,f,f -arg '-splicing 10',,,,,,,,,,,,, -nastring . -vcfinput -polish --maxgenethread 20 --thread 20

perl /pkg/biology/ANNOVAR/ANNOVAR_20191024/table_annovar.pl ${ID}_${DAY}_${Ver}_bwamem.filtered.Mutect2.vcf.gz /project/GP1/u3710062/AI_SHARE/reference/annovar_2016Feb01/humandb/ -buildver hg19 -out ${ID}_${DAY}_${Ver}_bwamem_Mutect2 -remove -protocol refGene,cytoBand,knownGene,ensGene,gnomad211_genome,avsnp150,TaiwanBiobank-official,gnomad211_exome,TWB_1497_joing_calling_AF,intervar_20180118,clinvar_20210123,cosmic_coding_GRCh37_v92,cosmic_noncoding_GRCh37_v92,icgc28,dbnsfp41a,cg69,kaviar_20150923,dbscsnv11,spidex,gwava -operation gx,r,gx,gx,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f -arg '-splicing 10',,,,,,,,,,,,,,,,,,, -nastring . -vcfinput -polish --maxgenethread 20 --thread 20

rm ${ID}_${DAY}_${Ver}_bwamem.sam ${ID}_${DAY}_${Ver}_bwamem.bam ${ID}_${DAY}_${Ver}_bwamem.marked.bam ${ID}_${DAY}_${Ver}_bwamem.marked.realigned.bam ${ID}_${DAY}_${Ver}_bwamem.marked.realigned.fixed.bam ${ID}_${DAY}_${Ver}_bwamem.marked.realigned.fixed.recal.bam ${ID}_${DAY}_${Ver}_bwamem.marked.realigned.fixed.recal.indexed.bam
