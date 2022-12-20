#!/usr/bin/python
#coding:utf-8
import os, sys


#============================================================#
memory = 'ngs96G' #48G/96G/192G/
path = "../BCR"
path1 = "/staging/reserve/paylong_ntu/Research/DF/DF_20220923_Research/"
data="_20220923"
tool = '_all'
align = '_bwamem'
threading ='16'
#=====================bed======================================#

bed ="/staging/reserve/paylong_ntu/AI_SHARE/Bed/DF/DF_v3/DF_v3.bed"
genelist ='/staging/reserve/paylong_ntu/AI_SHARE/Bed/DF/DF_v3/DF_v3_genelist.txt'

#==================SC(script)=================================#
#sh_gene_v1: origin
#sh_gene_v2:
# 1. modification of BWA    ( remove -M )  ( keep the splite read )
# 2. add  tool ( findSH and findBP  )
# 3. add tool ( bam_read )
#sh_gene_v2:
# 1. add gnomad in database
# 2. add  tools ( findBPv2_20170921.pl  and findSHv2_20170921.pl )

#============================================================#
### software direction
bwa="/home/u1151339/software/bwa-0.7.12/bwa"
bowtie2 ="/home/u1151339/software/bowtie2-2.2.5/bowtie2"
jar = 'java  -Xmx80g -jar'; #40g/80g
picard ="/home/u1151339/software/picard-tools-1.134/picard.jar"
GATK ="/home/u1151339/software/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar"
GATK4 ="/home/u1151339/software/GATK_v4.1.8.0/gatk-package-4.1.8.0-local.jar"
#python2="/home/u00scs00/software/Python-2.7.8/python"
python3="/home/u1151339/software/Python/Python-3.6.1/python"
#Platypus ="/home/u00scs00/software/Platypus_0.8.1/Platypus.py"
#pindel_1="source /pkg/biology/pindel/pindel_v0.2.5b8/env.sh"
#pindel="/pkg/biology/pindel/pindel_v0.2.5b8/pindel"
annovar ="/work/opt/ohpc/Taiwania3/pkg/biology/ANNOVAR/annovar_20210819/"#原20191024_新20210819_20200608
humandb="/staging/reserve/paylong_ntu/AI_SHARE/reference/annovar_2016Feb01/humandb/"
ANNOVAR = "/work/opt/ohpc/Taiwania3/pkg/biology/ANNOVAR/annovar_20210819/table_annovar.pl"#原20191024_新20210819_20200608
samtools_path="/home/u1151339/software/samtools-1.2/"
samtools="/home/u1151339/software/samtools-1.2/bin/samtools"
ref_ucsc="/home/u1151339/reference/ucsc.hg19.NC_012920.fasta"
ref_bowtie2="/home/u1151339/reference/ucsc.hg19"
snp="/home/u1151339/reference/dbsnp_137.hg19.rmMT.vcf"
findSH = "/home/u1151339/software/BP/findSHv2_20170921.pl"
findBP = "/home/u1151339/software/BP/findBPv2_20170921.pl"
bam_read="/home/u1151339/software/bam-readcount-master/bam-readcount"
bam_readcount="/staging/reserve/paylong_ntu/AI_SHARE/Bed/bamReadCount_rearrangement.py"
total_average = "/staging/reserve/paylong_ntu/AI_SHARE/Bed/total_average.py"
read_average_v2 = "/staging/reserve/paylong_ntu/AI_SHARE/Bed/read_average_v2.py"

#=============================================================#

# Open a file
#建立路徑
# 製造資料夾(filter)
#os.mkdir(path1, mode=0o777)


#dirs (檔案名稱含副檔名)
dirs = os.listdir( path1 )



#將dirs中只含有fastq 檔案
dirName =[]
for i in dirs:
    if i.find("fastq.gz") == -1:
        continue
    else:
        dirName.append(i)

# dirName2 ( 去除檔案副檔名)
#dirName3( **_R1_fastq  只留下 前面 sample name)
#dirName4(去除相同 sample name)
dirName2=[]
dirName3=[]

for i in dirName:
    a=os.path.splitext(i)[0]
    
    b=a.index("_")
    
    c=a[0:b]
    
    dirName2.append(a)
    dirName3.append(c)
    
dirName4=list(set(dirName3))

#print(dirName2)
#print(dirName3)
#產生 sample name file
for i in dirName4:
    os.system('mkdir'+' '+i)
    file1 =open(path1+i+data+align+tool+'.sh', 'w', encoding='utf-8')


    
#========================================================#

#生成 sh
    file1.write( '#!/bin/bash'+'\n' )#使用的shell名稱
    file1.write('#SBATCH -p '+ memory +'\n')#queue_name=partition
    file1.write('#SBATCH -c 28'+'\n')#48=14/96=28Core數_參考queue設定
    file1.write('#SBATCH --mem=92g'+'\n')#48=46g/96=92g記憶體量_參考queue設定
    file1.write('#SBATCH -A MST109178' +'\n')#project_ID
    file1.write('#SBATCH -J '+i+'\n' )#job_name
    file1.write('#SBATCH -o '+path1+i+'/'+i+data+align+tool+'_out.txt'+'\n' )
    file1.write('#SBATCH -e '+path1+i+'/'+i+data+align+tool+'_err.txt'+'\n' )
    file1.write('#SBATCH --mail-user=sharonchiu0104@gmail.com'+'\n' )#E-mail_address
    file1.write('#SBATCH --mail-type=FAIL,END'+'\n')#mail_events
    file1.write('\n' )
###threading  ###
    
###cp data from home to all individual path ###
    fastq_r1 = path1+i+'_*R1*.fastq.gz'
    fastq_r2 =path1+i+'_*R2*.fastq.gz'
    
### bwa_mem ###
    file1.write(bwa+' mem -t '+threading+' -R '+'\'@RG\\tID:'+i+data+align+'\\tLB:'+i+data+align+'\\tSM:'+i+data+align+'\\tPL:ILLUMINA\\\' '+
                ref_ucsc+' '+fastq_r1+' '+fastq_r2+' > '+path1+i+'/'+i+data+align+'.sam'+' '+'&&'+'\n' )
    file1.write('\n' )
### picard_SortSam ### 
    file1.write(jar+' '+picard+' ' +'SortSam INPUT='+path1+i+'/'+i+data+align+'.sam'+' '+'OUTPUT='+path1+i+'/'+i+data+align+'.bam'+
                ' '+'SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true'+' '+'&&'+'\n') 
    file1.write('\n' )
    
### picard_MarkDuplicates ###
    file1.write(jar+' '+picard+' '+'MarkDuplicates'+' '+'INPUT='+path1+i+'/'+i+data+align+'.bam'+' '+'OUTPUT='+path1+i+'/'+i+data+align+'.marked.bam'+' '+'METRICS_FILE='+path1+i+'/'+i+data+align+'_metrics'+' '+'VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true'+' '+'&&'+'\n')
    file1.write('\n' )

### GATK_RealignerTargetCreator ###
    file1.write(jar+' '+GATK+' '+'-T RealignerTargetCreator -R'+' '+ref_ucsc+' '+'-I'+' '+path1+i+'/'+i+data+align+'.marked.bam'+' '+'-o'+' '+path1+i+'/'+i+data+align+'.marked.bam.intervals'+' '+'-nt'+' '+threading+' '+'&&'+'\n')
    file1.write('\n' )

### GATK_IndelRealigner ###
    file1.write(jar+' '+GATK+' '+'-T IndelRealigner -R'+' '+ref_ucsc+' '+'-I'+' '+path1+i+'/'+i+data+align+'.marked.bam'+' '+'-targetIntervals'+' '+path1+i+'/'+i+data+align+'.marked.bam.intervals'+' '+'-o'+' '+path1+i+'/'+i+data+align+'.marked.realigned.bam'+' '+'&&'+'\n')
    file1.write('\n' )

### picard_FixMateInformation ###
    file1.write(jar+' '+picard+' '+'FixMateInformation'+' '+'INPUT='+path1+i+'/'+i+data+align+'.marked.realigned.bam'+' '+'OUTPUT='+path1+i+'/'+i+data+align+'.marked.realigned.fixed.bam'+' '+'SO=coordinate VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true'+' '+'&&'+'\n')
    file1.write('\n' )

### GATK_BaseRecalibrator ###
    file1.write(jar+' '+GATK+' '+'-T BaseRecalibrator -R'+' '+ref_ucsc+' '+'-knownSites'+' '+snp+' '+'-I'+' '+path1+i+'/'+i+data+align+'.marked.realigned.fixed.bam'+' '+'-o'+' '+path1+i+'/'+i+data+align+'.marked.realigned.fixed.bam.recal_data.grp'+' '+'-rf BadCigar'+' '+'-nct'+' '+threading+' '+'&&'+'\n')
    file1.write('\n' )

### GATK_PrintReads ###
    file1.write(jar+' '+GATK+' '+'-T PrintReads -R'+' '+ref_ucsc+' '+'-I'+' '+path1+i+'/'+i+data+align+'.marked.realigned.fixed.bam'+' '+'-BQSR'+' '+path1+i+'/'+i+data+align+'.marked.realigned.fixed.bam.recal_data.grp'+' '+'-o'+' '+path1+i+'/'+i+data+align+'.marked.realigned.fixed.recal.bam'+' '+'-nct'+' '+threading+' '+'&&'+'\n')
    file1.write('\n' )

### sh09_picard_SortSam ###
    file1.write(jar+' '+picard+' '+'SortSam'+' '+'INPUT='+path1+i+'/'+i+data+align+'.marked.realigned.fixed.recal.bam'+' '+'OUTPUT='+path1+i+'/'+i+data+align+'.marked.realigned.fixed.recal.indexed.bam'+' '+'SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true'+' '+'&&'+'\n')
    file1.write('\n' )
    
### Pindel ###
    file1.write('cp'+' '+path1+i+'/'+i+data+align+'.marked.realigned.fixed.recal.indexed.bai'+' '+path1+i+'/'+i+data+align+'.marked.realigned.fixed.recal.indexed.bam.bai'+'\n' )
    file1.write('\n' )

### FindBP ###
    file1.write('perl'+' '+findSH+' '+path1+i+'/'+i+data+align+'.sam'+' '+path1+i+'/'+i+data+align+'_FindSH.txt'+ '\n')
    file1.write('\n' )
    file1.write('perl'+' '+findBP +' '+path1+i+'/'+i+data+align+'_FindSH.txt'+' '+path1+i+'/'+i+data+align+'_findBP_original.txt'+' '+path1+i+'/'+i+data+align+'_findBP_merged.txt'+'\n')
    file1.write('\n' )

### bam_readcount ###
    file1.write(bam_read+' '+'-f'+' '+ref_ucsc+' '+'-l'+' '+bed+' '+'-w'+' '+'0'+' '+path1+i+'/'+i+data+align+'.marked.realigned.fixed.recal.indexed.bam'+' '+'>'+' '+path1+i+'/'+i+data+align+'_nucleotide_composition.txt'+'\n')
    file1.write('\n' )
    file1.write('cp'+' '+bam_readcount+' '+path1+i+'/'+'\n' )
    file1.write('\n' )
    file1.write(python3 +' '+path1+i+'/'+'bamReadCount_rearrangement.py'+'\n'  )
    file1.write('\n' )
    
### read_average ###
    file1.write('cp'+' '+genelist+' '+path1+i+'/'+'\n' )
    file1.write('\n' )
    file1.write('cp'+' '+total_average +' '+path1+i+'/'+'\n' )
    file1.write('\n' )
    file1.write('cp'+' '+read_average_v2 +' '+path1+i+'/'+'\n' )
    file1.write('\n' )
    file1.write(python3 +' '+path1+i+'/'+'total_average.py'+'\n'  )
    file1.write('\n' )
    file1.write(python3 +' '+path1+i+'/'+'read_average_v2.py'+'\n'  )
    file1.write('\n' )

### GATK_HaplotypeCaller_gvcf ###
    file1.write(jar+' '+GATK+' '+'-T HaplotypeCaller  -l INFO -R'+' '+ref_ucsc+' '+'-I'+' '+path1+i+'/'+i+data+align+'.marked.realigned.fixed.recal.indexed.bam'+' '+'--emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000'+' '
                +' '+'--dbsnp'+' '+snp+' '+'--max_alternate_alleles 30'+' '+'-o'+' '+path1+i+'/'+i+data+align+'.haplotype.SnpIndel.g.vcf.gz'+' '+'-nct'+' '+threading+' '+'\n')
    file1.write('\n' )

### GATK_HaplotypeCaller ###
    file1.write(jar+' '+GATK+' '+'-T GenotypeGVCFs -R'+' '+ref_ucsc+' '+'-V'+' '+path1+i+'/'+i+data+align+'.haplotype.SnpIndel.g.vcf.gz'+' '+'-o'+' '+path1+i+'/'+i+data+align+'.haplotype.SnpIndel.vcf.gz'+' '+'\n')
    file1.write('\n' )
    
### sh11-1_GATK_VariantFiltration ###
    file1.write(jar+' '+GATK+' '+'-T VariantFiltration -R'+' '+ref_ucsc+' '+'--variant'+' '+path1+i+'/'+i+data+align+'.haplotype.SnpIndel.vcf.gz'+' '+'-o'+' '+path1+i+'/'+i+data+align+'.filtered.haplotype.SnpIndel.vcf.gz'+' '+'--clusterWindowSize 10 --filterExpression "DP < 5" --filterName "LowCoverage" --filterExpression "QUAL < 30.0" --filterName "VeryLowQual" --filterExpression "QUAL > 30.0 && QUAL < 50.0" --filterName "LowQual" --filterExpression "QD < 1.5" --filterName "LowQD"'+' '+'\n' ) 
    file1.write('\n' )
    
### annovar_Table_annovar_HaplotypeCaller ###
    file1.write('perl'+' '+ANNOVAR+' '+path1+i+'/'+i+data+align+'.filtered.haplotype.SnpIndel.vcf.gz'+' '+humandb+' '+'-buildver hg19 -out'+' '+path1+i+'/'+i+data+align+'_haplotype'+' '+'-remove'+' '+'-protocol'+' '+'refGene,cytoBand,knownGene,ensGene,gnomad211_genome,avsnp150,TaiwanBiobank-official,gnomad211_exome,TWB_1497_joing_calling_AF,intervar_20180118,clinvar_20210501,cosmic_coding_GRCh37_v92,cosmic_noncoding_GRCh37_v92,icgc28,dbnsfp41a,cg69,kaviar_20150923,dbscsnv11,spidex,gwava,wgRna,targetScanS -operation gx,r,gx,gx,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f -arg \'-splicing 10\',,,,,,,,,,,,,,,,,,,,, -nastring . -vcfinput -polish --maxgenethread 20 --thread 20'+' '+'\n' )
    file1.write('\n' )

### GATK_Mutect2 ###
    file1.write(jar+' '+GATK4+' '+'Mutect2 -R'+' '+ref_ucsc+' '+'-I'+' '+path1+i+'/'+i+data+align+'.marked.realigned.fixed.recal.indexed.bam'+' '+'-O'+' '+path1+i+'/'+i+data+align+'.Mutect2.vcf.gz'+' '+'\n')
    file1.write('\n' )
    
## sh11-1_GATK_FilterMutectCalls ###
    file1.write(jar+' '+GATK4+' '+'FilterMutectCalls -R'+' '+ref_ucsc+' '+'--variant'+' '+path1+i+'/'+i+data+align+'.Mutect2.vcf.gz'+' '+'-O'+' '+path1+i+'/'+i+data+align+'.filtered.Mutect2.vcf.gz'+' '+'\n' ) 
    file1.write('\n' )

### annovar_Table_annovar_Mutect2 ###
    file1.write('perl'+' '+ANNOVAR+' '+path1+i+'/'+i+data+align+'.filtered.Mutect2.vcf.gz'+' '+humandb+' '+'-buildver hg19 -out'+' '+path1+i+'/'+i+data+align+'_Mutect2'+' '+'-remove'+' '+'-protocol'+' '+'refGene,cytoBand,knownGene,ensGene,gnomad211_genome,avsnp150,TaiwanBiobank-official,gnomad211_exome,TWB_1497_joing_calling_AF,intervar_20180118,clinvar_20210501,cosmic_coding_GRCh37_v92,cosmic_noncoding_GRCh37_v92,icgc28,dbnsfp41a,cg69,kaviar_20150923,dbscsnv11,spidex,gwava,wgRna,targetScanS -operation gx,r,gx,gx,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f -arg \'-splicing 10\',,,,,,,,,,,,,,,,,,,,, -nastring . -vcfinput -polish --maxgenethread 20 --thread 20'+' '+'\n' )
    file1.write('\n' )    
    file1.close()
