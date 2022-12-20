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

JOBDIR=""  #working directory
list_of_gvcfs=
SampleName= #output prefix

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
#fastq_folder="/staging/reserve/paylong_ntu/AI_SHARE/temp/SPB03"
#fastq_1="${fastq_folder}/${SampleName}/${SampleName}*R1*.gz" #NGS1_20170305A.R1.fastq.gz
#fastq_2="${fastq_folder}/${SampleName}/${SampleName}*R2*.gz"  #If using Illumina paired data
#sample="SM_"${SampleName}
#group="GP_"${SampleName}
platform="ILLUMINA"

# Update with the location of the reference data files
ref_dir="/staging/reserve/paylong_ntu/AI_SHARE/reference/GATK_bundle/2.8/hg19"

fasta="${ref_dir}/ucsc.hg19.fasta"
dbsnp="${ref_dir}/dbsnp_138.hg19.vcf"
known_Mills_indels="${ref_dir}/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"
known_1000G_indels="${ref_dir}/1000G_phase1.indels.hg19.sites.vcf"


# Determine whether Variant Quality Score Recalibration will be run
# VQSR should only be run when there are sufficient variants called
run_vqsr="yes"
#fastq_1="${fastq_folder}/${SampleName}/${SampleName}*R1*.gz" #NGS1_20170305A.R1.fastq.gz
#fastq_2="${fastq_folder}/${SampleName}/${SampleName}*R2*.gz"  #If using Illumina paired data
#sample="SM_"${SampleName}
#group="GP_"${SampleName}
platform="ILLUMINA"


# Update with the location of the reference data files
ref_dir="/staging/reserve/paylong_ntu/AI_SHARE/reference/GATK_bundle/2.8/hg19"

fasta="${ref_dir}/ucsc.hg19.fasta"
dbsnp="${ref_dir}/dbsnp_138.hg19.vcf"
known_Mills_indels="${ref_dir}/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"
known_1000G_indels="${ref_dir}/1000G_phase1.indels.hg19.sites.vcf"

# Determine whether Variant Quality Score Recalibration will be run
# VQSR should only be run when there are sufficient variants called
run_vqsr="yes"

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
workdir=${JOBDIR} #Determine where the output files will be stored

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

# GenotypeGVCFs & joint calling
gvcf_argument=""
while read -r line; do
 gvcf_argument=$gvcf_argument" -v $line"
done < ${list_of_gvcfs}

$release_dir/bin/sentieon driver -r $fasta -t $nt --algo GVCFtyper ${SampleName}.output-hc.vcf.gz $gvcf_argument

# ******************************************
# 7b. Variant Recalibration
# ******************************************
if [ "$run_vqsr" = "yes" ]; then
  #for SNP
  #create the resource argument
  resource_text="--resource $vqsr_1000G_phase1 --resource_param 1000G,known=false,training=true,truth=false,prior=10.0 "
  resource_text="$resource_text --resource $vqsr_1000G_omni --resource_param omni,known=false,training=true,truth=true,prior=12.0 "
  resource_text="$resource_text --resource $vqsr_dbsnp --resource_param dbsnp,known=true,training=false,truth=false,prior=2.0 "
  resource_text="$resource_text --resource $vqsr_hapmap --resource_param hapmap,known=false,training=true,truth=true,prior=15.0"

  #create the annotation argument
  annotation_array="QD MQ MQRankSum ReadPosRankSum FS"

  #Initial annotate_text variable
  annotate_text=""
  for annotation in $annotation_array; do
    annotate_text="$annotate_text --annotation $annotation"
  done

  #Run the VQSR
  tranches="--tranche 100.0 --tranche 99.9 --tranche 99.0 --tranche 98.0 --tranche 97.0 --tranche 96.0 --tranche 95.0 --tranche 94.0 --tranche 93.0 --tranche 92.0 --tranche 91.0 --tranche 90.0"
  $release_dir/bin/sentieon driver -r $fasta -t $nt --algo VarCal -v ${SampleName}.output-hc.vcf.gz $resource_text $annotate_text --var_type SNP --plot_file ${SampleName}.vqsr_SNP.hc.plot_file.txt --max_gaussians 8 --srand 47382911 --tranches_file ${SampleName}.vqsr_SNP.hc.tranches ${SampleName}.vqsr_SNP.hc.recal $tranches

  #apply the VQSR
  $release_dir/bin/sentieon driver -r $fasta -t $nt --algo ApplyVarCal -v ${SampleName}.output-hc.vcf.gz --var_type SNP --recal ${SampleName}.vqsr_SNP.hc.recal --tranches_file ${SampleName}.vqsr_SNP.hc.tranches --sensitivity 99.5 ${SampleName}.vqsr_SNP.hc.recaled.vcf.gz

  #plot the report
  $release_dir/bin/sentieon plot vqsr -o ${SampleName}.vqsr_SNP.VQSR.pdf ${SampleName}.vqsr_SNP.hc.plot_file.txt

  #for indels after SNPs
  #create the resource argument
  resource_text="--resource $vqsr_1000G_phase1_indel --resource_param 1000G,known=false,training=true,truth=false,prior=10.0 "
  resource_text="$resource_text --resource $vqsr_Mill --resource_param Mills,known=false,training=true,truth=true,prior=12.0 "
  resource_text="$resource_text --resource $vqsr_dbsnp --resource_param dbsnp,known=true,training=false,truth=false,prior=2.0 "

  #create the annotation argument
  annotation_array="QD MQ ReadPosRankSum FS"
  annotate_text=""
  for annotation in $annotation_array; do
    annotate_text="$annotate_text --annotation $annotation"
  done

  #Run the VQSR
  tranches="--tranche 100.0 --tranche 99.9 --tranche 99.0 --tranche 98.0 --tranche 97.0 --tranche 96.0 --tranche 95.0 --tranche 94.0 --tranche 93.0 --tranche 92.0 --tranche 91.0 --tranche 90.0"
  $release_dir/bin/sentieon driver -r $fasta -t $nt --algo VarCal -v ${SampleName}.vqsr_SNP.hc.recaled.vcf.gz $resource_text $annotate_text --var_type INDEL --plot_file ${SampleName}.vqsr_SNP_INDEL.hc.plot_file.txt --max_gaussians 4 --srand 47382911 --tranches_file ${SampleName}.vqsr_SNP_INDEL.hc.tranches ${SampleName}.vqsr_SNP_INDEL.hc.recal $tranches

  #apply the VQSR
  tranches="--tranche 100.0 --tranche 99.9 --tranche 99.0 --tranche 98.0 --tranche 97.0 --tranche 96.0 --tranche 95.0 --tranche 94.0 --tranche 93.0 --tranche 92.0 --tranche 91.0 --tranche 90.0"
  $release_dir/bin/sentieon driver -r $fasta -t $nt --algo VarCal -v ${SampleName}.vqsr_SNP.hc.recaled.vcf.gz $resource_text $annotate_text --var_type INDEL --plot_file ${SampleName}.vqsr_SNP_INDEL.hc.plot_file.txt --max_gaussians 4 --srand 47382911 --tranches_file ${SampleName}.vqsr_SNP_INDEL.hc.tranches ${SampleName}.vqsr_SNP_INDEL.hc.recal $tranches

  #apply the VQSR
  $release_dir/bin/sentieon driver -r $fasta -t $nt --algo ApplyVarCal -v ${SampleName}.vqsr_SNP.hc.recaled.vcf.gz --var_type INDEL --recal ${SampleName}.vqsr_SNP_INDEL.hc.recal --tranches_file ${SampleName}.vqsr_SNP_INDEL.hc.tranches --sensitivity 99.5 ${SampleName}.vqsr_SNP_INDEL.hc.recaled.vcf.gz

  #plot the report
  $release_dir/bin/sentieon plot vqsr -o ${SampleName}.vqsr_SNP_INDEL.VQSR.pdf ${SampleName}.vqsr_SNP_INDEL.hc.plot_file.txt
fi


set +x
exec >&3 2>&4
exec 3<&- 4<&-

printf "#############################################################################\n"
printf "###                  Work completed: $(date +%Y-%m-%d:%H:%M:%S)                  ###\n"
printf "#############################################################################\n"
