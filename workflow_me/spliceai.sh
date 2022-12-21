#!/bin/bash
#SBATCH -p ngs48G
#SBATCH -c 14
#SBATCH --mem=46g
#SBATCH -A MST109178
#SBATCH -J SAMPLE_ID
#SBATCH -e /work/u3003390/err 
#SBATCH -o /work/u3003390/out
#SBATCH --mail-user=k32650805@gmail.com
#SBATCH --mail-type=FAIL,END

# This code is used for NGS germline variant calling (SNV + Indels)

echo "processes directories"
# before start, you should copy the fastq to fastq_folder
# fold to be set
# change the following fold as you needed
sample_id=SAMPLE_ID 
user_dir="/work/u3003390"
fastq_dir="/work/u3003390/FASTQ" # fold of fastq
temp_dir="/work/u3003390/TEMP" # fold for temporary file
release_dir="/work/u3003390/RESULT" # fold for major result file
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


#splice ai
module load biology/Tensorflow/2.7.1
module load biology/SpliceAI/1.3

#load python (**must be the last)
module load biology/Python/3.9.5 


spliceai \
-I ${release_dir}/${sample_id}.${Date}.annotate.hg19_multianno.vcf \
-O ${release_dir}/${sample_id}.${Date}.annotate.hg19_multianno.spliceai.vcf \
-R ${fasta} \
-A grch37
#
