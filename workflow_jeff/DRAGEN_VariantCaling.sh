#!/usr/bin/bash
Samplelist=/staging/yuting0413/sample.list

while read -r sample
do
mkdir -p ${sample}
output_dir=/staging/yuting0413/${sample}
r1=/staging/yuting0413/${sample}*R1*.fastq.gz
r2=/staging/yuting0413/${sample}*R2*.fastq.gz
ref_dir="/staging/human/reference/hg38_graph_v4.0"
output_prefix=${sample}_dragen_v4.0.3_hs38DH_graph
RGID=@A00361
RGSM=${sample}
TIME=`date +%Y%m%d%H%M`
logfile=./${sample}/${TIME}_${sample}_run.log
exec 3<&1 4<&2
exec >$logfile 2>&1
set -euo pipefail
# 20220816
dragen --output-dir ${output_dir} --output-file-prefix ${output_prefix} -1 ${r1} -2 ${r2} --ref-dir ${ref_dir} --enable-map-align-output true --RGID ${RGID} --RGSM ${RGSM} --output-format CRAM --enable-variant-caller true --vc-enable-vcf-output true --vc-emit-ref-confidence GVCF --enable-sv true --enable-cnv true --cnv-enable-ref-calls true --cnv-enable-self-normalization true --cnv-enable-gcbias-correction true --repeat-genotype-enable true --enable-duplicate-marking true --remove-duplicates true --enable-hla true --enable-cyp2d6 true --enable-cyp2b6 true --enable-gba true --enable-star-allele true --enable-smn true --enable-pgx true

done <${Samplelist}
