#!/usr/bin/sh
IME=`date +%Y%m%d%H%M`
logfile=/staging/yuting0413/20221014_WGS_reJC/${TIME}_55WGS_HQ_JC_run.log
exec 3<&1 4<&2
exec >$logfile 2>&1
set -euo pipefail


dragen -f \
-r /staging/human/reference/hg38_graph_v4.0 \
--enable-joint-genotyping true \
--output-directory /staging/yuting0413/20221014_WGS_reJC \
--output-file-prefix 20221014_55WGS_TWBHQ164_dragen_graph38_JC \
--variant /staging/yuting0413/snv_vcf/[sampleID]_dragen_v4.0.3_hs38DH_graph.hard-filtered.gvcf.gz \
--variant /staging/yuting0413/snv_vcf/[sampleID]_dragen_v4.0.3_hs38DH_graph.hard-filtered.gvcf.gz \
