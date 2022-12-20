#!/usr/bin/sh
workdir=/home/u3003390/work/FASTQ
SampleList=/home/u3003390/work/FASTQ/list.text
PIPELINE=/home/u3003390/Script/nchc_script/workflow_me/NGS.sh
DAY=`date +%Y%m%d`

while read -r ID; 
	do
	rsync ${PIPELINE} ${workdir}/${DAY}_${ID}.sh 
	sed -i "s|SAMPLE_NAME|${ID}|g" ${workdir}/${DAY}_${ID}.sh 
	sleep 3s
	done<${SampleList}
