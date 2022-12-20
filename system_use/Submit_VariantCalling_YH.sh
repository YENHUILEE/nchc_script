#!/usr/bin/sh
SampleList =/home/u3003390/work/FASTQ/list.text
PIPELINE=/home/u3003390/Script/nchc_script/workflow_me/NGS.sh
DAY=`date +%Y%m%d`

while read -r ID; 
	do
	rsync ${PIPELINE} ./${DAY}_${PIPELINE}_${ID}.sh 
	sed -i "s|SAMPLE_NAME|${ID}|g" ./${DAY}_${PIPELINE}_${ID}.sh 
	sleep 3s
	done<${SampleList}
