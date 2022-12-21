#!/usr/bin/sh
workdir=/home/u3003390/work/FASTQ
SampleList=/home/u3003390/work/FASTQ/list.text
#PIPELINE=/home/u3003390/Script/nchc_script/workflow_me/NGS.sh
PIPELINE=/home/u3003390/Script/nchc_script/workflow_me/spliceai.sh
DAY=`date +%Y%m%d`

while read -r ID; 
	do
	rsync ${PIPELINE} ${workdir}/${DAY}_${ID}.sh 
	sed -i "s|SAMPLE_ID|${ID}|g" ${workdir}/${DAY}_${ID}.sh 
    echo "create sh"
	sleep 0.5s
	done<${SampleList}

