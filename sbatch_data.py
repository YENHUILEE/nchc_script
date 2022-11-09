#!/usr/bin/python
#coding:utf-8
import os, sys


#============================================================#
memory = '16G' #4G/16G/48G/128G/128G-16/256G-32
path1 = os.path.split(os.path.realpath(__file__))[0]

#============================================================#


# Open a file
#建立路徑
# 製造資料夾(filter)
#os.mkdir(path1, mode=0o777)


#dirs (檔案名稱含副檔名)
dirs = os.listdir( path1)


#os.system('bsub'+' '+'<'+'*.sh' )


#將dirs中只含有fastq 檔案
dirName =[]
for i in dirs:
    if i.find('data.sh') == -1:
        continue
    else:
        dirName.append(i)

for i in dirName:
    os.system('sbatch'+' '+i )
