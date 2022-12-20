#!/usr/bin/python
#coding:utf-8
import os, sys


#============================================================#

path1 = os.path.split(os.path.realpath(__file__))[0]

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
    if i.find("rearrangment.txt") == -1:
        continue
    else:
        dirName.append(i)
#print(dirName)

# dirName2 ( 去除檔案副檔名)
#dirName3( **_R1_fastq  只留下 前面 sample name)
#dirName4(去除相同 sample name)
dirName2=[]
dirName3=[]

for i in dirName:
    a=os.path.splitext(i)[0]
    dirName2.append(a)
    

#產生 sample name file
for i in dirName2:
    
    file0 =open(path1+'/'+i+'.txt','r', encoding='utf-8')
    file1 =open(path1+'/'+i+'_total_average.txt', 'w', encoding='utf-8')

    total_read=[]

    for j in file0.readlines()[1:]:
        data = j.strip().split("\t")
        
        total_read.append(int(data[3]))
        
    sum_read=sum( total_read)
    number_read=len(total_read)

    average_read=sum_read / number_read
    
    file1.write(i+'\t'+'average'+'\t'+str(round(average_read,2))+'\n' )
    file1.close()

    

