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



#將dirs中只含有rearrangment.txt 檔案
dirName =[]
for i in dirs:
    if i.find("rearrangment.txt") == -1:
        continue
    else:
        dirName.append(i)

#=============================================================#        
dirName2=[]


for i in dirName:
    a=os.path.splitext(i)[0]
    dirName2.append(a)
#產生 sample name file

#=============================================================#


dirName3 =[]
for i in dirs:
    if i.find("genelist.txt") == -1:
        continue
    else:
        dirName3.append(i)



#=============================================================#
#file1 = genelist.txt(chr1	11994723	11994912	PLOD1_NM_000302_exon1)
#file0=**nucleotide_composition_rearrangment.txt
#(chr	position	nucleotide	total_count	A	C   	G	T	N) 

for l in dirName3:
    for z in dirName2:
        
        file2 =open(path1+'/'+z+'_read_average.txt', 'w', encoding='utf-8')
        file1 =open(path1+'/'+l,'r', encoding='utf-8')
        file2.write('target'+'\t'+'chr'+'\t'+'start'+'\t'+'end'+'\t'+'read_average'+'\n' )
        for j in file1.readlines():
            line = j.strip().split("\t")
            total_read=[]
            for k in dirName2:
                file0 =open(path1+'/'+k+'.txt','r', encoding='utf-8')
                for g in file0.readlines()[1:]:
                    data = g.strip().split("\t")

                    if data[0] in  line[0] :
                        if int(data[1]) in  range ( int(line[1]),int(line[2])+1) :
                            total_read.append(int(data[3]))
                            

                file0.close()

            sum_read=sum( total_read)
            number_read = len(total_read)


            if sum_read == 0 :
                file2.write(line[3]+'\t'+'0'+'\n' )

            else :
                average_read=sum_read / number_read
                file2.write( line[3]+'\t'+line[0]+'\t'+line[1]+'\t'+line[2]+'\t'+str(round(average_read,2))+'\n' )


