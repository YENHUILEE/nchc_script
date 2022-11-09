#!/usr/bin/python
#coding:utf-8
import os, sys

#==========================================================#
path = os.path.split(os.path.realpath(__file__))[0]
#==========================================================#
#dirs (檔案名稱含副檔名)
dirs = os.listdir( path )

#將dirs中只含有nucleotide_composition.txt 檔案
dirName =[]
for i in dirs:
    if i.find("nucleotide_composition.txt") == -1:
        continue
    else:
        dirName.append(i)

# dirName2 ( 去除檔案副檔名)
#dirName3( **nucleotide_composition.txt  只留下 前面 sample name)
#dirName4(去除相同 sample name)
dirName2=[]

for i in dirName:
    a=os.path.splitext(i)[0]
    
    dirName2.append(a)
    

for j in dirName2:
    file1 = open(path+'/'+j+'.txt',"r",encoding = 'utf-8')
    
    file2 = open(path+'/'+j+"_rearrangment.txt","w",encoding = 'utf-8')
    file2.write("chr"+"\t"+"position"+"\t"+"nucleotide"+"\t"+"total_count"+"\t"+"A"+"\t"+"C"+"\t"+"G"+"\t"+"T"+"\t"+"N"+"\t"+"\n")

    for i in file1.readlines():
        new_rearrange = []
        old_data_list = i.strip().split("\t")
        new_rearrange.append(old_data_list[0])
        new_rearrange.append(old_data_list[1])
        new_rearrange.append(old_data_list[2])
        new_rearrange.append(old_data_list[3])
        A_count_list = old_data_list[5].strip().split(":")
        A_count = str(A_count_list[1])
        new_rearrange.append(A_count)
        C_count_list = old_data_list[6].strip().split(":")
        C_count = str(C_count_list[1])
        new_rearrange.append(C_count)
        G_count_list = old_data_list[7].strip().split(":")
        G_count = str(G_count_list[1])
        new_rearrange.append(G_count)
        T_count_list = old_data_list[8].strip().split(":")
        T_count = str(T_count_list[1])
        new_rearrange.append(T_count)
        N_count_list = old_data_list[9].strip().split(":")
        N_count = str(N_count_list[1])
        new_rearrange.append(N_count)
        write_file_data = "\t".join(new_rearrange)
        file2.write(write_file_data+"\n")
    file1.close()
    file2.close()

print("done")
