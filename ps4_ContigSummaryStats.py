# -*- coding: utf-8 -*-
"""
PS4 De Novo Transcriptome Assembly

Bi623

August 30, 2016 

@author: jsmith

Summary info on De Novo transcript assembly using Trinity 

"""

#Using Python regular expressions
import re 

#read in the fasta files from the trinity assembly
#contigFile = "/home12/jsmith16/bi623/160829_PS4/block2/trinityNorm/pipeNorm_trinity.fasta"
contigFile = "/Users/jsmith/Documents/bi623/160829_class11_PS4/pipeNorm_trinity.fasta"
contig = open(contigFile, "r")

#define a function to identify the mean contig length 
def median(List):    
    List.sort()        #sort the list of values  
    lenArray = len(List)      #find the total number of entries in the list
    if  lenArray % 2 == 1: #if the number of entries is an odd number
            position1 = int(((lenArray - 1)/2)) #median is the number in the middle of the list
            median = List[position1]             
            return median   
    else: #if the number of entries in the list is even 
            position1 = int((lenArray/2) - 1)              
            position2 = int((lenArray/2))            
            median = float(((List[position1] + List[position2]) / 2)) #the median is the average of the two middle positions in the list
            return median              
            
#initialize empty arrays for the length of kmers
lengthArray = []


#extract k-mer length of each contig (in red below). In addition, extract the k-mer coverage for the contig
counter = 0
for line in contig:    
    line = line.strip()
    if line[0] == ">":   #use this if statement to parse out header lines in the python script  
        counter += 1         
        header = re.split(" ", line) #this splits each line of the contigHeader file into seperate arrays    
        length = header[1]
        length = int(length[4:]) #the length is the 5th charated to the end because the header format is  len=000 
        lengthArray.append(length)

  
#Calculate the number of contigs  
numContigs = counter
print("the total number of contigs in this dataset is", numContigs)

#the maximum contig length
maxContigLen = max(lengthArray)
print("the maximum length of contigs in this dataset is", maxContigLen)

#find the total length of all contrifgs in the assmebly
totalLength = sum(lengthArray)

# the mean contig length
meanContigLen = totalLength / numContigs
print("The mean contig lenght is", meanContigLen) 

#Use the median function to find the median contig length 
medContig = median(lengthArray)
print("The median contig length is", medContig)




