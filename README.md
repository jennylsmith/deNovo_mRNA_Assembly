# deNovo_mRNA_Assembly

 PS4 De Novo Transcriptome Assembly 
 August 29, 2016


 Block 1


1. Preprocessing reads of a fastq file. Trimming raw reads using process_shortreads 
from the Stacks module.
 
 Files used:
 /research/bi610/data/rnaseq/SscoPE_R1.fastq 
 /research/bi610/data/rnaseq/SscoPE_R2.fastq  

PBS Job script for trimming of adaptor sequences:
 !/usr/bin/env bash
 PBS -N process_shortReads
 PBS -q fatnodes
 PBS -l nodes=1:ppn=4 
 PBS -d /home12/jsmith16/bi623/160829_PS4

cd /home12/jsmith16/bi623/160829_PS4

module load stacks

process_shortreads -1 SscoPE_R1.fastq -2 SscoPE_R2.fastq \
 -i fastq -o /home12/jsmith16/bi623/160829_PS4/trimmedFastq -y fastq -c -q --adapter_mm 2 \
 --adapter_1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC --adapter_2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT 



2. We want to quantify how much of our data was trimmed. Using a single UNIX command 
(including pipes), calculate the distribution of read lengths (one for R1s, one for R2s) 
in your remaining data. Plot these distributions in R and turn in with the rest of the
 assignment.
 
  unix command to create distribution of read lengths. 
 cat SscoPE_R2.2.fq | grep -A 1 "^@2_" | grep -v "\-" | grep -v "^@" | awk '{print length($0)}' \
 | sort -n | uniq -c | sort -nk 1 > R2_LengthDist.txt
 
 cat SscoPE_R2.2.fq | grep -A 1 "^@2_" | grep -v "\-" | grep -v "^@" | awk '{print length($0)}' \
 | sort -n | uniq -c | sort -nk 1 > R1_LengthDist.txt



3. Filter the clean data for rare k-mers using kmer_filter,also part of the Stacks package. 
Perform rare k-mer filtering ONLY, and rename the new files SscoPEcleanfil_R1.fastq and 
SscoPEcleanfil_R2.fastq. Use the defaults for k-mer size and --min_lim. Explain what this 
is doing to your data. 


PBS Job script for filtering rare kmers:
 !/usr/bin/env bash
 PBS -N kmer_filter_rare
 PBS -q fatnodes
 PBS -l nodes=1:ppn=4 
 PBS -d /home12/jsmith16/bi623/160829_PS4

cd /home12/jsmith16/bi623/160829_PS4/trimmedFastq

module load stacks

kmer_filter -1 SscoPE_R1.1.fq -2 SscoPE_R2.2.fq -i fastq \
-o /home12/jsmith16/bi623/160829_PS4/filterFastq -y fastq --rare
 

The rare kmer filter from the Stacks module will kmerize the reads, with a default
size of a 15mer. The reads are  discarded if the rare kmers are consecutive and at least
80% of kmers fall below a median kmer coverage, which is the default minimum 
limit set by stacks. It is important to remove rare kmers, because these can represent 
sequencing errors in the reads. It is important to have the minimum limit set for 
transcriptomic datasets, since there are some biologically rare transcripts that should not 
be discarded. 



4. Runkmer_filter again via the ACISS queueing system,to normalize the cleaned,
 filtered data to 2 different coverages (40x and 20x). Use the defaul kmer size of 15. 
 Briefly explain in a sentence or two what this is doing to your data.

PBS Job script for kmer normalization: 
 !/usr/bin/env bash
 PBS -N kmer_filter_cov
 PBS -q fatnodes
 PBS -l nodes=1:ppn=4 
 PBS -d /home12/jsmith16/bi623/160829_PS4

cd /home12/jsmith16/bi623/160829_PS4/filterFastq

module load stacks

 kmer_filter -1 SscoPE_R1.1.fil.fq -2 SscoPE_R2.2.fil.fq -i fastq \
 -o /home12/jsmith16/bi623/160829_PS4/filterFastqCov/40x \
 -y fastq --normalize 40


kmer_filter -1 SscoPE_R1.1.fil.fq -2 SscoPE_R2.2.fil.fq -i fastq \
-o /home12/jsmith16/bi623/160829_PS4/filterFastqCov/20x \
-y fastq --normalize 20


Kmer normalization using Stacks will normalize read depth according to kmer coverage. First,
the program will count the number of occurrences of each kmer in a sequence read. Then 
it will look up the count for each of the kmers that occurred in the entire dataset, and 
finally, calculate a median kmer coverage for each read. A read is discared if the median 
for that read is higher than the threshhold. In this case, we will discard reads which have 
coverage greater than 40x and 20x. 


5. At this point you should have 4 sets of processed paired-end reads: 
Count and report the number of total read pairs for each set.

- a set that has been cleaned(trimmed)

for file in *.fastq; do readPairs=$(grep -A 1 "@2_" $file | grep -v "\_" | grep -v "\-" | wc -l);\
echo $file":" $readPairs; done
SscoPEclean_R1.fastq: 48,984,846
SscoPEclean_R2.fastq: 48,984,846

- a set that has be cleaned and k-mer filtered

for file in *.fq; do readPairs=$(grep -A 1 "@2_" $file | grep -v "\_" | grep -v "\-" | wc -l);\
echo $file":" $readPairs; done
SscoPE_R1.1.filRare.fq: 40,423,875
SscoPE_R2.2.filRare.fq: 40,423,875


- a set that has been cleaned, k-mer filtered,and normalized to 20X
for file in *.fq ; do readPairs=$( grep -A 1 "@2_" $file | grep -v "\_" | grep -v "\-" | wc -l); \
echo $file":" $readPairs ; done
SscoPE_R1.1.fil.norm20x.fq: 9,401,682
SscoPE_R2.2.fil.norm20x.fq: 9,401,682	


- a set that has been cleaned, k-mer filtered, and normalized to 40X

For file in *.fq ; do readPairs=$( grep -A 1 "@2_" $file | grep -v "\_" | grep -v "\-" | wc -l); \
echo $file":" $readPairs ; done
SscoPE_R1.1.fil.norm40x.fq: 14326197
SscoPE_R2.2.fil.norm40x.fq: 14,326,197


6. Now run kmer_filter 3 times,on the cleaned set,the cleaned/filtered set, and the
cleaned/filtered/20X normalized, set, using just the --k_dist flag, and nothing else.
Make sure to direct the output to a separate file each time. This will generate a 
distribution of k-mer frequencies for each set. Plot the distributions in R and include
the figures with the rest of the assignment. How do the 3 distributions differ?

PBS Job script for kmer distribution: 
 !/usr/bin/env bash
 PBS -N kmerFilter_Dist_trimmed
 PBS -q fatnodes
 PBS -l nodes=1:ppn=12 
 PBS -d /home12/jsmith16/bi623/160829_PS4/trimmedFastq/kmerDist
 PBS -o /home12/jsmith16/bi623/160829_PS4/trimmedFastq/kmerDist
 PBS -e /home12/jsmith16/bi623/160829_PS4/trimmedFastq/kmerDist

module load stacks

 cleaned fastq
kmer_filter -1 /home12/jsmith16/bi623/160829_PS4/trimmedFastq/SscoPEclean_R1.fastq \
-2  /home12/jsmith16/bi623/160829_PS4/trimmedFastq/SscoPEclean_R2.fastq  \
-o /home12/jsmith16/bi623/160829_PS4/trimmedFastq/kmerDist  --k_dist

 cleaned, filtered 
kmer_filter -1 /home12/jsmith16/bi623/160829_PS4/filterFastq/SscoPE_R1.1.filRare.fq  \
-2 /home12/jsmith16/bi623/160829_PS4/filterFastq/SscoPE_R2.2.filRare.fq  \
-o /home12/jsmith16/bi623/160829_PS4/filterFastq/kmerDist  --k_dist 

 cleaned, filtered, 20x normalized
kmer_filter -1 /home12/jsmith16/bi623/160829_PS4/filterFastqCov/20x/SscoPE_R1.1.fil.norm20x.fq \
-2  /home12/jsmith16/bi623/160829_PS4/filterFastqCov/20x/SscoPE_R2.2.fil.norm20x.fq \
-o /home12/jsmith16/bi623/160829_PS4/filterFastqCov/20x/kmerDist   --k_dist 


The distributions of the clean and clean filtered fastq files have similar distributions,
but the cleaned fastq does have a higher frequency of  abundant kmers, which are lost 
after trimming. The clean, filtered, and 20x normalized distribution is much more narrow,
and shows much fewer reads overall, as well as lower frequency of abundant kmers. 



 Block 2 

7. Copy these files into a working directory on the ACISS cluster.

Files used:
/research/bi610/data/reads/pipe.fil.1.fastq
/research/bi610/data/reads/pipe.fil.2.fastq
/research/bi610/data/reads/pipe.filnorm.1.fastq
/research/bi610/data/reads/pipe.filnorm.2.fastq



8. Using UNIX commands, count and report the number of reads in each dataset file.

 unix command to pull out header lines and count the headers to determine number of 
 read pairs in the files 

for file in *.fastq; do grep -c -H "^@2_" $file; done
pipe.fil.1.fastq:40,423,875
pipe.fil.2.fastq:40,423,875
pipe.filnorm.1.fastq:14,326,197
pipe.filnorm.2.fastq:14,326,197   

There are 40,423,875 read pairs in the unnormalized dataset and 14,326,197 read pairs
in the normalized dataset. 



9. Use Trinity to assemble both the normalized and unnormalized RNA-seq datasets.
Set --JM 50, --CPU 10, --min_contig_length 300, -- min_kmer_cov 3, and 
--group_pairs_distance 800.

PBS Job script for De Novo transcriptome assembly using Trinity with normalized reads: 
 !/usr/bin/env bash
 PBS -N Trinity
 PBS -q fatnodes
 PBS -l nodes=1:ppn=10 
 PBS -d /home12/jsmith16/bi623/160829_PS4

cd /home12/jsmith16/bi623/160829_PS4/block2

module load trinity/131110

 Normalized RNA-seq dataset
Trinity.pl --seqType fq --JM 50G --min_contig_length 300  --min_kmer_cov 3 --CPU 10 \
 --group_pairs_distance 800 --left pipe.filnorm.1.fastq  --right pipe.filnorm.2.fastq \
 --output /home12/jsmith16/bi623/160829_PS4/block2/trinityNorm 


PBS Job script for De Novo transcriptome assembly using Trinity with unnormalized reads: 
 !/usr/bin/env bash
 PBS -N Trinity
 PBS -q longfatnode
 PBS -l nodes=1:ppn=10 
 PBS -d /home12/jsmith16/bi623/160829_PS4

cd /home12/jsmith16/bi623/160829_PS4/block2

module load trinity/131110

 Unnormalized RNA-seq dataset
Trinity.pl --seqType fq --JM 50G --min_contig_length 300  --min_kmer_cov 3 --CPU 10 \
 --group_pairs_distance 800 --left pipe.fil.1.fastq --right pipe.fil.2.fastq    \
 --output /home12/jsmith16/bi623/160829_PS4/block2/trinityUnnorm



10. Using information in the final assembly files (named Trinity.fasta, but please 
rename accordingly), calculate the number of transcripts, the maximum transcript length, 
the minimum transcript length, the mean transcript length, and the median transcript
length. Plot the contig length distributions for each assembly using R hist() and 
boxplot() functions. Make sure your plots are clean and labeled appropriately.

Summary statistics for the normalized dataset:
>python3.5 ps4_ContigSummaryStats.py
	the total number of contigs in this dataset is 139042.
	the maximum length of contigs in this dataset is 22465.
	The minimun contig length is 301.
	The mean contig lenght is 2556.7385682024137.
	The median contig length is 2116.0.


Summary statistics for the normalized dataset:
>python3.5 ps4_ContigSummaryStats.py 
	the total number of contigs in this dataset is 113427.
	the maximum length of contigs in this dataset is 18857.
	The minimun contig length is 301.
	The mean contig lenght is 2478.1113491496735.
	The median contig length is 1963.



11. Based on your assembly statistics and what you know about transcripts, 
comment on whether there is a clear difference in quality between the two Trinity 
assemblies. Also comment on differences in the total number of transcripts and transcript 
groups (“loci”) for the two assemblies. 


 Unix Commands to find groups of transcripts in then normalized dataset. This will pull out
 all transcripts from the header of the fasta file, and then pull out transcripts 
 with counts of more than 1. 
cat pipeNorm_trinity.fasta | grep ">" | \
sed -E 's/>(comp[0-9]+\_[c][0-9]+)\_seq[0-9]+ len.*/\1/' | sort | uniq -c | sort -nk1 |\
grep -v -E "[ ]{1,}[1]{1} comp[0-9]+.*" | wc -l 
10777

 unix commands to find groups of transcripts in the unnormalized dataset. 
cat pipeUnnorm_trinity.fasta | grep ">" | \
sed -E 's/>(comp[0-9]+\_[c][0-9]+)\_seq[0-9]+ len.*/\1/' | sort | uniq -c | sort -nk1 | \
grep -v -E "[ ]{1,}[1]{1} comp[0-9]+.*" | wc -l 
10401
 
The normalized assembly had a better output due to the fact that it had a longer
maximum contig length, 22,465bp versus 18,857bp in the unnormalized assembly. Also the mean 
and median contig lengths are longer, at 2,556bp and 2,116bp  indicating more efficient 
assembly. In addition, the boxplots indicate that the normalized distribution of contig
lengths are less variable, since the inter-quartile range is smaller, and the upper
whisker is shorter. When looking at "groups" of transcripts, in which the transcript
name is identical, but have more than one sequence assembled, you can see that more 
"groups" were found in the normalized dataset, with 10,777 identified. While the unnormalized
had a total of 10,401 groups. These groups indicate the existence of transcript variants 
that the assembler was able to identify, and it appears the normalized dataset increased
resolution of these transcript variants. 



 Block 3

12.Based on the assembly profiles, choose your best Trinity assembly to use in 
this assignment.

Dataset used:
/home12/jsmith16/bi623/160829_PS4/block2/trinityNorm/pipeNorm_trinity.fasta

This is the normalized dataset produced in Block 2. 



13.Create a BLAST database from your assembled transcriptome.
 We don't need to create a blast database of the assembled transcriptome. 

14.BLAST your pipefish transcripts against threespine stickleback protein sequences.
a. Since your pipefish data are nucleotides, and the stickle back data are amino acids, 
you will need to use a translating BLAST program. Only consider hits with an e-value 
less than 1x10-5, and retain only the best stickleback hit for each transcript.

PBS Job script for Blastx with transcript assmebly query agains a Stickleback protein database:
 PBS -N blastx 
 PBS -q generic
 PBS -l nodes=1:ppn=12
 PBS -d /home12/jsmith16/ 
 PBS -e /home12/jsmith16/bi623/160829_PS4/block2
 PBS -o /home12/jsmith16/bi623/160829_PS4/block2
 PBS -M jsmith16@uoregon.edu
 PBS -m ae 


cd /home12/jsmith16/bi623/160829_PS4/block2/trinityNorm 

module load blast

blastx -query /home12/jsmith16/bi623/160829_PS4/block2/trinityNorm/pipeNorm_trinity.fasta \
-db /home12/jsmith16/bi623/160818_PS1/stickleback -outfmt 6 -out NormalizedTrinity_blast.txt \
-num_threads 12 -evalue 0.00001 -max_hsps_per_subject 1 -max_target_seqs 1




15.Download the threespine stickleback GTF file from Ensembl.

File used:
/research/bi610/data/rnaseq/Gasterosteus_aculeatus.BROADS1.85.gtf



16. For each top stickleback BLAST hit you identify, look up the human readable gene name
(AKA external gene ID) from the stickleback GTF file. 

Please see "ps4_blastHits_final.py". 


17.Produce a table of the top BLAST hits you identified for each pipefish transcript, 
including the transcript IDs, e-value, and, if available, external ID.

Please see attached "Final_blast_outfile.txt" 


18. Read about cegma

http://korflab.ucdavis.edu/datasets/cegma/README


19. Run CEGMA on the normailzed trinity assembly to see how many of the conserved 248 
genes are found in your transcriptome assembly. How many of the 248 CEGs are “complete” in 
the pipefish transcriptome?

PBS Job script for CEGMA:
 !/usr/bin/env bash
 PBS -N cegma
 PBS -q fatnodes
 PBS -l nodes=1:ppn=12 
 PBS -d /home12/jsmith16/bi623/160829_PS4/block2
 PBS -o /home12/jsmith16/bi623/160829_PS4/block2 
 PBS -e /home12/jsmith16/bi623/160829_PS4/block2 


module load cegma/2.4

cegma -g /home12/jsmith16/bi623/160829_PS4/block2/trinityNorm/pipeNorm_trinity.fasta \
-T 12 -o pipefish

The results of the completeness report indicate that this De Novo transcriptome assembly
is 97.98% complete, since 243 or 248 highly conserved proteins (core eukaryotic proteins,
also called CEGs) were identified. Partial CEGs were identified at 247 out of the 248 
reference CEGs. However, there may be a large amount of redundancy in the assembly 
because the average number of orthologs is much greater than 1. There is approximately
3 or more or orthologs on average in the de novo assembly. 




















