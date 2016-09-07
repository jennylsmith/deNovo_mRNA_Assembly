#PS4 De Novo Transcriptome Assembly
#August 30, 2016 
#script to plot length distribution of sequence reads after trimming
#script to plot kmer distribution of trimmed, rare filtered, and normalized reads


#set working directory 
getwd()
setwd("/Users/jsmith/Documents/bi623/160829_class11_PS4/kmerDist/")

#read in the text files with the read length distributions
R1 <- read.table("R1_lengthDist.txt", header = FALSE, sep = "")
R2 <- read.table("R2_lengthDist.txt", header = FALSE, sep = "")

#define a function to add a header to the text files. 
add_header <- function(x){
  names(x) <- c("Frequency", "Length of Read")
  x
}

#Log2 transfrom the frequencies in R1 to help with scaling when plotted. 
logR1 <- R1
logR1[, 1] <- log2(R1[, 1])

#Add a header to the dataset 
logR1 <- add_header(logR1)
head(logR1)

#Log2 transform the frequencies in R2. 
logR2 <-R2
logR2[,1] <- log2(R2[,1])

#add a header to the dataset
logR2 <- add_header(logR2)
head(logR2)

#create a pdf file for the frequency plots.
pdf(file = "ReadLengthDistributions.pdf")
#set parameters to have both plots in one figure 
op <- par(mfrow = c(2, 1))
#plot the read length frequencies. Read length on the x-axis and frequency on the y-axis.
plot(logR1[,2], logR1[,1], xlim = c(20,101), ylab = "Frequency (log2)", xlab = "Read Length (bp)", 
     col = "darkmagenta", pch = 19, type = "h", main = "Distribution of Read Lengths \n after Trimming R1") 
plot(logR2[,2], logR1[,1], xlim = c(20,101), ylab = "Frequency (log2)", xlab = "Read Length (bp)", 
     col = "darkred", pch = 19, type = "h", main = "Distribution of Read Lengths \n after Trimming R2")
#close the device
dev.off()

#Plot the distribution of Kmer Frequencies after Trimming, Kmer Filtering, and Normalization

#set working directory
setwd("/Users/jsmith/Documents/bi623/160829_class11_PS4/kmerDist/")

#read in the kmer distribution text files as a data frame. 
clean <- read.table("kmerDist_trimmed.txt", header = TRUE, sep = "")
filtered <- read.table("kmerDist_CleanFiltered.txt", header = TRUE, sep = "")
normalized <- read.table("kmerDist_CleanFil20x.txt", header = TRUE, sep = "")


#create a pdf file for the frequency plots.
pdf( file = "/Users/jsmith/Documents/bi623/160829_class11_PS4/kmerDist/KmerFrequencyDistributions.pdf")

#set parameters to have both plots in one figure 
op <- par(mfrow = c(3, 1))

#use plot function to plot kmer frequency against kmer count for each of the pre-processed fastq datasets. 
#type = "h" is for histogram like vertical lines
#log = "y" log transforms the y axis for better scaling
plot(clean[,1], clean[,2], log = "y", col = "purple", main = "Kmer Frequency Distribution \n for Trimmed Reads",
     type = "h", xlab = "Kmer Frequency", xlim = c(0,100000), ylab = "Number of Occurances (log10)")
plot(filtered[,1],filtered[,2], log = "y", col = "blue", main = "Kmer Frequency Distribution \n for Trimmed and Filtered Reads",
     xlab = "Kmer Frequency", xlim = c(0,100000), ylab = "Number of Occurances (log10)", type = "h")
plot(normalized[,1], normalized[,2], log = "y", col = "red", 
     main = "Kmer Frequency Distribution \n for Trimmed, Filtered, and Normalized Reads",
     type = "h", xlab = "Kmer Frequency", xlim = c(0,100000), ylab = "Number of Occurances (log10)")

#turn off the device
dev.off()




