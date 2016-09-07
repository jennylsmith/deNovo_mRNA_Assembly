#PS4 De Novo Transcriptome Assembly
#August 30, 2016 
#script to plot the contig length distributions for each assembly using R hist() and boxplot() functions. 



#set working directory
getwd()
setwd("//Users/jsmith/Documents/bi623/160829_class11_PS4/block2_trinityAssembly/")

#read in the contig lengths table. 
norm <- read.table("contigLengths_Norm.txt", header = FALSE, sep = "")
head(norm)
dim(norm)
max(norm)

#read in the contig lengths table. 
unnorm <- read.table("contigLengths_Unnorm.txt", header = FALSE, sep = "")
head(unnorm)
class(unnorm)
dim(unnorm)
max(unnorm)

#create a pdf
pdf(file = "Contig_distributions.pdf")
#set parameters
op <- par(mfrow = c(2, 2))
#plot the distrubutions using histogram function.
hist(norm[,1], main = "Distribution of Contig Lengths\n for a Normalized Dataset", xlab = "Contig Length (bp)",
     col = "steelblue", xlim = c(0,22500))
hist(norm[,1], main = "Distribution of Contig Lengths\n for an Unnormalized Dataset", xlab = "Contig Length (bp)",
     col = "slateblue4", xlim = c(0,22500))

#Plot the distributions with a boxplot function.
boxplot(norm[,1], main = "Distribution of Contig Lengths\n for a Normalized Dataset",
        col = "steelblue", ylab = "Contig Length (bp)", cex = 0.25)
boxplot(unnorm[,1], main = "Distribution of Contig Lengths\n for an Unnormalized Dataset",
        col = "slateblue4", ylab = "Contig Length (bp)", cex = 0.25)
dev.off()




