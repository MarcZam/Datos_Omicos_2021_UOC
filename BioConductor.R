# Part 1 Introduction

setwd("D:/Docs/Masteres/Master en Bioinformática y Bioestadística/Tercer Semestre/M0.157 - Análisis de datos ómicos/Moudlo02/Debate")
# set working directory

info <-readLines("GSE58435_series_matrix.txt", n=70)
# reads the first 70 lines in the file where the information is located

rows2read <- 54743 -66 -2
# sets the number of rows to be read to be 54743 minus the first 66 and the last 2
x <-read.table("GSE58435_series_matrix.txt", skip=66, header=TRUE, sep="\t",row.names=1, nrows = rows2read)
# reads the table of microarray results

dim(x)
# we can see that now the dimensions of our array are 54675 observations upon 10 variables of which the first 5
# are of Turner syndrome and the remaining 5 are control variables

colnames(x)<- c(paste("Turner", 1:5, sep="_"), paste("Control",1:5, sep="_"))
colnames(x)
# we change the names of our columns to match the group of each of the samples

head(x)

round(apply(x,2, summary),3)
# we show some basic statistics of our data rounded to 3 decimal points

boxplot(x, main = "Expression values for 5 Turner and 5 Control samples", xlab = "Slides", ylab = "Expression", col = rep(c("red", "green"),each=5))
# we show a boxplot to verify symmetry of results

pcX <-prcomp(t(x), scale=TRUE) 
# prcomp - performs a principal components analysis
# we use the t - transpose function to turn over the matrix of results
loads<- round(pcX$sdev^2/sum(pcX$sdev^2)*100,1)
# obtains the relative standard deviation of each sample in reference to the rest of the samples in the batch

xlab<-c(paste("PC1",loads[1],"%"))
ylab<-c(paste("PC2",loads[2],"%"))
plot(pcX$x[,1:2],xlab=xlab,ylab=ylab, xlim=c(-150, 150))
title("Principal components (PCA)")
text(pcX$x[,1],pcX$x[,2],colnames(x), pos=4)

clust.euclid.average <- hclust(dist(t(x)),method="ward.D2")
plot(clust.euclid.average, hang=-1)

# Part 2 Bioconductor

require(Biobase)

expressionValues <- matrix(rnorm (300), nrow=30)
colnames(expressionValues) <- paste0("sample",1:10)
head(expressionValues)

targets <- data.frame(sampleNames = paste0("sample",1:10),
                      group=c(paste0("CTL",1:5),paste0("TR",1:5)),
                      age = rpois(10, 30), 
                      sex=as.factor(sample(c("Male", "Female"),10,replace=TRUE)),
                      row.names=1)
head(targets, n=10)

myGenes <-  paste0("gene",1:30)

pcs <- prcomp(expressionValues)
names(pcs)
barplot(pcs$sdev)
