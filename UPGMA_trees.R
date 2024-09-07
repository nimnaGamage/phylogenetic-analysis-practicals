#Phylogenetic Analysis
#Practical 01
#BUILDING UPGMA TREES USING R

#set working directory
setwd("D:/LMS 3rd Year - Sem 2 - Bioinformatics/BT 3167 - Phylo. Analysis/Practicals/Prac1/R")

#Install the necessary packages
install.packages("ape")
install.packages("BiocManager")
BiocManager::install("msa")

#load the packages into the R session
library(msa)
library(ape)

#read in a fasta file
dna = readDNAStringSet("sequences.fasta")

#perform multiple sequence alignment using clustalW AND MUSCLE
clustalw_aln = msa(dna,method = "ClustalW")
muscle_aln = msa(dna, method = "Muscle")

#Convert the MSA object to a DNAbin object for tree construction
clustal_bin = msaConvert(clustalw_aln,type = "ape::DNAbin")
muscle_bin = msaConvert(muscle_aln,type = "ape::DNAbin")

#calculate the pairwise distance matrix using the kimura 2-parameter model
clustalw_dists = dist.dna(clustal_bin)
clustalw_dists_GG95 = dist.dna(clustal_bin, model = "GG95")
clustalw_dists_indel = dist.dna(clustal_bin, model = "indel")
muscle_dists = dist.dna(muscle_bin, model = "K80")
muscle_dists_JC69 = dist.dna(muscle_bin, model = "JC69")
muscle_dists_F81 = dist.dna(muscle_bin, model = "F81")

summary(clustalw_dists)
clustalw_dists
summary(clustalw_dists_GG95)
clustalw_dists_GG95
summary(clustalw_dists_indel)
clustalw_dists_indel
summary(muscle_dists)
muscle_dists
summary(muscle_dists_JC69)
muscle_dists_JC69
summary(muscle_dists_F81)
muscle_dists_F81

#generate the upgma plot for the alignments
clustalw_upgma = hclust(clustalw_dists, method = "average")
muscle_upgma = hclust(muscle_dists, method = "average")

#extra step to make the upgma trees display properly when plotting
clustalw_upgma = as.dendrogram(clustalw_upgma)
muscle_upgma = as.dendrogram(muscle_upgma)

#plotting the trees
png(file="clustalw.png",1200,900)
par(cex=1.2,mar=c(5,1,2,4))
plot(clustalw_upgma, horiz = TRUE)
title("Clustalw UPGMA Tree")
dev.off()
png(file="MUSCLE.png",1200,900)
par(cex=1.2,mar=c(5,1,2,4))
plot(muscle_upgma, horiz = TRUE)
title("MUSCLE UPGMA Tree")
dev.off()
