#Phylogenetic Analysis
#Practical 02
#Date: 16/01/2023
#BUILDING NEIGHBOR JOINING TREES FROM MOLECULAR DATA

#install packages
install.packages("BiocManager")
BiocManager::install("msa")
library(msa)
install.packages("phangorn")
#BiocManager::install("phangorn")
library(phangorn)

#set working directory
setwd("D:/LMS 3rd Year - Sem 2 - Bioinformatics/BT 3167 - Phylo. Analysis/Practicals/Prac2")

# Read in a FASTA file
dna = readDNAStringSet( "PR2_sequences.fasta" )

# Perform multiple sequence alignment using ClustalW and MUSCLE
clustalw_aln = msa(dna,method= "ClustalW" )
muscle_aln = msa(dna, method = "Muscle" )

# Convert the MSA object to a DNAbin object for tree construction
#convert function also can be used
clustal_bin = as.DNAbin(clustalw_aln) 
muscle_bin = as.DNAbin(muscle_aln)

# Calculate the pairwise distance matrix - K80
#clustalw-K80
clustalw_dists_K80 = dist.dna(clustal_bin, model = "K80" )
summary(clustalw_dists_K80)
clustalw_dists_K80

#muscle-K80
muscle_dists_K80 = dist.dna(muscle_bin, model = "K80" )
summary(muscle_dists_K80)
muscle_dists_K80

# Calculate the pairwise distance matrix - JC69
#clustalw-JC69
clustalw_dists_JC69 = dist.dna(clustal_bin, model = "JC69" )
summary(clustalw_dists_JC69)
clustalw_dists_JC69

#muscle-JC69
muscle_dists_JC69 = dist.dna(muscle_bin, model = "JC69" )
summary(muscle_dists_JC69)
muscle_dists_JC69

# Generate the neighbor joining plots for the alignments
clustalw_tree_K80 = nj(clustalw_dists_K80)
muscle_tree_K80 = nj(muscle_dists_K80)
clustalw_tree_JC69 = nj(clustalw_dists_JC69)
muscle_tree_JC69 = nj(muscle_dists_JC69)

# Convert the alignment data into character-based form
clustalw_data = as.phyDat(clustalw_aln)
muscle_data = as.phyDat(muscle_aln)

# Display the parsimony scores for the generated trees
parsimony(clustalw_tree_K80,clustalw_data)
parsimony(muscle_tree_K80,muscle_data)
parsimony(clustalw_tree_JC69,clustalw_data)
parsimony(muscle_tree_JC69,muscle_data)

# Parsimony optimization for the constructed trees
clustalw_opt_K80 = optim.parsimony(clustalw_tree_K80,clustalw_data)
muscle_opt_K80 = optim.parsimony(muscle_tree_K80,muscle_data)
clustalw_opt_JC69 = optim.parsimony(clustalw_tree_JC69,clustalw_data)
muscle_opt_JC69 = optim.parsimony(muscle_tree_JC69,muscle_data)

# Display the parsimony scores for the optimize trees
parsimony(clustalw_opt_K80,clustalw_data)
parsimony(muscle_opt_K80,muscle_data)
parsimony(clustalw_opt_JC69,clustalw_data)
parsimony(muscle_opt_JC69,muscle_data)

# plot the optimized NJ trees
par(cex=0.85, mar = c(0.1,0.1,1.5,0.1))
plot(clustalw_opt_K80, type = "unrooted" )
title("ClustalW optimized Neighbor joining tree-K80" )

plot(muscle_opt_K80, type = "unrooted" )
title("MUSCLE optimized Neighbor joining tree-K80" )

plot(clustalw_opt_JC69, type = "unrooted" )
title("ClustalW optimized Neighbor joining tree-JC69" )

plot(muscle_opt_JC69, type = "unrooted" )
title("MUSCLE optimized Neighbor joining tree-JC69" )

# plot the unoptimized NJ trees
par(cex=0.85, mar = c(0.1,0.1,1.5,0.1))
plot(clustalw_tree_K80, type = "unrooted" )
title("ClustalW unoptimized Neighbor joining tree-K80" )

plot(muscle_tree_K80, type = "unrooted" )
title("MUSCLE unoptimized Neighbor joining tree-K80" )

plot(clustalw_tree_JC69, type = "unrooted" )
title("ClustalW unoptimized Neighbor joining tree-JC69" )

plot(muscle_tree_JC69, type = "unrooted" )
title("MUSCLE unoptimized Neighbor joining tree-JC69" )
