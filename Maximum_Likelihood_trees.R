#Phylogenetic Analysis
#Practical 04
#BUILDING MAXIMUM LIKELIHOOD TREES FROM MOLECULAR DATA

#Set working directory
setwd("D:/LMS 3rd Year - Sem 2 - Bioinformatics/BT 3167 - Phylo. Analysis/Practicals/Prac4/R")

# Load the packages into the R session
library(msa)
library(phangorn)

# Read the sequence file
fasta = readDNAStringSet("PR4_sequences_with_outgroup.fasta", format = "fasta")

# Perform the multiple sequence alignment
aln = msa(fasta, method = "ClustalW")

# Converting the alignment to other formats for use later
phyd = as.phyDat(aln)
dnabin = as.DNAbin(aln)

# calculating distance matrix and building NJ tree
dist_mat = dist.dna(dnabin, model = "RAW")
tree = nj(dist_mat)
plot(tree)

# Compute the likelihood of a phylogenetic tree given an alignment and an initial tree
mltree = pml(tree,phyd, k=4)
AIC(tree, k=log(nobs(object)))
AIC(mltree)

# Optimize the ml tree
opt_pml = optim.pml(mltree, optNni = TRUE, optBf = TRUE, optQ = TRUE, optGamma =
                      TRUE )
AIC(opt_pml)

# Root the tree using outgroup rooting
rooted_tree = root(opt_pml$tree, "Marchantia polymorpha" )
plot(rooted_tree)
plot(midpoint(opt_pml$tree))

# A function to make the plotted tree look neater
tr = ladderize(rooted_tree)
plot(tr)

# Bootstrapping ml trees
ts = bootstrap.pml(opt_pml,bs = 100, trees = T , multicore = T )
plotBS(rooted_tree,ts,type = "phylogram" )
