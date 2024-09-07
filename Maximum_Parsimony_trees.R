# Phylogenetic Analysis
# Practical 03
#BUILDING MAXIMUM PARSIMONY TREES AND CARRYING OUT BOOTSTRAP 
ANALYSIS IN R 

#Set working directory
setwd("D:/LMS 3rd Year - Sem 2 - Bioinformatics/BT 3167 - Phylo. Analysis/Practicals/Prac3/R")

# Load the packages into the R session
library(msa)
library(phangorn)

# Read the sequence file
fasta = readAAStringSet("PR3_sequences.fasta", format = "fasta")

# Perform the multiple sequence alignment
aln = msa(fasta, method = "ClustalW")

# Converting the alignment to other formats for use later
phyd = as.phyDat(aln)
aabin = as.AAbin(aln)

# calculating distance matrix and building NJ tree
dist_mat = dist.ml(aabin, model = "JTT")
tree = nj(dist_mat)
plot(tree)

# Building the MP tree using the parsimony ratchet for 100 iterations minimum
ptree = pratchet(phyd,minit = 100, trace = 0)

#ploting the output tree
plotBS(ptree, bs.col = "red", type = "phylogram")

# Get the CI and RI values for the tree - ptree
CI(ptree,phyd)
RI(ptree,phyd)

# Carrying out bootstrap analysis using a given tree and an alignment
# Here B is the number of iterations (1000) to bootstrap and FUN is the function
# that was used to generate the tree from the original data
bstrees <- boot.phylo(tree, aabin, FUN = 
                        function(e) nj(dist.ml(e,model = "JTT")),
                      B = 1000, trees = T)$trees

btree = tree

# Get the CI and RI values for the tree - bstraptree
CI(btree,phyd)
RI(btree,phyd)

# Plot the MP tree with bootstrap values
par(cex=0.7, mar = c(0.1,.01,0.1,0.1))
plotBS(tree, bstrees, bs.col = "blue ", frame = "rect", type = "phylogram" )

# The exhaustive search tree
ex_tree = bab(phyd, trace = 0)
plot(ex_tree)
