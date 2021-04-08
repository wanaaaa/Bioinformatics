##################################
## Question 1
library(seqinr)
library(Biostrings);
ts1r3_human <- read.fasta(file = "TS1R3_HUMAN.fastA")
ts1r3_mouse <- read.fasta(file = "TS1R3_MOUSE.fastA")

ts1r3_humanSeq <- ts1r3_human[[1]]
ts1r3_mouseSeq <- ts1r3_mouse[[1]]

humanSeqForDot <- ts1r3_humanSeq
mouseSeqForDot <- ts1r3_mouseSeq

ts1r3_humanStr <- c2s(ts1r3_humanSeq)
ts1r3_mouseStr <- c2s(ts1r3_mouseSeq)

ts1r3_humanStr <- toupper(ts1r3_humanStr)
ts1r3_mouseStr <- toupper(ts1r3_mouseStr)

save(ts1r3_humanStr, file = "./Data/ts1r3_humanStr.data")
save(ts1r3_mouseStr, file = "./Data/ts1r3_mouseStr.data")

data(BLOSUM50)

# Pairwise local alignment of protein sequences using the Smith-Waterman algorithm
localAlignSmithWaterman <- pairwiseAlignment(ts1r3_humanStr, ts1r3_mouseStr, substitutionMatrix = BLOSUM50, gapOpening = -2, gapExtension = -8, scoreOnly = FALSE, type="local")

localAlignSmithWaterman # Print out the optimal local alignment and its score

##---
source("seq_print_utils.R")
#printPairwiseAlignment(localAlignSmithWaterman, 60)
###
source("calculateMatchMismatch.R")

calMatchMisMatch(localAlignSmithWaterman, 60)
##---
# We will use a simulation to calculate a p-value.
# Generate random (null hypothesis) sequences
randomseqs <- generateSeqsWithMultinomialModel(ts1r3_mouseStr,1000)
randomseqs[1:10] # Print out the first 10 random sequences

# Align ts1r3_humanStr to the first of the 1000 random sequences:

# get corresponding score
pairwiseAlignment(ts1r3_humanStr, randomseqs[1], substitutionMatrix = "BLOSUM50", gapOpening = -2,gapExtension = -8, scoreOnly = TRUE)

# get score against non-random sequence
pairwiseAlignment(ts1r3_humanStr, ts1r3_humanStr, substitutionMatrix = "BLOSUM50", gapOpening = -2,gapExtension = -8, scoreOnly = TRUE)

# repeat this 1000 times to get a null distribution
randomscores <- double(1000) 
# Create a numeric vector with 1000 elements
for (i in 1:1000) {
  score <- pairwiseAlignment(ts1r3_humanStr, randomseqs[i], substitutionMatrix = "BLOSUM50",gapOpening = -2, gapExtension = -8, scoreOnly = TRUE)
  randomscores[i] <- score
}

randomscores

hist(randomscores, col="red") # plot a histogram


sum(randomscores >= 4180)/1000
#######################################################
##Question 2
dotPlot(humanSeqForDot, mouseSeqForDot)






