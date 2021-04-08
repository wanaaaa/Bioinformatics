library("bio3d"); library("msa"); library(seqinr); library(Biostrings); 
library("ape"); library("seqinr"); library(phangorn)

#Q11

alnForTree <- read.phyDat("sars.fasta", format = "fasta", type = "DNA")
dmDomain <- dist.ml(alnForTree)
treeNJ <- NJ(dmDomain)
virusTree <- root(treeNJ, resolve.root = T, outgroup = "Himalayan palm civet sars cov, complete genome")

plot(virusTree)
#################################################################
#Q12
#i)
library(ROCR); library(class)

allData<-load("ALL_example.RData")
load("ALL_example.RData")

predProject <- knn(train = train_expr, test = test_expr, cl = trainLabels, k = 1)

errRate <- mean(testLabels != predProject)
errRate

#ii)
combined_expr <- rbind(train_expr, test_expr)

trainLabelCharc <- as.character(trainLabels); testLabelCharc <- as.character(testLabels)
combinedLabelFactor <- c(trainLabelCharc, testLabelCharc)
combinedLabelFactor <- factor(combinedLabelFactor)

predCombined <- knn(train = combined_expr, test = combined_expr, cl = combinedLabelFactor, k = 1)

errRateCombined <- mean( combinedLabelFactor != predCombined)
errRateCombined
#####################################################################
#Q13
library("HMM")

States = c("Exon", "base1", "base2", "base3", "base4", "Intron")
Symbols = c("A","C","G","T")
transProbs = matrix(c('EE'=0.9,'E1'=0.1, 'E2'=0 , 'E3'=0, 'E4'=0, 'EI'=0,  
                      '1E'=0,  '11'=0,   '12'= 1, '13'=0, '14'=0, '1I'=0,
                      '2E'=0,  '21'=0,   '22'=0,  '23'=1, '24'=0, '2I'=0, 
                      '3E'=0,  '31'=0,   '32'=0,  '33'=0, '34'=1, '3I'=0,
                      '4E'=0,  '41'=0,   '42'=0,  '43'=0, '44'=0, '4I'=1,
                      'IE'=0,  'I1'=0,   'I2'=0,  'I3'=0, 'I4'=0, 'II'=1), 
                    c(length(States), length(States)), byrow = TRUE)

emissionProbs = matrix(c('A'=0.200, 'C'=0.300, 'G'=0.300, 'T'=0.200, 
                         'A'=0.997, 'C'=0.001, 'G'=0.001, 'T'=0.001,
                         'A'=0.001, 'C'=0.001, 'G'=0.997, 'T'=0.001,
                         'A'=0.001, 'C'=0.001, 'G'=0.997, 'T'=0.001,
                         'A'=0.001, 'C'=0.001, 'G'=0.001, 'T'=0.997,
                         'A'=0.150, 'C'=0.350, 'G'=0.350, 'T'=0.150), 
                       c(length(States), length(Symbols)), byrow = TRUE)

hg38_seq <- s2c('ACGAGGCGTTCATCGAGGAGGGCACATTCCTTTTCACCTCAGAGTCGGTCGGGGAAGGCCACCCAGGTGAGGGGACGGCCTGAAGCGAAGCGTGGGGCGGGGCAGAAGGCAGCGCCAAGGTCCGGCTGGCTGCGGCCGGCCGGTGGTGGGGCCCGC')

hmm <- initHMM(States, Symbols, startProbs = c(1,0,0,0), transProbs = transProbs, 
               emissionProbs = emissionProbs)

vit <- viterbi(hmm, hg38_seq)

sim <- simHMM(hmm, length(hg38_seq)) 

post <- posterior(hmm, hg38_seq)
