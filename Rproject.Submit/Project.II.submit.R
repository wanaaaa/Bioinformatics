#########################################################
#Q5
library("bio3d"); library("msa"); library(seqinr);library(Biostrings)

load("./Data/ts1r3_humanStr.data")
load("./Data/ts1r3_mouseStr.data")

resProjet <- hmmer(ts1r3_humanStr, type="phmmer", db="swissprot")
save(resProjet, file = "./Data/resProject.data")
load("./Data/resProject.data")
resProjet

#######################################################################
#Q6
famProject <- hmmer(ts1r3_humanStr, type="hmmscan", db="pfam")

#save(famProject, file =  "./Data/famProject.data")
#load("./Data/famProject.data")
famProject$hit.tbl  # This has all the domain hits

####################################################################################
#Q7
library(msa)
flySequences <- readAAStringSet("./fly_protein_unaligned.fasta", format = "fasta")
flyAlignment <- msa(flySequences)
flyAlignment
print(flyAlignment, show="complete")

seq219_244 <- "FYCDLCEYQTNTRYSFLRHKKSKKHQ"
fam219_244 <- hmmer(seq219_244, type="hmmscan", db="pfam")
fam219_244$hit.tbl
pfam219_244.aln <- pfam(fam219_244$hit.tbl$acc[1])
pfam219_244.aln

hmm219_244 <- hmmer(pfam219_244.aln, type="hmmsearch", db="pdb")
hmm219_244$hit.tbl
hmm219_244$hit.tbl$species
hmm219_244$hit.tbl[3, ]
####################################################################################
#Q8

library("HMM")
library(seqinr)

States = c("Exon", "5site",  "SecBaseIntron" ,"Intron")

Symbols = c("A","C","G","T")

transProbs = matrix(c('EE'=0.9,'E5'=0.1, 'ES'=0 ,'EI'=0, 
                      '5E'=0, '55'=0, '5S'= 1 ,'5I'=0,
                      'SE'=0, 'S5'=0, 'SS'=0, 'SI'=1.0,
                      'IE'=0, 'I5'=0, 'IS'= 0 ,'II'=1.0), 
                    c(length(States), length(States)), byrow = TRUE)

emissionProbs = matrix(c('A'=0.25,'C'=0.25,'G'=0.25,'T'=0.25, 
                         'A'=0.05,'C'=0.0,'G'=0.95,'T'=0.0,
                         'A'=0.05, 'C'=0.05, 'G'= 0.1, 'T'=0.8,
                         'A'=0.4,'C'=0.1,'G'=0.1,'T'=0.4), 
                          c(length(States), length(Symbols)), byrow = TRUE)

hg38_seq <- s2c('ACGAGGCGTTCATCGAGGAGGGCACATTCCTTTTCACCTCAGAGTCGGTCGGGGAAGGCCACCCAGGTGAGGGGACGGCCTGAAGCGAAGCGTGGGGCGGGGCAGAAGGCAGCGCCAAGGTCCGGCTGGCTGCGGCCGGCCGGTGGTGGGGCCCGC')
observed_seq <- s2c('CTTCATGTGAAAGCAGACGTAAGTCA')

hmm <- initHMM(States, Symbols, startProbs = c(1,0,0,0), transProbs = transProbs, 
               emissionProbs = emissionProbs)

vit <- viterbi(hmm, hg38_seq)
vitObserved <- viterbi(hmm, observed_seq)

sim <- simHMM(hmm, length(hg38_seq)) 
#sim <- simHMM(hmm, length(observed_seq)) 

post <- posterior(hmm, hg38_seq)
#post <- posterior(hmm, observed_seq)
############################################################################################
#Q9
library("bio3d"); library("msa"); library(seqinr); library(Biostrings); 
library("ape"); library("seqinr"); library(phangorn)

alnANF_receptor <- pfam("ANF_receptor", alignment = "seed")
#alnANF_receptor <- pfam("PF01094.27", alignment = "seed")
#save(alnANF_receptor, file = "./Data/alnANF_receptor.data")
load("./Data/alnANF_receptor.data")

speciesVec <- alnANF_receptor$id
speciesVec

iA = 1;
fastaStr = "";
cat(NULL, file = "domains.fasta")
for (specie in speciesVec) {
  fastaStr <- paste(">", specie, sep = "")
  
  write(fastaStr, file = "domains.fasta", append = TRUE)
  seqA <- alnANF_receptor$ali[iA, ]
  seqA <- paste(seqA, collapse = "")
  
  write(seqA, file = "domains.fasta", append = TRUE)
  iA <- iA + 1

}

alnForTree <- read.phyDat("domains.fasta", format = "fasta", type = "AA")
dmDomain <- dist.ml(alnForTree)
treeUPGMAall <- upgma(dmDomain)
treeNJall <- NJ(dmDomain)

### Distance based methods
plot(treeUPGMAall, main="UPGMA")
plot(treeNJall, "unrooted", main="NJ")

### Parsimony
#get parsimony score

parsimony(treeUPGMAall, alnForTree, method = "sankoff")
parsimony(treeNJall, alnForTree,method = "sankoff")

treePars <- optim.parsimony(treeUPGMAall, alnForTree, method = "sankoff")
treePars <- optim.parsimony(treeNJall, alnForTree, method = "sankoff")

plot(treePars, "unrooted", main="parsimony")

### Maximum likelihood 
#fit = pml(treeNJall, data=alnForTree, model="Jukes-Cantor")
fit = pml(treeNJall, alnForTree)
fit
methods(class="pml")
logLik(fit)

fitJC <- optim.pml(fit, TRUE)
logLik(fitJC)

bs <- bootstrap.pml(fitJC, bs=100, optNni=TRUE, multicore=FALSE, control = pml.control(trace=0))

save(bs, file = "./Data/bs.data")
load("./Data/bs.data")

plotBS(midpoint(fitJC$tree), bs, p = 50, type="p")
write.tree(bs, file="bsDomains.tre")
############################################################
#Q10
toyTree <- read.tree("toy_tree.tre")
plot(toyTree, type="unrooted", main="Phylogenetic Tree of toy_tree")

toySeq <- readAAStringSet("toy_MSA.fasta")
toySeq

#L1 and L2
atVec <- c(0.09, 0.103, 0.12)
for(at in atVec) {
  E <- (1/4)*(1+3*exp(-4*at))
  D <- (1/4)*(1-exp(-4*at))
  L1 <- (E^4*D + 2*E^3*D^2 + 4*E^2*D^3 + 2*E*D^4 + 7*D^5)/4
  L2 <- (2*E^3*D^2 + 2*E^2*D^3 + 10*E*D^4 + 2*D^5)/4
  lnL1 <- log(L1)
  lnL2 <- log(L2)
  resultStr <- paste("at=", at, " ", "E=", E, " " ,"D=",D, " ", "L1=", 
         L1, " ", "L2=", L2, " ", "lnL1=", lnL1, " ", "lnL2=", lnL2, " ", sep = " ")
  print(resultStr)
}
#-------
toyDNA <- read.phyDat("toy_MSA.fasta", format = "fasta", type = "DNA")
toyDNA

fit = pml(toyTree, data = toyDNA, model = "Jukes-Cantor")
fit
methods(class="pml")
logLik(fit)

#fitJC <- optim.pml(fit, TRUE)
#logLik(fitJC)


