printPairwiseAlignment <- function(alignment, chunksize=60, returnlist=FALSE) {
   require(Biostrings)           # This function requires the Biostrings package
   seq1aln <- pattern(alignment) # Get the alignment for the first sequence
   seq2aln <- subject(alignment) # Get the alignment for the second sequence
   alnlen  <- nchar(seq1aln)     # Find the number of columns in the alignment
   starts  <- seq(1, alnlen, by=chunksize)
   n       <- length(starts)
   seq1alnresidues <- 0
   seq2alnresidues <- 0
   for (i in 1:n) {
      chunkseq1aln <- substring(seq1aln, starts[i], starts[i]+chunksize-1)
      chunkseq2aln <- substring(seq2aln, starts[i], starts[i]+chunksize-1)
      #print(chunkseq1aln[i])
      # Find out how many gaps there are in chunkseq1aln:
      gaps1 <- countPattern("-",chunkseq1aln) # countPattern() is from Biostrings package
      # Find out how many gaps there are in chunkseq2aln:
      gaps2 <- countPattern("-",chunkseq2aln) # countPattern() is from Biostrings package
      # Calculate how many residues of the first sequence we have printed so far in the alignment:
      seq1alnresidues <- seq1alnresidues + chunksize - gaps1
      # Calculate how many residues of the second sequence we have printed so far in the alignment:
      seq2alnresidues <- seq2alnresidues + chunksize - gaps2
      if (returnlist == 'FALSE') {
         print(paste(chunkseq1aln,seq1alnresidues))
         print(paste(chunkseq2aln,seq2alnresidues))     
         #print("-----------------")
         print(paste(' '))
         oneSeq = paste(chunkseq1aln,seq1alnresidues)
         twoSeq = paste(chunkseq2aln,seq2alnresidues)
         for( j in 1:chunksize) {
            oneStr = substr(oneSeq, j, j)
            twoStr = substr(twoSeq, j, j)
            threeStr = paste(oneStr, " ", twoStr)
            if(oneStr == twoStr) {
              print(threeStr)
            }
            
         }

      }
   }
   if (returnlist == 'TRUE'){
    vector1 <- s2c(substring(seq1aln, 1, nchar(seq1aln)))
    vector2 <- s2c(substring(seq2aln, 1, nchar(seq2aln)))
    mylist <- list(vector1, vector2)
    return(mylist)
   }
}


generateSeqsWithMultinomialModel <- function(inputsequence, X)
{
  # Change the input sequence into a vector of letters
  require("seqinr") # This function requires the SeqinR package.
  inputsequencevector <- s2c(inputsequence)
  # Find the frequencies of the letters in the input sequence "inputsequencevector":
  mylength <- length(inputsequencevector)
  mytable <- table(inputsequencevector)
  # Find the names of the letters in the sequence
  letters <- rownames(mytable)
  numletters <- length(letters)
  probabilities <- numeric() # Make a vector to store the probabilities of letters
  for (i in 1:numletters)
  {
    letter <- letters[i]
    count <- mytable[[i]]
    probabilities[i] <- count/mylength
  }
  # Make X random sequences using the multinomial model with probabilities "probabilities"
  seqs <- numeric(X)
  for (j in 1:X)
  {
    seq <- sample(letters, mylength, rep=TRUE, prob=probabilities) # Sample with replacement
    seq <- c2s(seq)
    seqs[j] <- seq
  }
  # Return the vector of random sequences
  return(seqs)
}