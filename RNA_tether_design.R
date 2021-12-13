
# Initalize required libraries
library(Biostrings)
library(tidyverse)

# Make empty character vector
nucs <- as.character()

# Determine all unique trinucleotide combinations
for(a in seq(1,4)){
  
  x <- DNA_BASES[a]
  
  for (b in seq(1, 4)) {
    
    y <- paste0(x, DNA_BASES[b])
      
      y <- paste0(y, DNA_BASES) 
      nucs <- append(nucs, y)
  }
  
}

nucs <- paste0(nucs, nucs)
nucs.all <- nucs

# Remove redundancy
remove.overlap <- function(nucs.position){

x <- nucs[nucs.position]
y1 <- narrow(x, start = 1, end = 3)
y2 <- narrow(x, start = 2, end = 4)
y3 <- narrow(x, start = 3, end = 5)

z <- nucs[-nucs.position]
z.1 <- grepl(y1, z)
z.1 <- z.1 == FALSE
z <- z[z.1]

z.2 <- grepl(y2, z)
z.2 <- z.2 == FALSE
z <- z[z.2]

z.3 <- grepl(y3, z)
z.3 <- z.3 == FALSE
z <- z[z.3]
z <- append(x, z)
return(z)}

nucs <- nucs.all

nucs <- remove.overlap(nucs.position = 1)
nucs <- remove.overlap(nucs.position = 2)
nucs <- remove.overlap(nucs.position = 3)
nucs <- remove.overlap(nucs.position = 4)
nucs <- remove.overlap(nucs.position = 5)
nucs <- remove.overlap(nucs.position = 6)
nucs <- remove.overlap(nucs.position = 7)
nucs <- remove.overlap(nucs.position = 8)
nucs <- remove.overlap(nucs.position = 9)
nucs <- remove.overlap(nucs.position = 10)
nucs <- remove.overlap(nucs.position = 11)
nucs <- remove.overlap(nucs.position = 12)
nucs <- remove.overlap(nucs.position = 13)
nucs <- remove.overlap(nucs.position = 14)
nucs <- remove.overlap(nucs.position = 15)
nucs <- remove.overlap(nucs.position = 16)
nucs <- remove.overlap(nucs.position = 17)
nucs <- remove.overlap(nucs.position = 18)
nucs <- remove.overlap(nucs.position = 19)
nucs <- remove.overlap(nucs.position = 20)
nucs <- remove.overlap(nucs.position = 21)
nucs <- remove.overlap(nucs.position = 22)
nucs <- remove.overlap(nucs.position = 23)
nucs <- remove.overlap(nucs.position = 24)

## Double check that all possible trinucleotide combinations are possible
v <- 0
for(i in seq(1, length(nucs.all))){
  
  if(length(grep(narrow(nucs.all[i], start = 1, end = 3), nucs)) >= 1) {print(i) ; v <- v+1}
  
}


nucs <- sort(nucs)
nucs

# Make list of RNA sequences, without a stem loop
rna.1 <- as.character(RNAStringSet(DNAStringSet(nucs)))
rna.x <- c()
for (i in seq(1, length(rna.1))){
  x <- strsplit(rna.1[[i]], "")
  rna.x <- append(rna.x, paste0("5'-FAM-AU", paste(x[[1]], collapse = ""), "AU-IABkFQ-3'"))
}
rna.1 <- rna.x

# Make list of RNA sequences, with a 4bp stemloop
rna.2 <- as.character(RNAStringSet(DNAStringSet(nucs)))
rna.x <- c()
for (i in seq(1, length(rna.2))){
  x <- strsplit(rna.2[[i]], "")
  y.1 <- x[[1]] == "A"
  y.2 <- x[[1]] == "U"
  y.1[y.1 == FALSE] <- 0 ; y.1[y.1 == TRUE] <- 1
  y.2[y.2 == FALSE] <- 0 ; y.2[y.2 == TRUE] <- 1
  y <- sum(y.1) + sum(y.2)
  if(y >= 4) {rna.x <- append(rna.x, paste0("5'-FAM-CGCG", paste(x[[1]], collapse = ""), "CGCG-IABkFQ-3'"))} 
  if(y <= 3) {rna.x <- append(rna.x, paste0("5'-FAM-AUAU", paste(x[[1]], collapse = ""), "AUAU-IABkFQ-3'"))} 
}
rna.2 <- rna.x

# Write out list
order.list <- data.frame(c(rna.1, rna.2))


