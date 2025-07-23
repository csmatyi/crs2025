### EXERCISE #1 ###

# Manipulating a DNA sequence as a string

install.packages("stringr")
library(stringr)

dna <- "GTCGTCGTAGTAGGTTTATTATTCAA"

# calculate GC%
c_count <- stringr::str_count(dna, "C")
g_count <- stringr::str_count(dna, "G")
dna_len <- stringr::str_length(dna)
gc <- (c_count + g_count)*100/dna_len
gc

stringr::str_to_lower(dna)

# mutations
# substitution
sequence_chars <- unlist(strsplit(dna, split = ""))
sequence_chars
sequence_chars[4] <- 'A'
dna_mut <- paste(sequence_chars, collapse = "")
dna_mut

# OPTIONAL: view DNA sequences via Biostrings
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Biostrings")
BiocManager::install("msa")
library(msa)

Biostrings::DNAString(dna)
Biostrings::DNAString(dna_mut)
# notice the difference in color at base #4

# run multiple alignment
seqs <- DNAStringSet(c(dna, dna_mut))
alignment <- msa(seqs, method = "ClustalW")
alignment


# inversion
paste(
  rev(
    strsplit(dna,"")[[1]]
    ),
  collapse="")

# insertion of AAAAA after base #7
dna2 <- paste0(
  stringr::str_sub(dna,1,7),
  "AAAAA",
  stringr::str_sub(dna,8,25),
  collapse=""
  )
dna2

# deletion of 5 bp at position #5
dna3 <- paste0(
  stringr::str_sub(dna,1,7),
  stringr::str_sub(dna,12,25),
  collapse=""
  )
dna3



# Go back to Bioinfo_R_tutorial.R

### EXERCISE #2 ###

# Write function for Hardy-Weinberg equilibrium
hardy_weinberg <- function(p) {
  # a little bit of error-handling
  try(if ((p < 0) || (p > 1)) stop("p must be between 0 and 1!"))
  q <- 1 - p
  homozygous <- 2*p*q / (p*p + 2*p*q + q*q)
  return(homozygous)
}

p1 <- 0.21
p1_homoz <- hardy_weinberg(p1)
print(paste0("The proportion of homozygots is ",p1_homoz))



# Go back to Bioinfo_R_tutorial.R

### EXERCISE #3 ###

# RNA-seq experiment

install.packages("dplyr")
library(dplyr)

rna_exp <- data.frame(
  genes = c("TUB", "TP53", "WNT2", "CYCD1"),
  tumor1 = c(5.9,22.7,0.1,50.9),
  tumor2 = c(5.7,24.8,0.2,48.6),
  tumor3 = c(5.9,28.1,0.1,61.1),
  control1 = c(5.6,11.1,0.3,99.8),
  control2 = c(6.0,10.7,0.0,116.7)
)
rna_exp

# calculate log2 fold change between tumor and control
log2fc <- function(vec,n=3) {
  avg_trt <- mean(vec[1:n])
  m <- n+1

  avg_ctrl <- mean(vec[m:length(vec)])
  l2fc <- avg_trt/avg_ctrl
  return(l2fc)
}

# calculate values
l2fcs <- apply(as.matrix(rna_exp[,2:6]),1,log2fc)
rna_exp <- cbind(rna_exp,log2fc=l2fcs)
rna_exp

# get those genes which are:
# a. induced (log2 fc >= 2)
rna_exp %>% filter(log2fc >= 2) %>% select(genes)

# b. repressed (log2 fc <= 0.5)
rna_exp %>% filter(log2fc <= 0.5) %>% select(genes)

# c. stay the same (0.5 <= log2 fc <= 2)
rna_exp %>% filter(log2fc < 2 & log2fc > 0.5) %>% select(genes)



# Go back to Bioinfo_R_tutorial.R

### EXERCISE #4 ###
# Using the hybridogram package on slightly a larger data set

install.packages("hybridogram")
# OR
install.packages("hybridogram_0.3.2.tar.gz",repos=NULL,dependencies=TRUE,type="source")
library("hybridogram")
installed.packages()

V1 <- c("Phoca largha","Phoca largha","Phoca caspica","Odobenus rosmarus")
V2 <- c("Phoca vitulina","Phoca caspica","Pusa hispida","Odobenus sp.")
V3 <- c(2,3,3,3)
hybrid_data <- data.frame(V1,V2,V3)

C1 <- c(1,2,3)
C2 <- c("No hybrid","Hybrid with same 3rd species","Documented hybrid")
codes <- data.frame(C1,C2)

jpeg("hybridogram.jpg")
hybridogram(hybrid_data, codes)
dev.off()



### EXERCISE #5 ###

# Use the matrixcut package to find clusters in matrixes

install.packages("matrixcut")
library("matrixcut")

# R packages may come with their own datasets
primates <- matrixcut::primates
primates

# use matrixcut to get clusters
matrixcut::matrixcut(primates,0.9)

# create matrixcut image
jpeg("primate_components.jpg")
matrixcut::componentplot(primates,0,1)
dev.off()
