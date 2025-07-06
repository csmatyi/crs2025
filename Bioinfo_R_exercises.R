# Exercise #1
# Manipulating a DNA sequence as a string

library(stringr)
dna <- "GTCGTCGTAGTAGGTTTATTATTCG"
c_count <- stringr::str_count(dna, "C")
g_count <- stringr::str_count(dna, "G")
dna_len <- stringr::str_length(dna)
gc <- (c_count + g_count)/dna_len
gc

stringr::str_to_lower(dna)

# mutations
# substitution
sequence_chars <- unlist(strsplit(dna, split = ""))
sequence_chars[4] <- 'A'
dna <- paste(sequence_chars, collapse = "")
dna

# inversion
paste(rev(strsplit(dna,"")[[1]]),collapse="")

# insertion
dna2 <- paste0(stringr::str_sub(dna,1,7),"AAAAA",stringr::str_sub(dna,8,25),collapse="")
dna2

# deletion
dna3 <- paste0(stringr::str_sub(dna,1,7),stringr::str_sub(dna,12,25),collapse="")
dna3



# Exercise #2
# Write function for Hardy-Weinberg equilibrium
hardy_weinberg <- function(p) {
  # a little bit of error-handling
  try(if ((p < 0) || (p > 1)) stop("p must be between 0 and 1!"))
  q <- 1 - p
  homozygous <- 2*p*q/(p*p+2*p*q+q*q)
  return(homozygous)
}
p1 <- 0.21
p1_homoz <- hardy_weinberg(p1)
print(paste0("The proportion of homozygots is ",p1_homoz))



# Exercise #3
# RNA-seq experiment

library(dplyr)

rna_exp <- data.frame(
  genes = c("TUB", "TP53", "WNT2", "CYCD1"),
  tumor1 = c(5.9,22.7,0.1,50.9),
  tumor2 = c(5.7,24.8,0.2,48.6),
  tumor3 = c(5.9,28.1,0.1,61.1),
  control1 = c(5.6,11.1,0.3,99.8),
  control2 = c(6.0,10.7,0.0,116.7)
)

log2fc <- function(vec,n=3) {
  avg_trt <- mean(vec[1:n])
  m <- n+1
  avg_ctrl <- mean(vec[m:length(vec)])
  l2fc <- avg_trt/avg_ctrl
  return(l2fc)
}

l2fcs <- apply(as.matrix(rna_exp[,2:6]),1,log2fc)
rna_exp <- cbind(rna_exp,log2fc=l2fcs)
rna_exp

# get those genes which are:
# a. induced
rna_exp %>% filter(log2fc >= 2) %>% select(genes)
# b. repressed
rna_exp %>% filter(log2fc <= 0.5) %>% select(genes)
# c. stay the same
rna_exp %>% filter(log2fc < 2 & log2fc > 0.5) %>% select(genes)



# Exercise #4
# Using the hybridogram package on slightly a larger data set

install.packages("hybridogram")
library("hybridogram")

V1 <- c("Phoca largha","Phoca largha","Phoca caspica","Odobenus rosmarus")
V2 <- c("Phoca vitulina","Phoca caspica","Pusa hispida","Odobenus sp.")
V3 <- c(2,3,3,3)
hybrid_data <- data.frame(V1,V2,V3)
C1 <- c(1,2,3)
C2 <- c("No hybrid","Hybrid with same 3rd species","Documented hybrid")
codes <- data.frame(C1,C2)
hybridogram(hybrid_data, codes)



# Exercise #5
# Use the matrixcut package to find clusters in matrixes

install.packages("matrixcut")
library("matrixcut")

# R packages may come with their own datasets
primates <- matrixcut::primates
primates

matrixcut::matrixcut(primates,0.9)
jpeg("primate_components.jpg") # save as an image
matrixcut::componentplot(primates,0,1)
dev.off()