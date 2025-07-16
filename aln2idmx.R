library("seqinr")

# Run these commands only if running at the command line, otherwise skip down if you are in R Studio.
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=3) {
  stop("Wrong parameters. Run as: Rscript aln2idmx.R <alignment file> <sequence id matrix> <alignment format>", call.=FALSE)
}

# Create identity matrix from MAFFT alignment
# args[3] can be fasta or clustal

# 1. read in the alignment file from MEGA or MAFFT
aln <- read.alignment(args[1],args[3])
# 2. run the dist.alignment method from seqinr to create an 'identity' matrix (same as seq sim matrix)
aln.idmx <- 1-as.matrix(dist.alignment(aln,"identity"))
# 3. write the matrix to a file as input for the next R script
write.table(aln.idmx,args[2],quote=F,sep="\t")
