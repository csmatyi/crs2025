# Author: Matthew Cserhati
# Date: July 1, 2025
# Email: csmatyi1@gmail.com

# Initialize: remove files, include libraries
options(warn=-1)
# set this to your own directory
invisible(install.packages("ggplot2",repos = "http://cran.us.r-project.org"))

install.packages(c("ggplot2","cluster","factoextra","ape","seqinr","msa"), dependencies=TRUE,repos = "http://cran.us.r-project.org")
library("ggplot2")
library("cluster")
library("factoextra")
library("ape")
library("seqinr")
library("msa")

# Transform clustalw alignment to sequence identity matrix
aln <- read.alignment(choose.files(filters = Filters[c("zip", "All"),]), "clustal")
n_species <- length(aln$nam)
idmx <- 1 - dist.alignment(aln,"identity")
idmx2 <- as.matrix(idmx,nrow=n_species,ncol=n_species)
mx <- idmx2

# OR:
# Read sequence identity matrix
mx <- as.matrix(read.table(choose.files(filters = Filters[c("zip", "All"),]),header=T,row.names=1,sep="\t",check.names=FALSE))

# Hopkins index
# 0 - 0.5: bad clustering
# 0.5 - 0.75: acceptable clustering
# 0.75 - 1: good clustering
result <- get_clust_tendency(mx, n = nrow(mx)-1, graph = FALSE)
print(result$hopkins_stat)

######################## Plots ###########################

message("Drawing plots...")

# Betweenness and withinness plots

### 1. Silhouette plot
jpeg("Silhouette.jpg")
fviz_nbclust(mx, kmeans, method = "silhouette", k.max=20) + theme_classic() # pam|kmeans|hcut|clara
dev.off()

### 2. Elbow plot
jpeg("Elbow.jpg")
fviz_nbclust(mx, kmeans, method = "wss", k.max=20) + theme_classic() # kmeans|pam (kmeans)
dev.off()

### 3. Heatmap
species <- row.names(mx)
myBreaks <- c(seq(0,1,by=0.01))
ceyy = cexx = 0.75

# normalize mx to 0-1
mx_hm <- mx
mx2 <- (mx_hm - min(mx_hm))/(max(mx_hm) - min(mx_hm))
mx_hm <- mx2

# color palette
clr = colorRampPalette(c("yellow","green","darkgreen"))(100) # plant studies
clr = colorRampPalette(c("green","white","red"))(100)
clr = colorRampPalette(c("white","yellow","red"))(100) # mitochondrial studies
clr = colorRampPalette(c("white","yellow","orange","red"))(100)
clr = gray.colors(100) # grayscale

# n is the estimated number of clusters/kinds/baramins
# if you change the estimate, restart from this point!
n <- 7

clusmeth="ward.D2" # other clustering methods: ward.D ward.D2 single median average mcquitty complete centroid
row.clusters = hclust(dist(mx_hm),method=clusmeth)
col.clusters = hclust(dist(mx_hm),method=clusmeth)
ctk <- cutree(row.clusters,k=n)
clus_clrs <- rainbow(length(as.matrix(unique(ctk))))
clus_clrs_vec <- clus_clrs[ctk[species[sort(row.clusters$order)]]]

### Heatmap code ###
# provide a name for the heatmap
heatmap_name=paste("heatmap_",clusmeth,"_",n,"_250724.jpg",sep="")
jpeg(filename = heatmap_name, height = 3000, width = 3000, units = "px", res=400) # topo.colors(100) 5500, 5000
h <- heatmap(mx_hm, symkey =F, symbreaks=F, scale="none", dendrogram = F, Rowv=T, Colv=T, col = clr, RowSideColors=clus_clrs_vec, breaks = myBreaks, border_color=NA, na.color="white", margin = c(15,15), # gray.colors(100)
        cexRow=cexx,cexCol=ceyy, key=T, trace="none", lmat=rbind( c(4, 3), c(2,1), c(0,0) ), lhei=c(1.8,6.5,1), hclustfun = function(d) hclust(d, method=clusmeth), # dendrogram="none",
        labCol=as.expression(lapply(colnames(mx), function(a) bquote(italic(.(a))))),labRow=as.expression(lapply(rownames(mx), function(a) bquote(italic(.(a))))))
invisible(dev.off())

######################### Stats files #########################

# stats and clusters files
# number of clusters
n <- 7

row.clusters = hclust(dist(mx),method=clusmeth)
ctk <- cutree(row.clusters,k=n)

# clusters from ctk
write.table(ctk, file="kclusters.txt", col.names=F, quote=F, sep="\t")

header = "cluster\tspecies\tmin\tmean\tmax\tSEM\tp-value\tneglog"
write(header, file="kstats.txt", sep="\t", append=T)

cluster_sizes = as.vector(table(ctk))
for (n_cluster in 1:n) {
  csize = cluster_sizes[n_cluster]
  if (csize >= 3) {
    m1 = as.matrix(mx[ctk == n_cluster,ctk == n_cluster])

    x = m1[upper.tri(m1)]
    ll = dim(m1)[1]

    m2 = as.matrix(cbind(mx[ctk != n_cluster,ctk == n_cluster],t(mx[ctk == n_cluster,ctk != n_cluster])))
    m2b = m2[!duplicated(colnames(m2))]

    t = t.test(x,m2b)
    pval = t$p.value
    nglog = -log10(pval)
    min = min(x)
    max = max(x)

    mean2 = sprintf("%.3f", mean(x))
    sem2 = sprintf("%.3f", sd(x)/sqrt(csize))
    min2 = sprintf("%.3f", min)
    max2 = sprintf("%.3f", max)
    pval2 = sprintf("%.3f", pval)
    nglog2 = sprintf("%.3f", nglog)

    stats = paste(n_cluster, ll, min2, mean2, max2, sem2, pval, nglog2, sep="\t")
    stats2 = gsub("\n","\t",stats)
    write(stats, file="kstats.txt", sep="\t", append = T)
  }
}
