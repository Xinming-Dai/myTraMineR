library(TraMineR)
library(tidyverse)
library(seriation)
data(mvad)
# frequency plot
mvad.alphab <- c("employment", "FE", "HE", "joblessness", "school", "training")
mvad.seq <- seqdef(mvad, 17:86, xtstep = 6, alphabet = mvad.alphab)
seqfplot(mvad.seq, idxs = 1:20)
method <- "HAM"
method <- "OM"
method <- "LCS"
seqdistance <- seqdist(mvad.seq, method=method)
n <- dim(seqdistance)[1]

#=====reorder for image====
method <- "PCA" # best
method <-  "AOE" # best
method <- "Identity"
method <- "LLE"

order_list <- seriate(seqdistance, method=method) 
order_vector <- get_order(order_list)
pimage(seqdistance, order_list, main = paste("Reordered by seriation with method", method))
# reorder seqdistance matrix
reorder <- order(order_vector, decreasing = TRUE)
seqdistance_reordered <- seqdistance[reorder, reorder]
head(seqdistance_reordered)

#======complete linkage clustering======
seqcluster <- function(seqdistance, h, cmethod){
  cluster_object <- hclust(as.dist(seqdistance), method=cmethod)
  cluster_cut <- cutree(cluster_object, h = h)
  clusters <- data.frame(index=as.numeric(names(cluster_cut)),
                         cluster=cluster_cut)
  clusters <- clusters[order(clusters$cluster),]
  Freq <- as.data.frame(table(clusters$cluster))
  clusters["Freq"] <- rep(Freq$Freq, Freq$Freq)
  clusters <- clusters[order(clusters$Freq, decreasing=T),]
  return(clusters)
}

h <- 1.5
cmethod <- 'complete'
clusters <- seqcluster(seqdistance, h, cmethod)

# reorder the original distance matrix
cluster_seqdistance <- seqdistance[clusters$index, clusters$index]
pimage(cluster_seqdistance, main="Clustered using complete linkage clustering")
pimage(cluster_seqdistance[1:50, 1:50], main="Clustered using complete linkage clustering (first 50 data)")

# first most frequent sequences
frequency <- sort(table(clusters$cluster), decreasing=TRUE)/n
sum(frequency[1:20])

# ===========function 1 frequency plot==========
# 1. replace all elements with their medoids; build a new dataframe, use rep-->frequency; make a plot use seqfplot
# use seqtab
# improvement, if I can change the seqtab data directly. then, I can plot it directly.
# users may also want to give stslist data instead of a dataframe, so I can write a generic function for both classes
seqrep_replace <- function(data, clusters, var=NULL, alphabet, idxs=1:10, ...){
  # replace the first a few largest sequence clusters with their representatives
  # return a state sequence object—stslist
  if(is.null(idxs) || 
     (length(idxs) > 1 && min(idxs) < 1) || 
     any(idxs < 0)){
    stop("idxs should be a non negative integer or a strictly positive vector.")
  }
  seqdata <- seqdef(data, var=var, alphabet=alphabet, ...)
  cluster.idx <- unique(clusters[,c("cluster", "Freq")])
  nbuseq <- nrow(cluster.idx)
  if (idxs[1] == 0 || max(idxs) > nbuseq) {
    idxs <- 1:nbuseq
  }
  cluster.idx <- cluster.idx[idxs,]
  
  for(i in idxs){
    condition <- cluster.idx[i,]
    seq.idx <- clusters[clusters$cluster == condition$cluster & clusters$Freq == condition$Freq ,]$index
    fseq <- seqtab(seqdata[seq.idx,], idxs=1)[1,]
    data[seq.idx, var] <- fseq # replace the original dataframe with the most frequent sequences
  }
  
  new.seq <- seqdef(data, var=var, alphabet=alphabet, ...)
  return(new.seq)
}

idxs <- 1:20
var <- 17:86
alphabet <- c("employment", "FE", "HE", "joblessness", "school", "training")

mvad.new.seq <- seqrep_replace(data, clusters, var, alphabet, idxs, xtsetp=6)
seqfplot(mvad.new.seq, idxs=idxs)

# 2. If you want to make a plot for original data, I want to use seqiplot, 
# just make a plot on the sorted dataframe according to cluster size. use sortv in seqiplot

# =========== function 2: cut height for the top p percentile most frequent sequences========
getheight <- function(seqdistance, cmethod, idxs, fq){
  h <- 0.5
  clusters <- seqcluster(seqdistance, h, cmethod)
  frequency <- sort(table(clusters$cluster), decreasing=TRUE)/n
  frequency <- sum(frequency[idxs])
  if(fq>1|fq<0) stop("frequency should be in range [0,1]")
  while(frequency<fq){
    h <- h + 1
    clusters <- seqcluster(seqdistance, h, cmethod)
    frequency <- sort(table(clusters$cluster), decreasing=TRUE)/n
    frequency <- sum(frequency[idxs])
  }
  h.fq <- list(h=h,frequency=frequency)
  return(h.fq)
}
fq <- 0.3
getheight(seqdistance, cmethod, idxs, fq)
