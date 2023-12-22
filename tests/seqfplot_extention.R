library(TraMineR)
library(seriation)

seqcluster <- function(seqdistance, h=1.5, cmethod='complete'){
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

