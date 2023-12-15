#' Replace sequences with their representatives
#'
#' Similar as `seqdef()`. It returns a `stslist` object but with elements replaced by their cluster representatives.
#'
#' @param data 	A data frame, matrix, or character string vector containing sequence data (tibble will be converted with as.data.frame).
#' @param clusters A data frame, containing three columns—the index of original data, group memberships, and group size. Or the cluster data frame returned by `seqcluster` function.
#' @param var The list of columns containing the sequences. Default is NULL, i.e. all the columns. The function detects automatically whether the sequences are in the compressed (successive states in a character string) or extended format.
#' @param alphabet Optional vector containing the alphabet (the list of all possible states). Use this option if some states in the alphabet don't appear in the data or if you want to reorder the states. The specified vector MUST contain AT LEAST all the states appearing in the data. It may possibly contain additional states not appearing in the data. If NULL, the alphabet is set to the distinct states appearing in the data as returned by the seqstatl function. See details.
#' @param idxs A integer or an array of integers. The Default is 1:10, meaning replacing sequences in the 10 largest clusters with their representatives. If idxs=0, then all clusters will be replaces.
#' @param ... options passed to the `seqdef` function for handling input data that is not in STS format.
#'
#' @return An object of class `stslist`.
#' @export
#'
#' @examples
#' data(mvad)
#' mvad.alphab <- c("employment", "FE", "HE", "joblessness", "school", "training")
#' mvad.seq <- seqdef(mvad, 17:86, xtstep = 6, alphabet = mvad.alphab)
#' seqdistance <- seqdist(mvad.seq, method="HAM")
#' clusters <- seqcluster(seqdistance, h=1.5, cmethod='complete')
#' mvad.new.seq <- seqrep_replace(mvad, clusters, var=17:86, alphabet, idxs=1:20, xtsetp=6)
seqrep_replace <- function(data, clusters, var=NULL, alphabet, idxs=1:10, ...){
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
    data[seq.idx, var] <- fseq
  }

  new.seq <- seqdef(data, var=var, alphabet=alphabet, ...)
  return(new.seq)
}
