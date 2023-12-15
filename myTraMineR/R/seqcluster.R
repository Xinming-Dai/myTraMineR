#' Cluster sequences
#'
#' This function clusters sequences using hierarchical cluster and returns group memberships of sequences.
#'
#' @param seqdistance A distance matrix or a distance array returned by TraMineR::seqdist() function.
#' @param h Numeric scalar or vector with heights where the tree should be cut. The default is 1.5.
#' @param cmethod The agglomeration method to be used. This should be (an unambiguous abbreviation of) one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC). The default clustering method is "complete", complete-linkage clustering.
#'
#' @return `seqcluster` returns a data frame with three columns—the index of original data, group memberships, and group size.
#' @export
#'
#' @examples
#' data(mvad)
#' mvad.alphab <- c("employment", "FE", "HE", "joblessness", "school", "training")
#' mvad.seq <- seqdef(mvad, 17:86, xtstep = 6, alphabet = mvad.alphab)
#' seqdistance <- seqdist(mvad.seq, method="HAM")
#' clusters <- seqcluster(seqdistance, h=1.5, cmethod='complete')
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
