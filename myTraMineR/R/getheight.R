#' Get the cut height for cutting a hierarchical clustering tree
#'
#' The returned cut height `h` will satisfy that the first a few most frequent groups, specified by `idxs`, will represent at least `fq` of the total.
#'
#' @param seqdistance A distance matrix or a distance array returned by TraMineR::seqdist() function.
#' @param cmethod The agglomeration method to be used. This should be (an unambiguous abbreviation of) one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC). The default clustering method is "complete", complete-linkage clustering.
#' @param idxs A integer or an array of integers. The Default is 1:10, meaning replacing sequences in the 10 largest clusters with their representatives. If idxs=0, then all clusters will be replaces.
#' @param fq A numeric number within the range \[0,1\].
#'
#' @return The cut height—a numeric number and the frequency, meaning how much these clusters will represent for the total.
#' @export
#'
#' @examples
#' getheight(seqdistance, fq=0.3)
#' getheight(seqdistance, fq=0.3, idxs=0)
getheight <- function(seqdistance, cmethod="complete", idxs=1:10, fq){
  h <- 0.5
  n <- nrow(seqdistance)
  if(is.null(idxs) ||
     (length(idxs) > 1 && min(idxs) < 1) ||
     any(idxs < 0)){
    stop("idxs should be a non negative integer or a strictly positive vector.")
  }
  clusters <- seqcluster(seqdistance, h, cmethod)
  frequency <- sort(table(clusters$cluster), decreasing=TRUE)/n
  nbuseq <- length(frequency)
  if (idxs[1] == 0 || max(idxs) > nbuseq) {
    idxs <- 1:nbuseq
  }
  frequency <- sum(frequency[idxs])
  if(fq>1|fq<0) stop("frequency should be in range [0,1]")
  while(frequency<fq){
    h <- h + 1
    clusters <- seqcluster(seqdistance, h, cmethod)
    frequency <- sort(table(clusters$cluster), decreasing=TRUE)/n
    nbuseq <- length(frequency)
    if (idxs[1] == 0 || max(idxs) > nbuseq) {
      idxs <- 1:nbuseq
    }
    frequency <- sum(frequency[idxs])
  }
  h.fq <- list(h=h,frequency=frequency)
  if (h==0.5){
    message("There is no need to cluster.")
    h.fq$h <- 0
    return(h.fq)
  }else{
    return(h.fq)
  }
}
