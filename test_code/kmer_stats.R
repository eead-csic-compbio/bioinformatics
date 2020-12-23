#' @author Najla Nksouri <najlaksouri@gmail.com> 
#' Revised by Jacques van Helden <Jacques.van-Helden@univ-amu.fr>
#' @description Compare the statistics between query and background k_mer counts
#' @param query.oligo.file file containing k-mer counts computed from the query sequences, as produced by RSAT oligo-analysis.
#' This file must contain 
#' @param background.oligo.file file containing background k-mer occurrences
kmer.stats <- function(query.oligo.file,
                       background.oligo.file) {
  
  ## Check existence of input files
  if (!file.exists(query.oligo.file)) {
    stop("kmer.stats() error: query.oligo.file does not exist: ", query.oligo.file)
  }
  if (!file.exists(background.oligo.file)) {
    stop("kmer.stats() error: background.oligo.file does not exist: ", background.oligo.file)
  }
  
  ## Read tables and rename firts column
  oligo.query.freq.table <- read.table(query.oligo.file, sep = "\t", comment.char = ";", header = TRUE)
  names(oligo.query.freq.table)[1] <- "seq"
  
  oligo.genome.freq.table <- read.table(background.oligo.file, sep = "\t", comment.char = ";", header = TRUE)
  names(oligo.genome.freq.table)[1] <- "seq"
  
  oligo.table <- merge(
    oligo.query.freq.table,
    oligo.genome.freq.table,
    by=c("seq", "id"), suffixes=c("_query", "_background")
  )
  
  ## Compute the p-value of the observed occurrences in query frequencies
  nb.oligos <- nrow(oligo.table)
  N <- sum(oligo.table$occ_query) 
  
  ## Compute th number of occurrences expected in random background fragments
  ## of the same size as the query sequences.
  oligo.table$exp.occ <- oligo.table$obs_freq_background * N
  
  ## Test the over-representation: righ-tailed test
  oligo.table$pval.over <- pbinom(q = oligo.table$occ_query -1, size = N, 
                                  prob = oligo.table$obs_freq_background, lower.tail=FALSE)
  
  ## Compute log p-value to avoid the floating point limitation to 1e-320
  oligo.table$log.pval.over <- pbinom(q = oligo.table$occ_query -1, size = N, 
                                      prob = oligo.table$obs_freq_background, lower.tail=FALSE, log=TRUE)
  
  ## Test under-representation: left-tailed test
  oligo.table$pval.under <- pbinom(q = oligo.table$occ_query, size = N, 
                                   prob = oligo.table$obs_freq_background, lower.tail=TRUE)
  oligo.table$log.pval.under <- pbinom(q = oligo.table$occ_query, size = N, 
                                       prob = oligo.table$obs_freq_background, lower.tail=TRUE, log=TRUE)
  
  ## Two-sided test
  oligo.table$pval <- pmin(oligo.table$pval.under, oligo.table$pval.over) * 2
  oligo.table$log.pval <- pmin(oligo.table$log.pval.under, oligo.table$log.pval.over) * 2
  
  ## Compute the E-value (multiple testing correction).`
  ##  E-value indicates the expected nuimber of false positives 
  ## (i.e. how many oligos would expected to be be declared positive under 
  ## the null hypothesis)
  oligo.table$eval <- oligo.table$pval * nb.oligos
  
  ## Compute the significance
  ## sig = -log10(eval)
  oligo.table$sig <- -log10(nb.oligos) - oligo.table$log.pval/log(10)
  
  ## Compute maximal frequency for both files
  max.axis <- max(c(oligo.table$obs_freq_query, oligo.table$obs_freq_background))
  
  ## Log2 ratio observed/expected
  oligo.table$freq.ratio <- oligo.table$obs_freq_query/oligo.table$obs_freq_background
  oligo.table$log2.ratio <- log2(oligo.table$freq.ratio)
  
  return(oligo.table)
}
  
