# Copyright (c) 2025 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit



#' Rarefy OTU counts.
#' 
#' Sub-sample OTU observations such that all samples have an equal number.
#' If called on data with non-integer abundances, values will be re-scaled to 
#' integers between 1 and `depth` such that they sum to `depth`.
#' 
#' @inherit documentation
#' 
#' @param depth   How many observations to keep per sample. When 
#'        `0 < depth < 1`, it is taken as the minimum percentage of the 
#'        dataset's observations to keep. Ignored when `n_samples` is 
#'        specified. Default: `0.1`
#'
#' @param n_samples   The number of samples to keep. When `0 < n_samples < 1`, 
#'        it is taken as the percentage of samples to keep. If negative, that 
#'        number of samples is dropped. If `0`, all samples are kept. 
#'        If `NULL`, then `depth` is used instead. Default: `NULL`
#'     
#' @param seed   An integer seed for randomizing which observations to keep or 
#'        drop. If you need to create different random rarefactions of the same 
#'        data, set the seed to a different number each time.
#'     
#' @param times   How many independent rarefactions to perform. If set, 
#'        `rarefy()` will return a list of matrices. The seeds for each matrix
#'        will be sequential, starting from `seed`.
#' 
#' 
#' @return An integer matrix.
#' 
#' @export
#' @examples
#'     # Create an OTU matrix with 4 samples (A-D) and 5 OTUs.
#'     counts <- matrix(
#'       data     = c(4,0,3,2,6,0,8,0,0,5,0,9,0,0,7,0,10,0,0,1),
#'       nrow     = 5,
#'       dimnames = list(paste0('OTU', 1:5), LETTERS[1:4]) )
#'     counts
#'     colSums(counts)
#'     
#'     counts <- rarefy(counts, depth = 14)
#'     counts
#'     colSums(counts)
#' 
rarefy <- function (
    counts, 
    depth     = 0.1, 
    n_samples = NULL, 
    seed      = 0, 
    times     = NULL,
    cpus      = n_cpus() ) {
  
  validate_args()
  
  if (any(counts %% 1 > 0))
    stop('Non-integer counts cannot be rarefied.')
  
  n_rows <- nrow(counts)
  n_cols <- ncol(counts)
  dnames <- dimnames(counts)
  
  counts <- matrix(
    data     = as.integer(counts), 
    nrow     = n_rows, 
    ncol     = n_cols,
    dimnames = dnames )
  
  
  # Set target depth according to number/pct of samples to keep/drop.
  if (!is.null(n_samples)) {
    if (n_samples == 0)     n_samples <- n_cols             # Keep all
    if (abs(n_samples) < 1) n_samples <- n_samples * n_cols # Keep/drop percentage
    if (n_samples <= -1)    n_samples <- n_cols + n_samples # Drop n_samples
    n_samples <- max(1, floor(n_samples))                   # Keep at least one
    target    <- rev(sort(colSums(counts)))[[n_samples]]
  }
  
  # Depth is given as minimum percent of observations to keep.
  else if (depth < 1) {
    depths <- colSums(counts)   # observations per sample
    target <- (sum(depths) * depth) / length(depths)
    target <- min(depths[depths >= target])
  }
  
  # Depth is given as observations per sample to keep.
  else {
    target <- depth
  }
  
  target <- as.integer(target)
  
  
  if (is.null(times)) {
    
    result_mtx <- matrix(0L, n_rows, n_cols, FALSE, dnames)
    .Call(C_rarefy, counts, target, seed, cpus, result_mtx)
    
  } else {
    
    seeds <- ((seed + 2**31 - 1 + seq_len(times)) %% 2**32) - 2**31
    
    lapply(seeds, function (seed) {
      result_mtx <- matrix(0L, n_rows, n_cols, FALSE, dnames)
      .Call(C_rarefy, counts, target, seed, cpus, result_mtx)
    })
  }
  
}


