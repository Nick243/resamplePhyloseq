#' Perform bootstrap resampling of an existing phyloseq object.
#' @description The function assumes the phyloseq object has at least sample_data and otu_table elements and that taxa are rows. If your taxa are columns, transpose before running function.
#'
#' @param ps_object An existing phyloseq object
#' @param resamples The number of bootstrap samples to generate (default = 1000)
#' @param cores Number of available cores
#' @param future_seed Seed for reproducing results (passed to future.seed)
#'
#' @return A list of bootstrap resampled phyloseq objects
#' @export
#'
#' @examples
resample_phyloseq <- function(ps_object = ps, resamples = 1000, cores = 1, future_seed = 123){
  future::plan("multisession", workers = cores)
  boot_ps_list <- future.apply::future_replicate(n = resamples, expr = ps_boot(ps_object = ps_object), simplify = FALSE, future.seed = future_seed)
  return(boot_ps_list)
}
