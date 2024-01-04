resample_phyloseq <- function(ps_object = ps, resamples = 1000, cores = 1, future_seed = 123){
  future::plan("multisession", workers = cores)
  boot_ps_list <- future.apply::future_replicate(n = resamples, expr = ps_boot(ps_object = ps_object), simplify = FALSE, future.seed = future_seed)
  return(boot_ps_list)
}
