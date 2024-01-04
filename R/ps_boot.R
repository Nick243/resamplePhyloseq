ps_boot <- function(ps_object){

  '%>%' <- magrittr::'%>%'

  ps_ids <- phyloseq::sample_names(ps_object)
  boot_ids <- sample(ps_ids, size = length(ps_ids), replace = TRUE)
  boot_df <- data.frame(seq_id = boot_ids)

  meta_df <- data.frame(phyloseq::sample_data(ps_object)) %>%
    tibble::rownames_to_column(var = "seq_id")

  join_df <- dplyr::left_join(boot_df, meta_df, by = "seq_id")

  otu_df <- data.frame(t(phyloseq::otu_table(ps_object))) %>%
    tibble::rownames_to_column(var = "seq_id")

  join_otu <- dplyr::left_join(boot_df, otu_df, by = "seq_id")


  join_df <- join_df %>%
    dplyr::mutate(boot_id = paste("S", dplyr::row_number(), sep = "")) %>%
    dplyr::select(-seq_id) %>%
    tibble::column_to_rownames(var = "boot_id")

  join_otu <- join_otu %>%
    dplyr::mutate(boot_id = paste("S", dplyr::row_number(), sep = "")) %>%
    dplyr::select(-seq_id) %>%
    tibble::column_to_rownames(var = "boot_id") %>%
    t() %>%
    data.frame()

  boot_ps <- phyloseq::phyloseq(phyloseq::sample_data(join_df),
                                phyloseq::otu_table(join_otu, taxa_are_rows = TRUE))

  return(boot_ps)
}
