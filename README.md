
<!-- README.md is generated from README.Rmd. Please edit that file -->

# resamplePhyloseq

<!-- badges: start -->
<!-- badges: end -->

The goal of **resamplePhyloseq** is to perform bootstrap resampling
(i.e., sampling with replacement to the same number of samples as the
original object) of an existing phyloseq object. A list of bootstrap
resampled phyloseq objects is returned. **Users can then iterate
additional functions over the list elements** to obtain bootstrap
resampled quantities of interest. Parallel processing is implemented via
the future and future.apply packages.

## Installation

You can install the development version of resamplePhyloseq from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("Nick243/resamplePhyloseq")
```

## Example

This basic example highlights the core functionality of the
resamplePhyloseq package.

First, an example phyloseq object is obtained from the phyloseq package.
Below we a phyloseq object with 31 samples, 217 taxa (where taxa are
rows), and 9 sample variables is generated.

``` r
library(resamplePhyloseq)
library(phyloseq)
library(dplyr)
library(ggplot2)
library(Maaslin2)
data(enterotype)

ps <- phyloseq::subset_samples(enterotype, ClinicalStatus %in% c("healthy", "obese"))
ps <- phyloseq::subset_taxa(ps, phyloseq::taxa_sums(ps) > 0)
phyloseq::sample_data(ps)$ClinicalStatus <- factor(phyloseq::sample_data(ps)$ClinicalStatus, levels = c("healthy", "obese"))

ps <- subset_taxa(ps, Genus != "Bacteria")
ps <- subset_taxa(ps, Genus != "-1")
ps <- transform_sample_counts(ps, function(x) x / sum(x))
ps
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 217 taxa and 31 samples ]
#> sample_data() Sample Data:       [ 31 samples by 9 sample variables ]
#> tax_table()   Taxonomy Table:    [ 217 taxa by 1 taxonomic ranks ]
```

<br>

Then we run the *resample_phyloseq()* function to generate a list of
bootstrap resampled phyloseq objects. Here I only generate n=100
resamples to speed up the computations, but **larger values should be
considered** for many/most applications. Run parallel::detectCores() to
determine the number of available cores. For projects with a large
number of samples, or when many resamples are needed, consider using
using parallel::detectCores() - 1.

We now have a list of resampled phyloseq objects and the number of times
a given sample shows up across resamples varies (i.e., we are using the
sample data to approximate the population of interest).

``` r
ps_boot_list <- resample_phyloseq(ps_object = ps, resamples = 100, cores = 7, future_seed = 243)
ps_boot_list[1:2]
#> [[1]]
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 217 taxa and 31 samples ]
#> sample_data() Sample Data:       [ 31 samples by 9 sample variables ]
#> 
#> [[2]]
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 217 taxa and 31 samples ]
#> sample_data() Sample Data:       [ 31 samples by 9 sample variables ]

table(sample_data(ps_boot_list[[1]])$Sample_ID)
#> 
#>   AM.AD.2 AM.F10.T2   DA.AD.1   DA.AD.2   DA.AD.4   ES.AD.2   ES.AD.4   FR.AD.4 
#>         1         4         1         2         1         1         2         1 
#>   FR.AD.5   FR.AD.6   FR.AD.8   JP.AD.1   JP.AD.2   JP.AD.5   JP.AD.6   JP.AD.7 
#>         1         2         1         3         1         3         2         2 
#>   JP.IN.4 
#>         3
table(sample_data(ps_boot_list[[2]])$Sample_ID)
#> 
#>   AM.AD.1 AM.F10.T1   DA.AD.1   DA.AD.2   DA.AD.3   DA.AD.4   ES.AD.2   FR.AD.1 
#>         1         1         3         3         2         2         1         1 
#>   FR.AD.2   FR.AD.4   FR.AD.5   FR.AD.6   FR.AD.7   JP.AD.2   JP.AD.4   JP.AD.5 
#>         1         3         2         2         2         2         2         1 
#>   JP.AD.8   JP.IN.3 
#>         1         1
```

<br>

Finally, an example is provided below where an additional user-written
function is developed and then iterated over the resampled phyloseq
objects to obtain bootstrap resampled estimates of uncertainty. Here I
use the default settings (other than the TSS normalization since we
already have relative abundance data/proportions) in the Maaslin2
package to generate bootstrap estimates of the differential rankings
based on the FDR corrected p-values. For a more detailed description of
differential ranking see this [blog
post](https://www.nicholas-ollberding.com/post/bootstrap-resampling-for-ranking-differentially-abundant-taxa/).
**Any appropriate function could be applied here where bootstrap
estimates are needed.** The variation seen here is a bit extreme given
the example used. I would not expect such wide 95% CIs in most settings
based on my experience using deeper sequencing and larger sample sizes.

``` r
get_ranks <- function(ps_object){
  m <- Maaslin2(
    input_data = data.frame(otu_table(ps_object)),
    input_metadata = data.frame(sample_data(ps_object)),
    output = "/Users/olljt2/desktop/maaslin2_output",
    min_abundance = 0.0,
    min_prevalence = 0.0,
    normalization = "NONE",
    transform = "LOG",
    analysis_method = "LM",
    fixed_effects = c("ClinicalStatus"),
    correction = "BH",
    standardize = FALSE,
    cores = 1,
    plot_heatmap = FALSE,
    plot_scatter = FALSE)

  mas_out <- m$results |>
    dplyr::filter(metadata == "ClinicalStatus") |>
    dplyr::select(-qval, -name, -metadata, -N, -N.not.zero) |>
    dplyr::arrange(pval) |>
    dplyr::mutate(qval = p.adjust(pval, method = "fdr"),
                  p_rank = rank(pval))
  return(mas_out)
}

boot_ranks <- lapply(ps_boot_list, get_ranks) 
boot_ranks_df <- bind_rows(boot_ranks, .id = "column_label")
```

``` r
pval_ranks_df <- boot_ranks_df |>
  dplyr::group_by(feature) |>
  dplyr::summarise(median_p_rank = round(quantile(p_rank, prob = (0.5)), 0),
                   lcl_p_rank = round(quantile(p_rank, prob = (0.025)), 0),
                   ucl_p_rank = round(quantile(p_rank, prob = (0.975)), 0),
                   min_p_rank = min(p_rank),
                   max_p_rank = max(p_rank),
                   p_rank_95_ci = paste("(", lcl_p_rank, "; ", ucl_p_rank, ")", sep = "")) |>
  dplyr::select(feature, median_p_rank, p_rank_95_ci, min_p_rank, max_p_rank) |>
  dplyr::ungroup() |>
  dplyr::arrange(median_p_rank, min_p_rank)

head(pval_ranks_df)
#> # A tibble: 6 × 5
#>   feature         median_p_rank p_rank_95_ci min_p_rank max_p_rank
#>   <chr>                   <dbl> <chr>             <dbl>      <dbl>
#> 1 Pyramidobacter              3 (1; 53)             1         88  
#> 2 Macrococcus                 7 (1; 44)             1         89  
#> 3 Actinomycetales             9 (1; 67)             1         94  
#> 4 Rhodococcus                 9 (1; 91)             1        104  
#> 5 Geobacter                  13 (2; 72)             1         86.5
#> 6 Thermincola                14 (3; 75)             2.5       97
```
