#!/usr/bin/env Rscript

library(argparse)
library(furrr)
library(tidyverse)
source("/bms-dms/src/model_utils.R")

plan(multicore, workers = 40)

parser <- ArgumentParser()
parser$add_argument("-i", "--input", required=TRUE, help="Path to the input summary statistics")
parser$add_argument("-o", "--output", required=TRUE, help="Path to the output midpoints file")
args <- parser$parse_args()

output_file <- args$output
marginals <- read_tsv(args$input)

weights <- marginals %>%
    group_by(chunk, pos, aa) %>%
    mutate(bin = case_when(condition_conc == 25 ~ 0.125,
                           condition_conc == 50 ~ 0.375,
                           condition_conc == 75 ~ 0.625,
                           TRUE ~ 0.875))

weights_nest <- weights %>% nest(data = c(-chunk, -pos, -aa))
metadata_df <- weights_nest %>% select(-data) %>% transpose()

midpoints <- future_map2_dfr(.x = weights_nest$data,
                             .y = metadata_df,
                             ~generate_resamples(.x$log2Marginal,
                                                 .x$log2MarginalError,
                                                 .y,
                                                 num = 1000))

wt_scores <- midpoints %>%
    filter(aa == "WT") %>% 
    group_by(chunk) %>%
    summarize("WT score" = median(score_mean),
              "WT score standard error" = median(score_se))

midpoints_test <- midpoints %>%
    filter(aa != "WT") %>%
    left_join(wt_scores, by = "chunk") %>%
    mutate(midpoint_shift = score_mean - `WT score`,
           midpoint_shift_se = sqrt(score_se^2 + `WT score standard error`^2),
           statistic = midpoint_shift/midpoint_shift_se,
           p.value = (1 - pnorm(abs(statistic))) * 2,
           p.adj = p.adjust(p.value, method = "BH"))

write_tsv(midpoints_test,
    output_file)