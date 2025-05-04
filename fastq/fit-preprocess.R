#!/usr/bin/env Rscript

## Usage:
## Rscript fit-preprocess.R <out_file> <oligo_map> <samp_prop> <bc_file_dir>

library(magrittr)
library(data.table)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
stop_aliases <-  c("*", "X", "Stop", "stop", "x")

out_file <- args[1]
oligo_map <- fread(args[2], col.names = c("barcode", "oligo"))
samp_prop <- read_tsv(args[3], col_names = TRUE)
bc_file_dir <- args[4]

bc_files <- list.files(bc_file_dir, pattern = ".rna-bcs.tsv", full.names = TRUE)
names(bc_files) <- gsub(str_c(bc_file_dir, "/|.rna-bcs.tsv"), "", bc_files, perl = TRUE)
bcs <- bc_files %>%
    map_dfr(fread, .id = "sample",
    col.names = c("count", "barcode"))

# Remove barcodes represented more than once
setkey(oligo_map, barcode)
oligo_map <- oligo_map[, if (.N == 1) .SD, by = key(oligo_map)]

# Join map to counts
setkey(bcs, barcode)
bc_oligo_join <- data.table::merge.data.table(bcs,
    oligo_map,
    by = "barcode",
    all.x = TRUE)
bc_oligo_join <- data.table::merge.data.table(bc_oligo_join,
    samp_prop,
    by = "sample",
    all.x = TRUE)

# Format and write to out_file
mapped_counts <- bc_oligo_join %>%
    filter(!is.na(oligo)) %>%
    separate(oligo, c("lib", "chunk", "wt_aa", "pos",
        "mut_aa", "wt_codon", "mut_codon"), "_") %>%
    group_by(sample, chunk) %>%
    mutate(condition = str_c(condition, "_", condition_conc),
        mut_aa = if_else(wt_aa == mut_aa | is.na(mut_aa), "WT", mut_aa),
        stop_counts = log(sum(count[which(mut_aa %in% stop_aliases)])))

fwrite(mapped_counts, out_file, quote = F, sep = "\t", na = "NA")
