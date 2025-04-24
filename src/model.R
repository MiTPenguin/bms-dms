#!/usr/bin/env Rscript
# Negative Binomial Mixed Modeling of DMS Variant Effects

library(argparse)
library(tidyverse)
source("/bms-dms/src/model_utils.R")

RNGkind("L'Ecuyer-CMRG")
set.seed(1)

# Command line arguments
parser <- ArgumentParser()
parser$add_argument("-f", "--file", type = "character",
    help = "Mapped Counts File", metavar = "file")
parser$add_argument("-o", "--outfile", type = "character",
    help = "Output File", metavar = "outfile")
parser$add_argument("-m", "--model", type = "character",
    help = "Model Type", metavar = "model")
parser$add_argument("-s", "--stops", type = "character",
    help = "Stop Handling: agg or nonagg", metavar = "stops")
parser$add_argument("-n", "--nworkers", type = "numeric", default = 35,
    help = "Number of workers to use for model fitting", metavar = "nworkers")

# Argument parsing and I/O
args <- parser$parse_args()
mapped_counts_file <- args$file
model_type <- args$model
stops <- args$stops
nworkers <- args$nworkers

outfile <- args$outfile
model_output_path <- str_replace(args$outfile, ".tsv", "_model_objects/")
marginals_outfile <- str_replace(args$outfile, ".tsv", ".marginals.tsv")

dir.create(model_output_path)

# Read VERSION and stop if absent
if (!file.exists("VERSION")) {
    stop("VERSION file not found")
} else {
    version <- readLines("VERSION")
}

# Determine offset scope and format mapped counts accordingly
mapped_counts <- read_tsv(mapped_counts_file) %>%
    mutate(condition_conc = as.factor(condition_conc),
           condition = relevel(as.factor(condition), ref = "DMSO_0"),
           mut_aa = relevel(as.factor(mut_aa), ref = "WT"))

# Aggregate stops to their own position if desired
if (stops == "agg") {
    mapped_counts <- mapped_counts %>%
        mutate(pos = if_else(mut_aa %in% stop_aliases, "-1", pos))
} else if (stops == "nonagg") {
    NULL
} else {
    stop("Invalid stop handling option")
}

# Set up formula and run model
nestvars <- "chunk"
form <- generate_formula(model_type)
nested_sumstats <- rand_effect_wrap(mapped_counts, form, nestvars, model_output_path, FALSE, nworkers)

# Extract and format coefficients and marginals
nested_coef <- nested_sumstats %>% select(-marginals) %>% unnest(coefs)
nested_marginals <- nested_sumstats %>% select(-coefs) %>% unnest(marginals)

# Add VERSION
nested_coef$version <- version

# Write summary statistics
write_tsv(nested_coef, outfile)
write_tsv(nested_marginals, marginals_outfile)
