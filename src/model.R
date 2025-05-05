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
parser$add_argument("-n", "--nworkers", type = "numeric", default = 35,
    help = "Number of workers to use for model fitting", metavar = "nworkers")

# Argument parsing and I/O
args <- parser$parse_args()
mapped_counts_file <- args$file
nworkers <- args$nworkers

outfile <- args$outfile
model_output_path <- str_replace(args$outfile, ".sumstats.tsv", "_model_objects/")
marginals_outfile <- str_replace(args$outfile, ".sumstats.tsv", ".marginals.tsv")

dir.create(model_output_path)

mapped_counts <- read_tsv(mapped_counts_file) %>%
    mutate(condition_conc = as.factor(condition_conc),
           condition = relevel(as.factor(condition), ref = "DMSO_0"),
           mut_aa = relevel(as.factor(mut_aa), ref = "WT"))

# Set up formula and run model
form <- generate_formula(model_type)
nested_sumstats <- rand_effect_wrap(mapped_counts, form, "chunk", model_output_path, FALSE, nworkers)

# Extract and format coefficients and marginals
nested_coef <- nested_sumstats %>%
    select(-marginals) %>%
    unnest(coefs) %>%
    select(-c("effect", "component", "group", "dispersion")) %>%
    filter(str_detect(term, "mut")) %>%
    separate(term, into = c("condition", "aa"), sep = ":") %>%
    mutate(condition = str_remove(condition, "condition"), aa = str_remove(aa, "mut_aa"))

nested_marginals <- nested_sumstats %>%
    select(-coefs) %>%
    unnest(marginals) %>%
    select(-c("df", "statistic", "p.value")) %>%
    rename("aa" = "mut_aa")

# Write summary statistics
write_tsv(nested_coef, outfile)
write_tsv(nested_marginals, marginals_outfile)
