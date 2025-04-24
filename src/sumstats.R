library(argparse)
library(tidyverse)

parser <- ArgumentParser()
parser$add_argument("-s", "--sumstats", type = "character",
    help = "Summary statistics file", metavar = "sumstats")
parser$add_argument("-m", "--marginals", type = "character",
    help = "Marginals file", metavar = "marginals")

# Argument parsing and I/O
args <- parser$parse_args()
sumstats <- read_tsv(args$sumstats)
marginals <- read_tsv(args$marginals)
sumstats_out <- str_replace(args$sumstats, ".sumstats.tsv", "-cleaned.sumstats.tsv")
marginals_out <- str_replace(args$marginals, ".sumstats.marginals.tsv", "-cleaned.sumstats.marginals.tsv")

df_sumstats_clean <- sumstats %>%
    select(-c("effect", "component", "group", "dispersion")) %>%
    filter(str_detect(term, "mut")) %>%
    separate(term, into = c("condition", "aa"), sep = ":") %>%
    mutate(condition = str_remove(condition, "condition"), aa = str_remove(aa, "mut_aa"))

df_marginals_clean <- marginals %>%
    select(-df, -statistic, -p.value) %>%
    rename("aa" = "mut_aa")

write_tsv(df_sumstats_clean, sumstats_out)
write_tsv(df_marginals_clean, marginals_out)
