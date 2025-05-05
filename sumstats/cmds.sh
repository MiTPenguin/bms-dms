#!/usr/bin/env bash

Rscript --vanilla src/model.R \
    --file mapped_counts/TYK2-run1.mapped-counts.tsv \
    --outfile sumstats/TYK2-run1.sumstats.tsv \
    --nworkers 35 \
    --model global

Rscript --vanilla src/model.R \
    --file mapped_counts/TYK2-run2.mapped-counts.tsv \
    --outfile sumstats/TYK2-run2.sumstats.tsv \
    --nworkers 35 \
    --model global

Rscript --vanilla src/model.R \
    --file mapped_counts/TYK2-run3.mapped-counts.tsv \
    --outfile sumstats/TYK2-run3.sumstats.tsv \
    --nworkers 35 \
    --model global

Rscript --vanilla src/model.R \
    --file mapped_counts/TYK2-vamp2.mapped-counts.tsv\
    --outfile sumstats/TYK2-vamp2.sumstats.tsv \
    --nworkers 35 \
    --model flow
