#!/usr/bin/env bash

find . -name "*.fastq.gz" | \
    cut -d . -f 1 | \
    parallel "./count-barcodes.sh {}.fastq.gz {}.rna-bcs.tsv"

Rscript --vanilla fit-preprocess.R \
    ../mapped_counts/TYK2-run1.mapped-counts.tsv \
    ../barcode_maps/tyk2.bcmap-final.tsv \
    run1/sample-properties-rna.tsv \
    run1 &

Rscript --vanilla fit-preprocess.R \
    ../mapped_counts/TYK2-run2.mapped-counts.tsv \
    ../barcode_maps/tyk2.bcmap-final.tsv \
    run2/sample-properties-rna.tsv \
    run2 &

Rscript --vanilla fit-preprocess.R \
    ../mapped_counts/TYK2-run3.mapped-counts.tsv \
    ../barcode_maps/tyk2.bcmap-final.tsv \
    run3/sample-properties-rna.tsv \
    run3 &

Rscript --vanilla fit-preprocess.R \
    ../mapped_counts/TYK2-vamp2.mapped-counts.tsv \
    ../barcode_maps/tyk2-vamp.bcmap-final.tsv \
    vamp2/sample-properties-rna.tsv \
    vamp2 &
