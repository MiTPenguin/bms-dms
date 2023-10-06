TYK2 Deep Mutational Scanning Summary Statistics

There are two TSV files representing two comparisons, each with the same number of rows (23647, ignoring the header) and columns (9). Each row contains summary statistics for a single variant, and quantifies the change in the Mutant vs Wild-Type activity difference between the indicated two conditions. In other words, it is a difference of differences: the gap between Mutant and Wild-Type activity is quantified in each condition, and then the effects here represent the change in that gap between conditions.

The columns are:
    - Position: TYK2 amino acid position, 1-1187
    - Chunk: TYK2 library segment, 1-17
    - AA: Amino acid, 20 expected per position (19 non-WT plus 1 stop)
    - Log2FoldChange: Estimated effect size (see above)
    - StandardError: Standard error of the Log2FoldChange
    - Statistic: Ratio of the Log2FoldChange to the StandardError
    - Pvalue: Un-adjusted p-value against the hypothesis that the Log2FoldChange is zero
    - Padj: Benjamini-Hochberg-adjusted p-values
    - Comparison: Conditions being compared at each variant

The directories titled vN are previous versions of the summary statistics sent to BMS. The current version is v2.