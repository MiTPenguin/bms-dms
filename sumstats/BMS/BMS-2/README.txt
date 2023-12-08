TYK2 DMS and FlowDMS Rescaled Scores

There are two TSV files, one each for DMS and FlowDMS summary statistics which have been re-scaled into 0-1 scores. Each row contains summary statistics for a single variant, and can be joined on the first three columns noted below.

For DMS, we use the normalized summary statistics from run 7, IFN-alpha 100 U/mL minus Untreated. For FlowDMS, we use the most recent midpoints computed from run 2.For both datasets, a score of zero corresponds to the average stop, and a score of one corresponds to wild-type.

NOTE: these scores are not bounded by 0 and 1, but rather those are the reference points for Stop and Wild-Type, respectively.

The columns in DMS_IFNalpha100vsUntreated_rescaled.tsv are:
    - Chunk: TYK2 library segment
    - Position: TYK2 amino acid position
    - AA: Amino acid
    - Log2FoldChange: Estimated effect size, log2 scale
    - Log2FoldChangeError: Standard error of the Log2FoldChange
    - FDR: Benjamini-Hochberg-adjusted p-values
    - ScaledScore: 0-1 rescaled score
    - ScaledScore_Lower: Lower 95% confidence boundary for the ScaledScore
    - ScaledScore_Upper: Upper 95% confidence boundary for the ScaledScore

The columns in FlowDMS_midpoints_rescaled.tsv.tsv are:
    - Chunk: TYK2 library segment
    - Position: TYK2 amino acid position
    - AA: Amino acid
    - Midpoint: Estimated sorting bin midpoint
    - Midpoint_Lower: Lower 95% confidence boundary for the Midpoint
    - Midpoint_Upper: Upper 95% confidence boundary for the Midpoint
    - ScaledScore: 0-1 rescaled score
    - ScaledScore_Lower: Lower 95% confidence boundary for the ScaledScore
    - ScaledScore_Upper: Upper 95% confidence boundary for the ScaledScore
