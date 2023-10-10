### TYK2 Assay Comparisons Across Cell Counts

With the first mini-DMS dataset, we identified candidate partial loss-of-function variants [here](https://github.com/octantbio/dms/blob/master/notebooks/TYK2-Partial-LoF.md) based on whether the change in mutation effect between IFNa and None was both significantly different from zero, and significantly different from the aggregated stop effect. Here, we extend and generalize that analysis to incorporate the second mini-DMS assay dataset. This dataset includes an additional condition for evaluation as a compositional normalized, IFNb, and three separate cell counts. For this analysis, we have the following relevant comparisons:

| Assay | Cell Count | Comparison | Replicate Count |
| :--- | :--- | :--- | :--- |
| #1 | 14 million | IFNa vs None | 3 vs 3 |
| #2 | 7 million | IFNa vs None | 3 vs 3 |
| #2 | 3 million | IFNa vs None | 3 vs 3 |
| #2 | 1 million | IFNa vs None | 3 vs 3 |

We also want to consider the effects of using IFNb as a normalization control, pretending that we have _only_ IFNb and "forgetting" about None. This happens because the way our model shares barcode information means that data points from the None condition can influence the results from comparing IFNa vs IFNb, in particular for barcodes that None shares with IFNa but not IFNb (or vice versa). So, we also have the following comparisons, computed without using None barcode information in any way:

| Assay | Cell Count | Comparison | Replicate Count |
| :--- | :--- | :--- | :--- |
| #2 | 7 million | IFNa vs IFNb | 3 vs 3 |
| #2 | 3 million | IFNa vs IFNb | 3 vs 3 |
| #2 | 1 million | IFNa vs IFNb | 3 vs 3 |

For each of these comparisons, we have a set of summary statistics. In each set, each mutation has an estimated effect size, denoted here as the log2 fold change of IFNa vs None or IFNa vs IFNb, and an associated standard error. The ratio of the effect size to the standard error defins the z-statistic, which is used to compute nominal and Benjamini-Hochberg adjusted p-values. In general, we show error bars as +/- 3 standard errors and define FDR cutoffs at 1%.

#### Distributions of Mutation Effects

In general, higher cell densities result in data with higher numbers of successfully quantified barcodes, leading to a decrease in standard error. We also observe a similar effect with fold changes, where lower cell densities have slightly exaggerated effect sizes along with elevated variance.


    
![png](TYK2-Assay-Comparisons_files/TYK2-Assay-Comparisons_5_0.png)
    


For each of these comparisons, we can plot the long-form confidence intervals ordered by log2 fold change. While the distributions shown in the two columns below are the same, the left side highlights stop mutations in red, while the right side highlights alanine mutations: 


    
![png](TYK2-Assay-Comparisons_files/TYK2-Assay-Comparisons_8_0.png)
    


If we aggregate the stop effects into a single value, we can compute a confidence interval on the global average stop effect. This is _a priori_ reaosnable, since we expect a stop codon to have a similar effect on protein function regardless of where in chunk 10 it is located.


    
![png](TYK2-Assay-Comparisons_files/TYK2-Assay-Comparisons_11_0.png)
    


For each comparison, we can visually observe candidate partial loss-of-function variants as those whose gray confidence intervals do not overlap either the black or red horizontal lines. These are more obvious in comparisons derived from data with larger cell counts.

#### Comparisons Across Assay Parameters

One straightforward way to see whether functional variant effects (whether LoF or partial LoF) are reproducible across assay conditions is to view a scatterplot with all pairwise comparisons. This shows that the correlation is generally very clear for all pairs, but the benefit of higher cell counts is clear. The most strongly correlated effects are observed between 14 million and 7 million cells, normalized to Untreated:


    
![png](TYK2-Assay-Comparisons_files/TYK2-Assay-Comparisons_15_0.png)
    


We can more compactly visualize the correlation matrix of the effect sizes (i.e. Log2 Fold Changes). This summarizes and quantifies the relationships already visible above:


    
![png](TYK2-Assay-Comparisons_files/TYK2-Assay-Comparisons_17_0.png)
    


This is all encouraging, and we know from the first section that the stop codons are all highly significant and highly negative in all comparisons. The next step is to identify partial loss-of-function variants in each condition. To do this, we consider the Z-statistic (defined as the effect size or log2 fold change divided by the standard error) of each variant relative to WT, and also relative to the combined Stop effect. The plots below show this score for the comparison to WT on the x-axis, and for the comparison to the combined Stop on the y-axis:


    
![png](TYK2-Assay-Comparisons_files/TYK2-Assay-Comparisons_20_0.png)
    


The red points indicate variants that are significant at a 1% FDR against both WT and the combined Stop effect. Visually, this is intuitive, as those are the points that have the most negative x-axis score (most negative relative to WT) while simulataneously having a high y-axis score (most positive relative to the combined Stop effect).

Using the partial loss-of-function variant sets defined above, we can create a higher-resolution comparison scatterplot between pairs of comparisons, now highlighting variants that are significantly partial LoF in one condition, the other, or both. These plots are generally too dense to show in a grid, so a directory of them (one file per pairwise comparison) is located [here](TYK2-Assay-Comparisons-pairwise-scatterplots), and two examples are shown below that describe the extremes: 14 million using untreated vs 7 million using untreated is the most reproducible, while 1 million using untreated vs 3 million using IFNb has the poorest correlation:


    
![png](TYK2-Assay-Comparisons_files/TYK2-Assay-Comparisons_25_0.png)
    

