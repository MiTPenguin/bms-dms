## OCNT-VAMPLIB-1 Mini Assay Run #1: Chunks 3, 10, 14, 15

For the first mini-VAMP-seq experiment, we have 16 samples with the following structure:

| ID | Bins | Barcodes Per Variant | Cells Per Barcode |
| --- | ----------- | --- | --- |
| A | 25, 50, 75, 100 | 15 | 15 | 
| B | 25, 50, 75, 100 |15 | 15 | 
| C | 25, 50, 75, 100 | 30 | 15 |
| D | 25, 50, 75, 100 | 15 | 30 |

A and B are replicates of each other, while C and D are one replicate each of two different experimental procedures. All samples use the same library, OCNT-VAMPLIB-1, and there are not yet distinct "conditions" such as drug treatments to consider. Each of four samples (A - D) was sorted into four approximately equal bins, typically denoted as 25 through 100.

1. [Barcode Sequencing Distributions](#part1)
2. [Inference and Stop Codon Effects](#part2)
3. [Visualizations](#part3)
4. [Condition Comparions](#part4)
5. [Effect of Barcode Map Reduction](#part5)

### Barcode Sequencing Distributions <a name="part1"></a>

The plots below show the number of unique barcodes per variant in each sample on a linear (left) and log10 (right) scale. The red lines indicate 15 and 30 barcodes per variant, and samples are either separated by bin (top) or shown in aggregate (bottom). To compute the aggregate distributions, we used the union of all barcodes detected across all four bins for the given sample.


    
![png](OCNT-VAMPLIB-1-assay-run1_files/OCNT-VAMPLIB-1-assay-run1_4_0.png)
    



```R
bc_counts_aa %>% ungroup() %>% group_by(sample, chunk) %>% summarize(median(n))
```

    [1m[22m`summarise()` has grouped output by 'sample'. You can override using the
    `.groups` argument.



<table class="dataframe">
<caption>A grouped_df: 12 × 3</caption>
<thead>
	<tr><th scope=col>sample</th><th scope=col>chunk</th><th scope=col>median(n)</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>A</td><td> 3</td><td>17.0</td></tr>
	<tr><td>A</td><td>10</td><td>36.0</td></tr>
	<tr><td>A</td><td>14</td><td>43.0</td></tr>
	<tr><td>B</td><td> 3</td><td>15.0</td></tr>
	<tr><td>B</td><td>10</td><td>22.0</td></tr>
	<tr><td>B</td><td>14</td><td>25.0</td></tr>
	<tr><td>C</td><td> 3</td><td>29.0</td></tr>
	<tr><td>C</td><td>10</td><td>98.0</td></tr>
	<tr><td>C</td><td>14</td><td>75.0</td></tr>
	<tr><td>D</td><td> 3</td><td>29.0</td></tr>
	<tr><td>D</td><td>10</td><td>58.0</td></tr>
	<tr><td>D</td><td>14</td><td>45.5</td></tr>
</tbody>
</table>




```R
mapped_counts %>%
```


    
![png](OCNT-VAMPLIB-1-assay-run1_files/OCNT-VAMPLIB-1-assay-run1_7_0.png)
    


Counting the number of unique barcodes per chunk in each sample shows the expected depletion of chunk 15 relative to the others, and shows that samples C and D generally have more quantified barcodes as expected (with the exception of sample A, bin 100).


    
![png](OCNT-VAMPLIB-1-assay-run1_files/OCNT-VAMPLIB-1-assay-run1_9_0.png)
    


As before, we can plot the individual amino acid coverage per position, generating a large facet plot for each sample. One of these is shown below, and the rest are located [here](./coverage-plots).


    
![png](OCNT-VAMPLIB-1-assay-run1_files/OCNT-VAMPLIB-1-assay-run1_13_0.png)
    


### Inference and Stop Codon Effects <a name="part2"></a>

We can quantify variant effects in the VAMP-seq context in two ways. First, we can consider each sample and bin (e.g. C75) individually, and estimate the mutant vs wild-type effect. Alternatively, we can consider all four bins for a sample jointly, and model the bin itself as a numerical predictor (e.g. 0/1/2/3).

In the first model, we obtain four summary statistics for each variant (one per bin). In the second model, we obtain only one summary statistic for each variant (the slope of the mutant vs WT effect _across_ bins). Most of the time, the "slope" model is the important result we care about, but the "per-bin" model is useful as well to visualize particular positions in high detail.

To get started we extract and plot Stop effects for each sample group, considering A and B jointly as replicates:


    
![png](OCNT-VAMPLIB-1-assay-run1_files/OCNT-VAMPLIB-1-assay-run1_17_0.png)
    


We can also examine directly the number of variants across all comparisons that are significant at a 1% FDR:


    
    
    |source | GoF (1% FDR)| LoF (1% FDR)| Non-Significant|
    |:------|------------:|------------:|---------------:|
    |10-AB  |            6|          492|             902|
    |10-C   |           13|          503|             884|
    |10-D   |           14|          548|             838|
    |14-AB  |           13|          458|             928|
    |14-C   |            8|          418|             973|
    |14-D   |           11|          500|             888|
    |3-AB   |           15|          134|            1251|
    |3-C    |            8|          133|            1259|
    |3-D    |            7|          178|            1215|


### Visualizations <a name="part3"></a>

As before, we can show the same Log2FoldChanges and Z-statistics as heatmaps:


    
![png](OCNT-VAMPLIB-1-assay-run1_files/OCNT-VAMPLIB-1-assay-run1_22_0.png)
    



    
![png](OCNT-VAMPLIB-1-assay-run1_files/OCNT-VAMPLIB-1-assay-run1_23_0.png)
    


To get a better intuition for what "GoF" and "LoF" mean in the context of slopes across bins, we can extract the most significant GoF and LoF variants across all the plotted variants above, and display their _per-bin_ mutant vs WT summary statistics. First, let's get the top 5 GoF and LoF:


    
    
    |source | pos|chunk |aa |  estimate| std.error| statistic|  p.value|    p.adj|
    |:------|---:|:-----|:--|---------:|---------:|---------:|--------:|--------:|
    |3-AB   | 167|3     |D  | 0.7564435| 0.1596624|  4.737769| 2.20e-06| 1.34e-05|
    |10-C   | 696|10    |L  | 0.5261324| 0.1117830|  4.706729| 2.50e-06| 1.55e-05|
    |3-D    | 192|3     |S  | 0.7359611| 0.1664602|  4.421244| 9.80e-06| 5.56e-05|
    |14-AB  | 939|14    |S  | 0.4734079| 0.1086648|  4.356589| 1.32e-05| 7.33e-05|
    |14-AB  | 932|14    |W  | 0.5186452| 0.1194197|  4.343044| 1.41e-05| 7.76e-05|



    
    
    |source | pos|chunk |aa |   estimate| std.error| statistic| p.value| p.adj|
    |:------|---:|:-----|:--|----------:|---------:|---------:|-------:|-----:|
    |10-D   | 669|10    |P  | -1.1952460| 0.0427962| -27.92879|       0|     0|
    |10-AB  | 669|10    |P  | -0.9647506| 0.0382167| -25.24423|       0|     0|
    |10-C   | 669|10    |P  | -0.7186952| 0.0321140| -22.37947|       0|     0|
    |10-D   | 684|10    |S  | -0.9383830| 0.0449853| -20.85976|       0|     0|
    |3-D    | 154|3     |*  | -1.1078474| 0.0610024| -18.16072|       0|     0|


We can plot the strongest GoF and LoF variant in the above tables, across all samples and bins. For GoF variants this slope should be positive, while for LoF variants it should be negative:


    
![png](OCNT-VAMPLIB-1-assay-run1_files/OCNT-VAMPLIB-1-assay-run1_28_0.png)
    


### Condition Comparisons <a name="part4"></a>

Considering the slope analysis, we can evaluate the correlation between scores from each condition set (AB, C, and D):


    
![png](OCNT-VAMPLIB-1-assay-run1_files/OCNT-VAMPLIB-1-assay-run1_32_0.png)
    


Stops are less correlated since they are more consistently negative and thus have a smaller range. We can examine the significance boundaries for each of the three pairs more closely:


    
![png](OCNT-VAMPLIB-1-assay-run1_files/OCNT-VAMPLIB-1-assay-run1_35_0.png)
    



    
![png](OCNT-VAMPLIB-1-assay-run1_files/OCNT-VAMPLIB-1-assay-run1_35_1.png)
    



    
![png](OCNT-VAMPLIB-1-assay-run1_files/OCNT-VAMPLIB-1-assay-run1_35_2.png)
    


### Effect of Barcode Map Reduction <a name="part5"></a>

If we separate the barcode map into A and B subcomponents, but leave the rest of the pipeline the same, we get the following numbers of entries in the barcode maps for each chunk:


    
    
    | chunk|   total|       A|       B|
    |-----:|-------:|-------:|-------:|
    |    10| 1509803|  647418|  862385|
    |    14| 2288729| 1006554| 1282175|
    |     3| 1073748|  360452|  713296|


To determine the effects of only integrating A or B but not both, we generated the barcode map using either only A or only B but with all other filters (e.g. across sequencing replicates) intact. From here, we can extract the median number of barcodes per variant under each map and in each condition:


    
    
    |sample | 3A| 3B| 3AB| 10A| 10B| 10AB| 14A| 14B| 14AB|
    |:------|--:|--:|---:|---:|---:|----:|---:|---:|----:|
    |A      |  9|  8|  17|   7|  29|   36|   9|  35|   43|
    |B      |  6|  8|  15|   7|  15|   22|   9|  16|   25|
    |C      | 19| 10|  29|  10|  88|   98|  11|  64|   75|
    |D      | 14| 14|  29|  13|  45|   58|  15|  31|   46|


The above table indicates that, for many samples, the number of barcodes detected in a sample from a given A or B pool is more imbalanced than would be expected given the contribution of each pool to the barcode map. To evaluate the extent to which this propagates to the final summary statistics, we run the model using either A or B separately, and compare it to the original model where we considered the union of A and B. The most useful way to summarize the result is via the standard errors:


    
![png](OCNT-VAMPLIB-1-assay-run1_files/OCNT-VAMPLIB-1-assay-run1_46_0.png)
    


In all cases, both A and B have higher standard errors than Full due to having fewer barcodes. However, the standard error increase is generally much larger when only pool A is used. This is consistent with (generally) a much larger number of detected barcodes originating from pool B compared to pool A. This is true despite the fact that both pool A and pool B contained similar numbers of barcodes during barcode mapping. Thus, taken together, this suggests some imbalance when pooling A and B for the assay.
