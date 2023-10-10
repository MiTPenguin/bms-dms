## TYK2 VAMP-seq Midpoint Inference

For the first mini-VAMP-seq experiment, we have four activity bins ordered low to high and with approximately (but not exactly) equal widths. This document walks through a prototype modeling procedure for these data using position 669 as a test case. This is a useful position since it contains a spike-in (669P) as well as the usual stop effect. We include all non-WT variants at this position, as well as a sampling (1000 barcodes per sample) of the WT barcodes.

1. [Background and Descriptive Profiling](#part1)
2. [Adapting the DMS Model to VAMP-seq](#part2)
3. [Rescaling Coefficients to Weights](#part3)
4. [Inferring Differences from WT](#part4)
5. [Extension to Chunk 10](#part5)

### Background and Descriptive Profiling <a name="part1"></a>

Let's start by looking at what we got using the original DMS model with VAMP-seq data. For this model, we estimated the mutant vs WT effects within each bin individually:


    
![png](TYK2-VAMPseq-Midpoint_files/TYK2-VAMPseq-Midpoint_4_0.png)
    


The shared pattern is increased abundance in low-activity bins, and reduced abundance in high-activity bins, relative to WT. This is why the linear slope across these intervals do a decent, but non-optimal, job of capturing GoF or LoF patterns. To examine the underlying model in more detail, let's regenerate it and extract the WT and mutant averages separately:


    
![png](TYK2-VAMPseq-Midpoint_files/TYK2-VAMPseq-Midpoint_8_0.png)
    


This is unusual: it looks like bin 50 is markedly higher than the others, particularly for WT. We can examine the underlying raw counts and count the total read counts, the number of unique barcodes, and the actual distributions directly:


    
    
    |sample | total read count| unique barcode count|
    |:------|----------------:|--------------------:|
    |C25    |         14.98547|                69178|
    |C50    |         15.38236|                29886|
    |C75    |         15.48106|               104692|
    |C100   |         15.41413|                67064|



    
![png](TYK2-VAMPseq-Midpoint_files/TYK2-VAMPseq-Midpoint_11_0.png)
    


Indeed, it does look like there are fewer unique barcodes in bin 50 and they generally have higher  read counts compared to the other three bins. This poses problems even for a simple ratio approach - if we took each barcode and computed the proportion of reads from a given barcode in that bin, many barcodes would appear "elevated" in bin 50. Let's see what we can do about that.

### Adapting the DMS Model to VAMP-seq <a name="part2"></a>

Recall that the main DMS model is specified as:

$$count \sim bin + bin:AA + (1|barcode) + offset(stopCount)$$

The issues described above are "normalized out" by comparison with WT, since WT provides the "within-bin" reference against which we compute variant effects. However, we would like (for many reasons) to be able to estimate WT indiviually as well, despite bin 50. One plausible solution could be the _offset_. In the main model, the offset is total stop counts, which is not appropriate given that stops will be very unevenly distributed across bins. Instead, the true offset here should be the _average_ read count instead of the general total or stop total. This results in a new model:

$$count \sim bin + bin:AA + (1|barcode) + offset(averageLogCount)$$

where `averageLogCount` is defined as the `mean(log(count))` for all counts from a given sample. If we run this model and extract the same marginals, with no other changes, we obtain:

    NOTE: A nesting structure was detected in the fitted model:
        mut_aa %in% condition_conc
    



    
![png](TYK2-VAMPseq-Midpoint_files/TYK2-VAMPseq-Midpoint_14_1.png)
    


This looks much better! WT is more clearly elevated to the right, and the two LoF variants are elevated to the left with the proline more negative than the stop (which is plausible given the structure of this assay). The wide standard errors are to be expected given that this is only sample C, and the generally weaker effect of stops here. We can extract and plot the mutant vs WT effects to generate a plot identical to that from the original model but with this change:


    
![png](TYK2-VAMPseq-Midpoint_files/TYK2-VAMPseq-Midpoint_16_0.png)
    


We generally see the same pattern, but shifted down a bit and with slightly smaller effect sizes _and_ standard errors (so similar significance results). However, the variant effects _across bins_ are much more left-shifted and (from the previous plot) have a reasonable WT profile.

### Rescaling Coefficients to Weights <a name="part3"></a>

Now, all this modeling so far relates a quantity shown on the y-axis (count/abundance or log2FoldChange) to the sample IDs on the x-axis. However, if we consider the x-axis bins as representing quantitative activity levels, we would like to compute a "midpoint" on a 0-1 scale of all activity from a given variant across bins.

To do this, we define a `score` computed as follows:

$$score = b_1w_1 + b_2w_2 + b_3w_3 + b_4w_4$$

where bins are ordered 1 through 4 in order of increasing activity. Here, $b_i$ indicates the activity level of the $ith$ bin, and $w_i$ indicates the `weight` of that bin. For our purposes, we set the $b_i$ equal to the midpoint of the bin, resulting in values of: `0.125`, `0.375`, `0.625`, and `0.875`.

However, that still leaves the question of determining appropriate values for the $w_i$. In the previous section, we noted that we can view a model in terms of the abundances of mutant and WT separately, or their log fold change. Let's consider one variant (669P) along with WT and how we might use regression summary statistics to compute weights:


    
![png](TYK2-VAMPseq-Midpoint_files/TYK2-VAMPseq-Midpoint_19_0.png)
    


Next, we want to only consider _relative_ changes in mean count rather than the absolute change - if a variant shifts all-up or all-down relative to WT but with the same pattern, that would generally reflect a compositional issue rather than a biological effect. Thus, we _normalize_ these values by i) subtracting the smallest value and ii) dividing by the sum of the remaining values. 


    
![png](TYK2-VAMPseq-Midpoint_files/TYK2-VAMPseq-Midpoint_21_0.png)
    


This is great! Now we have weights we can use to compute the score for each mutant and WT:


<table class="dataframe">
<caption>A tibble: 2 × 6</caption>
<thead>
	<tr><th scope=col>mut_aa</th><th scope=col>25</th><th scope=col>50</th><th scope=col>75</th><th scope=col>100</th><th scope=col>score</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>WT</td><td>0.1950280</td><td>0</td><td>0.3677614</td><td>0.43721054</td><td>0.649125</td></tr>
	<tr><td>P </td><td>0.7220757</td><td>0</td><td>0.2038817</td><td>0.07404262</td><td>0.281125</td></tr>
</tbody>
</table>



Now, we would like to test for a difference between mutant and WT. To do this, we need some estimate of the error on these scores, which is very tricky given the post-processing of the regression coefficients we had to do.

### Inferring Differences from WT <a name="part4"></a>

Recall that the model provieds us a set of normalized log abundance estimates with standard errors for each bin and each variant, including WT. These errors are the only source of error or variance at all in the final score, but propagating the error from the regression outputs to that score is tricky and non-intuitive. An alternative is based in sampling - our regression outputs are Z-distributed (and tested with the Z-statistic), and so we can sample from them using the appropriate normal distribution.

To determine how variability in our abundance estimates propagates to our weights and final score, we can sample many sets of abundances to see how random changes (scaled by their standard errors) influence the final score. Let's walk through step-by-step:

First, we have the abundance estimates from the model:


    
![png](TYK2-VAMPseq-Midpoint_files/TYK2-VAMPseq-Midpoint_27_0.png)
    


Next, we take the means and standard errors from this model and sample from them many times - below, we shown an example with a handful of points to illustrate:


    
![png](TYK2-VAMPseq-Midpoint_files/TYK2-VAMPseq-Midpoint_33_0.png)
    


Then, we can compute the score for each simulated set of points after following the same procedure as before:

- Subtract the smallest value
- Divide by the sum of the remaining three values to obtain the weights $w_i$
- Compute the score as $\sum_iw_is_i$ for $s_i \in {0.125, 0.375, 0.625, 0.875}$

If we do this many more times (here, 10000) we obtain the following distributions:

    [1m[22m`summarise()` has grouped output by 'n'. You can override using the `.groups`
    argument.



    
![png](TYK2-VAMPseq-Midpoint_files/TYK2-VAMPseq-Midpoint_35_1.png)
    


We can do hypothesis testing and other types of inference on these distributions, but the main quantity we are interested in is their mean and standard deviation, the latter corresponding to the standard error of our estimated quantity. Here, we can extract the mean and SD of these distributions and test for a difference between P and WT:


    
    
    | estimate_P| estimate_WT| std.error_P| std.error_WT|   estimate| std.error| statistic|
    |----------:|-----------:|-----------:|------------:|----------:|---------:|---------:|
    |  0.2806714|   0.6372446|   0.0367323|    0.0110475| -0.3565732| 0.0383576| -9.296014|


### Extension to Chunk 10 <a name="part5"></a>

We set the bin values to their midpoints of 0.125, 0.375, 0.625, and 0.875. Then, we compute the weights and midpoints by resampling from the marginals extracted from the negative binomial model. Notably, we get a value of WT for each position (and all these values are extremely similar but not identical) so we specify the WT value as the median mean and median standard error across all positions. Doing so returns the following WT estimates:


    
    
    |  WT score| WT score standard error|
    |---------:|-----------------------:|
    | 0.6756111|               0.0058894|


Similarly, we can take the mean and standard deviation of the scores computed for each variant across all simulations, and plot either the subtractive difference from WT (the contrast effect size) or the Z-statistic of the difference. Using the Z-statistic results in a plot that is quite similar to those for DMS, including the pronounced influence of the spike-ins:


    
![png](TYK2-VAMPseq-Midpoint_files/TYK2-VAMPseq-Midpoint_46_0.png)
    



    
![png](TYK2-VAMPseq-Midpoint_files/TYK2-VAMPseq-Midpoint_47_0.png)
    



    
![png](TYK2-VAMPseq-Midpoint_files/TYK2-VAMPseq-Midpoint_47_1.png)
    

