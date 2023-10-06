## OCNT-DMSLIB-1 Assay, Chunk 10

This report documents the processing of the first complete DMS assay covering chunk 10 of TYK2. The new sequencing run with expressed barcodes is `221110_VH00964_9_AAAW57KHV` which includes 3 replicates each of 4 sample groups, for a total of 12 samples. The four sample groups correspond to two clones (c1 and c8) in each of two conditions (0 or 100 U/mL cytokine treatment). The barcode map is from sequencing run `220824_VH00964_6_AAAMJF3HV` and is the combination of 10C and 10D after removing any barcodes detected in both (experimentally, they should be non-overlapping).

1. [Barcode Sequencing Quality Control](#part1)
2. [Statistical Models](#part2)
3. [Summary Statistics](#part3)
    - [Within Condition](#part3a)
    - [Between Conditions](#part3b)
4. [Visualization and Interpretation](#part4)

### Barcode Sequencing Quality Control <a name="part1"></a>

To associate barcodes with oligos, we need to use a barcode map sequenced in a previous 2x150 run. We use the combination of the maps from `10C` and `10D`, after subtracting any barcodes that are detected in both. This returns approximately 5.5 million unique, oligo-barcode pairs that map to perfect oligos.

The expressed barcode 1x26 samples are sequenced to 65-80 million reads per sample. In all samples, the large majority of barcodes were previously observed in a barcode mapping run, though frequently not associated with a perfect oligo (as expected). In general, approximately 60-65% of the raw reads in each sample were matched to a perfect, usable oligo. Here, "subthreshold" means that the barcode failed one of the filters applied during barcode mapping (>3 reads with >75% purity in both replicates). These categories are **non-overlapping** and account for all observed barcodes in each sample:


```R
mapped_counts <- fread("../pipeline-tmp/DMSRUN1.mapped-counts.tsv") %>%
    separate(oligo, c("lib", "chunk", "wt_aa", "pos",
        "mut_aa", "wt_codon", "mut_codon"), "_") %>%
    group_by(sample) %>%
    mutate(total_counts = sum(count[which(mut_aa == "*")]) / 1000000)

mapped_counts$mut_aa[which(mapped_counts$wt_aa == mapped_counts$mut_aa | is.na(mapped_counts$mut_aa))] <- "WT" #nolint
mapped_counts$mut_aa <- relevel(as.factor(mapped_counts$mut_aa), ref = "WT")
mapped_counts$condition <- relevel(as.factor(mapped_counts$condition), ref = "none")
```


    
    
    |sample | total depth|  perfect| perfect, subthreshold| imperfect, previously observed| not previously observed|
    |:------|-----------:|--------:|---------------------:|------------------------------:|-----------------------:|
    |1A     |    81260580| 49350272|              12302795|                       17073741|                 2533772|
    |1B     |    70139871| 42512705|              10643354|                       14790333|                 2193479|
    |1C     |    66728935| 40365852|               9997156|                       14104873|                 2261054|
    |2A     |    60768685| 38803335|               9527227|                       10582745|                 1855378|
    |2B     |    76426470| 48605224|              11978344|                       13287914|                 2554988|
    |2C     |    84482619| 53752398|              13234312|                       14822459|                 2673450|
    |3A     |    79167010| 48283274|              11939808|                       16368307|                 2575621|
    |3B     |    68929027| 42121217|              10326971|                       14269903|                 2210936|
    |3C     |    79821325| 48698703|              12096126|                       16445310|                 2581186|
    |4A     |    72174838| 45944997|              11361745|                       12489972|                 2378124|
    |4B     |    74707942| 47538600|              11785181|                       12944377|                 2439784|
    |4C     |    71630461| 45643352|              11313965|                       12365913|                 2307231|


Below, the same table expressed as percentages:


    
    
    |sample | total depth|  perfect| perfect, subthreshold| imperfect, previously observed| not previously observed|
    |:------|-----------:|--------:|---------------------:|------------------------------:|-----------------------:|
    |1A     |         100| 60.73089|              15.13993|                       21.01110|                3.118083|
    |1B     |         100| 60.61132|              15.17447|                       21.08691|                3.127293|
    |1C     |         100| 60.49228|              14.98174|                       21.13757|                3.388416|
    |2A     |         100| 63.85416|              15.67786|                       17.41480|                3.053181|
    |2B     |         100| 63.59737|              15.67303|                       17.38653|                3.343067|
    |2C     |         100| 63.62539|              15.66513|                       17.54498|                3.164497|
    |3A     |         100| 60.98913|              15.08180|                       20.67567|                3.253402|
    |3B     |         100| 61.10810|              14.98204|                       20.70231|                3.207554|
    |3C     |         100| 61.00964|              15.15400|                       20.60265|                3.233705|
    |4A     |         100| 63.65792|              15.74198|                       17.30516|                3.294949|
    |4B     |         100| 63.63259|              15.77500|                       17.32664|                3.265763|
    |4C     |         100| 63.72059|              15.79491|                       17.26348|                3.221019|


This is generally consistent with the previously estimated "sequencing + synthesis" error rate of 35%, along with a small fraction (~3.3%) of barcodes which are observed for the first time in the new data. These are probably sequencing errors which occurred in the barcode sequencing data, but were not observed in any prior run.

Now, we consider only those barcodes which are in the `perfect` category in the above tables. Each of these barcodes is attached to a particular oligo, so we can view the spread of how many barcodes are attached to each codon or each amino-acid level sequence:


    
![png](OCNT-DMSLIB-1-assay-run1_files/OCNT-DMSLIB-1-assay-run1_14_0.png)
    



    
![png](OCNT-DMSLIB-1-assay-run1_files/OCNT-DMSLIB-1-assay-run1_14_1.png)
    


Each of the above points is a particular codon variant or residue variant, and the count is the number of unique barcodes for that variant. We can break out the clone 8 distributions more directly as densities for presentation:


    
![png](OCNT-DMSLIB-1-assay-run1_files/OCNT-DMSLIB-1-assay-run1_16_0.png)
    



    
![png](OCNT-DMSLIB-1-assay-run1_files/OCNT-DMSLIB-1-assay-run1_17_0.png)
    


We can also plot the actual read counts of each barcode across samples. This shows an expected spike at zero for barcodes with exactly one read count:


    
![png](OCNT-DMSLIB-1-assay-run1_files/OCNT-DMSLIB-1-assay-run1_19_0.png)
    


Note that while these barcodes have only one read count in a given sample in the 1x26 dataset, the fact that it is in the barcode map means it was detected in at least two sequencing replicates with at least three reads each in the barcode mapping data. Nonetheless, we can verify that we still would have enough barcode coverage even if we removed "lowly expressed" barcodes - again, we will not do this for the actual analysis, but we can just check:


    
![png](OCNT-DMSLIB-1-assay-run1_files/OCNT-DMSLIB-1-assay-run1_22_0.png)
    



    
![png](OCNT-DMSLIB-1-assay-run1_files/OCNT-DMSLIB-1-assay-run1_22_1.png)
    


We still observe a median of >30 barcodes per codon level variant, even requiring that a barcode have at least four reads to be included for a given sample. We do not need to apply this filter given how our regression model works, but even if we did our depth would still be quite good.

We can also assess coverage by plotting the numbers of barcodes per variant across each amino acid or codon and position. Since the samples are fairly similar


    
![png](OCNT-DMSLIB-1-assay-run1_files/OCNT-DMSLIB-1-assay-run1_25_0.png)
    



    
![png](OCNT-DMSLIB-1-assay-run1_files/OCNT-DMSLIB-1-assay-run1_25_1.png)
    


We can break out combined and polished versions of these plots for external presentation with just clone 8:


    
![png](OCNT-DMSLIB-1-assay-run1_files/OCNT-DMSLIB-1-assay-run1_28_0.png)
    


### Statistical Models <a name="part2"></a>

We use a hierarchical negative binomial mixed model to quantify mutation vs wild-type (WT) effects, and to test whether they are different from zero or across conditions. We apply two slight variations on the same underlying model: one that shares barcode information across condition for the inference of _across condition_ effects, and another that limits barcode sharing to within condition for the inference of _within condition_ effects. The architecture of the underlying model is shown below, with the difference being how information is shared between the colored circles at the top:


    
![png](OCNT-DMSLIB-1-assay-run1_files/OCNT-DMSLIB-1-assay-run1_30_0.png)
    


In the separated model, we share information within each condition only. In the schematic below, the dashed lines mark the boundaries of information sharing:


    
![png](OCNT-DMSLIB-1-assay-run1_files/OCNT-DMSLIB-1-assay-run1_32_0.png)
    


In contrast, the separated model pools information across both conditions to maximise power for detecting within-barcode changes across condition:


    
![png](OCNT-DMSLIB-1-assay-run1_files/OCNT-DMSLIB-1-assay-run1_34_0.png)
    


We obtain three sets of summary statistics: the WT vs Mutant effects in each of Untreated and IFNa separately, as well as the difference between those two effects. We share barcode information between conditions only when the test is for a difference between conditions. Applying a Benjamini-Hochberg correction and a 1% FDR threshold, we obtain the following counts of significant mutations in each comparison:


    
    
    |comparison         |model    |clone | not significant| significant|
    |:------------------|:--------|:-----|---------------:|-----------:|
    |none               |separate |c1    |            1400|           0|
    |none               |separate |c8    |            1400|           0|
    |IFNa2a             |separate |c1    |            1110|         290|
    |IFNa2a             |separate |c8    |            1196|         204|
    |IFNa2a - Untreated |combined |c1    |            1090|         310|
    |IFNa2a - Untreated |combined |c8    |            1149|         251|


When Untreated, we observe no significant effects at a 1% FDR in either clone, while when treated with IFNa we observe either 290 or 204 mutations with significantly different activity from wild-type without sharing barcode information. If we allow the sharing of barcode information and test for a difference between conditions, the number of significant mutations increases to 310 and 251 in each clone, reflecting the expected increase in power. However, note that while more powerful this is technically a different comparison (e.g. evaluating whether two Mut vs WT effects are different from each other, rather than whether one Mut vs WT effect is different from zero).

A natural next question is whether the same or different variants are being detected between clones and comparisons, which we can visualize via an upset plot:


    
![png](OCNT-DMSLIB-1-assay-run1_files/OCNT-DMSLIB-1-assay-run1_42_0.png)
    


The large majority of mutations with a signiicant effect are significant in both IFNa considered alone, as well when testing the difference in effect from Untreated. However, a subset of mutations become significant only in the "IFNa - Untreated" comparison due to increased power, and a very small number become non-significant. In this last case, these are variants which display a similar significant difference from wild-type in both Untreated and IFNa.

To see the concrete improvement in power that comes with barcode information sharing, we can compare the `IFNa - Untreated` comparison above to the same comparison, but instead calculate it by subtracting the separately-estimated `IFNa` and `Untreated` values of using the combined model.

Below, we show the Z-statistic of the final (Mut vs WT) vs (Mut vs WT) comparison with either information sharing (x-axis) or no sharing (y-axis), and colored by stop codon. Generally, the observed loss-of-function variants including stops are much more significant (i.e. much more negative) in both clones when combining barcode information rather than keeping it separate.


    
![png](OCNT-DMSLIB-1-assay-run1_files/OCNT-DMSLIB-1-assay-run1_45_0.png)
    


### Visualization and Interpretation <a name="part3"></a>

A simple way to visualize summary statistics across a large variant set is with volcano plots, which show the log fold change on the x-axis and the -log10 p-value on the y-axis. Here, stop codons are still colored red:

    Warning message:
    “Removed 8 rows containing missing values (geom_point).”



    
![png](OCNT-DMSLIB-1-assay-run1_files/OCNT-DMSLIB-1-assay-run1_47_1.png)
    


We can extract and highlight only the difference plot for clone 8:

    Warning message:
    “Removed 2 rows containing missing values (geom_point).”



    
![png](OCNT-DMSLIB-1-assay-run1_files/OCNT-DMSLIB-1-assay-run1_49_1.png)
    


Alternatively, we can use heatmaps which show the Z-statistic for each comparison. The Z-statistic is the estimated effect size (here, it is the log fold change) divided by the standard error of that estimate. Note the scale - the limits are -30 and +10, which are _extremely_ significant effects:


    
![png](OCNT-DMSLIB-1-assay-run1_files/OCNT-DMSLIB-1-assay-run1_51_0.png)
    


We can also make the same plot using the log fold change values, for a different perspective on the same pattern:


    
![png](OCNT-DMSLIB-1-assay-run1_files/OCNT-DMSLIB-1-assay-run1_53_0.png)
    


This is broadly consistent with our expectations from the model, namely that we capture substantial additional power by using random intercepts, but that the underlying pattern is consistent even if we completely each dataset by treatment condition and then compare only after the fact.

We can also grab the conservation scores at the residue level from ConSurf, for comparison in this chunk:


    
![png](OCNT-DMSLIB-1-assay-run1_files/OCNT-DMSLIB-1-assay-run1_58_0.png)
    

