## DMS Barcode Resampling Simulation Prototype

0. [Core Conclusions](#part0)
1. [Background](#part1)
2. [Data Subset](#part2)
3. [Sampling Across Clones and Replicates](#part3)
4. [Downsampling Results and Interpretation](#part4)

### Core Conclusions <a name="part0"></a>

NOTE: An important caveat is that the following analysis considers only all barcodes from the TYK2 chunk 10 assay that are either WT or contain a variant at position 638. This analysis could be expanded to all positions and variants.

- We observe that relative library composition at the variant level remains relatively constant throughout various stages in the process.
    -  For example, variants that are 10X more abundant than the rest of the library at the chip design stage are 10X more abundant at the barcode mapping _and_ barcode expression stages.
-  A large majority of mapped barcodes are "lost" between mapping and expression, presumably during integration.
-  Most barcodes in each RNA sample are detected in more than one replicate, and a defined subset of barcodes is specific to each clone. Clone 1 has more clone-specific barcodes than Clone 8, which likely explains its slightly increased power.
-  Resampling of expression counts under all weighting schemes demonstrates a straightforward tradeoff between barcodes per variant and power very similar to that observed previously with MC4R.

### Background <a name="part1"></a>

In the initial TYK2 DMS across chunk 10, two clone pools were generated after barcode mapping from the same library (named c1 and c8). After creation and expansion, they are frozen down, stored, and re-expanded in triplicate for the barcode expression component of the experiment. The schematic below depicts this procedure:


    
![png](Barcode-Sampling-Proto_files/Barcode-Sampling-Proto_4_0.png)
    


There are several points at which we sample a subpopulation of variants from a larger group. However, we have two readouts - the barcode mapping procedure and the barcode expression assay. While barcode mapping provides some quantitative information about barcode abundance (e.g. the number of reads for each barcode), it is primarily used for determining barcode-oligo associations. The molecular processing and storage after barcode mapping may change barcode frequencies relative to those displayed during barcode mapping, and thus read counts at barcode mapping may not provide quantitative information about sampling probabilities.

### Data Subset <a name="part2"></a>

The "unit of regression" in DMS analysis includes all wild-type barcodes and all barcodes mapped to a variant at a specific, individual position. With the benefit of our spike-ins, we can choose a position with an enriched variant to capture resampling effects at different scales. We select position 638, which includes the enriched stop variant 638X. Retrieving all barcodes in either barcode mapping or in barcode expression corresponding to WT or 638-MUT constructs returns a tractable dataset for resampling:

- 475102 barcodes identified during barcode mapping
- 62244 out of 475102 barcodes were detected in at least one replicate during barcode expression


    
![png](Barcode-Sampling-Proto_files/Barcode-Sampling-Proto_10_0.png)
    


These look pretty well correlated, both across replicates and between the barcode map and the replicates. These are small sets of values so we can compute the correlation matrix easily:


    
![png](Barcode-Sampling-Proto_files/Barcode-Sampling-Proto_12_0.png)
    


The barcodes detected per residue are remarkably consistent. There are two clearly outlying points though, which are the WT barcodes and 638X barcodes. We can remove those and plot this again to get a clearer view, and also to double check that these values are indeed different across samples:


    
![png](Barcode-Sampling-Proto_files/Barcode-Sampling-Proto_14_0.png)
    


We could imagine subsampling either the barcode map (followed by removing and counts from barcode expression which were from a removed barcode) or the barcode expression data. Given the ambiguity about the extent to which barcode counts in barcode mapping are quantitatively useful, we start by simply subsampling barcodes from each of the three replicates of clone 1. We can do this either in a _weighted_ manner, where the read counts themselves provide the relative probability of selection, or _unweighted_ where those probabilities are all equal.

### Sampling Across Clones and Replicates <a name="part3"></a>

Above, we describe the numbers of barcodes per variant at this position compared to WT. However, we can more directly examing the overlap of barcodes between each sample and the barcode map using an upset plot. To see how to interpret these plots, lets make one of these plots that includes just three replicates from clone 1:


    
![png](Barcode-Sampling-Proto_files/Barcode-Sampling-Proto_19_0.png)
    


Considering samples 1A, 1B, and 1C, there are seven options for each barcode. It might be in all samples, only one sample, or more than one but not all. We can count up the number of barcodes in each of these categories, where each sample is a row whose sum is shown in the horizontal left bars, while each column is a non-overlapping barcode set detected in those samples with black dots.

Now, let's see what it looks like for all samples, including or excluding the barcode map:


    
![png](Barcode-Sampling-Proto_files/Barcode-Sampling-Proto_22_0.png)
    


This is hard to see, so we can remove the "bcmap" set and just view the rest:


    
![png](Barcode-Sampling-Proto_files/Barcode-Sampling-Proto_24_0.png)
    


It immediately shows that most variants from barcode mapping are never seen again, but those that are seen again are a mix of those shared across all samples (both clones), those restricted to only one clone, and those that are specific to one or a handful of samples:

### Downsampling Results and Interpretation <a name="part4"></a>

We can re-sample this test dataset in a very large number of ways:

- [Unweighted Resampling](#partA)
- [Expression-Weighted Resampling](#partB)
- [Mapping-Weighted Resampling](#partC)

The easiest way to understand how this works is a simple diagram. Imagine there are a few reads that have the following counts in barcode mapping and expression:

| Barcode    | Mapping Count | Expression Count |
| ---------- | ----------- | ----------- |  
| AACAAGTCTATATTCAATATT   | 5       | 15 | 
| TTTGTGTTTGTGTTTATCAAA   | 10        | 5 |  
| ... | ... | ... |

  -  In Unweighted resampling, we randomly select with replacement from the `Expression Count` column some number of times, where each row has the same probability of selection.
  -  In Mapping-Weighted resampling, each row has a probability of selection proportional to the `Mapping Count`
  -  In Expression-Weighted resampling, each row has a probability of selection proportional to the `Expression Count`


#### Unweighted Resampling <a name="partA"></a>

We consider the vector of barcode expression counts across all samples, and resample these counts evenly such that any count has an equal probability of being sampled, with no other modifications (yet). We then compute summary statistics and plot their distributions stratified by resampling size:


    
![png](Barcode-Sampling-Proto_files/Barcode-Sampling-Proto_28_0.png)
    


#### Expression-Weighted Resampling <a name="partB"></a>

The simplest way to weight each barcode count is according to the normalized count frequency itself. So, barcode counts that are larger are proportionally more likely to be sampled than those that are smaller. This has a visible but relatively small effect on the resulting summary statistics:


    
![png](Barcode-Sampling-Proto_files/Barcode-Sampling-Proto_30_0.png)
    


#### Mapping-Weighted Resampling <a name="partC"></a>

Alternatively, we can specify the probability of each barcode count being included as proportional to its frequency during _barcode mapping_ rather than expression itself. This may better model the experimental schematic shown at the start of this document, where many of the bottlenecking and cell sampling steps are in between barcode mapping and expression. We can retrieve these frequencies from the barcode mapping data and examine them:


    
![png](Barcode-Sampling-Proto_files/Barcode-Sampling-Proto_33_0.png)
    



    
![png](Barcode-Sampling-Proto_files/Barcode-Sampling-Proto_33_1.png)
    


They look pretty uncorrelated - to see what the sampling effects of using the barcode mapping count are likely to be, we can look at the one-dimensional distributions. The red lines below indicate read counts of 8 and 16. Those are the boundaries of the middle 50% of barcode frequencies from mapping for each sample, and there are no absurd outliers at the barcode-mapping level. So, we might expect this to be pretty similar to uniform sampling:


    
![png](Barcode-Sampling-Proto_files/Barcode-Sampling-Proto_35_0.png)
    



    
![png](Barcode-Sampling-Proto_files/Barcode-Sampling-Proto_35_1.png)
    


Indeed, the summary statistic distributions by resampling size look pretty similar to those that came before:


    
![png](Barcode-Sampling-Proto_files/Barcode-Sampling-Proto_37_0.png)
    

