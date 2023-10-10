## Sequencing and Synthesis Errors in the TYK2 DMS Library


0. [Executive Summary](#part0)
1. [Global Error Rate Estimates](#part1)
2. [Distinguishing Synthesis from Sequencing Errors](#part2)
3. [Error Distributions by Position](#part3)
4. [Error Distributions by Base Identity](#part4)

### Executive Summary <a name="part0"></a>

- This analysis of sequencing error rates uses barcode mapping data from TYK2 Chunks 3-6, sequencing run `VH00964_8_AAAVVT7HV`
- The median proportion of reads which exactly match an expected oligo is 66.5% across barcode mapping samples
    - Thus, the median combined error rate (including synthesis, PCR, and sequencing) is 33.5%
- Substitutions are highly correlated across all samples in the same chunk, while indels are much more correlated within rather than across A/B group
- Errors (including substitutions and indels) are uniformly distributed along the length of the oligo
- Most substitutions are likely sequencing errors (very low base calling quality score)
  - For example, if calling genetic variants in a human, base qualities this low would not be used
  - However, the most common substitution in _both_ high and low quality bases is G>T
  - In high quality basecalls only, the second most abundant substitution is C>T

### Global Error Rate Estimates <a name="part1"></a>

For the main DMS experiment, we want to use only error-free, valid oligo-barcode associations from these data and we apply a series of filters to define what "valid" means. However, it is also useful to characterize the actual content of the removed barcode mapping data. Each barcode mapping sample generates a set of paired-end sequences, which are classified into the bins shown below. These bins are _non-overlapping_ such that `total fragments` equals the sum of all other fields for a given sample:


    
    
    |sample | total fragments| uncombined| empty| no assignment| ambiguous assignment| unique assignment with errors| unique assignment with no errors|
    |:------|---------------:|----------:|-----:|-------------:|--------------------:|-----------------------------:|--------------------------------:|
    |3A_1   |        56995569|    1467442|     0|       1088960|              3550136|                      11830789|                         39058242|
    |3A_2   |        57329625|    1503360|     2|       1215396|              3629955|                      10894684|                         40086228|
    |3B_1   |        48504137|    1266729|     1|        931425|              2998564|                       6165102|                         37142316|
    |3B_2   |        41862423|    1109183|     0|        895508|              2692060|                      12800230|                         24365442|
    |4A_1   |        56866293|    1462627|     2|       1049978|              3019290|                      11464951|                         39869445|
    |4A_2   |        57735109|    1393315|     4|       1061637|              3054369|                      11344057|                         40881727|
    |4B_1   |        57669587|    1456096|     0|       1080222|              3067813|                      12322100|                         39743356|
    |4B_2   |        51393103|    1323227|     1|        946596|              2712739|                      12885545|                         33524995|
    |5A_1   |        49946515|    1312957|     0|       1375675|              3471511|                       9931246|                         33855126|
    |5A_2   |        50238182|    1282210|     0|       1327774|              3404602|                      12855552|                         31368044|
    |5B_1   |        50381602|    1102999|     0|       1244888|              3276753|                      12155314|                         32601648|
    |5B_2   |        43417816|     863565|     1|        932892|              2761650|                      13042989|                         25816719|
    |6A_1   |        43311834|    1239324|     1|        921196|              2939123|                      10093106|                         28119084|
    |6A_2   |        41659294|    1211553|     1|        914082|              2872745|                      12638482|                         24022431|
    |6B_1   |        26367894|     711296|     1|        567254|              1806061|                      10015693|                         13267589|
    |6B_2   |        49243114|    1394180|     0|       1079673|              3346455|                       9638255|                         33784551|


And the same table as a rounded percentage:


    
    
    |sample |total fragments |uncombined |empty     |no assignment |ambiguous assignment |unique assignment with errors |unique assignment with no errors |
    |:------|:---------------|:----------|:---------|:-------------|:--------------------|:-----------------------------|:--------------------------------|
    |3A_1   |100             |2.575      |0         |1.911         |6.229                |20.76                         |68.53                            |
    |3A_2   |100             |2.622      |3.489e-06 |2.12          |6.332                |19                            |69.92                            |
    |3B_1   |100             |2.612      |2.062e-06 |1.92          |6.182                |12.71                         |76.58                            |
    |3B_2   |100             |2.65       |0         |2.139         |6.431                |30.58                         |58.2                             |
    |4A_1   |100             |2.572      |3.517e-06 |1.846         |5.309                |20.16                         |70.11                            |
    |4A_2   |100             |2.413      |6.928e-06 |1.839         |5.29                 |19.65                         |70.81                            |
    |4B_1   |100             |2.525      |0         |1.873         |5.32                 |21.37                         |68.92                            |
    |4B_2   |100             |2.575      |1.946e-06 |1.842         |5.278                |25.07                         |65.23                            |
    |5A_1   |100             |2.629      |0         |2.754         |6.95                 |19.88                         |67.78                            |
    |5A_2   |100             |2.552      |0         |2.643         |6.777                |25.59                         |62.44                            |
    |5B_1   |100             |2.189      |0         |2.471         |6.504                |24.13                         |64.71                            |
    |5B_2   |100             |1.989      |2.303e-06 |2.149         |6.361                |30.04                         |59.46                            |
    |6A_1   |100             |2.861      |2.309e-06 |2.127         |6.786                |23.3                          |64.92                            |
    |6A_2   |100             |2.908      |2.4e-06   |2.194         |6.896                |30.34                         |57.66                            |
    |6B_1   |100             |2.698      |3.792e-06 |2.151         |6.849                |37.98                         |50.32                            |
    |6B_2   |100             |2.831      |0         |2.193         |6.796                |19.57                         |68.61                            |


Under this framework, the fragments that are `uncombined`, `empty`, or have `no assignment` are the remaining dark matter of the dataset - they would require further processing to characterize. This is definitely possible - some manual exploration suggests that virtually all reads (>99%) can be aligned to TYK2, but simply cannot be meaningfully resolved into our expected oligos.

The remaining reads can be uniquely identified as originating from a particular oligo, and are either imperfect or exact matches. The median exact match frequncy is 66.5%, meaning that approx. 2/3 of the raw read count from a given paired end data set are exact copies of an expected oligo sequence.

### Distinguishing Synthesis from Sequencing Errors <a name="part2"></a>

While useful, the above does not distinguish between squencing and synthesis errors. Fortunately, we have replicates which should differ in their sequencing, but not synthesis, errors although the difference is probabilistic rather than absolute. One measure of similarity is the Jaccard index, which quantifies the overlap between two groups of errors. It is defined as the intersection of the errors detected in two samples divided by the union of those errors - a value of 1 means that all detected errors are shared, while a value of zero means that no detected variants are shared. A "detected" variant means at least one read - we aren't yet considering depth. We can show this index as a heatmap, separated by indels and substitutions:


    
![png](TYK2-Synth-Seq-Errors_files/TYK2-Synth-Seq-Errors_14_0.png)
    



    
![png](TYK2-Synth-Seq-Errors_files/TYK2-Synth-Seq-Errors_15_0.png)
    


Samples within the same chunk share 75-80% of substitution errors regardless of whether they are group A or B. In contrast, indel error sharing is much higher within A/B group in each chunk (50-55% vs ~35%).

Another informative view that incorporates the abundance of each error is the Pearson correlation of the log counts across errors in each sample. For replicate pairs, the observed errors (both indels and substitutions) are very well correlated:

    Using sample as id variables
    



    
![png](TYK2-Synth-Seq-Errors_files/TYK2-Synth-Seq-Errors_19_0.png)
    



    
![png](TYK2-Synth-Seq-Errors_files/TYK2-Synth-Seq-Errors_21_0.png)
    


### Error Distributions by Position <a name="part3"></a>

In the above, we considered overlap between indels and substitutions in two samples, regardless of the frequency of each modification in either sample. Alternatively, we can instead ignore the particular oligo of origin, and pile up the numbers of reads with a substitution or indel at each position along the 210bp oligo length. This shows no obvious patterns or positional enrichment for any error category:


    
![png](TYK2-Synth-Seq-Errors_files/TYK2-Synth-Seq-Errors_23_0.png)
    



    
![png](TYK2-Synth-Seq-Errors_files/TYK2-Synth-Seq-Errors_24_0.png)
    


Note that in the indel plot, there is a gap towards the ends followed by a spike. This is induced by the aligner - it does not effectively resolve indels located in the first or last 6 bases of the fragment, and so just chops off the end entirely. In those cases, the reads are represented by the spike at each end, depending on whether it was a 5-prime or 3-prime truncation. In all other cases, the shown indel position is the **right-most** position of the deleted or inserted region, which is why the right-side gap is a bit more filled in than the left-side.

### Error Distributions by Base Identity <a name="part4"></a>

Since there is not much positional specificity, we can instead look at the base composition stratified by position. For indels this is tricky since indels have variable length and involve multiple bases, but for substitutions it is straightforward to count how many reads contain each of the possible single base changes:


    
![png](TYK2-Synth-Seq-Errors_files/TYK2-Synth-Seq-Errors_28_0.png)
    

In all samples the most frequent substitution is G>T or C>A, and presumably this contains at least some of both sequencing and synthesis errors. An additional piece of information is the base quality score which quantifies the probability of a basecalling error. The distribution of base qualities for substituted bases in each sample is fairly consistent, and shows that (in principle) most of these substitutions are explicable as low-quality base calls:

    
![png](TYK2-Synth-Seq-Errors_files/TYK2-Synth-Seq-Errors_31_0.png)
    


We can use the base quality to separate out substitutions that are likely sequencing or basecalling errors from substitutions which is unlikely to be a sequencing error (according to the basecalling algorithm). This changes the error spectrum substantially:


    
![png](TYK2-Synth-Seq-Errors_files/TYK2-Synth-Seq-Errors_34_0.png)
    

