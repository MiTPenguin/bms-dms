## OCNT-DMSLIB-0-CRE DMS Assay Run #1: MC4R, Gs (DMS5)

| ID | Condition | Concentration | Replicates
| --- | ----------- | --- | ----------- |
| 1 | Forskolin | 2.5e-5 | 4 |
| 2 | None | 0 | 4 |
| 3 | aMSH | 5e-10 | 4 |
| 4 | aMSH | 2e-8 | 4 |
| 5 | THIQ | 4e-10 | 4 |
| 6 | THIQ | 1.2e-8 | 4 |
| 7 | aMSH | 5e-9 | 3 |
| 8 | THIQ | 4e-9 | 4 |

This report generates summary statistics and plots for MC4R, Gs (DMS5) identically to that of TYK2. Some sections that are TYK2-specific (e.g. drug resistence GoF, spike-ins) are omitted, and the ClinVar comparison is added (since, for MC4R, that association is much stronger than for TYK2).

1. [Barcode Sequencing Distributions](#part1)
2. [Inference and Stop Codon Effects](#part2)
3. [Visualizations](#part3)
4. [ClinVar](#part4)
5. [Unnormalized Summary Statistics](#part5)
6. [Comparing aMSH to THIQ](#part6)
7. [Published FunctionaL Variants](#part7)

### Barcode Sequencing Distributions <a name="part1"></a>

To start, we plot the distribution of unique barcodes per variant across samples:


    
![png](OCNT-DMSLIB-0-CRE-assay-run1_files/OCNT-DMSLIB-0-CRE-assay-run1_5_0.png)
    



    
![png](OCNT-DMSLIB-0-CRE-assay-run1_files/OCNT-DMSLIB-0-CRE-assay-run1_5_1.png)
    



    
![png](OCNT-DMSLIB-0-CRE-assay-run1_files/OCNT-DMSLIB-0-CRE-assay-run1_6_0.png)
    


We can stratify this by residue, but at the expense of generating very large plots. All full coverage plots are stored [here](./coverage-plots) for all samples, and a representative sample (1A) is shown below:


    
![png](OCNT-DMSLIB-0-CRE-assay-run1_files/OCNT-DMSLIB-0-CRE-assay-run1_9_0.png)
    



    
![png](OCNT-DMSLIB-0-CRE-assay-run1_files/OCNT-DMSLIB-0-CRE-assay-run1_10_0.png)
    



    
![png](OCNT-DMSLIB-0-CRE-assay-run1_files/OCNT-DMSLIB-0-CRE-assay-run1_11_0.png)
    


### Inference and Stop Codon Effects <a name="part2"></a>

Next, we can compute the normalizing contrasts against Forskolin or None, and plot the results. First, we can quickly check if there are any positions missing, and see that there are two partial positional dropouts at position 8 and 293 (meaning total coverage is still approx. 99.89%):


    
    
    |aa | missing at position|
    |:--|-------------------:|
    |C  |                   8|
    |F  |                   8|
    |W  |                   8|
    |*  |                   8|
    |Y  |                   8|
    |M  |                 293|
    |W  |                 293|


Now, we compute the normalizing contrasts and evaluate stop codon effects across each chunk:


```R
write_tsv(sumstats_all, "../sumstats/MC4R/MC4R-DMS5-Gs.tsv")
```


    
![png](OCNT-DMSLIB-0-CRE-assay-run1_files/OCNT-DMSLIB-0-CRE-assay-run1_18_0.png)
    



    
![png](OCNT-DMSLIB-0-CRE-assay-run1_files/OCNT-DMSLIB-0-CRE-assay-run1_18_1.png)
    


Finally, we can count up all significant effects at a 1% FDR in all comparisons. This verifies what is visually clear above, namely the relatively high basal activity (meaning significant variant effects in the absence of agonist): 


    
    
    |condition                  | Non-Significant| Significant (FDR < 0.01)|
    |:--------------------------|---------------:|------------------------:|
    |Forsk2.5e-05 - None0       |            5466|                     1167|
    |aMSH5e-10 - None0          |            4951|                     1682|
    |aMSH5e-09 - None0          |            5689|                      944|
    |aMSH2e-08 - None0          |            5698|                      935|
    |THIQ4e-10 - None0          |            5478|                     1155|
    |THIQ4e-09 - None0          |            5813|                      820|
    |THIQ1.2e-08 - None0        |            5790|                      843|
    |None0 - Forsk2.5e-05       |            5466|                     1167|
    |aMSH5e-10 - Forsk2.5e-05   |            4338|                     2295|
    |aMSH5e-09 - Forsk2.5e-05   |            5283|                     1350|
    |aMSH2e-08 - Forsk2.5e-05   |            5490|                     1143|
    |THIQ4e-10 - Forsk2.5e-05   |            4934|                     1699|
    |THIQ4e-09 - Forsk2.5e-05   |            5780|                      853|
    |THIQ1.2e-08 - Forsk2.5e-05 |            5938|                      695|


### Visualizations <a name="part3"></a>

Next, we generate the standard visualizations for these data, namely heatmaps (using log2FoldChange or Z-statistic) and volcano plots:


    
![png](OCNT-DMSLIB-0-CRE-assay-run1_files/OCNT-DMSLIB-0-CRE-assay-run1_23_0.png)
    



    
![png](OCNT-DMSLIB-0-CRE-assay-run1_files/OCNT-DMSLIB-0-CRE-assay-run1_23_1.png)
    



    
![png](OCNT-DMSLIB-0-CRE-assay-run1_files/OCNT-DMSLIB-0-CRE-assay-run1_24_0.png)
    



    
![png](OCNT-DMSLIB-0-CRE-assay-run1_files/OCNT-DMSLIB-0-CRE-assay-run1_24_1.png)
    


### ClinVar <a name="part4"></a>

For MC4R, we know the association between ClinVar pathogenicity status and the observed summary statistics are quite strong for the Gs (DMS5) dataset. Since we have it on hand, we can extract just those variants which are annotated in ClinVar under various categories, and plot them across all conditions. We could do this for any variant subset, and so show a few selections below:


    
![png](OCNT-DMSLIB-0-CRE-assay-run1_files/OCNT-DMSLIB-0-CRE-assay-run1_27_0.png)
    



    
![png](OCNT-DMSLIB-0-CRE-assay-run1_files/OCNT-DMSLIB-0-CRE-assay-run1_28_0.png)
    



    
![png](OCNT-DMSLIB-0-CRE-assay-run1_files/OCNT-DMSLIB-0-CRE-assay-run1_31_0.png)
    



    
![png](OCNT-DMSLIB-0-CRE-assay-run1_files/OCNT-DMSLIB-0-CRE-assay-run1_32_0.png)
    


### Unnormalized Summary Statistics <a name="part5"></a>

The summary statistics without any contrast or computed difference are located in the same directory as the normalized summary statistics [here](../sumstats/MC4R).

### Comparing aMSH to THIQ <a name="part6"></a>

We have three concentrations each for aMSH and THIQ. If we match these conditions together as "Low", "Medium", and "High" we can compute the pairwise difference in mutation effets between aMSH and THIQ at each level. This generates a new set of summary statistics, which are shown below. 


    
![png](OCNT-DMSLIB-0-CRE-assay-run1_files/OCNT-DMSLIB-0-CRE-assay-run1_39_0.png)
    



    
![png](OCNT-DMSLIB-0-CRE-assay-run1_files/OCNT-DMSLIB-0-CRE-assay-run1_39_1.png)
    



    
![png](OCNT-DMSLIB-0-CRE-assay-run1_files/OCNT-DMSLIB-0-CRE-assay-run1_39_2.png)
    


Most of these effects are negative, including a faint band of stops. This means that many variant effects are _more negative_ in aMSH compared to THIQ. This may reflect systematic differences, e.g. in the dosage of each compound. However, we can check which variants reverse this trend, namely those that are more negative in THIQ (instead of aMSH). Doing so, and imposing a FDR of 0.01, we find a resulting set of "THIQ-inhibiting" variants:


    
    
    | pos|chunk |aa |  estimate| std.error|     p.adj|condition             |
    |---:|:-----|:--|---------:|---------:|---------:|:---------------------|
    |  48|1     |D  | 0.8316088| 0.1657876| 0.0000664|aMSH5e-10 - THIQ4e-10 |
    | 129|3     |H  | 1.3480205| 0.3010955| 0.0005952|aMSH5e-10 - THIQ4e-10 |
    | 129|3     |N  | 1.0304938| 0.2719743| 0.0072895|aMSH5e-10 - THIQ4e-10 |
    | 129|3     |S  | 1.2465230| 0.1724862| 0.0000000|aMSH5e-10 - THIQ4e-10 |
    | 129|3     |T  | 1.0552811| 0.2210029| 0.0001882|aMSH5e-10 - THIQ4e-10 |
    | 284|6     |Y  | 0.6867885| 0.1752546| 0.0046356|aMSH5e-10 - THIQ4e-10 |
    | 129|3     |D  | 1.1923011| 0.2948011| 0.0029155|aMSH5e-09 - THIQ4e-09 |


These clearly consolidate to residue 129, which is a known interactor of THIQ based on a co-complex crystal structure. Residue 284 is also represented, but with a weaker effect compared to the residue 129 hits. It is difficult to identify variants which are "aMSH-inhibiting" in the same way, because _many_ variants appear "inhibitory" in aMSH if compared to THIQ. However, if we are conservative and look in the reverse direction, we resolve a very precise signal.


    
![png](OCNT-DMSLIB-0-CRE-assay-run1_files/OCNT-DMSLIB-0-CRE-assay-run1_44_0.png)
    



    
![png](OCNT-DMSLIB-0-CRE-assay-run1_files/OCNT-DMSLIB-0-CRE-assay-run1_46_0.png)
    


### Published Functional Variants <a name="part7"></a>

Using the normalized summary statistics, we can compare our summary statistics to those from Huang et al (2017):


    
![png](OCNT-DMSLIB-0-CRE-assay-run1_files/OCNT-DMSLIB-0-CRE-assay-run1_50_0.png)
    



    
![png](OCNT-DMSLIB-0-CRE-assay-run1_files/OCNT-DMSLIB-0-CRE-assay-run1_50_1.png)
    

